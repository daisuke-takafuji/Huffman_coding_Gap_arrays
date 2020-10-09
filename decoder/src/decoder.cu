/*
 *  Decode files using Huffman encoding with gap arrays
 *  Usage: ./ghuffman_decoder inputfile outputfile 
 *  Copyright (C) 2020  Daisuke Takafuji
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/cu_timer.h"
#include "../include/constants.hpp"
#include "../include/decoder.cuh"

#define CUERROR  {\
cudaError_t cuError;\
	if( cudaSuccess != (cuError = cudaGetLastError()) ){\
		printf("Error: %s : at %s line %d\n", cudaGetErrorString(cuError), __FILE__, __LINE__);\
		exit(EXIT_FAILURE);\
	} \
};

__device__ __host__ void show_bits( unsigned int bits, unsigned char length ){
	unsigned int mask = 1<<(length-1);
	for(int i=0; i<length; i++){
		if( mask&bits ) printf("1");
		else printf("0");
		mask >>= 1;
	}
	printf("\n");
}

__global__ void nullkernel(){
	return;
}

__global__ void cu_make_table(
		struct Symbol *symbols,
		struct TableInfo *tableinfo,
		int symbol_count,
		int prefix_bit,
		void *decode_table ){

	__shared__ uint4 length_sum[MAX_CODE_NUM];

	if( threadIdx.x < symbol_count ){
		struct Symbol this_symbol = symbols[threadIdx.x];
		int codelength = this_symbol.length;
		int shift = 8*(codelength - 1);
		int pos = shift/32;
		shift -= pos*32;
		uint4 length = make_uint4(0,0,0,0);
		unsigned int *this_length;
		this_length = (unsigned int *)&length;

		this_length[pos] = 1<<(shift);
		length_sum[threadIdx.x] = length;
		
		unsigned int code = 0;

		__syncthreads();

		for(int i=1; i<MAX_CODE_NUM; i<<=1){
			if( (int)threadIdx.x-i>=0 ){
				length_sum[threadIdx.x].x += length_sum[threadIdx.x-i].x;
				length_sum[threadIdx.x].y += length_sum[threadIdx.x-i].y;
				length_sum[threadIdx.x].z += length_sum[threadIdx.x-i].z;
				length_sum[threadIdx.x].w += length_sum[threadIdx.x-i].w;
			}
			__syncthreads();
		}
		
		if( threadIdx.x > 0 ){
			int mul=1;
			length = length_sum[threadIdx.x - 1];
			while( pos >= 0 ){
				unsigned int current_value = this_length[pos];
				while(shift>=0){
					unsigned int mask = (0xFF)<<shift;
					 code += ((mask&current_value)>>shift)*mul;
					 mul<<=1;
					 shift -= 8;
				}
				shift = 24;
				pos--;
			}
		}
	}
}

void cu_make_table_call(
		struct Symbol *symbols,
		struct TableInfo *tableinfo,
		int symbol_count,
		int prefix_bit,
		void *decode_table ){

	double nullsum = 0,
		   testsum = 0;
	for(int i=0; i<LOOP*2; i++){
	TIMER_START(null)
		nullkernel<<<1,256>>>();
	TIMER_STOP(null)

	TIMER_START(test)
	cu_make_table<<<1, 256>>>(
				symbols,
				tableinfo, 
				symbol_count,
				prefix_bit,
				decode_table);
	TIMER_STOP(test)
		if( i>=LOOP ){
			nullsum += timenull;
			testsum += timetest;
		}
	}

	printf("code_gentime, %lf, %lf, %lf\n", nullsum*1000/LOOP, testsum*1000/LOOP, (testsum - nullsum)*1000/LOOP);
	CUERROR;
}

// decode with threads 
__global__ void gpu_dec_l1_l2_multi(
		unsigned int* gpu_input, 
		int inputfilesize, 
		void* dectable,
		volatile unsigned long long int* output_pos, 
		unsigned int* gaparray, 
		int gap_element_num,
		int gaparray_size,
		unsigned int* gpu_output,
		int outputfilesize, 
		int* counter,
		int tablesize,
		unsigned int prefix_bit,
		unsigned int symbol_count,
		struct TableInfo tableinfo){

	unsigned int* thread_input;
	int tinIdx = 0;
	int inputSize = SEGMENT_SIZE/MAX_BITS;

	int code_count = 0;
	int bits_count = 0;
	int at = 0;
	int first_at = 0;
	
	extern __shared__ unsigned char table[];
	__shared__ unsigned long long int exclusive_sum;
	__shared__ unsigned char local_gap[NUM_THREADS*LOCAL_SEGMENT_NUM];
	__shared__ unsigned int local_output_pos[NUM_THREADS*LOCAL_SEGMENT_NUM];
	__shared__ unsigned int blockprefixsum[NUM_THREADS/WARP_SIZE];

	unsigned char *dectable_tmp;
	unsigned int *ptr_table;
	unsigned char *l1_table, *l2_table, *length_table;
	const unsigned int boundary_code = tableinfo.l1table_size;

	unsigned int local_bits_count=0;
	unsigned int local_gap_idx = 1;

	dectable_tmp = (unsigned char *)dectable;
	ptr_table = (unsigned int *)table;
	length_table = (unsigned char *)table + tableinfo.ptrtable_size*sizeof(unsigned int);
	l1_table = length_table + MAX_CODE_NUM;
	l2_table = l1_table + tableinfo.l1table_size;
	for( int i=threadIdx.x; i<tablesize; i+=blockDim.x ){
		table[i] = dectable_tmp[i];
	}
	__shared__ int block_counter;
	int block_pos;
	if(threadIdx.x == 0) block_counter = atomicAdd(counter, 1);
	__syncthreads();

	block_pos = block_counter;

	int input_ptr = block_pos*blockDim.x + threadIdx.x;
	if( input_ptr < gap_element_num ){
		thread_input = (gpu_input + input_ptr*inputSize + gaparray_size);
		if( input_ptr != 0 ){
			first_at = at = ((gaparray[(input_ptr-1)/GAP_FAC_NUM]>>(((input_ptr-1)&7)*GAP_LENGTH_MAX)))&0xF;
		}
		local_gap[threadIdx.x*LOCAL_SEGMENT_NUM] = first_at;
		input_ptr *= inputSize;

		const unsigned int prefix_mask = ~(((unsigned int) (0) - 1) >> prefix_bit);
		const unsigned int prefix_shift = MAX_BITS - prefix_bit;

		unsigned int window = 0;
		unsigned int next = 0;
		unsigned int copy_next = 0;

		unsigned int taken=0;

		bits_count = at;
		local_bits_count = bits_count;

		window = thread_input[0];
		next = thread_input[1];
		copy_next = next;
		copy_next >>= MAX_BITS - at;
		next <<= at;
		window <<= at;
		window |= copy_next;
		while( bits_count < SEGMENT_SIZE && input_ptr < inputfilesize ){
			while( at < MAX_BITS && bits_count < SEGMENT_SIZE ){
				unsigned int deccode = (window&prefix_mask)>>prefix_shift;
				if( deccode < boundary_code ){
					unsigned char current_symbol = l1_table[deccode];
					taken = length_table[current_symbol];
				}else{
					unsigned int ptrIdx = deccode - boundary_code;
					unsigned int table_Idx = ptr_table[ptrIdx];
					unsigned int table_width = ((table_Idx>>16)&0xFFFF);
					table_Idx &= 0xFFFF;
					unsigned int l2_mask = (~(((unsigned int) (0) - 1) >> table_width))>>prefix_bit;
					unsigned int l2_shift = MAX_BITS - prefix_bit - table_width;

					deccode = (window&l2_mask)>>l2_shift;
					unsigned char current_symbol = l2_table[ table_Idx + deccode ];
					taken = length_table[current_symbol];
				}

				// store local gaps
				if( local_bits_count >= LOCAL_SEGMENT_SIZE ){
					local_bits_count -= LOCAL_SEGMENT_SIZE;
					local_gap[threadIdx.x*LOCAL_SEGMENT_NUM + local_gap_idx] = at&(MAX_CODEWORD_LENGTH-1);
					local_output_pos[threadIdx.x*LOCAL_SEGMENT_NUM + local_gap_idx] = code_count;
					local_gap_idx++;
				}

				// shift by decoded code bits
				copy_next = next;
				copy_next >>= MAX_BITS - taken;

				next <<= taken;
				window <<= taken;
				at += taken;
				bits_count += taken;
				local_bits_count += taken;
				window |= copy_next;
				code_count++;
			}
			input_ptr++;
			tinIdx++;
			window = thread_input[tinIdx];
			next = thread_input[tinIdx+1];

			if( at == MAX_BITS ){
				at = 0;
			}else{
				at -= MAX_BITS;
				copy_next = next;
				window <<= at;
				next <<= at;

				copy_next >>= MAX_BITS - at;
				window |= copy_next;
			}
		}
	}

	// compute prefixsums
	int warpLane = threadIdx.x&(WARP_SIZE-1);
	int warpId = threadIdx.x/WARP_SIZE;
	unsigned int tmp_count = code_count;
	unsigned int tmp_value = 0;
	for( unsigned int delta = 1; delta < WARP_SIZE; delta <<= 1 ){
		tmp_value = __shfl_up_sync(0xFFFFFFFF, tmp_count, delta, WARP_SIZE);
		if( warpLane >= delta ) tmp_count += tmp_value;
	}
	// compute prefixsums in blocks
	if( warpLane == WARP_SIZE-1 ){
		blockprefixsum[warpId] = tmp_count;
	}
	__syncthreads();

	if( threadIdx.x < BLOCK_LAST ){
		tmp_value = blockprefixsum[threadIdx.x];
		const unsigned int shfl_mask = ~( (~0) << BLOCK_LAST );
		for( unsigned int delta=1; delta < BLOCK_LAST; delta <<=1 ){
			unsigned int tmp = __shfl_up_sync( shfl_mask, tmp_value, delta, BLOCK_LAST);
			if( threadIdx.x >= delta ) tmp_value += tmp;
		}
		blockprefixsum[threadIdx.x] = tmp_value;
	}
	exclusive_sum = 0;
	__syncthreads();

	// One warp looks back to 32 block 
	if( warpId == 0 ){
		int posIdx = block_pos - warpLane;
		unsigned long long int local_sum;
		unsigned long long int before_sum = 0;
		if( warpLane == 0 ){
			local_sum = blockprefixsum[NUM_THREADS/WARP_SIZE-1];
			if( block_pos == 0 ){
				output_pos[block_pos] = local_sum | FLAG_P;
			}else{
				output_pos[block_pos] = local_sum | FLAG_A;
			}
		}

		if( block_pos > 0 ){
			while( 1 ){
				while( posIdx > 0 && (before_sum == 0)){
					before_sum = output_pos[posIdx-1];
				}
				unsigned long long int tmp_sum = 0;
				for(unsigned int delta = 1; delta < WARP_SIZE; delta <<= 1){
					tmp_sum = __shfl_down_sync(0xFFFFFFFF, before_sum, delta, WARP_SIZE);
					if( warpLane < WARP_SIZE - delta && ((before_sum&FLAG_P)==0) ) before_sum += tmp_sum;
				}
				local_sum += (before_sum&(~FLAG_MASK));
				before_sum = __shfl_sync(0xFFFFFFFF, before_sum, 0);
	
				if( before_sum&FLAG_P ){
					break;
				}
				posIdx -= WARP_SIZE;
				before_sum = 0;
			}
			if( warpLane == 0 ){
				output_pos[block_pos] = local_sum | FLAG_P;
				exclusive_sum = local_sum - blockprefixsum[BLOCK_LAST - 1];
			}
		}
	}

	__syncthreads();

	unsigned long long int before_sum = 0;
	int tmp_count_plus = 0;
	if(warpLane == 0){
		before_sum = exclusive_sum;
	}
	if(warpId > 0 && warpLane == 0){
		before_sum = exclusive_sum;
		tmp_count_plus = blockprefixsum[warpId-1];
	}
	before_sum = __shfl_sync(0xFFFFFFFF, before_sum, 0);
	tmp_count += __shfl_sync(0xFFFFFFFF, tmp_count_plus, 0);

	tmp_count_plus = before_sum + tmp_count - code_count;
	local_output_pos[threadIdx.x*LOCAL_SEGMENT_NUM] = tmp_count_plus;
	for( int i=1; i<LOCAL_SEGMENT_NUM; i++ ){
		local_output_pos[threadIdx.x*LOCAL_SEGMENT_NUM + i] += tmp_count_plus;
	}
	__syncthreads();

	input_ptr = block_pos*blockDim.x*inputSize;
	int local_segIdx = 0;
	unsigned int local_inputSize = LOCAL_SEGMENT_SIZE/MAX_BITS;
	unsigned int local_input_ptr = input_ptr + threadIdx.x*local_inputSize;

	while( local_input_ptr < inputfilesize && local_segIdx < LOCAL_SEGMENT_NUM ){
		thread_input = (gpu_input + local_input_ptr + gaparray_size);
		unsigned int output_ptr = local_output_pos[threadIdx.x + NUM_THREADS*local_segIdx];
		unsigned int output_ptr_idx = output_ptr&3;
		output_ptr >>= 2;
		const unsigned int prefix_mask = ~(((unsigned int) (0) - 1) >> prefix_bit);
		const unsigned int prefix_shift = MAX_BITS - prefix_bit;

		unsigned int window = 0;
		unsigned int next = 0;
		unsigned int copy_next = 0;

		unsigned int taken=0;

		at = local_gap[threadIdx.x + NUM_THREADS*local_segIdx];
		bits_count = at;
		tinIdx = 0;

		window = thread_input[0];
		next = thread_input[1];
		copy_next = next;
		copy_next >>= MAX_BITS - at;
		next <<= at;
		window <<= at;
		window |= copy_next;

		unsigned int output_symbols = 0;
		int first_flag=1;
		unsigned char decsymbol;
		while( bits_count < LOCAL_SEGMENT_SIZE && local_input_ptr < inputfilesize ){
			while( at < MAX_BITS && bits_count < LOCAL_SEGMENT_SIZE ){
				unsigned int deccode = (window&prefix_mask)>>prefix_shift;
				if( deccode < boundary_code ){
					decsymbol = l1_table[deccode];
					taken = length_table[decsymbol];
				}else{
					unsigned int ptrIdx = deccode - boundary_code;
					unsigned int table_Idx = ptr_table[ptrIdx];
					unsigned int table_width = ((table_Idx>>16)&0xFFFF);
					table_Idx &= 0xFFFF;
					unsigned int l2_mask = (~(((unsigned int) (0) - 1) >> table_width))>>prefix_bit;
					unsigned int l2_shift = MAX_BITS - prefix_bit - table_width;

					deccode = (window&l2_mask)>>l2_shift;
					decsymbol = l2_table[ table_Idx + deccode ];
					taken = length_table[decsymbol];
				}

				copy_next = next;
				copy_next >>= MAX_BITS - taken;
	
				next <<= taken;
				window <<= taken;
				at += taken;
				bits_count += taken;
				window |= copy_next;

				// output decoded symbols 
				UINT_OUT(output_symbols, decsymbol, output_ptr_idx);
				output_ptr_idx++;
				if( (output_ptr_idx&3)== 0 ){
					if( first_flag ){
						first_flag = 0;
						atomicOr( &gpu_output[output_ptr++], output_symbols );
					}else{
						gpu_output[output_ptr++] = output_symbols;
					}
					output_symbols=0;
					output_ptr_idx=0;
				}
			}
			local_input_ptr++;
			tinIdx++;

			window = thread_input[tinIdx];
			next = thread_input[tinIdx+1];

			if( at == MAX_BITS ){
				at = 0;
			}else{
				at -= MAX_BITS;
				copy_next = next;
				window <<= at;
				next <<= at;

				copy_next >>= MAX_BITS - at;
				window |= copy_next;
			}
		}
		// output remained symbols
		if( output_ptr_idx != 0 ){
			atomicOr( &gpu_output[output_ptr], output_symbols );
		}
		local_segIdx++;
		local_input_ptr = input_ptr + (threadIdx.x + NUM_THREADS*local_segIdx)*local_inputSize;
	}
}

__global__ void gpu_dec_l1_l2(
		unsigned int* gpu_input, 
		int inputfilesize, 
		void* dectable,
		volatile unsigned long long int* output_pos, 
		unsigned int* gaparray, 
		int gap_element_num,
		int gaparray_size,
		unsigned int* gpu_output,
		int outputfilesize, 
		int* counter,
		int tablesize,
		unsigned int prefix_bit,
		unsigned int symbol_count,
		struct TableInfo tableinfo){

	unsigned int* thread_input;
	int tinIdx = 0;
	int inputSize = SEGMENT_SIZE/MAX_BITS;

	int code_count = 0;
	int bits_count = 0;
	int at = 0;
	int first_at = 0;

	extern __shared__ unsigned char table[];
	__shared__ unsigned long long int exclusive_sum;

	unsigned char *dectable_tmp;
	unsigned int *ptr_table;
	unsigned char *l1_table, *l2_table, *length_table;
	const unsigned int boundary_code = tableinfo.l1table_size;
	dectable_tmp = (unsigned char *)dectable;
	ptr_table = (unsigned int *)table;
	length_table = (unsigned char *)table + tableinfo.ptrtable_size*sizeof(unsigned int);
	l1_table = length_table + MAX_CODE_NUM;
	l2_table = l1_table + tableinfo.l1table_size;
	for( int i=threadIdx.x; i<tablesize; i+=blockDim.x ){
		table[i] = dectable_tmp[i];
	}
	__shared__ int block_counter;
	int block_pos;
	if(threadIdx.x == 0) block_counter = atomicAdd(counter, 1);
	__syncthreads();

	block_pos = block_counter;

	int input_ptr = block_pos*blockDim.x + threadIdx.x;

	if( input_ptr < gap_element_num ){
		thread_input = (gpu_input + input_ptr*inputSize + gaparray_size);
		if( input_ptr != 0 ){
			first_at = at = ((gaparray[(input_ptr-1)/GAP_FAC_NUM]>>(((input_ptr-1)&7)*GAP_LENGTH_MAX)))&0xF;
		}
		input_ptr *= inputSize;

		const unsigned int prefix_mask = ~(((unsigned int) (0) - 1) >> prefix_bit);
		const unsigned int prefix_shift = MAX_BITS - prefix_bit;

		unsigned int window = 0;
		unsigned int next = 0;
		unsigned int copy_next = 0;

		unsigned int taken=0;

		bits_count = at;

		window = thread_input[0];
		next = thread_input[1];
		copy_next = next;
		copy_next >>= MAX_BITS - at;
		next <<= at;
		window <<= at;
		window |= copy_next;

		while( bits_count < SEGMENT_SIZE && input_ptr < inputfilesize ){
			while( at < MAX_BITS && bits_count < SEGMENT_SIZE ){
				unsigned int deccode = (window&prefix_mask)>>prefix_shift;
				if( deccode < boundary_code ) {
					unsigned char current_symbol = l1_table[deccode];
					taken = length_table[current_symbol];
				}else{
					unsigned int ptrIdx = deccode - boundary_code;
					unsigned int table_Idx = ptr_table[ptrIdx];
					unsigned int table_width = ((table_Idx>>16)&0xFFFF);
					table_Idx &= 0xFFFF;
					unsigned int l2_mask = (~(((unsigned int) (0) - 1) >> table_width))>>prefix_bit;
					unsigned int l2_shift = MAX_BITS - prefix_bit - table_width;
	
					deccode = (window&l2_mask)>>l2_shift;
					unsigned char current_symbol = l2_table[ table_Idx + deccode ];
					taken = length_table[current_symbol];
				}
				copy_next = next;
				copy_next >>= MAX_BITS - taken;
	
				next <<= taken;
				window <<= taken;
				at += taken;
				bits_count += taken;
				window |= copy_next;
				code_count++;
			}
			input_ptr++;
			tinIdx++;
			window = thread_input[tinIdx];
			next = thread_input[tinIdx+1];

			at -= MAX_BITS;
			copy_next = next;
			window <<= at;
			next <<= at;

			copy_next >>= MAX_BITS - at;
			window |= copy_next;
		}

		// compute prefixsums in warps
		int warpLane = threadIdx.x&(WARP_SIZE-1);
		int warpId = threadIdx.x/WARP_SIZE;
		unsigned int tmp_count = code_count;
		unsigned int tmp_value = 0;
		for( unsigned int delta = 1; delta < WARP_SIZE; delta <<= 1 ){
			tmp_value = __shfl_up_sync(0xFFFFFFFF, tmp_count, delta, WARP_SIZE);
			if( warpLane >= delta ) tmp_count += tmp_value;
		}

		// compute prefixsums in blocks
		__shared__ unsigned int blockprefixsum[NUM_THREADS/WARP_SIZE];
		if( warpLane == WARP_SIZE-1 ){
			blockprefixsum[threadIdx.x/WARP_SIZE] = tmp_count;
		}
		__syncthreads();

		// The first warp compute prefixsums
		if( threadIdx.x < BLOCK_LAST ){
			tmp_value = blockprefixsum[threadIdx.x];
			const unsigned int shfl_mask = ~( (~0) << BLOCK_LAST );
			for( unsigned int delta=1; delta < BLOCK_LAST; delta <<=1 ){
				unsigned int tmp = __shfl_up_sync( shfl_mask, tmp_value, delta, BLOCK_LAST);
				if( threadIdx.x >= delta ) tmp_value += tmp;
			}
			blockprefixsum[threadIdx.x] = tmp_value;
		}
		exclusive_sum = 0;
		__syncthreads();
	
		if( warpId == 0 ){
			int posIdx = block_pos - warpLane;
			unsigned long long int local_sum;
			unsigned long long int before_sum = 0;
			if( warpLane == 0 ){
				local_sum = blockprefixsum[NUM_THREADS/WARP_SIZE-1];
				if( block_pos == 0 ){
					output_pos[block_pos] = local_sum | FLAG_P;
				}else{
					output_pos[block_pos] = local_sum | FLAG_A;
				}
			}

			if( block_pos > 0 ){
				while( 1 ){
					while( posIdx > 0 && (before_sum == 0)){
						before_sum = output_pos[posIdx-1];
					}
					unsigned long long int tmp_sum = 0;
					for(unsigned int delta = 1; delta < WARP_SIZE; delta <<= 1){
						tmp_sum = __shfl_down_sync(0xFFFFFFFF, before_sum, delta, WARP_SIZE);
						if( warpLane < WARP_SIZE - delta && ((before_sum&FLAG_P)==0) ) before_sum += tmp_sum;
					}
					local_sum += (before_sum&(~FLAG_MASK));
					before_sum = __shfl_sync(0xFFFFFFFF, before_sum, 0);
		
					if( before_sum&FLAG_P ){
						break;
					}
					posIdx -= WARP_SIZE;
					before_sum = 0;
				}
				if( warpLane == 0 ){
					output_pos[block_pos] = local_sum | FLAG_P;
					exclusive_sum = local_sum - blockprefixsum[BLOCK_LAST - 1];
				}
			}
		}
		__syncthreads();
		unsigned long long int before_sum = 0;
		int tmp_count_plus = 0;
		if(warpLane == 0){
			before_sum = exclusive_sum;
		}
		if(warpId > 0 && warpLane == 0){
			before_sum = exclusive_sum;
			tmp_count_plus = blockprefixsum[warpId-1];
		}
		before_sum = __shfl_sync(0xFFFFFFFF, before_sum, 0);
		tmp_count += __shfl_sync(0xFFFFFFFF, tmp_count_plus, 0);

		unsigned int output_ptr = (before_sum + tmp_count - code_count)/4;
		unsigned int output_ptr_idx = (before_sum + tmp_count - code_count)&3;

		at = first_at;
		bits_count = at;
		input_ptr = (block_pos*blockDim.x + threadIdx.x)*inputSize;
		tinIdx = 0;

		window = thread_input[0];
		next = thread_input[1];
		copy_next = next;
		copy_next >>= MAX_BITS - at;
		next <<= at;
		window <<= at;
		window |= copy_next;

		unsigned int output_symbols = 0;
		int first_flag=1;
		unsigned char decsymbol;

		while( bits_count < SEGMENT_SIZE && input_ptr < inputfilesize ){
			while( at < MAX_BITS && bits_count < SEGMENT_SIZE ){
				unsigned int deccode = (window&prefix_mask)>>prefix_shift;
				if( deccode < boundary_code ){
					decsymbol = l1_table[deccode];
					taken = length_table[decsymbol];
				}else{
					unsigned int ptrIdx = deccode - boundary_code;
					unsigned int table_Idx = ptr_table[ptrIdx];
					unsigned int table_width = ((table_Idx>>16)&0xFFFF);
					table_Idx &= 0xFFFF;
					unsigned int l2_mask = (~(((unsigned int) (0) - 1) >> table_width))>>prefix_bit;
					unsigned int l2_shift = MAX_BITS - prefix_bit - table_width;

					deccode = (window&l2_mask)>>l2_shift;
					decsymbol = l2_table[ table_Idx + deccode ];
					taken = length_table[decsymbol];
				}
				copy_next = next;
				copy_next >>= MAX_BITS - taken;
	
				next <<= taken;
				window <<= taken;
				at += taken;
				bits_count += taken;
				window |= copy_next;

				UINT_OUT(output_symbols, decsymbol, output_ptr_idx);
				output_ptr_idx++;
				if( (output_ptr_idx&3)== 0 ){
					if( first_flag ){
						first_flag = 0;
						atomicOr( &gpu_output[output_ptr++], output_symbols );
					}else{
						gpu_output[output_ptr++] = output_symbols;
					}
					output_symbols=0;
					output_ptr_idx=0;
				}
			}
			input_ptr++;
			tinIdx++;

			window = thread_input[tinIdx];
			next = thread_input[tinIdx+1];

			at -= MAX_BITS;
			copy_next = next;
			window <<= at;
			next <<= at;

			copy_next >>= MAX_BITS - at;
			window |= copy_next;
		}
		if( output_ptr_idx != 0 ){
			atomicOr( &gpu_output[output_ptr], output_symbols );
		}
	}
}

void decoder_l1_l2(
		unsigned int *input,
		int inputfilesize,
		unsigned int *output,
		int outputfilesize,
		int gap_elements_num,
		void *dectable,
		int tablesize,
		unsigned int prefix_bit,
		unsigned int symbol_count,
		struct TableInfo tableinfo
		){

	unsigned int* gpu_input;
	unsigned int* gpu_output;
	long long unsigned int* gpu_outputpos;
	int* gBlock_counter;
	void *d_stage_table;

	int gaparray_size = (gap_elements_num + GAP_FAC_NUM-1)/GAP_FAC_NUM;
	cudaMalloc( (void **) &gpu_input, sizeof(int)*(gaparray_size+inputfilesize));
	cudaMalloc( (void **) &gpu_output, sizeof(int)*((outputfilesize+3)/4));
	cudaMalloc( (void **) &gpu_outputpos, sizeof(long long int)*(gap_elements_num+1) );
	cudaMalloc( (void **) &gBlock_counter, sizeof(int) );
	cudaMalloc( (void **) &d_stage_table, tablesize );

	// H2D
	float HtDsum=0, decsum=0, DtHsum=0;
	for(int i=0; i<LOOP*2; i++){
		cudaMemset( gBlock_counter, 0, sizeof(int) );
		cudaMemset( gpu_outputpos, 0, sizeof(long long int)*(gap_elements_num+1) );
		cudaMemset( gpu_output, 0, sizeof(int)*(gaparray_size+inputfilesize) );
	
		TIMER_START(H2D)
		cudaMemcpy(gpu_input, input, sizeof(int)*(gaparray_size+inputfilesize), cudaMemcpyHostToDevice);
		cudaMemcpy( d_stage_table, dectable, tablesize, cudaMemcpyHostToDevice );
		TIMER_STOP(H2D)
		CUERROR
	
		int numThreads = NUM_THREADS;
		int numBlocks = (gap_elements_num + numThreads-1)/numThreads;
	
		TIMER_START(dec)
		gpu_dec_l1_l2<<<numBlocks, numThreads, tablesize>>>
		( gpu_input,
		  inputfilesize,
		  d_stage_table,
		  gpu_outputpos,
		  gpu_input,
		  gap_elements_num, 
		  gaparray_size, 
		  gpu_output, 
		  outputfilesize, 
		  gBlock_counter, 
		  tablesize, 
		  prefix_bit,
		  symbol_count,
		  tableinfo);
		TIMER_STOP(dec)
		CUERROR
		// D2H
		TIMER_START(D2H)
		cudaMemcpy(output, gpu_output, sizeof(char)*outputfilesize, cudaMemcpyDeviceToHost);
		TIMER_STOP(D2H)
	
		if(i>=LOOP){
			HtDsum += timeH2D;
			decsum += timedec;
			DtHsum += timeD2H;
		}
	}
	printf("HtoD,%f, dec,%f, DtoH,%f", HtDsum/LOOP, decsum/LOOP, DtHsum/LOOP);

        printf(",SEGMENTSIZE,%d",SEGMENT_SIZE);
        printf(",THREAD_NUM,%d",NUM_THREADS);
        printf(",LOCAL_SEGMENT_NUM,%d",LOCAL_SEGMENT_NUM);
        printf("\n");

	cudaFree( gpu_input );
	cudaFree( gpu_input );
	cudaFree( gpu_outputpos );
	cudaFree( gBlock_counter );
	cudaFree( d_stage_table );
}
