#include <stdio.h>
#include <stdlib.h>

#include "../include/cu_timer.h"
#include "../include/encoder.cuh"
#include "../include/constants.hpp"

#define HIST_THREADS 192
#define WARP_SIZE 32
#define WARP_COUNT (HIST_THREADS/WARP_SIZE)
#define HIST_SIZE 256
#define S_HIST_SIZE (WARP_COUNT*HIST_SIZE)

#define HIST_BLOCK 240

inline __device__ void addByte(
		unsigned int *s_WarpHist,
		unsigned int data){

	atomicAdd(s_WarpHist + data, 1);
}

inline __device__ void addWord(
		unsigned int *s_WarpHist,
		unsigned int data ){

	addByte(s_WarpHist, (data >>  0) & 0xFFU );
	addByte(s_WarpHist, (data >>  8) & 0xFFU );
	addByte(s_WarpHist, (data >> 16) & 0xFFU );
	addByte(s_WarpHist, (data >> 24) & 0xFFU );
}

__global__ void cu_histgram(
		unsigned int *d_PartialHistograms,
		unsigned int *d_Data,
		unsigned int dataCount,
		unsigned int byteCount
		){

	__shared__ unsigned int s_Hist[S_HIST_SIZE];
	unsigned int *s_WarpHist = s_Hist + (threadIdx.x >> 5) * HIST_SIZE;
	unsigned int warpLane = threadIdx.x&31;

	// initialize an array in shared memory
	for(unsigned int i=warpLane; i<HIST_SIZE; i += WARP_SIZE){
		s_WarpHist[i] = 0;
	}
	__syncthreads();

	unsigned int pos=0;
	for( pos = (blockIdx.x*blockDim.x) + threadIdx.x; pos<dataCount-1; pos += (blockDim.x*gridDim.x) ){
		unsigned int data = d_Data[pos];
		addWord(s_WarpHist, data);
	}

	if(pos == dataCount-1){
		unsigned int data = d_Data[pos];
		switch(byteCount&3){
			case 1:
				addByte(s_WarpHist, (data >>  0) & 0xFFU );
				break;
			case 2:
				addByte(s_WarpHist, (data >>  0) & 0xFFU );
				addByte(s_WarpHist, (data >>  8) & 0xFFU );
				break;
			case 3:
				addByte(s_WarpHist, (data >>  0) & 0xFFU );
				addByte(s_WarpHist, (data >>  8) & 0xFFU );
				addByte(s_WarpHist, (data >> 16) & 0xFFU );
				break;
			default:
				addByte(s_WarpHist, (data >>  0) & 0xFFU );
				addByte(s_WarpHist, (data >>  8) & 0xFFU );
				addByte(s_WarpHist, (data >> 16) & 0xFFU );
				addByte(s_WarpHist, (data >> 24) & 0xFFU );
		}
	}

	__syncthreads();

	//
	for(unsigned int bin = threadIdx.x; bin<HIST_SIZE; bin += HIST_THREADS){
		unsigned int sum = 0;
		for(unsigned int i = 0; i < WARP_COUNT; i++){
			sum += s_Hist[bin + i * HIST_SIZE];
		}
		d_PartialHistograms[blockIdx.x * HIST_SIZE + bin] = sum;
	}
}

#define MERGE_THREADBLOCK_SIZE 256
__global__ void mergeHistogram(
		unsigned int *d_Histogram,
		unsigned int *d_PartialHistograms,
		unsigned int histogramCount){

	unsigned int sum = 0;

	for( unsigned int i = threadIdx.x; i < histogramCount; i+= MERGE_THREADBLOCK_SIZE){
		sum += d_PartialHistograms[blockIdx.x + i * HIST_SIZE];
	}

	__shared__ unsigned int data[MERGE_THREADBLOCK_SIZE];
	data[threadIdx.x] = sum;

	for(unsigned int stride = MERGE_THREADBLOCK_SIZE/2; stride > 0; stride >>= 1){
		__syncthreads();
		if( threadIdx.x < stride ){
			data[threadIdx.x] += data[threadIdx.x + stride];
		}
	}

	if( threadIdx.x == 0 ){
		d_Histogram[blockIdx.x] = data[0];
	}
}

extern "C" void histgram(
		unsigned int *d_Histogram,
		unsigned int *d_Data,
		unsigned int byteCount){

	unsigned int dataCount = (byteCount + 3)/4;
	unsigned int *d_PartialHistograms;
	cudaMalloc((void **)&d_PartialHistograms, sizeof(int) * HIST_BLOCK * HIST_SIZE);
	
	cu_histgram<<<HIST_BLOCK, HIST_THREADS>>>(
			d_PartialHistograms,
			d_Data,
			dataCount,
			byteCount);

	mergeHistogram<<<HIST_SIZE, MERGE_THREADBLOCK_SIZE>>>(
			d_Histogram,
			d_PartialHistograms,
			HIST_BLOCK );

	cudaFree(d_PartialHistograms);
	return;
}

__global__ void cuencoder(
		unsigned int* outputfile, 
		unsigned int outputfilesize, 
		unsigned int* inputfile, 
		unsigned int inputfilesize,
		struct Codetable* codetable,
		volatile unsigned long long int *inclusive_sum, 
		unsigned int* gap_array_bytes, 
		unsigned int gap_array_elements_num,
		unsigned int *counter){


	unsigned int* threadInput;
	unsigned int threadInput_idx = 0;
	unsigned int block_idx = 0;
	unsigned int blockNum = inputfilesize/(THREAD_ELEMENT*THREAD_NUM) + (inputfilesize%(THREAD_ELEMENT*THREAD_NUM)!=0);
	__shared__ struct Codetable shared_codetable[MAX_CODE_NUM];
	__shared__ unsigned long long int shared_exclusive_sum;
	__shared__ unsigned int shared_block_idx;

	for( int i=threadIdx.x; i<MAX_CODE_NUM; i+=blockDim.x ){
		shared_codetable[i] = codetable[i];
	}
	if( threadIdx.x == 0 ){
		shared_block_idx = atomicAdd( counter, 1 );
	}
	__syncthreads();
	block_idx = shared_block_idx;
	threadInput_idx = (block_idx*blockDim.x + threadIdx.x)*THREAD_ELEMENT;
	
	while( block_idx < blockNum ){
		unsigned int window = 0;
		unsigned int window_pos = 0;
		unsigned int output_pos = 0;
		unsigned int input_pos = 0;
		unsigned int output_code_bits = 0;
		unsigned int input = 0;
		threadInput = inputfile+(threadInput_idx/4);

		// Get length per thread input
		// /------------------------------------------------------------/
		input = threadInput[0];
		while( input_pos < THREAD_ELEMENT && threadInput_idx + input_pos < inputfilesize ){
			const struct Codetable code = shared_codetable[ GET_CHAR( input, input_pos&3) ];
			output_code_bits += code.length;
			input_pos++;
			if( (input_pos&3) == 0 ) input = threadInput[input_pos/4];
		}
		// /------------------------------------------------------------/

		// // compute prefixsums
		// // prefixsums in warps
		int warpLane = threadIdx.x&(WARP_SIZE-1);
		int warpIdx = threadIdx.x/WARP_SIZE;
		unsigned int tmp_output_code_bits = output_code_bits;
		unsigned int tmp_value = 0;
		for(int delta = 1; delta < WARP_SIZE; delta <<= 1){
			tmp_value = __shfl_up_sync( 0xFFFFFFFF, tmp_output_code_bits, delta, WARP_SIZE );
			if( warpLane >= delta ) tmp_output_code_bits += tmp_value;
		}

		// prefixsums in blocks
		__shared__ unsigned int shared_output_bits[THREAD_NUM/WARP_SIZE];
		if( warpLane == WARP_SIZE-1 ){
			shared_output_bits[warpIdx] = tmp_output_code_bits;
		}
		__syncthreads();

		// Threads in the first warps compute prefixsums 
		if( threadIdx.x < TNUM_DIV_WSIZE ){
			tmp_value = shared_output_bits[threadIdx.x];
			const unsigned int shfl_mask = ~( (~0) << TNUM_DIV_WSIZE );
			for( int delta=1; delta < TNUM_DIV_WSIZE; delta <<= 1 ){
				unsigned int tmp = __shfl_up_sync( shfl_mask, tmp_value, delta, TNUM_DIV_WSIZE );
				if( threadIdx.x >= delta ) tmp_value += tmp;
			}
			shared_output_bits[threadIdx.x] = tmp_value;
		}
		__syncthreads();

		// The first block looks back 32 blocks simultaneously.
		if( warpIdx == 0 ){
			int posIdx = block_idx - warpLane;
			unsigned long long int local_inclusive_sum = 0;
			unsigned long long int exclusive_sum = 0;
			if( warpLane == 0 ){
				local_inclusive_sum = shared_output_bits[TNUM_DIV_WSIZE-1];
				if( block_idx == 0 ){
					inclusive_sum[block_idx] = local_inclusive_sum | FLAG_P;
					shared_exclusive_sum = 0;
				}
				else{
					inclusive_sum[block_idx] = local_inclusive_sum | FLAG_A;
				}
			}

			if( block_idx > 0 ){
				while( 1 ){
					while( posIdx > 0 && (exclusive_sum == 0)){
						exclusive_sum = inclusive_sum[posIdx-1];
					}
					unsigned long long int tmp_sum = 0;
					for(unsigned int delta = 1; delta < WARP_SIZE; delta <<= 1){
						tmp_sum = __shfl_down_sync(0xFFFFFFFF, exclusive_sum, delta, WARP_SIZE);
						if( warpLane < (WARP_SIZE - delta) && ((exclusive_sum&FLAG_P)==0) ) exclusive_sum += tmp_sum;
					}
					local_inclusive_sum += (exclusive_sum&(~FLAG_MASK));
					exclusive_sum = __shfl_sync(0xFFFFFFFF, exclusive_sum, 0);
		
					if( exclusive_sum&FLAG_P ){
						break;
					}
		
					posIdx -= WARP_SIZE;
					exclusive_sum = 0;
				}
				if( warpLane == 0 ){
					inclusive_sum[block_idx] = ((local_inclusive_sum&(~FLAG_MASK)) | FLAG_P);
					shared_exclusive_sum = local_inclusive_sum - shared_output_bits[TNUM_DIV_WSIZE-1];
				}
			}
		}
		__syncthreads();

		unsigned long long int exclusive_sum = 0;
		unsigned int tmp_count_plus = 0;
		if(warpLane == 0){
			exclusive_sum = shared_exclusive_sum;
		}
		if( warpIdx > 0 && warpLane == 0 ){
			tmp_count_plus = shared_output_bits[warpIdx-1];
		}
		exclusive_sum = __shfl_sync( 0xFFFFFFFF, exclusive_sum, 0 );
		tmp_output_code_bits += __shfl_sync(0xFFFFFFFF, tmp_count_plus, 0);

		exclusive_sum = exclusive_sum + tmp_output_code_bits - output_code_bits;
		// /------------------------------------------------------------/
		// Output encoded data
		// /------------------------------------------------------------/
		unsigned int output_bits = (exclusive_sum&(SEGMENTSIZE-1));
		window_pos = (exclusive_sum&(MAX_BITS-1));
		output_pos = (exclusive_sum/MAX_BITS);

		window = 0;
		input_pos = 0;
		int first_flag = 1;
		int last_out_flag = 0;
		input = threadInput[0];
		while( input_pos < THREAD_ELEMENT && threadInput_idx + input_pos < inputfilesize ){
			struct Codetable code = shared_codetable[ GET_CHAR( input, input_pos&3 ) ];
			input_pos++;
			if( (input_pos&3)==0 ) input = threadInput[input_pos/4];
			while( window_pos + code.length < MAX_BITS && threadInput_idx + input_pos < inputfilesize && input_pos < THREAD_ELEMENT ){
				window <<= code.length;
				window += code.code;
				window_pos += code.length;

				output_bits += code.length;
				if( threadInput_idx + input_pos < inputfilesize && input_pos < THREAD_ELEMENT ){
					code = shared_codetable[ GET_CHAR( input, input_pos&3 ) ];
					input_pos++;
					if( (input_pos&3) == 0 ) input = threadInput[input_pos/4];
				}
			}

			output_bits += code.length;
			if( output_bits/SEGMENTSIZE != (output_bits-code.length)/SEGMENTSIZE ){
				const int gap_pos = output_pos/(SEGMENTSIZE/MAX_BITS);
				unsigned int gap_elements = output_bits&(MAX_CODEWORD_LENGTH-1);
				gap_array_bytes[gap_pos] = gap_elements;
			}

			const int diff = window_pos + code.length - MAX_BITS;
			last_out_flag = diff;

			if( diff >= 0 ){
				window <<= code.length - diff;
				window += (code.code >> diff);

				if( first_flag ){
					atomicOr( &outputfile[output_pos++], window );
					first_flag = 0;
				}
				else{
					outputfile[output_pos++] = window;
				}

				window = code.code & ~(~0 << diff);
				window_pos = diff;
			}
			else{
				window <<= code.length;
				window |= code.code;

				const int shift = MAX_BITS - (window_pos+code.length);
				window <<= shift;
				atomicOr( &outputfile[output_pos++], window );
				window_pos = 0;
				last_out_flag = 0;
			}
		}
		// Output remained bits
		if( last_out_flag != 0 ){
			window <<= (MAX_BITS - last_out_flag);
			atomicOr( &outputfile[output_pos++], window );
		}
		// assign segments to blocks
		if( threadIdx.x == 0 ) shared_block_idx = atomicAdd( counter, 1 );
		__syncthreads();
		block_idx = shared_block_idx;
		threadInput_idx = (block_idx*blockDim.x + threadIdx.x)*THREAD_ELEMENT;
		// /------------------------------------------------------------/
	}
}

// Kernel for gathering elements for gap array
__global__ void cu_get_gaparray( unsigned int *gap_array_bytes, unsigned int *gap_array, unsigned int gap_array_elements_num, unsigned int gap_array_size ){
	__shared__ unsigned int shared_gap_array[GAP_BLOCK_RANGE];

	int block_start_pos = blockIdx.x*GAP_BLOCK_RANGE;
	int thread_start_pos = threadIdx.x*GAP_ELEMENTS_NUM;
	int gap_array_idx = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int gap_element=0;

	for(int i=threadIdx.x; i<GAP_BLOCK_RANGE && block_start_pos+i < gap_array_elements_num; i+=blockDim.x){
		shared_gap_array[i] = gap_array_bytes[ block_start_pos + i ];
	}
	__syncthreads();

	for(int i=0; i<GAP_ELEMENTS_NUM && block_start_pos+thread_start_pos+i < gap_array_elements_num; i++){
		const int gap_shift = GAP_LENGTH_MAX*i;
		gap_element |= (shared_gap_array[ thread_start_pos + i ]<<gap_shift);

	}
	if( gap_array_idx < gap_array_size )
		gap_array[gap_array_idx] = gap_element;

}

// launch Kernels for encode
void encode(
		unsigned int* outputfile, 
		unsigned int outputfilesize, 
		unsigned int* d_inputfile, 
		unsigned int inputfilesize, 
		unsigned int gap_array_elements_num,
		struct Codetable* codetable ){

	unsigned int
				 *d_outputfile,
				 *d_gap_array_bytes,
				 *d_gap_array,
				 *d_counter;
	unsigned long long int *d_inclusive_sum;
	struct Codetable *d_codetable;

	unsigned int *gap_array;
	gap_array = outputfile + outputfilesize;
	unsigned int gap_array_size=0;
	gap_array_size = (gap_array_elements_num + GAP_ELEMENTS_NUM-1)/GAP_ELEMENTS_NUM;

	unsigned int blockNum = (inputfilesize + THREAD_ELEMENT*THREAD_NUM-1)/(THREAD_ELEMENT*THREAD_NUM);
	cudaMalloc( &d_outputfile, sizeof(int)*outputfilesize );
	cudaMalloc( &d_codetable, sizeof(struct Codetable)*MAX_CODE_NUM );
	cudaMalloc( &d_inclusive_sum, sizeof(unsigned long long int)*blockNum );
	cudaMalloc( &d_gap_array_bytes, sizeof(unsigned int)*gap_array_elements_num );
	cudaMalloc( &d_gap_array, sizeof(unsigned int)*gap_array_size );
	cudaMalloc( &d_counter, sizeof(int) );
	CUERROR;

	cudaMemset( d_outputfile, 0, sizeof(int)*outputfilesize );
	cudaMemset( d_inclusive_sum, 0, sizeof(unsigned long long int)*blockNum );
	cudaMemset( d_gap_array_bytes, 0, sizeof(int)*gap_array_elements_num );
	cudaMemset( d_counter, 0, sizeof(int) );

	TIMER_START(HtD)
	cudaMemcpy( d_codetable, codetable, sizeof(struct Codetable)*MAX_CODE_NUM, cudaMemcpyHostToDevice );
	TIMER_STOP(HtD)
	CUERROR;

	TIMER_START(cuenc)
	cuencoder<<<BLOCKNUM, THREAD_NUM>>>
		( d_outputfile,
		  outputfilesize, 
		  d_inputfile,
		  inputfilesize,
		  d_codetable,
		  d_inclusive_sum,
		  d_gap_array_bytes,
		  gap_array_elements_num,
		  d_counter);
	cudaDeviceSynchronize();
	TIMER_STOP(cuenc)
	CUERROR;


	blockNum = gap_array_size;
	blockNum = (blockNum+GAP_THREADS-1)/GAP_THREADS;

	TIMER_START(gapget)
	cu_get_gaparray<<<blockNum, GAP_THREADS>>>(d_gap_array_bytes, d_gap_array, gap_array_elements_num, gap_array_size);
	TIMER_STOP(gapget)
	CUERROR;

	TIMER_START(DtH)
	cudaMemcpy(outputfile, d_outputfile, sizeof(int)*outputfilesize, cudaMemcpyDeviceToHost);
	cudaMemcpy(gap_array, d_gap_array, sizeof(int)*gap_array_size, cudaMemcpyDeviceToHost);
	TIMER_STOP(DtH)

	cudaFree(d_inputfile);
	cudaFree(d_outputfile);
	cudaFree(d_codetable);
	cudaFree(d_inclusive_sum);
	cudaFree(d_gap_array_bytes);
	cudaFree(d_gap_array);
}
