// Usage: bin/encoder input output
#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <unistd.h>
#include "../include/huffman.h"

void fatal(const char* str){
	printf("%s", str);
	fflush(stdout);
	exit(EXIT_FAILURE);
}

int compare( const void *p, const void *q ){
	struct Symbol *P = (struct Symbol *) p;
	struct Symbol *Q = (struct Symbol *) q;
	return P->num - Q->num;
}

int main(int argc, char **argv){
	if( argc != 3 ){
		printf("Please enter : bin/encoder input output\n");
		exit(EXIT_FAILURE);
	}

	FILE *input, *output;
	unsigned int *inputfile;
	unsigned int *outputfile;
	unsigned int *gap_array;
	unsigned int inputfilesize=0;
	unsigned int outputfilesize=0;
	unsigned long long int outputfilesize_bits=0;
	unsigned int gap_array_size = 0;
	unsigned int gap_array_elements_num = 0;
	unsigned int *num_of_symbols;
	num_of_symbols = (unsigned int *)malloc(sizeof(int)*MAX_CODE_NUM);
	struct Symbol symbols[MAX_CODE_NUM] = {};
	struct Codetable *codetable;

	input = fopen(argv[1], "rb");
	output = fopen(argv[2], "wb");

	fseek(input, 0L, SEEK_END);
	inputfilesize = ftell(input);
	fseek(input, 0L, SEEK_SET);

	cudaMallocHost( &codetable, sizeof(struct Codetable)*MAX_CODE_NUM );
	cudaMallocHost( &inputfile, sizeof(int)*((inputfilesize+3)/4) );

	fread(inputfile, sizeof(char), inputfilesize, input);
	fsync(input->_fileno);

	int symbol_count = 0;

	unsigned int *d_num_of_symbols;
	unsigned int *d_inputfile;
	cudaMalloc((void **)&d_inputfile, sizeof(int)*((inputfilesize+3)/4 + 1 ));
	cudaMalloc((void **)&d_num_of_symbols, sizeof(int)*MAX_CODE_NUM);

	TIMER_START(HtD1)
	cudaMemcpy(d_inputfile, inputfile, sizeof(int)*((inputfilesize+3)/4), cudaMemcpyHostToDevice);
	CUERROR
	TIMER_STOP(HtD1)

	// compute the histogram
	TIMER_START(hist)
	histgram(d_num_of_symbols, d_inputfile, inputfilesize);
	cudaMemcpy(num_of_symbols, d_num_of_symbols, sizeof(int)*MAX_CODE_NUM, cudaMemcpyDeviceToHost);
	TIMER_STOP(hist)

	cudaFree(d_num_of_symbols);

	// make symbols
	CPU_TIMER_START(dic)
	symbol_count = store_symbols(num_of_symbols, symbols);
	qsort(symbols, symbol_count, sizeof(struct Symbol) ,compare);

	boundary_PM(symbols, symbol_count, codetable);
	CPU_TIMER_STOP(dic)

	outputfilesize_bits = get_outputfilesize(symbols, symbol_count);
	gap_array_elements_num = (outputfilesize_bits + SEGMENTSIZE-1)/SEGMENTSIZE;
	gap_array_size = gap_array_elements_num/GAP_ELEMENTS_NUM + ((gap_array_elements_num%GAP_ELEMENTS_NUM) != 0);

	outputfilesize = (outputfilesize_bits+MAX_BITS-1)/MAX_BITS;

	cudaMallocHost( &outputfile, sizeof(unsigned int)*(outputfilesize+gap_array_size) );
	gap_array = outputfile + outputfilesize;

	encode(
		outputfile, 
		outputfilesize, 
		d_inputfile, 
		inputfilesize, 
		gap_array_elements_num,
		codetable );

	CUERROR;

        printf("file,%s",argv[1]);
        printf(",SEGMENTSIZE,%d",SEGMENTSIZE);
        printf(",THREAD_NUM,%d",THREAD_NUM);
        printf(",Bytes_per_THREAD,%d",THREAD_ELEMENT);
        printf("\n");

	//------------------------------------------------------------ 

	size_t tmp_symbol_count = symbol_count;
	CPU_TIMER_START(write)
	fwrite( &tmp_symbol_count, sizeof(size_t), 1, output );
	for(int i=symbol_count-1; i>=0; i--){
		unsigned char tmpsymbol = symbols[i].symbol;
		unsigned char tmplength = symbols[i].length;
		fwrite( &tmpsymbol, sizeof(tmpsymbol), 1, output );
		fwrite( &tmplength, sizeof(tmplength), 1, output );
	}

	fwrite( &inputfilesize, sizeof(inputfilesize), 1, output );
	fwrite( &outputfilesize, sizeof(outputfilesize), 1, output );
	fwrite( &gap_array_elements_num, sizeof(gap_array_elements_num), 1, output );
	fwrite( gap_array, sizeof(int), gap_array_size, output );
	fwrite( outputfile, sizeof(unsigned int), outputfilesize, output);

	fdatasync(output->_fileno);
	CPU_TIMER_STOP(write)

		printf("running time (ms) (including CPU time), %lf\n", sum_of_time);
	int outsize = ftell(output);
	printf("ratio=outputfile/inputfile,%lf", (double)outsize/inputfilesize);
	printf(",outputfilesize,%d,inputfilesize,%d\n", outsize, inputfilesize);
	printf("\n");

	fclose(input);
	fclose(output);

	cudaFreeHost(inputfile);
	cudaFreeHost(outputfile);
	cudaFreeHost(codetable);

	return 0;
}

