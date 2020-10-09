#ifndef ENCODER_
#define ENCODER_

#define CUERROR  {\
cudaError_t cuError;\
	if( cudaSuccess != (cuError = cudaGetLastError()) ){\
		printf("Error: %s : at %s line %d\n", cudaGetErrorString(cuError), __FILE__, __LINE__);\
		exit(EXIT_FAILURE);\
	} \
}

void encode(
		unsigned int* outputfile, 
		unsigned int outputfilesize, 
		unsigned int* inputfile, 
		unsigned int inputfilesize, 
		unsigned int gap_array_elements_num,
		struct Codetable* codetable );

extern "C" void histgram(
		unsigned int *d_Histogram,
		unsigned int *d_Data,
		unsigned int dataCount);

#define CUERROR  {\
cudaError_t cuError;\
	if( cudaSuccess != (cuError = cudaGetLastError()) ){\
		printf("Error: %s : at %s line %d\n", cudaGetErrorString(cuError), __FILE__, __LINE__);\
		exit(EXIT_FAILURE);\
	} \
}

#endif
