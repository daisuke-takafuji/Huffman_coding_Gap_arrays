#ifndef ENCODER_
#define ENCODER_

void encode(
		unsigned int* outputfile, 
		unsigned int outputfilesize, 
		unsigned int* inputfile, 
		unsigned int inputfilesize, 
		unsigned int gap_array_elements_num,
		struct Codetable* codetable );

extern "C" void histgram(
		unsigned int *d_Histogram,
		void *d_Data,
		unsigned int dataCount);


#endif
