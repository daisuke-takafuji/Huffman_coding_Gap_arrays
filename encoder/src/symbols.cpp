#include "../include/symbols.hpp"
#include <stdio.h>

void get_symbols_num(
		unsigned int *inputfile, 
		unsigned int inputfilesize,
		unsigned int *num_of_symbols)
{
	unsigned int input_idx=0;
	while( input_idx<inputfilesize ){
		num_of_symbols[ GET_CHAR( inputfile[input_idx>>2], input_idx&3 ) ]++;
		input_idx++;
	}
	return;
}


void get_symbols_num(
		unsigned char *inputfile,
		unsigned int inputfilesize,
		unsigned int *num_of_symbols)
{
	for(unsigned int i=0; i<inputfilesize; i++){
		num_of_symbols[inputfile[i]]++;
	}
	return;
}

int store_symbols(
		unsigned int *num_of_symbols,
		struct Symbol *symbols )
{

	int symbol_count = 0;
	for(int i=0; i<MAX_CODE_NUM; i++){
		if( num_of_symbols[i] != 0 ){
			symbols[symbol_count].symbol = i;
			symbols[symbol_count].num = num_of_symbols[i];
			symbol_count++;
		}
	}
	return symbol_count;
}

unsigned long long int get_outputfilesize(
		struct Symbol *symbols,
		unsigned int symbol_count)
{
	unsigned long long int outputfilesize=0;
	for(unsigned int i=0; i<symbol_count; i++){
		struct Symbol symbol = symbols[i];
		outputfilesize += symbol.length*symbol.num;
	}
	return outputfilesize;
}
