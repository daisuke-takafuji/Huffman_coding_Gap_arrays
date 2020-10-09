#ifndef SYMBOLS_
#define SYMBOLS_
#include "constants.hpp"

void get_symbols_num(
		unsigned char *inputfile,
		unsigned int inputfilesize,
		unsigned int *num_of_symbols);

void get_symbols_num(
		unsigned int *inputfile,
		unsigned int inputfilesize,
		unsigned int *num_of_symbols);

int store_symbols(
		unsigned int *num_of_symbols,
		struct Symbol *symbols);


unsigned long long int get_outputfilesize(
		struct Symbol *symbols,
		unsigned int symbol_count);

#endif
