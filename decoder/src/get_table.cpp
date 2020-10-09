#include "../include/get_table.h"

unsigned int get_table_info(
		struct Symbol* symbols,
		int symbol_count,
		int prefix_bit,
		struct TableInfo &tableinfo){

	unsigned int code=0;
	unsigned int previous_prefix_code=0;
	unsigned int prefix_code = 0;
	unsigned int ptrtable_size=0;
	unsigned int l2table_size=0;
	unsigned int l1table_size=0;
	unsigned int prefix_mask = ~ ( (~(0))<<prefix_bit );

	for(int i=0; i<symbol_count; i++){
		int current_length = symbols[i].length;
		int next_length = (i+1 == symbol_count) ? current_length : symbols[i+1].length;
		code = (code+1) << (next_length - current_length);

		if( current_length > prefix_bit ){
			unsigned int prefix_shift = current_length - prefix_bit;
			prefix_code = ((code>>prefix_shift)&prefix_mask);
			if( l1table_size == 0 ){
				l1table_size = (code>>prefix_shift);
			}

			if( previous_prefix_code != prefix_code ){
				l2table_size += 1<<prefix_shift;
				ptrtable_size++;
			}
			previous_prefix_code = prefix_code;
		}
	}

	tableinfo.l1table_size = l1table_size;
	tableinfo.l2table_size = l2table_size;
	tableinfo.ptrtable_size = ptrtable_size;
	if( l1table_size == 0 ){
		int bit = symbols[symbol_count-1].length;
		tableinfo.l1table_size = (1<<bit);
		return bit;
	}
	return prefix_bit;
}

unsigned int get_twolevel_table(
		void *stage_table,
		int prefix_bit,
		struct Symbol* symbols, 
		int symbol_count,
		struct TableInfo tableinfo){
	
	unsigned int code = 0;
	unsigned int boundary_code = 0;
	unsigned char *l1_table;
	unsigned char *l2_table;
	unsigned char *length_table;
	unsigned int *ptr_table;

	unsigned int ptrIdx = 0;
	unsigned char ptr_counter=0;

	ptr_table = (unsigned int *)stage_table;
	length_table = (unsigned char *)stage_table + tableinfo.ptrtable_size*sizeof(unsigned int);
	l1_table = length_table + MAX_CODE_NUM;
	l2_table = l1_table + tableinfo.l1table_size;

	for( int i=0; i<symbol_count; ){
		int current_length = symbols[i].length;
		unsigned char current_symbol = symbols[i].symbol;
		length_table[current_symbol] = current_length;

		if( current_length <= prefix_bit ){
			const int prefix_shift = prefix_bit - current_length;
			const int num_elements = 1<<prefix_shift;
			for( int elementsIdx = 0; elementsIdx < num_elements; elementsIdx++ ){
				l1_table[ (code<<prefix_shift) + elementsIdx ] = current_symbol;
			}
			int next_length = (i==symbol_count-1) ? current_length : symbols[i+1].length;
			code = ((code+1) << (next_length - current_length));
			i++;
		}else{
			int prefix_shift = current_length - prefix_bit;
			int prefix_shift_tmp = prefix_shift;
			unsigned int prefix_mask = ~( (~0)<<prefix_bit );

			unsigned int code_tmp = code;
			unsigned int last_length=current_length - prefix_bit;
			int counter = 0;

			if( boundary_code == 0 ){
				boundary_code = (code>>prefix_shift);
			}

			ptr_counter++;

			code_tmp = code;
			while( ((code_tmp>>prefix_shift_tmp)&prefix_mask) == ((code>>prefix_shift)&prefix_mask) ){
				if( i+counter==symbol_count ) break;
				counter++;
				current_length = symbols[i + counter-1].length;
				int next_length = (i+counter==symbol_count-1) ? current_length : symbols[i+counter].length;
				prefix_shift_tmp = next_length - prefix_bit;
				last_length = current_length - prefix_bit;
				code_tmp = ((code_tmp+1) << (next_length - current_length));
			}
			code_tmp = code;

			const unsigned int table_width = last_length;
			unsigned char* this_table = l2_table + ptrIdx;
			for(int symbolIdx=0; symbolIdx<counter; symbolIdx++){
				current_length = symbols[ i + symbolIdx ].length;
				current_symbol = symbols[ i + symbolIdx ].symbol;
				length_table[current_symbol] = current_length;
				const unsigned int suffix_shift = current_length - prefix_bit;
				const unsigned int suffix_mask = ~ ( (~(0)) << suffix_shift );
				const unsigned int suffix_code = (code_tmp&suffix_mask)<<(table_width-suffix_shift);
				const int elementsNum = 1<<(table_width-suffix_shift);


				for(int elementsIdx=0; elementsIdx<elementsNum; elementsIdx++){
					this_table[ suffix_code + elementsIdx ] = current_symbol;
				}
				int next_length = (i+symbolIdx==symbol_count-1) ? current_length : symbols[i+symbolIdx+1].length;
				code_tmp = ((code_tmp+1) << (next_length - current_length));
			}
			
			unsigned int ptr_value = (table_width<<16) | ptrIdx;
			ptr_table[ ptr_counter-1 ] = ptr_value;

			i += counter;
			ptrIdx += (1<<table_width);
			code = code_tmp;
		}
	}
	return boundary_code;
}


int get_suffixnum(
		struct Dectable* dectable,
		struct Symbol* symbols,
		int symbol_count ){

	unsigned int code = 0;
	unsigned int current_length = 0;
	unsigned int next_length = 0;

	enum PREFIX_LENGTH {
		P_9,
		P_10,
		P_11,
		P_12
	};
	// prefix  9-bit 0
	// prefix 10-bit 1
	// prefix 11-bit 2
	unsigned int previous_prefix[4] = {};
	int sum[4]={};

	current_length = symbols[0].length;
	for( int i = 0; i < symbol_count; i++ ){

		if( current_length > SEPARATE_MIN ){
			int pidx = current_length - SEPARATE_MIN - 1;
			unsigned int prefix_idx[4] = {};
			prefix_idx[P_9] = current_length - SEPARATE_MIN - P_9;
			prefix_idx[P_10] = current_length - SEPARATE_MIN - P_10;
			prefix_idx[P_11] = current_length - SEPARATE_MIN - P_11;
			prefix_idx[P_12] = current_length - SEPARATE_MIN - P_12;
	
			switch( pidx ){
				case P_9:
					if( previous_prefix[P_9] != (code>>prefix_idx[P_9]) ){
						sum[P_9]++;
						previous_prefix[P_9] = (code>>prefix_idx[P_9]);
					}
					break;
				case P_10:
					if( previous_prefix[P_9] != (code>>prefix_idx[P_9]) ){
						sum[P_9]++;
						previous_prefix[P_9] = (code>>prefix_idx[P_9]);
					}
					if( previous_prefix[P_10] != (code>>prefix_idx[P_10]) ){
						sum[P_10]++;
						previous_prefix[P_10] = (code>>prefix_idx[P_10]);
					}
					break;
				case P_11:
					if( previous_prefix[P_9] != (code>>prefix_idx[P_9]) ){
						sum[P_9]++;
						previous_prefix[P_9] = (code>>prefix_idx[P_9]);
					}
					if( previous_prefix[P_10] != (code>>prefix_idx[P_10]) ){
						sum[P_10]++;
						previous_prefix[P_10] = (code>>prefix_idx[P_10]);
					}
					if( previous_prefix[P_11] != (code>>prefix_idx[P_11]) ){
						sum[P_11]++;
						previous_prefix[P_11] = (code>>prefix_idx[P_11]);
					}
					break;
				default:
					if( previous_prefix[P_9] != (code>>prefix_idx[P_9]) ){
						sum[P_9]++;
						previous_prefix[P_9] = (code>>prefix_idx[P_9]);
					}
					if( previous_prefix[P_10] != (code>>prefix_idx[P_10]) ){
						sum[P_10]++;
						previous_prefix[P_10] = (code>>prefix_idx[P_10]);
					}
					if( previous_prefix[P_11] != (code>>prefix_idx[P_11]) ){
						sum[P_11]++;
						previous_prefix[P_11] = (code>>prefix_idx[P_11]);
					}
					if( previous_prefix[P_12] != (code>>prefix_idx[P_12]) ){
						sum[P_12]++;
						previous_prefix[P_12] = (code>>prefix_idx[P_12]);
					}
					break;
			}
		}
		next_length = (i==symbol_count-1) ? current_length : symbols[i+1].length;
		code = (code+1) << ( next_length - current_length );
		current_length = next_length;
	}

	int tablesize = (1<<(P_9 + SEPARATE_MIN)) + sum[P_9]*(1<<(MAX_CODEWORD_LENGTH-SEPARATE_MIN));
	int tablesize_tmp = (1<<(P_9 + SEPARATE_MIN)) + sum[P_9]*(1<<(MAX_CODEWORD_LENGTH-SEPARATE_MIN));
	int min_value=0;
	for( int i=0; i<4; i++ ){
		tablesize_tmp = (1<<(P_9 + SEPARATE_MIN + i)) + sum[P_9+i]*(1<<(MAX_CODEWORD_LENGTH - SEPARATE_MIN - i));
		if( tablesize_tmp < tablesize ){
			min_value = i;
			tablesize = tablesize_tmp;
		} 
	}
	return sum[min_value] + ((min_value+SEPARATE_MIN)<<16);
}


void get_dectable( struct Dectable* dectable, struct Symbol* symbols, int symbol_count ){
	unsigned int code = 0;
	unsigned int current_length = 0;
	unsigned int next_length = 0;

	current_length = symbols[0].length;
	for( int i = 0; i < symbol_count; i++ ){
		const unsigned char current_symbol = symbols[i].symbol;

		const int shift = MAX_CODEWORD_LENGTH - current_length;
		const int num_elements = 1<<shift;

		for( int elementsIdx = 0; elementsIdx < num_elements; elementsIdx++ ){
			dectable[ (code<<shift) + elementsIdx ].symbol = current_symbol;
			dectable[ (code<<shift) + elementsIdx ].length = current_length;
		}

		next_length = (i==symbol_count-1) ? current_length : symbols[i+1].length;
		code = (code+1) << ( next_length - current_length );
		current_length = next_length;
	}
}

void get_prefix_suffix_table(
		struct Dectable* prefixtable,
		struct Dectable *suffixtable, 
		struct Symbol* symbols,
		int symbol_count,
		int prefix_bit){

	unsigned int code = 0;
	unsigned int prefix_code = 0;
	unsigned int suffix_code = 0;
	unsigned int previous_prefix_code = 0;
	int current_length = 0;
	unsigned int next_length = 0;
	int suffix_pos = -1;
	unsigned int suffix_bit = MAX_CODEWORD_LENGTH-prefix_bit;

	for( int i=0; i<symbol_count; i++ ){
		current_length = symbols[i].length;
		const unsigned char current_symbol = symbols[i].symbol;

		if( current_length <= prefix_bit ){	// only prefix_code
			const int prefix_shift = prefix_bit - current_length;
			const int num_elements = 1<<prefix_shift;
			for( int elementsIdx = 0; elementsIdx < num_elements; elementsIdx++ ){
				prefixtable[ (code<<prefix_shift) + elementsIdx ].symbol = current_symbol;
				prefixtable[ (code<<prefix_shift) + elementsIdx ].length = current_length;
			}
		}
		else{	// with prefix_code and suffix_code
			const int prefix_shift = current_length - prefix_bit;
			const int suffix_shift = suffix_bit - prefix_shift;
			const int suffix_mask = ~( (~(0)) << suffix_bit );

			prefix_code = code>>prefix_shift;
			suffix_code = ((code<<suffix_shift)&suffix_mask);

			if( prefix_code != previous_prefix_code ) suffix_pos++;
			previous_prefix_code = prefix_code;
			prefixtable[ prefix_code ].length = 0xFF;
			prefixtable[ prefix_code ].symbol = suffix_pos;

			const int num_elements = 1<<suffix_shift;
			int suffixIdx = suffix_pos*(1<<suffix_bit);
			for( int elementsIdx=0; elementsIdx < num_elements; elementsIdx++ ){
				suffixtable[ suffixIdx + (suffix_code) + elementsIdx ].symbol = current_symbol;
				suffixtable[ suffixIdx + (suffix_code) + elementsIdx ].length = current_length;
			}
		}
		next_length = (i==symbol_count-1) ? current_length : symbols[i+1].length;
		code = (code+1) << (next_length - current_length);
		current_length = next_length;
	}
}
