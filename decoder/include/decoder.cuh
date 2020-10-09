#ifndef DECODER_
#define DECODER_

void decoder_l1_l2(
		unsigned int *input,
		int inputfilesize,
		unsigned int *output,
		int outputfilesize,
		int gap_element_num,
		void *dectable,
		int tablesize,
		unsigned int prefix_bit,
		unsigned int symbol_count,
		struct TableInfo tableinfo
		);

void cu_make_table_call(
		struct Symbol *symbols,
		struct TableInfo *tableinfo,
		int symbol_count,
		int prefix_bit,
		void *decode_table );

#endif
