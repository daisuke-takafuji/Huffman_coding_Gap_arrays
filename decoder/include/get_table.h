#ifndef GET_TABLE_
#define GET_TABLE_
#include "constants.hpp"

unsigned int get_table_info(
		struct Symbol* symbols,
		int symbol_count,
		int prefix_bit,
		struct TableInfo &tableinfo);

unsigned int get_twolevel_table(
		void *stage_table,
		int prefix_bit,
		struct Symbol* symbols, 
		int symbol_count,
		struct TableInfo tableinfo);

int get_suffixnum(
		struct Dectable* dectable,
		struct Symbol* symbols, 
		int symbol_count );

void get_prefix_suffix_table(
		struct Dectable* prefixtable,
		struct Dectable *suffixtable, 
		struct Symbol* symbols,
		int symbol_count,
		int prefix_bit);

#endif
