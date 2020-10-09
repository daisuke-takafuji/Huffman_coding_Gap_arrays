#ifndef PACKAGE_MERGE_
#define PACKAGE_MERGE_

#include "constants.hpp"

#define ull unsigned long long int

struct NodePack{
	ull left_weight;
	ull right_weight;
	ull weight;
	int counter;
	int LR;
	int pack_pos;
	struct NodePack *chain;
};

struct NodePackList{
	struct NodePack *list;
	int current_pos[MAX_CODEWORD_LENGTH];
};

void boundary_PM(
		struct Symbol *symbols,
		unsigned int symbol_count,
		struct Codetable *codetable);


void add_node(
		struct Symbol *symbols, 
		int symbol_count, 
		struct NodePackList &pack, 
		int listpos,
		int symbol_length);


void read_chain(
		struct Symbol *symbols,
		struct NodePack *list);
#endif
