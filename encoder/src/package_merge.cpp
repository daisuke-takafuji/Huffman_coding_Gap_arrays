#include <stdio.h>
#include <stdlib.h>

#include "../include/package_merge.hpp"
#include "../include/constants.hpp"

#define LEFT 0
#define RIGHT 1
#define ull unsigned long long int
#define LIST_SIZE 2

void add_node(
		struct Symbol *symbols,
		int symbol_count, 
		struct NodePackList &pack, 
		int listpos,
		int symbols_length){

	struct NodePack *next_pack;
	struct NodePack *current_pack;
	struct NodePack *before_pack;
	// struct NodePack *save_pack;
	int maxlistLength = 2*(symbol_count-1);
	int above_pack_pos = 0;
	if( listpos > 0 ){
		above_pack_pos = pack.current_pos[listpos-1];
	}
	current_pack = &pack.list[ maxlistLength*listpos + pack.current_pos[listpos] ];
	next_pack = &pack.list[ maxlistLength*listpos + pack.current_pos[listpos]+1 ];
	int symbolsIdx = current_pack->counter;
	next_pack->counter = current_pack->counter;
	if(listpos > 0){
		// left
		{
			above_pack_pos = pack.current_pos[listpos-1];
			before_pack = &pack.list[ maxlistLength*(listpos-1) + above_pack_pos];
			if( before_pack->weight <= symbols[symbolsIdx].num || symbolsIdx >= symbols_length){
				next_pack->left_weight = before_pack->weight;
				next_pack->chain = before_pack;
				add_node( symbols, symbol_count, pack, listpos-1 , symbols_length);
			}
			else{
				next_pack->left_weight = symbols[symbolsIdx++].num;
				next_pack->counter++;
				next_pack->chain = current_pack->chain;
			}
		}
		// right
		{
			above_pack_pos = pack.current_pos[listpos-1];
			before_pack = &pack.list[ maxlistLength*(listpos-1) + above_pack_pos];
			if( before_pack->weight <= symbols[symbolsIdx].num || symbolsIdx >= symbols_length ){
				next_pack->right_weight = before_pack->weight;
				next_pack->chain = before_pack;
				add_node( symbols, symbol_count, pack, listpos-1, symbols_length );
			}
			else{
				next_pack->right_weight = symbols[symbolsIdx++].num;
				next_pack->counter++;
			}
		}
	}
	else{
		if(symbolsIdx < symbols_length){
			next_pack->left_weight = symbols[symbolsIdx++].num;
			next_pack->counter++;
		}
		if(symbolsIdx < symbols_length){
			next_pack->right_weight = symbols[symbolsIdx++].num;
			next_pack->counter++;
		}
		else{
			next_pack->right_weight = 0;
		}
		next_pack->chain = NULL;
	}

	pack.current_pos[listpos] += 1;
	next_pack->weight = next_pack->left_weight + next_pack->right_weight;
}

void read_chain(
		struct Symbol *symbols,
		struct NodePack *list){

	struct NodePack *chain_node;
	int listpos = MAX_CODEWORD_LENGTH-2;
	struct NodePack *before_pack, *before_pack_save, *save_pack;
	before_pack = &list[(listpos-1)*LIST_SIZE + 1];
	before_pack_save = before_pack - 1;


	while( (chain_node = before_pack->chain) != NULL ){
		save_pack = &list[MAX_CODEWORD_LENGTH*LIST_SIZE + listpos];
		*before_pack_save = *before_pack;
		*save_pack = *before_pack;
		before_pack = chain_node;
		before_pack_save = before_pack - 1;
		listpos--;
	}
	if( chain_node == NULL ){
		printf("READ listpos %d\n\n", listpos);
	}
}

// Boudary package merge algorithm
void boundary_PM(
		struct Symbol *symbols,
		unsigned int symbol_count,
		struct Codetable *codetable){

	struct NodePackList pack;
	int maxlistLength = 2*(symbol_count-1);
	pack.list = (struct NodePack *)malloc( sizeof(struct NodePack)*maxlistLength*MAX_CODEWORD_LENGTH ) ;
	for(int i=0; i<MAX_CODEWORD_LENGTH; i++){
		pack.current_pos[i] = 0;
	}

	int num_of_symbols = symbol_count;
	int looptime = (num_of_symbols-1)*2;

	// initialize NodePack
	for( int i=0; i<MAX_CODEWORD_LENGTH; i++ ){
		struct NodePack init_node;
		init_node.left_weight = symbols[0].num;
		init_node.right_weight = symbols[1].num;
		init_node.weight = init_node.left_weight + init_node.right_weight;
		init_node.counter = 2;
		init_node.LR = LEFT;
		init_node.chain = NULL;
		init_node.pack_pos = 0;
		pack.list[i*maxlistLength] = init_node;
	}

	int listpos = MAX_CODEWORD_LENGTH-1;
	ull last_counter = 2;
	struct NodePack *before_pack;
	struct NodePack *chain_head;
	chain_head = NULL;
	before_pack = &pack.list[ maxlistLength*(listpos-1)];
	for( int i=2; i<looptime; i++ ){
		if( last_counter < symbol_count ){
			if( before_pack->weight <= symbols[last_counter].num ){
				add_node(symbols, symbol_count, pack, listpos-1, symbol_count);
				chain_head = before_pack;
				before_pack++;
			}
			else{
				last_counter++;
			}
		}
		else{
			add_node(symbols, symbol_count, pack, listpos-1, symbol_count);
			chain_head = before_pack;
			before_pack++;
		}
	}

	for(unsigned int i=0; i<last_counter; i++) symbols[i].length++;
	struct NodePack *chain;
	chain = chain_head;
	while( chain != NULL ){
		int counter = chain->counter;
		chain = chain->chain;
		for( int j=0; j<counter; j++ ) symbols[j].length++;
	}

	unsigned int code=0;
	unsigned int current_length=0;
	unsigned int next_length=0;

	current_length = symbols[symbol_count-1].length;
	for(int i=symbol_count-1; i>=0; i--){
		codetable[symbols[i].symbol].code = code;
		codetable[symbols[i].symbol].length = current_length;

		next_length = (i==0) ? current_length : symbols[i-1].length;

		code = (code + 1) << (next_length - current_length);
		current_length = next_length;
	}
}
