#ifndef CONSTANTS_
#define CONSTANTS_

#define MAX_BITS 32
#define MAX_CODEWORD_LENGTH 16
#define MAX_CODE_NUM 256
#define PACK_SIZE 512

#define GAP_LENGTH_MAX 4
#define GAP_ELEMENTS_NUM 8

#define SEGMENTSIZE 128

// For Tesla V100, the number of CUDA blocks per SM is 8,
// the number of CUDA blocks is BLOCKNUM=8*80 
#define BLOCKNUM 8*80
#define THREAD_NUM 256
#define THREAD_ELEMENT 32
#define TNUM_DIV_WSIZE THREAD_NUM/WARP_SIZE
#define WARP_SIZE 32
#define FLAG_A    0x0100000000000000
#define FLAG_P    0x8000000000000000
#define FLAG_MASK 0xFF00000000000000

#define GAP_THREADS 256
#define GAP_BLOCK_RANGE GAP_THREADS*GAP_ELEMENTS_NUM

#define GET_CHAR(value, shift) ((value>>((shift)*8))&0xFF)

#define LEFT 0
#define RIGHT 1
#define ull unsigned long long int

struct Symbol{
	unsigned char symbol;
	unsigned char length;
	unsigned int num;
};

struct Codetable{
	unsigned int code;
	unsigned int length;
};

#endif
