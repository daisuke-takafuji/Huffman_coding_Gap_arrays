#ifndef CONSTANTS_
#define CONSTANTS_

#define MAX_BITS 32
#define MAX_CODEWORD_LENGTH 16
#define MAX_CODE_NUM 256
#define PACK_SIZE 512

#define SEGMENT_SIZE 128
// #define SEGMENT_SIZE 256

// LOCAL_SEGMENT_NUM <= 4
#define LOCAL_SEGMENT_NUM 2
#define LOCAL_SEGMENT_SIZE SEGMENT_SIZE/LOCAL_SEGMENT_NUM

#define NUM_THREADS 128
#define WARP_SIZE 32
#define BLOCK_LAST NUM_THREADS/WARP_SIZE

#define GAP_LENGTH_MAX 4
#define GAP_FAC_NUM 8

#define SEPARATE_MIN 9

#define FLAG_A    0x0100000000000000
#define FLAG_P    0x8000000000000000
#define FLAG_MASK 0xFF00000000000000

#define X 0
#define Y 1
#define Z 2
#define W 3
#define O 4

#define LOOP 100

#define UINT_OUT( symbols, symbol, pos ) \
	symbols = (symbols|(symbol<<( pos*8 )));

#define FIXED_PREFIX_BIT 10

struct Symbol{
	unsigned char symbol;
	unsigned char length;
	unsigned int num;
};

struct Dectable{
	unsigned char symbol;
	unsigned char length;
};

struct TableInfo{
	unsigned int l1table_size;
	unsigned int l2table_size;
	unsigned int ptrtable_size;
};

#endif
