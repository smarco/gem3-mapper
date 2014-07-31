//#include "myers-common.h"
#include <stdint.h>

#define		NUM_BITS			4
#define		NUM_BASES			5
#define		NUM_BASES_ENTRY		128
#define     SIZE_GPU_HW_WORD   	32
#define		NUM_SUB_ENTRIES 	(NUM_BASES_ENTRY / SIZE_GPU_HW_WORD)

typedef struct {
	uint32_t bitmap[NUM_BASES][NUM_SUB_ENTRIES];
} qryEntry_t;

typedef struct {
	uint32_t column;
	uint32_t score;
} resEntry_t;

typedef struct {
	uint32_t query;
	uint64_t position;
} candInfo_t;

typedef struct {
	uint32_t posEntry;
	uint32_t size;
} qryInfo_t;


//Obtain Buffers
inline qryEntry_t* getPEQBuffer	   	    (void *myersBuffer);
inline candInfo_t* getCandidatesBuffer  (void *myersBuffer);
inline qryInfo_t*  getPEQInfoBuffer 	(void *myersBuffer);
inline resEntry_t* getResultsBuffer 	(void *myersBuffer);

//Get elements
inline uint32_t getMaxPEQEntries  (void *myersBuffer);
inline uint32_t getMaxCandidates  (void *myersBuffer);
inline uint32_t getMaxQueries  (void *myersBuffer);

//Main functions
int initMyers			(void ***myersBuffer, uint32_t numBuffers, const char* referenceFileName, 
			  			 int32_t averageQuerySize, int32_t candidatesPerQuery);
void sendMyersBuffer	(void *myersBuffer, float distance, uint32_t numPEQEntries, uint32_t numQueries, uint32_t numCandidates);
void receiveMyersBuffer	(void *myersBuffer);
void endMyers			(void ***myersBuffer);

