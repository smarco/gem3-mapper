#include <stdint.h>

#define		NUM_BITS			4
#define		NUM_BASES			5
#define		NUM_BASES_ENTRY		128
#define     SIZE_GPU_HW_WORD   	32
#define		NUM_SUB_ENTRIES 	(NUM_BASES_ENTRY / SIZE_GPU_HW_WORD)
#define     MAX_VALUE       	0xFFFFFFFF

#define CUDA_ERROR(error) (CudaError(error, __FILE__, __LINE__ ))
#define MYERS_ERROR(error) (MyersError(error, __FILE__, __LINE__ ))

typedef enum
{
    SUCCESS,
    E_OPENING_FILE,
    E_READING_FILE,
	E_INSUFFICIENT_MEM_GPU,
	E_ALLOCATE_MEM,
	E_INCOMPATIBLE_GPU
} _config_error;

/*Error type for the Myers API */
typedef _config_error myersError_t;


/*****************************
Common types for Device & Host
*****************************/

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

typedef struct {
	uint32_t bitmap[NUM_BASES][NUM_SUB_ENTRIES];
} qryEntry_t;


/*****************************
Internal Objects
*****************************/

typedef struct {
	uint32_t size;
	uint32_t numEntries; 
	uint32_t *h_reference;
	uint32_t *d_reference;
} reference_buffer_t;

typedef struct {
	uint32_t numResults;
 	uint32_t numReorderedResults;
	resEntry_t* h_results;
	resEntry_t* h_reorderResults;
	resEntry_t* d_reorderResults;
} results_buffer_t;

typedef struct {
	uint32_t	numBuckets;
	uint32_t	candidatesPerBufffer;
	uint32_t	numWarps;
	uint32_t   	*h_reorderBuffer;
	uint32_t   	*d_reorderBuffer;
	uint32_t 	*h_initPosPerBucket;
	uint32_t 	*h_initWarpPerBucket;
	uint32_t 	*d_initPosPerBucket;
	uint32_t 	*d_initWarpPerBucket;
} reorder_buffer_t;

typedef struct {
	uint32_t 	numCandidates;
	candInfo_t 	*h_candidates;
	candInfo_t 	*d_candidates;
} candidates_buffer_t;



/*************************************
Specific types for the Devices (GPUs)
**************************************/

typedef struct {
	uint4 bitmap[NUM_BASES];
} d_qryEntry_t;

typedef struct {
	uint32_t 	totalQueriesEntries;
	uint32_t 	numQueries;
	d_qryEntry_t 	*h_queries;
	d_qryEntry_t 	*d_queries;
	qryInfo_t 		*h_qinfo;
	qryInfo_t 		*d_qinfo;
} d_queries_buffer_t;



/*************************************
Specific types for the Host (CPU)
**************************************/

typedef struct {
	uint32_t 	totalQueriesEntries;
	uint32_t 	numQueries;
	float		distance;
	qryEntry_t 	*h_queries;
	qryEntry_t 	*d_queries;
	qryInfo_t 	*h_qinfo;
	qryInfo_t 	*d_qinfo;
} queries_buffer_t;



/*************************************
Interface Objects
**************************************/

typedef struct {
	uint32_t			numBuffers;
	uint32_t 			idBuffer;
	uint32_t 			maxPEQEntries;
	uint32_t			maxCandidates;
	uint32_t			maxQueries;
	uint32_t			maxReorderBuffer;
	uint32_t			idDevice;
	uint32_t			majorCC;
	uint32_t			minorCC;
	cudaStream_t		idStream;
	reference_buffer_t 	*reference;
	queries_buffer_t 	*queries;
	candidates_buffer_t *candidates;
	reorder_buffer_t 	*reorderBuffer;
	results_buffer_t 	*results;
} buffer_t;

