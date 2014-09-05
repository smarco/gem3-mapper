#include "myers-interface.h"

/********************************
Common constants for Device & Host
*********************************/

#define UINT32_ZEROS			0x00000000u
#define UINT32_ONES  			0xFFFFFFFFu
#define UINT32_ONE_LAST_MASK  	0x80000000u
#define	UINT32_SIZE				4

#define	CUDA_ERROR(error)		(CudaError(error, __FILE__, __LINE__ ))
#define	MYERS_ERROR(error)		(MyersError(error, __FILE__, __LINE__ ))
#define	MYERS_INLINE			inline

/* Defines related to Reference representation */
#define	REFERENCE_CHAR_LENGTH		4
#define	REFERENCE_CHARS_PER_ENTRY	(BMP_GPU_UINT32_LENGTH / REFERENCE_CHAR_LENGTH)
#define REFERENCE_END_PADDING		1250

/* Defines related to GPU Architecture */
#define	WARP_SIZE						32
#ifndef CUDA_NUM_THREADS
		#define		CUDA_NUM_THREADS	128
#endif

/* Defines to distribute and balance the work in GPU*/
#define	PEQ_LENGTH_PER_CUDA_THREAD		128
#define	NUM_BUCKETS_FOR_BINNING			(WARP_SIZE + 1)
#define MIN_CANDIDATES_PER_BUFFER		1000

/* Functions inline */
#define DIV_CEIL(NUMERATOR,DENOMINATOR) (((NUMERATOR)+((DENOMINATOR)-1))/(DENOMINATOR))
#define ROUND(NUM) ((int)((NUM) < 0 ? ((NUM) - 0.5) : ((NUM) + 0.5)))
#define MIN(NUM_A,NUM_B) (((NUM_A) < (NUM_B)) ? (NUM_A) : (NUM_B))

/* Conversion utils */
#define CONVERT_B_TO_KB(number) ((number)/(1024))
#define CONVERT_B_TO_MB(number) ((number)/(1024*1024))
#define CONVERT_B_TO_GB(number) ((number)/(1024*1024*1024))
#define CONVERT_MB_TO_B(number) ((number)*1024*1024)

/* System */
#define FILE_SIZE_LINES				250

/* Encoded DNA Nucleotides */
#define ENC_DNA_CHAR_A 0
#define ENC_DNA_CHAR_C 1
#define ENC_DNA_CHAR_G 2
#define ENC_DNA_CHAR_T 3

#define ENC_DNA_CHAR_N    4
#define ENC_DNA_CHAR_SEP  5
#define ENC_DNA_CHAR_JUMP 6

typedef enum
{
    SUCCESS,
    E_OPENING_FILE,
    E_READING_FILE,
	E_INSUFFICIENT_MEM_GPU,
	E_ALLOCATE_MEM,
	E_INCOMPATIBLE_GPU,
	E_NO_SUPPORTED_GPUS,
	E_REFERENCE_CODING
} _config_error;

/*Error type for the Myers API */
typedef _config_error myersError_t;


/*****************************
Internal Objects
*****************************/

typedef struct {
	uint64_t	size;
	uint64_t	numEntries;
	uint32_t	*h_reference;
	uint32_t	**d_reference;
} reference_buffer_t;

typedef struct {
	uint32_t	numResults;
 	uint32_t	numReorderedResults;
	bpm_gpu_res_entry_t	*h_results;
	bpm_gpu_res_entry_t	*h_reorderResults;
	bpm_gpu_res_entry_t	*d_reorderResults;
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
	bpm_gpu_cand_info_t 	*h_candidates;
	bpm_gpu_cand_info_t 	*d_candidates;
} candidates_buffer_t;


typedef struct {
	uint32_t			numDevices;
	uint32_t			numSupportedDevices;
	uint32_t			idDevice;
	uint32_t			idSupportedDevice;
	bpm_gpu_dev_arch_t	architecture;
	uint32_t			cudaCores;
	uint32_t			frequency;   /*Mhz*/
	float				relativePerformance;
} device_info_t;


/*************************************
Specific types for the Devices (GPUs)
**************************************/

typedef struct {
	uint4 bitmap[BMP_GPU_PEQ_ALPHABET_SIZE];
} d_qryEntry_t;

typedef struct {
	uint32_t		totalQueriesEntries;
	uint32_t		numQueries;
	d_qryEntry_t	*h_queries;
	d_qryEntry_t	*d_queries;
	bpm_gpu_qry_info_t		*h_qinfo;
	bpm_gpu_qry_info_t		*d_qinfo;
} d_queries_buffer_t;


/*************************************
Specific types for the Host (CPU)
**************************************/

typedef struct {
	uint32_t 				totalQueriesEntries;
	uint32_t 				numQueries;
	bpm_gpu_qry_entry_t 	*h_queries;
	bpm_gpu_qry_entry_t 	*d_queries;
	bpm_gpu_qry_info_t 				*h_qinfo;
	bpm_gpu_qry_info_t 				*d_qinfo;
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
	cudaStream_t		idStream;
	device_info_t		*device;
	reference_buffer_t 	*reference;
	queries_buffer_t 	*queries;
	candidates_buffer_t *candidates;
	reorder_buffer_t 	*reorderBuffer;
	results_buffer_t 	*results;
} buffer_t;


/*************************************
GPU Side defines (ASM instructions)
**************************************/

// output temporal carry in internal register
#define UADD__CARRY_OUT(c, a, b) \
     asm volatile("add.cc.u32 %0, %1, %2;" : "=r"(c) : "r"(a) , "r"(b));

// add & output with temporal carry of internal register
#define UADD__IN_CARRY_OUT(c, a, b) \
     asm volatile("addc.cc.u32 %0, %1, %2;" : "=r"(c) : "r"(a) , "r"(b));

// add with temporal carry of internal register
#define UADD__IN_CARRY(c, a, b) \
     asm volatile("addc.u32 %0, %1, %2;" : "=r"(c) : "r"(a) , "r"(b));
