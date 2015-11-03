#ifndef GPU_REFERENCE_H_
#define GPU_REFERENCE_H_

#include "gpu_commons.h"

/* Defines related to Reference representation */
#define	GPU_REFERENCE_CHAR_LENGTH		4
#define GPU_REFERENCE_CHARS_PER_UINT1	(GPU_UINT32_LENGTH / GPU_REFERENCE_CHAR_LENGTH)
#define GPU_REFERENCE_CHARS_PER_UINT2	(GPU_REFERENCE_CHARS_PER_UINT1 * 2)
#define GPU_REFERENCE_CHARS_PER_UINT4	(GPU_REFERENCE_CHARS_PER_UINT1 * 4)
#define	GPU_REFERENCE_CHARS_PER_ENTRY	GPU_REFERENCE_CHARS_PER_UINT2
#define GPU_REFERENCE_BYTES_PER_ENTRY	GPU_UINT64_SIZE
#define GPU_REFERENCE_END_PADDING		625

/*****************************
Internal Objects (General)
*****************************/

typedef struct {
	uint64_t			size;
	uint64_t			numEntries;
	uint64_t			*h_reference;
	uint64_t			**d_reference;
	memory_alloc_t		*memorySpace;
} gpu_reference_buffer_t;


#endif /* GPU_REFERENCE_H_ */
