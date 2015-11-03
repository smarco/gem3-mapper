
#ifndef GPU_FMI_CORE_H_
#define GPU_FMI_CORE_H_

#include "gpu_resources.h"
#include "gpu_commons.h"
#include "gpu_scheduler.h"

inline __device__ uint32_t count_bitmap(uint32_t bitmap, int32_t shift, uint32_t idxCounterGroup)
{
	uint32_t mask = GPU_UINT32_ONES << (GPU_UINT32_LENGTH - shift);

    mask = (shift > GPU_UINT32_LENGTH) ? GPU_UINT32_ONES : mask;
    mask = (shift > 0) ? mask : GPU_UINT32_ZEROS;

	mask = (idxCounterGroup) ? ~mask : mask;
    return (__popc(bitmap & mask));
}

inline __device__ uint32_t compute_bitmaps(uint4 bitmap, uint32_t bitmapPosition, uint32_t bit0, uint32_t bit1, uint32_t localEntryThreadIdx, uint32_t idxCounterGroup)
{
	const int32_t  relativePosition = bitmapPosition - (localEntryThreadIdx * GPU_UINT32_LENGTH);
	uint32_t resultBitmaps, bmpCollapsed, numCaracters;

	bitmap.x = bit0 ? bitmap.x : ~bitmap.x;
	bitmap.y = bit1 ? bitmap.y : ~bitmap.y;

	bmpCollapsed  = bitmap.x & bitmap.y & bitmap.z;
	resultBitmaps = count_bitmap(bmpCollapsed, relativePosition, idxCounterGroup);
	numCaracters  = reduce_add_warp(resultBitmaps, GPU_FMI_THREADS_PER_ENTRY);

	return (numCaracters);
}

inline __device__ uint64_t select_counter(uint4 vectorCounters, uint32_t indexBase)
{
	const uint2 vectorCountersA = {vectorCounters.x, vectorCounters.y};
	const uint2 vectorCountersB = {vectorCounters.z, vectorCounters.w};
	uint64_t scalarCountersA, scalarCountersB;

	V2S_B64(vectorCountersA, scalarCountersA);
	V2S_B64(vectorCountersB, scalarCountersB);

	return((indexBase == 0) ? scalarCountersA : scalarCountersB);
}

inline __device__ void gather_bitmaps(uint4 loadData, volatile gpu_fmi_exch_bmp_mem_t * exchBMP, uint32_t localEntryThreadIdx)
{
	const uint32_t idBMP = (localEntryThreadIdx == 0) ?  GPU_FMI_THREADS_PER_ENTRY - 1 : localEntryThreadIdx - 1;
	exchBMP->v[idBMP] = loadData.w;
}

inline __device__ uint64_t LF_Mapping(uint4 loadEntry, gpu_fmi_exch_bmp_mem_t * exchBMP, uint32_t missedEntry,
									  uint32_t localEntryThreadIdx, uint32_t bitmapPosition, uint32_t bit1, uint32_t bit0)
{
	// Select the counter candidate
    const uint64_t resultCounters = select_counter(loadEntry, bit0);

    // Reorganize entry layout between threads (send bitmaps to thread 0 using shared mem)
    gather_bitmaps(loadEntry, exchBMP, localEntryThreadIdx);
    if(localEntryThreadIdx == 0) loadEntry = exchBMP->s;

    // Count the number of occ in the bitmap
	const uint32_t resultBitmaps  = compute_bitmaps(loadEntry, bitmapPosition, bit0, bit1, localEntryThreadIdx, missedEntry);

	// Compute interval alternate counters
	const uint64_t interval = (missedEntry) ? resultCounters - resultBitmaps : resultCounters + resultBitmaps;

	return(interval);
}

#endif /* GPU_FMI_CORE_H_ */

