/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <nmmintrin.h>
#include <time.h>
#include <sys/time.h>
#include <limits.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#include <cuda.h>
#include <cuda_runtime.h>

#define	BWT_CHAR_LENGTH				    3 			    // 3 bits / character

#define FMI_NUM_COUNTERS			    4 			    // 4 virtual counters
#define FMI_ALTERNATE_COUNTERS		2 			    // 2
#define FMI_ENTRY_SIZE				    128			    // 128 bases / FMI entry
#define FMI_BITS_LOAD_PER_THREAD	128			    // accesses of 16 Bytes / thread

#define	UINT32_LENGTH				      32 			    // 32 bits
#define	UINT64_LENGTH				      64 			    // 64 bits
#define	UINT64_MAX_VALUE			    ULONG_MAX   // 64 bits
#define MASK_ONES					        0xFFFFFFFF
#define MASK_ONE					        0x00000001
#define MASK_ZEROS					      0x00000000

#define WARP_SIZE					        32
#define MAX_THREADS_PER_SM			  128

#define FMI_COUNTERS_PER_ENTRY		(FMI_NUM_COUNTERS / FMI_ALTERNATE_COUNTERS)									                    //   2 physical counters / FMI entry
#define FMI_BITMAPS_PER_ENTRY		  (FMI_ENTRY_SIZE * BWT_CHAR_LENGTH / UINT32_LENGTH)							                //   12 physical bitmaps / FMI entry
#define FMI_ENTRY_LENGTH			    ((FMI_COUNTERS_PER_ENTRY * UINT64_LENGTH) + (FMI_ENTRY_SIZE * BWT_CHAR_LENGTH)) //   512 bits
#define FMI_THREADS_PER_ENTRY		  (FMI_ENTRY_LENGTH / FMI_BITS_LOAD_PER_THREAD)  								                  //   4 threads
#define FMI_ENTRIES_PER_WARP		  (WARP_SIZE / FMI_THREADS_PER_ENTRY)
#define FMI_ENTRIES_PER_BLOCK		  (MAX_THREADS_PER_SM / FMI_THREADS_PER_ENTRY)

#define FMI_ENTRIES_PER_DECODE		2 											       // Necessary FMI entries for unaligned FMI threads
#define FMI_THREADS_PER_COUNTERS	2 											       // Necessary threads to load all the counters
#define DECODE_THREADS_PER_ENTRY	(FMI_THREADS_PER_ENTRY * FMI_ENTRIES_PER_DECODE)   // FMI_THREADS_PER_ENTRY + 1 (rounded to pow of 2)
#define DECODE_THREADS_PER_LOAD		(FMI_THREADS_PER_ENTRY + FMI_THREADS_PER_COUNTERS) // Threads to load an entry for the decode primitive

#define V2S_B64(v,s)	asm("mov.b64 %0, {%1,%2};" : "=l"(s) : "r"(v.x), "r"(v.y))
#define S2V_B64(s,v)	asm("mov.b64 {%0,%1}, %2;" : "=r"(v.x), "=r"(v.y) : "l"(s))

#define MIN(NUM_A, NUM_B) 				      ((NUM_A < NUM_B) ? NUM_A : NUM_B)
#define CATCH_ERROR(error) 				      {{if (error) { fprintf(stderr, "%s\n", processError(error)); exit(EXIT_FAILURE); }}}
#define DIV_CEIL(NUMERATOR,DENOMINATOR)	(((NUMERATOR)+((DENOMINATOR)-1))/(DENOMINATOR))


/* Encoded DNA Nucleotides */
#define ENC_DNA_CHAR_A    0
#define ENC_DNA_CHAR_C    1
#define ENC_DNA_CHAR_G    2
#define ENC_DNA_CHAR_T    3


// FMI Entry (64 Bytes) using:
//    4 letters: Alternate counters   (2  uint64_t)
//    5 letters: 12 Bitmaps x 32 bits (3 uint128_t)
typedef struct {
	uint64_t counters[FMI_COUNTERS_PER_ENTRY];
	uint32_t bitmaps[FMI_ENTRY_SIZE * BWT_CHAR_LENGTH / UINT32_LENGTH];
} scalar_fmi_entry_t;

typedef uint4 vector_fmi_entry_t;

typedef union {
  scalar_fmi_entry_t s;
  vector_fmi_entry_t v[FMI_THREADS_PER_ENTRY];
} fmi_entry_t;

typedef union {
  uint4 s;
  uint32_t v[FMI_BITS_LOAD_PER_THREAD / UINT32_LENGTH];
} exch_bmp_mem_t;

typedef struct {
	uint32_t        numDecodings;
	uint64_t       *h_initBWTPos;
	uint64_t       *d_initBWTPos;
	uint64_t       *h_endBWTPos;
	uint32_t       *h_steps;
} gem_profile_t;

typedef struct {
	uint32_t        numDecodings;
	ulonglong2     *h_BWTPos;
	ulonglong2     *d_BWTPos;
	ulonglong2     *h_BWTPos_host;
} dec_results_t;

typedef struct {
	uint64_t        bwtSize;
	uint64_t        numEntries;
	fmi_entry_t    *h_fmi;
	fmi_entry_t    *d_fmi;
} fmi_t;

inline double sampleTime()
{
    struct timespec tv;
    #ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
        clock_serv_t cclock;
        mach_timespec_t mts;
        host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
        clock_get_time(cclock, &mts);
        mach_port_deallocate(mach_task_self(), cclock);
        tv.tv_sec = mts.tv_sec;
        tv.tv_nsec = mts.tv_nsec;
    #else
        clock_gettime(CLOCK_REALTIME, &tv);
    #endif
    return((tv.tv_sec+tv.tv_nsec/1000000000.0));
}

#define CUDA_ERROR(error) (HandleError(error, __FILE__, __LINE__ ))

static void HandleError( cudaError_t err, const char *file,  int32_t line ) {
  if (err != cudaSuccess) {
    printf( "%s in %s at line %d\n", cudaGetErrorString(err),  file, line );
    exit( EXIT_FAILURE );
  }
}

const char * processError(uint32_t e){
  switch(e) {
    case 0:  return "No error";
    case 30: return "Cannot open reference file";
    case 31: return "Cannot allocate reference";
    case 32: return "Reference file isn't multifasta format";
    case 37: return "Cannot open reference file on write mode";
    case 42: return "Cannot open queries file";
    case 43: return "Cannot allocate queries";
    case 45: return "Cannot allocate results";
    case 47: return "Cannot open results file for save intervals";
    case 48: return "Cannot open results file for load intervals";
    case 99: return "Not implemented";
    default: return "Unknown error";
  }
}

uint64_t charToBinASCII(unsigned char base)
{
	switch(base)
	{
    case 'A':
    case 'a':
      return(ENC_DNA_CHAR_A);
    case 'C':
    case 'c':
      return(ENC_DNA_CHAR_C);
    case 'G':
    case 'g':
      return(ENC_DNA_CHAR_G);
    case 'T':
    case 't':
      return(ENC_DNA_CHAR_T);
    default :
      return(ENC_DNA_CHAR_A);
	}
}

inline __device__ uint32_t countBitmap(uint32_t bitmap, int32_t shift, uint32_t idxCounterGroup)
{
  uint32_t mask = MASK_ONES << (UINT32_LENGTH - shift);

  mask = (shift > UINT32_LENGTH) ? MASK_ONES : mask;
  mask = (shift > 0) ? mask : MASK_ZEROS;

  mask = (idxCounterGroup) ? ~mask : mask;
  return (__popc(bitmap & mask));
}

inline __device__ uint32_t reduceEntry(uint32_t resultBitmaps)
{
	for (int32_t i = 1; i < FMI_THREADS_PER_ENTRY; i *= 2){
		int32_t n = __shfl_down((int) resultBitmaps, i, WARP_SIZE);
		resultBitmaps += n;
	}
	return(resultBitmaps);
}

inline __device__ uint32_t computeBitmaps(uint4 bitmap, uint32_t bitmapPosition, uint32_t bit0, uint32_t bit1, uint32_t localEntryThreadIdx, uint32_t idxCounterGroup)
{
	const int32_t  relativePosition = bitmapPosition - (localEntryThreadIdx * UINT32_LENGTH);
	uint32_t resultBitmaps, bmpCollapsed, numCaracters;

	bitmap.x = bit0 ? bitmap.x : ~bitmap.x;
	bitmap.y = bit1 ? bitmap.y : ~bitmap.y;

	bmpCollapsed  = bitmap.x & bitmap.y & bitmap.z;
	resultBitmaps = countBitmap(bmpCollapsed, relativePosition, idxCounterGroup);
	numCaracters  = reduceEntry(resultBitmaps);

	return (numCaracters);
}

inline __device__ uint64_t selectCounter(uint4 vectorCounters, uint32_t indexBase)
{
	const uint2 vectorCountersA = {vectorCounters.x, vectorCounters.y};
	const uint2 vectorCountersB = {vectorCounters.z, vectorCounters.w};
	uint64_t scalarCountersA, scalarCountersB;
	
	V2S_B64(vectorCountersA, scalarCountersA);
	V2S_B64(vectorCountersB, scalarCountersB);

	return((indexBase == 0) ? scalarCountersA : scalarCountersB);
}

inline __device__ uint64_t shfl_64(uint64_t scalarValue, int lane)
{
	uint2 vectorValue;
	S2V_B64(scalarValue, vectorValue);
	vectorValue.x = __shfl(vectorValue.x, lane);
	vectorValue.y = __shfl(vectorValue.y, lane);
	V2S_B64(vectorValue, scalarValue);
	return (scalarValue);
}

inline __device__ void gatherBitmaps(uint4 loadData, volatile exch_bmp_mem_t * exchBMP, uint32_t localEntryThreadIdx)
{
	const uint32_t idBMP = (localEntryThreadIdx == 0) ?  FMI_THREADS_PER_ENTRY - 1 : localEntryThreadIdx - 1;
	exchBMP->v[idBMP] = loadData.w;
}

inline uint32_t countBitmapCPU(uint32_t bitmap, int32_t shift, uint32_t idxCounterGroup)
{
  uint32_t mask = MASK_ONES << (UINT32_LENGTH - shift);

  mask = (shift > UINT32_LENGTH) ? MASK_ONES : mask;
  mask = (shift > 0) ? mask : MASK_ZEROS;

  mask = (idxCounterGroup) ? ~mask : mask;
  return (_mm_popcnt_u32(bitmap & mask));
}

inline uint32_t computeBitmapsCPU(uint3 vbitmap, uint32_t bitmapPosition, uint32_t bit0, uint32_t bit1, uint32_t missedEntry, uint32_t idBitmap)
{
	const int32_t  relativePosition = bitmapPosition - (idBitmap * UINT32_LENGTH);
	uint32_t bmpCollapsed, numCaracters;

	vbitmap.x = bit0 ? vbitmap.x : ~vbitmap.x;
	vbitmap.y = bit1 ? vbitmap.y : ~vbitmap.y;

	bmpCollapsed  = vbitmap.x & vbitmap.y & vbitmap.z;
	numCaracters = countBitmapCPU(bmpCollapsed, relativePosition, missedEntry);

	return (numCaracters);
}

inline uint32_t gatherBaseFromBWTCPU(uint3 vbitmap, uint32_t bitmapPosition, uint32_t * bit1, uint32_t * bit0, uint32_t * bit2 /*, uint32_t idDecoding*/)
{
	const int32_t  relativePosition = (UINT32_LENGTH - (bitmapPosition % UINT32_LENGTH)) - 1;
	(* bit0) = (vbitmap.x >> relativePosition) & 1;
	(* bit1) = (vbitmap.y >> relativePosition) & 1;
	(* bit2) = (vbitmap.z >> relativePosition) & 1;
	return (0);
}

char decodeBase(uint32_t bit0, uint32_t bit1)
{
	if(bit1 == 0 && bit0 == 0) return ('A');
	if(bit1 == 0 && bit0 == 1) return ('C');
	if(bit1 == 1 && bit0 == 0) return ('G');
	if(bit1 == 1 && bit0 == 1) return ('T');	
	return ('A');
}

void decodingPositionsFMICPUKernel(fmi_entry_t *fmi, uint64_t bwtSize, uint32_t numDecodings, uint64_t *h_initBWTPos, ulonglong2 *h_BWTPos, uint32_t samplingRate)
{
	const uint32_t NUM_BITMAPS = FMI_ENTRY_SIZE / UINT32_LENGTH; 
	const uint32_t LUT[12] = {3,7,11,0,1,2,4,5,6,8,9,10};

	for(uint32_t idDecoding = 0; idDecoding < numDecodings; ++idDecoding){
    uint64_t interval = h_initBWTPos[idDecoding];
    uint32_t idStep = 0, bit0, bit1, bit2;

		while(interval % samplingRate){
			const uint64_t entryIdx	      = interval / FMI_ENTRY_SIZE;
			const uint32_t bitmapPosition = interval % FMI_ENTRY_SIZE;
			      uint32_t numCharacters  = 0;

			// Gathering the base of the decodification from the BWT
			const uint32_t basePositionBMP = (bitmapPosition / UINT32_LENGTH) * BWT_CHAR_LENGTH;
			uint3 vbitmapBase  = { fmi[entryIdx].s.bitmaps[LUT[basePositionBMP]], fmi[entryIdx].s.bitmaps[LUT[basePositionBMP + 1]], fmi[entryIdx].s.bitmaps[LUT[basePositionBMP + 2]]};
			gatherBaseFromBWTCPU(vbitmapBase, bitmapPosition, &bit1, &bit0, &bit2);

			const uint32_t missedEntry    = (entryIdx % FMI_ALTERNATE_COUNTERS == bit1) ? 0 : 1;
			const uint64_t bigCounter     = fmi[entryIdx + missedEntry].s.counters[bit0];
			// reorder bitmaps layout for low flag
			for(uint32_t idBitmap = 0; idBitmap < NUM_BITMAPS; ++idBitmap){
			 	const uint32_t initBitmap  = idBitmap * BWT_CHAR_LENGTH;
			 	uint3 vbitmap  = {fmi[entryIdx].s.bitmaps[LUT[initBitmap]], fmi[entryIdx].s.bitmaps[LUT[initBitmap + 1]], fmi[entryIdx].s.bitmaps[LUT[initBitmap + 2]]};
			 	numCharacters += computeBitmapsCPU(vbitmap, bitmapPosition, bit0, bit1, missedEntry, idBitmap);
			}

			interval = (missedEntry) ? bigCounter - numCharacters : bigCounter + numCharacters;
			interval = (bit2 == 0) ? UINT64_MAX_VALUE : interval;
			idStep++;
		}

		h_BWTPos[idDecoding].x = interval;
	}
}

inline __device__ uint64_t LF_Mapping(uint4 loadEntry, exch_bmp_mem_t * exchBMP, uint32_t missedEntry, 
									  uint32_t localEntryThreadIdx, uint32_t bitmapPosition, uint32_t bit1, uint32_t bit0)
{
  // Select the counter candidate
  const uint64_t resultCounters = selectCounter(loadEntry, bit0);

  // Reorganize entry layout between threads (send bitmaps to thread 0 using shared mem)
  gatherBitmaps(loadEntry, exchBMP, localEntryThreadIdx);
  if(localEntryThreadIdx == 0) loadEntry = exchBMP->s;

  // Count the number of occ in the bitmap
  const uint32_t resultBitmaps  = computeBitmaps(loadEntry, bitmapPosition, bit0, bit1, localEntryThreadIdx, missedEntry);

  // Compute interval alternate counters
  const uint64_t interval = (missedEntry) ? resultCounters - resultBitmaps : resultCounters + resultBitmaps;

  return(interval);
}


inline __device__ void gatherBaseFromBWT(uint4 vbitmap, exch_bmp_mem_t * exchBMP, uint32_t bitmapPosition, 
										 uint32_t fmiEntryThreadIdx, uint32_t decodeEntryThreadIdx, uint32_t decodeEntryIdx,
										 uint32_t * bit1, uint32_t * bit0, uint32_t * bit2 /*, uint32_t idDecoding*/)
{
  uint32_t localBit0, localBit1, localBit2;

  gatherBitmaps(vbitmap, exchBMP, fmiEntryThreadIdx);
  if(fmiEntryThreadIdx == 0) vbitmap = exchBMP->s;

  const int32_t  relativePosition = (UINT32_LENGTH - (bitmapPosition % UINT32_LENGTH)) - 1;
  localBit0 = (vbitmap.x >> relativePosition) & 1;
  localBit1 = (vbitmap.y >> relativePosition) & 1;
  localBit2 = (vbitmap.z >> relativePosition) & 1;

  const uint32_t lane = (decodeEntryIdx * DECODE_THREADS_PER_ENTRY) + (bitmapPosition / UINT32_LENGTH);
  (* bit0) = __shfl(localBit0, lane);
  (* bit1) = __shfl(localBit1, lane);
  (* bit2) = __shfl(localBit2, lane);
}

inline __device__ uint4 gatherCountersFMI(uint4 loadEntry, uint32_t missedEntry, uint32_t decodeEntryIdx, uint32_t decodeEntryThreadIdx)
{
	uint4 auxData;
	const uint32_t idThreadCounter = (decodeEntryIdx * DECODE_THREADS_PER_ENTRY) + FMI_THREADS_PER_ENTRY;
	const uint32_t lane = (missedEntry) ? idThreadCounter + 1 : idThreadCounter;

 	auxData.x = __shfl(loadEntry.x, lane);
 	auxData.y = __shfl(loadEntry.y, lane);
 	auxData.z = __shfl(loadEntry.z, lane);
 	auxData.w = __shfl(loadEntry.w, lane);

	if(decodeEntryThreadIdx == 0) loadEntry = auxData;

	return(loadEntry);
}

inline __device__ uint32_t generateLoadThreadIdx(uint64_t decodeEntryThreadIdx)
{
	uint32_t localThreadIdx = (decodeEntryThreadIdx  > (FMI_THREADS_PER_ENTRY + 1)) ? 0                     : decodeEntryThreadIdx;
	         localThreadIdx = (decodeEntryThreadIdx == (FMI_THREADS_PER_ENTRY + 1)) ? FMI_THREADS_PER_ENTRY : localThreadIdx;
	         localThreadIdx = (decodeEntryThreadIdx ==  FMI_THREADS_PER_ENTRY)      ? 0                     : localThreadIdx;
  return(localThreadIdx);
}

void __global__ decodingPositionsFMIGPUKernel(fmi_entry_t *fmi, uint64_t bwtSize, uint32_t numDecodings, uint64_t *d_initBWTPos, ulonglong2 *d_BWTPos, uint32_t samplingRate)
{
	const uint32_t globalThreadIdx      = blockIdx.x * MAX_THREADS_PER_SM + threadIdx.x;
	const uint32_t localWarpThreadIdx   = globalThreadIdx % WARP_SIZE;
	const uint32_t idDecoding 	 	    = globalThreadIdx / DECODE_THREADS_PER_ENTRY;

	if((threadIdx.x < MAX_THREADS_PER_SM) && (idDecoding < numDecodings)){

		const uint32_t decodeEntryIdx        = localWarpThreadIdx  / DECODE_THREADS_PER_ENTRY;
		const uint32_t decodeEntryThreadIdx  = localWarpThreadIdx  % DECODE_THREADS_PER_ENTRY;
		const uint32_t fmiEntryThreadIdx     = localWarpThreadIdx  % FMI_THREADS_PER_ENTRY;
		const uint32_t loadFMIEntryThreadIdx = generateLoadThreadIdx(decodeEntryThreadIdx);

          uint4    loadEntry;
          uint64_t interval = d_initBWTPos[idDecoding], idStep = 0;
          uint32_t foundBaseN = 0;

		__shared__ exch_bmp_mem_t   exchBMP[FMI_ENTRIES_PER_BLOCK];
				       exch_bmp_mem_t * decExchBMP = &exchBMP[threadIdx.x / FMI_THREADS_PER_ENTRY];

		while((interval % samplingRate) && (foundBaseN == 0)){			
			const uint64_t entryIdx    	  =  interval / FMI_ENTRY_SIZE;
			const uint32_t bitmapPosition =  interval % FMI_ENTRY_SIZE;
			      uint32_t bit1, bit0, bit2;

			// Loading FM-index entry in thread cooperative way
			if(( 0 < decodeEntryThreadIdx) && (decodeEntryThreadIdx < DECODE_THREADS_PER_LOAD)) 
				loadEntry = fmi[entryIdx].v[loadFMIEntryThreadIdx];

			// Gathering the base and sharing it with the rest of the threads  
			gatherBaseFromBWT(loadEntry, decExchBMP, bitmapPosition, fmiEntryThreadIdx, decodeEntryThreadIdx, decodeEntryIdx, &bit1, &bit0, &bit2);

			// Gathering the counters
			const uint32_t missedEntry = (entryIdx % FMI_ALTERNATE_COUNTERS != bit1) ? 1 : 0;
 			loadEntry = gatherCountersFMI(loadEntry, missedEntry, decodeEntryIdx, decodeEntryThreadIdx);

			// Compute LF-Mapping th0 of each group contain the result)
 			interval = LF_Mapping(loadEntry, decExchBMP, missedEntry, fmiEntryThreadIdx, bitmapPosition, bit1, bit0);

			// Share interval along the thread group
 			const uint32_t lane = decodeEntryIdx * DECODE_THREADS_PER_ENTRY;
			interval 		        = shfl_64(interval, lane);

			// Exit condition
			if(bit2 == 0) foundBaseN = 1;

			// Increment for the next Backward-Search
			idStep++;
		}

		ulonglong2 BWTPos = {interval, idStep};
		if(foundBaseN) BWTPos = make_ulonglong2(UINT64_MAX_VALUE, UINT64_MAX_VALUE);
		
		// Save intervals
		if(decodeEntryThreadIdx  == 0) d_BWTPos[idDecoding] = BWTPos;
	}
}

uint32_t loadFMI(const char *fn, fmi_t *fmi)
{
  FILE *fp = NULL;

  fp = fopen(fn, "rb");
  if (fp == NULL) return (8);

  fread(&fmi->numEntries, sizeof(uint64_t), 1, fp);
  fread(&fmi->bwtSize, sizeof(uint64_t), 1, fp);

  fmi->h_fmi = (fmi_entry_t *) malloc(fmi->numEntries * sizeof(fmi_entry_t));
  fread(fmi->h_fmi, sizeof(fmi_entry_t), fmi->numEntries, fp);
  fclose(fp);

  return (0);
}

uint32_t loadGEMProfile(const char *fn, gem_profile_t *profRegions)
{
  FILE *fp = NULL;

  uint32_t steps;
  uint64_t initBWTPos, endBWTPos;
  uint32_t bookmark, numAs, numCs, numGs, numTs, numXs;

  fp = fopen(fn, "r");
  if (fp == NULL) return (8);

  //allocate the memory for the structures on profRegions
  profRegions->h_initBWTPos = (uint64_t *) malloc(profRegions->numDecodings * sizeof(uint64_t));
  profRegions->h_endBWTPos  = (uint64_t *) malloc(profRegions->numDecodings * sizeof(uint64_t));
  profRegions->h_steps      = (uint32_t *) malloc(profRegions->numDecodings * sizeof(uint32_t));

  for(uint32_t idDecoding = 0; idDecoding < profRegions->numDecodings; ++idDecoding){
    fscanf(fp, "%llu %llu %u %c %u %u %u %u %u",
           &initBWTPos, &endBWTPos, &steps, &bookmark,
           &numAs, &numCs, &numGs, &numTs, &numXs);
    profRegions->h_initBWTPos[idDecoding] = initBWTPos;
    profRegions->h_endBWTPos[idDecoding]  = endBWTPos;
    profRegions->h_steps[idDecoding]      = steps;
  }

  fclose(fp);
  return (0);
}

uint32_t allocateResults(dec_results_t *dec, uint32_t numDecodings)
{
	dec->numDecodings   = numDecodings;
  dec->h_BWTPos       = (ulonglong2 *) malloc(dec->numDecodings * sizeof(ulonglong2));
  dec->h_BWTPos_host  = (ulonglong2 *) malloc(dec->numDecodings * sizeof(ulonglong2));
  return (0);
}


uint32_t inspectGEMProfile(const char *fn, gem_profile_t *profRegions)
{
  FILE *fp = NULL;
  const uint32_t NUM_ELEMENTS_PER_LINE = 9;
        uint32_t parsedElements = NUM_ELEMENTS_PER_LINE;

  uint32_t steps, numDecodings = 0;
  uint32_t bookmark, numAs, numCs, numGs, numTs, numXs;
  uint64_t initBWTPos, endBWTPos;

  fp = fopen(fn, "r");
  if (fp == NULL) return (8);

  while(parsedElements == NUM_ELEMENTS_PER_LINE){
		parsedElements = fscanf(fp, "%llu %llu %u %c %u %u %u %u %u", 
		                        &initBWTPos, &endBWTPos, &steps, &bookmark, &numAs, &numCs, &numGs, &numTs, &numXs);
		if(parsedElements == NUM_ELEMENTS_PER_LINE) numDecodings++;
	}

	profRegions->numDecodings = numDecodings;

  fclose(fp);
  return (0);
}

inline char BinASCIItoChar(uint32_t base)
{
	switch(base)
	{
    case ENC_DNA_CHAR_A:
      return('A');
    case ENC_DNA_CHAR_C:
      return('C');
    case ENC_DNA_CHAR_G:
      return('G');
    case ENC_DNA_CHAR_T:
      return('T');
    default :
      return('X');
	}
}

inline uint32_t printGEMProfile(gem_profile_t *profRegions, const uint32_t maxDecodings)
{
  for(uint32_t idDecoding = 0; idDecoding < MIN(profRegions->numDecodings, maxDecodings); ++idDecoding){
    printf("[%d] init_pos=%llu   \t end_pos=%llu \t steps=%u \n", idDecoding, profRegions->h_initBWTPos[idDecoding],
           profRegions->h_endBWTPos[idDecoding], profRegions->h_steps[idDecoding]);
  }
  return (0);
}

uint32_t checkIntervalsGPU(gem_profile_t *profRegions, dec_results_t *results)
{
	const uint32_t maxDecodings = 10; // Just check and print the first results
	      uint32_t missMatches  = 0;

  for(uint32_t idDecoding = 0; idDecoding < profRegions->numDecodings; ++idDecoding){
    if(profRegions->h_endBWTPos[idDecoding] != results->h_BWTPos[idDecoding].x){
      if(missMatches < maxDecodings){
      printf("[%d] init_pos=%llu   \t end_pos=%llu \t steps=%u \t\t (Results) end_pos=%llu steps=%llu \n",
        idDecoding, profRegions->h_initBWTPos[idDecoding], profRegions->h_endBWTPos[idDecoding],
        profRegions->h_steps[idDecoding], results->h_BWTPos[idDecoding].x, results->h_BWTPos[idDecoding].y);
      }
      missMatches++;
    }
	}
  return (0);
}

uint32_t checkIntervalsCPU(gem_profile_t *profRegions, dec_results_t *results)
{
	const uint32_t maxDecodings = 10; // Just check and print the first results
	      uint32_t missMatches  = 0, unResolvedDecoding = 0;

  for(uint32_t idDecoding = 0; idDecoding < profRegions->numDecodings; ++idDecoding){
    if(profRegions->h_endBWTPos[idDecoding] != results->h_BWTPos_host[idDecoding].x
      && profRegions->h_steps[idDecoding] != results->h_BWTPos_host[idDecoding].y){
      if(results->h_BWTPos_host[idDecoding].x == 0) unResolvedDecoding++;
      if(missMatches < maxDecodings){
        printf("[%d] init_pos=%llu   \t end_pos=%llu \t steps=%u \t\t (Results) end_pos=%llu steps=%llu \n",
              idDecoding, profRegions->h_initBWTPos[idDecoding], profRegions->h_endBWTPos[idDecoding],
              profRegions->h_steps[idDecoding], results->h_BWTPos_host[idDecoding].x, results->h_BWTPos_host[idDecoding].y);
      }
      missMatches++;
    }
	}

	printf("=> Non solved decodings: %u \n", unResolvedDecoding);
	return (0);
}

uint32_t initProfileGEM(gem_profile_t **profRegions)
{
  gem_profile_t *prof  = (gem_profile_t *) malloc(sizeof(gem_profile_t));
	prof->numDecodings   = 0;

	prof->h_initBWTPos   = NULL;
	prof->d_initBWTPos   = NULL;
	prof->h_endBWTPos    = NULL;
	prof->h_steps        = NULL;

  (* profRegions) = prof;
  return (0);
}

uint32_t initFMI(fmi_t **fm_index)
{
  fmi_t *fmi  	  = (fmi_t *) malloc(sizeof(fmi_t));
	fmi->bwtSize    = 0;
	fmi->numEntries = 0;

	fmi->h_fmi      = NULL;
	fmi->d_fmi      = NULL;

  (* fm_index)    = fmi;
  return (0);
}

uint32_t initResults(dec_results_t **decoding)
{
  dec_results_t *dec = (dec_results_t *) malloc(sizeof(dec_results_t));
	dec->numDecodings  = 0;

	dec->h_BWTPos      = NULL;
	dec->d_BWTPos      = NULL;
	dec->h_BWTPos_host = NULL;

  (* decoding) = dec;
  return (0);
}

uint32_t freeFMI(fmi_t **fmi)
{   
  if((* fmi)->h_fmi != NULL){
    free((* fmi)->h_fmi);
    (* fmi)->h_fmi = NULL;
  }

  if((* fmi)->d_fmi != NULL){
    cudaFree((* fmi)->d_fmi);
    (* fmi)->d_fmi = NULL;
  }

  if((* fmi) != NULL){
    free(* fmi);
    (* fmi) = NULL;
  }
  return(0);
}


uint32_t freeResults(dec_results_t **dec)
{   
  if((* dec)->h_BWTPos_host != NULL){
    free((* dec)->h_BWTPos_host);
    (* dec)->h_BWTPos_host = NULL;
  }

  if((* dec)->h_BWTPos != NULL){
    free((* dec)->h_BWTPos);
    (* dec)->h_BWTPos = NULL;
  }

  if((* dec)->d_BWTPos != NULL){
    cudaFree((* dec)->d_BWTPos);
    (* dec)->d_BWTPos = NULL;
  }

  if((* dec) != NULL){
    free(* dec);
    (* dec) = NULL;
  }
  return(0);
}

uint32_t freeProfileGEM(gem_profile_t **prof)
{   
  if((* prof)->h_initBWTPos != NULL){
    free((* prof)->h_initBWTPos);
    (* prof)->h_initBWTPos = NULL;
  }

  if((* prof)->d_initBWTPos != NULL){
    cudaFree((* prof)->d_initBWTPos);
    (* prof)->d_initBWTPos = NULL;
  }

  if((* prof)->h_endBWTPos != NULL){
    free((* prof)->h_endBWTPos);
    (* prof)->h_endBWTPos = NULL;
  }

  if((* prof)->h_steps != NULL){
    free((* prof)->h_steps);
    (* prof)->h_steps = NULL;
  }

  if((* prof) != NULL){
    free(* prof);
    (* prof) = NULL;
  }
  return(0);
}

uint32_t decodingPositionsFMICPU(fmi_entry_t *h_fmi, uint64_t bwtSize, uint32_t numDecodings, uint64_t *h_initBWTPos, ulonglong2 *h_BWTPos)
{
	const uint32_t samplingRate = 4;
	const uint32_t nreps = 10;
	double start, stop;
    start = sampleTime() * 1000;

		for(uint32_t iteration = 0; iteration < nreps; ++iteration)
			decodingPositionsFMICPUKernel(h_fmi, bwtSize, numDecodings, h_initBWTPos, h_BWTPos, samplingRate);

    stop = sampleTime() * 1000;
	printf("\t Time Kernel CPU:  %8.2f ms\n", (stop - start) / nreps);

	return(0);
}


uint32_t decodingPositionsFMIGPU(fmi_entry_t *d_fmi, uint64_t bwtSize, uint32_t numDecodings, uint64_t *d_initBWTPos, ulonglong2 *d_BWTPos)
{
	const uint32_t samplingRate = 4;
	const uint32_t threads = 128;
	const uint32_t blocks  = DIV_CEIL(numDecodings * DECODE_THREADS_PER_ENTRY, threads);
	const uint32_t nreps   = 10;

	printf("GPU numDecodings=%u \n", numDecodings);

	float elapsed_time_ms = 0.0f;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

		for(uint32_t iteration = 0; iteration < nreps; ++iteration)
			decodingPositionsFMIGPUKernel<<<blocks,threads>>>(d_fmi, bwtSize, numDecodings, d_initBWTPos, d_BWTPos, samplingRate);

	cudaEventRecord(stop, 0);
	cudaThreadSynchronize();
	cudaEventElapsedTime(&elapsed_time_ms, start, stop);
	elapsed_time_ms /= nreps;

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	printf("\t Time Kernel GPU:  %8.2f ms\n", elapsed_time_ms);

	return(0);
}

int32_t transferCPUtoGPU(gem_profile_t *profile, fmi_t *fmi, dec_results_t *results)
{
	//allocate & transfer Seeds to GPU
	CUDA_ERROR(cudaMalloc((void**)&profile->d_initBWTPos, profile->numDecodings * sizeof(uint64_t)));
	CUDA_ERROR(cudaMemcpy(profile->d_initBWTPos, profile->h_initBWTPos, profile->numDecodings * sizeof(uint64_t), cudaMemcpyHostToDevice));

	//allocate & transfer FMIndex to GPU
	CUDA_ERROR(cudaMalloc((void**)&fmi->d_fmi, fmi->numEntries * sizeof(fmi_entry_t)));
	CUDA_ERROR(cudaMemcpy(fmi->d_fmi, fmi->h_fmi, fmi->numEntries * sizeof(fmi_entry_t), cudaMemcpyHostToDevice));

	//allocate & initialize Results (SA intervals)
 	CUDA_ERROR(cudaMalloc((void**)&results->d_BWTPos, results->numDecodings * sizeof(ulonglong2)));
 	CUDA_ERROR(cudaMemset(results->d_BWTPos, 0, results->numDecodings * sizeof(ulonglong2)));

	return (0);
}

int32_t transferGPUtoCPU(dec_results_t *results)
{	
	CUDA_ERROR(cudaMemcpy(results->h_BWTPos, results->d_BWTPos, results->numDecodings * sizeof(ulonglong2), cudaMemcpyDeviceToHost));
	return (0);
}

int32_t main(int argc, char *argv[])
{
  fmi_t         *fmIndex     = NULL;
  dec_results_t *results     = NULL;
  gem_profile_t *profRegions = NULL;

  char *fmiFile     = argv[1];
  char *profileFile = argv[2];

	CATCH_ERROR(initFMI(&fmIndex));
	CATCH_ERROR(initProfileGEM(&profRegions));  
	CATCH_ERROR(initResults(&results));    

	printf("=> Loading FMI ... \n");
	CATCH_ERROR(loadFMI(fmiFile, fmIndex));
	printf("\t BWT size: %llu, Number of entries: %llu\n", fmIndex->bwtSize, fmIndex->numEntries);
	
	printf("=> Loading GEM Profile ... \n");
	CATCH_ERROR(inspectGEMProfile(profileFile, profRegions));
	CATCH_ERROR(loadGEMProfile(profileFile, profRegions));

	printf("=> Initializing and Allocating Results ... \n");
	CATCH_ERROR(allocateResults(results, profRegions->numDecodings));

	printf("=> Print sample of decodings (num=%d) ... \n", profRegions->numDecodings);
	CATCH_ERROR(printGEMProfile(profRegions, 4));

	printf("=> Launching FMI Kernel (CPU) ... \n");
	CATCH_ERROR(decodingPositionsFMICPU(fmIndex->h_fmi, fmIndex->bwtSize, profRegions->numDecodings, 
									                    profRegions->h_initBWTPos, results->h_BWTPos_host));

	printf("=> Allocating and sending memory to GPU ... \n");
	CATCH_ERROR(transferCPUtoGPU(profRegions, fmIndex, results));

	printf("=> Launching FMI Kernel (GPU) ... \n");
	CATCH_ERROR(decodingPositionsFMIGPU(fmIndex->d_fmi, fmIndex->bwtSize, profRegions->numDecodings, 
	 								                    profRegions->d_initBWTPos, results->d_BWTPos));

	printf("=> Gathering results from GPU ... \n");
	CATCH_ERROR(transferGPUtoCPU(results));

	printf("=> Checking Intervals CPU ... \n");
  CATCH_ERROR(checkIntervalsCPU(profRegions, results));

	printf("=> Checking Intervals GPU ... \n");
	CATCH_ERROR(checkIntervalsGPU(profRegions, results));

	// printf("=> Saving Intervals ... \n");
	// CATCH_ERROR(saveIntervals(fmiFile, results));

	CATCH_ERROR(freeFMI(&fmIndex));
	CATCH_ERROR(freeResults(&results));
	CATCH_ERROR(freeProfileGEM(&profRegions));
	printf("=> Done! \n");

  return (0);
}







