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

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#include <cuda.h>
#include <cuda_runtime.h>

#define	BWT_CHAR_LENGTH				    3 			// 3 bits / character

#define FMI_NUM_COUNTERS			    4 			// 4 virtual counters
#define FMI_ALTERNATE_COUNTERS		2 			// 2 
#define FMI_ENTRY_SIZE				    128			// 128 bases / FMI entry
#define FMI_BITS_LOAD_PER_THREAD	128			// accesses of 16 Bytes / thread
#define FMI_ENTRIES_PER_SEED		  2

#define	UINT32_LENGTH				      32 			// 32 bits
#define	UINT64_LENGTH				      64 			// 64 bits
#define MASK_ONES					        0xFFFFFFFF
#define MASK_ONE					        0x00000001
#define MASK_ZEROS					      0x00000000

#define WARP_SIZE					        32
#define MAX_THREADS_PER_SM			  128

#define FMI_COUNTERS_PER_ENTRY		(FMI_NUM_COUNTERS / FMI_ALTERNATE_COUNTERS)									    //   2 physical counters / FMI entry
#define FMI_BITMAPS_PER_ENTRY		  (FMI_ENTRY_SIZE * BWT_CHAR_LENGTH / UINT32_LENGTH)							    //   12 physical bitmaps / FMI entry
#define FMI_ENTRY_LENGTH			    ((FMI_COUNTERS_PER_ENTRY * UINT64_LENGTH) + (FMI_ENTRY_SIZE * BWT_CHAR_LENGTH)) //   512 bits
#define FMI_THREADS_PER_ENTRY		  (FMI_ENTRY_LENGTH / FMI_BITS_LOAD_PER_THREAD)  								    //   4 threads
#define FMI_ENTRIES_PER_WARP		  (WARP_SIZE / FMI_THREADS_PER_ENTRY)
#define FMI_ENTRIES_PER_BLOCK		  (MAX_THREADS_PER_SM / FMI_THREADS_PER_ENTRY)

//#define FMI_ENTRIES_PER_BLOCK		(MAX_THREADS_PER_SM / FMI_THREADS_PER_ENTRY)

#define SEED_ENTRY_LENGTH			    128			// 128 bits
#define SEED_FIELD_SIZE				    8			  // 8 bits
#define SEED_MAX_CHARS				    60			// 60 bases
#define SEED_CHAR_LENGTH			    2 		  // 2 bits

#define SEED_THREADS_PER_ENTRY		(FMI_THREADS_PER_ENTRY * FMI_ENTRIES_PER_SEED)
#define SEED_ENTRIES_PER_WARP		  (WARP_SIZE / SEED_THREADS_PER_ENTRY)
#define SEED_BASES_PER_ENTRY		  (UINT64_LENGTH / SEED_CHAR_LENGTH)

#define V2S_B64(v,s)	asm("mov.b64 %0, {%1,%2};" : "=l"(s) : "r"(v.x), "r"(v.y))
#define S2V_B64(s,v)	asm("mov.b64 {%0,%1}, %2;" : "=r"(v.x), "=r"(v.y) : "l"(s))

#define MIN(NUM_A, NUM_B) 				        ((NUM_A < NUM_B) ? NUM_A : NUM_B)
#define CATCH_ERROR(error)                {{if (error) { fprintf(stderr, "%s\n", processError(error)); exit(EXIT_FAILURE); }}}
#define DIV_CEIL(NUMERATOR,DENOMINATOR)		(((NUMERATOR)+((DENOMINATOR)-1))/(DENOMINATOR))


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
	uint32_t    numSeeds;
	ulonglong2 *d_seeds;
	ulonglong2 *h_seeds;
	ulonglong2 *h_intervalsGEM;
	uint32_t   *h_stepsGEM;
} gem_profile_t;

typedef struct {
	uint32_t    numIntervals;
	ulonglong2 *h_intervals;
	ulonglong2 *h_intervals_host;
	ulonglong2 *d_intervals;
} intervals_t;

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

char decodeBase(uint32_t bit0, uint32_t bit1)
{
	if(bit1 == 0 && bit0 == 0) return ('A');
	if(bit1 == 0 && bit0 == 1) return ('C');
	if(bit1 == 1 && bit0 == 0) return ('G');
	if(bit1 == 1 && bit0 == 1) return ('T');	
	return ('A');
}

void backwardSearchFMICPUKernel(fmi_entry_t *fmi, uint64_t bwtSize, uint32_t numSeeds, ulonglong2 *seeds, ulonglong2 *resIntervals)
{
	const uint32_t NUM_BITMAPS = FMI_ENTRY_SIZE / UINT32_LENGTH; 
	const uint32_t LUT[12] = {3,7,11,0,1,2,4,5,6,8,9,10};

	for(uint32_t idSeed = 0; idSeed < numSeeds; ++idSeed){
		const uint32_t   seedSize = seeds[idSeed].y >> (UINT64_LENGTH - SEED_FIELD_SIZE);
		const ulonglong2 seed = seeds[idSeed];
		      uint64_t   interval_lo = 0, interval_hi = bwtSize;
		      uint64_t   currentSeed = seed.x;
		      uint32_t   idStep = 0, foundSeed = 0;

		while((idStep < seedSize) && (foundSeed == 0)){
			const uint64_t entryIdx_lo	     = interval_lo / FMI_ENTRY_SIZE;
			const uint64_t entryIdx_hi 	     = interval_hi / FMI_ENTRY_SIZE;
			const uint32_t bitmapPosition_lo = interval_lo % FMI_ENTRY_SIZE;
			const uint32_t bitmapPosition_hi = interval_hi % FMI_ENTRY_SIZE;
				  uint32_t numCharacters_lo  = 0, numCharacters_hi = 0;

			// Gathering the base of the seed
			currentSeed = (idStep == SEED_BASES_PER_ENTRY) ? seed.y : currentSeed;
			const uint32_t bit0 =  currentSeed & 0x1L;
			const uint32_t bit1 = (currentSeed & 0x2L) >> 1;
			currentSeed >>= SEED_CHAR_LENGTH;

			const uint32_t missedEntry_lo    = (entryIdx_lo % FMI_ALTERNATE_COUNTERS == bit1) ? 0 : 1;
			const uint64_t bigCounter_lo     = fmi[entryIdx_lo + missedEntry_lo].s.counters[bit0];
			// reorder bitmaps layout for low flag
			for(uint32_t idBitmap = 0; idBitmap < NUM_BITMAPS; ++idBitmap){
			 	const uint32_t initBitmap  = idBitmap * BWT_CHAR_LENGTH;
			 	uint3 vbitmap     = {fmi[entryIdx_lo].s.bitmaps[LUT[initBitmap]], fmi[entryIdx_lo].s.bitmaps[LUT[initBitmap + 1]], fmi[entryIdx_lo].s.bitmaps[LUT[initBitmap + 2]]};
			 	numCharacters_lo += computeBitmapsCPU(vbitmap, bitmapPosition_lo, bit0, bit1, missedEntry_lo, idBitmap);
			}

			const uint32_t missedEntry_hi    = (entryIdx_hi % FMI_ALTERNATE_COUNTERS == bit1) ? 0 : 1;
			const uint64_t bigCounter_hi     = fmi[entryIdx_hi + missedEntry_hi].s.counters[bit0];
			//reorder bitmaps layout for high flag
			for(uint32_t idBitmap = 0; idBitmap < NUM_BITMAPS; ++idBitmap){
				const uint32_t initBitmap  = idBitmap * BWT_CHAR_LENGTH;
			 	uint3 vbitmap     = {fmi[entryIdx_hi].s.bitmaps[LUT[initBitmap]], fmi[entryIdx_hi].s.bitmaps[LUT[initBitmap + 1]], fmi[entryIdx_hi].s.bitmaps[LUT[initBitmap + 2]]};
			 	numCharacters_hi += computeBitmapsCPU(vbitmap, bitmapPosition_hi, bit0, bit1, missedEntry_hi, idBitmap);
			}

			interval_lo = (missedEntry_lo) ? bigCounter_lo - numCharacters_lo : bigCounter_lo + numCharacters_lo;
			interval_hi = (missedEntry_hi) ? bigCounter_hi - numCharacters_hi : bigCounter_hi + numCharacters_hi;
			foundSeed = (interval_lo == interval_hi) ? 1 : 0;
			idStep++;
		}

		resIntervals[idSeed] = make_ulonglong2(interval_lo, interval_hi);
	}
}

void __global__ backwardSearchFMIGPUKernel(fmi_entry_t *fmi, uint64_t bwtSize, uint32_t numSeeds, ulonglong2 *seeds, ulonglong2 *resIntervals)
{
	const uint32_t globalThreadIdx     = blockIdx.x * MAX_THREADS_PER_SM + threadIdx.x;
	const uint32_t localWarpThreadIdx  = globalThreadIdx % WARP_SIZE;
	const uint32_t localEntryIdx       = localWarpThreadIdx  / FMI_THREADS_PER_ENTRY;
	const uint32_t localEntryThreadIdx = localWarpThreadIdx  % FMI_THREADS_PER_ENTRY;
	const uint32_t idSeed 	 		       = globalThreadIdx / SEED_THREADS_PER_ENTRY;

	if ( (threadIdx.x < MAX_THREADS_PER_SM) && (globalThreadIdx < (numSeeds * SEED_THREADS_PER_ENTRY)) ){

		const uint32_t   localIdSeed   = idSeed % SEED_ENTRIES_PER_WARP;
		const ulonglong2 seed 	       = seeds[idSeed];
		const uint32_t   seedSize      = seed.y >> (UINT64_LENGTH - SEED_FIELD_SIZE);
			    uint64_t   currentSeed   = seed.x;
			    uint64_t   interval      = (localEntryIdx % FMI_ENTRIES_PER_SEED) ? bwtSize : 0;
			    uint32_t   resultBitmaps, idStep = 0, foundSeed = 0;
			    uint64_t   sharedInterval, resultCounters;

		__shared__ exch_bmp_mem_t   exchBMP[FMI_ENTRIES_PER_BLOCK];
				   exch_bmp_mem_t * seedExchBMP = &exchBMP[threadIdx.x / FMI_THREADS_PER_ENTRY];

		while((idStep < seedSize) && (foundSeed == 0)){			
			const uint64_t entryIdx    	  =  interval / FMI_ENTRY_SIZE;
			const uint32_t bitmapPosition =  interval % FMI_ENTRY_SIZE;

			// Gathering the base of the seed
			currentSeed = (idStep == SEED_BASES_PER_ENTRY) ? seed.y : currentSeed;
			const uint32_t bit0 =  currentSeed & 0x1L;
			const uint32_t bit1 = (currentSeed & 0x2L) >> 1;
			currentSeed >>= SEED_CHAR_LENGTH;
 
			// Loading FM-index entry in thread cooperative way
			const uint32_t missedEntry   = (entryIdx % FMI_ALTERNATE_COUNTERS != bit1) ? 1 : 0;
			const uint64_t entryIdxFixed = (localEntryThreadIdx == 0) ? entryIdx + missedEntry : entryIdx;
			uint4 loadEntry              = fmi[entryIdxFixed].v[localEntryThreadIdx];

			// Entry computation
      resultCounters = selectCounter(loadEntry, bit0);
      // Reorganize entry layout between threads (send bitmaps to thread 0)
      gatherBitmaps(loadEntry, seedExchBMP, localEntryThreadIdx);
      if(localEntryThreadIdx == 0) loadEntry = seedExchBMP->s;
			resultBitmaps  = computeBitmaps(loadEntry, bitmapPosition, bit0, bit1, localEntryThreadIdx, missedEntry);

			// Shared results
			const uint32_t lane = (localEntryThreadIdx == 0) ? localWarpThreadIdx + FMI_THREADS_PER_ENTRY : localEntryIdx * FMI_THREADS_PER_ENTRY;
			interval      		  = (missedEntry) ? resultCounters - resultBitmaps : resultCounters + resultBitmaps;
			sharedInterval      = shfl_64(interval, lane);

			// Early exit condition
			foundSeed           = (__ballot(interval == sharedInterval) >> (localIdSeed * SEED_THREADS_PER_ENTRY)) & MASK_ONE;

			// Update interval for bitmap threads
			if (localEntryThreadIdx) interval = sharedInterval;

			// Increment for the next Backward-Search
			idStep++;
		}
		// Save intervals
		if((localWarpThreadIdx % SEED_THREADS_PER_ENTRY) == 0)
			resIntervals[idSeed] = make_ulonglong2(interval, sharedInterval);
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

inline void tranformSeed(unsigned long long * bitmap1, unsigned long long * bitmap0, const char * seedASCII, const uint64_t seedSize)
{
	uint64_t localBitmap[2] = {MASK_ZEROS, MASK_ZEROS};
	uint32_t index = 0;

	for(int32_t idBase = seedSize - 1; idBase >= 0; --idBase){
		const char     base = seedASCII[idBase];
		const uint32_t localIdBase = ((seedSize - 1) - idBase) % UINT32_LENGTH;
		localBitmap[index] |= charToBinASCII(base) << (localIdBase * SEED_CHAR_LENGTH);
		if(localIdBase == (UINT32_LENGTH - 1)) index++;
	}
	localBitmap[1] |= (seedSize << (UINT64_LENGTH - SEED_FIELD_SIZE));

	(* bitmap1) = localBitmap[1];
	(* bitmap0) = localBitmap[0];
}

uint32_t loadGEMProfile(const char *fn, gem_profile_t *profRegions)
{
  FILE *fp = NULL;
  const uint32_t SEED_MAX_SIZE = 20000;

  uint32_t steps;
  char 	 seedASCII[SEED_MAX_SIZE];
  uint64_t seedSize, lo, hi;

  fp = fopen(fn, "r");
  if (fp == NULL) return (8);

  //allocate the memory for the structures on profRegions
  profRegions->h_seeds	    = (ulonglong2 *) malloc(profRegions->numSeeds * sizeof(ulonglong2));
  profRegions->h_intervalsGEM = (ulonglong2 *) malloc(profRegions->numSeeds * sizeof(ulonglong2));
  profRegions->h_stepsGEM     = (uint32_t *)   malloc(profRegions->numSeeds * sizeof(uint32_t));

  for(uint32_t idSeed = 0; idSeed < profRegions->numSeeds; ++idSeed){
    fscanf(fp, "%llu %s %llu %llu %u", &seedSize, seedASCII, &lo, &hi, &steps);
    tranformSeed(&profRegions->h_seeds[idSeed].y, &profRegions->h_seeds[idSeed].x, seedASCII, seedSize);
    profRegions->h_stepsGEM[idSeed]       = steps;
    profRegions->h_intervalsGEM[idSeed].x = lo;
    profRegions->h_intervalsGEM[idSeed].y = hi;
	}

  fclose(fp);
  return (0);
}

uint32_t allocateResults(intervals_t *results, uint32_t numSeeds)
{
	results->numIntervals      = numSeeds;
  results->h_intervals       = (ulonglong2 *) malloc(results->numIntervals * sizeof(ulonglong2));
  results->h_intervals_host  = (ulonglong2 *) malloc(results->numIntervals * sizeof(ulonglong2));
  return (0);
}

uint32_t inspectGEMProfile(const char *fn, gem_profile_t *profRegions)
{
  FILE *fp = NULL;
  const uint32_t NUM_ELEMENTS_PER_LINE = 5;
  const uint32_t SEED_MAX_SIZE = 20000;
        uint32_t parsedElements = NUM_ELEMENTS_PER_LINE;

  uint32_t steps, numSeeds = 0;
  char 	   seedASCII[SEED_MAX_SIZE];
  uint64_t seedSize, lo, hi;

  fp = fopen(fn, "r");
  if (fp == NULL) return (8);

  while(parsedElements == NUM_ELEMENTS_PER_LINE){
    parsedElements = fscanf(fp, "%llu %s %llu %llu %u", &seedSize, seedASCII, &lo, &hi, &steps);
    if(parsedElements == NUM_ELEMENTS_PER_LINE) numSeeds++;
	}

	profRegions->numSeeds = numSeeds;
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

inline uint32_t printSeed(ulonglong2 seed, uint32_t seedSize)
{
	uint64_t bitmap = seed.x;
	for(uint32_t idBase = 0; idBase < seedSize; ++idBase){
		uint32_t base = bitmap & 0x3;
		printf("%c", BinASCIItoChar(base));
		bitmap >>= SEED_CHAR_LENGTH;
		if(idBase == UINT32_LENGTH) bitmap = seed.y;
	}
	return(0);
}

inline uint32_t printGEMProfile(gem_profile_t *profRegions, const uint32_t maxSeeds)
{
  for(uint32_t idSeed = 0; idSeed < MIN(profRegions->numSeeds, maxSeeds); ++idSeed){
    if(idSeed == 0) printf("seed_0.x: %lld  seed_0.y: %lld\n", profRegions->h_seeds[idSeed].x, profRegions->h_seeds[idSeed].y);
    const uint32_t seedSize = profRegions->h_seeds[idSeed].y >> (UINT64_LENGTH - SEED_FIELD_SIZE);
    printf("[%d] seed=", idSeed);
    printSeed(profRegions->h_seeds[idSeed], seedSize);
    printf("\t size=%d \t lo=%llu \t hi=%llu \t steps=%u \n", seedSize, profRegions->h_intervalsGEM[idSeed].x,
          profRegions->h_intervalsGEM[idSeed].y, profRegions->h_stepsGEM[idSeed]);
  }

  return (0);
}

uint32_t checkIntervalsGPU(gem_profile_t *profRegions, intervals_t *intervals)
{
	const uint32_t maxSeeds    = 10; // Just check and print the first results
	      uint32_t missMatches = 0;

  for(uint32_t idSeed = 0; idSeed < profRegions->numSeeds; ++idSeed){
    if((profRegions->h_intervalsGEM[idSeed].y != intervals->h_intervals[idSeed].y) ||
       (profRegions->h_intervalsGEM[idSeed].x != intervals->h_intervals[idSeed].x)){
      if(missMatches < maxSeeds){
        const uint32_t seedSize = profRegions->h_seeds[idSeed].y >> (UINT64_LENGTH - SEED_FIELD_SIZE);
        printf("[%d] seed=", idSeed);
        printSeed(profRegions->h_seeds[idSeed], seedSize);
        printf("\t size=%d \t lo=%llu \t hi=%llu \t steps=%u (GPU) lo=%llu \t hi=%llu \n",
              seedSize, profRegions->h_intervalsGEM[idSeed].x,
              profRegions->h_intervalsGEM[idSeed].y, profRegions->h_stepsGEM[idSeed],
              intervals->h_intervals[idSeed].x, intervals->h_intervals[idSeed].y);
      }
      missMatches++;
    }
	}
  return (0);
}

uint32_t checkIntervalsCPU(gem_profile_t *profRegions, intervals_t *intervals)
{
	const uint32_t maxSeeds    = 10; // Just check and print the first results
	      uint32_t missMatches = 0;

  for(uint32_t idSeed = 0; idSeed < profRegions->numSeeds; ++idSeed){
    if((profRegions->h_intervalsGEM[idSeed].y != intervals->h_intervals_host[idSeed].y) ||
       (profRegions->h_intervalsGEM[idSeed].x != intervals->h_intervals_host[idSeed].x)){
      if(missMatches < maxSeeds){
        const uint32_t seedSize = profRegions->h_seeds[idSeed].y >> (UINT64_LENGTH - SEED_FIELD_SIZE);
        printf("[%d] seed=", idSeed);
        printSeed(profRegions->h_seeds[idSeed], seedSize);
        printf("\t size=%d \t lo=%llu \t hi=%llu \t steps=%u (GPU) lo=%llu \t hi=%llu \n",
               seedSize, profRegions->h_intervalsGEM[idSeed].x,
               profRegions->h_intervalsGEM[idSeed].y, profRegions->h_stepsGEM[idSeed],
               intervals->h_intervals_host[idSeed].x, intervals->h_intervals_host[idSeed].y);
      }
      missMatches++;
    }
	}
  return (0);
}

uint32_t initProfileGEM(gem_profile_t **profRegions)
{
  gem_profile_t *prof  = (gem_profile_t *) malloc(sizeof(gem_profile_t));
	prof->numSeeds        = 0;
	prof->d_seeds         = NULL;
	prof->h_seeds         = NULL;
	prof->h_intervalsGEM  = NULL;
	prof->h_stepsGEM      = NULL;

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

uint32_t initResults(intervals_t **intervals)
{
  intervals_t *inter      = (intervals_t *) malloc(sizeof(intervals_t));
  inter->numIntervals     = 0;
  inter->h_intervals      = NULL;
  inter->h_intervals_host = NULL;
  inter->d_intervals      = NULL;

  (* intervals) = inter;
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

uint32_t freeResults(intervals_t **res)
{   
  if((* res)->h_intervals != NULL){
    free((* res)->h_intervals);
    (* res)->h_intervals = NULL;
  }

  if((* res)->d_intervals != NULL){
    cudaFree((* res)->d_intervals);
    (* res)->d_intervals = NULL;
  }

  if((* res) != NULL){
    free(* res);
    (* res) = NULL;
  }
  return(0);
}

uint32_t freeProfileGEM(gem_profile_t **prof)
{   
  if((* prof)->h_seeds != NULL){
    free((* prof)->h_seeds);
    (* prof)->h_seeds = NULL;
  }

  if((* prof)->d_seeds != NULL){
    cudaFree((* prof)->d_seeds);
    (* prof)->d_seeds = NULL;
  }

  if((* prof)->h_intervalsGEM != NULL){
    free((* prof)->h_intervalsGEM);
    (* prof)->h_intervalsGEM = NULL;
  }

  if((* prof)->h_stepsGEM != NULL){
    free((* prof)->h_stepsGEM);
    (* prof)->h_stepsGEM = NULL;
  }

  if((* prof) != NULL){
    free(* prof);
    (* prof) = NULL;
  }
  return(0);
}

uint32_t backwardSearchFMICPU(fmi_entry_t *h_fmi, uint64_t bwtSize, uint32_t numSeeds, ulonglong2 *h_seeds, ulonglong2 *h_intervals_host)
{
  const uint32_t nreps = 10;
  double start, stop;
  start = sampleTime() * 1000;

  for(uint32_t iteration = 0; iteration < nreps; ++iteration)
    backwardSearchFMICPUKernel(h_fmi, bwtSize, numSeeds, h_seeds, h_intervals_host);

  stop = sampleTime() * 1000;
  printf("\t Time Kernel CPU:  %8.2f ms\n", (stop - start) / nreps);
	return(0);
}


uint32_t backwardSearchFMIGPU(fmi_entry_t *d_fmi, uint64_t bwtSize, uint32_t numSeeds, ulonglong2 *d_seeds, ulonglong2 *d_intervals)
{
	const uint32_t threads = 128;
	const uint32_t blocks  = DIV_CEIL(numSeeds * SEED_THREADS_PER_ENTRY, threads);
	const uint32_t nreps   = 10;

	float elapsed_time_ms = 0.0f;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

		for(uint32_t iteration = 0; iteration < nreps; ++iteration)
			backwardSearchFMIGPUKernel<<<blocks,threads>>>(d_fmi, bwtSize, numSeeds, d_seeds, d_intervals);

	cudaEventRecord(stop, 0);
	cudaThreadSynchronize();
	cudaEventElapsedTime(&elapsed_time_ms, start, stop);
	elapsed_time_ms /= nreps;

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	printf("\t Time Kernel GPU:  %8.2f ms\n", elapsed_time_ms);

	return(0);
}

int32_t transferCPUtoGPU(gem_profile_t *profile, fmi_t *fmi, intervals_t *results)
{
	//allocate & transfer Seeds to GPU
	CUDA_ERROR(cudaMalloc((void**)&profile->d_seeds, profile->numSeeds * sizeof(ulonglong2)));
	CUDA_ERROR(cudaMemcpy(profile->d_seeds, profile->h_seeds, profile->numSeeds * sizeof(ulonglong2), cudaMemcpyHostToDevice));

	//allocate & transfer FMIndex to GPU
	CUDA_ERROR(cudaMalloc((void**)&fmi->d_fmi, fmi->numEntries * sizeof(fmi_entry_t)));
	CUDA_ERROR(cudaMemcpy(fmi->d_fmi, fmi->h_fmi, fmi->numEntries * sizeof(fmi_entry_t), cudaMemcpyHostToDevice));

	//allocate & initialize Results (SA intervals)
 	CUDA_ERROR(cudaMalloc((void**)&results->d_intervals, results->numIntervals * sizeof(ulonglong2)));
 	CUDA_ERROR(cudaMemset(results->d_intervals, 0, results->numIntervals * sizeof(ulonglong2)));

	return (0);
}

int32_t transferGPUtoCPU(intervals_t *results)
{	
	CUDA_ERROR(cudaMemcpy(results->h_intervals, results->d_intervals, results->numIntervals * sizeof(ulonglong2), cudaMemcpyDeviceToHost));
	return (0);
}

int32_t main(int argc, char *argv[])
{
  fmi_t         *fmIndex     = NULL;
  intervals_t   *results     = NULL;
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
	CATCH_ERROR(allocateResults(results, profRegions->numSeeds));

	printf("=> Print sample of seeds (num=%d) ... \n", profRegions->numSeeds);
	CATCH_ERROR(printGEMProfile(profRegions, 4));

	printf("=> Launching FMI Kernel (CPU) ... \n");
	CATCH_ERROR(backwardSearchFMICPU(fmIndex->h_fmi, fmIndex->bwtSize, profRegions->numSeeds, 
									                 profRegions->h_seeds, results->h_intervals_host));

	printf("=> Allocating and sending memory to GPU ... \n");
	CATCH_ERROR(transferCPUtoGPU(profRegions, fmIndex, results));

	printf("=> Launching FMI Kernel (GPU) ... \n");
	CATCH_ERROR(backwardSearchFMIGPU(fmIndex->d_fmi, fmIndex->bwtSize, profRegions->numSeeds, 
	 						                     profRegions->d_seeds, results->d_intervals));

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







