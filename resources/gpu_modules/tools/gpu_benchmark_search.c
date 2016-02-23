/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#include "../gpu_interface.h"
#include "gpu_benchmark_commons.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>


typedef struct {
	/* Example test fields */
	uint32_t 					        numSeeds;
	uint32_t 					        numIntervals;
	gpu_fmi_search_seed_t		  *seeds;
	gpu_fmi_search_sa_inter_t	*intervals;
	gpu_fmi_search_sa_inter_t	*intervalsGEM;
} test_t;

inline void tranformSeed(uint64_t * bitmap1, uint64_t * bitmap0, const char * seedASCII, const uint64_t seedSize)
{
	uint64_t localBitmap[2] = {MASK_ZEROS, MASK_ZEROS};
	uint32_t index = 0;
	int32_t idBase;

	for(idBase = seedSize - 1; idBase >= 0; --idBase){
		const char     base = seedASCII[idBase];
		const uint32_t localIdBase = ((seedSize - 1) - idBase) % UINT32_LENGTH;
		localBitmap[index] |= charToBinASCII(base) << (localIdBase * SEED_CHAR_LENGTH);
		if(localIdBase == (UINT32_LENGTH - 1)) index++;
	}
	localBitmap[1] |= (seedSize << (UINT64_LENGTH - SEED_FIELD_SIZE));

	(* bitmap1) = localBitmap[1];
	(* bitmap0) = localBitmap[0];
}

uint32_t loadGEMProfile(const char *fn, test_t *profRegions)
{
  FILE *fp = NULL;
  const uint32_t SEED_MAX_SIZE = 20000;

  uint32_t steps, idSeed;
  char 	 seedASCII[SEED_MAX_SIZE];
  uint64_t seedSize, lo, hi;

  fp = fopen(fn, "r");
  if (fp == NULL) return (8);

  for(idSeed = 0; idSeed < profRegions->numSeeds; ++idSeed){
    fscanf(fp, "%llu %s %llu %llu %u", &seedSize, seedASCII, &lo, &hi, &steps);
    tranformSeed(&profRegions->seeds[idSeed].low, &profRegions->seeds[idSeed].hi, seedASCII, seedSize);
    profRegions->intervalsGEM[idSeed].low = lo;
    profRegions->intervalsGEM[idSeed].hi  = hi;
	}

  fclose(fp);
  return (0);
}

uint32_t inspectGEMProfile(const char *fn, uint32_t *totalNumSeeds)
{
  FILE *fp = NULL;
  const uint32_t NUM_ELEMENTS_PER_LINE = 5;
  const uint32_t SEED_MAX_SIZE = 20000;
        uint32_t parsedElements = NUM_ELEMENTS_PER_LINE;

  uint32_t steps, numSeeds = 0;
  char 	 seedASCII[SEED_MAX_SIZE];
  uint64_t seedSize, lo, hi;

  fp = fopen(fn, "r");
  if (fp == NULL) return (8);

  while(parsedElements == NUM_ELEMENTS_PER_LINE){
    parsedElements = fscanf(fp, "%llu %s %llu %llu %u", &seedSize, seedASCII, &lo, &hi, &steps);
    if(parsedElements == NUM_ELEMENTS_PER_LINE) numSeeds++;
  }

  (* totalNumSeeds) = numSeeds;

  fclose(fp);
  return (0);
}

uint32_t loadTestData(char *seedsFile, test_t *testData)
{
	int32_t error, totalNumSeeds;

	error = inspectGEMProfile(seedsFile, &totalNumSeeds);
	if(error != 0){fprintf(stderr, "Error %d, inspecting seeds \n", error); exit(EXIT_FAILURE);}

	testData->numSeeds     = totalNumSeeds;
	testData->numIntervals = totalNumSeeds;

	//Reserve memory
	testData->seeds = (gpu_fmi_search_seed_t *) malloc(testData->numSeeds * sizeof(gpu_fmi_search_seed_t));
	if (testData->seeds == NULL){fprintf(stderr, "Error allocating seeds \n"); exit(EXIT_FAILURE);}
	testData->intervals = (gpu_fmi_search_sa_inter_t *) malloc(testData->numIntervals * sizeof(gpu_fmi_search_sa_inter_t));
	if (testData->intervals == NULL){fprintf(stderr, "Error allocating results \n"); exit(EXIT_FAILURE);}
	testData->intervalsGEM = (gpu_fmi_search_sa_inter_t *) malloc(testData->numIntervals * sizeof(gpu_fmi_search_sa_inter_t));
	if (testData->intervalsGEM == NULL){fprintf(stderr, "Error allocating GEM results \n"); exit(EXIT_FAILURE);}

	error = loadGEMProfile(seedsFile, testData);
	if(error != 0){fprintf(stderr, "Error %d, loading seeds \n", error); exit(EXIT_FAILURE);}

	return (0);
}

uint32_t saveResults(const char *fn, test_t *testData, uint32_t numBuffers)
{
	FILE *fp = NULL;
	char resFileOut[512];
	char cadena[512];
	uint32_t error, idTotalInterval = 0, idBuffer, idInterval;

	for(idBuffer = 0; idBuffer < numBuffers; idBuffer++){
		sprintf(resFileOut, "%s.buffer%d.res", fn, idBuffer);

		fp = fopen(resFileOut, "wb");
		if (fp == NULL) {fprintf(stderr, "Error opening file \n"); exit(EXIT_FAILURE);}

		sprintf(cadena, "%u\n", testData[idBuffer].numIntervals);
		fputs(cadena, fp);

		idTotalInterval = 0;
		for(idInterval = 0; idInterval < testData[idBuffer].numIntervals; ++idInterval){
			sprintf(cadena, "%u %llu %llu\n", idTotalInterval, testData[idBuffer].intervals[idInterval].low, testData[idBuffer].intervals[idInterval].hi);
			fputs(cadena, fp);
			idTotalInterval++;
		}
		fclose(fp);
	}

	return (0);
}

uint32_t putIntoBuffer(void *buffer, test_t *testData)
{
	uint32_t i;
	gpu_fmi_search_seed_t *seed_buffer = gpu_fmi_search_buffer_get_seeds_(buffer);

	testData->numSeeds 		= MIN(testData->numSeeds, gpu_fmi_search_buffer_get_max_seeds_(buffer));
	testData->numIntervals 	= testData->numSeeds;

	for(i = 0; i < testData->numSeeds; ++i){
		seed_buffer[i].hi = testData->seeds[i].hi;
		seed_buffer[i].low = testData->seeds[i].low;
	}

	return (0);
}

uint32_t getFromBuffer(void *buffer, test_t *testData)
{
	uint32_t i;
	gpu_fmi_search_sa_inter_t *intervals_buffer = gpu_fmi_search_buffer_get_sa_intervals_(buffer);

	for(i = 0; i < testData->numIntervals; i++){
		testData->intervals[i].hi = intervals_buffer[i].hi;
		testData->intervals[i].low = intervals_buffer[i].low;
	}

	return (0);
}


uint32_t freeTestData(test_t *testData, uint32_t numThreads)
{
	uint32_t threadID;
	for(threadID = 0; threadID < numThreads; ++threadID){
		free(testData[threadID].seeds);
		free(testData[threadID].intervals);
		free(testData[threadID].intervalsGEM);
	}
	return (0);
}

inline uint32_t printSeed(gpu_fmi_search_seed_t seed, uint32_t seedSize)
{
	uint64_t bitmap = seed.low;
	uint32_t idBase;

	for(idBase = 0; idBase < seedSize; ++idBase){
		uint32_t base = bitmap & 0x3;
		printf("%c", BinASCIItoChar(base));
		bitmap >>= SEED_CHAR_LENGTH;
		if(idBase == UINT32_LENGTH) bitmap = seed.hi;
	}
	return(0);
}

uint32_t checkIntervalsGPU(test_t *profRegions)
{
	const uint32_t maxSeeds    = 10; // Just check and print the first results
	      uint32_t missMatches = 0;
	      uint32_t idSeed;

  for(idSeed = 0; idSeed < profRegions->numSeeds; ++idSeed){
    if((profRegions->intervalsGEM[idSeed].hi != profRegions->intervals[idSeed].hi) ||
       (profRegions->intervalsGEM[idSeed].low != profRegions->intervals[idSeed].low)){
      if(missMatches < maxSeeds){
        const uint32_t seedSize = profRegions->seeds[idSeed].hi >> (UINT64_LENGTH - SEED_FIELD_SIZE);
      printf("[%d] seed=", idSeed);
      printSeed(profRegions->seeds[idSeed], seedSize);
      printf("\t size=%d \t (CPU) lo=%llu \t hi=%llu \t (GPU) lo=%llu \t hi=%llu \n",
          seedSize, profRegions->intervalsGEM[idSeed].low, profRegions->intervalsGEM[idSeed].hi,
          profRegions->intervals[idSeed].low, profRegions->intervals[idSeed].hi);
    }
    missMatches++;
    }
	}

  return (0);
}

double processSearchFMI(char *fmiFile, char *seedsFile, uint32_t numBuffers, uint32_t numThreads, float maxMbPerBuffer, uint32_t numTasks)
{
	test_t				testData[numThreads];
	uint32_t 			error, idBuffer, iteration, threadID, idTask;
	double				ts, ts1;

	gpu_buffers_dto_t 	buff  = {.buffer 		 		        = NULL,
								               .numBuffers 	 		      = numBuffers,
								               .maxMbPerBuffer 		    = maxMbPerBuffer,
								               .activeModules  		    = GPU_FMI_EXACT_SEARCH};
	gpu_index_dto_t 	  index = {.fmi 			 		        = fmiFile,
								               .indexCoding 	 		    = GPU_INDEX_PROFILE_FILE,
								               .bwtSize 		 		      = 0};
	gpu_reference_dto_t ref   = {.reference 	 		      = NULL,
								               .refCoding 			      = GPU_REF_NONE,
								               .refSize				        = 0};
  gpu_info_dto_t      sys   = {.selectedArchitectures = GPU_ARCH_SUPPORTED,
                               .userAllocOption       = GPU_LOCAL_OR_REMOTE_DATA,
                               .activatedModules      = GPU_NONE_MODULES,
                               .allocatedStructures   = GPU_NONE_MODULES
                               .verbose               = verbose};

	for(threadID = 0; threadID < numThreads; ++threadID)
		loadTestData(seedsFile, &testData[threadID]);

	// Initialize the systems and buffers
	gpu_init_buffers_(&buff, &index, &ref, &sys);

	// Master thread initialize all the buffers
	// Better each thread initialize self buffers
	for(idBuffer = 0; idBuffer < numBuffers; ++idBuffer)
		gpu_alloc_buffer_(buff.buffer[idBuffer]);

	ts = sample_time();

  #pragma omp parallel for num_threads(numThreads) schedule(static,1) private(idTask, idBuffer)
  for(idTask = 0; idTask < numTasks; ++idTask){
    idBuffer = idTask % numBuffers;

    //Trace the jobs
    printf("Host thread %d \t sent job %d \t to buffer %d \t in device %d \n", omp_get_thread_num(), idTask, idBuffer, gpu_buffer_get_id_device_(buff.buffer[idBuffer]));

    //Fill the buffer (generate work)
    gpu_fmi_search_init_buffer_(buff.buffer[idBuffer]);
    putIntoBuffer(buff.buffer[idBuffer], &testData[idBuffer]);
    gpu_fmi_search_send_buffer_(buff.buffer[idBuffer],testData[idBuffer].numSeeds);
    gpu_fmi_search_receive_buffer_(buff.buffer[idBuffer]);

    //Get the results from the buffer (consume results)
    getFromBuffer(buff.buffer[idBuffer], &testData[idBuffer]);
  }

	ts1 = sample_time();
	gpu_destroy_buffers_(&buff);

	for(threadID = 0; threadID < numThreads; ++threadID)
		checkIntervalsGPU(&testData[threadID]);

	//Save last state of each buffer
	error = saveResults(seedsFile, &testData[0], numBuffers);
	if(error != 0){fprintf(stderr, "Error %d, saving results \n", error); exit(EXIT_FAILURE);}

	freeTestData(&testData[0], numThreads);
	return(ts1-ts);
}


uint32_t main(int argc, char *argv[])
{
	void 		 	**buffer;
	uint32_t		numBuffers, numThreads, idBuffer, maxMbPerBuffer, numTasks, error;
	char 			*fmiFile, *seedsFile;
	double 			timeElapsed;

	fmiFile    = argv[2];
	seedsFile  = argv[3];

	/* FMI Backward-Search GPU */
	/* Execution parameters:
	 * 1) 1 
	 * 2) Reference file path:	example: ~/Myers/data/human_g1k_v37.cleaned.fasta.3101804739.4bits.ref
	 * 3) Queries paths: 		example: ~/Myers/data/1000.qry
	 * 4) Number of buffers		example: 8
	 * 5) Maximum buffer size	example: 32  (in MBytes)
	 * 6) Number of tasks		example: 120
	 */
	if(atoi(argv[1]) == 0){
		numBuffers 		  = atoi(argv[4]);
		numThreads 		  = atoi(argv[4]);
		maxMbPerBuffer 	= atoi(argv[5]);
		numTasks 		    = atoi(argv[6]);
	 	timeElapsed 	  = processSearchFMI(fmiFile, seedsFile, numBuffers, numThreads, (float) maxMbPerBuffer, numTasks);
	}

	/* Myers-GPU */
	if(atoi(argv[1]) > 0){
		printf("Option for FMI backward search not implemented: exiting \n");
		exit(EXIT_FAILURE);		
	}	

	printf("TOTAL TIME: \t %f \n", timeElapsed);

    return (0);
}
