/*
 * PROJECT: Bit-Parallel Myers on GPU
 * FILE: myers-interface.h
 * DATE: 4/7/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Code example of using the BMP GPU interface 
 */

#include "../gpu_interface.h"
#include "gpu_benchmark_commons.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

typedef struct {
	/* Example test fields */
	uint32_t 					numInitDecodings;
	uint32_t 					numEndDecodings;
	uint32_t					samplingRate;
	gpu_fmi_decode_init_pos_t	*initDecodings;
	gpu_fmi_decode_end_pos_t	*endDecodings;
	gpu_fmi_decode_end_pos_t	*endDecodingsGEM;
} test_t;


uint32_t loadGEMProfile(const char *fn, test_t *profRegions)
{
    FILE *fp = NULL;

    uint32_t idDecoding, steps;
    uint64_t initBWTPos, endBWTPos;
    uint32_t bookmark, numAs, numCs, numGs, numTs, numXs;

    fp = fopen(fn, "r");
    if (fp == NULL) return (8);

    for(idDecoding = 0; idDecoding < profRegions->numInitDecodings; ++idDecoding){
		fscanf(fp, "%llu %llu %u %c %u %u %u %u %u",
				&initBWTPos, &endBWTPos, &steps, &bookmark,
				&numAs, &numCs, &numGs, &numTs, &numXs);
		profRegions->initDecodings[idDecoding]             = initBWTPos;
		profRegions->endDecodingsGEM[idDecoding].interval  = endBWTPos;
		profRegions->endDecodingsGEM[idDecoding].steps     = steps;
	}

    fclose(fp);
    return (0);
}

uint32_t inspectGEMProfile(const char *fn, uint32_t *numTotalDecodings)
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

    (* numTotalDecodings) = numDecodings;

    fclose(fp);
    return (0);
}

uint32_t loadTestData(char *decodeFile, test_t *testData)
{
	int32_t error, totalNumDecodings;

	error = inspectGEMProfile(decodeFile, &totalNumDecodings);
	if(error != 0){fprintf(stderr, "Error %d, inspecting seeds \n", error); exit(EXIT_FAILURE);}

	testData->numInitDecodings = totalNumDecodings;
	testData->numEndDecodings = totalNumDecodings;


    //allocate the memory for the structures on testData
	testData->initDecodings    = (gpu_fmi_decode_init_pos_t *) malloc(testData->numInitDecodings * sizeof(gpu_fmi_decode_init_pos_t));
	if (testData->initDecodings == NULL){fprintf(stderr, "Error allocating initialized decodings \n"); exit(EXIT_FAILURE);}
	testData->endDecodings     = (gpu_fmi_decode_end_pos_t *)  malloc(testData->numEndDecodings  * sizeof(gpu_fmi_decode_end_pos_t));
	if (testData->endDecodings == NULL){fprintf(stderr, "Error allocating end decodings \n"); exit(EXIT_FAILURE);}
	testData->endDecodingsGEM  = (gpu_fmi_decode_end_pos_t *)  malloc(testData->numEndDecodings  * sizeof(gpu_fmi_decode_end_pos_t));
	if (testData->endDecodingsGEM == NULL){fprintf(stderr, "Error allocating end GEM decodings \n"); exit(EXIT_FAILURE);}

	error = loadGEMProfile(decodeFile, testData);
	if(error != 0){fprintf(stderr, "Error %d, loading decodings input \n", error); exit(EXIT_FAILURE);}

	return (0);
}

uint32_t saveResults(const char *fn, test_t *testData, uint32_t numBuffers)
{
	FILE *fp = NULL;
	char resFileOut[512];
	char cadena[512];
	uint32_t error, idTotalDecodings = 0, idBuffer, idDecoding;

	for(idBuffer = 0; idBuffer < numBuffers; idBuffer++){
		sprintf(resFileOut, "%s.buffer%d.res", fn, idBuffer);

		fp = fopen(resFileOut, "wb");
		if (fp == NULL) {fprintf(stderr, "Error opening file \n"); exit(EXIT_FAILURE);}

		sprintf(cadena, "%u\n", testData[idBuffer].numInitDecodings);
		fputs(cadena, fp);

		idTotalDecodings = 0;
		for(idDecoding = 0; idDecoding < testData[idBuffer].numInitDecodings; ++idDecoding){
			sprintf(cadena, "%u %llu %llu %llu\n", idTotalDecodings, testData[idBuffer].initDecodings[idDecoding], testData[idBuffer].endDecodings[idDecoding].interval, testData[idBuffer].endDecodings[idDecoding].steps);
			fputs(cadena, fp);
			idTotalDecodings++;
		}
		fclose(fp);
	}

	return (0);
}

uint32_t putIntoBuffer(void *buffer, test_t *testData)
{
	uint32_t i;
	gpu_fmi_decode_init_pos_t *initDec_buffer = gpu_fmi_decode_buffer_get_init_pos_(buffer);

	testData->numInitDecodings 	= MIN(testData->numInitDecodings, gpu_fmi_decode_buffer_get_max_positions_(buffer));
	testData->numEndDecodings 	= testData->numInitDecodings;

	for(i = 0; i < testData->numInitDecodings; ++i){
		initDec_buffer[i]  = testData->initDecodings[i];
	}

	return (0);
}

uint32_t getFromBuffer(void *buffer, test_t *testData)
{
	uint32_t i;
	gpu_fmi_decode_end_pos_t *endDec_buffer = gpu_fmi_decode_buffer_get_end_pos_(buffer);

	for(i = 0; i < testData->numEndDecodings; i++){
		testData->endDecodings[i].interval = endDec_buffer[i].interval;
		testData->endDecodings[i].steps = endDec_buffer[i].steps;
	}

	return (0);
}


uint32_t freeTestData(test_t *testData, uint32_t numThreads)
{
	uint32_t threadID;
	for(threadID = 0; threadID < numThreads; ++threadID){
		free(testData[threadID].initDecodings);
		free(testData[threadID].endDecodings);
		free(testData[threadID].endDecodingsGEM);
	}
	return (0);
}

uint32_t checkDecodingsGPU(test_t *profRegions)
{
	const uint32_t maxDecodings  = 10; // Just check and print the first results
	      uint32_t missMatches 	 = 0;
	      uint32_t idDecoding;

    for(idDecoding = 0; idDecoding < profRegions->numEndDecodings; ++idDecoding){
    	if((profRegions->endDecodingsGEM[idDecoding].interval != profRegions->endDecodings[idDecoding].interval) ||
    	   (profRegions->endDecodingsGEM[idDecoding].steps 	  != profRegions->endDecodings[idDecoding].steps)){
    		if(missMatches < maxDecodings){
				printf("[%d] Decoding=", idDecoding);
				printf("\t (CPU) pos=%llu \t steps=%llu \t (GPU) posv=%llu \t steps=%llu \n",
						profRegions->endDecodingsGEM[idDecoding].interval, profRegions->endDecodingsGEM[idDecoding].steps,
						profRegions->endDecodings[idDecoding].interval, profRegions->endDecodings[idDecoding].steps);
			}
			missMatches++;
    	}
	}

    return (0);
}

double processDecodeFMI(char *fmiFile, char *decodeFile, uint32_t numBuffers, uint32_t numThreads, float maxMbPerBuffer, uint32_t numTasks)
{
	void 	 	**buffer;
	test_t		testData[numThreads];

	const uint32_t SAMPLING_RATE = 4;
	uint32_t 	error, idBuffer, iteration, threadID, idTask;
	double		ts, ts1;

	for(threadID = 0; threadID < numThreads; ++threadID){
		testData[threadID].samplingRate = SAMPLING_RATE;
		loadTestData(decodeFile, &testData[threadID]);
	}

	gpu_init_buffers_(&buffer, numBuffers, maxMbPerBuffer,
			 	 	  NULL, GPU_NONE_DATA, 0,
					  fmiFile, GPU_REF_PROFILE_FILE, 0,
					  GPU_FMI_DECODE_POS, GPU_ARCH_SUPPORTED, GPU_LOCAL_DATA, 0);

	// Master thread initialize all the buffers
	// Better each thread initialize self buffers
	for(idBuffer = 0; idBuffer < numBuffers; ++idBuffer)
		gpu_alloc_buffer_(buffer[idBuffer]);

	ts = sample_time();

		#pragma omp parallel for num_threads(numThreads) schedule(static,1) private(idTask, idBuffer)
		for(idTask = 0; idTask < numTasks; ++idTask){
				idBuffer = idTask % numBuffers;

				//Trace the jobs
				printf("Host thread %d \t sent job %d \t to buffer %d \t in device %d \n", omp_get_thread_num(), idTask, idBuffer, gpu_buffer_get_id_device_(buffer[idBuffer]));

				//Fill the buffer (generate work)
				gpu_fmi_decode_init_buffer_(buffer[idBuffer]);
				putIntoBuffer(buffer[idBuffer], &testData[idBuffer]);
				gpu_fmi_decode_send_buffer_(buffer[idBuffer],testData[idBuffer].numInitDecodings, testData[idBuffer].samplingRate);
				gpu_fmi_decode_receive_buffer_(buffer[idBuffer]);

				//Get the results from the buffer (consume results)
				getFromBuffer(buffer[idBuffer], &testData[idBuffer]);
		}

	ts1 = sample_time();
	gpu_destroy_buffers_(&buffer);

	for(threadID = 0; threadID < numThreads; ++threadID)
		checkDecodingsGPU(&testData[threadID]);

	//Save last state of each buffer
	error = saveResults(decodeFile, &testData[0], numBuffers);
		if(error != 0){fprintf(stderr, "Error %d, saving results \n", error); exit(EXIT_FAILURE);}

	freeTestData(&testData[0], numThreads);

	return(ts1-ts);
}


uint32_t main(int argc, char *argv[])
{
	void 		 	**buffer;
	uint32_t		numBuffers, numThreads, idBuffer, maxMbPerBuffer, numTasks, error;
	char 			*fmiFile, *decodeFile;
	double 			timeElapsed;

	fmiFile     = argv[2];
	decodeFile  = argv[3];

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
		numBuffers 		= atoi(argv[4]);
		numThreads 		= atoi(argv[4]);
		maxMbPerBuffer 	= atoi(argv[5]);
		numTasks 		= atoi(argv[6]);

	 	timeElapsed 	= processDecodeFMI(fmiFile, decodeFile, numBuffers, numThreads, (float) maxMbPerBuffer, numTasks);
	}

	/* Myers-GPU */
	if(atoi(argv[1]) > 0){
		printf("Option for FMI decode primitive not implemented: exiting \n");
		exit(EXIT_FAILURE);		
	}	

	printf("TOTAL TIME: \t %f \n", timeElapsed);

    return (0);
}
