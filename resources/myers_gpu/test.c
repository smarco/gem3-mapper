#include "myers-interface.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#ifndef MIN
	#define MIN(_a, _b) (((_a) < (_b)) ? (_a) : (_b))
#endif

typedef struct {
	uint32_t totalQueriesEntries;
	uint32_t numQueries;
	uint32_t numCandidates;
	uint32_t numResults;
	qryEntry_t	*queries;
	candInfo_t	*candidates;
	qryInfo_t	*qinfo;
	resEntry_t	*results;
} test_t;

double sampleTime()
{
	struct timespec tv;
	clock_gettime(CLOCK_REALTIME, &tv);
	return((tv.tv_sec+tv.tv_nsec/1000000000.0));
}

uint32_t loadQueries(const char *fn, test_t *testData, uint32_t *averageQuerySize, uint32_t * averageCandidatesPerQuery)
{
	FILE *fp = NULL;
	size_t result;
	uint32_t idQuery, dummy;

	fp = fopen(fn, "rb");
	if (fp == NULL) return (10);

	fread(&dummy, sizeof(uint32_t), 1, fp);
	fread(&testData->totalQueriesEntries, sizeof(uint32_t), 1, fp);
	fread(&dummy, sizeof(uint32_t), 1, fp);
	fread(&testData->numQueries, sizeof(uint32_t), 1, fp);
	fread(&testData->numCandidates, sizeof(uint32_t), 1, fp);
	fread(averageCandidatesPerQuery, sizeof(uint32_t), 1, fp);
	fread(averageQuerySize, sizeof(uint32_t), 1, fp);

	testData->queries = NULL;
	testData->candidates = NULL;
	testData->qinfo = NULL;

	testData->queries = (qryEntry_t *) malloc(testData->totalQueriesEntries * sizeof(qryEntry_t));
		if (testData->queries == NULL) return (11);
	result = fread(testData->queries, sizeof(qryEntry_t), testData->totalQueriesEntries, fp);
		if (result != testData->totalQueriesEntries) return (12);

	testData->candidates = (candInfo_t *) malloc(testData->numCandidates * sizeof(candInfo_t));
		if (testData->candidates == NULL) return (13);
	result = fread(testData->candidates , sizeof(candInfo_t), testData->numCandidates, fp);
		if (result != testData->numCandidates) return (14);

	testData->qinfo = (qryInfo_t *) malloc(testData->numQueries * sizeof(qryInfo_t));
		if (testData->qinfo == NULL) return (15);
	result = fread(testData->qinfo , sizeof(qryInfo_t), testData->numQueries, fp);
		if (result != testData->numQueries) return (16);

	fclose(fp);
	return (0);
}

uint32_t saveResults(const char *fn, test_t *testData, uint32_t numBuffers)
{
	FILE *fp = NULL;
	char resFileOut[512];
	char cadena[512];
	uint32_t error, idTotalCandidate = 0, idBuffer, idCandidate;
	uint32_t totalResult = 0;

	for(idBuffer = 0; idBuffer < numBuffers; idBuffer++){
		sprintf(resFileOut, "%s.buffer%d.res", fn, idBuffer);

		fp = fopen(resFileOut, "wb");
			if (fp == NULL) return (47);

		sprintf(cadena, "%u\n", testData[idBuffer].numResults);
		fputs(cadena, fp);

		idTotalCandidate = 0;
		for(idCandidate = 0; idCandidate < testData[idBuffer].numResults; idCandidate++){
			sprintf(cadena, "%u %u %u\n", idTotalCandidate, testData[idBuffer].results[idCandidate].score, testData[idBuffer].results[idCandidate].column);
			fputs(cadena, fp);
			idTotalCandidate++;
		}
		fclose(fp);
	}

	return (0);
}

uint32_t putIntoBuffer(void *buffer, test_t *testData)
{
	uint32_t i;

	qryEntry_t * 	queries_buffer 		= getPEQBuffer(buffer);
	candInfo_t *	candidates_buffer 	= getCandidatesBuffer(buffer);
	qryInfo_t *		queryInfo_buffer	= getPEQInfoBuffer(buffer);

	testData->totalQueriesEntries = MIN(testData->totalQueriesEntries, getMaxPEQEntries(buffer));
	testData->numCandidates 	  = MIN(testData->numCandidates, getMaxCandidates(buffer));
	testData->numQueries 		  = MIN(testData->numQueries, getMaxQueries(buffer));
	testData->numResults 		  = testData->numCandidates;

	for(i = 0; i < testData->totalQueriesEntries; i++)
		queries_buffer[i] = testData->queries[i];
	for(i = 0; i < testData->numQueries; i++)
		queryInfo_buffer[i] = testData->qinfo[i];
	for(i = 0; i < testData->numCandidates; i++)
		candidates_buffer[i] = testData->candidates[i];

	return (0);
}

uint32_t getFromBuffer(void *buffer, test_t *testData)
{
	uint32_t i;
	resEntry_t * results_buffer = getResultsBuffer(buffer);

	for(i = 0; i < testData->numResults; i++)
		testData->results[i] = results_buffer[i];

	return (0);
}

uint32_t loadTestData(char *qryFile, test_t *testData, uint32_t *averageQuerySize, uint32_t *averageCandidatesPerQuery)
{
	int32_t error;

	error = loadQueries(qryFile, testData, averageQuerySize, averageCandidatesPerQuery);
		if(error != 0){fprintf(stderr, "Error %d, loading queries \n", error); exit(EXIT_FAILURE);}
	testData->results = (resEntry_t *) malloc(testData->numCandidates * sizeof(resEntry_t));
		if (testData->queries == NULL){fprintf(stderr, "Error allocating results \n"); exit(EXIT_FAILURE);}
	return (0);
}

uint32_t freeTestData(test_t *testData, uint32_t numThreads)
{
	uint32_t threadID;
	for(threadID = 0; threadID < numThreads; ++threadID){
		free(testData[threadID].queries);
		free(testData[threadID].candidates);
		free(testData[threadID].qinfo);
		free(testData[threadID].results);
	}
	return (0);
}

/* Execution parameters:
 * 1) Reference file path:	example: ~/Myers/data/human_g1k_v37.cleaned.fasta.3101804739.4bits.ref
 * 2) Queries paths: 		example: ~/Myers/data/1000.qry
 * 3) Number of buffers		example: 8
 * 4) Maximum buffer size	example: 32  (in MBytes)
 * 5) Number of tasks		example: 120
 */

uint32_t main(int argc, char *argv[])
{
	void 		 	**buffer;
	uint32_t		numBuffers = atoi(argv[3]), numThreads = atoi(argv[3]), error;
	uint32_t 		idBuffer, iteration, threadID, maxMbPerBuffer = atoi(argv[4]), idTask, numTasks = atoi(argv[5]);
	uint32_t		averageQuerySize, averageCandidatesPerQuery;
	char 			*refFile  = argv[1];
	double			ts, ts1;

	test_t			testData[numThreads]; //simulate the gen. work threads
	char 			*qryFile = argv[2];


	for(threadID = 0; threadID < numThreads; ++threadID)
		loadTestData(qryFile, &testData[threadID], &averageQuerySize, &averageCandidatesPerQuery);

	initMyers(&buffer, numBuffers, maxMbPerBuffer, refFile, PROFILE_REFERENCE_FILE, 0, averageQuerySize,
				averageCandidatesPerQuery, ARCH_SUPPORTED);

	ts = sampleTime();

		#pragma omp parallel for num_threads(numThreads) schedule(static,1) private(idTask, idBuffer)
		for(idTask = 0; idTask < numTasks; ++idTask){
				idBuffer = idTask % numBuffers;

				//Trace the jobs
				printf("Host thread %d \t sent job %d \t to buffer %d \t in device %d \n", omp_get_thread_num(), idTask, idBuffer, getIdDeviceBuffer(buffer[idBuffer]));

				//Fill the buffer (generate work)
				putIntoBuffer(buffer[idBuffer], &testData[idBuffer]);

				sendMyersBuffer(buffer[idBuffer], testData[idBuffer].totalQueriesEntries, testData[idBuffer].numQueries, testData[idBuffer].numCandidates);
				receiveMyersBuffer(buffer[idBuffer]);

				//Get the results from the buffer (consume results)
				getFromBuffer(buffer[idBuffer], &testData[idBuffer]);
		}

	ts1 = sampleTime();
	endMyers(&buffer);

	printf("TOTAL TIME: \t %f \n", ts1-ts);

	//Save last state of each buffer
	error = saveResults(qryFile, &testData[0], numBuffers);
		if(error != 0){fprintf(stderr, "Error %d, saving results \n", error); exit(EXIT_FAILURE);}

	freeTestData(&testData[0], numThreads);

    return (0);
}
