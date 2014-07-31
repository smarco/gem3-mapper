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

int loadQueries(const char *fn, test_t *testData)
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

	testData->queries = NULL;
	testData->candidates = NULL;
	testData->qinfo == NULL;

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

int saveResults(const char *fn, test_t *testData, uint32_t numBuffers)
{
	FILE *fp = NULL;
	char resFileOut[512];
	char cadena[512];
	uint32_t error, idTotalCandidate = 0, idBuffer, idCandidate;
	uint32_t totalResult = 0;

	sprintf(resFileOut, "%s.gem-int.pad.gpu", fn);

	fp = fopen(resFileOut, "wb");
		if (fp == NULL) return (47);

	for(idBuffer = 0; idBuffer < numBuffers; idBuffer++){
		totalResult += testData[idBuffer].numResults;
	}

	sprintf(cadena, "%u\n", totalResult);
	fputs(cadena, fp);

	for(idBuffer = 0; idBuffer < numBuffers; idBuffer++){
		for(idCandidate = 0; idCandidate < testData[idBuffer].numResults; idCandidate++){
			sprintf(cadena, "%u %u %u\n", idTotalCandidate, testData[idBuffer].results[idCandidate].score, testData[idBuffer].results[idCandidate].column);
			fputs(cadena, fp);
			idTotalCandidate++;
		}
	}

	fclose(fp);
	return (0);
}

int putIntoBuffer(void *buffer, test_t *testData)
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

int getFromBuffer(void *buffer, test_t *testData)
{
	uint32_t i;
	resEntry_t * results_buffer = getResultsBuffer(buffer);

	for(i = 0; i < testData->numResults; i++)
		testData->results[i] = results_buffer[i];

	return (0);
}

int loadTestData(unsigned char *qryFile, test_t *testData)
{
	int32_t error;

	error = loadQueries(qryFile, testData);
		if(error != 0){fprintf(stderr, "Error %d, loading queries \n", error); exit(EXIT_FAILURE);}
	testData->results = (resEntry_t *) malloc(testData->numCandidates * sizeof(resEntry_t));
		if (testData->queries == NULL){fprintf(stderr, "Error allocating results \n"); exit(EXIT_FAILURE);}
	return (0);
}

int freeTestData(test_t *testData)
{
	free(testData->queries);
	free(testData->candidates);
	free(testData->qinfo);
	free(testData->results);
	return (0);
}

/* Execution parameters:
 * 1)   Distance (overlap): 	0.2
 * 2)   Reference file path:	~/Myers/data/human_g1k_v37.cleaned.fasta.3101804739.4bits.ref
 * 3-N) Queries paths: 			 ~/Myers/data/1000.qry
 */

int main(int argc, char *argv[])
{
  void 		 		**buffer;
	int32_t		   idBuffer, numBuffers = argc - 3, error;
  float 		    distance = atof(argv[1]);
  unsigned char *refFile  = argv[2];
  double 		   ts, ts1;

	test_t			testData[numBuffers];
  unsigned char 	*qryFile[numBuffers];

	printf("numBuffers: %d \n", numBuffers);


	for(idBuffer = 0; idBuffer < numBuffers; idBuffer++){
		qryFile[idBuffer] = argv[idBuffer + 3];
		printf("qryFile: %s \n", qryFile[idBuffer]);
		loadTestData(qryFile[idBuffer], &testData[idBuffer]);
	}

	initMyers(&buffer, numBuffers, refFile, 1000, 10);
	ts = sampleTime();

	//Parallel region
	#pragma omp parallel for private(idBuffer)
	for(idBuffer = 0; idBuffer < numBuffers; idBuffer++){
		//Fill the buffer
		putIntoBuffer(buffer[idBuffer], &testData[idBuffer]);
		sendMyersBuffer(buffer[idBuffer], distance, testData[idBuffer].totalQueriesEntries, testData[idBuffer].numQueries, testData[idBuffer].numCandidates);
		receiveMyersBuffer(buffer[idBuffer]);

		//Get the results from the buffer
		getFromBuffer(buffer[idBuffer], &testData[idBuffer]);
	}

	ts1 = sampleTime();
	endMyers(&buffer);

	printf("TOTAL TIME: \t %f \n", ts1-ts);

  error = saveResults(qryFile[0], &testData[0], numBuffers);
		if(error != 0){fprintf(stderr, "Error %d, saving results \n", error); exit(EXIT_FAILURE);}

	for(idBuffer = 0; idBuffer < numBuffers; idBuffer++)
    	freeTestData(&testData[idBuffer]);

    return (0);
}
