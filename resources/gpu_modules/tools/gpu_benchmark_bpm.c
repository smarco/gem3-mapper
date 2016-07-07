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
	uint32_t 			      totalSizeQueries;
	uint32_t 			      totalQueriesEntries;
	uint32_t 			      numQueries;
	uint32_t 			      numCandidates;
	uint32_t 			      numResults;
	gpu_bpm_qry_entry_t	*queries;
	gpu_bpm_cand_info_t	*candidates;
	gpu_bpm_qry_info_t	*qinfo;
	gpu_bpm_alg_entry_t	*results;
	/* Debug fields */ 
	char 				        *raw_queries;
	uint32_t 			      *pos_raw_queries;
	uint32_t			      *GEM_score;
} test_t;

uint32_t loadQueries(const char *fn, test_t *testData, uint32_t *averageQuerySize, uint32_t *averageCandidatesPerQuery)
{
	FILE *fp = NULL;
	size_t result;
	uint32_t idQuery, dummy;

	fp = fopen(fn, "rb");
	if (fp == NULL) return (10);

	fread(&testData->totalSizeQueries, sizeof(uint32_t), 1, fp);
	fread(&testData->totalQueriesEntries, sizeof(uint32_t), 1, fp);
	fread(&dummy, sizeof(uint32_t), 1, fp);
	fread(&testData->numQueries, sizeof(uint32_t), 1, fp);
	fread(&testData->numCandidates, sizeof(uint32_t), 1, fp);
	fread(averageCandidatesPerQuery, sizeof(uint32_t), 1, fp);
	fread(averageQuerySize, sizeof(uint32_t), 1, fp);

	testData->numResults = testData->numCandidates;

	testData->queries = NULL;
	testData->candidates = NULL;
	testData->qinfo = NULL;

	testData->queries = (gpu_bpm_qry_entry_t *) malloc(testData->totalQueriesEntries * sizeof(gpu_bpm_qry_entry_t));
		if (testData->queries == NULL) {fprintf(stderr, "Error allocating queries \n"); exit(EXIT_FAILURE);}
	result = fread(testData->queries, sizeof(gpu_bpm_qry_entry_t), testData->totalQueriesEntries, fp);
		if (result != testData->totalQueriesEntries) {fprintf(stderr, "Error reading queries \n"); exit(EXIT_FAILURE);}

	testData->candidates = (gpu_bpm_cand_info_t *) malloc(testData->numCandidates * sizeof(gpu_bpm_cand_info_t));
		if (testData->candidates == NULL) {fprintf(stderr, "Error allocating candidates \n"); exit(EXIT_FAILURE);}
	result = fread(testData->candidates , sizeof(gpu_bpm_cand_info_t), testData->numCandidates, fp);
		if (result != testData->numCandidates) {fprintf(stderr, "Error reading candidates \n"); exit(EXIT_FAILURE);}

	testData->qinfo = (gpu_bpm_qry_info_t *) malloc(testData->numQueries * sizeof(gpu_bpm_qry_info_t));
		if (testData->qinfo == NULL) {fprintf(stderr, "Error allocating information for candidates \n"); exit(EXIT_FAILURE);}
	result = fread(testData->qinfo , sizeof(gpu_bpm_qry_info_t), testData->numQueries, fp);
		if (result != testData->numQueries) {fprintf(stderr, "Error reading information for candidates \n"); exit(EXIT_FAILURE);}

	testData->GEM_score = (uint32_t *) malloc(testData->numCandidates * sizeof(uint32_t));
		if (testData->GEM_score == NULL) {fprintf(stderr, "Error allocating GEM score results \n"); exit(EXIT_FAILURE);}
	result = fread(testData->GEM_score, sizeof(uint32_t), testData->numCandidates, fp);
		if (result != testData->numCandidates) {fprintf(stderr, "Error reading GEM score results \n"); exit(EXIT_FAILURE);}

	testData->raw_queries = (char *) malloc(testData->totalSizeQueries * sizeof(char));
		if (testData->raw_queries == NULL) {fprintf(stderr, "Error allocating RAW queries \n"); exit(EXIT_FAILURE);}
	result = fread(testData->raw_queries, sizeof(char), testData->totalSizeQueries, fp);
		if (result != testData->totalSizeQueries) {fprintf(stderr, "Error reading RAW queries \n"); exit(EXIT_FAILURE);}

	testData->pos_raw_queries = (uint32_t *) malloc(testData->numQueries * sizeof(uint32_t));
		if (testData->pos_raw_queries == NULL) {fprintf(stderr, "Error allocating information for RAW queries \n"); exit(EXIT_FAILURE);}
	result = fread(testData->pos_raw_queries , sizeof(uint32_t), testData->numQueries, fp);
		if (result != testData->numQueries) {fprintf(stderr, "Error reading information for RAW queries \n"); exit(EXIT_FAILURE);}


	fclose(fp);
	return (0);
}

uint32_t loadReference(char *fn, char **reference, uint32_t *size)
{
	FILE *fp = NULL;
	char lineFile[250], *tmp_reference = NULL;
	uint64_t sizeFile = 0, position = 0;
	int32_t charsRead = 0;

	fp = fopen(fn, "rb");
	if (fp == NULL) {fprintf(stderr, "Error opening file \n"); exit(EXIT_FAILURE);}

	fseek(fp, 0L, SEEK_END);
	sizeFile = ftell(fp);
	rewind(fp);

	tmp_reference = (char*) malloc(sizeFile * sizeof(char));
	if (tmp_reference == NULL) {fprintf(stderr, "Error allocating reference \n"); exit(EXIT_FAILURE);}

	if ((fgets(lineFile, 250, fp) == NULL) || (lineFile[0] != '>'))
	{fprintf(stderr, "Error reference file is not a Multifasta file \n"); exit(EXIT_FAILURE);}
	
	while((!feof(fp)) && (fgets(lineFile, 250, fp) != NULL)){
		if (lineFile[0] != '>'){
			charsRead = strlen(lineFile);
			if(charsRead) charsRead--;
			memcpy((tmp_reference + position), lineFile, charsRead);
			position +=  charsRead;
		}
	}

	(* size) = position;
	(* reference) = tmp_reference;

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
		if (fp == NULL) {fprintf(stderr, "Error opening file \n"); exit(EXIT_FAILURE);}

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

	gpu_bpm_qry_entry_t * 	queries_buffer 		= gpu_bpm_buffer_get_peq_entries_(buffer);
	gpu_bpm_cand_info_t *	candidates_buffer 	= gpu_bpm_buffer_get_candidates_(buffer);
	gpu_bpm_qry_info_t *	queryInfo_buffer	= gpu_bpm_buffer_get_peq_info_(buffer);

	testData->totalQueriesEntries 	= MIN(testData->totalQueriesEntries, gpu_bpm_buffer_get_max_peq_entries_(buffer));
	testData->numCandidates 		= MIN(testData->numCandidates, gpu_bpm_buffer_get_max_candidates_(buffer));
	testData->numQueries 			= MIN(testData->numQueries, gpu_bpm_buffer_get_max_queries_(buffer));
	testData->numResults 			= testData->numCandidates;

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
	gpu_bpm_alg_entry_t *results_buffer = gpu_bpm_buffer_get_alignments_(buffer);

	for(i = 0; i < testData->numResults; i++)
		testData->results[i] = results_buffer[i];

	return (0);
}

uint32_t loadTestData(char *qryFile, test_t *testData, uint32_t *averageQuerySize, uint32_t *averageCandidatesPerQuery)
{
	int32_t error;

	error = loadQueries(qryFile, testData, averageQuerySize, averageCandidatesPerQuery);
		if(error != 0){fprintf(stderr, "Error %d, loading queries \n", error); exit(EXIT_FAILURE);}
	testData->results = (gpu_bpm_alg_entry_t *) malloc(testData->numCandidates * sizeof(gpu_bpm_alg_entry_t));
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

void computeSW(char *query, char *candidate, uint32_t idCandidate, uint32_t sizeQuery, uint32_t sizeCandidate,
			  uint32_t sizeRef, uint32_t positionRef, uint32_t *matrix, gpu_bpm_alg_entry_t *results)
{
	int32_t sizeColumn = sizeQuery + 1, numColumns = sizeCandidate + 1;
	int32_t idColumn, i, y, j;
	int32_t cellLeft, cellUpper, cellDiagonal, delta;
	int32_t score, minScore = sizeQuery;
	int32_t minColumn = 0;
	char base;

	if((positionRef < sizeRef) && (sizeRef - positionRef) > sizeCandidate){

		// Horizontal initialization 
		for(i = 0; i < numColumns; i++)
			matrix[i] = 0;

		// Vertical initialization 
		for(i = 0; i < sizeColumn; i++)
			matrix[i * numColumns] = i;

		// Compute SW-MATRIX
		for(idColumn = 1; idColumn < numColumns; idColumn++){

			base = candidate[idColumn - 1];
			for(y = 1; y < sizeColumn; y++){

				delta = (base != query[y - 1]);
				cellLeft = matrix[y * numColumns + (idColumn - 1)] + 1;
				cellUpper = matrix[(y - 1) * numColumns + idColumn] + 1;
				cellDiagonal = matrix[(y - 1) * numColumns + (idColumn - 1)] + delta;
				
				score = MIN(cellDiagonal, MIN(cellLeft, cellUpper));
				matrix[y * numColumns + idColumn] = score;
			}

			if(score < minScore){
				minScore = score;
				minColumn = idColumn - 1;
			}
		}
    	results[idCandidate].column = minColumn;
    	results[idCandidate].score = minScore;
	}
}

double processMyersGPU(char *refFile, char *qryFile, uint32_t numBuffers, uint32_t numThreads, float maxMbPerBuffer, uint32_t numTasks)
{
	test_t				testData[numThreads];
	uint32_t 			averageQuerySize, averageCandidatesPerQuery;
	uint32_t 			error, idBuffer, iteration, threadID, idTask;
	double				ts, ts1;

	gpu_buffers_dto_t 	buff  = {.buffer 		 		        = NULL,
								               .numBuffers 	 		      = numBuffers,
								               .maxMbPerBuffer 		    = maxMbPerBuffer,
								               .activeModules			    = GPU_BPM};
	gpu_index_dto_t 	  index = {.fmi 			 		        = NULL,
								               .indexCoding 	 		    = GPU_INDEX_NONE,
								               .bwtSize 		 		      = 0};
	gpu_reference_dto_t ref   = {.reference 	 		      = refFile,
								               .refCoding 			      = GPU_REF_PROFILE_FILE,
								               .refSize				        = 0};
	gpu_info_dto_t 		  sys	  = {.selectedArchitectures = GPU_ARCH_SUPPORTED,
								               .userAllocOption   	  = GPU_LOCAL_OR_REMOTE_DATA,
	                             .activatedModules      = GPU_NONE_MODULES,
	                             .allocatedStructures   = GPU_NONE_MODULES
	                             .verbose               = verbose};

	for(threadID = 0; threadID < numThreads; ++threadID)
		loadTestData(qryFile, &testData[threadID], &averageQuerySize, &averageCandidatesPerQuery);

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
				gpu_bpm_init_buffer_(buff.buffer[idBuffer], averageQuerySize, averageCandidatesPerQuery);
				putIntoBuffer(buff.buffer[idBuffer], &testData[idBuffer]);
				gpu_bpm_send_buffer_(buff.buffer[idBuffer], testData[idBuffer].totalQueriesEntries, testData[idBuffer].numQueries, testData[idBuffer].numCandidates, 0);
				gpu_bpm_receive_buffer_(buff.buffer[idBuffer]);

				//Get the results from the buffer (consume results)
				getFromBuffer(buff.buffer[idBuffer], &testData[idBuffer]);
		}

	ts1 = sample_time();
	gpu_destroy_buffers_(&buff);

	//Save last state of each buffer
	error = saveResults(qryFile, &testData[0], numBuffers);
		if(error != 0){fprintf(stderr, "Error %d, saving results \n", error); exit(EXIT_FAILURE);}

	freeTestData(&testData[0], numThreads);

	return(ts1-ts);
}

double processSmithWatermanCPU(char *refFile, char *qryFile, uint32_t numThreads, uint32_t numBuffers)
{
	test_t		testData;
	double		ts, ts1;
	uint32_t 	**matrix;
	char 	 	*query = NULL, *candidate = NULL, *reference = NULL;
	uint32_t 	sizeRef, error, tID;
	uint32_t 	averageQuerySize, averageCandidatesPerQuery;

	loadTestData(qryFile, &testData, &averageQuerySize, &averageCandidatesPerQuery);
	loadReference(refFile, &reference, &sizeRef);	

	matrix = (uint32_t **) malloc(numThreads * sizeof(uint32_t *));
	for(tID = 0; tID < numThreads; tID++){
		matrix[tID] = (uint32_t *) malloc((averageQuerySize * 2) * (averageQuerySize * 4) * sizeof(uint32_t));
		if (matrix[tID] == NULL) {fprintf(stderr, "Error allocating SW matrix \n"); exit(EXIT_FAILURE);}
	}

	ts = sample_time();

	#pragma omp parallel num_threads(numThreads)

	{
		uint32_t threadID = omp_get_thread_num();
		uint32_t idQuery, idCandidate, sizeQuery, sizeCandidate, positionCandidate;
		char 	 *query = NULL, *candidate = NULL;

		#pragma omp for schedule(static) 
		for (idCandidate = 0; idCandidate < testData.numCandidates; idCandidate++){

			idQuery 			= testData.candidates[idCandidate].query;
			query				= testData.raw_queries + testData.pos_raw_queries[idQuery];
			candidate			= reference + testData.candidates[idCandidate].position;
			sizeQuery 			= testData.qinfo[idQuery].size;
			sizeCandidate 		= testData.candidates[idCandidate].size;
			positionCandidate 	= testData.candidates[idCandidate].position;

			computeSW(query, candidate, idCandidate, sizeQuery, sizeCandidate, sizeRef, positionCandidate, matrix[threadID], testData.results);
		}
	}

	ts1 = sample_time();

	error = saveResults(qryFile, &testData, numBuffers);
	if(error != 0) {fprintf(stderr, "Error %d, saving results \n", error); exit(EXIT_FAILURE);}

	freeTestData(&testData, numBuffers);

	return(ts1-ts);
}



uint32_t main(int argc, char *argv[])
{
	void 		 	**buffer;
	uint32_t		numBuffers, numThreads, idBuffer, maxMbPerBuffer, numTasks, error;
	char 			*refFile, *qryFile;
	double 			timeElapsed;

	refFile  = argv[2];
	qryFile  = argv[3];
	
	/* Smith and Waterman */
	/* Execution parameters:
	 * 1) 0 
	 * 2) Reference file path:	example: ~/Myers/data/human_g1k_v37.cleaned.fasta (in Multi-Fasta)
	 * 3) Queries paths: 		example: ~/Myers/data/1000.qry
	 * 4) Number threads		example: 8
	 */
	if(atoi(argv[1]) == 0){
		numBuffers 	= 1;
		numThreads 	= atoi(argv[4]);
		//numThreads 	= 1;
		numTasks 	= 1;
	 	timeElapsed = processSmithWatermanCPU(refFile, qryFile, numThreads, numBuffers);
	}	

	/* Myers-GPU */
	/* Execution parameters:
	 * 1) 1 
	 * 2) Reference file path:	example: ~/Myers/data/human_g1k_v37.cleaned.fasta.3101804739.4bits.ref
	 * 3) Queries paths: 		example: ~/Myers/data/1000.qry
	 * 4) Number of buffers		example: 8
	 * 5) Maximum buffer size	example: 32  (in MBytes)
	 * 6) Number of tasks		example: 120
	 */
	if(atoi(argv[1]) == 1){
		numBuffers 		= atoi(argv[4]);
		numThreads 		= atoi(argv[4]);
		maxMbPerBuffer 	= atoi(argv[5]);
		numTasks 		= atoi(argv[6]);

	 	timeElapsed 	= processMyersGPU(refFile, qryFile, numBuffers, numThreads, (float) maxMbPerBuffer, numTasks);
	}

	/* Myers-GPU */
	if(atoi(argv[1]) > 1){
		printf("Option for dynamic programing not implemented: exiting \n");
		exit(EXIT_FAILURE);		
	}	

	printf("TOTAL TIME: \t %f \n", timeElapsed);

    return (0);
}
