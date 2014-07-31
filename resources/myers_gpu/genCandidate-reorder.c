#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define		NUM_BASES		5
#define		MAX_SIZE_LINE	20000
#define		NUM_BASES_ENTRY	128
#define		SIZE_HW_WORD	32
#define		NUM_SUB_ENTRIES (NUM_BASES_ENTRY / SIZE_HW_WORD)

#define CATCH_ERROR(error) {{if (error) { fprintf(stderr, "%s\n", processError(error)); exit(EXIT_FAILURE); }}}

#ifndef MAX
	#define MAX(_a, _b) (((_a) > (_b)) ? (_a) : (_b))
#endif 

#ifndef MIN
	#define MIN(_a, _b) (((_a) < (_b)) ? (_a) : (_b))
#endif


typedef struct {
	uint bitmap[NUM_BASES][NUM_SUB_ENTRIES];
} qryEntry_t;

typedef struct {
	uint posEntry;
	uint size;
} qryInfo_t;

typedef struct {
	uint query;
	uint position;
} candInfo_t;

typedef struct {
	uint query;
	unsigned long long int position;
} candInfo_GEM_t;

typedef struct {
	uint totalSizeQueries;
	uint totalQueriesEntries;
	uint sizeQueries;
	uint numQueries;
	uint numCandidates;

	unsigned char *char_queries;
	qryEntry_t *h_queries;
	qryEntry_t *d_queries;
	candInfo_t *h_candidates;
	candInfo_t *d_candidates;
	qryInfo_t *h_infoQueries;
	qryInfo_t *d_infoQueries;
	uint *h_results;	
} qry_t;

int countCandidates(unsigned char* cadena, uint *numberOfCandidates, uint *querySize)
{
	int localNumCandidates = 0, i = 0, flagQueryMesured = 0;

	while(cadena[i] != '\n'){
		if(cadena[i] == '\t'){
			if(flagQueryMesured == 0){
				(* querySize) = i;
				flagQueryMesured = 1;
			}
			localNumCandidates++;
		}
		i++;
	}
	(* numberOfCandidates) = localNumCandidates;
	return (0);
}

int countQueries(FILE *fp, void *queries)
{   //Read all the file and return the number of candidates and queries   
	qry_t *qry = (qry_t *) queries;
	unsigned char cadena[MAX_SIZE_LINE];
	uint numberOfCandidates = 0, numberOfQueries = 0, sizeAllQueries = 0;
	uint localCandidates = 0;
	uint currentQuerySize = 0, lastQuerySize = 0;

	rewind(fp);
	
	if (fgets(cadena, MAX_SIZE_LINE, fp) == NULL) 
		return (32);

	countCandidates(cadena, &localCandidates, &currentQuerySize);
	lastQuerySize = currentQuerySize;
	numberOfCandidates += localCandidates;
	sizeAllQueries += currentQuerySize;
	numberOfQueries++;

	while((!feof(fp)) && (fgets(cadena, MAX_SIZE_LINE, fp) != NULL)){
		countCandidates(cadena, &localCandidates, &currentQuerySize); 

		//if(lastQuerySize != currentQuerySize)
		//	lastQuerySize = 0;

		numberOfCandidates += localCandidates;
		sizeAllQueries += currentQuerySize;
		numberOfQueries++;
	}

	qry->sizeQueries = 0;
	qry->numQueries = numberOfQueries;
	qry->numCandidates = numberOfCandidates;
	qry->totalSizeQueries = sizeAllQueries;

	return (0);
}

int processCandidate(char *candidate, uint *pos, unsigned int *res)
{
	unsigned char area[50];
	unsigned char symbol[10];
	unsigned char position[15];
	unsigned char result[500];
	int tokenCount;

	tokenCount = sscanf(candidate,"%49[^:]:%9[^:]:%14[^:]:%s", area, symbol, position, result);
	if (tokenCount != 4) 
	{ fprintf(stderr, "processCandidate: Something went wrong\n"); exit(42); }

	(* pos) = (uint) atoi(position);
	(* res) = (uint) atoi(result);
		
	return(0);
}

int processQuery(uint queryNumber, unsigned char *textLine, 
				 unsigned char *queries, candInfo_t *listCandidates, uint *listResults, 
				 uint *retSizeQueries, uint *retNumCandidates)
{
	uint position,
		 numCandidates = 0,
		 sizeQuery,
		 result,
		 tokenCount;

	int	 sizeLastCandidates,
		 sizeCurrentCandidate;

	unsigned char *pLastCandidates;

	unsigned char query[MAX_SIZE_LINE],
	  	      lastCandidates[MAX_SIZE_LINE],
	  	      currentCandidate[500];

	pLastCandidates = &lastCandidates[0];

  	sscanf(textLine, "%s\t%[^\n]", query, pLastCandidates);
	sizeQuery = strlen(query);
	sizeLastCandidates = strlen(pLastCandidates);
	memcpy(queries, query, sizeQuery);
	
	while (sizeLastCandidates > 0){
		tokenCount = sscanf(pLastCandidates,"%s", currentCandidate);
		if (tokenCount < 1) { fprintf(stderr, "processQuery: Something went wrong\n"); exit(42); }

		processCandidate(currentCandidate, &position, &result);
		listCandidates[numCandidates].query = queryNumber;
		listCandidates[numCandidates].position = position;
		listResults[numCandidates] = result;
		numCandidates++;

		//update values next iteration
		sizeCurrentCandidate = strlen(currentCandidate) + 1;
		pLastCandidates += sizeCurrentCandidate;
		sizeLastCandidates -= sizeCurrentCandidate;
	}

	//return the number of data used
	(* retSizeQueries) 	  = sizeQuery;
	(* retNumCandidates)  =	numCandidates;

	return (0);
}

int loadQueries(const unsigned char *fn, void *queries)
{      
	qry_t *qry = (qry_t *) queries;
	FILE *fp = NULL;
	unsigned char textLine[MAX_SIZE_LINE];
	uint queryNumber = 0,
		 sizeCurrentQuery, 
		 numQueryCandidates,
		 processedCandidates = 0, 
		 processedBases = 0,
		 lastEntriesPerQuery = 0;

	fp = fopen(fn, "rb");
	if (fp==NULL) return (30);

	countQueries(fp, (void *) qry);
	rewind(fp);

	qry->char_queries 	= (unsigned char *) malloc(qry->totalSizeQueries * sizeof(unsigned char));
		if (qry->char_queries == NULL) return (31);
	qry->h_candidates 	= (candInfo_t *) malloc(qry->numCandidates * sizeof(candInfo_t));
		if (qry->h_candidates == NULL) return (31);
	qry->h_infoQueries 	= (qryInfo_t *) malloc(qry->numQueries * sizeof(qryInfo_t));
		if (qry->h_infoQueries == NULL) return (31);
	qry->h_results 		= (uint *) malloc(qry->numCandidates * sizeof(uint));
		if (qry->h_results == NULL) return (31);
	
	if (fgets(textLine, MAX_SIZE_LINE, fp) == NULL) 
		return (32);

	processQuery(queryNumber, textLine, 
				 qry->char_queries, qry->h_candidates, qry->h_results, 
				 &sizeCurrentQuery, &numQueryCandidates);
	qry->h_infoQueries[queryNumber].size = sizeCurrentQuery;
	processedBases = sizeCurrentQuery;
	qry->h_infoQueries[queryNumber].posEntry = lastEntriesPerQuery;
	lastEntriesPerQuery += (sizeCurrentQuery / NUM_BASES_ENTRY) + ((sizeCurrentQuery % NUM_BASES_ENTRY) ? 1 : 0);
	processedCandidates += numQueryCandidates;
	queryNumber++;

	while((!feof(fp)) && (fgets(textLine, MAX_SIZE_LINE, fp) != NULL)){
		processQuery(queryNumber, textLine,
					 qry->char_queries + processedBases, qry->h_candidates + processedCandidates, qry->h_results + processedCandidates, 
				     &sizeCurrentQuery, &numQueryCandidates);
		qry->h_infoQueries[queryNumber].size = sizeCurrentQuery;
		processedBases += sizeCurrentQuery;
		qry->h_infoQueries[queryNumber].posEntry += lastEntriesPerQuery;
		lastEntriesPerQuery += (sizeCurrentQuery / NUM_BASES_ENTRY) + ((sizeCurrentQuery % NUM_BASES_ENTRY) ? 1 : 0);
		processedCandidates += numQueryCandidates;
		queryNumber++;
	}

	fclose(fp);
	return (0);
}

uint base2number(unsigned char base)
{
	switch(base)
	{
    	case 'A':
    	case 'a':
    	    return(1);
    	case 'C':
    	case 'c':
    	    return(2);
    	case 'G':
    	case 'g':
    	    return(4);
    	case 'T':
    	case 't':
    	    return(8);
    	default :
    	    return(16);
	}
}

int char2bin(unsigned char *query, uint subBMP, qryEntry_t *binaryQuery, int numBases)
{
	int i; 
	uint indexBase;
	uint bitmapA 	= 0;
	uint bitmapC 	= 0;
	uint bitmapG 	= 0;
	uint bitmapT 	= 0;
	uint bitmapN 	= 0;

	for(i = 0; i < numBases; i++){
		indexBase = base2number(query[i]);
		bitmapA |= ((indexBase & 1 )      ) << i;
		bitmapC |= ((indexBase & 2 )  >> 1) << i;
		bitmapG |= ((indexBase & 4 )  >> 2) << i;
		bitmapT |= ((indexBase & 8 )  >> 3) << i;
		bitmapN |= ((indexBase & 16)  >> 4) << i;
	}

	binaryQuery->bitmap[0][subBMP] = bitmapA;
	binaryQuery->bitmap[1][subBMP] = bitmapC;
	binaryQuery->bitmap[2][subBMP] = bitmapG;
	binaryQuery->bitmap[3][subBMP] = bitmapT;
	binaryQuery->bitmap[4][subBMP] = bitmapN;

	return (0);
}

int transformQueries(void *queries)
{
	qry_t *qry = (qry_t *) queries;
	uint numEntriesPerQuery;
	uint pos, idQuery, intraQuery;
	uint processedBases = 0;
	uint startQueryEntry, intraBMP;
	int intraBasesProcessed, sizeQuery, numBases;

	for(idQuery = 0; idQuery < qry->numQueries; idQuery++)
		qry->totalQueriesEntries += (qry->h_infoQueries[idQuery].size / NUM_BASES_ENTRY) + ((qry->h_infoQueries[idQuery].size % NUM_BASES_ENTRY) ? 1 : 0);

	qry->h_queries 	= (qryEntry_t *) malloc(qry->totalQueriesEntries * sizeof(qryEntry_t));
		if (qry->h_queries == NULL) return (31);

	for(idQuery = 0; idQuery < qry->numQueries; idQuery++){
		intraQuery = 0;
		startQueryEntry = qry->h_infoQueries[idQuery].posEntry;
		sizeQuery = qry->h_infoQueries[idQuery].size;
		for(pos = 0; pos < sizeQuery; pos += NUM_BASES_ENTRY){
			for(intraBMP = 0; intraBMP < NUM_SUB_ENTRIES; intraBMP++){
				intraBasesProcessed = SIZE_HW_WORD * intraBMP;
				numBases = MIN(sizeQuery - (int)(pos + intraBasesProcessed), SIZE_HW_WORD);
				numBases = (numBases < 0) ? 0 : numBases; 
				char2bin(qry->char_queries + processedBases + pos + intraBasesProcessed, intraBMP, qry->h_queries + startQueryEntry + intraQuery, numBases);
			}
			intraQuery++;
		}
		processedBases += sizeQuery;
	}

	return (0);
}

int stadistics(void *queries)
{
	qry_t *qry = (qry_t *) queries;
	uint idHistogram, idQuery;
	uint queryHistogram[10000];
	uint numelements = 10000;

	uint minQueryLong = qry->h_infoQueries[0].size;
	uint maxQueryLong = qry->h_infoQueries[0].size;

	for(idHistogram = 0; idHistogram < numelements; idHistogram++)
		queryHistogram[idHistogram] = 0;

	for(idQuery = 0; idQuery < qry->numQueries; idQuery++){
		maxQueryLong = MAX(qry->h_infoQueries[idQuery].size, maxQueryLong);
		minQueryLong = MIN(qry->h_infoQueries[idQuery].size, minQueryLong);
	}

	for(idQuery = 0; idQuery < qry->numQueries; idQuery++)
		queryHistogram[qry->h_infoQueries[idQuery].size]++;

	printf("Number of queries: %d\n", qry->numQueries);
	printf("Number of candidates: %d\n", qry->numCandidates);
	printf("Ratio: %.2f\n", (float)qry->numCandidates/(float)qry->numQueries);
	printf("Minimum Query Longitud %d - Maximum Query Longitud %d\n", minQueryLong, maxQueryLong);

	for(idHistogram = 0; idHistogram < maxQueryLong; idHistogram++){
		printf("%d\n",queryHistogram[idHistogram]);
	}

	return (0);
}

int saveQueries(const char *fn, void *query)
{
    qry_t *qry = (qry_t *) query;
	candInfo_GEM_t *candidates_GEM = NULL; 
	uint idCandidate;

    char qryFileOut[512];
    FILE *fp = NULL;

    sprintf(qryFileOut, "%s.%u.%u.%u.gem.qry", fn, qry->numQueries, qry->numCandidates / qry->numQueries, qry->sizeQueries);
    
    fp = fopen(qryFileOut, "wb");
    	if (fp == NULL) return (8);

	//Temporal Fix
	candidates_GEM = (candInfo_GEM_t *) malloc(qry->numCandidates * sizeof(candInfo_GEM_t));
		if (candidates_GEM == NULL) return (31);
	for(idCandidate = 0; idCandidate < qry->numCandidates; idCandidate++){
		candidates_GEM[idCandidate].position = (unsigned long long int) qry->h_candidates[idCandidate].position;
		candidates_GEM[idCandidate].query = qry->h_candidates[idCandidate].query;
	}

    fwrite(&qry->totalSizeQueries, 		sizeof(uint), 1, fp);
    fwrite(&qry->totalQueriesEntries, 	sizeof(uint), 1, fp);
    fwrite(&qry->sizeQueries, 			sizeof(uint), 1, fp);
    fwrite(&qry->numQueries, 			sizeof(uint), 1, fp);
    fwrite(&qry->numCandidates, 		sizeof(uint), 1, fp);

    fwrite(qry->h_queries, 		sizeof(qryEntry_t), 	qry->totalQueriesEntries, 	fp);
    fwrite(candidates_GEM, 		sizeof(candInfo_GEM_t),	qry->numCandidates, 		fp);
    fwrite(qry->h_infoQueries, 	sizeof(qryInfo_t), 		qry->numQueries, 			fp);
    fwrite(qry->h_results, 		sizeof(uint), 			qry->numCandidates, 		fp);
    fwrite(qry->char_queries, 	sizeof(char),		 	qry->totalSizeQueries, 		fp);

	free(candidates_GEM);
    fclose(fp);
    return (0);
}

int freeQueries(void *queries)
{   
    qry_t *qry = (qry_t *) queries;  
	
    if(qry->char_queries != NULL){
        free(qry->char_queries);
        qry->char_queries = NULL;
    }  

    if(qry->h_queries != NULL){
        free(qry->h_queries);
        qry->h_queries = NULL;
    }

    if(qry->h_candidates != NULL){
        free(qry->h_candidates);
        qry->h_candidates = NULL;
    }  

    if(qry->h_infoQueries != NULL){
        free(qry->h_infoQueries);
        qry->h_infoQueries = NULL;
    }  
 
     if(qry->h_results != NULL){
        free(qry->h_results);
        qry->h_results = NULL;
    }
	
    return(0);
}

char *processError(int e){ 
    switch(e) {
        case 0:  return "No error"; break; 
        case 30: return "Cannot open reference file"; break;
        case 31: return "Cannot allocate reference"; break;
        case 32: return "Reference file isn't multifasta format"; break;
        case 37: return "Cannot open reference file on write mode"; break;
        case 42: return "Cannot open queries file"; break;
        case 43: return "Cannot allocate queries"; break;
        case 45: return "Cannot allocate results"; break;
        case 47: return "Cannot open results file for save intervals"; break;
        case 48: return "Cannot open results file for load intervals"; break;
        case 99: return "Not implemented"; break;
        default: return "Unknown error";
    }   
}

int initQueries(void **queries)
{
    qry_t *qry = (qry_t *) malloc(sizeof(qry_t));
	qry->totalSizeQueries = 0;
	qry->totalQueriesEntries = 0;
	qry->sizeQueries = 0;
	qry->numQueries = 0;
	qry->numCandidates = 0;

	qry->char_queries = NULL;
	qry->h_queries = NULL;
	qry->d_queries = NULL;
	qry->h_infoQueries = NULL;
	qry->d_infoQueries = NULL;
	qry->h_candidates = NULL;
	qry->d_candidates = NULL;
	qry->h_results = NULL;

    (*queries) = qry;
    return (0);
}

int main(int argc, char *argv[])
{
	if(argc < 2) {
		printf("Usage: gcandidates regionsfile.prof\
		     \n example: bin/gcandidates input/regions-score.1M.prof \n");
		exit(0);
	}

    void *queries;
    unsigned char *qryFile = argv[1];
    int error;

	error = initQueries(&queries);    
	CATCH_ERROR(error);

	error = loadQueries(qryFile, queries);
    CATCH_ERROR(error);
	
	error = transformQueries(queries);
    CATCH_ERROR(error);

    error = saveQueries(qryFile, queries);
    CATCH_ERROR(error);
	
	//error = stadistics(queries);
    //CATCH_ERROR(error);
    
    error = freeQueries(queries);
    CATCH_ERROR(error);

    return (0);
}

