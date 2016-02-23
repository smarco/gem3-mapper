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

#define	REFERENCE_CHAR_LENGTH		  4
#define	UINT32_LENGTH				      32
#define	UINT64_LENGTH				      64
#define	REFERENCE_CHARS_PER_ENTRY	(UINT64_LENGTH / REFERENCE_CHAR_LENGTH)
#define REFERENCE_END_PADDING		  625
#define FILE_SIZE_LINES				    250

#define CATCH_ERROR(error)              {{if (error) { fprintf(stderr, "%s\n", processError(error)); exit(EXIT_FAILURE); }}}
#define	DIV_CEIL(NUMERATOR,DENOMINATOR) (((NUMERATOR)+((DENOMINATOR)-1))/(DENOMINATOR))
#define MIN(NUM_A, NUM_B)               ((NUM_A < NUM_B) ? NUM_A : NUM_B)

/* Encoded DNA Nucleotides */
#define ENC_DNA_CHAR_A    0LL
#define ENC_DNA_CHAR_C    1LL
#define ENC_DNA_CHAR_G    2LL
#define ENC_DNA_CHAR_T    3LL

#define ENC_DNA_CHAR_N    4LL
#define ENC_DNA_CHAR_SEP  5LL
#define ENC_DNA_CHAR_JUMP 6LL

typedef struct {
	uint64_t size;
	uint64_t numEntries;
	uint64_t *h_reference;
	uint64_t *d_reference;
	unsigned char *char_reference;
} reference_buffer_t;


uint64_t charToBinASCII(unsigned char base)
{
  switch(base)
  {
    case 'A':
    case 'a':
      return(ENC_DNA_CHAR_A);
    case 'C':
    case 'c':
      return(ENC_DNA_CHAR_C << (UINT64_LENGTH - REFERENCE_CHAR_LENGTH));
    case 'G':
    case 'g':
      return(ENC_DNA_CHAR_G << (UINT64_LENGTH - REFERENCE_CHAR_LENGTH));
    case 'T':
    case 't':
      return(ENC_DNA_CHAR_T << (UINT64_LENGTH - REFERENCE_CHAR_LENGTH));
    default :
      return(ENC_DNA_CHAR_N << (UINT64_LENGTH - REFERENCE_CHAR_LENGTH));
  }
}

uint32_t transformReferenceASCII(const char *referenceASCII, reference_buffer_t *reference)
{
	uint64_t indexBase, bitmap;
	uint64_t idEntry, i, referencePosition;
	unsigned char referenceChar;

	for(idEntry = 0; idEntry < reference->numEntries; ++idEntry){
		bitmap = 0;
		for(i = 0; i < REFERENCE_CHARS_PER_ENTRY; i++){
			referencePosition = (idEntry * REFERENCE_CHARS_PER_ENTRY) + i;
			if (referencePosition < reference->size) referenceChar = referenceASCII[referencePosition];
				else referenceChar = 'N'; //filling reference padding
			indexBase = charToBinASCII(referenceChar);
			bitmap = (bitmap >> REFERENCE_CHAR_LENGTH) | indexBase;
		}
		reference->h_reference[referencePosition / REFERENCE_CHARS_PER_ENTRY] = bitmap;
	}
	return(0);
}


uint32_t loadReferenceMFASTA(const char *fn, void *reference)
{
	reference_buffer_t *ref = (reference_buffer_t *) reference;
	FILE *fp = NULL;
	char lineFile[FILE_SIZE_LINES], *tmp_reference;
	uint64_t sizeFile = 0, position = 0;
	int32_t charsRead = 0;

	fp = fopen(fn, "rb");
	if (fp == NULL) return (30);

	fseek(fp, 0L, SEEK_END);
	sizeFile = ftell(fp);
	rewind(fp);

	tmp_reference = (char*) malloc(sizeFile * sizeof(char));
	if (ref == NULL) return (31);

	if ((fgets(lineFile, FILE_SIZE_LINES, fp) == NULL) || (lineFile[0] != '>'))
		return (32);

	while((!feof(fp)) && (fgets(lineFile, FILE_SIZE_LINES, fp) != NULL)){
		if (lineFile[0] != '>'){
			charsRead = strlen(lineFile);
			if(charsRead) charsRead--;
			memcpy((tmp_reference + position), lineFile, charsRead);
			position +=  charsRead;
		}
	}

	ref->size = position;
	ref->numEntries = DIV_CEIL(ref->size, REFERENCE_CHARS_PER_ENTRY) + REFERENCE_END_PADDING;
	printf("Reference size: %llu, Number of entries: %llu\n", ref->size, ref->numEntries);

	ref->h_reference = (uint64_t *) malloc(ref->numEntries * sizeof(uint64_t));
	if (ref->h_reference == NULL) return (31);

	transformReferenceASCII(tmp_reference, ref);

	fclose(fp);
	free(tmp_reference);
	return (0);
}

uint32_t saveRef(const char *fn, void *reference)
{
  reference_buffer_t *ref = (reference_buffer_t *) reference;

  char fileName[512];
  FILE *fp = NULL;
  uint64_t i;
  uint32_t error;

  sprintf(fileName, "%s.%llu.%ubits.ref", fn, ref->size, REFERENCE_CHAR_LENGTH);

  fp = fopen(fileName, "wb");
  if (fp == NULL) return (8);

  fwrite(&ref->numEntries, sizeof(uint64_t), 1, fp);
  fwrite(&ref->size, sizeof(uint64_t), 1, fp);

  fwrite(ref->h_reference, sizeof(uint64_t), ref->numEntries, fp);
  fclose(fp);
  return (0);
}

uint32_t freeReference(void *reference)
{   
  reference_buffer_t *ref = (reference_buffer_t *) reference;

  if(ref->h_reference != NULL){
    free(ref->h_reference);
    ref->h_reference = NULL;
  }

  if(ref->char_reference != NULL){
    free(ref->char_reference);
    ref->char_reference = NULL;
  }

  return(0);
}

char *processError(uint32_t e)
{
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

uint32_t initReference(void **reference)
{
  reference_buffer_t *ref = (reference_buffer_t *) malloc(sizeof(reference_buffer_t));
  ref->size = 0;
  ref->numEntries = 0;

  ref->h_reference = NULL;
  ref->d_reference = NULL;
  ref->char_reference = NULL;

  (*reference) = ref;
  return (0);
}

int32_t main(int argc, char *argv[])
{
  void *reference;
  char *refFile = argv[1];
  int error;

  error = initReference(&reference);
	CATCH_ERROR(error);

	error = loadReferenceMFASTA(refFile, reference);
  CATCH_ERROR(error);

  error = saveRef(refFile, reference);
  CATCH_ERROR(error);
    
  error = freeReference(reference);
  CATCH_ERROR(error);

  return (0);
}

