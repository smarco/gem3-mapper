/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define	BWT_CHAR_LENGTH				  3
#define FMI_NUM_COUNTERS			  4
#define FMI_ALTERNATE_COUNTERS	2
#define FMI_ENTRY_SIZE				  128

#define	UINT32_LENGTH				    32
#define	UINT64_LENGTH				    64
#define FILE_SIZE_LINES				  250

#define CATCH_ERROR(error)              {{if (error) { fprintf(stderr, "%s\n", processError(error)); exit(EXIT_FAILURE); }}}
#define	DIV_CEIL(NUMERATOR,DENOMINATOR) (((NUMERATOR)+((DENOMINATOR)-1))/(DENOMINATOR))
#define MIN(NUM_A, NUM_B) 				      ((NUM_A < NUM_B) ? NUM_A : NUM_B)

/* Encoded DNA Nucleotides */
#define ENC_DNA_CHAR_A          0
#define ENC_DNA_CHAR_C          1
#define ENC_DNA_CHAR_G          2
#define ENC_DNA_CHAR_T          3

//whatever different than the above bases
#define ENC_DNA_CHAR_X          4

// FMI Entry (64 Bytes) using:
//    4 letters: Alternate counters   (2  uint64_t)
//    5 letters: 12 Bitmaps x 32 bits (3 uint128_t)
typedef struct {
	uint64_t counters[FMI_NUM_COUNTERS / FMI_ALTERNATE_COUNTERS];
	uint32_t bitmaps[FMI_ENTRY_SIZE * BWT_CHAR_LENGTH / UINT32_LENGTH];
} fmi_entry_t;

typedef struct {
	uint32_t bitmaps[BWT_CHAR_LENGTH];
} bwt_entry_t;

typedef struct {
	uint64_t counters[FMI_NUM_COUNTERS];
} fmi_counters_t;

typedef struct {
	uint64_t        size;
	uint64_t        numEntries;
	fmi_entry_t    *h_fmi;
	fmi_counters_t *h_countersFMI;
	bwt_entry_t    *h_bitmapBWT;
	unsigned char  *h_asciiBWT;
} fmi_buffer_t;

inline
uint32_t charToBinASCII(unsigned char base)
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
      return(ENC_DNA_CHAR_X);
  }
}

inline char decodeBase(uint32_t indexBase)
{
  switch(indexBase)
  {
    case ENC_DNA_CHAR_A:
      return('A');
    case ENC_DNA_CHAR_C:
      return('C');
    case ENC_DNA_CHAR_G:
      return('G');
    case ENC_DNA_CHAR_T:
      return('T');
    default:
    return('N');
  }
}

void encodeBWTtoEntryPEQ(bwt_entry_t *bwt_entry, unsigned char base, uint32_t size)
{
	uint32_t idBase, binBase = charToBinASCII(base);
	for(idBase = 0; idBase < BWT_CHAR_LENGTH; ++idBase){
		bwt_entry->bitmaps[idBase] |= ((binBase >> idBase) & 0x1) << (UINT32_LENGTH - size - 1);
	}
}

inline 
void countBasesASCII(fmi_counters_t *counterEntry, unsigned char base)
{
	switch(base){
    case 'A':
    case 'a':
      counterEntry->counters[0]++;
      break;
    case 'C':
    case 'c':
      counterEntry->counters[1]++;
      break;
    case 'G':
    case 'g':
      counterEntry->counters[2]++;
      break;
    case 'T':
    case 't':
      counterEntry->counters[3]++;
    break;
	}
}

uint32_t transformBWTtoPEQ(fmi_buffer_t *fmi)
{
	uint64_t idEntry, i, bwtPosition;
	unsigned char bwtChar;

	//Padded to FMI_ENTRY_SIZE (128 bases)
	uint32_t bwtNumEntries = fmi->numEntries * (FMI_ENTRY_SIZE / UINT32_LENGTH);

	fmi->h_bitmapBWT = (bwt_entry_t *) malloc(bwtNumEntries * sizeof(bwt_entry_t));
	if (fmi->h_bitmapBWT == NULL) return (31);

	//printf("[READING BWT]: \n");
	for(idEntry = 0; idEntry < bwtNumEntries; ++idEntry){
		//printf("[%2lld] ", idEntry);
		for(i = 0; i < UINT32_LENGTH; ++i){
			bwtPosition = (idEntry * UINT32_LENGTH) + i;
			if (bwtPosition < fmi->size) bwtChar = fmi->h_asciiBWT[bwtPosition];
				else bwtChar = 'N'; //filling bwt padding
			//printf("%c", bwtChar);
			encodeBWTtoEntryPEQ(&fmi->h_bitmapBWT[idEntry], bwtChar, i);
		}
		//printf("\n ");
	}
	return(0);
}

uint32_t buildCounters(fmi_buffer_t *fmi)
{
	uint64_t idEntry, i, bwtPosition;
	uint32_t idBase, idCounter;
	unsigned char bwtChar;

	fmi_counters_t  localCounters;
	uint32_t        countersNumEntries = fmi->numEntries * (FMI_ENTRY_SIZE / UINT32_LENGTH);

	fmi->h_countersFMI = (fmi_counters_t *) malloc(countersNumEntries * sizeof(fmi_counters_t));
	if (fmi->h_countersFMI == NULL) return (31);

	// Init first local BWT entry
	for(idEntry = 0; idEntry < countersNumEntries; ++idEntry){
		for(idBase = 0; idBase < FMI_NUM_COUNTERS; ++idBase){
			fmi->h_countersFMI[idEntry].counters[idBase] = 0;
		}
	}

	// Accumulate values locally for each BWT entry
	for(idEntry = 0; idEntry < countersNumEntries - 1; ++idEntry){
		for(i = 0; i < UINT32_LENGTH; ++i){
			bwtPosition = (idEntry * UINT32_LENGTH) + i;
			if (bwtPosition < fmi->size) bwtChar = fmi->h_asciiBWT[bwtPosition];
				else bwtChar = 'N'; //filling bwt padding (will be N)
			countBasesASCII(&fmi->h_countersFMI[idEntry + 1], bwtChar);
		}
	}

	// Accumulate values globally for all the BWT
	for(idEntry = 1; idEntry < countersNumEntries; ++idEntry){
		for(idBase = 0; idBase < FMI_NUM_COUNTERS; ++idBase){
			fmi->h_countersFMI[idEntry].counters[idBase] += fmi->h_countersFMI[idEntry - 1].counters[idBase];
		}
	}

	// Prepare the local counters with the previous accumulative letters
	localCounters.counters[0] = 0;
	for(idBase = 1; idBase < FMI_NUM_COUNTERS; ++idBase){
		localCounters.counters[idBase] = localCounters.counters[idBase - 1] + fmi->h_countersFMI[countersNumEntries - 4].counters[idBase - 1];
	}

	// Accumulate the previous alphabet letters to the global counters
	for(idEntry = 0; idEntry < countersNumEntries; ++idEntry){
		for(idBase = 1; idBase < FMI_NUM_COUNTERS; ++idBase){
			fmi->h_countersFMI[idEntry].counters[idBase] += localCounters.counters[idBase];
		}
	}

	return (0);
}

void setCounterLayout(fmi_counters_t *h_countersFMI, fmi_entry_t *h_fmi, const uint64_t idEntry)
{
	uint32_t idCounter;
	const uint32_t FMI_COUNTERS_PER_ENTRY = FMI_NUM_COUNTERS / FMI_ALTERNATE_COUNTERS; // (4 letters) => 2 counters / fmi_entry

	for(idCounter = 0; idCounter < FMI_COUNTERS_PER_ENTRY; ++idCounter){
		const uint32_t FMI_SET_ALT_COUNTER = idEntry % FMI_COUNTERS_PER_ENTRY;
		h_fmi->counters[idCounter] = h_countersFMI->counters[(FMI_SET_ALT_COUNTER * FMI_COUNTERS_PER_ENTRY) + idCounter];
	}
}

void setBitmapLayout(bwt_entry_t *h_bitmapBWT, fmi_entry_t *h_fmi)
{
	const uint32_t LUT[12] = {3,7,11,0,1,2,4,5,6,8,9,10};  									  // 1 1 (1) 0 2 2 (2) 0 3 3 (3) (0)
	uint32_t idPacket, idFMIBucket, padding = 0;
	const uint32_t FMI_NUM_BITMAPS        = FMI_ENTRY_SIZE * BWT_CHAR_LENGTH / UINT32_LENGTH; // 12 words             (4 FMI entries x 3 bits)
	const uint32_t NUM_PACKET_BMP_ENTRIES = FMI_ENTRY_SIZE  / UINT32_LENGTH;                  // 4 BMP entries        (128 bases / 32 bits) 
	const uint32_t NUM_PACKET_FMI_ENTRIES = FMI_NUM_BITMAPS / NUM_PACKET_BMP_ENTRIES;         // 3 FMI bitmap entries (12 bitmaps / 4 packets) 

	for(idFMIBucket = 0; idFMIBucket < NUM_PACKET_BMP_ENTRIES; ++idFMIBucket)                 // Iterate over BMP entries (4)
		h_bitmapBWT[idFMIBucket].bitmaps[NUM_PACKET_FMI_ENTRIES - 1] = ~ h_bitmapBWT[idFMIBucket].bitmaps[NUM_PACKET_FMI_ENTRIES - 1];

		for(idPacket = 0; idPacket < NUM_PACKET_BMP_ENTRIES; ++idPacket){          		      // Iterate over BMP entries (4)
			for(idFMIBucket = 0; idFMIBucket < NUM_PACKET_FMI_ENTRIES; ++idFMIBucket){        // Iterate over FMI bitmaps (3)
				h_fmi->bitmaps[LUT[idPacket * NUM_PACKET_FMI_ENTRIES + idFMIBucket]] =  h_bitmapBWT[idPacket].bitmaps[idFMIBucket];
		}
	}
}


uint32_t buildFMI(fmi_buffer_t *fmi)
{
	const uint32_t BITMAPS_PER_FMI = FMI_ENTRY_SIZE / UINT32_LENGTH; // 4 FMI bitmap entries
	uint64_t idEntry; 

	fmi->h_fmi = (fmi_entry_t *) malloc(fmi->numEntries * sizeof(fmi_entry_t));
	if (fmi->h_fmi == NULL) return (31);

	for(idEntry = 0; idEntry < fmi->numEntries; ++idEntry){
		setCounterLayout(&fmi->h_countersFMI[idEntry * BITMAPS_PER_FMI], &fmi->h_fmi[idEntry], idEntry);
		setBitmapLayout(&fmi->h_bitmapBWT[idEntry * BITMAPS_PER_FMI], &fmi->h_fmi[idEntry]);
	}

	return(0);
}

char decodeFMIBase(uint32_t bit0, uint32_t bit1, uint32_t bit2)
{
	if(bit1 == 0 && bit0 == 0 && bit2 == 1) return ('A');
	if(bit1 == 0 && bit0 == 1 && bit2 == 1) return ('C');
	if(bit1 == 1 && bit0 == 0 && bit2 == 1) return ('G');
	if(bit1 == 1 && bit0 == 1 && bit2 == 1) return ('T');
	return ('N');
}

void printCountersPerEntry(fmi_buffer_t *fmi, uint64_t idEntry)
{
	const uint32_t FMI_COUNTERS_PER_ENTRY = FMI_NUM_COUNTERS / FMI_ALTERNATE_COUNTERS; // (4 letters) => 2 counters / fmi_entry
	      uint32_t idCounter;

	for(idCounter = 0; idCounter < FMI_COUNTERS_PER_ENTRY; ++idCounter){
		printf("%c=[%3lld] ", decodeBase((idEntry % FMI_COUNTERS_PER_ENTRY) * FMI_COUNTERS_PER_ENTRY + idCounter), fmi->h_fmi[idEntry].counters[idCounter]);
	}
}

uint32_t checkingFMI(fmi_buffer_t *fmi)
{
	const uint32_t LUT[12] = {3,7,11,0,1,2,4,5,6,8,9,10};
	const uint32_t FMI_NUM_BITMAPS = FMI_ENTRY_SIZE * BWT_CHAR_LENGTH / UINT32_LENGTH; // 12 words             (4 FMI entries x 3 bits)
	const uint32_t BITMAPS_PER_FMI = FMI_ENTRY_SIZE / UINT32_LENGTH;                   // 4 FMI bitmap entries
	const uint32_t PACKET_FMI_ENTRIES = FMI_NUM_BITMAPS / BITMAPS_PER_FMI;             // 3 FMI bitmap entries (12 bitmaps / 4 packets) 

	uint64_t idEntry;
	uint32_t idBitmap, idBase;
	
	printf("\n ");
	for(idEntry = 0; idEntry < fmi->numEntries; ++idEntry){
		printf("[%2lld] ", idEntry);
		printCountersPerEntry(fmi, idEntry);

		for(idBitmap = 0; idBitmap < BITMAPS_PER_FMI; ++idBitmap){
			uint32_t bitmap0 = fmi->h_fmi[idEntry].bitmaps[LUT[idBitmap * PACKET_FMI_ENTRIES]];
			uint32_t bitmap1 = fmi->h_fmi[idEntry].bitmaps[LUT[idBitmap * PACKET_FMI_ENTRIES + 1]];
			uint32_t bitmap2 = fmi->h_fmi[idEntry].bitmaps[LUT[idBitmap * PACKET_FMI_ENTRIES + 2]];

			for(idBase = 0; idBase < UINT32_LENGTH; ++idBase){
				const uint32_t bitmap_alignment = UINT32_LENGTH - idBase - 1;
				const uint32_t bit0 = (bitmap0 >> bitmap_alignment) & 1;
				const uint32_t bit1 = (bitmap1 >> bitmap_alignment) & 1;
				const uint32_t bit2 = (bitmap2 >> bitmap_alignment) & 1;

				printf("%c", decodeFMIBase(bit0,bit1,bit2));
			}
			printf(" | ");
		}
		printf("\n ");
	}

	return(0);
}


uint32_t loadBwtMFASTA(const char *fn, fmi_buffer_t *fmi)
{
	FILE *fp = NULL;
	char lineFile[FILE_SIZE_LINES], *tmp_reference;
	uint64_t sizeFile = 0, position = 0;
	int32_t charsRead = 0;

	fp = fopen(fn, "rb");
	if (fp == NULL) return (30);

	fseek(fp, 0L, SEEK_END);
	sizeFile = ftell(fp);
	rewind(fp);

	fmi->h_asciiBWT = (unsigned char*) malloc(sizeFile * sizeof(unsigned char));
	if (fmi->h_asciiBWT == NULL) return (31);

	while((!feof(fp)) && (fgets(lineFile, FILE_SIZE_LINES, fp) != NULL)){
		if (lineFile[0] != '>'){
			charsRead = strlen(lineFile);
			if(charsRead) charsRead--;
			memcpy((fmi->h_asciiBWT + position), lineFile, charsRead);
			position +=  charsRead;
		}
	}

	fmi->size       = position;
	fmi->numEntries = DIV_CEIL(fmi->size, FMI_ENTRY_SIZE) + 1;

	fclose(fp);
	return (0);
}

uint32_t saveFMI(const char *fn, fmi_buffer_t *fmi)
{
  char fileName[512];
  FILE *fp = NULL;

  sprintf(fileName, "%s.%llu.%u.fmi", fn, fmi->size, FMI_ENTRY_SIZE);

  fp = fopen(fileName, "wb");
  if (fp == NULL) return (8);

  fwrite(&fmi->numEntries, sizeof(uint64_t), 1, fp);
  fwrite(&fmi->size, sizeof(uint64_t), 1, fp);
  fwrite(fmi->h_fmi, sizeof(fmi_entry_t), fmi->numEntries, fp);

  printf("FMI SIZE: %lld, FMI ENTRIES: %lld \n", fmi->size, fmi->numEntries);

  fclose(fp);
  return (0);
}

uint32_t freeFMI(fmi_buffer_t **fmi)
{   
  if((* fmi)->h_fmi != NULL){
    free((* fmi)->h_fmi);
    (* fmi)->h_fmi = NULL;
  }

  if((* fmi)->h_countersFMI != NULL){
    free((* fmi)->h_countersFMI);
    (* fmi)->h_countersFMI = NULL;
  }

  if((* fmi)->h_bitmapBWT != NULL){
    free((* fmi)->h_bitmapBWT);
    (* fmi)->h_bitmapBWT = NULL;
  }

  if((* fmi)->h_asciiBWT != NULL){
    free((* fmi)->h_asciiBWT);
    (* fmi)->h_asciiBWT = NULL;
  }

  if((* fmi) != NULL){
    free(* fmi);
    (* fmi) = NULL;
  }

  return(0);
}

inline char *processError(uint32_t e){
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

uint32_t initFMI(fmi_buffer_t **fm_index)
{
  fmi_buffer_t *fmi  = (fmi_buffer_t *) malloc(sizeof(fmi_buffer_t));
	fmi->size          = 0;
	fmi->numEntries    = 0;

	fmi->h_fmi         = NULL;
	fmi->h_countersFMI = NULL;
	fmi->h_bitmapBWT   = NULL;
	fmi->h_asciiBWT    = NULL;

  (* fm_index) = fmi;
  return (0);
}

int32_t main(int argc, char *argv[])
{
  fmi_buffer_t *fmi;
  char *refFile = argv[1];
  int error;

	error = initFMI(&fmi);    
	CATCH_ERROR(error);

	printf("=> Loading MultiFasta BWT ... \n");
	error = loadBwtMFASTA(refFile, fmi);
  CATCH_ERROR(error);

	printf("=> BWT size: %llu, Number of entries: %llu\n", fmi->size, fmi->numEntries);
	printf("=> Transforming RAW BWT to bitmap BWT ... \n");
	error = transformBWTtoPEQ(fmi);
  CATCH_ERROR(error);

	printf("=> Building FMI counters ... \n");
	error = buildCounters(fmi);
  CATCH_ERROR(error);

	printf("=> Packing FMI structure ... \n");
 	error = buildFMI(fmi);
  CATCH_ERROR(error);

 // printf("=> Checking FMI structure ... \n");
 // error = checkingFMI(fmi);
 // CATCH_ERROR(error);

	printf("=> Saving FMI ... \n");
  error = saveFMI(refFile, fmi);
  CATCH_ERROR(error);
    
  printf("=> Done! \n");
  error = freeFMI(&fmi);
  CATCH_ERROR(error);

  return (0);
}

