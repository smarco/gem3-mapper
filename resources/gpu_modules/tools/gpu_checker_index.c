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

#define	BWT_CHAR_LENGTH				  3
#define FMI_NUM_COUNTERS			  4
#define FMI_ALTERNATE_COUNTERS	2
#define FMI_ENTRY_SIZE				  128

#define	UINT32_LENGTH				    32
#define	UINT64_LENGTH				    64

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
	uint64_t        bwtSize;
	uint64_t        numEntries;
	fmi_entry_t    *h_fmi;
} fmi_t;


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

inline char decodeFMIBase(const uint32_t bit0, const uint32_t bit1, const uint32_t bit2)
{
	if(bit1 == 0 && bit0 == 0 && bit2 == 1) return ('A');
	if(bit1 == 0 && bit0 == 1 && bit2 == 1) return ('C');
	if(bit1 == 1 && bit0 == 0 && bit2 == 1) return ('G');
	if(bit1 == 1 && bit0 == 1 && bit2 == 1) return ('T');
	return ('N');
}

inline void printCountersPerEntry(const fmi_t *fmi, uint64_t idEntry)
{
	const uint32_t FMI_COUNTERS_PER_ENTRY = FMI_NUM_COUNTERS / FMI_ALTERNATE_COUNTERS; // (4 letters) => 2 counters / fmi_entry
	      uint32_t idCounter;

	for(idCounter = 0; idCounter < FMI_COUNTERS_PER_ENTRY; ++idCounter){
		printf("%c=[%3lld] ", decodeBase((idEntry % FMI_COUNTERS_PER_ENTRY) * FMI_COUNTERS_PER_ENTRY + idCounter), fmi->h_fmi[idEntry].counters[idCounter]);
	}
}

/*
static char *binrep(uint32_t val, char *buff, int32_t sz)
{
  char *pbuff = buff;

  if (sz < 1) return NULL;
  if (val == 0) {
    *pbuff++ = '0';
    *pbuff = '\0';
    return buff;
  }

  pbuff += sz;
  *pbuff-- = '\0';
  while (val != 0) {
    if (sz-- == 0) return NULL;
    *pbuff-- = ((val & 1) == 1) ? '1' : '0';
    val >>= 1;
  }
  return pbuff+1;
}*/


static char *binrep(uint32_t val, char *buff, int32_t sz)
{
	int32_t idBit;
  char *pbuff = buff;

  pbuff   += sz;
  *pbuff-- = '\0';
  for(idBit = 0; idBit < sz; ++idBit){
    *pbuff-- = ((val & 1) == 1) ? '1' : '0';
    val >>= 1;
  }
	return pbuff+1;
}

inline void printBinary(const uint32_t num)
{
	char buffer [UINT32_LENGTH + 1], *pbuff = &buffer[0];
	pbuff = binrep (num, pbuff, UINT32_LENGTH);
  printf ("%s",pbuff);
}


inline uint32_t printEntryFMI(char* TAG, const fmi_t *fmi, const uint64_t idEntry)
{
	const uint32_t LUT[12] = {3,7,11,0,1,2,4,5,6,8,9,10};
	const uint32_t FMI_NUM_BITMAPS    = FMI_ENTRY_SIZE * BWT_CHAR_LENGTH / UINT32_LENGTH; // 12 words             (4 FMI entries x 3 bits)
	const uint32_t BITMAPS_PER_FMI    = FMI_ENTRY_SIZE / UINT32_LENGTH;                   // 4 FMI bitmap entries
	const uint32_t PACKET_FMI_ENTRIES = FMI_NUM_BITMAPS / BITMAPS_PER_FMI;                // 3 FMI bitmap entries (12 bitmaps / 4 packets) 
	const uint32_t enableRawBitmaps   = 1, enableBinary = 1;

	uint32_t idBitmap, idBase;

	printf("%s, [%2lld] ", TAG, idEntry);
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
	printf("\n");

	if(enableRawBitmaps){
		uint32_t idRow;
		for(idRow = 0; idRow < BWT_CHAR_LENGTH; ++idRow){
			printf("\t\t\t BITMAP %d :", idRow);
			for(idBitmap = 0; idBitmap < BITMAPS_PER_FMI; ++idBitmap){
				const uint32_t bitmap = fmi->h_fmi[idEntry].bitmaps[LUT[idBitmap * PACKET_FMI_ENTRIES + idRow]];
				if(enableBinary){
					printBinary(bitmap);
					printf(" | ");
				}else{
					printf("%8x | ", bitmap);
				}
			}
			printf("\n");
		}
	}

	return(0);
}

inline uint32_t initFMI(fmi_t **fm_index)
{
  fmi_t *fmi  = (fmi_t *) malloc(sizeof(fmi_t));
	fmi->bwtSize       = 0;
	fmi->numEntries    = 0;
	fmi->h_fmi         = NULL;

  (* fm_index) = fmi;
  return (0);
}

inline uint32_t loadFMI(const char *fn, fmi_t *fmi)
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

inline uint32_t freeFMI(fmi_t **fmi)
{   
  if((* fmi)->h_fmi != NULL){
    free((* fmi)->h_fmi);
    (* fmi)->h_fmi = NULL;
  }

  if((* fmi) != NULL){
    free(* fmi);
    (* fmi) = NULL;
  }
  return(0);
}

inline uint32_t equalEntry(const fmi_t *fmi_A, const fmi_t *fmi_B, const uint64_t idEntry)
{
	const uint32_t FMI_COUNTERS_PER_ENTRY = FMI_NUM_COUNTERS / FMI_ALTERNATE_COUNTERS; // (4 letters) => 2 counters / fmi_entry
	const uint32_t FMI_BITMAPS_PER_ENTRY  = FMI_ENTRY_SIZE * BWT_CHAR_LENGTH / UINT32_LENGTH;
	uint32_t sameEntry = 1, idBitmap, idCounter;

	for(idCounter = 0; idCounter < FMI_COUNTERS_PER_ENTRY; ++idCounter)
		if(fmi_A->h_fmi[idEntry].counters[idCounter] != fmi_B->h_fmi[idEntry].counters[idCounter]) sameEntry = 0;

	for(idBitmap = 0; idBitmap < FMI_BITMAPS_PER_ENTRY; ++idBitmap)
		if(fmi_A->h_fmi[idEntry].bitmaps[idBitmap] != fmi_B->h_fmi[idEntry].bitmaps[idBitmap]) sameEntry = 0;

	return(sameEntry);
}

inline uint32_t checkingFMI(const fmi_t *fmi_A, const fmi_t *fmi_B)
{
	const uint64_t maxEntries	    = 10; // Just check and print the first results
	const uint64_t numEntries	    = MIN(fmi_A->numEntries, fmi_B->numEntries);
	uint64_t idEntry, missMatches = 0;

  for(idEntry = 0; idEntry < numEntries; ++idEntry){
    if(!equalEntry(fmi_A, fmi_B, idEntry)){
      missMatches++;
      if(missMatches < maxEntries){
        CATCH_ERROR(printEntryFMI("FMI_A", fmi_A, idEntry));
        CATCH_ERROR(printEntryFMI("FMI_B", fmi_B, idEntry));
        printf("\n");
      }
    }
  }

	printf("=> Divergent FMI entries %d (%.2f \%) \n", missMatches, ((float)missMatches/(float)numEntries)*100.0);
	return(0);
}

int32_t main(int argc, char *argv[])
{
  fmi_t *fmIndex_A     = NULL;
  fmi_t *fmIndex_B     = NULL;

  char *fmiFile_A     = argv[1];
  char *fmiFile_B     = argv[2];
  int error;

	error = initFMI(&fmIndex_A);    
	CATCH_ERROR(error);
	error = initFMI(&fmIndex_B);    
	CATCH_ERROR(error);

	printf("=> Loading FMI_A: %s \n", fmiFile_A);
	CATCH_ERROR(loadFMI(fmiFile_A, fmIndex_A));
	printf("\t BWT size: %llu, Number of entries: %llu\n", fmIndex_A->bwtSize, fmIndex_A->numEntries);
	printf("=> Loading FMI_B: %s \n", fmiFile_B);
	CATCH_ERROR(loadFMI(fmiFile_B, fmIndex_B));
	printf("\t BWT size: %llu, Number of entries: %llu\n", fmIndex_B->bwtSize, fmIndex_B->numEntries);

 	printf("=> Checking FMI structure ... \n");
 	error = checkingFMI(fmIndex_A, fmIndex_B);
 	CATCH_ERROR(error);

  printf("=> Done! \n");
  error = freeFMI(&fmIndex_A);
  CATCH_ERROR(error);
  error = freeFMI(&fmIndex_B);
  CATCH_ERROR(error);

  return (0);
}

