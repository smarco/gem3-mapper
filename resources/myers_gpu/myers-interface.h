/*
 * PROJECT: GEMMapper
 * FILE: gem-mapper.c
 * DATE:5/12/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Genomic Read Mapper
 */

#include <stdint.h>

/*
 * Constants
 */
#define PEQ_ALPHABET_SIZE 5
#define PEQ_ENTRY_LENGTH  128
#define UINT32_LENGTH     32
#define PEQ_SUBENTRIES    (PEQ_ENTRY_LENGTH / UINT32_LENGTH)
#define UINT32_ONE_MASK   0x00000001u

/*
 * Enum types for Device & Host
 */
typedef enum
{
 MFASTA_FILE,
 PROFILE_REFERENCE_FILE,
 ASCII,
 GEM
} _reference_coding;
typedef _reference_coding refCoding_t;
typedef enum
{
  ARCH_TESLA     = UINT32_ONE_MASK << 0,
  ARCH_FERMI_1G  = UINT32_ONE_MASK << 1,
  ARCH_FERMI_2G  = UINT32_ONE_MASK << 2,
  ARCH_KEPLER_1G = UINT32_ONE_MASK << 3,
 ARCH_KEPLER_2G  = UINT32_ONE_MASK << 4,
 ARCH_MAXWELL    = UINT32_ONE_MASK << 5,
 ARCH_NEWGEN     = UINT32_ONE_MASK << 31,
 ARCH_SUPPORTED  = ARCH_FERMI_1G | ARCH_FERMI_2G | ARCH_KEPLER_1G | ARCH_KEPLER_2G | ARCH_MAXWELL | ARCH_NEWGEN
} _dev_arch;
typedef _dev_arch devArch_t;

/*
 * Common types for Device & Host
 */
typedef struct { /* each row 1 PEQ Entry (128bits) */
 uint32_t bitmap[PEQ_ALPHABET_SIZE][PEQ_SUBENTRIES];
} qryEntry_t;
typedef struct {
 uint32_t column;
 uint32_t score;
} resEntry_t;
typedef struct {
 uint64_t position;
 uint32_t query;
 uint32_t size;
} candInfo_t;
typedef struct {
 uint32_t posEntry;
 uint32_t size;
} qryInfo_t;

/*
 * Obtain Buffers
 */
inline qryEntry_t* getPEQBuffer(void* myersBuffer);
inline candInfo_t* getCandidatesBuffer(void* myersBuffer);
inline qryInfo_t* getPEQInfoBuffer(void* myersBuffer);
inline resEntry_t* getResultsBuffer(void* myersBuffer);

/*
 * Get elements
 */
inline uint32_t getMaxPEQEntries(void* myersBuffer);
inline uint32_t getMaxCandidates(void* myersBuffer);
inline uint32_t getMaxQueries(void* myersBuffer);
inline uint32_t getIdDeviceBuffer(void* myersBuffer);

/*
 * Main functions
 */
void initMyers(void*** myersBuffer,uint32_t numBuffers,uint32_t maxMbPerBuffer,
    const char* referenceRaw,refCoding_t refCoding,const uint64_t refSize,
    uint32_t averageQuerySize,uint32_t candidatesPerQuery,devArch_t selectedArchitectures);
void sendMyersBuffer(void* myersBuffer,uint32_t numPEQEntries,uint32_t numQueries,uint32_t numCandidates);
void receiveMyersBuffer(void* myersBuffer);
void endMyers(void*** myersBuffer);

