/*
 * PROJECT: FM-index on GPU (backward-search and decode primitives)
 * FILE: fmi-interface.h
 * DATE: 1/10/2015
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 * DESCRIPTION: Interface for FMI on GPU
 */

#ifndef GPU_FMI_INTERFACE_H_
#define GPU_FMI_INTERFACE_H_

/*
 * Constants
 */
/* Defines related to BWT representation */
#define	GPU_FMI_BWT_CHAR_LENGTH     3       // 3 bits / character

#define GPU_FMI_NUM_COUNTERS        4       // 4 virtual counters
#define GPU_FMI_ALTERNATE_COUNTERS  2       // 2
#define GPU_FMI_ENTRY_SIZE          128     // 128 bases / FMI entry

#define GPU_FMI_COUNTERS_PER_ENTRY  (GPU_FMI_NUM_COUNTERS / GPU_FMI_ALTERNATE_COUNTERS)                 //   2 physical counters / FMI entry
#define GPU_FMI_BITMAPS_PER_ENTRY   (GPU_FMI_ENTRY_SIZE * GPU_FMI_BWT_CHAR_LENGTH / GPU_UINT32_LENGTH)  //   12 physical bitmaps / FMI entry

/*
 * Enum types for Device & Host
 */
typedef enum
{
  GPU_INDEX_NONE,
  GPU_INDEX_MFASTA_FILE,
  GPU_INDEX_PROFILE_FILE,
  GPU_INDEX_ASCII,
  GPU_INDEX_GEM_FULL,
  GPU_INDEX_GEM_FILE
} gpu_index_coding_t;

/*
 * Common types for Device & Host
 */
/* SA data structures */
typedef uint64_t  gpu_sa_decode_text_pos_t;
typedef uint64_t  gpu_sa_entry_t;

typedef struct {
  uint64_t*         sa;
  uint64_t          sa_sampling;
  uint64_t          sa_length;
  uint64_t          index_coding;
} gpu_gem_sa_dto_t;

typedef struct {
  gpu_sa_entry_t      *h_sa;
  uint64_t            numEntries;
  uint32_t            samplingRate;
  gpu_index_coding_t  indexCoding;
} gpu_sa_dto_t;

/* FMI data structures */
typedef struct {                                  // FMI Entry (64 Bytes) using:
  uint64_t counters[GPU_FMI_COUNTERS_PER_ENTRY];	// 4 letters: Alternate counters   (2  uint64_t)
  uint32_t bitmaps[GPU_FMI_BITMAPS_PER_ENTRY];		// 5 letters: 12 Bitmaps x 32 bits (3 uint128_t)
} gpu_fmi_entry_t;

typedef uint64_t	gpu_fmi_decode_init_pos_t;
typedef uint64_t  gpu_fmi_decode_text_pos_t;

typedef struct{
  uint64_t low;
  uint64_t hi;
} gpu_sa_search_inter_t;

typedef struct{
  uint64_t hi;
  uint64_t low;
} gpu_fmi_search_seed_t;

typedef struct{
  uint64_t interval;
  uint64_t steps;
} gpu_fmi_decode_end_pos_t;

typedef struct {
  uint64_t*          c;               // Occurrences of each character
  uint64_t*          C;               // The cumulative occurrences ("ranks") of the symbols of the string
  uint64_t*          mayor_counters;  // Pointer to the Mayor Counters (Rank)
  uint64_t*          bwt_mem;         // Pointer to the BWT structure in memory
  uint64_t           bwt_length;      // Length of the BWT
  gpu_index_coding_t index_coding;
} gpu_gem_fmi_dto_t;

typedef struct {
  gpu_fmi_entry_t     *h_fmi;
  uint64_t            bwtSize;
  gpu_index_coding_t  indexCoding;
} gpu_fmi_dto_t;

/*
 * Index General Structures
 */
typedef struct {
  char          *filename;
  gpu_fmi_dto_t fmi;
  gpu_sa_dto_t  sa;
} gpu_index_dto_t;

/*
 * Obtain Buffers
 */
gpu_fmi_search_seed_t* 	   gpu_fmi_search_buffer_get_seeds_(const void* const fmiBuffer);
gpu_sa_search_inter_t*     gpu_fmi_search_buffer_get_sa_intervals_(const void* const fmiBuffer);
gpu_fmi_decode_init_pos_t* gpu_fmi_decode_buffer_get_init_pos_(const void* const fmiBuffer);
gpu_fmi_decode_end_pos_t*  gpu_fmi_decode_buffer_get_end_pos_(const void* const fmiBuffer);
gpu_sa_decode_text_pos_t*  gpu_sa_decode_buffer_get_ref_pos_(const void* const fmiBuffer);

/*
 * Get elements
 */
uint32_t gpu_fmi_search_buffer_get_max_seeds_(const void* const fmiBuffer);
uint32_t gpu_fmi_decode_buffer_get_max_positions_(const void* const fmiBuffer);

/*
 * Main functions
 */
void gpu_fmi_search_init_buffer_(void* const fmiBuffer);
void gpu_fmi_search_send_buffer_(void* const fmiBuffer, const uint32_t numSeeds);
void gpu_fmi_search_receive_buffer_(const void* const fmiBuffer);
void gpu_fmi_search_init_and_realloc_buffer_(void* const fmiBuffer, const uint32_t numSeeds);


void gpu_fmi_decode_init_buffer_(void* const fmiBuffer);
void gpu_fmi_decode_send_buffer_(void* const fmiBuffer, const uint32_t numDecodings, const uint32_t samplingRate);
void gpu_fmi_decode_receive_buffer_(const void* const fmiBuffer);
void gpu_fmi_decode_init_and_realloc_buffer_(void* const fmiBuffer, const uint32_t numDecodings);

#endif /* GPU_FMI_INTERFACE_H_ */

