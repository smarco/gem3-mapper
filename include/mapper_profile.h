/*
 * PROJECT: GEMMapper
 * FILE: mapper_profile.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Basically, these counters/timers focus on the measurement of:
 *     - Time spent in certain functions/blocks-code (and number of calls/executions)
 *     - Functional profile of execution flow branches (counting)
 */

#ifndef MAPPER_PROFILE_H_
#define MAPPER_PROFILE_H_

#include "essentials.h"

/*
 * Mapper
 */
#define GP_MAPPER_ALL             0
#define GP_MAPPER_LOAD_INDEX      1
#define GP_OUTPUT_MAP_SE          2
#define GP_OUTPUT_SAM_SE          3

#define GP_MAPPER_CUDA                4
#define GP_MAPPER_CUDA_PUT_BPM_BUFFER 5

#define GP_INPUT_FASTA_PARSE_SEQUENCE 6

/*
 * Archive Search
 */
#define GP_ARCHIVE_SEARCH_SE                  10
#define GP_ARCHIVE_SELECT_MATCHES             11
#define GP_ARCHIVE_SEARCH_GENERATE_CANDIDATES 12
#define GP_ARCHIVE_SEARCH_COPY_CANDIDATES     13
#define GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES 14

/*
 * I/O
 */
#define GP_INPUT_FILL_BUFFER                            20
#define GP_BUFFERED_INPUT_RELOAD                        21
#define GP_BUFFERED_INPUT_BUFFER_SIZE                   22
#define GP_OUTPUT_WRITE_BUFFER                          23
#define GP_OUTPUT_BYTES_WRITTEN                         24
#define GP_OUTPUT_BUFFER_EXTENSIONS                     25
#define GP_OUTPUT_BUFFER_REQUESTS                       26
#define GP_OUTPUT_BUFFER_REQUESTS_STALLS                27
#define GP_OUTPUT_BUFFER_REQUESTS_STALLS_BUSY           28
#define GP_OUTPUT_BUFFER_REQUESTS_STALLS_NOT_PRIORITY   29

/*
 * Approximate search
 */
#define GP_AS_MAIN             30
#define GP_AS_READ_RECOVERY    31
#define GP_AS_EXACT_SEARCH     32

#define GP_AS_ADAPTIVE_MCS     40

#define GP_AS_FILTER_REGIONS                       60
#define GP_AS_FILTER_REGIONS_NUM_ELEGIBLE_REGIONS  61
#define GP_AS_FILTER_REGIONS_SEARCH_D2             62
#define GP_AS_FILTER_REGIONS_SEARCH_D2_HIT         63
#define GP_AS_FILTER_REGIONS_SEARCH_D1             64
#define GP_AS_FILTER_REGIONS_SEARCH_D1_HIT         65
#define GP_AS_FILTER_REGIONS_SEARCH_D0             66
#define GP_AS_FILTER_REGIONS_PROCESSED             67
#define GP_AS_FILTER_REGIONS_SKIPPED               68

/*
 * Region Profile
 */
#define GP_REGION_PROFILE_SOFT          80
#define GP_REGION_PROFILE_ADAPTIVE      81
#define GP_REGION_PROFILE_QUIT_PROFILE  82

#define GP_REGION_PROFILE_SOFT_NUM_REGIONS       83
#define GP_REGION_PROFILE_SOFT_REGION_LENGTH     84
#define GP_REGION_PROFILE_SOFT_REGION_CANDIDATES 85
#define GP_REGION_PROFILE_SOFT_TOTAL_CANDIDATES  86

/*
 * Filtering Candidates (Verifying)
 */
#define GP_FC_VERIFY                     100
#define GP_FC_DECODE_POSITIONS           101
#define GP_FC_VERIFY_CANDIDATE_REGIONS   102
#define GP_FC_REALIGN_CANDIDATE_REGIONS  103
#define GP_FC_RETRIEVE_CANDIDATE_REGIONS 104

#define GP_FC_NUM_CANDIDATE_POSITIONS         110
#define GP_FC_NUM_CANDIDATE_REGIONS           111
#define GP_FC_NUM_ACCEPTED_REGIONS            117
#define GP_FC_CANDIDATE_REGIONS_COVERAGE      112
#define GP_FC_CANDIDATE_REGIONS_EXT_COVERAGE  113
#define GP_FC_LEVENSHTEIN_BPM                 114
#define GP_FC_LEVENSHTEIN_ACCEPTED            115
#define GP_FC_LEVENSHTEIN_BPM_QUICK_ABANDON   116

#define GP_FC_KMER_COUNTER_FILTER             117
#define GP_FC_KMER_COUNTER_FILTER_DISCARDED   118

/*
 * FM-Index
 */
#define GP_FMIDX_LOOKUP_DIST    120

/*
 * Archive search-group Dispatcher
 */
#define GP_SGDISPATCHER_REQUESTS_GENERATING                            140
#define GP_SGDISPATCHER_REQUESTS_GENERATING_EXTENSION                  141
#define GP_SGDISPATCHER_REQUESTS_GENERATING_EXTENSION_VOLUNTEER        142
#define GP_SGDISPATCHER_REQUESTS_GENERATING_STALLS                     143
#define GP_SGDISPATCHER_REQUESTS_GENERATING_STALLS_BUSY                144
#define GP_SGDISPATCHER_REQUESTS_GENERATING_STALLS_NOT_PRIORITY        145
#define GP_SGDISPATCHER_REQUESTS_SELECTING                             146
#define GP_SGDISPATCHER_REQUESTS_SELECTING_EXTENSION                   147
#define GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS                      148
#define GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS_IDLE                 149
#define GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS_NO_SINGLE_GROUPS     150
#define GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS_EXTENSION_NOT_READY  151

/*
 * BPM-GPU
 */
#define GP_BPM_GPU_BUFFER_SEND              160
#define GP_BPM_GPU_BUFFER_USAGE_CANDIDATES  161
#define GP_BPM_GPU_BUFFER_USAGE_QUERIES     162
#define GP_BPM_GPU_BUFFER_USAGE_PEQ_ENTRIES 163
#define GP_BPM_GPU_BUFFER_CHECK_TIME        164

/*
 * System
 */
GEM_INLINE void mapper_profile_print_io(FILE* const stream);
GEM_INLINE void mapper_profile_print_mem_structs(FILE* const stream);
/*
 * Global Mapper
 */
GEM_INLINE void mapper_profile_print_mapper_adaptive(FILE* const stream);
GEM_INLINE void mapper_profile_print_mapper_adaptive_ranks(FILE* const stream);
GEM_INLINE void mapper_profile_print_mapper_cuda_adaptive(FILE* const stream);

GEM_INLINE void mapper_profile_print_mapper_efficiency_ratios(FILE* const stream);
/*
 * Approximate string search
 */
GEM_INLINE void mapper_profile_print_approximate_search(FILE* const stream);
GEM_INLINE void mapper_profile_print_approximate_search_ranks(FILE* const stream);
/*
 * Region Profile
 */
GEM_INLINE void mapper_profile_print_region_profile_soft(FILE* const stream);
/*
 * Filtering Generating
 */
GEM_INLINE void mapper_profile_print_filtering_generating(FILE* const stream);
GEM_INLINE void mapper_profile_print_filtering_generating_ranks(FILE* const stream);
/*
 * Filtering Verification
 */
GEM_INLINE void mapper_profile_print_filtering_verifying(FILE* const stream,const bool verification_profile);
GEM_INLINE void mapper_profile_print_filtering_verifying_ranks(FILE* const stream);
/*
 * Neighborhood Search
 */
GEM_INLINE void mapper_profile_print_neighborhood_search(FILE* const stream);
GEM_INLINE void mapper_profile_print_neighborhood_search_ranks(FILE* const stream);
/*
 * Archive Search
 */
GEM_INLINE void mapper_profile_print_archive_search(FILE* const stream);
/*
 * Archive Select
 */
GEM_INLINE void mapper_profile_print_archive_select(FILE* const stream);
/*
 * Archive Search-Group (Dispatcher, BMP-Buffers, ...)
 */
GEM_INLINE void mapper_profile_print_archive_search_group(FILE* const stream);

/*
 * Error Messages
 */
//#define GEM_ERROR_MAPPER_PROFILE ""

#endif /* MAPPER_PROFILE_H_ */
