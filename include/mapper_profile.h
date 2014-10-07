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
#define GP_ARCHIVE_SEARCH_SE      2
#define GP_ARCHIVE_SELECT_MATCHES 3
#define GP_OUTPUT_MAP_SE          4
#define GP_OUTPUT_SAM_SE          6

/*
 * I/O
 */
#define GP_INPUT_FILL_BUFFER                            10
#define GP_BUFFERED_INPUT_RELOAD                        11
#define GP_BUFFERED_INPUT_BUFFER_SIZE                   12
#define GP_OUTPUT_WRITE_BUFFER                          13
#define GP_OUTPUT_BYTES_WRITTEN                         14
#define GP_OUTPUT_BUFFER_EXTENSIONS                     15
#define GP_OUTPUT_BUFFER_REQUESTS                       16
#define GP_OUTPUT_BUFFER_REQUESTS_STALLS                17
#define GP_OUTPUT_BUFFER_REQUESTS_STALLS_BUSY           18
#define GP_OUTPUT_BUFFER_REQUESTS_STALLS_NOT_PRIORITY   19

/*
 * Approximate search
 */
#define GP_AS_MAIN             30
#define GP_AS_READ_RECOVERY    31
#define GP_AS_EXACT_SEARCH     32

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

/*
 * Filtering Candidates (Verifying)
 */
#define GP_FC_VERIFY  100
#define GP_FC_DECODE  101
#define GP_FC_CHECK   102

/*
 * FM-Index
 */
#define GP_FMIDX_LOOKUP_DIST    120

/*
 * Archive search-group Dispatcher
 */
#define GP_SGDISPATCHER_REQUESTS_GENERATING                            140
#define GP_SGDISPATCHER_REQUESTS_GENERATING_EXTENSION                  141
#define GP_SGDISPATCHER_REQUESTS_GENERATING_STALLS                     142
#define GP_SGDISPATCHER_REQUESTS_GENERATING_STALLS_BUSY                143
#define GP_SGDISPATCHER_REQUESTS_GENERATING_STALLS_NOT_PRIORITY        144
#define GP_SGDISPATCHER_REQUESTS_SELECTING                             145
#define GP_SGDISPATCHER_REQUESTS_SELECTING_EXTENSION                   146
#define GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS                      147
#define GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS_IDLE                 148
#define GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS_NO_SINGLE_GROUPS     149
#define GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS_EXTENSION_NOT_READY  150

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
GEM_INLINE void mapper_profile_print_system_info(FILE* const stream);
GEM_INLINE void mapper_profile_print_io(FILE* const stream);
/*
 * Global Mapper
 */
GEM_INLINE void mapper_profile_print_mapper_adaptive(FILE* const stream);
GEM_INLINE void mapper_profile_print_mapper_efficiency_ratios(FILE* const stream);
/*
 * Approximate string search
 */
GEM_INLINE void mapper_profile_print_approximate_search(FILE* const stream);
GEM_INLINE void mapper_profile_print_approximate_search_ranks(FILE* const stream);
/*
 * Filtering Generating
 */
GEM_INLINE void mapper_profile_print_filtering_generating(FILE* const stream);
GEM_INLINE void mapper_profile_print_filtering_generating_ranks(FILE* const stream);
/*
 * Filtering Verification
 */
GEM_INLINE void mapper_profile_print_filtering_verifying(FILE* const stream);
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
