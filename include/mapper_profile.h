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
#define GP_MAPPER_ALL                                         0
#define GP_MAPPER_LOAD_INDEX                                  1
#define GP_OUTPUT_MAP_SE                                      2
#define GP_OUTPUT_MAP_PE                                      3
#define GP_OUTPUT_SAM_SE                                      4
#define GP_OUTPUT_SAM_PE                                      5

#define GP_MAPPER_CUDA_THREAD                                 6
#define GP_MAPPER_CUDA_THREAD_GENERATING                      7
#define GP_MAPPER_CUDA_THREAD_SELECTING                       8
#define GP_MAPPER_CUDA_THREAD_RESTART_UNFIT                   9

#define GP_MAPPER_NUM_READS                                  10

/*
 * Archive Search
 */
#define GP_ARCHIVE_SEARCH_SE                                110
#define GP_ARCHIVE_SEARCH_SE_INIT                           111
#define GP_ARCHIVE_SEARCH_PREPARE_SEQUENCE                  112
#define GP_ARCHIVE_SEARCH_SE_GENERATE_CANDIDATES            113
#define GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES              114
#define GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH                  115
#define GP_ARCHIVE_SELECT_SE_MATCHES                        116
#define GP_ARCHIVE_SCORE_SE_MATCHES                         117

// Copy/Retrieve candidates
#define GP_ARCHIVE_SEARCH_COPY_CANDIDATES                   118
#define GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES               119

// Search Group
#define GP_ARCHIVE_SEARCH_GROUP_BUFFERS_USED                120
#define GP_ARCHIVE_SEARCH_GROUP_RETRIEVE_CANDIDATES_DELAY   121

// Paired Mode
#define GP_ARCHIVE_SEARCH_PE                                122

#define GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT             123
#define GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT_SUCCESS     124
#define GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY             125
#define GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES              126
#define GP_ARCHIVE_SEARCH_PE_GENERATE_CANDIDATES            139
#define GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH                  140

#define GP_ARCHIVE_SELECT_PE_MATCHES                        142
#define GP_ARCHIVE_SCORE_PE_MATCHES                         143
#define GP_PAIRED_MATCHES_FIND_PAIRS                        144

/*
 * I/O
 */
#define GP_INPUT_FASTA_PARSE_SEQUENCE                       145
#define GP_INPUT_FILL_BUFFER                                146
#define GP_BUFFERED_INPUT_RELOAD                            147
#define GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED             148
#define GP_BUFFERED_INPUT_BUFFER_SIZE                       149
#define GP_BUFFERED_OUTPUT_DUMP                             150
#define GP_OUTPUT_WRITE_BUFFER                              151
#define GP_OUTPUT_BYTES_WRITTEN                             152
#define GP_OUTPUT_BUFFER_EXTENSIONS                         153
#define GP_OUTPUT_BUFFER_REQUESTS                           154
#define GP_OUTPUT_BUFFER_REQUESTS_STALLS                    155
#define GP_OUTPUT_BUFFER_REQUESTS_STALLS_BUSY               156
#define GP_OUTPUT_BUFFER_REQUESTS_STALLS_NOT_PRIORITY       157

/*
 * Approximate search
 */
#define GP_AS_MAIN                                          160
#define GP_AS_READ_RECOVERY                                 162
#define GP_AS_EXACT_SEARCH                                  163
#define GP_AS_NEIGHBORHOOD_SEARCH                           164

#define GP_AS_ADAPTIVE_MCS                                  170

#define GP_AS_FILTERING_EXACT                               171
#define GP_AS_FILTERING_EXACT_BOOST                         172
#define GP_AS_FILTERING_INEXACT                             173
#define GP_AS_FILTERING_UNBOUNDED_ALIGN                     174

#define GP_AS_FILTERING_EXACT_MAPPED                        175
#define GP_AS_FILTERING_EXACT_BOOST_MAPPED                  176
#define GP_AS_FILTERING_INEXACT_MAPPED                      177
#define GP_AS_FILTERING_UNBOUNDED_ALIGN_MAPPED              178

#define GP_AS_FILTERING_EXACT_MCS                           179
#define GP_AS_FILTERING_INEXACT_MCS                         181

#define GP_AS_GENERATE_CANDIDATES                           185
#define GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS      186
#define GP_AS_GENERATE_CANDIDATES_SEARCH_D2                 187
#define GP_AS_GENERATE_CANDIDATES_SEARCH_D2_HIT             188
#define GP_AS_GENERATE_CANDIDATES_SEARCH_D2_HIT_CANDIDATES  189
#define GP_AS_GENERATE_CANDIDATES_SEARCH_D1                 190
#define GP_AS_GENERATE_CANDIDATES_SEARCH_D1_HIT             191
#define GP_AS_GENERATE_CANDIDATES_SEARCH_D1_HIT_CANDIDATES  192
#define GP_AS_GENERATE_CANDIDATES_SEARCH_D0_HIT             193
#define GP_AS_GENERATE_CANDIDATES_SEARCH_D0_HIT_CANDIDATES  194
#define GP_AS_GENERATE_CANDIDATES_PROCESSED                 195
#define GP_AS_GENERATE_CANDIDATES_SKIPPED                   196
#define GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING         197

/*
 * Region Profile
 */
#define GP_REGION_PROFILE_ADAPTIVE                          200
#define GP_REGION_PROFILE_QUIT_PROFILE                      201

#define GP_REGION_PROFILE_LIGHTWEIGHT                       205
#define GP_REGION_PROFILE_LIGHTWEIGHT_NUM_REGIONS           206
#define GP_REGION_PROFILE_LIGHTWEIGHT_NUM_REGIONS_UNIQUE    207
#define GP_REGION_PROFILE_LIGHTWEIGHT_NUM_REGIONS_STANDARD  208
#define GP_REGION_PROFILE_LIGHTWEIGHT_REGION_LENGTH         209
#define GP_REGION_PROFILE_LIGHTWEIGHT_REGION_CANDIDATES     210
#define GP_REGION_PROFILE_LIGHTWEIGHT_TOTAL_CANDIDATES      211

#define GP_REGION_PROFILE_BOOST                             213
#define GP_REGION_PROFILE_BOOST_NUM_REGIONS                 214
#define GP_REGION_PROFILE_BOOST_NUM_REGIONS_UNIQUE          215
#define GP_REGION_PROFILE_BOOST_NUM_REGIONS_STANDARD        216
#define GP_REGION_PROFILE_BOOST_REGION_LENGTH               217
#define GP_REGION_PROFILE_BOOST_REGION_CANDIDATES           218
#define GP_REGION_PROFILE_BOOST_TOTAL_CANDIDATES            219

#define GP_REGION_PROFILE_DELIMIT                           221
#define GP_REGION_PROFILE_DELIMIT_NUM_REGIONS               222
#define GP_REGION_PROFILE_DELIMIT_NUM_REGIONS_UNIQUE        223
#define GP_REGION_PROFILE_DELIMIT_NUM_REGIONS_STANDARD      224
#define GP_REGION_PROFILE_DELIMIT_REGION_LENGTH             225
#define GP_REGION_PROFILE_DELIMIT_REGION_CANDIDATES         226
#define GP_REGION_PROFILE_DELIMIT_TOTAL_CANDIDATES          227

/*
 * Filtering Candidates (Verifying)
 */
#define GP_FC_VERIFICATION                                  330
#define GP_FC_PROCESS_CANDIDATES                            331
#define GP_FC_DECODE_POSITIONS                              332
#define GP_FC_VERIFY_CANDIDATE_REGIONS                      333
#define GP_FC_RETRIEVE_BPM_BUFFER_CANDIDATE_REGIONS         334
#define GP_FC_REALIGN_CANDIDATE_REGIONS                     335
#define GP_FC_UNBOUNDED_ALIGNMENT                           336
#define GP_FC_REALIGN_BPM_BUFFER_CANDIDATE_REGIONS          337
#define GP_FC_RETRIEVE_CANDIDATE_REGIONS                    338
#define GP_FC_COMPOSE_REGIONS                               339

#define GP_FC_KMER_COUNTER_FILTER                           340
#define GP_FC_KMER_COUNTER_FILTER_DISCARDED                 341

#define GP_FC_EXTEND_MATCH                                  350
#define GP_FC_EXTEND_RETRIEVE_CANDIDATE_REGIONS             351
#define GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS               352
#define GP_FC_EXTEND_VERIFY_CANDIDATE_LENGTH                353
#define GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS              354

#define GP_CANDIDATE_POSITIONS                              400
#define GP_CANDIDATE_REGIONS                                402
#define GP_CANDIDATE_REGIONS_DUPLICATED                     403
#define GP_CANDIDATE_REGION_LENGTH                          404
#define GP_KMER_COUNTER_FILTER                              405
#define GP_KMER_COUNTER_FILTER_DISCARDED                    406
#define GP_LEVENSHTEIN_ACCEPTED                             407

#define GP_ACCEPTED_REGIONS                                 425
#define GP_ACCEPTED_EXACT                                   426
#define GP_ACCEPTED_INEXACT                                 427
#define GP_ACCEPTED_REGIONS_LENGTH                          428

#define GP_MATCHING_REGIONS_CHAIN                           450
#define GP_MATCHING_REGIONS_CHAIN_SUCCESS                   451
#define GP_MATCHING_REGIONS_CHAIN_COVERAGE                  452
#define GP_MATCHING_REGIONS_EXTEND                          453
#define GP_MATCHING_REGIONS_EXTEND_COVERAGE                 454
#define GP_MATCHING_REGIONS_EDIT_SCAFFOLD                   455
#define GP_MATCHING_REGIONS_EDIT_SCAFFOLDED                 456
#define GP_MATCHING_REGIONS_EDIT_SCAFFOLD_COVERAGE          457
#define GP_MATCHING_REGIONS_SWG_SCAFFOLD                    458
#define GP_MATCHING_REGIONS_SWG_SCAFFOLDED                  459

#define GP_BPM_TILED                                        470
#define GP_BMP_TILED_NUM_TILES                              471
#define GP_BMP_TILED_NUM_TILES_VERIFIED                     472
#define GP_BPM_QUICK_ABANDON                                473
#define GP_BPM_ALL                                          474
#define GP_BPM_ALL_QUICK_ABANDON                            475
#define GP_BPM_ALL_MATCHES_FOUND                            476

/*
 * BPM-GPU
 */
#define GP_BPM_GPU_INIT                                     501
#define GP_BPM_GPU_BUFFER_INIT                              502
#define GP_BPM_GPU_BUFFER_SEND                              503
#define GP_BPM_GPU_BUFFER_CHECK_TIME                        504
#define GP_BPM_GPU_BUFFER_NUM_CANDIDATES                    505
#define GP_BPM_GPU_BUFFER_CANDIDATES_LENGTH                 506
#define GP_BPM_GPU_BUFFER_USAGE_CANDIDATES                  507
#define GP_BPM_GPU_BUFFER_USAGE_QUERIES                     508
#define GP_BPM_GPU_BUFFER_USAGE_PEQ_ENTRIES                 509

/*
 * FM-Index
 */
#define GP_FMIDX_LOOKUP_DIST                                550

/*
 * Matches Align
 */
#define GP_MATCHES_ALIGN_EXACT                              580
#define GP_MATCHES_ALIGN_HAMMING                            581
#define GP_MATCHES_ALIGN_LEVENSHTEIN                        582
#define GP_MATCHES_ALIGN_SWG                                583
#define GP_MATCHES_ALIGN_LOCAL_SWG                          584

/*
 * SWG
 */
#define GP_SWG_ALIGN_FULL                                   600
#define GP_SWG_ALIGN_FULL_LENGTH                            601
#define GP_SWG_ALIGN_BANDED                                 602
#define GP_SWG_ALIGN_BANDED_LENGTH                          603

/*
 * Neighborhood Search
 */
#define GP_NS_BEST_MATCH                                    620

#define GP_NS_NODES_EXPLORED                                621
#define GP_NS_NODES_EXPLORED_MTABLE                         622
#define GP_NS_NODE_SUCCESS                                  623
#define GP_NS_FAILED_OPT                                    624
#define GP_NS_NODE_CLOSED                                   625
#define GP_NS_NODE_CLOSED_DEPTH                             626

/*
 * Checks
 */
#define GP_CHECK_NUM_READS                                  650
#define GP_CHECK_NUM_MAPS                                   651
#define GP_CHECK_INCORRECT                                  652
#define GP_CHECK_SUBOPTIMAL                                 653
#define GP_CHECK_SUBOPTIMAL_SUBDOMINANT                     654
#define GP_CHECK_SUBOPTIMAL_DIFF                            655
#define GP_CHECK_SUBOPTIMAL_SCORE                           656
#define GP_CHECK_SUBOPTIMAL_DISTANCE                        657

/*
 * Dummy
 */
#define GP_DUMMY1                                           700
#define GP_DUMMY2                                           701

/*
 * Mapper SE
 */
void mapper_profile_print_mapper_single_end(
    FILE* const stream,const bool map_output,const uint64_t num_threads);
void mapper_profile_print_mapper_single_end_cuda(
    FILE* const stream,const bool map_output,const uint64_t num_threads);

/*
 * Mapper PE
 */
void mapper_profile_print_mapper_paired_end(
    FILE* const stream,const bool map_output,const uint64_t num_threads);
void mapper_profile_print_mapper_paired_end_cuda(
    FILE* const stream,const bool map_output,const uint64_t num_threads);

/*
 * Error Messages
 */
//#define GEM_ERROR_MAPPER_PROFILE ""

#endif /* MAPPER_PROFILE_H_ */
