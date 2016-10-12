/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Counters Codes. Basically, these counters/timers focus on the measurement of:
 *     - Time spent in certain functions/blocks-code (and number of calls/executions)
 *     - Functional profile of execution flow branches (counting)
 */

#ifndef MAPPER_PROFILE_COUNTERS_H_
#define MAPPER_PROFILE_COUNTERS_H_

/*
 * Mapper
 */
#define GP_MAPPER_ALL                                         0
#define GP_MAPPER_LOAD_INDEX                                  1
#define GP_MAPPER_NUM_READS                                   2

#define GP_OUTPUT_MAP_SE                                      3
#define GP_OUTPUT_MAP_PE                                      4
#define GP_OUTPUT_SAM_SE                                      5
#define GP_OUTPUT_SAM_PE                                      6

/*
 * Archive Search SE
 */
#define GP_ARCHIVE_SEARCH_SE                                 50
#define GP_ARCHIVE_SEARCH_SE_INIT                            51
#define GP_ARCHIVE_SEARCH_PE_INIT                            52
#define GP_ARCHIVE_SEARCH_SE_PREPARE_SEQUENCE                53
#define GP_ARCHIVE_SEARCH_SE_GENERATE_CANDIDATES             54
#define GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES               55
#define GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH                   56

// StepWise
#define GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_GENERATE         60
#define GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_COPY             61
#define GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_RETRIEVE         62
#define GP_ARCHIVE_SEARCH_SE_DECODE_CANDIDATES_GENERATE      63
#define GP_ARCHIVE_SEARCH_SE_DECODE_CANDIDATES_COPY          64
#define GP_ARCHIVE_SEARCH_SE_DECODE_CANDIDATES_RETRIEVE      65
#define GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES_GENERATE      66
#define GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES_COPY          67
#define GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES_RETRIEVE      68

#define GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_GENERATE         70
#define GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_COPY             71
#define GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_RETRIEVE         72
#define GP_ARCHIVE_SEARCH_PE_DECODE_CANDIDATES_GENERATE      73
#define GP_ARCHIVE_SEARCH_PE_DECODE_CANDIDATES_COPY          74
#define GP_ARCHIVE_SEARCH_PE_DECODE_CANDIDATES_RETRIEVE      75
#define GP_ARCHIVE_SEARCH_PE_VERIFY_CANDIDATES_GENERATE      76
#define GP_ARCHIVE_SEARCH_PE_VERIFY_CANDIDATES_COPY          77
#define GP_ARCHIVE_SEARCH_PE_VERIFY_CANDIDATES_RETRIEVE      78

// Post-processing
#define GP_ARCHIVE_SELECT_SE_MATCHES                        116
#define GP_ARCHIVE_SCORE_SE_MATCHES                         117

/*
 * Archive Search PE
 */
#define GP_ARCHIVE_SEARCH_PE                                122
#define GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT             123
#define GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT_TOTAL       124
#define GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT_SUCCESS     125
#define GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY             126
#define GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY_TOTAL       127
#define GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY_SUCCESS     128
#define GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES              129
#define GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_TOTAL        130
#define GP_ARCHIVE_SEARCH_PE_EXTEND_NUM_MATCHES             131
#define GP_ARCHIVE_SEARCH_PE_GENERATE_CANDIDATES            132
#define GP_ARCHIVE_SEARCH_PE_FIND_PAIRS                     133
#define GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH                  134

// Post-processing
#define GP_ARCHIVE_SELECT_PE_MATCHES                        142
#define GP_ARCHIVE_SCORE_PE_MATCHES                         143

/*
 * I/O
 */
#define GP_INPUT_FASTA_PARSE_SEQUENCE                       145
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
#define GP_AS_MCS                                           161

#define GP_AS_FILTERING_ADATIVE                             165
#define GP_AS_FILTERING_ADATIVE_CALL                        166
#define GP_AS_FILTERING_ADATIVE_MCS                         167
#define GP_AS_NEIGHBORHOOD_SEARCH                           168
#define GP_AS_NEIGHBORHOOD_SEARCH_MCS                       169
#define GP_AS_NEIGHBORHOOD_SEARCH_CALL                      170
#define GP_AS_LOCAL_ALIGN                                   171
#define GP_AS_LOCAL_ALIGN_CALL                              172

#define GP_AS_GENERATE_CANDIDATES                           205
#define GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS      206
#define GP_AS_GENERATE_CANDIDATES_PROCESSED                 205
#define GP_AS_GENERATE_CANDIDATES_SKIPPED                   206
#define GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING         207

#define GP_ASSW_REGION_PROFILE_UNSUCCESSFUL                 220
#define GP_ASSW_DECODE_CANDIDATES_COPIED                    221
#define GP_ASSW_DECODE_CANDIDATES_RETRIVED                  222
#define GP_ASSW_VERIFY_CANDIDATES_TILES_COPIED              223
#define GP_ASSW_VERIFY_CANDIDATES_TILES_RETRIVED            224

/*
 * Region Profile
 */
#define GP_REGION_PROFILE_ADAPTIVE                          250
#define GP_REGION_PROFILE_QUIT_PROFILE                      251

#define GP_REGION_PROFILE                                   259
#define GP_REGION_PROFILE_NUM_REGIONS                       260
#define GP_REGION_PROFILE_NUM_REGIONS_UNIQUE                261
#define GP_REGION_PROFILE_NUM_REGIONS_STANDARD              262
#define GP_REGION_PROFILE_REGION_LENGTH                     263
#define GP_REGION_PROFILE_REGION_CANDIDATES                 264
#define GP_REGION_PROFILE_TOTAL_CANDIDATES                  265

/*
 * Filtering Candidates (Verifying)
 */
#define GP_FC_VERIFICATION                                  330
#define GP_FC_PROCESS_CANDIDATES                            331
#define GP_FC_COMPOSE_REGIONS                               332
#define GP_FC_DECODE_POSITIONS                              333
#define GP_FC_DECODE_CANDIDATES_BUFFERED                    334
#define GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL       335
#define GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL_TOTAL 336
#define GP_FC_VERIFY_CANDIDATES                             337
#define GP_FC_VERIFY_CANDIDATES_REGION                      338
#define GP_FC_VERIFY_CANDIDATES_BUFFERED                    339
#define GP_FC_VERIFY_CANDIDATES_BUFFERED_DDIFF              340
#define GP_FC_REALIGN_CANDIDATE_REGIONS                     341
#define GP_FC_REALIGN_LOCAL_CANDIDATE_REGIONS               342

#define GP_FC_KMER_COUNTER_FILTER                           350
#define GP_FC_KMER_COUNTER_FILTER_NA                        351
#define GP_FC_KMER_COUNTER_FILTER_DISCARDED                 352
#define GP_FC_KMER_COUNTER_FILTER_ACCEPTED                  353

#define GP_FC_EXTEND_MATCH                                  360
#define GP_FC_EXTEND_RETRIEVE_CANDIDATE_REGIONS             361
#define GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS               362
#define GP_FC_EXTEND_VERIFY_CANDIDATES_LENGTH               363
#define GP_FC_EXTEND_VERIFY_CANDIDATES_FOUND                364
#define GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS              365

#define GP_FC_CACHE_SEARCH_HIT                              372
#define GP_FC_SELECT_PRUNE_HIT                              373

#define GP_CANDIDATE_POSITIONS                              400
#define GP_CANDIDATE_REGIONS                                402
#define GP_CANDIDATE_REGIONS_DUPLICATED                     403
#define GP_CANDIDATE_REGION_LENGTH                          404
#define GP_CANDIDATE_REGION_ALIGNMENT_REGIONS_TOTAL         405
#define GP_CANDIDATE_REGION_ALIGNMENT_COVERAGE              406
#define GP_CANDIDATE_REGION_LOCAL                           407
#define GP_CANDIDATE_REGION_LOCAL_ALIGNED                   408

#define GP_ACCEPTED_REGIONS                                 410
#define GP_DISCARDED_REGIONS                                411

#define GP_ALIGNED_REGIONS                                  415
#define GP_ALIGNED_EXACT                                    416
#define GP_ALIGNED_INEXACT                                  417
#define GP_ALIGNED_REGIONS_LENGTH                           418
#define GP_ALIGNED_ACCEPTED                                 419
#define GP_ALIGNED_DISCARDED                                420
#define GP_ALIGNED_DISCARDED_SWG                            421
#define GP_ALIGNED_DISCARDED_SWG_THRESHOLD                  422
#define GP_ALIGNED_DISCARDED_MATCHING_BASES                 423
#define GP_ALIGNED_DISCARDED_SCAFFOLD                       424

#define GP_MATCH_SCAFFOLD_ALIGNMENT                         450
#define GP_MATCH_SCAFFOLD_ALIGNMENT_REGIONS                 452
#define GP_MATCH_SCAFFOLD_ALIGNMENT_COVERAGE                453
#define GP_MATCH_SCAFFOLD_ALIGNMENT_ADAPTIVE                451
#define GP_MATCH_SCAFFOLD_CHAIN_REGIONS                     454
#define GP_MATCH_SCAFFOLD_CHAIN_REGIONS_SCAFFOLDS           455
#define GP_MATCH_SCAFFOLD_CHAIN_REGIONS_COVERAGE            456
#define GP_MATCH_SCAFFOLD_EDIT                              457
#define GP_MATCH_SCAFFOLD_EDIT_SCAFFOLDS                    458
#define GP_MATCH_SCAFFOLD_EDIT_COVERAGE                     459
#define GP_MATCH_SCAFFOLD_EDIT_TILES_TOTAL                  460
#define GP_MATCH_SCAFFOLD_EDIT_TILES_ALIGN                  461
#define GP_MATCH_SCAFFOLD_EDIT_TILES_SKIPPED                462
#define GP_MATCH_SCAFFOLD_EDIT_CELLS                        463

#define GP_BPM_DISTANCE                                     470
#define GP_BMP_DISTANCE_NUM_TILES                           471
#define GP_BMP_DISTANCE_NUM_TILES_VERIFIED                  472
#define GP_BPM_DISTANCE_QUICK_ABANDON                       473
#define GP_BPM_DISTANCE_TEXT_LENGTH                         474
#define GP_BPM_DISTANCE_KEY_LENGTH                          475
#define GP_BPM_DISTANCE_CELLS                               476

#define GP_BPM_ALL                                          480
#define GP_BPM_ALL_QUICK_ABANDON                            481
#define GP_BPM_ALL_MATCHES_FOUND                            482


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
 * Matches
 */
#define GP_MATCHES_MAPS_ADDED                               590
#define GP_MATCHES_MAPS_DUP                                 591

#define GP_MATCHES_ACCURACY_CASE_CALLS                      600
#define GP_MATCHES_ACCURACY_CASE_HIT                        601
#define GP_MATCHES_ACCURACY_CASE_MISS                       602

#define GP_MATCHES_ACCURACY_CASE_MAX_MATCHES                603
#define GP_MATCHES_ACCURACY_CASE_UNMAPPED                   604
#define GP_MATCHES_ACCURACY_CASE_TIE                        605
#define GP_MATCHES_ACCURACY_CASE_MMAP_D1                    606
#define GP_MATCHES_ACCURACY_CASE_MMAP                       607
#define GP_MATCHES_ACCURACY_CASE_UNIQUE                     608
#define GP_MATCHES_ACCURACY_CASE_MAP_INCOMPLETE             609
#define GP_MATCHES_ACCURACY_CASE_MAP_FRONTIER               610

/*
 * SWG
 */
#define GP_SWG_ALIGN_FULL                                   620
#define GP_SWG_ALIGN_FULL_LENGTH                            621
#define GP_SWG_ALIGN_BANDED                                 622
#define GP_SWG_ALIGN_BANDED_LENGTH                          623
#define GP_SWG_ALIGN_BANDED_CELLS                           624

/*
 * Neighborhood Search
 */
#define GP_NS_GENERATION                                    630
#define GP_NS_VERIFICATION                                  631
#define GP_NS_ALIGN                                         632

#define GP_NS_NODES                                         635
#define GP_NS_NODES_SUCCESS                                 636
#define GP_NS_NODES_FAIL                                    637
#define GP_NS_NODES_DP_PREVENT                              638
#define GP_NS_DP_CELLS_COMPUTED                             640

#define GP_NS_SEARCH_DEPTH                                  645
#define GP_NS_SEARCH_CANDIDATES_GENERATED                   646
#define GP_NS_BRANCH_CANDIDATES_GENERATED                   647

/*
 * Checks
 */
#define GP_CHECK_NUM_READS                                  650
#define GP_CHECK_NUM_MAPS                                   651
#define GP_CHECK_INCORRECT                                  652
#define GP_CHECK_PRIMARY_SUBOPTIMAL                         653
#define GP_CHECK_PRIMARY_SUBOPTIMAL_FAIL                    654
#define GP_CHECK_PRIMARY_SUBOPTIMAL_SCORE                   656
#define GP_CHECK_PRIMARY_SUBOPTIMAL_DIFF                    655
#define GP_CHECK_PRIMARY_SUBOPTIMAL_DISTANCE                657
#define GP_CHECK_SUBDOMINANT_SUBOPTIMAL                     658
#define GP_CHECK_SUBDOMINANT_SUBOPTIMAL_FAIL                659
#define GP_CHECK_SUBDOMINANT_SUBOPTIMAL_SCORE               660
#define GP_CHECK_SUBDOMINANT_SUBOPTIMAL_DIFF                661
#define GP_CHECK_SUBDOMINANT_SUBOPTIMAL_DISTANCE            662

/*
 * Dummy
 */
#define GP_DUMMY1                                           700
#define GP_DUMMY2                                           701

/*
 * CUDA SE Mapper
 */
#define GP_MAPPER_CUDA_SE                                    710
#define GP_MAPPER_CUDA_SE_REGION_PROFILE                     711
#define GP_MAPPER_CUDA_SE_DECODE_CANDIDATES                  712
#define GP_MAPPER_CUDA_SE_VERIFY_CANDIDATES                  713
#define GP_MAPPER_CUDA_SE_FINISH_SEARCH                      714

#define GP_SEARCH_STAGE_REGION_PROFILE_BUFFERS_USED          720
#define GP_SEARCH_STAGE_REGION_PROFILE_SEARCHES_IN_BUFFER    721
#define GP_SEARCH_STAGE_DECODE_CANDIDATES_BUFFERS_USED       722
#define GP_SEARCH_STAGE_DECODE_CANDIDATES_SEARCHES_IN_BUFFER 723
#define GP_SEARCH_STAGE_VERIFY_CANDIDATES_BUFFERS_USED       724
#define GP_SEARCH_STAGE_VERIFY_CANDIDATES_SEARCHES_IN_BUFFER 725

/*
 * CUDA PE Mapper
 */
#define GP_MAPPER_CUDA_PE                                    750
#define GP_MAPPER_CUDA_PE_REGION_PROFILE                     751
#define GP_MAPPER_CUDA_PE_DECODE_CANDIDATES                  752
#define GP_MAPPER_CUDA_PE_VERIFY_CANDIDATES                  753
#define GP_MAPPER_CUDA_PE_FINISH_SEARCH                      754

/*
 * GPU BUFFERS
 */
#define GP_GPU_BUFFER_COLLECTION_INIT                        790

#define GP_GPU_BUFFER_FMI_SEARCH_ALLOC                       800
#define GP_GPU_BUFFER_FMI_SEARCH_NUM_QUERIES                 801
#define GP_GPU_BUFFER_FMI_SEARCH_SEND                        802
#define GP_GPU_BUFFER_FMI_SEARCH_RECEIVE                     803
#define GP_GPU_BUFFER_FMI_SEARCH_USAGE_QUERIES               804
#define GP_GPU_BUFFER_FMI_SEARCH_DUTY_CYCLE                  805

#define GP_GPU_BUFFER_FMI_DECODE_ALLOC                       810
#define GP_GPU_BUFFER_FMI_DECODE_NUM_QUERIES                 811
#define GP_GPU_BUFFER_FMI_DECODE_SEND                        812
#define GP_GPU_BUFFER_FMI_DECODE_RECEIVE                     813
#define GP_GPU_BUFFER_FMI_DECODE_USAGE_CANDIDATES            814
#define GP_GPU_BUFFER_FMI_DECODE_DUTY_CYCLE                  815

#define GP_GPU_BUFFER_ALIGN_BPM_ALLOC                        820
#define GP_GPU_BUFFER_ALIGN_BPM_NUM_QUERIES                  821
#define GP_GPU_BUFFER_ALIGN_BPM_CANDIDATE_LENGTH             822
#define GP_GPU_BUFFER_ALIGN_BPM_CANDIDATE_PER_TILE           823
#define GP_GPU_BUFFER_ALIGN_BPM_CELLS                        824
#define GP_GPU_BUFFER_ALIGN_BPM_SEND                         825
#define GP_GPU_BUFFER_ALIGN_BPM_RECEIVE                      826
#define GP_GPU_BUFFER_ALIGN_BPM_DUTY_CYCLE                   827
#define GP_GPU_BUFFER_ALIGN_BPM_USAGE_CANDIDATES             828
#define GP_GPU_BUFFER_ALIGN_BPM_USAGE_QUERIES                829
#define GP_GPU_BUFFER_ALIGN_BPM_USAGE_PEQ_ENTRIES            830

#endif /* MAPPER_PROFILE_COUNTERS_H_ */
