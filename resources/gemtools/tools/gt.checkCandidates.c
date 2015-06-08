/*
 * PROJECT: GEM-Tools library
 * FILE: gt.checkCandidates.c
 * DATE: 08/10/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include <gem_tools.h>
#include <getopt.h>
#include <gt_buffered_input_file.h>
#include <gt_buffered_output_file.h>
#include <gt_commons.h>
#include <gt_dna_string.h>
#include <gt_error.h>
#include <gt_input_file.h>
#include <gt_input_parser.h>
#include <gt_map.h>
#include <gt_map_align.h>
#include <gt_map_align_bpm.h>
#include <gt_map_align_bpm_simd.h>
#include <gt_mm.h>
#include <gt_options_menu.h>
#include <gt_output_file.h>
#include <gt_profiler.h>
#include <gt_sequence_archive.h>
#include <gt_string.h>
#include <gt_vector.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

/*
 * Constants
 */
#define GT_CC_EXPANSION_FACTOR 0.20

#define GT_GPROF_THREAD_ALL 0
#define GT_GPROF_NUM_CANDIDATES 1
#define GT_GPROF_RETRIEVE_TEXT 2
#define GT_GPROF_DP_DISTANCE 3

#define GT_GPROF_BASE_BPM64_DISTANCE 4
#define GT_GPROF_BASE_BPM64_CUTOFF_DISTANCE 5

#define GT_GPROF_BPM8_DISTANCE 6
#define GT_GPROF_BPM8_CUTOFF_DISTANCE 7
#define GT_GPROF_BPM16_DISTANCE 8
#define GT_GPROF_BPM16_CUTOFF_DISTANCE 9
#define GT_GPROF_BPM32_DISTANCE 10
#define GT_GPROF_BPM32_CUTOFF_DISTANCE 11
#define GT_GPROF_BPM64_DISTANCE 12
#define GT_GPROF_BPM64_CUTOFF_DISTANCE 13
#define GT_GPROF_BPM128_DISTANCE 14
#define GT_GPROF_BPM128_CUTOFF_DISTANCE 15

#define GT_GPROF_GLOBAL 16

#define GT_GPROF_BPM_SIMD128_DISTANCE 17
#define GT_GPROF_BPM_SIMD128_UNFOLDED_DISTANCE 18

/*
 * Profile
 */
gt_profile** tprof;

/*
 * Parameters
 */
typedef struct {
  /* I/O */
  bool paired_end;
  char* name_input_file;
  char* name_output_file;
  char* name_reference_file;
  char* name_gem_index_file;
  float error;
  /* Misc */
  uint64_t num_threads;
  bool verbose;
} gt_check_candidates_args;
gt_check_candidates_args parameters = {
  /* I/O */
  .name_input_file=NULL,
  .paired_end=false,
  .name_output_file=NULL,
  .name_reference_file=NULL,
  .name_gem_index_file=NULL,
  .error = 0.20,
  /* Misc */
  .num_threads=1,
  .verbose=false,
};
/*
 * Candidates Job
 */
typedef struct {
  gt_string* reference_read;
  gt_vector* candidates; // (gt_map*)
  uint64_t job_id;
} gt_candidates_job;

GT_INLINE gt_candidates_job* gt_candidates_job_new() {
  gt_candidates_job* const candidates_job = gt_alloc(gt_candidates_job);
  candidates_job->reference_read = gt_string_new(1000);
  candidates_job->candidates = gt_vector_new(200,sizeof(gt_map*));
  candidates_job->job_id = 0;
  return candidates_job;
}
GT_INLINE void gt_candidates_job_clean(gt_candidates_job* const candidates_job) {
  gt_string_clear(candidates_job->reference_read);
  GT_VECTOR_ITERATE(candidates_job->candidates,candidate,position,gt_map*) {
    gt_map_delete(*candidate);
  }
  gt_vector_clear(candidates_job->candidates);
  candidates_job->job_id++;
}
GT_INLINE void gt_candidates_job_delete(gt_candidates_job* const candidates_job) {
  gt_candidates_job_clean(candidates_job);
  gt_string_delete(candidates_job->reference_read);
  gt_vector_delete(candidates_job->candidates);
}
GT_INLINE void gt_candidates_job_add_candidate(gt_candidates_job* const candidates_job,gt_map* const map) {
  gt_vector_insert(candidates_job->candidates,map,gt_map*);
}
/*
 * Cheching candidates
 */
#define gt_map_block_bpm_get_distance__cutoff_uint128_t gt_map_block_bpm_get_distance__cutoff___int128
#define gt_map_block_bpm_get_distance_uint128_t gt_map_block_bpm_get_distance___int128

#define GT_CHECK_CANDIDATES_CALL_BMP(VL,true_position,true_distance) { \
  /* Resulting Position & Distance */ \
  uint64_t __bpm_position, __bpm_distance; \
  /* Calculate distance using BPM (Myers) */ \
  GPROF_INC_COUNTER(tprof[thread_id],GT_GPROF_BPM##VL##_DISTANCE); \
  GPROF_START_TIMER(tprof[thread_id],GT_GPROF_BPM##VL##_DISTANCE); \
  gt_map_block_bpm_get_distance_uint##VL##_t(bpm##VL##_pattern, \
      gt_string_get_string(candidates_sequence),gt_string_get_length(candidates_sequence), \
      &__bpm_position,&__bpm_distance); \
  GPROF_STOP_TIMER(tprof[thread_id],GT_GPROF_BPM##VL##_DISTANCE); \
  /* Check results */ \
  gt_cond_error_msg(true_position!=__bpm_position || true_distance!=__bpm_distance, \
      "BitParalell algorithms error at BMP-%d implementation",VL); \
  /* Calculate distance using BPM-cutoff (Myers) */  \
  GPROF_INC_COUNTER(tprof[thread_id],GT_GPROF_BPM##VL##_CUTOFF_DISTANCE); \
  GPROF_START_TIMER(tprof[thread_id],GT_GPROF_BPM##VL##_CUTOFF_DISTANCE); \
  gt_map_block_bpm_get_distance__cutoff_uint##VL##_t(bpm##VL##_pattern, \
      gt_string_get_string(candidates_sequence),gt_string_get_length(candidates_sequence), \
      &__bpm_position,&__bpm_distance,max_distance); \
  GPROF_STOP_TIMER(tprof[thread_id],GT_GPROF_BPM##VL##_CUTOFF_DISTANCE); \
  /* Check results */ \
  gt_cond_error_msg((true_distance<=max_distance) && \
      (true_position!=__bpm_position || true_distance!=__bpm_distance), \
      "BitParalell algorithms error at BMP-%d-CUTOFF implementation",VL); \
}
GT_INLINE void gt_check_candidates_align(
    const uint64_t thread_id,gt_sequence_archive* const sequence_archive,
    gt_candidates_job* const candidates_job,gt_vector* const dp_buffer,
    gt_bpm_pattern* const bpm8_pattern,gt_bpm_pattern* const bpm16_pattern,
    gt_bpm_pattern* const bpm32_pattern,gt_bpm_pattern* const bpm64_pattern,
    gt_bpm_pattern* const bpm128_pattern) {
  // Prepare placeholder for reference sequences
  const uint64_t reference_length = gt_string_get_length(candidates_job->reference_read);
  const uint64_t extra_length =
      (uint64_t) ceil(GT_CC_EXPANSION_FACTOR*(float)reference_length);
  gt_string* const candidates_sequence = gt_string_new(reference_length+2*extra_length+1);
  // Compile pattern for Myers
  const uint64_t max_distance = (uint64_t)((float)gt_string_get_length(candidates_job->reference_read)*parameters.error);
  gt_map_bpm_pattern_compile(bpm8_pattern,1,
      gt_string_get_string(candidates_job->reference_read),
      gt_string_get_length(candidates_job->reference_read));
  gt_map_bpm_pattern_compile(bpm16_pattern,2,
      gt_string_get_string(candidates_job->reference_read),
      gt_string_get_length(candidates_job->reference_read));
  gt_map_bpm_pattern_compile(bpm32_pattern,4,
      gt_string_get_string(candidates_job->reference_read),
      gt_string_get_length(candidates_job->reference_read));
  gt_map_bpm_pattern_compile(bpm64_pattern,8,
      gt_string_get_string(candidates_job->reference_read),
      gt_string_get_length(candidates_job->reference_read));
  gt_map_bpm_pattern_compile(bpm128_pattern,16,
      gt_string_get_string(candidates_job->reference_read),
      gt_string_get_length(candidates_job->reference_read));
  // Align all candidates
  GT_VECTOR_ITERATE(candidates_job->candidates,candidate_ptr,pos,gt_map*) {
    GPROF_INC_COUNTER(tprof[thread_id],GT_GPROF_NUM_CANDIDATES);
    gt_map* const candidate = *candidate_ptr;

    // Retrieve reference sequence (candidate text)
    GPROF_INC_COUNTER(tprof[thread_id],GT_GPROF_RETRIEVE_TEXT);
    GPROF_START_TIMER(tprof[thread_id],GT_GPROF_RETRIEVE_TEXT);
    if (gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
        gt_map_get_seq_name(candidate),FORWARD,gt_map_get_position(candidate),
        reference_length+extra_length,extra_length,candidates_sequence)) {
      continue; // Ignore negative positions due to gem-mapper alignment offset
    }
    GPROF_STOP_TIMER(tprof[thread_id],GT_GPROF_RETRIEVE_TEXT);

     /* Calculate distance using DP-levenshtein (reference test) */
    GPROF_INC_COUNTER(tprof[thread_id],GT_GPROF_DP_DISTANCE);
    GPROF_START_TIMER(tprof[thread_id],GT_GPROF_DP_DISTANCE);
    uint64_t dp_position, dp_distance;
    gt_map_realign_levenshtein_get_distance(
        gt_string_get_string(candidates_job->reference_read),gt_string_get_length(candidates_job->reference_read),
        gt_string_get_string(candidates_sequence),gt_string_get_length(candidates_sequence),
        true,&dp_position,&dp_distance,dp_buffer);
    GPROF_STOP_TIMER(tprof[thread_id],GT_GPROF_DP_DISTANCE);

    /* Calculate distance using base implementation of BPM (Myers) */
    GPROF_INC_COUNTER(tprof[thread_id],GT_GPROF_BASE_BPM64_DISTANCE);
    GPROF_START_TIMER(tprof[thread_id],GT_GPROF_BASE_BPM64_DISTANCE);
    uint64_t reference_position, reference_distance;
    gt_map_block_bpm_get_distance(bpm64_pattern,
        gt_string_get_string(candidates_sequence),gt_string_get_length(candidates_sequence),
        &reference_position,&reference_distance,max_distance);
    GPROF_STOP_TIMER(tprof[thread_id],GT_GPROF_BASE_BPM64_DISTANCE);
    /* Calculate distance using base implementation of BPM (Myers) + CUTTOFF */
    GPROF_INC_COUNTER(tprof[thread_id],GT_GPROF_BASE_BPM64_CUTOFF_DISTANCE);
    GPROF_START_TIMER(tprof[thread_id],GT_GPROF_BASE_BPM64_CUTOFF_DISTANCE);
    uint64_t reference_position_cutoff, reference_distance_cutoff;
    gt_map_block_bpm_get_distance__cutoff(bpm64_pattern,
        gt_string_get_string(candidates_sequence),gt_string_get_length(candidates_sequence),
        &reference_position_cutoff,&reference_distance_cutoff,max_distance);
    GPROF_STOP_TIMER(tprof[thread_id],GT_GPROF_BASE_BPM64_CUTOFF_DISTANCE);
    /* Check consistency */
    if (reference_position!=dp_position || reference_distance!=dp_distance ||
        ((dp_distance<=max_distance) &&
         (reference_position_cutoff!=dp_position || reference_distance_cutoff!=dp_distance))) {
      gt_error_msg("DynamicPrograming approach returns different results");
    }

    // Calculate distance using generic optimized vector BPM-8 (Myers)
    GT_CHECK_CANDIDATES_CALL_BMP(8,reference_position,reference_distance);
    // Calculate distance using generic optimized vector BPM-16 (Myers)
    GT_CHECK_CANDIDATES_CALL_BMP(16,reference_position,reference_distance);
    // Calculate distance using generic optimized vector BPM-32 (Myers)
    GT_CHECK_CANDIDATES_CALL_BMP(32,reference_position,reference_distance);
    // Calculate distance using generic optimized vector BPM-64 (Myers)
    GT_CHECK_CANDIDATES_CALL_BMP(64,reference_position,reference_distance);
    // Calculate distance using generic optimized vector BPM-128 (Myers)
    GT_CHECK_CANDIDATES_CALL_BMP(128,reference_position,reference_distance);

    // SIMD Implementations
    uint64_t bpm_simd_position, bpm_simd_distance;
    /* Calculate distance using BPM-SIMD unfolded-64 */
    GPROF_INC_COUNTER(tprof[thread_id],GT_GPROF_BPM_SIMD128_UNFOLDED_DISTANCE);
    GPROF_START_TIMER(tprof[thread_id],GT_GPROF_BPM_SIMD128_UNFOLDED_DISTANCE);
    gt_map_block_bpm_get_distance_simd128_unfold64(bpm128_pattern,
        gt_string_get_string(candidates_sequence),gt_string_get_length(candidates_sequence),
        &bpm_simd_position,&bpm_simd_distance);
    GPROF_STOP_TIMER(tprof[thread_id],GT_GPROF_BPM_SIMD128_UNFOLDED_DISTANCE);
    gt_cond_error_msg(reference_position!=bpm_simd_position || reference_distance!=bpm_simd_distance,
        "BitParalell algorithms error at SIMD-128 (unfolded) implementation");
    /* Calculate distance using BPM-SIMD SSE4.2 */
    GPROF_INC_COUNTER(tprof[thread_id],GT_GPROF_BPM_SIMD128_DISTANCE);
    GPROF_START_TIMER(tprof[thread_id],GT_GPROF_BPM_SIMD128_DISTANCE);
    gt_map_block_bpm_get_distance_simd128(bpm128_pattern,
        gt_string_get_string(candidates_sequence),gt_string_get_length(candidates_sequence),
        &bpm_simd_position,&bpm_simd_distance);
    GPROF_STOP_TIMER(tprof[thread_id],GT_GPROF_BPM_SIMD128_DISTANCE);
    gt_cond_error_msg(reference_position!=bpm_simd_position || reference_distance!=bpm_simd_distance,
        "BitParalell algorithms error at SIMD-128 implementation");

  }
  // Free
  gt_string_delete(candidates_sequence);
}
/*
 * Display Profile
 */
#define GT_CHECK_CANDIDATES_SHOW_PROFILE_BMP(VL) \
  fprintf(stderr,"  --> BPM-" #VL ".distance  %2.3f (%2.3f%%)\t(%lu calls)\t(%2.3f ms/call)\n", \
      GPROF_GET_TIMER(tprof[thread_id],GT_GPROF_BPM##VL##_DISTANCE), \
      GPROF_TIME_PERCENTAGE(tprof[thread_id],GT_GPROF_BPM##VL##_DISTANCE,GT_GPROF_THREAD_ALL), \
      GPROF_GET_COUNTER(tprof[thread_id],GT_GPROF_BPM##VL##_DISTANCE), \
      GPROF_TIME_PER_CALL(tprof[thread_id],GT_GPROF_BPM##VL##_DISTANCE,GT_GPROF_BPM##VL##_DISTANCE)); \
  fprintf(stderr,"  --> BPM-" #VL ".CUTOFF.distance  %2.3f (%2.3f%%)\t(%lu calls)\t(%2.3f ms/call)\n", \
      GPROF_GET_TIMER(tprof[thread_id],GT_GPROF_BPM##VL##_CUTOFF_DISTANCE), \
      GPROF_TIME_PERCENTAGE(tprof[thread_id],GT_GPROF_BPM##VL##_CUTOFF_DISTANCE,GT_GPROF_THREAD_ALL), \
      GPROF_GET_COUNTER(tprof[thread_id],GT_GPROF_BPM##VL##_CUTOFF_DISTANCE), \
      GPROF_TIME_PER_CALL(tprof[thread_id],GT_GPROF_BPM##VL##_CUTOFF_DISTANCE,GT_GPROF_BPM##VL##_CUTOFF_DISTANCE))
GT_INLINE void gt_check_candidates_show_profile(const uint64_t thread_id) {
  fprintf(stderr,"[Thread-%lu]Total.Time %2.3f\n",thread_id,GPROF_GET_TIMER(tprof[thread_id],GT_GPROF_THREAD_ALL));
  fprintf(stderr,"  --> Retrieve.Text %2.3f (%2.3f%%)\t(%lu calls)\t(%2.3f ms/call)\n",
      GPROF_GET_TIMER(tprof[thread_id],GT_GPROF_RETRIEVE_TEXT),
      GPROF_TIME_PERCENTAGE(tprof[thread_id],GT_GPROF_RETRIEVE_TEXT,GT_GPROF_THREAD_ALL),
      GPROF_GET_COUNTER(tprof[thread_id],GT_GPROF_RETRIEVE_TEXT),
      GPROF_TIME_PER_CALL(tprof[thread_id],GT_GPROF_RETRIEVE_TEXT,GT_GPROF_RETRIEVE_TEXT));
  /* STD Dynamic Programming */
  fprintf(stderr,"  --> DP.distance   %2.3f (%2.3f%%)\t(%lu calls)\t(%2.3f ms/call)\n",
      GPROF_GET_TIMER(tprof[thread_id],GT_GPROF_DP_DISTANCE),
      GPROF_TIME_PERCENTAGE(tprof[thread_id],GT_GPROF_DP_DISTANCE,GT_GPROF_THREAD_ALL),
      GPROF_GET_COUNTER(tprof[thread_id],GT_GPROF_DP_DISTANCE),
      GPROF_TIME_PER_CALL(tprof[thread_id],GT_GPROF_DP_DISTANCE,GT_GPROF_DP_DISTANCE));
  /* STD Myers Base Implementation */
  fprintf(stderr,"  --> BPM.Base.distance   %2.3f (%2.3f%%)\t(%lu calls)\t(%2.3f ms/call)\n",
      GPROF_GET_TIMER(tprof[thread_id],GT_GPROF_BASE_BPM64_DISTANCE),
      GPROF_TIME_PERCENTAGE(tprof[thread_id],GT_GPROF_BASE_BPM64_DISTANCE,GT_GPROF_THREAD_ALL),
      GPROF_GET_COUNTER(tprof[thread_id],GT_GPROF_BASE_BPM64_DISTANCE),
      GPROF_TIME_PER_CALL(tprof[thread_id],GT_GPROF_BASE_BPM64_DISTANCE,GT_GPROF_BASE_BPM64_DISTANCE));
  fprintf(stderr,"  --> BPM.Base.CUTOFF.distance   %2.3f (%2.3f%%)\t(%lu calls)\t(%2.3f ms/call)\n",
      GPROF_GET_TIMER(tprof[thread_id],GT_GPROF_BASE_BPM64_CUTOFF_DISTANCE),
      GPROF_TIME_PERCENTAGE(tprof[thread_id],GT_GPROF_BASE_BPM64_CUTOFF_DISTANCE,GT_GPROF_THREAD_ALL),
      GPROF_GET_COUNTER(tprof[thread_id],GT_GPROF_BASE_BPM64_CUTOFF_DISTANCE),
      GPROF_TIME_PER_CALL(tprof[thread_id],GT_GPROF_BASE_BPM64_CUTOFF_DISTANCE,GT_GPROF_BASE_BPM64_CUTOFF_DISTANCE));
  /* Generic opt. vector impl.*/
  GT_CHECK_CANDIDATES_SHOW_PROFILE_BMP(8);
  GT_CHECK_CANDIDATES_SHOW_PROFILE_BMP(16);
  GT_CHECK_CANDIDATES_SHOW_PROFILE_BMP(32);
  GT_CHECK_CANDIDATES_SHOW_PROFILE_BMP(64);
  GT_CHECK_CANDIDATES_SHOW_PROFILE_BMP(128);
  /* SIMD */
  fprintf(stderr,"  --> BPM.SIMD.unfolded64.distance   %2.3f (%2.3f%%)\t(%lu calls)\t(%2.3f ms/call)\n",
      GPROF_GET_TIMER(tprof[thread_id],GT_GPROF_BPM_SIMD128_UNFOLDED_DISTANCE),
      GPROF_TIME_PERCENTAGE(tprof[thread_id],GT_GPROF_BPM_SIMD128_UNFOLDED_DISTANCE,GT_GPROF_THREAD_ALL),
      GPROF_GET_COUNTER(tprof[thread_id],GT_GPROF_BPM_SIMD128_UNFOLDED_DISTANCE),
      GPROF_TIME_PER_CALL(tprof[thread_id],GT_GPROF_BPM_SIMD128_UNFOLDED_DISTANCE,GT_GPROF_BPM_SIMD128_UNFOLDED_DISTANCE));
  fprintf(stderr,"  --> BPM.SIMD.distance   %2.3f (%2.3f%%)\t(%lu calls)\t(%2.3f ms/call)\n",
      GPROF_GET_TIMER(tprof[thread_id],GT_GPROF_BPM_SIMD128_DISTANCE),
      GPROF_TIME_PERCENTAGE(tprof[thread_id],GT_GPROF_BPM_SIMD128_DISTANCE,GT_GPROF_THREAD_ALL),
      GPROF_GET_COUNTER(tprof[thread_id],GT_GPROF_BPM_SIMD128_DISTANCE),
      GPROF_TIME_PER_CALL(tprof[thread_id],GT_GPROF_BPM_SIMD128_DISTANCE,GT_GPROF_BPM_SIMD128_DISTANCE));
  /* */
  fprintf(stderr,"\n");
}
#define GT_CHECK_CANDIDATES_SHOW_GENERAL_PROFILE_BMP(VL) \
fprintf(stderr,"  --> BPM-" #VL ".distance   %2.3f (%2.3f%%)\t(%lu calls)\t(%2.3f ms/call)\n", \
    GPROF_GET_TIMER(tprof[0],GT_GPROF_BPM##VL##_DISTANCE), \
    GPROF_TIME_PERCENTAGE(tprof[0],GT_GPROF_BPM##VL##_DISTANCE,GT_GPROF_GLOBAL), \
    GPROF_GET_COUNTER(tprof[0],GT_GPROF_BPM##VL##_DISTANCE), \
    GPROF_TIME_PER_CALL(tprof[0],GT_GPROF_BPM##VL##_DISTANCE,GT_GPROF_BPM##VL##_DISTANCE)); \
fprintf(stderr,"  --> BPM-" #VL ".CUTOFF.distance   %2.3f (%2.3f%%)\t(%lu calls)\t(%2.3f ms/call)\n", \
    GPROF_GET_TIMER(tprof[0],GT_GPROF_BPM##VL##_CUTOFF_DISTANCE), \
    GPROF_TIME_PERCENTAGE(tprof[0],GT_GPROF_BPM##VL##_CUTOFF_DISTANCE,GT_GPROF_GLOBAL), \
    GPROF_GET_COUNTER(tprof[0],GT_GPROF_BPM##VL##_CUTOFF_DISTANCE), \
    GPROF_TIME_PER_CALL(tprof[0],GT_GPROF_BPM##VL##_CUTOFF_DISTANCE,GT_GPROF_BPM##VL##_CUTOFF_DISTANCE))
GT_INLINE void gt_check_candidates_show_general_profile() {
  GPROF_SUM_OVERLAP(tprof,parameters.num_threads);
  fprintf(stderr,"[GeneralProfile]\n");
  fprintf(stderr,"Total.Time %2.3f\t(%lu candidates)\n",
      GPROF_GET_TIMER(tprof[0],GT_GPROF_GLOBAL),
      GPROF_GET_COUNTER(tprof[0],GT_GPROF_NUM_CANDIDATES));
  fprintf(stderr,"  --> DP.distance   %2.3f (%2.3f%%)\t(%lu calls)\t(%2.3f ms/call)\n",
      GPROF_GET_TIMER(tprof[0],GT_GPROF_DP_DISTANCE),
      GPROF_TIME_PERCENTAGE(tprof[0],GT_GPROF_DP_DISTANCE,GT_GPROF_GLOBAL),
      GPROF_GET_COUNTER(tprof[0],GT_GPROF_DP_DISTANCE),
      GPROF_TIME_PER_CALL(tprof[0],GT_GPROF_DP_DISTANCE,GT_GPROF_DP_DISTANCE));
  /* STD Myers Base Implementation */
  fprintf(stderr,"  --> BPM.Base.distance   %2.3f (%2.3f%%)\t(%lu calls)\t(%2.3f ms/call)\n",
      GPROF_GET_TIMER(tprof[0],GT_GPROF_BASE_BPM64_DISTANCE),
      GPROF_TIME_PERCENTAGE(tprof[0],GT_GPROF_BASE_BPM64_DISTANCE,GT_GPROF_THREAD_ALL),
      GPROF_GET_COUNTER(tprof[0],GT_GPROF_BASE_BPM64_DISTANCE),
      GPROF_TIME_PER_CALL(tprof[0],GT_GPROF_BASE_BPM64_DISTANCE,GT_GPROF_BASE_BPM64_DISTANCE));
  fprintf(stderr,"  --> BPM.Base.CUTOFF.distance   %2.3f (%2.3f%%)\t(%lu calls)\t(%2.3f ms/call)\n",
      GPROF_GET_TIMER(tprof[0],GT_GPROF_BASE_BPM64_CUTOFF_DISTANCE),
      GPROF_TIME_PERCENTAGE(tprof[0],GT_GPROF_BASE_BPM64_CUTOFF_DISTANCE,GT_GPROF_THREAD_ALL),
      GPROF_GET_COUNTER(tprof[0],GT_GPROF_BASE_BPM64_CUTOFF_DISTANCE),
      GPROF_TIME_PER_CALL(tprof[0],GT_GPROF_BASE_BPM64_CUTOFF_DISTANCE,GT_GPROF_BASE_BPM64_CUTOFF_DISTANCE));
  /* Generic opt. vector impl.*/
  GT_CHECK_CANDIDATES_SHOW_GENERAL_PROFILE_BMP(8);
  GT_CHECK_CANDIDATES_SHOW_GENERAL_PROFILE_BMP(16);
  GT_CHECK_CANDIDATES_SHOW_GENERAL_PROFILE_BMP(32);
  GT_CHECK_CANDIDATES_SHOW_GENERAL_PROFILE_BMP(64);
  GT_CHECK_CANDIDATES_SHOW_GENERAL_PROFILE_BMP(128);
  /* SIMD */
  fprintf(stderr,"  --> BPM.SIMD.unfolded64.distance   %2.3f (%2.3f%%)\t(%lu calls)\t(%2.3f ms/call)\n",
      GPROF_GET_TIMER(tprof[0],GT_GPROF_BPM_SIMD128_UNFOLDED_DISTANCE),
      GPROF_TIME_PERCENTAGE(tprof[0],GT_GPROF_BPM_SIMD128_UNFOLDED_DISTANCE,GT_GPROF_THREAD_ALL),
      GPROF_GET_COUNTER(tprof[0],GT_GPROF_BPM_SIMD128_UNFOLDED_DISTANCE),
      GPROF_TIME_PER_CALL(tprof[0],GT_GPROF_BPM_SIMD128_UNFOLDED_DISTANCE,GT_GPROF_BPM_SIMD128_UNFOLDED_DISTANCE));
  fprintf(stderr,"  --> BPM.SIMD.distance   %2.3f (%2.3f%%)\t(%lu calls)\t(%2.3f ms/call)\n",
      GPROF_GET_TIMER(tprof[0],GT_GPROF_BPM_SIMD128_DISTANCE),
      GPROF_TIME_PERCENTAGE(tprof[0],GT_GPROF_BPM_SIMD128_DISTANCE,GT_GPROF_THREAD_ALL),
      GPROF_GET_COUNTER(tprof[0],GT_GPROF_BPM_SIMD128_DISTANCE),
      GPROF_TIME_PER_CALL(tprof[0],GT_GPROF_BPM_SIMD128_DISTANCE,GT_GPROF_BPM_SIMD128_DISTANCE));
}
/*
 * I/O Work Loop
 */
GT_INLINE void gt_check_candidates_parse_candidates_line(
      const char** const text_line,gt_candidates_job* const candidates_job) {
  // Reset candidates_job
  gt_candidates_job_clean(candidates_job);
  /*
   * Parse line:
   *   TCAGATGCATCG.....CGAACAG chr10:+:38880860:-1 chr10:+:42383932:-1
   */
  // Read
  gt_input_parse_field(text_line,TAB,candidates_job->reference_read);
  // Candidates
  while (!gt_input_parse_eol(text_line)) {
    gt_map* map = gt_map_new();
    gt_input_parse_field(text_line,COLON,map->seq_name); // chr10
    gt_input_parse_skip_chars(text_line,2); // Skip "+:"
    gt_input_parse_integer(text_line,(int64_t*)&map->position); // 38880860
    gt_input_parse_next_char(text_line); // ":"
    gt_input_parse_integer(text_line,(int64_t*)&map->gt_score); // -1
    gt_input_parse_field(text_line,TAB,NULL); // Skip the rest
    // Add the candidate
    gt_candidates_job_add_candidate(candidates_job,map);
  }
}
GT_INLINE gt_status gt_check_candidates_parse_candidates(
    gt_buffered_input_file* const buffered_input,gt_candidates_job* const candidates_job) {
  gt_status error_code;
  // Check the end_of_block. Reload buffer if needed
  if (gt_buffered_input_file_eob(buffered_input)) {
    if ((error_code=gt_buffered_input_file_reload(buffered_input,100))!=GT_INPUT_STATUS_OK) return error_code;
  }
  // Parse alignment
  gt_check_candidates_parse_candidates_line(gt_buffered_input_file_get_text_line(buffered_input),candidates_job);
  // Next record
  gt_buffered_input_file_skip_line(buffered_input);
  // OK
  return GT_INPUT_STATUS_OK;
}
void gt_check_candidates_read__write() {
  // Allocate profiles
  uint64_t i;
  tprof = gt_malloc(sizeof(gt_profile*)*parameters.num_threads);
  for (i=0;i<parameters.num_threads;++i) tprof[i] = GPROF_NEW(100);

  // Open I/O files
  gt_input_file* const input_file = gt_tools_open_input_file(parameters.name_input_file,GT_COMPRESSION_NONE);
  gt_output_file* const output_file = gt_tools_open_output_file(parameters.name_output_file,GT_COMPRESSION_NONE);

  // Open reference file
  gt_sequence_archive* sequence_archive =
      gt_tools_open_sequence_archive(parameters.name_gem_index_file,parameters.name_reference_file,true);

  // Parallel reading+process
  GPROF_START_TIMER(tprof[0],GT_GPROF_GLOBAL);
  #pragma omp parallel num_threads(parameters.num_threads)
  {
    // Thread ID
    const uint64_t thread_id = omp_get_thread_num();

    // Prepare IN/OUT buffers & printers
    gt_buffered_input_file* const buffered_input = gt_buffered_input_file_new(input_file);
    gt_buffered_output_file* const buffered_output = gt_buffered_output_file_new(output_file);
    gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_output);

    /*
     * READ + PROCCESS Loop
     */
    gt_vector* const buffer = gt_vector_new(1000,8);
    gt_bpm_pattern* const bpm8_pattern = gt_map_bpm_pattern_new();
    gt_bpm_pattern* const bpm16_pattern = gt_map_bpm_pattern_new();
    gt_bpm_pattern* const bpm32_pattern = gt_map_bpm_pattern_new();
    gt_bpm_pattern* const bpm64_pattern = gt_map_bpm_pattern_new();
    gt_bpm_pattern* const bpm128_pattern = gt_map_bpm_pattern_new();
    gt_candidates_job* const candidates_job = gt_candidates_job_new();
    while (gt_check_candidates_parse_candidates(buffered_input,candidates_job)) {

      // Do the job
      GPROF_START_TIMER(tprof[thread_id],GT_GPROF_THREAD_ALL);
      gt_check_candidates_align(thread_id,sequence_archive,candidates_job,buffer,
          bpm8_pattern,bpm16_pattern,bpm32_pattern,bpm64_pattern,bpm128_pattern);
      GPROF_STOP_TIMER(tprof[thread_id],GT_GPROF_THREAD_ALL);

    }

    // Clean
    gt_vector_delete(buffer);
    gt_map_bpm_pattern_delete(bpm8_pattern);
    gt_map_bpm_pattern_delete(bpm16_pattern);
    gt_map_bpm_pattern_delete(bpm32_pattern);
    gt_map_bpm_pattern_delete(bpm64_pattern);
    gt_map_bpm_pattern_delete(bpm128_pattern);
    gt_candidates_job_delete(candidates_job);
    gt_buffered_input_file_close(buffered_input);
    gt_buffered_output_file_close(buffered_output);
  }
  GPROF_STOP_TIMER(tprof[0],GT_GPROF_GLOBAL);

  // Global stats
  for (i=0;i<parameters.num_threads;++i) gt_check_candidates_show_profile(i);
  gt_check_candidates_show_general_profile();
  // Release archive & Clean
  gt_sequence_archive_delete(sequence_archive);
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
  for (i=0;i<parameters.num_threads;++i) GPROF_DELETE(tprof[i]);
  gt_free(tprof);
}
/*
 * Parse arguments
 */
gt_option gt_check_candidates_options[] = {
  /* I/O */
  { 'i', "input", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file>" , "" },
  { 'o', "output", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file>" , "" },
  { 'r', "reference", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file> (MultiFASTA/FASTA)" , "" },
  { 'I', "gem-index", GT_OPT_REQUIRED, GT_OPT_STRING, 2 , true, "<file> (GEM2-Index)" , "" },
  { 'e', "error", GT_OPT_REQUIRED, GT_OPT_FLOAT, 2 , true, "<error>" , "" },
  /* Misc */
  { 't', "threads", GT_OPT_REQUIRED, GT_OPT_INT, 3 , true, "" , "" },
  { 'h', "help", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , true, "" , "" },
  {  0, "", 0, 0, 0, false, "", ""}
};
char* gt_check_candidates_groups[] = {
  /*  0 */ "Null",
  /*  1 */ "Unclassified",
  /*  2 */ "I/O",
  /*  3 */ "Misc"
};
void parse_arguments(int argc,char** argv) {
  GT_OPTIONS_ITERATE_BEGIN(check_candidates,option) {
    /* I/O */
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    case 'r':
      parameters.name_reference_file = optarg;
      break;
    case 'I':
      parameters.name_gem_index_file = optarg;
      break;
    case 'e':
      parameters.error = atof(optarg);
      break;
    /* Misc */
    case 't': // threads
      parameters.num_threads = atol(optarg);
      break;
    case 'h': // help
      fprintf(stderr, "USE: ./gt.checkCandidates [ARGS]...\n");
      gt_options_fprint_menu(stderr,gt_check_candidates_options,gt_check_candidates_groups,false,false);
      exit(1);
    case '?':
    default:
      gt_fatal_error_msg("Option not recognized");
    }
  } GT_OPTIONS_ITERATE_END;
  /*
   * Parameters check
   */
  if (parameters.name_reference_file==NULL && parameters.name_gem_index_file==NULL) {
    gt_fatal_error_msg("Reference file required");
  }
}
/*
 * Main
 */
int main(int argc,char** argv) {
  // GT error handler
  gt_handle_error_signals();

  // Parsing command-line options
  parse_arguments(argc,argv);

  // Filter !!
  gt_check_candidates_read__write();

  return 0;
}

