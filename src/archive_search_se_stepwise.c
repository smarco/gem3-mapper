/*
 * PROJECT: GEMMapper
 * FILE: archive_search_se_stepwise.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_search_se_stepwise.h"
#include "archive_search_se.h"
#include "approximate_search_filtering_adaptive_stepwise.h"
#include "approximate_search_filtering_stages.h"
#include "filtering_candidates_bpm_buffer.h"

/*
 * Debug
 */
#define DEBUG_ARCHIVE_SEARCH_SE_STEPWISE GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * SE Archive Search Stepwise: Init Search
 */
GEM_INLINE void archive_search_se_stepwise_init_search(archive_search_t* const archive_search) {
  // Reset initial values (Prepare pattern(s), instantiate parameters values, ...)
  archive_search_reset(archive_search);
}
/*
 * SE Archive Search Stepwise: Region-Profile Generation
 */
GEM_INLINE void archive_search_se_stepwise_generate_region_profile_partition(
    archive_search_t* const archive_search) {
  // TODO TODO TODO TODO
}
GEM_INLINE void archive_search_se_stepwise_generate_region_profile_copy_partition(
    archive_search_t* const archive_search) {
  // TODO TODO TODO TODO
}
GEM_INLINE void archive_search_se_stepwise_generate_region_profile_retrieve_partition(
    archive_search_t* const archive_search) {
  // TODO TODO TODO TODO
}
/*
 * SE Archive Search Stepwise: Candidate Verification
 */
GEM_INLINE void archive_search_se_stepwise_generate_candidates(archive_search_t* const archive_search) {
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE.generate.candidates\n");
    tab_fprintf(gem_log_get_stream(),"  => Tag %s\n",archive_search->sequence.tag.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence %s\n",archive_search->sequence.read.buffer);
    tab_global_inc();
  }
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_GENERATE_CANDIDATES,PROFILE_LEVEL);
  // Run the search (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  approximate_search_filtering_adaptive_stepwise_generate_candidates(forward_asearch); // Forward search
  if (archive_search->emulate_rc_search) {
    // Run the search (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    approximate_search_filtering_adaptive_stepwise_generate_candidates(reverse_asearch); // Reverse emulated-search
  }
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_GENERATE_CANDIDATES,PROFILE_LEVEL);
}
GEM_INLINE void archive_search_se_stepwise_verify_candidates(
    archive_search_t* const archive_search,matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE.VerifyCPU.Candidates (stage=%s)\n",
        approximate_search_state_label[archive_search->forward_search_state.search_state]);
    tab_fprintf(gem_log_get_stream(),"  => Tag %s\n",archive_search->sequence.tag.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence %s\n",archive_search->sequence.read.buffer);
    tab_global_inc();
  }
  // Verify candidates (FORWARD)
  approximate_search_verify(&archive_search->forward_search_state,matches);
  if (archive_search->emulate_rc_search) {
    // Verify candidates (REVERSE)
    approximate_search_verify(&archive_search->reverse_search_state,matches);
  }
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES,PROFILE_LEVEL);
}
GEM_INLINE void archive_search_se_stepwise_copy_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  PROFILE_START(GP_ARCHIVE_SEARCH_COPY_CANDIDATES,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE.Copy.Candidates (stage=%s)\n",
        approximate_search_state_label[archive_search->forward_search_state.search_state]);
    tab_fprintf(gem_log_get_stream(),"  => Tag %s\n",archive_search->sequence.tag.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence %s\n",archive_search->sequence.read.buffer);
    tab_global_inc();
  }
  // Add candidates (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  forward_asearch->bpm_buffer_offset = bpm_gpu_buffer_get_num_candidates(bpm_gpu_buffer);
  forward_asearch->bpm_buffer_candidates = filtering_candidates_bpm_buffer_add(
      forward_asearch->filtering_candidates,&forward_asearch->pattern,bpm_gpu_buffer);
  if (archive_search->emulate_rc_search) {
    // Add candidates (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    reverse_asearch->bpm_buffer_offset = bpm_gpu_buffer_get_num_candidates(bpm_gpu_buffer);
    reverse_asearch->bpm_buffer_candidates = filtering_candidates_bpm_buffer_add(
        reverse_asearch->filtering_candidates,&reverse_asearch->pattern,bpm_gpu_buffer);
  }
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_COPY_CANDIDATES,PROFILE_LEVEL);
}
GEM_INLINE void archive_search_se_stepwise_retrieve_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer,
    matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE.Retrieve.Candidates (stage=%s)\n",
        approximate_search_state_label[archive_search->forward_search_state.search_state]);
    tab_fprintf(gem_log_get_stream(),"  => Tag %s\n",archive_search->sequence.tag.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence %s\n",archive_search->sequence.read.buffer);
    tab_global_inc();
  }
  // Verified candidates (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  approximate_search_verify_using_bpm_buffer(forward_asearch,bpm_gpu_buffer,forward_asearch->bpm_buffer_offset,
      forward_asearch->bpm_buffer_offset+forward_asearch->bpm_buffer_candidates,matches);
  if (archive_search->emulate_rc_search) {
    // Verified candidates (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    approximate_search_verify_using_bpm_buffer(reverse_asearch,bpm_gpu_buffer,reverse_asearch->bpm_buffer_offset,
        reverse_asearch->bpm_buffer_offset+reverse_asearch->bpm_buffer_candidates,matches);
  }
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES,PROFILE_LEVEL);
}
/*
 * SE Archive Search Stepwise: Finish Search
 */
GEM_INLINE void archive_search_se_stepwise_finish_search(archive_search_t* const archive_search,matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE.Finish (stage=%s)\n",
        approximate_search_state_label[archive_search->forward_search_state.search_state]);
    tab_fprintf(gem_log_get_stream(),"  => Tag %s\n",archive_search->sequence.tag.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence %s\n",archive_search->sequence.read.buffer);
    tab_global_inc();
  }
  // Run the search until the end (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  approximate_search_filtering_adaptive_stepwise_finish(forward_asearch,matches); // Forward search
  if (archive_search->emulate_rc_search) {
    // Run the search until the end (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    approximate_search_filtering_adaptive_stepwise_finish(reverse_asearch,matches); // Reverse emulated-search
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_global_inc();
    archive_search_se_print(gem_log_get_stream(),archive_search,matches);
    tab_global_dec();
    tab_global_dec();
  }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH,PROFILE_LEVEL);
}
