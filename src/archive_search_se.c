/*
 * PROJECT: GEMMapper
 * FILE: archive_search_se.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_search_se.h"
#include "archive_select.h"
#include "archive_score.h"
#include "matches_classify.h"
#include "approximate_search_filtering_stages.h"

/*
 * Debug
 */
#define DEBUG_ARCHIVE_SEARCH_SE GEM_DEEP_DEBUG

/*
 * Setup
 */
GEM_INLINE void archive_search_single_end_configure(archive_search_t* const archive_search,mm_search_t* const mm_search) {
  archive_search_configure(archive_search,single_end,mm_search);
}

/*
 * SE Archive Search buiding blocks
 */
GEM_INLINE void archive_search_continue(
    archive_search_t* const archive_search,const bool verify_candidates,matches_t* const matches) {
  // Run the search (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  approximate_search(forward_asearch,matches); // Forward search
  if (archive_search->emulate_rc_search) {
    // Run the search (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    approximate_search(reverse_asearch,matches); // Reverse emulated-search
  }
}
GEM_INLINE void archive_search_generate_candidates(archive_search_t* const archive_search) {
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.Generate.Candidates\n");
    tab_fprintf(gem_log_get_stream(),"  => Tag %s\n",archive_search->sequence.tag.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence %s\n",archive_search->sequence.read.buffer);
    tab_global_inc();
  }
  // Reset initial values (Prepare pattern(s), instantiate parameters values, ...)
  archive_search_reset(archive_search);
  // Run the search (stop before filtering)
  PROF_START(GP_ARCHIVE_SEARCH_SE_GENERATE_CANDIDATES);
  archive_search_continue(archive_search,false,NULL);
  PROF_STOP(GP_ARCHIVE_SEARCH_SE_GENERATE_CANDIDATES);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE) { tab_global_dec(); }
}
GEM_INLINE void archive_search_verify_candidates(archive_search_t* const archive_search,matches_t* const matches) {
  PROF_START(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES);
  // Verify candidates (FORWARD)
  approximate_search_verify(&archive_search->forward_search_state,matches);
  if (archive_search->emulate_rc_search) {
    // Verify candidates (REVERSE)
    approximate_search_verify(&archive_search->reverse_search_state,matches);
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES);
}
GEM_INLINE void archive_search_finish_search(archive_search_t* const archive_search,matches_t* const matches) {
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.Finish.Search\n");
  }
  // Run the search up to the end
  PROF_START(GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH);
  archive_search_continue(archive_search,true,matches);
  PROF_STOP(GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE) {
    tab_global_inc();
    archive_search_print(gem_log_get_stream(),archive_search,matches);
    tab_global_dec();
    tab_global_dec();
  }
}
GEM_INLINE void archive_search_copy_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  PROF_START(GP_ARCHIVE_SEARCH_COPY_CANDIDATES);
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
  PROF_STOP(GP_ARCHIVE_SEARCH_COPY_CANDIDATES);
}
GEM_INLINE void archive_search_retrieve_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer,matches_t* const matches) {
  PROF_START(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES);
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
  PROF_STOP(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES);
}
/*
 * Single-End Indexed Search (SE Online Approximate String Search)
 */
GEM_INLINE void archive_search_single_end(archive_search_t* const archive_search,matches_t* const matches) {
  PROF_START(GP_ARCHIVE_SEARCH_SE);
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.SE\n");
    tab_fprintf(gem_log_get_stream(),"  => Tag %s\n",archive_search->sequence.tag.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence %s\n",archive_search->sequence.read.buffer);
    tab_global_inc();
  }
  // Reset initial values (Prepare pattern(s), instantiate parameters values, ...)
  archive_search_reset(archive_search);
  // Search the pattern(s)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  if (!archive_search->emulate_rc_search) {
    // Compute the full search
    approximate_search(forward_asearch,matches); // Forward search
  } else {
    // Configure search stage to stop at
    const as_parameters_t* const actual_parameters = &archive_search->as_parameters;
    const bool lower_max_difference =
        actual_parameters->complete_strata_after_best_nominal < forward_asearch->max_complete_error;
    if (lower_max_difference && archive_search->probe_strand) forward_asearch->stop_before_neighborhood_search = true;
    // Run the search (FORWARD)
    approximate_search(forward_asearch,matches); // Forward search
    // Check the number of matches & keep searching
    if (!forward_asearch->max_matches_reached) {
      // Keep on searching
      approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
      // Run the search (REVERSE)
      approximate_search(reverse_asearch,matches); // Reverse emulated-search
      // Resume forward search (if not completed before)
      if (forward_asearch->search_state != asearch_end && !forward_asearch->max_matches_reached) {
        approximate_search(forward_asearch,matches);
      }
    }
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_SE);
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE) {
    tab_global_inc();
    archive_search_print(gem_log_get_stream(),archive_search,matches);
    tab_global_dec();
    tab_global_dec();
  }
}
/*
 * Compute Predictors
 */
GEM_INLINE void archive_search_compute_predictors(
    archive_search_t* const archive_search,matches_t* const matches,
    matches_predictors_t* const predictors) {
  const uint64_t read_length = sequence_get_length(&archive_search->sequence);
  const swg_penalties_t* const swg_penalties = &archive_search->as_parameters.search_parameters->swg_penalties;
  const uint64_t max_region_length = archive_search_get_max_region_length(archive_search);
  const uint64_t num_zero_regions = archive_search_get_num_zero_regions(archive_search);
  const uint64_t proper_length = fm_index_get_proper_length(archive_search->archive->fm_index);
  matches_classify_compute_predictors(matches,predictors,
      swg_penalties,read_length,max_region_length,proper_length,
      matches->max_complete_stratum==ALL ? 0 : matches->max_complete_stratum,num_zero_regions);
}
/*
 * Display
 */
GEM_INLINE void archive_search_print(
    FILE* const stream,archive_search_t* const archive_search,matches_t* const matches) {
  tab_fprintf(stream,"[GEM]>ArchiveSearch.SE\n");
  tab_global_inc();
  tab_fprintf(stream,"=> Search.Forward\n");
  tab_global_inc();
  approximate_search_print(stream,&archive_search->forward_search_state);
  tab_global_dec();
  if (archive_search->emulate_rc_search) {
    tab_fprintf(stream,"=> Search.Reverse\n");
    tab_global_inc();
    approximate_search_print(stream,&archive_search->reverse_search_state);
    tab_global_dec();
  }
  if (matches!=NULL) {
    tab_fprintf(stream,"=> Matches.end\n");
    tab_global_inc();
    matches_print(stream,matches);
    tab_global_dec();
  }
  tab_global_dec();
}
