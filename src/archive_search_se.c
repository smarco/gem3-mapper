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
#include "approximate_search_filtering_adaptive.h"
#include "approximate_search_filtering_stages.h"

/*
 * Debug
 */
#define DEBUG_ARCHIVE_SEARCH_SE GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * Setup
 */
GEM_INLINE void archive_search_se_configure(
    archive_search_t* const archive_search,mm_search_t* const mm_search) {
  archive_search_configure(archive_search,single_end,mm_search);
}
/*
 * Archive Search SE Continue
 */
GEM_INLINE void archive_search_se_continue(archive_search_t* const archive_search,matches_t* const matches) {
  // Run the search (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  approximate_search(forward_asearch,matches); // Forward search
  if (archive_search->emulate_rc_search) {
    // Run the search (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    approximate_search(reverse_asearch,matches); // Reverse emulated-search
  }
}
/*
 * Single-End Indexed Search (SE Online Approximate String Search)
 */
GEM_INLINE void archive_search_se(archive_search_t* const archive_search,matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE,PROFILE_LEVEL);
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
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE) {
    tab_global_inc();
    archive_search_se_print(gem_log_get_stream(),archive_search,matches);
    tab_global_dec();
    tab_global_dec();
  }
}
/*
 * Compute Predictors
 */
GEM_INLINE void archive_search_se_compute_predictors(
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
GEM_INLINE void archive_search_se_print(
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
