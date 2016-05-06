/*
 * PROJECT: GEMMapper
 * FILE: archive_search_se.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive/archive_search_se.h"
#include "archive/archive_select.h"
#include "archive/archive_score_se.h"
#include "archive/archive_check.h"
#include "approximate_search/approximate_search_filtering_adaptive.h"
#include "approximate_search/approximate_search_filtering_stages.h"
#include "filtering/region_profile_fixed.h"
#include "matches/matches_classify.h"

/*
 * Debug
 */
#define DEBUG_ARCHIVE_SEARCH_SE GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * Memory Injection (Support Data Structures)
 */
void archive_search_se_inject_mm(
    archive_search_t* const archive_search,
    mm_search_t* const mm_search) {
  archive_search_inject_mm_stack(archive_search,mm_search->mm_stack);
  archive_search_inject_mapper_stats(archive_search,mm_search->mapper_stats);
  archive_search_inject_interval_set(archive_search,&mm_search->interval_set);
  archive_search_inject_text_collection(archive_search,&mm_search->text_collection);
  archive_search_inject_filtering_candidates(archive_search,
      &mm_search->filtering_candidates_forward_end1,
      &mm_search->filtering_candidates_reverse_end1,
      &mm_search->text_collection,mm_search->mm_stack);
}
/*
 * Archive Search SE Continue
 */
void archive_search_se_continue(
    archive_search_t* const archive_search,
    matches_t* const matches) {
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
void archive_search_se(
    archive_search_t* const archive_search,
    matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.SE\n");
    tab_fprintf(gem_log_get_stream(),"  => Tag %s\n",archive_search->sequence.tag.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence %s\n",archive_search->sequence.read.buffer);
    tab_global_inc();
//    approximate_search_t* const search = &archive_search->forward_search_state;
//    region_profile_print_mappability(gem_log_get_stream(),search->archive->fm_index,
//        search->as_parameters->search_parameters->allowed_enc,search->pattern.key,
//        search->pattern.key_length,false,search->mm_stack);
  }
  // Reset initial values (Prepare pattern(s), instantiate parameters values, ...)
  archive_search_reset(archive_search);
  // Search the pattern(s)
  search_parameters_t* const search_parameters = &archive_search->search_parameters;
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  if (!archive_search->emulate_rc_search) {
    // Compute the full search
    approximate_search(forward_asearch,matches); // Forward search
  } else {
    // Configure search stage to stop at
    const bool lower_max_difference =
        search_parameters->complete_strata_after_best_nominal < forward_asearch->max_complete_error;
    if (lower_max_difference && archive_search->probe_strand) forward_asearch->stop_before_neighborhood_search = true;
    // Run the search (FORWARD)
    approximate_search(forward_asearch,matches); // Forward search
    // Keep on searching
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    // Run the search (REVERSE)
    approximate_search(reverse_asearch,matches); // Reverse emulated-search
    // Resume forward search (if not completed before)
    if (forward_asearch->search_stage != asearch_stage_end) {
      approximate_search(forward_asearch,matches);
    }
  }
  // Select Matches
  archive_select_se_matches(archive_search,&search_parameters->select_parameters_report,matches);
  // Select alignment-Model and process accordingly
  archive_score_matches_se(archive_search,matches);
  // Check matches
  if (search_parameters->check_type!=archive_check_nothing) {
    archive_check_se_matches(
        archive_search->archive,search_parameters->alignment_model,
        &search_parameters->swg_penalties,&archive_search->sequence,
        matches,search_parameters->check_type,archive_search->mm_stack);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE) {
    tab_global_inc();
    archive_search_se_print(gem_log_get_stream(),archive_search,matches);
    tab_global_dec();
    tab_global_dec();
  }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE,PROFILE_LEVEL);
}
/*
 * Compute Predictors
 */
void archive_search_se_compute_predictors(
    archive_search_t* const archive_search,
    matches_t* const matches,
    matches_predictors_t* const predictors) {
  const uint64_t max_complete_stratum = matches->max_complete_stratum==ALL ? 0 : matches->max_complete_stratum;
  matches_predictors_compute(matches,predictors,&archive_search->forward_search_state.metrics,max_complete_stratum);
}
/*
 * Display
 */
void archive_search_se_print(
    FILE* const stream,
    archive_search_t* const archive_search,
    matches_t* const matches) {
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
