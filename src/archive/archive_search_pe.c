/*
 * PROJECT: GEMMapper
 * FILE: archive_search_pe.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive/archive_search_pe.h"
#include "archive/archive_search_se.h"
#include "archive/archive_search_se_stepwise.h"
#include "archive/archive_select.h"
#include "archive/archive_score_se.h"
#include "archive/archive_score_pe.h"
#include "archive/archive_check.h"
#include "approximate_search/approximate_search_control.h"
#include "approximate_search/approximate_search_verify_candidates.h"
#include "matches/matches_classify_logit_models.h"

/*
 * Debug
 */
#define DEBUG_ARCHIVE_SEARCH_PE GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * Constants
 */
#define ARCHIVE_SEARCH_PE_EXTENSION_MAX_READ_LENGTH 500

/*
 * Memory Injection (Support Data Structures)
 */
void archive_search_pe_inject_mm(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    mm_search_t* const mm_search) {
  // Search end/1
  archive_search_inject_mm_stack(archive_search_end1,mm_search->mm_stack);
  archive_search_inject_mapper_stats(archive_search_end1,mm_search->mapper_stats);
  archive_search_inject_interval_set(archive_search_end1,&mm_search->interval_set);
  archive_search_inject_text_collection(archive_search_end1,&mm_search->text_collection);
  archive_search_inject_filtering_candidates(archive_search_end1,
      &mm_search->filtering_candidates_forward_end1,
      &mm_search->filtering_candidates_reverse_end1,
      &mm_search->text_collection,mm_search->mm_stack);
  // Search end/2
  archive_search_inject_mm_stack(archive_search_end2,mm_search->mm_stack);
  archive_search_inject_mapper_stats(archive_search_end2,mm_search->mapper_stats);
  archive_search_inject_interval_set(archive_search_end2,&mm_search->interval_set);
  archive_search_inject_text_collection(archive_search_end2,&mm_search->text_collection);
  archive_search_inject_filtering_candidates(archive_search_end2,
      &mm_search->filtering_candidates_forward_end2,
      &mm_search->filtering_candidates_reverse_end2,
      &mm_search->text_collection,mm_search->mm_stack);
}
/*
 * PE Extension Control
 */
bool archive_search_pe_is_extension_feasible(archive_search_t* const archive_search) {
  // Check the number of samples to derive the expected template size
  return mapper_stats_template_length_is_reliable(archive_search->mapper_stats);
}
bool archive_search_pe_use_shortcut_extension(
    archive_search_t* const archive_search_extended,
    archive_search_t* const archive_search_candidate,
    matches_t* const matches) {
  // Check extension enabled
  search_parameters_t* const search_parameters = &archive_search_extended->search_parameters;
  if (!search_parameters->search_paired_parameters.paired_end_extension_shortcut) return false;
  // Check key-length
  const uint64_t key_length = archive_search_candidate->forward_search_state.pattern.key_length;
  if (key_length > ARCHIVE_SEARCH_PE_EXTENSION_MAX_READ_LENGTH) return false;
  // Check the number of samples to derive the expected template size
  if (!mapper_stats_template_length_is_reliable(archive_search_extended->mapper_stats)) return false;
  // Test if the shortcut extension will provide a reliable mapping
  return matches->metrics.mapq >= 30;
}
bool archive_search_pe_use_recovery_extension(
    archive_search_t* const archive_search_extended,
    archive_search_t* const archive_search_candidate,
    matches_t* const matches_extended,
    matches_t* const matches_candidate) {
  // Check extension enabled
  search_parameters_t* const search_parameters = &archive_search_extended->search_parameters;
  if (!search_parameters->search_paired_parameters.paired_end_extension_recovery) return false;
  // Check key-length
  const uint64_t key_length = archive_search_candidate->forward_search_state.pattern.key_length;
  if (key_length > ARCHIVE_SEARCH_PE_EXTENSION_MAX_READ_LENGTH) return false;
  // Check the number of samples to derive the expected template size
  if (!mapper_stats_template_length_is_reliable(archive_search_extended->mapper_stats)) return false;
  // Test suitability for extension
  return (matches_extended->metrics.mapq == 0 || matches_candidate->metrics.mapq == 0);
}
/*
 * PE Extension
 */
uint64_t archive_search_pe_extend_matches(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,
    const sequence_end_t candidate_end) {
  PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_TOTAL);
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES,PROFILE_LEVEL);
  // Parameters
  search_parameters_t* const search_parameters = &archive_search_end1->search_parameters;
  search_paired_parameters_t* const search_paired_parameters = &search_parameters->search_paired_parameters;
  mapper_stats_t* const mapper_stats = archive_search_end1->mapper_stats;
  // Extend in all possible concordant orientations
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  /*
   * Configure extension
   *   All extensions are done against the forward strand. If the candidate is in the reverse,
   *   the reverse-candidate region is generated reverse-complementing the forward. Thus, all
   *   the extension is done using @forward_search_state (matches contain all the strand needed info)
   */
  archive_search_t* candidate_archive_search;
  matches_t* extended_matches, *candidate_matches;
  if (candidate_end==paired_end2) {
    // Extend towards end/2
    candidate_archive_search = archive_search_end2;
    extended_matches = matches_end1;
    candidate_matches = matches_end2;
  } else {
    // Extend towards end/1
    candidate_archive_search = archive_search_end1;
    extended_matches = matches_end2;
    candidate_matches = matches_end1;
  }
  filtering_candidates_t* const filtering_candidates =
      candidate_archive_search->forward_search_state.filtering_candidates;
  pattern_t* const candidate_pattern = &candidate_archive_search->forward_search_state.pattern;
  uint64_t total_matches_found = 0;
  // Iterate over all matches of the extended end
//  const uint64_t first_stratum_distance = matches_get_match_trace_buffer(extended_matches)[0].distance;
  VECTOR_ITERATE(extended_matches->position_matches,extended_match,en,match_trace_t) {
//    if (extended_match->distance > first_stratum_distance) break; // Break beyond first stratum
    if (extended_match->type == match_type_extended) continue;    // Skip matches retrieve from extension
    if (search_paired_parameters->pair_orientation[pair_orientation_FR] == pair_relation_concordant) {
      // Extend (filter nearby region)
      PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_NUM_MATCHES);
      total_matches_found += approximate_search_verify_extend_candidate(filtering_candidates,
          candidate_pattern,extended_match,mapper_stats,paired_matches,candidate_end);
    }
  }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES,PROFILE_LEVEL);
  // (Re)Score Matches
  if (total_matches_found > 0) {
    archive_score_matches_se(candidate_archive_search,candidate_matches);
  }
  // Return
  return total_matches_found;
}
/*
 * Archive Search PE Continue Search
 */
void archive_search_pe_continue(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE,PROFILE_LEVEL);
  // Parameters
  search_parameters_t* const search_parameters = &archive_search_end1->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters_report;
  search_paired_parameters_t* const search_paired_parameters = &search_parameters->search_paired_parameters;
  mapper_stats_t* const mapper_stats = archive_search_end1->mapper_stats;
  matches_t* const matches_end1 = (paired_matches!=NULL) ? paired_matches->matches_end1 : NULL;
  matches_t* const matches_end2 = (paired_matches!=NULL) ? paired_matches->matches_end2 : NULL;
  // Callback (switch to proper search stage)
  while (archive_search_end1->pe_search_state != archive_search_pe_end) {
    switch (archive_search_end1->pe_search_state) {
      case archive_search_pe_begin: { // Beginning of the search (Init)
        PROFILE_START(GP_ARCHIVE_SEARCH_PE_INIT,PROFILE_LEVEL);
        archive_search_end1->pair_searched = false;
        archive_search_end1->pair_extended = false;
        archive_search_end1->pair_extended_shortcut = false;
        archive_search_end2->pair_searched = false;
        archive_search_end2->pair_extended = false;
        archive_search_reset(archive_search_end1); // Init (End/1)
        archive_search_reset(archive_search_end2); // Init (End/2)
        PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_INIT,PROFILE_LEVEL);
      }
      // No break
      case archive_search_pe_search_end1:
        // Full Search (End/1)
        archive_search_se_stepwise_finish_search(archive_search_end1,matches_end1);
        archive_select_se_matches(archive_search_end1,select_parameters,matches_end1);
        archive_search_end1->pair_searched = true;
        // Test for extension of End/1 (Shortcut to avoid mapping end/2)
        archive_search_end1->pair_extended =
            archive_search_pe_use_shortcut_extension(archive_search_end1,archive_search_end2,matches_end1);
        if (archive_search_end1->pair_extended) {
          // Extend End/1
          archive_search_end1->pair_extended_shortcut = true; // Debug
          PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT_TOTAL);
          PROFILE_START(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT,PROFILE_LEVEL);
          const uint64_t num_matches_found = archive_search_pe_extend_matches(
              archive_search_end1,archive_search_end2,paired_matches,paired_end2);
          PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT,PROFILE_LEVEL);
          // Check results of the extension
          if (num_matches_found > 0) {
            PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT_SUCCESS);
            archive_search_end1->pe_search_state = archive_search_pe_find_pairs;
            break; // Callback
          }
          /*
           * Extension failed because
           *  (1) Insert size is beyond the expected distribution (Any insert size filtering/search must be discarded)
           *  (2) Sensitivity of the candidates/end1 search is not enough
           */
        }
        // No break
      case archive_search_pe_search_end2:
        // Full Search (End/2)
        archive_search_se_stepwise_finish_search(archive_search_end2,matches_end2);
        archive_select_se_matches(archive_search_end2,select_parameters,matches_end2);
        archive_search_end2->pair_searched = true;
        // No break
      case archive_search_pe_recovery:
        // Paired-end recovery by extension
        if (!archive_search_end1->pair_extended) {
          if (archive_search_pe_use_recovery_extension(
              archive_search_end1,archive_search_end2,matches_end1,matches_end2)) {
            // Extend End/1
#ifdef GEM_PROFILE
            PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY_TOTAL);
            PROFILE_START(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY,PROFILE_LEVEL);
            const uint64_t num_matches_found =
#endif
            archive_search_pe_extend_matches(
                archive_search_end1,archive_search_end2,paired_matches,paired_end2);
            PROF_ADD_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY_SUCCESS,num_matches_found>0);
            archive_search_end1->pair_extended = true;
            PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY,PROFILE_LEVEL);
          }
        }
        if (!archive_search_end2->pair_extended) {
          if (archive_search_pe_use_recovery_extension(
              archive_search_end2,archive_search_end1,matches_end2,matches_end1)) {
            // Extend End/2
#ifdef GEM_PROFILE
            PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY_TOTAL);
            PROFILE_START(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY,PROFILE_LEVEL);
            const uint64_t num_matches_found =
#endif
            archive_search_pe_extend_matches(
                archive_search_end1,archive_search_end2,paired_matches,paired_end1);
            PROF_ADD_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY_SUCCESS,num_matches_found>0);
            archive_search_end2->pair_extended = true;
            PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY,PROFILE_LEVEL);
          }
        }
        // No break
      case archive_search_pe_find_pairs: {
        // Pair matches (Cross-link matches from both ends)
        const uint64_t num_matches_end1 = matches_get_num_match_traces(paired_matches->matches_end1);
        const uint64_t num_matches_end2 = matches_get_num_match_traces(paired_matches->matches_end2);
        if (num_matches_end1 > 0 && num_matches_end2 > 0) {
          PROFILE_START(GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH,PROFILE_LEVEL);
          paired_matches_clear(paired_matches,false); // Clean sheet
          paired_matches_find_pairs(paired_matches,search_paired_parameters,mapper_stats);
          paired_matches_find_discordant_pairs(paired_matches,search_paired_parameters); // Find discordant (if required)
          PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH,PROFILE_LEVEL);
        }
        // Check number of paired-matches
        if (paired_matches_get_num_maps(paired_matches) == 0) {
          if (!archive_search_end2->pair_searched) {
            archive_search_end1->pe_search_state = archive_search_pe_search_end2;
            break; // Callback
          }
        }
        // Set MCS
        paired_matches->max_complete_stratum =
            ((matches_end1->max_complete_stratum!=ALL) ? matches_end1->max_complete_stratum : 0) +
            ((matches_end2->max_complete_stratum!=ALL) ? matches_end2->max_complete_stratum : 0);
        archive_search_end1->pe_search_state = archive_search_pe_end; // End of the workflow
        gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE) {
          archive_search_pe_print(stderr,archive_search_end1,archive_search_end2,paired_matches);
        }
        break;
      }
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE,PROFILE_LEVEL);
}
/*
 * Paired-End Indexed Search (PE Online Approximate String Search)
 */
void archive_search_pe(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.PE\n");
    tab_fprintf(gem_log_get_stream(),"  => Tag %s\n",archive_search_end1->sequence.tag.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence/1 %s\n",archive_search_end1->sequence.read.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence/2 %s\n",archive_search_end2->sequence.read.buffer);
    tab_global_inc();
  }
  // Parameters
  search_parameters_t* const search_parameters = &archive_search_end1->search_parameters;
  // Init
  archive_search_end1->pe_search_state = archive_search_pe_begin;
  // archive_search_end2->pe_search_state = archive_search_pe_begin;
  // PE search
  archive_search_pe_continue(archive_search_end1,archive_search_end2,paired_matches);
  // PE select
  archive_select_pe_matches(archive_search_end1,archive_search_end2,
      &search_parameters->select_parameters_report,paired_matches);
  // PE Score (Select alignment-Model and process accordingly)
  archive_score_matches_pe(archive_search_end1,archive_search_end2,paired_matches);
  // PE Check matches
  if (search_parameters->check_type!=archive_check_nothing) {
    archive_check_pe_matches(
        archive_search_end1->archive,search_parameters->alignment_model,
        &search_parameters->swg_penalties,&archive_search_end1->sequence,
        &archive_search_end2->sequence,paired_matches,
        search_parameters->check_type,archive_search_end1->mm_stack);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE) {
    tab_global_inc();
    archive_search_pe_print(gem_log_get_stream(),archive_search_end1,archive_search_end2,paired_matches);
    tab_global_dec();
    tab_global_dec();
  }
}
/*
 * Compute Predictors
 */
void archive_search_pe_compute_predictors(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,
    matches_predictors_t* const predictors) {
  // Compute predictors
  approximate_search_metrics_t* const search_end1 = &archive_search_end1->forward_search_state.metrics;
  approximate_search_metrics_t* const search_end2 = &archive_search_end2->forward_search_state.metrics;
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  const uint64_t mcs_end1 = (matches_end1->max_complete_stratum!=ALL) ? matches_end1->max_complete_stratum : 0;
  const uint64_t mcs_end2 = (matches_end2->max_complete_stratum!=ALL) ? matches_end2->max_complete_stratum : 0;
  paired_matches_predictors_compute(paired_matches,predictors,search_end1,search_end2,mcs_end1,mcs_end2);
}
/*
 * Display
 */
void archive_search_pe_print(
    FILE* const stream,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  tab_fprintf(stream,"[GEM]>ArchiveSearch.PE\n");
  tab_global_inc();
  tab_fprintf(stream,"=> PE.Input\n");
  tab_fprintf(stream,"=>   Tag %s\n",archive_search_end1->sequence.tag.buffer);
  tab_fprintf(stream,"=>   Sequence/1 %s\n",archive_search_end1->sequence.read.buffer);
  tab_fprintf(stream,"=>   Sequence/2 %s\n",archive_search_end2->sequence.read.buffer);
  tab_fprintf(stream,"=>   Search.State %s\n",archive_search_pe_state_label[archive_search_end1->pe_search_state]);
  tab_fprintf(stream,"=> ASM.Parameters\n");
  tab_global_inc();
  search_parameters_print(stream,&archive_search_end1->search_parameters);
  tab_global_dec();
  tab_fprintf(stream,"=> End/1\n");
  tab_fprintf(stream,"  => Searched %s\n",archive_search_end1->pair_searched ? "yes" : "no");
  tab_fprintf(stream,"  => Extended %s (shortcut-extension=%s)\n",
      archive_search_end1->pair_extended ? "yes" : "no",
      archive_search_end1->pair_extended_shortcut ? "yes" : "no");
  tab_fprintf(stream,"=> End/2\n");
  tab_fprintf(stream,"  => Searched %s\n",archive_search_end2->pair_searched ? "yes" : "no");
  tab_fprintf(stream,"  => Extended %s\n",archive_search_end2->pair_extended ? "yes" : "no");
  if (!mapper_stats_template_length_is_reliable(archive_search_end1->mapper_stats)) {
    tab_fprintf(stream,"=> Template-length 'n/a'\n");
  } else {
    const double template_length_mean = mapper_stats_template_length_get_mean(archive_search_end1->mapper_stats);
    const double template_length_stddev = mapper_stats_template_length_get_stddev(archive_search_end1->mapper_stats);
    const uint64_t template_length_expected_max = mapper_stats_template_length_get_expected_max(archive_search_end1->mapper_stats);
    const uint64_t template_length_expected_min = mapper_stats_template_length_get_expected_min(archive_search_end1->mapper_stats);
    tab_fprintf(stream,"=> Template-length {min=%lu,max=%lu,mean=%2.1f,stddev==%2.1f}\n",
        template_length_expected_min,template_length_expected_max,template_length_mean,template_length_stddev);
  }
  // Archive Search
//  tab_fprintf(stream,"=> Archive.Search.End/1\n");
//  tab_global_inc();
//  archive_search_print(stream,archive_search_end1,NULL);
//  tab_global_dec();
//  tab_fprintf(stream,"=> Archive.Search.End/2\n");
//  tab_global_inc();
//  archive_search_print(stream,archive_search_end2,NULL);
//  tab_global_dec();
  // Paired Matches
  tab_fprintf(stream,"=> Paired.Matches\n");
  tab_global_inc();
  paired_matches_print(stream,paired_matches);
  tab_global_dec();
  tab_global_dec();
}

