/*
 * PROJECT: GEMMapper
 * FILE: archive_search_pe_stages.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive/archive_search_pe_stages.h"
#include "archive/archive_search_pe.h"
#include "archive/archive_search_se_stepwise.h"
#include "archive/archive_select.h"
#include "archive/archive_score_se.h"
#include "approximate_search/approximate_search_verify_candidates.h"

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
 * PE Extension Control
 */
bool archive_search_pe_use_shortcut_extension(
    archive_search_t* const archive_search_extended,
    archive_search_t* const archive_search_candidate,
    matches_t* const matches) {
  // Check extension enabled
  search_parameters_t* const search_parameters = &archive_search_extended->search_parameters;
  if (!search_parameters->search_paired_parameters.paired_end_extension_shortcut) return false;
  // Check key-length
  const uint64_t key_length = archive_search_candidate->approximate_search.pattern.key_length;
  if (key_length==0 || key_length > ARCHIVE_SEARCH_PE_EXTENSION_MAX_READ_LENGTH) return false;
  // Check the number of samples to derive the expected template size
  if (!mapper_stats_template_length_is_reliable(archive_search_extended->mapper_stats)) return false;
  // Test if the shortcut extension will provide a reliable mapping
  return matches->metrics.mapq >= 30;
}
bool archive_search_pe_use_recovery_extension(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Check extension enabled
  search_parameters_t* const search_parameters = &archive_search_end1->search_parameters;
  if (!search_parameters->search_paired_parameters.paired_end_extension_recovery) return false;
  // Check key-length
  const uint64_t key_length_end1 = archive_search_end1->approximate_search.pattern.key_length;
  if (key_length_end1==0 || key_length_end1 > ARCHIVE_SEARCH_PE_EXTENSION_MAX_READ_LENGTH) return false;
  const uint64_t key_length_end2 = archive_search_end2->approximate_search.pattern.key_length;
  if (key_length_end2==0 || key_length_end2 > ARCHIVE_SEARCH_PE_EXTENSION_MAX_READ_LENGTH) return false;
  // Check the number of samples to derive the expected template size
  if (!mapper_stats_template_length_is_reliable(archive_search_end1->mapper_stats)) return false;
  // Test suitability for extension
  if (paired_matches->matches_end1->metrics.mapq < MAPQ_CONFIDENCE_SCORE_MIN ||
      paired_matches->matches_end2->metrics.mapq < MAPQ_CONFIDENCE_SCORE_MIN) {
    return true;
  } else if (paired_matches_get_num_maps(paired_matches) > 0) {
    paired_map_t* const primary_paired_map = paired_matches_get_maps(paired_matches);
    return primary_paired_map->match_trace_end1->mapq_score < MAPQ_CONFIDENCE_SCORE_MIN ||
           primary_paired_map->match_trace_end2->mapq_score < MAPQ_CONFIDENCE_SCORE_MIN;
  }
  return false;
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
      candidate_archive_search->approximate_search.filtering_candidates;
  pattern_t* const candidate_pattern = &candidate_archive_search->approximate_search.pattern;
  uint64_t total_matches_found = 0;
  // Iterate over all matches of the extended end
  const uint64_t num_extended_match_traces = matches_get_num_match_traces(extended_matches);
  match_trace_t** const extended_match_traces = matches_get_match_traces(extended_matches);
  uint64_t i;
  for (i=0;i<num_extended_match_traces;++i) {
    if (extended_match_traces[i]->type == match_type_extended) continue; // Skip matches retrieved from extension
    if (search_paired_parameters->pair_orientation[pair_orientation_FR] == pair_relation_concordant) {
      // Extend (filter nearby region)
      PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_NUM_MATCHES);
      total_matches_found += approximate_search_verify_extend_candidate(filtering_candidates,
          candidate_pattern,extended_match_traces[i],mapper_stats,paired_matches,candidate_end);
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
 * Archive Search PE Stages
 */
void archive_search_pe_begin(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_INIT,PROFILE_LEVEL);
  archive_search_end1->pair_searched = false;
  archive_search_end1->pair_extended = false;
  archive_search_end1->pair_extended_shortcut = false;
  archive_search_end2->pair_searched = false;
  archive_search_end2->pair_extended = false;
  archive_search_reset(archive_search_end1); // Init (End/1)
  archive_search_reset(archive_search_end2); // Init (End/2)
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_INIT,PROFILE_LEVEL);
  // Next State
  archive_search_end1->pe_search_state = archive_search_pe_state_search_end1;
}
void archive_search_pe_search_end1(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Parameters
  matches_t* const matches_end1 = paired_matches->matches_end1;
  // Full Search (End/1)
  archive_search_se_stepwise_finish_search(archive_search_end1,matches_end1,true);
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
      // Next State
      archive_search_end1->pe_search_state = archive_search_pe_state_find_pairs;
      return;
    }
    /*
     * Extension failed because
     *  (1) Insert size is beyond the expected distribution (Any insert size filtering/search must be discarded)
     *  (2) Sensitivity of the candidates/end1 search is not enough
     */
  }
  // Next State
  archive_search_end1->pe_search_state = archive_search_pe_state_search_end2;
}
void archive_search_pe_search_end2(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Full Search (End/2)
  archive_search_se_stepwise_finish_search(archive_search_end2,paired_matches->matches_end2,true);
  archive_search_end2->pair_searched = true;
  // Next State
  archive_search_end1->pe_search_state = archive_search_pe_state_recovery;
}
void archive_search_pe_recovery(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Paired-end recovery by extension
  if (!archive_search_end1->pair_extended) {
    if (archive_search_pe_use_recovery_extension(archive_search_end1,archive_search_end2,paired_matches)) {
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
    if (archive_search_pe_use_recovery_extension(archive_search_end1,archive_search_end2,paired_matches)) {
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
  // Next State
  archive_search_end1->pe_search_state = archive_search_pe_state_find_pairs;
}
void archive_search_pe_find_pairs(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Parameters
  search_parameters_t* const search_parameters = &archive_search_end1->search_parameters;
  search_paired_parameters_t* const search_paired_parameters = &search_parameters->search_paired_parameters;
  mapper_stats_t* const mapper_stats = archive_search_end1->mapper_stats;
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
  const uint64_t num_matches = paired_matches_get_num_maps(paired_matches);
  if (num_matches == 0) {
    if (!archive_search_end2->pair_searched) {
      archive_search_end1->pe_search_state = archive_search_pe_state_search_end2;
      return;
    }
  }
  // PE Select Matches
  archive_select_pe_matches(archive_search_end1,archive_search_end2,
      &search_parameters->select_parameters_report,paired_matches);
  // Check for subdominant ends & force extension
  if ((!archive_search_end1->pair_extended || !archive_search_end1->pair_extended) &&
      archive_search_pe_use_recovery_extension(archive_search_end1,archive_search_end2,paired_matches)) {
    archive_search_end1->pe_search_state = archive_search_pe_state_recovery;
    return;
  }
  // Next State
  archive_search_end1->pe_search_state = archive_search_pe_state_end;
}
void archive_search_pe_end(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Parameters
  search_parameters_t* const search_parameters = &archive_search_end1->search_parameters;
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  // Set metrics
  matches_metrics_set_proper_length(&paired_matches->metrics,
      fm_index_get_proper_length(archive_search_end1->archive->fm_index));
  matches_metrics_set_read_length(&paired_matches->metrics,
      archive_search_end1->approximate_search.pattern.key_length +
      archive_search_end2->approximate_search.pattern.key_length);
  matches_metrics_set_swg_match_score(&paired_matches->metrics,
      search_parameters->swg_penalties.generic_match_score);
  region_profile_t* const region_profile_end1 = &archive_search_end1->approximate_search.region_profile;
  region_profile_t* const region_profile_end2 = &archive_search_end2->approximate_search.region_profile;
  matches_metrics_set_max_region_length(&paired_matches->metrics,
      MAX(region_profile_end1->max_region_length,region_profile_end2->max_region_length));
  matches_metrics_set_kmer_frequency(&paired_matches->metrics,
      MAX(region_profile_end1->kmer_frequency,region_profile_end2->kmer_frequency));
  // Set MCS
  paired_matches->max_complete_stratum =
      ((matches_end1->max_complete_stratum!=ALL) ? matches_end1->max_complete_stratum : 0) +
      ((matches_end2->max_complete_stratum!=ALL) ? matches_end2->max_complete_stratum : 0);
  archive_search_end1->pe_search_state = archive_search_pe_state_end; // End of the workflow
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE) {
    archive_search_pe_print(stderr,archive_search_end1,archive_search_end2,paired_matches);
  }
}

