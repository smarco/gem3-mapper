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
 *   Archive-Search Paired-End module encapsulating basic
 *   PE-search stages
 */

#include "archive/search/archive_search_pe_stages.h"
#include "archive/search/archive_search_pe.h"
#include "archive/search/archive_search_se_stepwise.h"
#include "archive/search/archive_search_pe_extend.h"
#include "archive/search/archive_select.h"
#include "archive/score/archive_score_se.h"
#include "matches/classify/matches_classify.h"
#include "matches/classify/paired_matches_classify.h"

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
#define ARCHIVE_SEARCH_PE_EXTENSION_MAX_READ_LENGTH 300

/*
 * PE Search Begin
 */
void archive_search_pe_begin(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_INIT,PROFILE_LEVEL);
  archive_search_end1->searched = false;
  archive_search_end2->searched = false;
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_INIT,PROFILE_LEVEL);
  // Next State
  archive_search_end1->pe_search_state = archive_search_pe_state_search_end1;
}
/*
 * PE Search individual ends
 */
void archive_search_pe_search_end1(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Parameters
  matches_t* const matches_end1 = paired_matches->matches_end1;
  // Full Search (End/1)
  archive_search_se_stepwise_finish_search(archive_search_end1,matches_end1);
  archive_search_end1->searched = true;
  // Test for extension of End/1 (Shortcut to avoid mapping end/2)
  if (archive_search_pe_use_shortcut_extension(archive_search_end1,archive_search_end2,matches_end1)) {
    // Extend End/1
    PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT_TOTAL);
    PROFILE_START(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT,PROFILE_LEVEL);
    archive_search_pe_extend_matches(
        archive_search_end1,archive_search_end2,paired_matches,paired_end2);
    PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT,PROFILE_LEVEL);
    // Check results of the extension
    if (paired_matches_get_num_maps(paired_matches) > 0) {
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
  archive_search_se_stepwise_finish_search(archive_search_end2,paired_matches->matches_end2);
  archive_search_end2->searched = true;
  // Next State
  archive_search_end1->pe_search_state = archive_search_pe_state_find_pairs;
}
/*
 * PE Search extension (aligning found matches nearby)
 */
bool archive_search_pe_use_shortcut_extension(
    archive_search_t* const archive_search_extended,
    archive_search_t* const archive_search_candidate,
    matches_t* const matches) {
  return false;
//  // Check extension enabled
//  search_parameters_t* const search_parameters = &archive_search_extended->search_parameters;
//  if (!search_parameters->search_paired_parameters.paired_end_extension_shortcut) return false;
//  // Check key-length
//  const uint64_t key_length = archive_search_candidate->approximate_search.pattern.key_length;
//  if (key_length==0 || key_length > ARCHIVE_SEARCH_PE_EXTENSION_MAX_READ_LENGTH) return false;
//  // Check the number of samples to derive the expected template size
//  if (!mapper_stats_template_length_is_reliable(archive_search_extended->mapper_stats)) return false;
//  // Test if the shortcut extension will provide a reliable mapping
//  return matches->metrics.mapq >= 30;
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
  if (key_length_end1==0 || key_length_end1>ARCHIVE_SEARCH_PE_EXTENSION_MAX_READ_LENGTH) return false;
  const uint64_t key_length_end2 = archive_search_end2->approximate_search.pattern.key_length;
  if (key_length_end2==0 || key_length_end2>ARCHIVE_SEARCH_PE_EXTENSION_MAX_READ_LENGTH) return false;
//  // Test suitability for extension
//  if (paired_matches->matches_end1->metrics.mapq < MAPQ_CONFIDENCE_SCORE_MIN ||
//      paired_matches->matches_end2->metrics.mapq < MAPQ_CONFIDENCE_SCORE_MIN) {
//    return true;
//  } else if (paired_matches_get_num_maps(paired_matches) > 0) {
//    paired_map_t* const primary_paired_map = paired_matches_get_primary_map(paired_matches);
//    return primary_paired_map->match_trace_end1->mapq_score < MAPQ_CONFIDENCE_SCORE_MIN ||
//           primary_paired_map->match_trace_end2->mapq_score < MAPQ_CONFIDENCE_SCORE_MIN;
//  }
//  return false;
  return true;
}
void archive_search_pe_recovery(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Paired-end recovery by extension (End/1)
  if (!paired_matches->matches_end1->matches_extended) {
    if (archive_search_pe_use_recovery_extension(archive_search_end1,archive_search_end2,paired_matches)) {
      // Extend End/1
      PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY_TOTAL);
      PROFILE_START(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY,PROFILE_LEVEL);
      archive_search_pe_extend_matches(
          archive_search_end1,archive_search_end2,paired_matches,paired_end2);
      PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY,PROFILE_LEVEL);
    }
    // Check Accuracy Reached (Quick abandon condition)
    if (paired_matches_classify_search_accomplished(paired_matches)) return;
  }
  // Paired-end recovery by extension (End/2)
  if (!paired_matches->matches_end2->matches_extended) {
    if (archive_search_pe_use_recovery_extension(archive_search_end1,archive_search_end2,paired_matches)) {
      // Extend End/2
      PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY_TOTAL);
      PROFILE_START(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY,PROFILE_LEVEL);
      archive_search_pe_extend_matches(
          archive_search_end1,archive_search_end2,paired_matches,paired_end1);
      PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY,PROFILE_LEVEL);
    }
  }
}
/*
 * PE Search cross-matching valid pairs
 */
void archive_search_pe_cross_pair_ends(
    search_parameters_t* const search_parameters,
    mapper_stats_t* const mapper_stats,
    paired_matches_t* const paired_matches,
    mm_allocator_t* const mm_allocator) {
  // Pair matches (Cross-link matches from both ends)
  const uint64_t num_matches_end1 = matches_get_num_match_traces(paired_matches->matches_end1);
  const uint64_t num_matches_end2 = matches_get_num_match_traces(paired_matches->matches_end2);
  if (num_matches_end1 > 0 && num_matches_end2 > 0) {
    PROFILE_START(GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH,PROFILE_LEVEL);
    paired_matches_clear(paired_matches,false); // Clean sheet
    paired_matches_find_pairs(paired_matches,search_parameters,mapper_stats,mm_allocator);
    PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH,PROFILE_LEVEL);
  }
}
void archive_search_pe_find_pairs(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Parameters
  search_parameters_t* const search_parameters = &archive_search_end1->search_parameters;
  mapper_stats_t* const mapper_stats = archive_search_end1->mapper_stats;
  mm_allocator_t* const mm_allocator = archive_search_end1->mm_allocator;
  // Find pairs
  archive_search_pe_cross_pair_ends(search_parameters,mapper_stats,paired_matches,mm_allocator);
  // Check accuracy reached
  if (!paired_matches_classify_search_accomplished(paired_matches)) {
    // Check paired-matches (for shortcut extension failure)
    if (!archive_search_end2->searched) {
      archive_search_end1->pe_search_state = archive_search_pe_state_search_end2;
      return;
    }
    // Recovery by extension
    archive_search_pe_recovery(archive_search_end1,archive_search_end2,paired_matches);
  }
  // Find discordant (if required)
  paired_matches_find_discordant_pairs(paired_matches,search_parameters);
  // PE Select Matches
  archive_select_pe_matches(&search_parameters->select_parameters,mapper_stats,paired_matches);
  // Next State
  archive_search_end1->pe_search_state = archive_search_pe_state_end;
}
/*
 * PE Search End
 */
void archive_search_pe_end(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Parameters
  search_parameters_t* const search_parameters = &archive_search_end1->search_parameters;
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
  matches_metrics_set_region_profile_metrics(&paired_matches->metrics,
      (region_profile_end1->avg_region_length+region_profile_end2->avg_region_length)/2,
      MAX(region_profile_end1->max_region_length,region_profile_end2->max_region_length),
      MAX(region_profile_end1->kmer_frequency,region_profile_end2->kmer_frequency));
  archive_search_end1->pe_search_state = archive_search_pe_state_end; // End of the workflow
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE) {
    archive_search_pe_print(stderr,archive_search_end1,archive_search_end2,paired_matches);
  }
}

