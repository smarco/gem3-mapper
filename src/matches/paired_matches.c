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
 */

#include "matches/paired_matches.h"
#include "matches/paired_matches_test.h"

/*
 * Constants
 */
#define PAIRED_MATCHES_INIT_COUNTERS 200
#define PAIRED_MATCHES_INIT_MATCHES   50
#define PAIRED_MATCHES_INIT_MATCHES   50

/*
 * Paired-Matches Classes
 */
const char* paired_matches_class_label[] =
{
    [0] = "unmapped",
    [1] = "perfect-tie",
    [2] = "tie",
    [3] = "mmap-d1",
    [4] = "mmap",
    [5] = "unique",
};

/*
 * Setup
 */
paired_matches_t* paired_matches_new(void) {
  // Alloc
  paired_matches_t* const paired_matches = mm_alloc(paired_matches_t);
  // State
  paired_matches->paired_matches_class = paired_matches_class_unmapped;
  paired_matches->max_complete_stratum = ALL;
  // Matches Counters
  paired_matches->counters = matches_counters_new();
  // Single-End Matches
  paired_matches->matches_end1 = matches_new();
  paired_matches->matches_end2 = matches_new();
  paired_matches->extended_matches = vector_new(PAIRED_MATCHES_INIT_MATCHES,match_trace_t*);
  // Paired-End Matches
  paired_matches->paired_maps = vector_new(PAIRED_MATCHES_INIT_MATCHES,paired_map_t);
  paired_matches->discordant_paired_maps = vector_new(PAIRED_MATCHES_INIT_MATCHES,paired_map_t);
  // Paired-End Metrics
  matches_metrics_init(&paired_matches->metrics);
  // Return
  return paired_matches;
}
void paired_matches_clear(paired_matches_t* const paired_matches,const bool clear_matches) {
  // State
  paired_matches->paired_matches_class = paired_matches_class_unmapped;
  paired_matches->max_complete_stratum = ALL;
  // Matches Counters
  matches_counters_clear(paired_matches->counters);
  // Single-End Matches
  if (clear_matches) {
    matches_clear(paired_matches->matches_end1);
    matches_clear(paired_matches->matches_end2);
  }
  // Paired-End Matches
  vector_clear(paired_matches->paired_maps);
  vector_clear(paired_matches->discordant_paired_maps);
  // Paired-End Metrics
  matches_metrics_init(&paired_matches->metrics);
}
void paired_matches_delete(paired_matches_t* const paired_matches) {
  // Matches Counters
  matches_counters_delete(paired_matches->counters);
  // Single-End Matches
  matches_delete(paired_matches->matches_end1);
  matches_delete(paired_matches->matches_end2);
  vector_delete(paired_matches->extended_matches);
  // Paired-End Matches
  vector_delete(paired_matches->paired_maps);
  vector_delete(paired_matches->discordant_paired_maps);
  // Delete handler
  mm_free(paired_matches);
}
/*
 * Accessors
 */
bool paired_matches_is_mapped(const paired_matches_t* const paired_matches) {
  return vector_get_used(paired_matches->paired_maps) > 0;
}
uint64_t paired_matches_get_num_maps(const paired_matches_t* const paired_matches) {
  return vector_get_used(paired_matches->paired_maps);
}
paired_map_t* paired_matches_get_maps(paired_matches_t* const paired_matches) {
  return vector_get_mem(paired_matches->paired_maps,paired_map_t);
}
uint64_t paired_matches_counters_get_count(paired_matches_t* const paired_matches,const uint64_t distance) {
  return matches_counters_get_count(paired_matches->counters,distance);
}
uint64_t paired_matches_counters_get_total_count(paired_matches_t* const paired_matches) {
  return matches_counters_get_total_count(paired_matches->counters);
}
uint64_t paired_matches_get_first_stratum_matches(paired_matches_t* const paired_matches) {
  const uint64_t min_edit_distance = matches_metrics_get_min_edit_distance(&paired_matches->metrics);
  return (min_edit_distance==UINT32_MAX) ? 0 : paired_matches_counters_get_count(paired_matches,min_edit_distance);
}
/*
 * Adding Paired-Matches
 */
void paired_matches_add(
    paired_matches_t* const paired_matches,
    match_trace_t* const match_trace_end1,
    match_trace_t* const match_trace_end2,
    const pair_relation_t pair_relation,
    const pair_orientation_t pair_orientation,
    const pair_layout_t pair_layout,
    const uint64_t template_length,
    const double template_length_sigma) {
  // Alloc paired match & add counters
  const uint64_t pair_event_distance = paired_map_compute_event_distance(match_trace_end1,match_trace_end2);
  const uint64_t pair_edit_distance = paired_map_compute_edit_distance(match_trace_end1,match_trace_end2);
  const int32_t pair_swg_score = paired_map_compute_swg_score(match_trace_end1,match_trace_end2);
  paired_map_t* paired_map;
  switch (pair_relation) {
    case pair_relation_concordant:
      vector_alloc_new(paired_matches->paired_maps,paired_map_t,paired_map);
      matches_counters_add(paired_matches->counters,pair_edit_distance,1); // Update counters
      paired_map->pair_relation = pair_relation_concordant;
      paired_matches_metrics_update(&paired_matches->metrics,
          pair_event_distance,pair_edit_distance,pair_swg_score,template_length_sigma);
      break;
    case pair_relation_discordant:
      vector_alloc_new(paired_matches->discordant_paired_maps,paired_map_t,paired_map);
      paired_map->pair_relation = pair_relation_discordant;
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Setup
  paired_map->match_trace_end1 = match_trace_end1;
  paired_map->match_trace_end2 = match_trace_end2;
  paired_map->pair_orientation = pair_orientation;
  paired_map->pair_layout = pair_layout;
  paired_map->template_length = template_length;
  paired_map->template_length_sigma = template_length_sigma;
  paired_map->index_position = MIN(
      match_trace_end1->match_alignment.match_position,
      match_trace_end2->match_alignment.match_position);
  paired_map->event_distance = pair_event_distance;
  paired_map->edit_distance = pair_edit_distance;
  paired_map->swg_score = pair_swg_score;
}
/*
 * Compute Pairs Relation
 */
pair_orientation_t paired_matches_compute_orientation(
    const match_trace_t* const match_trace_end1,
    const match_trace_t* const match_trace_end2) {
  if (match_trace_end1->strand == Forward) {
    if (match_trace_end2->strand == Forward) {
      return pair_orientation_FF;
    } else { // match_trace_end2->strand == Reverse
      if (match_trace_end1->text_position <= match_trace_end2->text_position) {
        return pair_orientation_FR;
      } else {
        return pair_orientation_RF;
      }
    }
  } else { // match_trace_end1->strand == Reverse
    if (match_trace_end2->strand == Reverse) {
      return pair_orientation_RR;
    } else { // match_trace_end2->strand == Forward
      if (match_trace_end2->text_position <= match_trace_end1->text_position) {
        return pair_orientation_FR;
      } else {
        return pair_orientation_RF;
      }
    }
  }
}
pair_layout_t paired_matches_compute_layout(
    const match_trace_t* const match_trace_end1,
    const match_trace_t* const match_trace_end2) {
  // Get matches location
  const uint64_t begin_position_1 = match_trace_end1->text_position;
  const uint64_t end_position_1 = begin_position_1 + match_trace_end1->match_alignment.effective_length;
  const uint64_t begin_position_2 = match_trace_end2->text_position;
  const uint64_t end_position_2 = begin_position_2 + match_trace_end2->match_alignment.effective_length;
  // Compute Pair-Layout
  if (end_position_1 < begin_position_2 || end_position_2 < begin_position_1) {
    return pair_layout_separate;
  } else if (begin_position_1 <= end_position_2) {
    // End-1 at the left
    if (begin_position_1 <= begin_position_2 && end_position_2 <= end_position_1) {
      return pair_layout_contain;
    } else {
      return pair_layout_overlap;
    }
  } else {
    // End-2 at the left
    if (begin_position_2 <= begin_position_1 && end_position_1 <= end_position_2) {
      return pair_layout_contain;
    } else {
      return pair_layout_overlap;
    }
  }
}
uint64_t paired_matches_compute_template_length(
    const match_trace_t* const match_trace_end1,
    const match_trace_t* const match_trace_end2,
    const pair_orientation_t pair_orientation,
    const pair_layout_t pair_layout) {
  // Get matches location
  const uint64_t effective_length_match_end1 = match_trace_end1->match_alignment.effective_length;
  const uint64_t effective_length_match_end2 = match_trace_end2->match_alignment.effective_length;
  uint64_t begin_position_1, end_position_1;
  uint64_t begin_position_2, end_position_2;
  if (match_trace_end1->text_position <= match_trace_end2->text_position) {
    begin_position_1 = match_trace_end1->text_position;
    end_position_1 = begin_position_1 + effective_length_match_end1;
    begin_position_2 = match_trace_end2->text_position;
    end_position_2 = begin_position_2 + effective_length_match_end2;
  } else {
    begin_position_1 = match_trace_end2->text_position;
    end_position_1 = begin_position_1 + effective_length_match_end2;
    begin_position_2 = match_trace_end1->text_position;
    end_position_2 = begin_position_2 + effective_length_match_end1;
  }
  /*
   * Compute Template-Length
   *   TLEN: signed observed Template LENgth. If all segments are mapped to the same reference, the
   *         unsigned observed template length equals the number of bases from the leftmost mapped base
   *         to the rightmost mapped base. The leftmost segment has a plus sign and the rightmost has a
   *         minus sign. The sign of segments in the middle is undefined. It is set as 0 for single-segment
   *         template or when the information is unavailable
   */
  switch (pair_orientation) {
    case pair_orientation_FR:
      switch (pair_layout) {
        case pair_layout_separate: return end_position_2-begin_position_1; // SAM spec compliant
        case pair_layout_overlap:  return end_position_2-begin_position_1; // SAM spec compliant
        case pair_layout_contain:  return MAX(effective_length_match_end1,effective_length_match_end2); // ???
      }
      break;
    case pair_orientation_RF:
      switch (pair_layout) {
        case pair_layout_separate: return end_position_2-begin_position_1; // SAM spec compliant
        case pair_layout_overlap:  return end_position_1-begin_position_2; // BWA-MEM compliant (not SAM spec)
        case pair_layout_contain:  return MAX(effective_length_match_end1,effective_length_match_end2); // ???
      }
      break;
    case pair_orientation_FF:
    case pair_orientation_RR:
      return MAX(end_position_1,end_position_2)-MIN(begin_position_1,begin_position_2); // SAM spec compliant
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  return 0;
}
pair_relation_t paired_matches_compute_relation(
    const search_paired_parameters_t* const parameters,
    mapper_stats_t* const mapper_stats,
    match_trace_t* const match_trace_end1,
    match_trace_t* const match_trace_end2,
    pair_orientation_t* const pair_orientation,
    pair_layout_t* const pair_layout,
    uint64_t* const template_length,
    double* const template_length_sigma) {
  pair_relation_t pair_relation = pair_relation_concordant; // Init
  // Compute BS-orientation
  if (match_trace_end1->bs_strand!=match_trace_end2->bs_strand) return pair_relation_invalid;
  // Compute orientation
  *pair_orientation = paired_matches_compute_orientation(match_trace_end1,match_trace_end2);
  if (parameters->pair_orientation[*pair_orientation]==pair_relation_invalid) return pair_relation_invalid;
  if (parameters->pair_orientation[*pair_orientation]==pair_relation_discordant) pair_relation = pair_relation_discordant;
  // Compute layout
  *pair_layout = paired_matches_compute_layout(match_trace_end1,match_trace_end2);
  if (parameters->pair_layout[*pair_layout]==pair_relation_invalid) return pair_relation_invalid;
  if (parameters->pair_layout[*pair_layout]==pair_relation_discordant) pair_relation = pair_relation_discordant;
  // Compute template-length & limits
  *template_length = paired_matches_compute_template_length(match_trace_end1,match_trace_end2,*pair_orientation,*pair_layout);
  if (*template_length < parameters->min_template_length || *template_length > parameters->max_template_length) {
    *template_length_sigma = MAX_TEMPLATE_LENGTH_SIGMAS;
    return pair_relation_invalid;
  }
  if (mapper_stats_template_length_is_reliable(mapper_stats)) {
    const uint64_t tl_expected_max = mapper_stats_template_length_get_expected_max(mapper_stats);
    const uint64_t tl_expected_min = mapper_stats_template_length_get_expected_min(mapper_stats);
    *template_length_sigma = mapper_stats_template_length_get_sigma_dev(mapper_stats,*template_length);
    if ((*template_length > tl_expected_max) || (*template_length < tl_expected_min)) {
      pair_relation = pair_relation_discordant;
    }
  } else {
    *template_length_sigma = 0.0; // Don't know (positive assumption)
  }
  return pair_relation;
}
/*
 * Cross pair matches
 */
void paired_matches_cross_pair(
    paired_matches_t* const paired_matches,
    search_parameters_t* const search_parameters,
    mapper_stats_t* const mapper_stats,
    match_trace_t* const match_trace_end1,
    match_trace_t* const match_trace_end2) {
  // Parameters
  search_paired_parameters_t* const search_paired_parameters = &search_parameters->search_paired_parameters;
  const pair_discordant_search_t discordant_search = search_paired_parameters->pair_discordant_search;
  // Check layout
  double template_length_sigma;
  pair_orientation_t pair_orientation;
  pair_layout_t pair_layout;
  uint64_t template_length;
  const pair_relation_t pair_relation = paired_matches_compute_relation(
      search_paired_parameters,mapper_stats,match_trace_end1,match_trace_end2,
      &pair_orientation,&pair_layout,&template_length,&template_length_sigma);
  switch (pair_relation) {
    case pair_relation_invalid:
      break;
    case pair_relation_discordant:
      if (discordant_search == pair_discordant_search_never) break;
      paired_matches_add(paired_matches,match_trace_end1,match_trace_end2,
          pair_relation_discordant,pair_orientation,pair_layout,template_length,template_length_sigma);
      break;
    case pair_relation_concordant: {
      paired_matches_add(paired_matches,match_trace_end1,match_trace_end2,
          pair_relation_concordant,pair_orientation,pair_layout,template_length,template_length_sigma);
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
}
/*
 * Find Pairs
 */
uint64_t paired_matches_find_pairs_locate(
    match_trace_t** const match_trace_sorted,
    const uint64_t num_match_traces,
    const char* const sequence_name) {
  uint64_t i;
  for (i=0;i<num_match_traces;++i) {
    if (gem_streq(sequence_name,match_trace_sorted[i]->sequence_name)) return i;
  }
  return num_match_traces;
}
void paired_matches_find_pairs_sorted(
    paired_matches_t* const paired_matches,
    search_parameters_t* const search_parameters,
    mapper_stats_t* const mapper_stats,
    match_trace_t** const matches_sorted_end1,
    const uint64_t num_matches_end1,
    match_trace_t** const matches_sorted_end2,
    const uint64_t num_matches_end2) {
  /*
   * Traverse all matches from the first end
   *   Might be interesting to implement a binary search of closest valid
   *   position and quick abandon after exploring feasible pairs (TODO)
   *   Shortcomings:
   *     Not compatible with discordant pairs
   *     Not really useful for regular mapping of a few matches on both sides
   */
  uint64_t position_end1;
  for (position_end1=0;position_end1<num_matches_end1;++position_end1) {
    // Fetch match-trace end/1
    match_trace_t* const match_trace_end1 = matches_sorted_end1[position_end1];
    const char* sequence_name = match_trace_end1->sequence_name;
    uint64_t position_end2 = paired_matches_find_pairs_locate(matches_sorted_end2,num_matches_end2,sequence_name);
    while (position_end2 < num_matches_end2) {
      // Fetch match-trace end/2
      match_trace_t* const match_trace_end2 = matches_sorted_end2[position_end2];
      if (!gem_streq(sequence_name,match_trace_end2->sequence_name)) break;
      // Cross pair
      paired_matches_cross_pair(
          paired_matches,search_parameters,mapper_stats,
          match_trace_end1,match_trace_end2);
      ++position_end2; // Next
    }
    // Check matches found
    const uint64_t max_searched_paired_matches = search_parameters->select_parameters.max_searched_paired_matches;
    const uint64_t num_paired_map = paired_matches_get_num_maps(paired_matches);
    if (num_paired_map >= max_searched_paired_matches) {
      paired_matches_sort_by_swg_score(paired_matches); // Sort
      if (paired_matches_test_accuracy_reached(paired_matches,search_parameters)) return; // Quick abandon condition
    }
  }
}
void paired_matches_find_pairs(
    paired_matches_t* const paired_matches,
    search_parameters_t* const search_parameters,
    mapper_stats_t* const mapper_stats,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  const uint64_t num_matches_end1 = matches_get_num_match_traces(matches_end1);
  const uint64_t num_matches_end2 = matches_get_num_match_traces(matches_end2);
  // Allocate matches placeholders
  mm_allocator_push_state(mm_allocator);
  match_trace_t** const matches_sorted_end1 = mm_allocator_calloc(mm_allocator,num_matches_end1,match_trace_t*,false);
  match_trace_t** const matches_sorted_end2 = mm_allocator_calloc(mm_allocator,num_matches_end2,match_trace_t*,false);
  // Create matches placeholders
  match_trace_t** match_trace_end1 = matches_get_match_traces(matches_end1);
  match_trace_t** match_trace_end2 = matches_get_match_traces(matches_end2);
  uint64_t i;
  for (i=0;i<num_matches_end1;++i) matches_sorted_end1[i] = match_trace_end1[i];
  for (i=0;i<num_matches_end2;++i) matches_sorted_end2[i] = match_trace_end2[i];
  // Sort both ends by (chr_name,position)
  matches_traces_sort_by_genomic_position(matches_sorted_end1,num_matches_end1);
  matches_traces_sort_by_genomic_position(matches_sorted_end2,num_matches_end2);
  // Traverse all matches and cross-pair
  paired_matches_find_pairs_sorted(
      paired_matches,search_parameters,mapper_stats,
      matches_sorted_end1,num_matches_end1,
      matches_sorted_end2,num_matches_end2);
  // Update number of matches candidates
  const uint64_t num_paired_maps = vector_get_used(paired_matches->paired_maps);
  matches_metrics_set_accepted_candidates(&paired_matches->metrics,num_paired_maps);
  // Pop allocator state
  mm_allocator_pop_state(mm_allocator);
}
void paired_matches_find_discordant_pairs(
    paired_matches_t* const paired_matches,
    search_parameters_t* const search_parameters) {
  // Parameters
  search_paired_parameters_t* const search_paired_parameters = &search_parameters->search_paired_parameters;
  // Check number of discordant pairs
  if (vector_get_used(paired_matches->discordant_paired_maps) > 0) return;
  // Check discordant mode
  switch (search_paired_parameters->pair_discordant_search) {
    case pair_discordant_search_never: return;
    case pair_discordant_search_only_if_no_concordant:
      if (vector_get_used(paired_matches->paired_maps) > 0) return;
    // No break
    case pair_discordant_search_always: {
      // Merge discordant paired-matches
      const uint64_t num_discordant_matches = vector_get_used(paired_matches->discordant_paired_maps);
      vector_reserve_additional(paired_matches->paired_maps,num_discordant_matches);
      paired_map_t* const concordant_map = vector_get_free_elm(paired_matches->paired_maps,paired_map_t);
      VECTOR_ITERATE_CONST(paired_matches->discordant_paired_maps,discordant_map,dn,paired_map_t) {
        // Add the discordant match
        concordant_map[dn] = *discordant_map;
        // Update counters
        matches_counters_add(paired_matches->counters,discordant_map->edit_distance,1);
        paired_matches_metrics_update(&paired_matches->metrics,discordant_map->event_distance,
            discordant_map->edit_distance,discordant_map->swg_score,discordant_map->template_length_sigma);
      }
      // Add to used
      vector_add_used(paired_matches->paired_maps,num_discordant_matches);
      // Clear discordant
      vector_clear(paired_matches->discordant_paired_maps);
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
}
/*
 * Filters
 */
void paired_matches_filter_by_mapq(
    paired_matches_t* const paired_matches,
    const uint8_t mapq_threshold) {
  paired_map_t* paired_map_out = paired_matches_get_maps(paired_matches);
  VECTOR_ITERATE(paired_matches->paired_maps,paired_map_in,n,paired_map_t) {
    if (paired_map_in->mapq_score >= mapq_threshold) {
      *paired_map_out = *paired_map_in;
      ++paired_map_out;
    }
  }
  vector_update_used(paired_matches->paired_maps,paired_map_out);
}
/*
 * Sort
 */
int paired_matches_cmp_swg_score(const paired_map_t* const a,const paired_map_t* const b) {
  // Compare relation (concordant first)
  if (a->pair_relation != b->pair_relation) return (b->pair_relation==pair_relation_concordant) ? 1 : -1;
  // Compare SWG-score
  const int distance_swg = (int)b->swg_score - (int)a->swg_score;
  if (distance_swg) return distance_swg;
  // Compare template-length-sigmas
  const int template_length_sigmas_diff = (int)a->template_length_sigma - (int)b->template_length_sigma;
  if (template_length_sigmas_diff) return template_length_sigmas_diff;
  if (a->template_length_sigma==0.0 || a->template_length_sigma==MAX_TEMPLATE_LENGTH_SIGMAS) {
    const int template_length_diff = (int)a->template_length - (int)b->template_length;
    if (template_length_diff) return template_length_diff;
  }
  // Compare layout (separated first)
  if (a->pair_layout != b->pair_layout) return (int)a->pair_relation - (int)b->pair_relation;
  // Untie using position (helps to stabilize & cmp results)
  if (a->index_position < b->index_position) return -1;
  if (a->index_position > b->index_position) return  1;
  return 0;
}
void paired_matches_sort_by_swg_score(paired_matches_t* const paired_matches) {
  // Sort global matches (match_trace_t) wrt swg-score
  qsort(vector_get_mem(paired_matches->paired_maps,paired_map_t),
      vector_get_used(paired_matches->paired_maps),sizeof(paired_map_t),
      (int (*)(const void *,const void *))paired_matches_cmp_swg_score);
}
/*
 * Display
 */
void paired_matches_print(FILE* const stream,paired_matches_t* const paired_matches) {
  tab_fprintf(stream,"[GEM]>Paired.Matches\n");
  tab_global_inc();
  tab_fprintf(stream,"=> Counters\t");
  matches_counters_print(stream,paired_matches->counters,paired_matches->max_complete_stratum);
  fprintf(stream,"\n");
  tab_fprintf(stream,"=> Paired.Maps %lu\n",vector_get_used(paired_matches->paired_maps));
  tab_fprintf(stream,"=> Discordant.Paired.Maps %lu\n",vector_get_used(paired_matches->discordant_paired_maps));
  tab_fprintf(stream,"=> Matches.end/1\n");
  tab_global_inc();
  matches_print(stream,paired_matches->matches_end1);
  tab_global_dec();
  tab_fprintf(stream,"=> Matches.end/2\n");
  tab_global_inc();
  matches_print(stream,paired_matches->matches_end2);
  tab_global_dec();
  tab_fprintf(stream,"=> Metrics.PE\n");
  tab_global_inc();
  matches_metrics_print(stream,&paired_matches->metrics);
  tab_global_dec();
  tab_global_dec();
}

