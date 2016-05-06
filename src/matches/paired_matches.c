/*
 * PROJECT: GEMMapper
 * FILE: paired_matches.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "matches/paired_matches.h"

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
    [1] = "subdominant-end",
    [2] = "tie-d0",
    [3] = "tie-d1",
    [4] = "mmap",
    [5] = "unique",
    [6] = "high-quality-ends"
};

/*
 * Setup
 */
paired_matches_t* paired_matches_new() {
  // Alloc
  paired_matches_t* const paired_matches = mm_alloc(paired_matches_t);
  // State
  paired_matches->paired_matches_class = paired_matches_class_unmapped;
  paired_matches->max_complete_stratum = ALL;
  // Text Collection Buffer
  paired_matches->text_collection = NULL;
  // Matches Counters
  paired_matches->counters = matches_counters_new();
  // Single-End Matches
  paired_matches->matches_end1 = matches_new();
  paired_matches->matches_end2 = matches_new();
  // Paired-End Matches
  paired_matches->paired_maps = vector_new(PAIRED_MATCHES_INIT_MATCHES,paired_map_t);
  paired_matches->discordant_paired_maps = vector_new(PAIRED_MATCHES_INIT_MATCHES,paired_map_t);
  // Paired-End Metrics
  matches_metrics_init(&paired_matches->metrics);
  // Return
  return paired_matches;
}
void paired_matches_configure(paired_matches_t* const paired_matches,text_collection_t* const text_collection) {
  // Text Collection Buffer
  paired_matches->text_collection = text_collection;
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
  const uint64_t min_distance = matches_metrics_get_min_distance(&paired_matches->metrics);
  return (min_distance==UINT32_MAX) ? 0 : paired_matches_counters_get_count(paired_matches,min_distance);
}
match_trace_t* paired_map_get_match_end1(
    paired_matches_t* const paired_matches,
    const paired_map_t* const paired_map) {
  vector_t* const matches_end1 = paired_matches->matches_end1->position_matches;
  return vector_get_elm(matches_end1,paired_map->match_end1_offset,match_trace_t);
}
match_trace_t* paired_map_get_match_end2(
    paired_matches_t* const paired_matches,
    const paired_map_t* const paired_map) {
  vector_t* const matches_end2 = paired_matches->matches_end2->position_matches;
  return vector_get_elm(matches_end2,paired_map->match_end2_offset,match_trace_t);
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
  const uint64_t pair_distance = paired_map_compute_distance(match_trace_end1,match_trace_end2);
  const uint64_t pair_edit_distance = paired_map_compute_edit_distance(match_trace_end1,match_trace_end2);
  const int32_t pair_swg_score = paired_map_compute_swg_score(match_trace_end1,match_trace_end2);
  paired_map_t* paired_map;
  switch (pair_relation) {
    case pair_relation_concordant:
      vector_alloc_new(paired_matches->paired_maps,paired_map_t,paired_map);
      matches_counters_add(paired_matches->counters,pair_distance,1); // Update counters
      paired_map->pair_relation = pair_relation_concordant;
      paired_matches_metrics_update(&paired_matches->metrics,
          pair_distance,pair_edit_distance,pair_swg_score,template_length_sigma);
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
  vector_t* const position_matches_end1 = paired_matches->matches_end1->position_matches;
  vector_t* const position_matches_end2 = paired_matches->matches_end2->position_matches;
  paired_map->match_end1_offset = match_trace_end1 - vector_get_mem(position_matches_end1,match_trace_t);
  paired_map->match_end2_offset = match_trace_end2 - vector_get_mem(position_matches_end2,match_trace_t);
  paired_map->pair_orientation = pair_orientation;
  paired_map->pair_layout = pair_layout;
  paired_map->template_length = template_length;
  paired_map->template_length_sigma = template_length_sigma;
  paired_map->index_position = MIN(
      match_trace_end1->match_alignment.match_position,
      match_trace_end2->match_alignment.match_position);
  paired_map->distance = pair_distance;
  paired_map->edit_distance = pair_edit_distance;
  paired_map->swg_score = pair_swg_score;
}
/*
 * Finding Pairs
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
    paired_matches_t* const paired_matches,
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
match_trace_t* paired_matches_find_pairs_locate_by_sequence_name(
    matches_t* const matches,
    const char* const sequence_name) {
  const match_trace_t* const match_trace_sentinel =
      matches_get_match_trace_buffer(matches) + matches_get_num_match_traces(matches);
  match_trace_t* match_trace = matches_get_match_trace_buffer(matches);
  while (match_trace < match_trace_sentinel) {
    if (gem_streq(sequence_name,match_trace->sequence_name)) return match_trace;
    ++match_trace;
  }
  return NULL;
}
void paired_matches_find_pairs(
    paired_matches_t* const paired_matches,
    const search_paired_parameters_t* const search_paired_parameters,
    mapper_stats_t* const mapper_stats) {
  // Matches
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  // Sort both ends by (chr_name,position)
  matches_sort_by_sequence_name__position(matches_end1);
  matches_sort_by_sequence_name__position(matches_end2);
  // Traverse all matches from the first end
  uint64_t num_concordant_pair_matches = vector_get_used(paired_matches->paired_maps);
  uint64_t num_discordant_pair_matches = vector_get_used(paired_matches->discordant_paired_maps);
  match_trace_t* match_trace_end1 = vector_get_mem(matches_end1->position_matches,match_trace_t);
  const match_trace_t* const match_trace_end1_sentinel = match_trace_end1 + vector_get_used(matches_end1->position_matches);
  const match_trace_t* const match_trace_end2_sentinel =
      vector_get_mem(matches_end2->position_matches,match_trace_t) + vector_get_used(matches_end2->position_matches);
  const pair_discordant_search_t discordant_search = search_paired_parameters->pair_discordant_search;
  uint64_t template_length;
  double template_length_sigma;
  while (match_trace_end1 < match_trace_end1_sentinel) {
//    // Check number of total pair-matches found so far
//    if (num_concordant_pair_matches > max_paired_matches) break; // TODO Install quick-quit cond.
    // TODO Binary search of closest valid position and quick abandon after exploring feasible pairs
    // Traverse all possible pairs for @match_trace_end1
    const char* sequence_name = match_trace_end1->sequence_name;
    match_trace_t* match_trace_end2 = paired_matches_find_pairs_locate_by_sequence_name(matches_end2,sequence_name);
    if (match_trace_end2 != NULL) {
      while (match_trace_end2 < match_trace_end2_sentinel && gem_streq(sequence_name,match_trace_end2->sequence_name)) {
        pair_orientation_t pair_orientation;
        pair_layout_t pair_layout;
        const pair_relation_t pair_relation = paired_matches_compute_relation(
            paired_matches,search_paired_parameters,mapper_stats,match_trace_end1,match_trace_end2,
            &pair_orientation,&pair_layout,&template_length,&template_length_sigma);
        switch (pair_relation) {
          case pair_relation_invalid: break;
          case pair_relation_discordant:
            if (discordant_search == pair_discordant_search_never) break;
            paired_matches_add(paired_matches,match_trace_end1,match_trace_end2,pair_relation_discordant,
                pair_orientation,pair_layout,template_length,template_length_sigma);
            ++num_discordant_pair_matches;
            break;
          case pair_relation_concordant: {
            paired_matches_add(paired_matches,match_trace_end1,match_trace_end2,pair_relation_concordant,
                pair_orientation,pair_layout,template_length,template_length_sigma);
            ++num_concordant_pair_matches;
            break;
          }
          default:
            GEM_INVALID_CASE();
            break;
        }
        ++match_trace_end2;
      }
    }
    ++match_trace_end1;
  }
}
void paired_matches_find_discordant_pairs(
    paired_matches_t* const paired_matches,
    const search_paired_parameters_t* const search_paired_parameters) {
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
        matches_counters_add(paired_matches->counters,discordant_map->distance,1);
        paired_matches_metrics_update(&paired_matches->metrics,discordant_map->distance,
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
int paired_matches_cmp_distance(const paired_map_t* const a,const paired_map_t* const b) {
  // Compare relation (concordant first)
  if (a->pair_relation != b->pair_relation) return (b->pair_relation==pair_relation_concordant) ? 1 : -1;
  // Compare distance
  const int distance_diff = (int)a->distance - (int)b->distance;
  if (distance_diff) return distance_diff;
  // Compare SWG-score
  const int distance_swg = (int)b->swg_score - (int)a->swg_score;
  if (distance_swg) return distance_swg;
  // Compare Edit-distance
  const int distance_edit = (int)a->edit_distance - (int)b->edit_distance;
  if (distance_edit) return distance_edit;
  // Compare layout (separated first)
  if (a->pair_layout != b->pair_layout) return (int)a->pair_relation - (int)b->pair_relation;
  // Compare template-length-sigmas
  const int template_length_sigmas_diff = (int)a->template_length_sigma - (int)b->template_length_sigma;
  if (template_length_sigmas_diff) return template_length_sigmas_diff;
  // Untie using position (helps to stabilize & cmp results)
  return (int)a->index_position - (int)b->index_position;
}
void paired_matches_sort_by_distance(paired_matches_t* const paired_matches) {
  // Sort global matches (match_trace_t) wrt distance
  qsort(vector_get_mem(paired_matches->paired_maps,paired_map_t),
      vector_get_used(paired_matches->paired_maps),sizeof(paired_map_t),
      (int (*)(const void *,const void *))paired_matches_cmp_distance);
}
int paired_matches_cmp_mapq_score(const paired_map_t* const a,const paired_map_t* const b) {
  return b->mapq_score - a->mapq_score;
}
void paired_matches_sort_by_mapq_score(paired_matches_t* const paired_matches) {
  // Sort global matches (match_trace_t) wrt distance
  qsort(vector_get_mem(paired_matches->paired_maps,paired_map_t),
      vector_get_used(paired_matches->paired_maps),sizeof(paired_map_t),
      (int (*)(const void *,const void *))paired_matches_cmp_mapq_score);
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

