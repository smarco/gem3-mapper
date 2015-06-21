/*
 * PROJECT: GEMMapper
 * FILE: paired_matches.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "paired_matches.h"

/*
 * Constants
 */
#define PAIRED_MATCHES_INIT_COUNTERS 200
#define PAIRED_MATCHES_INIT_MATCHES   50
#define PAIRED_MATCHES_INIT_MATCHES   50

/*
 * Setup
 */
GEM_INLINE paired_matches_t* paired_matches_new() {
  // Alloc
  paired_matches_t* const paired_matches = mm_alloc(paired_matches_t);
  // State
  paired_matches->max_complete_stratum = ALL;
  // Text Collection Buffer
  paired_matches->text_collection = NULL;
  // Matches Counters
  paired_matches->counters = vector_new(PAIRED_MATCHES_INIT_COUNTERS,uint64_t);
  paired_matches->discordant_counters = vector_new(PAIRED_MATCHES_INIT_COUNTERS,uint64_t);
  // Single-End Matches
  paired_matches->matches_end1 = matches_new();
  paired_matches->matches_end2 = matches_new();
  // Paired-End Matches
  paired_matches->matches = vector_new(PAIRED_MATCHES_INIT_MATCHES,paired_match_t);
  paired_matches->discordant_matches = vector_new(PAIRED_MATCHES_INIT_MATCHES,paired_match_t);
  // Return
  return paired_matches;
}
GEM_INLINE void paired_matches_configure(paired_matches_t* const paired_matches,text_collection_t* const text_collection) {
  // Text Collection Buffer
  paired_matches->text_collection = text_collection;
}
GEM_INLINE void paired_matches_clear(paired_matches_t* const paired_matches) {
  // State
  paired_matches->max_complete_stratum = ALL;
  // Matches Counters
  vector_clear(paired_matches->counters);
  vector_clear(paired_matches->discordant_counters);
  // Single-End Matches
  matches_clear(paired_matches->matches_end1);
  matches_clear(paired_matches->matches_end2);
  // Paired-End Matches
  vector_clear(paired_matches->matches);
  vector_clear(paired_matches->discordant_matches);
}
GEM_INLINE void paired_matches_delete(paired_matches_t* const paired_matches) {
  // Matches Counters
  vector_delete(paired_matches->counters);
  vector_delete(paired_matches->discordant_counters);
  // Single-End Matches
  matches_delete(paired_matches->matches_end1);
  matches_delete(paired_matches->matches_end2);
  // Paired-End Matches
  vector_delete(paired_matches->matches);
  vector_delete(paired_matches->discordant_matches);
  // Delete handler
  mm_free(paired_matches);
}
/*
 * Accessors
 */
GEM_INLINE bool paired_matches_is_mapped(paired_matches_t* const paired_matches) {
  return vector_get_used(paired_matches->matches) > 0;
}
/*
 * Counters
 */
GEM_INLINE void paired_matches_counters_add(vector_t* const counters,const uint64_t pair_distance) {
  // Reserve Memory
  if (pair_distance >= vector_get_used(counters)) {
    vector_reserve(counters,pair_distance+1,true);
    vector_set_used(counters,pair_distance+1);
  }
  ++(*vector_get_elm(counters,pair_distance,uint64_t));
}
/*
 * Adding Paired-Matches
 */
GEM_INLINE void paired_matches_add(
    paired_matches_t* const paired_matches,match_trace_t* const match_trace_end1,
    match_trace_t* const match_trace_end2,const pair_orientation_t pair_orientation,const uint64_t template_length) {
  // Alloc paired match & add counters
  const uint64_t pair_distance = paired_match_calculate_distance(match_trace_end1,match_trace_end2);
  paired_match_t* paired_match;
  switch (pair_orientation) {
    case pair_orientation_concordant:
      // Alloc
      vector_alloc_new(paired_matches->matches,paired_match_t,paired_match);
      paired_matches_counters_add(paired_matches->counters,pair_distance); // Update counters
      paired_match->pair_orientation = pair_orientation_concordant;
      break;
    case pair_orientation_discordant:
      // Alloc
      vector_alloc_new(paired_matches->discordant_matches,paired_match_t,paired_match);
      paired_matches_counters_add(paired_matches->discordant_counters,pair_distance); // Update counters
      paired_match->pair_orientation = pair_orientation_discordant;
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Setup
  paired_match->match_end1 = match_trace_end1;
  paired_match->match_end2 = match_trace_end2;
  paired_match->template_length = template_length;
  paired_match->distance = pair_distance;
}
/*
 * Finding Pairs
 */
GEM_INLINE bool paired_matches_get_layout(
    paired_matches_t* const paired_matches,const match_trace_t* const match_trace_end1,
    const match_trace_t* const match_trace_end2,search_parameters_t* const search_parameters,
    mapper_stats_t* const mapper_stats,uint64_t* const template_length,
    pair_layout_t* const pair_layout,bool* const concordant_template_length) {
  // Get template observed length
  const uint64_t begin_position_1 = match_trace_end1->text_position;
  const uint64_t end_position_1 = begin_position_1 + match_trace_end1->match_alignment.effective_length;
  const uint64_t begin_position_2 = match_trace_end2->text_position;
  const uint64_t end_position_2 = begin_position_2 + match_trace_end2->match_alignment.effective_length;
  const uint64_t template_observed_length = paired_match_get_template_observed_length(
      begin_position_1,end_position_1,begin_position_2,end_position_2);
  *template_length = template_observed_length; // Assign
  // Check Template-Length limits
  if (template_observed_length < search_parameters->min_template_length) {
    return false;
  } else if (template_observed_length > search_parameters->max_template_length) {
    return false;
  }
  // Compute Template-Length constraints
  if (mapper_stats_template_length_estimation_within_ci(mapper_stats,MS_TEMPLATE_LENGTH_DEFAULT_MOE)) {
    *concordant_template_length =
        (template_observed_length <= mapper_stats_template_length_get_expected_max(mapper_stats)) &&
        (template_observed_length >= mapper_stats_template_length_get_expected_min(mapper_stats));
  } else {
    *concordant_template_length = true;
  }
  // Compute Pair-Layout constraints
  if (end_position_1 < begin_position_2 || end_position_2 < begin_position_1) {
    *pair_layout = pair_layout_separate;
  } else if (begin_position_1 <= end_position_2) {
    // End-1 at the left
    if (begin_position_1 <= begin_position_2 && end_position_2 <= end_position_1) {
      *pair_layout = pair_layout_contain;
    } else if (end_position_1 >= begin_position_2) {
      *pair_layout = pair_layout_overlap;
    } else {
      *pair_layout = pair_layout_dovetail;
    }
  } else {
    // End-2 at the left
    if (begin_position_2 <= begin_position_1 && end_position_1 <= end_position_2) {
      *pair_layout = pair_layout_contain;
    } else if (end_position_2 >= begin_position_1) {
      *pair_layout = pair_layout_overlap;
    } else {
      *pair_layout = pair_layout_dovetail;
    }
  }
  // Return allowed layout
  return true;
}
GEM_INLINE bool paired_matches_compute_layout(
    paired_matches_t* const paired_matches,match_trace_t* const match_trace_end1,
    match_trace_t* const match_trace_end2,search_parameters_t* const search_parameters,
    mapper_stats_t* const mapper_stats,uint64_t* const template_length,
    bool* const concordant_layout,bool* const concordant_template_length) {
  pair_layout_t pair_layout;
  const bool allowed_layout = paired_matches_get_layout(paired_matches,
      match_trace_end1,match_trace_end2,search_parameters,
      mapper_stats,template_length,&pair_layout,concordant_template_length);
  if (!allowed_layout) return false;
  switch (pair_layout) {
    case pair_layout_separate:
      *concordant_layout = search_parameters->pair_layout_separate;
      break;
    case pair_layout_overlap:
      *concordant_layout = search_parameters->pair_layout_overlap;
      break;
    case pair_layout_contain:
      *concordant_layout = search_parameters->pair_layout_contain;
      break;
    case pair_layout_dovetail:
      *concordant_layout = search_parameters->pair_layout_dovetail;
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  return true;
}
GEM_INLINE pair_orientation_t paired_matches_get_orientation(
    const match_trace_t* const match_trace_end1,const match_trace_t* const match_trace_end2,
    search_parameters_t* const search_parameters) {
  if (match_trace_end1->strand == Forward) {
    if (match_trace_end2->strand == Forward) {
      return search_parameters->pair_orientation_FF;
    } else { // match_trace_end2->strand == Reverse
      if (match_trace_end1->text_position <= match_trace_end2->text_position) {
        return search_parameters->pair_orientation_FR;
      } else {
        return search_parameters->pair_orientation_RF;
      }
    }
  } else { // match_trace_end1->strand == Reverse
    if (match_trace_end2->strand == Reverse) {
      return search_parameters->pair_orientation_RR;
    } else { // match_trace_end2->strand == Forward
      if (match_trace_end2->text_position <= match_trace_end1->text_position) {
        return search_parameters->pair_orientation_FR;
      } else {
        return search_parameters->pair_orientation_RF;
      }
    }
  }
}
GEM_INLINE match_trace_t* paired_matches_find_pairs_locate_by_sequence_name(
    matches_t* const matches,const char* const sequence_name) {
  const match_trace_t* const match_trace_sentinel =
      matches_get_match_traces(matches) + matches_get_num_match_traces(matches);
  match_trace_t* match_trace = matches_get_match_traces(matches);
  while (match_trace < match_trace_sentinel) {
    if (gem_streq(sequence_name,match_trace->sequence_name)) return match_trace;
    ++match_trace;
  }
  return NULL;
}
GEM_INLINE void paired_matches_pair_match_with_mates(
    paired_matches_t* const paired_matches,search_parameters_t* const search_parameters,
    mapper_stats_t* const mapper_stats,const pair_orientation_t pair_orientation,
    match_trace_t* const match_trace,const sequence_end_t mate_end,
    match_trace_t* const mates_array,const uint64_t num_mates_trace) {
  // Traverse all mates
  uint64_t mate_pos, template_length;
  bool concordant_template_length, concordant_layout;
  for (mate_pos=0;mate_pos<num_mates_trace;++mate_pos) {
    match_trace_t* const mate_trace = mates_array+mate_pos;
    // Check layout & add
    const bool allowed_layout = paired_matches_compute_layout(paired_matches,match_trace,mate_trace,
        search_parameters,mapper_stats,&template_length,&concordant_layout,&concordant_template_length);
    if (allowed_layout) {
      const pair_orientation_t pair_orientation = concordant_layout ? pair_orientation_concordant : pair_orientation_discordant;
      if (mate_end==paired_end2) {
        paired_matches_add(paired_matches,match_trace,mate_trace,pair_orientation,template_length);
      } else {
        paired_matches_add(paired_matches,mate_trace,match_trace,pair_orientation,template_length);
      }
    }
  }
}
GEM_INLINE void paired_matches_find_pairs(
    paired_matches_t* const paired_matches,search_parameters_t* const search_parameters,
    mapper_stats_t* const mapper_stats) {
  PROF_START(GP_PAIRED_MATCHES_FIND_PAIRS);
  // Matches
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  // Sort both ends by (chr_name,position)
  matches_sort_by_sequence_name__position(matches_end1);
  matches_sort_by_sequence_name__position(matches_end2);
  // Traverse all matches from the first end
  const uint64_t max_search_matches = search_parameters->max_search_matches;
  uint64_t num_concordant_pair_matches = vector_get_used(paired_matches->matches);
  uint64_t num_discordant_pair_matches = vector_get_used(paired_matches->discordant_matches);
  match_trace_t* match_trace_end1 = vector_get_mem(matches_end1->position_matches,match_trace_t);
  const match_trace_t* const match_trace_end1_sentinel = match_trace_end1 + vector_get_used(matches_end1->position_matches);
  const match_trace_t* const match_trace_end2_sentinel =
      vector_get_mem(matches_end2->position_matches,match_trace_t) + vector_get_used(matches_end2->position_matches);
  const pair_discordant_search_t discordant_search = search_parameters->pair_discordant_search;
  uint64_t template_length;
  bool concordant_template_length, allowed_layout, concordant_layout;
  while (match_trace_end1 < match_trace_end1_sentinel) {
    // Check number of total pair-matches found so far
    if (num_concordant_pair_matches > max_search_matches || num_discordant_pair_matches > max_search_matches) break;
    // TODO Binary search of closest valid position and quick abandon after exploring feasible pairs
    // Traverse all possible pairs for @match_trace_end1
    const char* sequence_name = match_trace_end1->sequence_name;
    match_trace_t* match_trace_end2 = paired_matches_find_pairs_locate_by_sequence_name(matches_end2,sequence_name);
    if (match_trace_end2 != NULL) {
      while (match_trace_end2 < match_trace_end2_sentinel && gem_streq(sequence_name,match_trace_end2->sequence_name)) {
        // Get orientation
        const pair_orientation_t pair_orientation =
            paired_matches_get_orientation(match_trace_end1,match_trace_end2,search_parameters);
        switch (pair_orientation) {
          case pair_orientation_invalid: break;
          case pair_orientation_discordant:
            if (discordant_search == pair_discordant_search_never) break;
            allowed_layout = paired_matches_compute_layout(paired_matches,match_trace_end1,match_trace_end2,
                search_parameters,mapper_stats,&template_length,&concordant_layout,&concordant_template_length);
            if (!allowed_layout) break;
            paired_matches_add(paired_matches,match_trace_end1,match_trace_end2,pair_orientation_discordant,template_length);
            ++num_discordant_pair_matches;
            break;
          case pair_orientation_concordant: {
            allowed_layout = paired_matches_compute_layout(paired_matches,match_trace_end1,match_trace_end2,
                search_parameters,mapper_stats,&template_length,&concordant_layout,&concordant_template_length);
            if (!allowed_layout) break;
            if (concordant_layout && concordant_template_length) {
              paired_matches_add(paired_matches,match_trace_end1,match_trace_end2,pair_orientation_concordant,template_length);
              ++num_concordant_pair_matches;
            } else if (discordant_search != pair_discordant_search_never) {
              paired_matches_add(paired_matches,match_trace_end1,match_trace_end2,pair_orientation_discordant,template_length);
              ++num_discordant_pair_matches;
            }
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
  // Update MCS
  PROF_STOP(GP_PAIRED_MATCHES_FIND_PAIRS);
}
GEM_INLINE void paired_matches_find_discordant_pairs(
    paired_matches_t* const paired_matches,search_parameters_t* const search_parameters) {
  const pair_discordant_search_t discordant_search = search_parameters->pair_discordant_search;
  switch (discordant_search) {
    case pair_discordant_search_never: break;
    case pair_discordant_search_only_if_no_concordant:
      if (vector_get_used(paired_matches->matches) > 0) break;
    // No break
    case pair_discordant_search_always:
      // Add discordant matches
      if (vector_get_used(paired_matches->discordant_matches) > 0) {
        // Merge discordant paired-matches
        const uint64_t num_discordant_matches = vector_get_used(paired_matches->discordant_matches);
        vector_reserve_additional(paired_matches->matches,num_discordant_matches);
        paired_match_t* const concordant_match = vector_get_free_elm(paired_matches->matches,paired_match_t);
        VECTOR_ITERATE_CONST(paired_matches->discordant_matches,discordant_match,dn,paired_match_t) {
          concordant_match[dn] = *discordant_match; // Add the discordant match
          paired_matches_counters_add(paired_matches->counters,discordant_match->distance); // Account the distance
        }
        vector_add_used(paired_matches->matches,num_discordant_matches);
        vector_clear(paired_matches->discordant_matches);
      }
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
/*
 * Utils
 */
int paired_matches_cmp_distance(const paired_match_t* const a,const paired_match_t* const b) {
  return a->distance - b->distance;
}
GEM_INLINE void paired_matches_sort_by_distance(paired_matches_t* const paired_matches) {
  // Sort global matches (match_trace_t) wrt distance
  qsort(vector_get_mem(paired_matches->matches,paired_match_t),
      vector_get_used(paired_matches->matches),sizeof(paired_match_t),
      (int (*)(const void *,const void *))paired_matches_cmp_distance);
}
int paired_matches_cmp_mapq_score(const paired_match_t* const a,const paired_match_t* const b) {
  return b->mapq_score - a->mapq_score;
}
GEM_INLINE void paired_matches_sort_by_mapq_score(paired_matches_t* const paired_matches) {
  // Sort global matches (match_trace_t) wrt distance
  qsort(vector_get_mem(paired_matches->matches,paired_match_t),
      vector_get_used(paired_matches->matches),sizeof(paired_match_t),
      (int (*)(const void *,const void *))paired_matches_cmp_mapq_score);
}
GEM_INLINE uint64_t paired_match_get_template_observed_length(
    const uint64_t begin_position_1,const uint64_t end_position_1,
    const uint64_t begin_position_2,const uint64_t end_position_2) {
  const uint64_t most_right = MAX(end_position_1,end_position_2);
  const uint64_t most_left = MIN(begin_position_1,begin_position_2);
  return most_right - most_left;
}
GEM_INLINE uint64_t paired_match_calculate_distance(
    const match_trace_t* const match_trace_end1,const match_trace_t* const match_trace_end2) {
  return match_trace_end1->distance+match_trace_end2->distance;
}

