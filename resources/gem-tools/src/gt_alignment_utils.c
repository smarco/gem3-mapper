/*
 * PROJECT: GEM-Tools library
 * FILE: gt_alignment_utils.c
 * DATE: 19/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_alignment_utils.h"

/*
 * Alignment basic tools
 */
GT_INLINE uint64_t gt_alignment_get_read_proportion(gt_alignment* const alignment,const float proportion) {
  GT_ALIGNMENT_CHECK(alignment);
  return gt_get_integer_proportion(proportion,gt_string_get_length(alignment->read));
}
/*
 * Alignment high-level insert (Update global state: counters, ...)
 */
GT_INLINE gt_map* gt_alignment_put_map(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,
    gt_map* const map,const bool check_duplicated,const bool replace_duplicated) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ALIGNMENT_CHECK(alignment); GT_MAP_CHECK(map);
  /*
   * Once the @map is given to @gt_alignment_put_map, it belong to the given @alignment.
   * So, in case of duplication, one of the maps will be deleted wrt @replace_duplicated
   */
  // Handle duplicates
  gt_map* found_map;
  uint64_t found_map_pos;
  if (gt_expect_false(check_duplicated && gt_alignment_find_map_fx(gt_map_cmp_fx,alignment,map,&found_map_pos,&found_map))) {
    if (gt_expect_true(replace_duplicated && gt_map_get_global_distance(map) < gt_map_get_global_distance(found_map))) {
      /* (gt_map_get_global_bases_aligned(map) <= gt_map_get_global_bases_aligned(found_map) ||
       *  gt_map_get_global_distance(map) < gt_map_get_global_distance(found_map)) */
      // Remove old map
      gt_alignment_dec_counter(alignment,gt_map_get_global_distance(found_map));
      gt_map_delete(found_map);
      // Replace old map
      gt_alignment_inc_counter(alignment,gt_map_get_global_distance(map));
      gt_alignment_set_map(alignment,map,found_map_pos);
      return map;
    } else {
      gt_map_delete(map);
      return found_map;
    }
  } else {
    // Add new map
    gt_alignment_inc_counter(alignment,gt_map_get_global_distance(map));
    gt_alignment_add_map(alignment,map);
    return map;
  }
}
GT_INLINE void gt_alignment_insert_map(gt_alignment* const alignment,gt_map* const map) {
  GT_ALIGNMENT_CHECK(alignment); GT_MAP_CHECK(map);
  gt_alignment_insert_map_fx(gt_map_cmp,alignment,map);
}
GT_INLINE void gt_alignment_insert_map_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,gt_map* const map) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ALIGNMENT_CHECK(alignment); GT_MAP_CHECK(map);
  // Prevent duplicates
  gt_alignment_put_map(gt_map_cmp_fx,alignment,map,true,true); /* TODO: Why replace_duplicated?? */
}
GT_INLINE void gt_alignment_insert_map_gt_vector(gt_alignment* const alignment,gt_vector* const map_vector) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_VECTOR_CHECK(map_vector);
  GT_VECTOR_ITERATE(map_vector,map_it,map_pos,gt_map*) {
    gt_alignment_insert_map(alignment,*map_it);
  }
}
GT_INLINE void gt_alignment_insert_map_fx_gt_vector(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,gt_vector* const map_vector) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ALIGNMENT_CHECK(alignment);
  GT_VECTOR_CHECK(map_vector);
  GT_VECTOR_ITERATE(map_vector,map_it,map_pos,gt_map*) {
    gt_alignment_insert_map_fx(gt_map_cmp_fx,alignment,*map_it);
  }
}
GT_INLINE bool gt_alignment_find_map_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,gt_map* const map,
    uint64_t* const found_map_pos,gt_map** const found_map) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ALIGNMENT_CHECK(alignment); GT_MAP_CHECK(map);
  GT_NULL_CHECK(found_map_pos); GT_NULL_CHECK(found_map);
  // Search for the map
  uint64_t pos = 0;
  GT_ALIGNMENT_ITERATE(alignment,map_it) {
    if (gt_map_cmp_fx(map_it,map)==0) {
      *found_map_pos = pos;
      *found_map = map_it;
      return true;
    }
    ++pos;
  }
  return false;
}
GT_INLINE bool gt_alignment_is_map_contained(gt_alignment* const alignment,gt_map* const map) {
  GT_ALIGNMENT_CHECK(alignment); GT_MAP_CHECK(map);
  gt_map* found_map;
  uint64_t found_map_pos;
  return gt_alignment_find_map_fx(gt_map_cmp,alignment,map,&found_map_pos,&found_map);
}
GT_INLINE bool gt_alignment_is_map_contained_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,gt_map* const map) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ALIGNMENT_CHECK(alignment); GT_MAP_CHECK(map);
  gt_map* found_map;
  uint64_t found_map_pos;
  return gt_alignment_find_map_fx(gt_map_cmp_fx,alignment,map,&found_map_pos,&found_map);
}
GT_INLINE void gt_alignment_reduce_maps(gt_alignment* const alignment,const uint64_t max_num_matches) {
  const uint64_t num_matches = gt_alignment_get_num_maps(alignment);
  if (max_num_matches < num_matches) {
    // Free unused maps
    GT_VECTOR_ITERATE_OFFSET(alignment->maps,map,map_pos,max_num_matches,gt_map*) {
      gt_map_delete(*map);
    }
    // Shrink maps vector
    gt_vector_set_used(alignment->maps,max_num_matches);
  }
}
/*
 * Alignment' Counters operators
 */
GT_INLINE bool gt_alignment_is_mapped(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  const bool unique_flag = gt_alignment_get_not_unique_flag(alignment);
  return unique_flag || gt_alignment_is_thresholded_mapped(alignment,UINT64_MAX);
}
GT_INLINE bool gt_alignment_is_thresholded_mapped(gt_alignment* const alignment,const int64_t max_allowed_strata) {
  GT_ALIGNMENT_CHECK(alignment);
  if (gt_alignment_get_not_unique_flag(alignment)) return true;
  gt_vector* vector = gt_alignment_get_counters_vector(alignment);
  if(gt_alignment_get_num_counters(alignment) == 0) return false;
  GT_VECTOR_ITERATE(vector,counter,counter_pos,uint64_t) {
    if (counter_pos>max_allowed_strata) return false;
    else if (*counter!=0) return true;
  }
  return false;
}
GT_INLINE void gt_alignment_recalculate_counters(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_vector_clear(gt_alignment_get_counters_vector(alignment));
  // Recalculate counters
  gt_alignment_map_iterator map_iterator;
  gt_alignment_new_map_iterator(alignment,&map_iterator);
  gt_map* map;
  while ((map=gt_alignment_next_map(&map_iterator))!=NULL) {
    gt_alignment_inc_counter(alignment,gt_map_get_global_distance(map));
  }
}
GT_INLINE void gt_alignment_recalculate_counters_no_splits(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  gt_vector_clear(gt_alignment_get_counters_vector(alignment));
  // Recalculate counters
  gt_alignment_map_iterator map_iterator;
  gt_alignment_new_map_iterator(alignment,&map_iterator);
  gt_map* map;
  while ((map=gt_alignment_next_map(&map_iterator))!=NULL) {
    gt_alignment_inc_counter(alignment,gt_map_get_no_split_distance(map));
  }
}
GT_INLINE int64_t gt_alignment_get_uniq_degree(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return gt_counters_get_uniq_degree(gt_alignment_get_counters_vector(alignment));
}
GT_INLINE int64_t gt_alignment_get_min_matching_strata(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  return gt_counters_get_min_matching_strata(gt_alignment_get_counters_vector(alignment));
}
GT_INLINE bool gt_alignment_get_next_matching_strata(
    gt_alignment* const alignment,const uint64_t begin_strata,
    uint64_t* const next_matching_strata,uint64_t* const num_maps) {
  GT_ALIGNMENT_CHECK(alignment);
  return gt_counters_get_next_matching_strata(
      gt_alignment_get_counters_vector(alignment),begin_strata,next_matching_strata,num_maps);
}
/*
 * Sort maps by score _> (int (*)(const void *,const void *))
 * i.e. (Sorting from smaller to bigger, <) := a-b
 *   cmp(a,b) := -n if (a<b)
 *                n if (a>b)
 *                0 if (a==b)
 */
int gt_alignment_cmp_distance__score(gt_map** const map_a,gt_map** const map_b) {
  // Sort by distance
  const int64_t distance_a = gt_map_get_global_distance(*map_a);
  const int64_t distance_b = gt_map_get_global_distance(*map_b);
  if (distance_a != distance_b) return distance_a-distance_b;
  // Sort by score (here we cannot do the trick as gt_score fills the whole uint64_t range)
  const uint64_t score_a = (*map_a)->gt_score;
  const uint64_t score_b = (*map_b)->gt_score;
  return (score_a > score_b) ? -1 : (score_a < score_b ? 1 : 0);
}
int gt_alignment_cmp_distance__score_no_split(gt_map** const map_a,gt_map** const map_b) {
  // Sort by distance
  const int64_t distance_a = gt_map_get_no_split_distance(*map_a);
  const int64_t distance_b = gt_map_get_no_split_distance(*map_b);
  if (distance_a != distance_b) return distance_a-distance_b;
  // Sort by score (here we cannot do the trick as gt_score fills the whole uint64_t range)
  const uint64_t score_a = (*map_a)->gt_score;
  const uint64_t score_b = (*map_b)->gt_score;
  return (score_a > score_b) ? -1 : (score_a < score_b ? 1 : 0);
}
GT_INLINE void gt_alignment_sort_by_distance__score(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  qsort(gt_vector_get_mem(alignment->maps,gt_map*),gt_vector_get_used(alignment->maps),
      sizeof(gt_map*),(int (*)(const void *,const void *))gt_alignment_cmp_distance__score);
}
GT_INLINE void gt_alignment_sort_by_distance__score_no_split(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  qsort(gt_vector_get_mem(alignment->maps,gt_map*),gt_vector_get_used(alignment->maps),
      sizeof(gt_map*),(int (*)(const void *,const void *))gt_alignment_cmp_distance__score_no_split);
}
/*
 * Alignment's Maps Utils
 */
GT_INLINE uint64_t gt_alignment_sum_mismatch_qualities(gt_alignment* const alignment,gt_map* const map) {
  GT_ALIGNMENT_CHECK(alignment);
  const char* const qualities = gt_alignment_get_qualities(alignment);
  uint64_t qv = 0;
  uint64_t split_map_offset = 0;
  GT_MAP_ITERATE_MAP_BLOCK(map, block){
    GT_MISMS_ITERATE(block, misms) {
      if (misms->misms_type == MISMS) {
        qv =+ qualities[misms->position + split_map_offset];
      }
    }
    split_map_offset += gt_map_get_base_length(block);
  }
  return qv;
}
GT_INLINE uint64_t gt_alignment_get_max_mismatch_quality(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  uint64_t max_qual = 0;
  GT_ALIGNMENT_ITERATE(alignment,map) {
    const uint64_t q = gt_alignment_sum_mismatch_qualities(alignment,map);
    if (q > max_qual) max_qual = q;
  }
  return max_qual;
}
/*
 * Alignment' Maps set-operators
 */
GT_INLINE void gt_alignment_merge_alignment_maps(gt_alignment* const alignment_dst,gt_alignment* const alignment_src) {
  GT_ALIGNMENT_CHECK(alignment_dst);
  GT_ALIGNMENT_CHECK(alignment_src);
  // Perform regular merge
  if (alignment_dst->alg_dictionary == NULL) {
    gt_alignment_merge_alignment_maps_fx(gt_map_cmp,alignment_dst,alignment_src);
    return;
  }
  // Perform fast iHash merge
  GT_ALIGNMENT_ITERATE(alignment_src,map_src) {
    gt_alignment_dictionary_element* alg_dicc_elem = NULL;
    gt_ihash_element* ihash_element_b = NULL;
    gt_ihash_element* ihash_element_e = NULL;
    const uint64_t vector_position = gt_vector_get_used(alignment_dst->maps);
    const uint64_t begin_position = gt_map_get_position(map_src)-gt_map_get_left_trim_length(map_src);
    const uint64_t end_position = gt_map_get_position(map_src)+gt_map_get_length(map_src);
    // Try add
    if (gt_alignment_dictionary_try_add(alignment_dst->alg_dictionary,map_src,
          begin_position,end_position,&alg_dicc_elem,&ihash_element_b,&ihash_element_e,vector_position)) {
      // (1) Ok, is new. We add it !
      gt_map* const map_src_cp = gt_map_copy(map_src);
      gt_alignment_inc_counter(alignment_dst,gt_map_get_global_distance(map_src_cp));
      gt_alignment_add_map(alignment_dst,map_src_cp);
    } else {
      // (2) One occurrence (could be a duplicate). Solve conflict
      uint64_t found_vector_position = 0;
      gt_map* map_found = NULL;
      bool found_candidate = false;
      if (ihash_element_b!=NULL) { // Check begin IDX ihash
        found_vector_position = *((uint64_t*)ihash_element_b->element);
        map_found = gt_alignment_get_map(alignment_dst,found_vector_position);
        if (gt_map_cmp(map_src,map_found)==0 && gt_map_less_than(map_src,map_found)) {
          found_candidate = true;
        }
      }
      if (!found_candidate && ihash_element_e!=NULL) { // Check end IDX ihash
        found_vector_position = *((uint64_t*)ihash_element_e->element);
        map_found = gt_alignment_get_map(alignment_dst,found_vector_position);
        if (gt_map_cmp(map_src,map_found)==0 && gt_map_less_than(map_src,map_found)) {
        	found_candidate = true;
        }
      }
      if (found_candidate) { // Is the same map !!
        gt_map* const map_src_cp = gt_map_copy(map_src);
        // Remove old map
        gt_alignment_dec_counter(alignment_dst,gt_map_get_global_distance(map_found));
        gt_map_delete(map_found);
        // Replace old map
        gt_alignment_inc_counter(alignment_dst,gt_map_get_global_distance(map_src_cp));
        gt_alignment_set_map(alignment_dst,map_src_cp,found_vector_position);
        // Record position at IDX iHash
        gt_alignment_dictionary_record_position(alignment_dst->alg_dictionary,begin_position,end_position,
            alg_dicc_elem,ihash_element_b,ihash_element_e,found_vector_position);
      } else {
        // (3) iHash won't solve the conflict. Resort to standard search
        if (gt_expect_false(gt_alignment_find_map_fx(gt_map_cmp,alignment_dst,map_src,&found_vector_position,&map_found))) {
          if (gt_expect_true(gt_map_less_than(map_src,map_found))) {
            gt_map* const map_src_cp = gt_map_copy(map_src);
            // Remove old map
            gt_alignment_dec_counter(alignment_dst,gt_map_get_global_distance(map_found));
            gt_map_delete(map_found);
            // Replace old map
            gt_alignment_inc_counter(alignment_dst,gt_map_get_global_distance(map_src_cp));
            gt_alignment_set_map(alignment_dst,map_src_cp,found_vector_position);
            // Record position at IDX iHash
            gt_alignment_dictionary_record_position(alignment_dst->alg_dictionary,begin_position,end_position,
                alg_dicc_elem,ihash_element_b,ihash_element_e,found_vector_position);
          }
        } else {
          // Add new map
          gt_map* const map_src_cp = gt_map_copy(map_src);
          gt_alignment_inc_counter(alignment_dst,gt_map_get_global_distance(map_src_cp));
          gt_alignment_add_map(alignment_dst,map_src_cp);
          // Record position at IDX iHash
          gt_alignment_dictionary_record_position(alignment_dst->alg_dictionary,begin_position,end_position,
              alg_dicc_elem,ihash_element_b,ihash_element_e,gt_vector_get_used(alignment_dst->maps)-1);
        }
      }
    }
  }
  gt_alignment_set_mcs(alignment_dst,GT_MIN(gt_alignment_get_mcs(alignment_dst),gt_alignment_get_mcs(alignment_src)));
}
GT_INLINE void gt_alignment_merge_alignment_maps_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_src) {
  GT_ALIGNMENT_CHECK(alignment_dst);
  GT_ALIGNMENT_CHECK(alignment_src);
  GT_ALIGNMENT_ITERATE(alignment_src,map_src) {
    gt_map* const map_src_cp = gt_map_copy(map_src);
    gt_alignment_put_map(gt_map_cmp_fx,alignment_dst,map_src_cp,true,true);
  }
  gt_alignment_set_mcs(alignment_dst,GT_MIN(gt_alignment_get_mcs(alignment_dst),gt_alignment_get_mcs(alignment_src)));
}
GT_INLINE gt_alignment* gt_alignment_union_alignment_maps_v(
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,va_list v_args) {
  GT_ZERO_CHECK(num_src_alignments);
  // Create new alignment
  gt_alignment* const alignment_union = gt_alignment_copy(alignment_src,false);
  alignment_union->alg_dictionary = gt_alignment_dictionary_new(alignment_union);
  gt_alignment_merge_alignment_maps(alignment_union,alignment_src);
  // Merge alignment sources into alignment_union
  uint64_t num_alg_merged = 1;
  while (num_alg_merged < num_src_alignments) {
    gt_alignment* alignment_target = va_arg(v_args,gt_alignment*);
    GT_ALIGNMENT_CHECK(alignment_target);
    gt_alignment_merge_alignment_maps(alignment_union,alignment_target);
    ++num_alg_merged;
  }
  // Clear
  gt_alignment_dictionary_delete(alignment_union->alg_dictionary);
  alignment_union->alg_dictionary = NULL;
  return alignment_union;
}
GT_INLINE gt_alignment* gt_alignment_union_alignment_maps_va(
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...) {
  GT_ZERO_CHECK(num_src_alignments);
  GT_ALIGNMENT_CHECK(alignment_src);
  va_list v_args;
  va_start(v_args,alignment_src);
  gt_alignment* const alignment_dst =
      gt_alignment_union_alignment_maps_v(num_src_alignments,alignment_src,v_args);
  va_end(v_args);
  return alignment_dst;
}
GT_INLINE gt_alignment* gt_alignment_union_alignment_maps_fx_v(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,va_list v_args) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ZERO_CHECK(num_src_alignments);
  // Create new alignment
  gt_alignment* const alignment_union = gt_alignment_copy(alignment_src,false);
  gt_alignment_merge_alignment_maps_fx(gt_map_cmp_fx,alignment_union,alignment_src);
  // Merge alignment sources into alignment_union
  uint64_t num_alg_merged = 1;
  while (num_alg_merged < num_src_alignments) {
    gt_alignment* alignment_target = va_arg(v_args,gt_alignment*);
    GT_ALIGNMENT_CHECK(alignment_target);
    gt_alignment_merge_alignment_maps_fx(gt_map_cmp_fx,alignment_union,alignment_target);
    ++num_alg_merged;
  }
  return alignment_union;
}
GT_INLINE gt_alignment* gt_alignment_union_alignment_maps_fx_va(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_ZERO_CHECK(num_src_alignments);
  GT_ALIGNMENT_CHECK(alignment_src);
  va_list v_args;
  va_start(v_args,alignment_src);
  gt_alignment* const alignment_dst =
      gt_alignment_union_alignment_maps_fx_v(gt_map_cmp_fx,num_src_alignments,alignment_src,v_args);
  va_end(v_args);
  return alignment_dst;
}
GT_INLINE gt_alignment* gt_alignment_subtract_alignment_maps_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_minuend,gt_alignment* const alignment_subtrahend) {
  GT_ALIGNMENT_CHECK(alignment_minuend);
  GT_ALIGNMENT_CHECK(alignment_subtrahend);
  // Create new alignment
  gt_alignment* const alignment_difference = gt_alignment_copy(alignment_minuend,false);
  // Copy not common maps
  GT_ALIGNMENT_ITERATE(alignment_minuend,map_minuend) {
    if (!gt_alignment_is_map_contained_fx(gt_map_cmp_fx,alignment_subtrahend,map_minuend)) { // TODO Improvement: Scheduled for v2.0
      gt_map* const map_cp = gt_map_copy(map_minuend);
      gt_alignment_put_map(gt_map_cmp_fx,alignment_difference,map_cp,true,false);
    }
  }
  return alignment_difference;
}
GT_INLINE gt_alignment* gt_alignment_subtract_alignment_maps(
    gt_alignment* const alignment_minuend,gt_alignment* const alignment_subtrahend) {
  GT_ALIGNMENT_CHECK(alignment_minuend);
  GT_ALIGNMENT_CHECK(alignment_subtrahend);
  return gt_alignment_subtract_alignment_maps_fx(gt_map_cmp,alignment_minuend,alignment_subtrahend);
}
GT_INLINE gt_alignment* gt_alignment_intersect_alignment_maps_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_src_A,gt_alignment* const alignment_src_B) {
  GT_ALIGNMENT_CHECK(alignment_src_A);
  GT_ALIGNMENT_CHECK(alignment_src_B);
  // Create new alignment
  gt_alignment* const alignment_intersection = gt_alignment_copy(alignment_src_A,false);
  // Copy common maps
  GT_ALIGNMENT_ITERATE(alignment_src_A,map_A) {
    if (gt_alignment_is_map_contained_fx(gt_map_cmp_fx,alignment_src_B,map_A)) { // TODO Improvement: Scheduled for v2.0
      gt_map* const map_cp = gt_map_copy(map_A);
      gt_alignment_put_map(gt_map_cmp_fx,alignment_intersection,map_cp,true,false);
    }
  }
  return alignment_intersection;
}
GT_INLINE gt_alignment* gt_alignment_intersect_alignment_maps(
    gt_alignment* const alignment_src_A,gt_alignment* const alignment_src_B) {
  GT_ALIGNMENT_CHECK(alignment_src_A);
  GT_ALIGNMENT_CHECK(alignment_src_B);
  return gt_alignment_intersect_alignment_maps_fx(gt_map_cmp,alignment_src_A,alignment_src_B);
}
/*
 * Alignment realignment
 */
GT_INLINE void gt_alignment_recover_mismatches(gt_alignment* const alignment,gt_sequence_archive* const sequence_archive) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  GT_ALIGNMENT_ITERATE(alignment,map) {
    gt_map_recover_mismatches_sa(map,alignment->read,sequence_archive);
  }
  gt_alignment_recalculate_counters(alignment);
}
GT_INLINE void gt_alignment_realign_hamming(gt_alignment* const alignment,gt_sequence_archive* const sequence_archive) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  GT_ALIGNMENT_ITERATE(alignment,map) {
    gt_map_realign_hamming_sa(map,alignment->read,sequence_archive);
  }
  gt_alignment_recalculate_counters(alignment);
}
GT_INLINE void gt_alignment_realign_levenshtein(gt_alignment* const alignment,gt_sequence_archive* const sequence_archive) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  GT_ALIGNMENT_ITERATE(alignment,map) {
    gt_map_realign_levenshtein_sa(map,alignment->read,sequence_archive);
  }
  gt_alignment_recalculate_counters(alignment);
}
GT_INLINE void gt_alignment_realign_weighted(
    gt_alignment* const alignment,gt_sequence_archive* const sequence_archive,int32_t (*gt_weigh_fx)(char*,char*)) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  GT_NULL_CHECK(gt_weigh_fx);
  GT_ALIGNMENT_ITERATE(alignment,map) {
    gt_map_realign_weighted_sa(map,alignment->read,sequence_archive,gt_weigh_fx);
  }
  gt_alignment_recalculate_counters(alignment);
}
/*
 * Alignment trimming
 */
GT_INLINE void gt_alignment_hard_trim(gt_alignment* const alignment,const uint64_t left,const uint64_t right) {
  GT_ALIGNMENT_CHECK(alignment);
  uint64_t read_length = gt_string_get_length(alignment->read);
  uint64_t qualities_length = gt_string_get_length(alignment->qualities);
  if (left+right >= read_length) return;
  /*
   * Apply left trim
   */
  if (left > 0) {
    gt_read_trim left_trim;
    left_trim.length = left;
    // Store trimmed read & apply trim
    left_trim.trimmed_read = gt_string_new(left+1);
    gt_string_set_nstring(left_trim.trimmed_read,gt_string_get_string(alignment->read),left);
    gt_string_trim_left(alignment->read,left);
    read_length-=left;
    if (qualities_length > 0) {
      // Store trimmed qualities & apply trim
      left_trim.trimmed_qualities = gt_string_new(left+1);
      gt_string_set_nstring(left_trim.trimmed_qualities,gt_string_get_string(alignment->qualities),left);
      gt_string_trim_left(alignment->qualities,left);
      qualities_length-=left;
    }
    // Annotate LEFT-trim
    gt_attributes_annotate_left_trim(alignment->attributes,&left_trim);
    // Apply LEFT-trim to all maps
    GT_ALIGNMENT_ITERATE(alignment,map) {
      gt_map_left_trim(map,left);
    }
  }
  /*
   * Apply right trim
   */
  if (right > 0) {
    gt_read_trim right_trim;
    right_trim.length = right;
    // Store trimmed read & apply trim
    right_trim.trimmed_read = gt_string_new(right+1);
    gt_string_set_nstring(right_trim.trimmed_read,gt_string_get_string(alignment->read),right);
    gt_string_trim_right(alignment->read,right);
    if (qualities_length > 0) {
      // Store trimmed qualities & apply trim
      right_trim.trimmed_qualities = gt_string_new(right+1);
      gt_string_set_nstring(right_trim.trimmed_qualities,gt_string_get_string(alignment->qualities),right);
      gt_string_trim_right(alignment->qualities,right);
    }
    // Annotate RIGHT-trim
    gt_attributes_annotate_right_trim(alignment->attributes,&right_trim);
    // Apply RIGHT-trim to all maps
    GT_ALIGNMENT_ITERATE(alignment,map) {
      gt_map_right_trim(map,right);
    }
  }
  // Recalculate counters
  gt_alignment_recalculate_counters(alignment);
}
GT_INLINE void gt_alignment_hard_trim_min_length(gt_alignment* const alignment,const uint64_t left,const uint64_t right,const uint64_t min_length) {
  GT_ALIGNMENT_CHECK(alignment);
  if (min_length+(left+right) > gt_string_get_length(alignment->read)) return;
  gt_alignment_hard_trim(alignment,left,right);
}
GT_INLINE void gt_alignment_quality_trim(gt_alignment* const alignment,const uint64_t quality_threshold,const uint64_t min_length) {
  GT_ALIGNMENT_CHECK(alignment);
  // TODO
}
GT_INLINE void gt_alignment_restore_trim(gt_alignment* const alignment) {
  GT_ALIGNMENT_CHECK(alignment);
  /*
   * Restore RIGHT-trim (if any)
   */
  gt_read_trim* const right_trim = gt_attributes_get_right_trim(alignment->attributes);
  if (right_trim!=NULL) {
    // Restore Read & Qualities
    if (right_trim->trimmed_read != NULL) gt_string_right_append_gt_string(alignment->read,right_trim->trimmed_read);
    if (right_trim->trimmed_qualities != NULL) gt_string_right_append_gt_string(alignment->qualities,right_trim->trimmed_qualities);
    // Restore Maps
    GT_ALIGNMENT_ITERATE(alignment,map) {
      gt_map_restore_right_trim(map,right_trim->length);
    }
    // Delete annotated trim
    gt_attributes_remove(alignment->attributes,GT_ATTR_ID_RIGHT_TRIM);
  }
  /*
   * Restore LEFT-trim (if any)
   */
  gt_read_trim* const left_trim = gt_attributes_get_left_trim(alignment->attributes);
  if (left_trim!=NULL) {
    // Restore Read & Qualities
    if (left_trim->trimmed_read != NULL) gt_string_left_append_gt_string(alignment->read,left_trim->trimmed_read);
    if (left_trim->trimmed_qualities != NULL) gt_string_left_append_gt_string(alignment->qualities,left_trim->trimmed_qualities);
    // Restore Maps
    GT_ALIGNMENT_ITERATE(alignment,map) {
      gt_map_restore_left_trim(map,left_trim->length);
    }
    // Delete annotated trim
    gt_attributes_remove(alignment->attributes,GT_ATTR_ID_LEFT_TRIM);
  }
  // Recalculate counters
  gt_alignment_recalculate_counters(alignment);
}
