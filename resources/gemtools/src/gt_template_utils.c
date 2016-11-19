/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template_utils.h
 * DATE: 19/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_template_utils.h"

/*
 * Template Pair Info Setup
 *   PE read => template.GT_ATTR_ID_TAG_PAIR = 0, template.end1.GT_ATTR_ID_TAG_PAIR = 1, template.end2.GT_ATTR_ID_TAG_PAIR = 2
 *   SE read => template.GT_ATTR_ID_TAG_PAIR = X, template.end1.GT_ATTR_ID_TAG_PAIR = X
 *              alignment.GT_ATTR_ID_TAG_PAIR = X
 *   All of them have the same CASAVA and EXTRA_TAG info
 */
GT_INLINE void gt_template_setup_pair_attributes_to_alignments(gt_template* const template,const bool copy_tags) {
  GT_TEMPLATE_CHECK(template);
  const uint64_t num_blocks = gt_template_get_num_blocks(template);
  int64_t i, p;
  // Setup all alignments' tags + attr
  for (i=0;i<num_blocks;i++) {
    gt_alignment* const alignment = gt_template_get_block(template,i);
    // Copy all attributes
    gt_attributes_copy(alignment->attributes,template->attributes);
    if (copy_tags) gt_string_copy(alignment->tag,template->tag);
    if (num_blocks>1) {
      p = i+1;
    } else {
      int64_t* const attr = gt_attributes_get(template->attributes,GT_ATTR_ID_TAG_PAIR);
      p = (attr!=NULL) ? *attr : 0;
    }
    gt_attributes_add(alignment->attributes,GT_ATTR_ID_TAG_PAIR,&p,int64_t);
  }
  // Clear template's pair info
  if (num_blocks > 1) {
    p = 0;
    gt_attributes_add(template->attributes,GT_ATTR_ID_TAG_PAIR,&p,int64_t);
  }
}
GT_INLINE uint64_t gt_template_get_read_proportion(gt_template* const template,const float proportion) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_read_proportion(alignment,proportion);
  } GT_TEMPLATE_END_REDUCTION;
  const uint64_t read_length_end1 = gt_string_get_length(gt_template_get_end1(template)->read);
  const uint64_t read_length_end2 = gt_string_get_length(gt_template_get_end2(template)->read);
  return gt_get_integer_proportion(proportion,read_length_end1+read_length_end2);
}

/*
 * Template's MMaps operators (Update global state: counters, ...)
 */
GT_INLINE void gt_template_alias_dup_mmap_members(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_template* const template,
    gt_map** const mmap,gt_map** const uniq_mmaps) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap); GT_NULL_CHECK(uniq_mmaps);
  // Resolve mmap
  const uint64_t num_blocks = gt_template_get_num_blocks(template);
  uint64_t i;
  for (i=0;i<num_blocks;++i) {
    GT_NULL_CHECK(mmap[i]);
    uniq_mmaps[i] = gt_alignment_put_map(gt_map_cmp_fx,gt_template_get_block(template,i),mmap[i],true,false);
  }
}
GT_INLINE gt_map** gt_template_raw_put_mmap(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_template* const template,
    gt_map** const mmap,gt_mmap_attributes* const mmap_attributes) {
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap);
  GT_NULL_CHECK(mmap_attributes);
  // Resolve mmap aliasing/insertion
  gt_map** uniq_mmaps = gt_calloc(gt_template_get_num_blocks(template),gt_map*,false);
  gt_template_alias_dup_mmap_members(gt_map_cmp_fx,template,mmap,uniq_mmaps);
  // Raw mmap insertion
  gt_template_add_mmap_array(template,uniq_mmaps,mmap_attributes);
  gt_map** const template_mmap = gt_template_get_mmap_array(template,gt_template_get_num_mmaps(template)-1,NULL);
  gt_free(uniq_mmaps); // Free auxiliary vector
  return template_mmap;
}
GT_INLINE gt_map** gt_template_put_mmap(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr,
    const bool check_duplicated,const bool replace_duplicated) {
  GT_NULL_CHECK(gt_mmap_cmp_fx); GT_NULL_CHECK(gt_map_cmp_fx);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap);
  GT_NULL_CHECK(mmap_attr);
  // Check mmap duplicates
  gt_map** found_mmap;
  gt_mmap_attributes* found_mmap_attributes=0;
  uint64_t found_mmap_pos=0;
  bool is_duplicated = false;
  if (check_duplicated) {
    is_duplicated = gt_expect_false(gt_template_find_mmap_fx(gt_mmap_cmp_fx,
          template,mmap,&found_mmap_pos,&found_mmap,&found_mmap_attributes));
  }
  gt_map** template_mmap;
  if (!check_duplicated || !is_duplicated || replace_duplicated) { // TODO: Chose which to replace (like alignment)
    // Resolve mmap aliasing/insertion
    gt_map** uniq_mmaps = gt_calloc(gt_template_get_num_blocks(template),gt_map*,false);
    gt_template_alias_dup_mmap_members(gt_map_cmp_fx,template,mmap,uniq_mmaps);
    // Insert mmap
    if (!check_duplicated || !is_duplicated) { // Add new mmap
      gt_template_inc_counter(template,mmap_attr->distance);
      gt_template_add_mmap_array(template,uniq_mmaps,mmap_attr);
      template_mmap=gt_template_get_mmap_array(template,gt_template_get_num_mmaps(template)-1,NULL);
    } else { // Replace mmap
      gt_template_dec_counter(template,found_mmap_attributes->distance); // Remove old mmap
      gt_template_set_mmap_array(template,found_mmap_pos,uniq_mmaps,mmap_attr); // Replace old mmap
      gt_template_inc_counter(template,mmap_attr->distance);
      template_mmap=gt_template_get_mmap_array(template,found_mmap_pos,NULL);
    }
    gt_free(uniq_mmaps); // Free auxiliary vector
  } else {
    // Delete mmap
    GT_MMAP_ITERATE_ENDS(mmap,gt_template_get_num_blocks(template),map,end_pos) {
      gt_map_delete(map);
    }
    template_mmap=gt_template_get_mmap_array(template,found_mmap_pos,NULL);
  }
  return template_mmap;
}
GT_INLINE void gt_template_insert_mmap(
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attributes) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap);
  GT_NULL_CHECK(mmap_attributes);
  gt_template_insert_mmap_fx(gt_mmap_cmp,template,mmap,mmap_attributes);
}
GT_INLINE void gt_template_insert_mmap_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attributes) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap);
  GT_NULL_CHECK(mmap_attributes);
  gt_template_put_mmap(gt_mmap_cmp_fx,gt_map_cmp,template,mmap,mmap_attributes,true,true); // TODO: Why replace, why?
}
GT_INLINE void gt_template_insert_mmap_gtvector(
    gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attributes) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap);
  GT_NULL_CHECK(mmap_attributes);
  gt_template_insert_mmap_gtvector_fx(gt_mmap_cmp,template,mmap,mmap_attributes);
}
GT_INLINE void gt_template_insert_mmap_gtvector_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attributes) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap);
  GT_NULL_CHECK(mmap_attributes);
  gt_check(gt_vector_get_used(mmap)!=gt_template_get_num_blocks(template),TEMPLATE_ADD_BAD_NUM_BLOCKS);
  gt_template_insert_mmap_fx(gt_mmap_cmp_fx,template,gt_vector_get_mem(mmap,gt_map*),mmap_attributes);
}
GT_INLINE bool gt_template_find_mmap_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap,
    uint64_t* const found_mmap_pos,gt_map*** const found_mmap,gt_mmap_attributes** const found_mmap_attributes) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap);
  GT_NULL_CHECK(found_mmap_pos);
  GT_NULL_CHECK(found_mmap);
  // Search for the mmap
  const uint64_t num_blocks = gt_template_get_num_blocks(template);
  uint64_t pos = 0;
  if (template->alg_dictionary == NULL || template->alg_dictionary->refs_dictionary == NULL) {
    GT_TEMPLATE_ITERATE_MMAP__ATTR_(template,template_mmap,mmap_attribute) {
      if (gt_mmap_cmp_fx(template_mmap,mmap,num_blocks)==0) {
        *found_mmap_pos = pos;
        *found_mmap = template_mmap;
        if (found_mmap_attributes) *found_mmap_attributes = mmap_attribute;
        return true;
      }
      ++pos;
    }
  } else {
    // indexed search only through other templates
    // with the first map on the same chromosome
    char* seq_name = gt_map_get_seq_name(mmap[0]);
    if (!gt_shash_is_contained(template->alg_dictionary->refs_dictionary, seq_name)) {
      return false;
    }
    gt_vector* dict_elements = gt_shash_get(template->alg_dictionary->refs_dictionary, seq_name, gt_vector);
    GT_VECTOR_ITERATE(dict_elements, e, c, gt_template_dictionary_map_element*) {
      gt_map** template_mmap = (*e)->mmap;
      if (gt_mmap_cmp_fx(template_mmap,mmap,num_blocks)==0) {
        *found_mmap_pos = pos;
        *found_mmap = template_mmap;
        if (found_mmap_attributes) *found_mmap_attributes = (*e)->mmap_attrs;
        return true;
      }
      ++pos;
    }

  }
  return false;
}
GT_INLINE bool gt_template_is_mmap_contained(gt_template* const template,gt_map** const mmap) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap);
  gt_map** found_mmap;
  uint64_t found_mmap_pos;
  return gt_template_find_mmap_fx(gt_mmap_cmp,template,mmap,
      &found_mmap_pos,&found_mmap,NULL);
}
GT_INLINE bool gt_template_is_mmap_contained_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap);
  gt_map** found_mmap;
  uint64_t found_mmap_pos;
  return gt_template_find_mmap_fx(gt_mmap_cmp_fx,template,mmap,
      &found_mmap_pos,&found_mmap,NULL);
}

GT_INLINE void gt_template_reduce_mmaps(gt_template* const template,const uint64_t max_num_matches) {
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_reduce_maps(alignment,max_num_matches);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  const uint64_t num_matches = gt_template_get_num_mmaps(template);
  if (max_num_matches < num_matches) {
    gt_vector_set_used(template->mmaps,max_num_matches);
    gt_template_recalculate_counters(template);
  }
}

/*
 * Insert size should be and estimate of the original fragment size that was sequenced.  We can get this from the following:
 *
 * position of rightmost read + no. bases in rightmost read - position of leftmost read
 * 
 * If the read has been trimmed from the start of the read then we can't get the original size, but this is a relatively
 * rare event.  Trimming from the end of the read does not effect the calculation as the rightmost read will be on the
 * negative strand so trimming x bases from the end of the read will shift the mapped position of the read by x.
 *
 * Split mappings require special handling in that we need to consider the number of bases read + the distance between the 
 * last block in each mapping as follows:
 *
 * Position of last block of rightmost read + no. bases in rightmost read - (position of last block of leftmost read + no. bases)
 * in all other blocks of leftmost read
 *
 * if start_x is non-zero then *start_x will be set to the start position of the left block before the insert
 *
 */
GT_INLINE int64_t gt_template_get_insert_size(gt_map** const mmap,gt_status* gt_error,uint64_t* start_x,gt_string** ctg) {
  // Get last block of each map
  gt_map *block[2]={0,0};
  uint64_t length[2]={0,0};
  int64_t x=0;
  *gt_error=GT_TEMPLATE_INSERT_SIZE_OK;
  GT_MAP_ITERATE(mmap[0],map_it) {
    block[0]=map_it;
    length[0]+=gt_map_get_base_length(block[0]);
  } 
  GT_MAP_ITERATE(mmap[1],map_it2) {
    block[1]=map_it2;
    length[1]+=gt_map_get_base_length(block[1]);
  } 
  if(gt_string_equals(block[0]->seq_name,block[1]->seq_name)) {
    if(block[0]->strand!=block[1]->strand) {
      if(block[0]->strand==FORWARD) {
      		x=1+block[1]->position+length[1]-(block[0]->position+length[0]-gt_map_get_base_length(block[0]));
      		if(start_x) *start_x=block[0]->position;
      } else {
      	x=1+block[0]->position+length[0]-(block[1]->position+length[1]-gt_map_get_base_length(block[1]));
    		if(start_x) *start_x=block[1]->position;
      }
      if(ctg) *ctg=block[0]->seq_name;
    } else {
      *gt_error=GT_TEMPLATE_INSERT_SIZE_SAME_STRAND;
    }
  } else {
    *gt_error=GT_TEMPLATE_INSERT_SIZE_DIFFERENT_CONTIGS;
  }
  return x;
}
/*
 * Template's Counters operators
 */
GT_INLINE bool gt_template_is_mapped(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  const bool unique_flag = gt_template_get_not_unique_flag(template);
  return unique_flag || gt_template_is_thresholded_mapped(template,UINT64_MAX);
}
GT_INLINE bool gt_template_is_thresholded_mapped(gt_template* const template,const uint64_t max_allowed_strata) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_is_thresholded_mapped(alignment,max_allowed_strata);
  } GT_TEMPLATE_END_REDUCTION;
  if (gt_template_get_not_unique_flag(template)) return true;
  if(gt_template_get_num_counters(template) == 0) return false;
  gt_vector* vector = gt_template_get_counters_vector(template);
  GT_VECTOR_ITERATE(vector,counter,counter_pos,uint64_t) {
    if (counter_pos>max_allowed_strata) return false;
    else if (*counter!=0) return true;
  }
  return false;
}

GT_INLINE void gt_template_recalculate_counters(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_recalculate_counters(alignment);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Clear previous counters
  gt_vector_clear(gt_template_get_counters_vector(template));
  // Recalculate counters
  const uint64_t num_blocks = gt_template_get_num_blocks(template);
  GT_TEMPLATE_ITERATE_MMAP__ATTR_(template,mmap,mmap_attr) {
    uint64_t i, total_distance = 0;
    for (i=0;i<num_blocks;++i) {
      total_distance+=gt_map_get_global_distance(mmap[i]);
    }
    mmap_attr->distance=total_distance;
    gt_template_inc_counter(template,total_distance);
  }
}

GT_INLINE void gt_template_recalculate_counters_no_splits(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_recalculate_counters_no_splits(alignment);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Clear previous counters
  gt_vector_clear(gt_template_get_counters_vector(template));
  // Recalculate counters
  const uint64_t num_blocks = gt_template_get_num_blocks(template);
  GT_TEMPLATE_ITERATE_MMAP__ATTR_(template,mmap,mmap_attr) {
    uint64_t i, total_distance = 0;
    for (i=0;i<num_blocks;++i) {
      total_distance+=gt_map_get_no_split_distance(mmap[i]);
    }
    mmap_attr->distance=total_distance;
    gt_template_inc_counter(template,total_distance);
  }
}


GT_INLINE int64_t gt_template_get_min_matching_strata(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_min_matching_strata(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return gt_counters_get_min_matching_strata(gt_template_get_counters_vector(template));
}
GT_INLINE int64_t gt_template_get_uniq_degree(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_uniq_degree(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return gt_counters_get_uniq_degree(gt_template_get_counters_vector(template));
}
GT_INLINE bool gt_template_get_next_matching_strata(
    gt_template* const template,const uint64_t begin_strata,
    uint64_t* const next_matching_strata,uint64_t* const num_maps) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_next_matching_strata(alignment,begin_strata,next_matching_strata,num_maps);
  } GT_TEMPLATE_END_REDUCTION;
  return gt_counters_get_next_matching_strata(gt_template_get_counters_vector(template),
      begin_strata,next_matching_strata,num_maps);
}

/*
 * Sort maps by score _> (int (*)(const void *,const void *))
 * i.e. (Sorting from smaller to bigger, <) := a-b
 *   cmp(a,b) := -n if (a<b)
 *                n if (a>b)
 *                0 if (a==b)
 */
typedef struct {
  gt_mmap_attributes attributes;
  uint64_t init_position;
  gt_map* end_1;
  gt_map* end_2;
} gt_mmap_placeholder;
int gt_mmap_cmp_distance__score(gt_mmap* const mmap_a,gt_mmap* const mmap_b) {
  // Sort by distance
  const int64_t distance_a = mmap_a->attributes.distance;
  const int64_t distance_b = mmap_b->attributes.distance;
  if (distance_a != distance_b) return distance_a-distance_b;
  // Sort by score (here we cannot do the trick as gt_score fills the whole uint64_t range)
  const uint64_t score_a = mmap_a->attributes.gt_score;
  const uint64_t score_b = mmap_b->attributes.gt_score;
  return (score_a > score_b) ? -1 : (score_a < score_b ? 1 : 0);
}
/*
 * Template's Maps Sorting
 */
GT_INLINE void gt_template_sort_by_distance__score(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_sort_by_distance__score(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  // Sort
  const uint64_t num_mmap = gt_template_get_num_mmaps(template);
  qsort(gt_vector_get_mem(template->mmaps,gt_mmap),num_mmap,sizeof(gt_mmap),
      (int (*)(const void *,const void *))gt_mmap_cmp_distance__score);
}

GT_INLINE void gt_template_sort_by_distance__score_no_split(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_sort_by_distance__score_no_split(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  // Sort
  const uint64_t num_mmap = gt_template_get_num_mmaps(template);
  qsort(gt_vector_get_mem(template->mmaps,gt_mmap),num_mmap,sizeof(gt_mmap),
      (int (*)(const void *,const void *))gt_mmap_cmp_distance__score);
}
/*
 * Template's MMaps Utils
 */
GT_INLINE uint64_t gt_template_sum_mismatch_qualities(gt_template* const template,gt_map** const mmap) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_sum_mismatch_qualities(alignment,mmap[0]);
  } GT_TEMPLATE_END_REDUCTION;
  return gt_alignment_sum_mismatch_qualities(gt_template_get_block(template,0),mmap[0]) +
         gt_alignment_sum_mismatch_qualities(gt_template_get_block(template,1),mmap[1]);
}
GT_INLINE uint64_t gt_template_get_max_mismatch_quality(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_max_mismatch_quality(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  uint64_t max_qual = 0;
  GT_TEMPLATE_ITERATE_MMAP_(template,mmap) {
    const uint64_t q = gt_template_sum_mismatch_qualities(template,mmap);
    if (q > max_qual) max_qual = q;
  }
  return max_qual;
}
/*
 * Template Set operators
 */
GT_INLINE void gt_template_merge_template_mmaps(gt_template* const template_dst,gt_template* const template_src) {
  GT_TEMPLATE_CHECK(template_dst);
  GT_TEMPLATE_CHECK(template_src);
  GT_TEMPLATE_COMMON_CONSISTENCY_ERROR(template_dst,template_src);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_src,alignment_src) {
    GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_dst,alignment_dst) {
      gt_alignment_merge_alignment_maps(alignment_dst,alignment_src);
    } GT_TEMPLATE_END_REDUCTION;
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Merge mmaps
  gt_template_merge_template_mmaps_fx(gt_mmap_cmp,gt_map_cmp,template_dst,template_src);
}
GT_INLINE void gt_template_merge_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template_dst,gt_template* const template_src) {
  GT_NULL_CHECK(gt_mmap_cmp_fx); GT_NULL_CHECK(gt_map_cmp_fx);
  GT_TEMPLATE_COMMON_CONSISTENCY_ERROR(template_dst,template_src);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_src,alignment_src) {
    GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_dst,alignment_dst) {
      gt_alignment_merge_alignment_maps_fx(gt_map_cmp_fx,alignment_dst,alignment_src);
    } GT_TEMPLATE_END_REDUCTION;
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Merge mmaps
  bool use_hash = gt_template_get_num_mmaps(template_src) > 100 || gt_template_get_num_mmaps(template_dst) > 100;
  if(use_hash){
    template_dst->alg_dictionary = gt_template_dictionary_new(template_dst);
    gt_template_dictionary_add_ref(template_dst->alg_dictionary, template_dst);
  }

  GT_TEMPLATE_CHECK(template_dst);
  GT_TEMPLATE_CHECK(template_src);
  GT_TEMPLATE_ITERATE_MMAP__ATTR(template_src,mmap,mmap_attr) {
    gt_map** mmap_copy = gt_mmap_array_copy(mmap,__mmap_num_blocks);
    gt_template_put_mmap(gt_mmap_cmp_fx,gt_map_cmp_fx,template_dst,mmap_copy,mmap_attr,true,true);
    gt_free(mmap_copy); // Free array handler
  }
  gt_template_set_mcs(template_dst,GT_MIN(gt_template_get_mcs(template_dst),gt_template_get_mcs(template_src)));
  if(use_hash){
    gt_template_dictionary_delete(template_dst->alg_dictionary);
    template_dst->alg_dictionary = NULL;
  }
}
GT_INLINE gt_template* gt_template_union_template_mmaps_v(
    const uint64_t num_src_templates,gt_template* const template_src,va_list v_args) {
  GT_ZERO_CHECK(num_src_templates);
  // Create new template
  gt_template* const template_union = gt_template_dup(template_src,false,false);
  gt_template_merge_template_mmaps(template_union,template_src);
  // Merge template sources into template_union
  uint64_t num_tmp_merged = 1;
  while (num_tmp_merged < num_src_templates) {
    gt_template* template_target = va_arg(v_args,gt_template*);
    GT_TEMPLATE_COMMON_CONSISTENCY_ERROR(template_union,template_target);
    GT_TEMPLATE_CHECK(template_target);
    gt_template_merge_template_mmaps(template_union,template_target);
    ++num_tmp_merged;
  }
  return template_union;
}
GT_INLINE gt_template* gt_template_union_template_mmaps_va(
    const uint64_t num_src_templates,gt_template* const template_src,...) {
  GT_ZERO_CHECK(num_src_templates);
  GT_TEMPLATE_CHECK(template_src);
  va_list v_args;
  va_start(v_args,template_src);
  gt_template* const template_union =
      gt_template_union_template_mmaps_v(num_src_templates,template_src,v_args);
  va_end(v_args);
  return template_union;
}
GT_INLINE gt_template* gt_template_union_template_mmaps_a(
    gt_template** const templates,const uint64_t num_src_templates) {
  GT_ZERO_CHECK(num_src_templates);
  // Create new template
  gt_template* const template_union = gt_template_dup(templates[0],false,false);
  gt_template_merge_template_mmaps(template_union,templates[0]);
  // Merge template sources into template_union
  uint64_t i;
  for (i=1;i<num_src_templates;++i) {
    GT_TEMPLATE_COMMON_CONSISTENCY_ERROR(template_union,templates[i]);
    GT_TEMPLATE_CHECK(templates[i]);
    gt_template_merge_template_mmaps(template_union,templates[i]);
  }
  return template_union;
}

GT_INLINE gt_template* gt_template_union_template_mmaps_fx_v(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    const uint64_t num_src_templates,gt_template* const template_src,va_list v_args) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_ZERO_CHECK(num_src_templates);
  // Create new template
  gt_template* const template_union = gt_template_dup(template_src,false,false);
  gt_template_merge_template_mmaps_fx(gt_mmap_cmp_fx,gt_map_cmp_fx,template_union,template_src);
  // Merge template sources into template_union
  uint64_t num_tmp_merged = 1;
  while (num_tmp_merged < num_src_templates) {
    gt_template* template_target = va_arg(v_args,gt_template*);
    GT_TEMPLATE_COMMON_CONSISTENCY_ERROR(template_union,template_target);
    GT_TEMPLATE_CHECK(template_target);
    gt_template_merge_template_mmaps_fx(gt_mmap_cmp_fx,gt_map_cmp_fx,template_union,template_target);
    ++num_tmp_merged;
  }
  return template_union;
}
GT_INLINE gt_template* gt_template_union_template_mmaps_fx_va(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    const uint64_t num_src_templates,gt_template* const template_src,...) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_ZERO_CHECK(num_src_templates);
  GT_TEMPLATE_CHECK(template_src);
  va_list v_args;
  va_start(v_args,template_src);
  gt_template* const template_union =
      gt_template_union_template_mmaps_fx_v(gt_mmap_cmp_fx,gt_map_cmp_fx,num_src_templates,template_src,v_args);
  va_end(v_args);
  return template_union;
}

GT_INLINE gt_template* gt_template_subtract_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template_minuend,gt_template* const template_subtrahend) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_TEMPLATE_COMMON_CONSISTENCY_ERROR(template_minuend,template_subtrahend);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_minuend,alignment_minuend) {
    GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_subtrahend,alignment_subtrahend) {
      gt_alignment* const alignment_difference =
          gt_alignment_subtract_alignment_maps_fx(gt_map_cmp_fx,alignment_minuend,alignment_subtrahend);
      gt_template* const template_difference = gt_template_new();
      gt_template_add_block(template_difference,alignment_difference);
      return template_difference;
    } GT_TEMPLATE_END_REDUCTION;
  } GT_TEMPLATE_END_REDUCTION;
  // Subtract
  GT_TEMPLATE_CHECK(template_minuend);
  GT_TEMPLATE_CHECK(template_subtrahend);
  gt_template* const template_difference = gt_template_dup(template_minuend,false,false);
  uint64_t found_mmap_pos;
  gt_map** found_mmap;
  gt_mmap_attributes* found_mmap_attr;
  GT_TEMPLATE_ITERATE_MMAP__ATTR(template_minuend,mmap,mmap_attr) {
    if (!gt_template_find_mmap_fx(gt_mmap_cmp_fx,template_subtrahend,mmap,&found_mmap_pos,&found_mmap,&found_mmap_attr)) {
      gt_map** mmap_copy = gt_mmap_array_copy(mmap,__mmap_num_blocks);
      gt_template_put_mmap(gt_mmap_cmp_fx,gt_map_cmp_fx,template_difference,mmap_copy,mmap_attr,true,false);
      gt_free(mmap_copy);
    }
  }
  return template_difference;
}
GT_INLINE gt_template* gt_template_subtract_template_mmaps(
    gt_template* const template_minuend,gt_template* const template_subtrahend) {
  GT_TEMPLATE_CHECK(template_minuend);
  GT_TEMPLATE_CHECK(template_subtrahend);
  return gt_template_subtract_template_mmaps_fx(gt_mmap_cmp,gt_map_cmp,template_minuend,template_subtrahend);
}

GT_INLINE gt_template* gt_template_intersect_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template_A,gt_template* const template_B) {
  GT_NULL_CHECK(gt_mmap_cmp_fx);
  GT_NULL_CHECK(gt_map_cmp_fx);
  GT_TEMPLATE_COMMON_CONSISTENCY_ERROR(template_A,template_B);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_A,alignment_A) {
    GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_B,alignment_B) {
      gt_alignment* const alignment_intersection =
          gt_alignment_intersect_alignment_maps_fx(gt_map_cmp_fx,alignment_A,alignment_B);
      gt_template* const template_intersection = gt_template_new();
      gt_template_add_block(template_intersection,alignment_intersection);
      return template_intersection;
    } GT_TEMPLATE_END_REDUCTION;
  } GT_TEMPLATE_END_REDUCTION;
  // Intersect
  GT_TEMPLATE_CHECK(template_A);
  GT_TEMPLATE_CHECK(template_B);
  gt_template* const template_intersection = gt_template_dup(template_A,false,false);
  uint64_t found_mmap_pos;
  gt_map** found_mmap;
  gt_mmap_attributes* found_mmap_attr;
  GT_TEMPLATE_ITERATE_MMAP__ATTR(template_A,mmap,mmap_attr) {
    if (gt_template_find_mmap_fx(gt_mmap_cmp_fx,template_B,mmap,&found_mmap_pos,&found_mmap,&found_mmap_attr)) {
      gt_map** mmap_copy = gt_mmap_array_copy(mmap,__mmap_num_blocks);
      gt_template_put_mmap(gt_mmap_cmp_fx,gt_map_cmp_fx,template_intersection,mmap_copy,mmap_attr,true,false);
      gt_free(mmap_copy);
    }
  }
  return template_intersection;
}
GT_INLINE gt_template* gt_template_intersect_template_mmaps(
    gt_template* const template_A,gt_template* const template_B) {
  GT_TEMPLATE_CHECK(template_A);
  GT_TEMPLATE_CHECK(template_B);
  return gt_template_intersect_template_mmaps_fx(gt_mmap_cmp,gt_map_cmp,template_A,template_B);
}

/*
 * Template realignment
 */
GT_INLINE void gt_template_recover_mismatches(gt_template* const template,gt_sequence_archive* const sequence_archive) {
  GT_TEMPLATE_CHECK(template);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
    gt_alignment_recover_mismatches(alignment,sequence_archive);
  }
  if (gt_template_get_num_blocks(template)>1) gt_template_recalculate_counters(template);
}
GT_INLINE void gt_template_realign_hamming(gt_template* const template,gt_sequence_archive* const sequence_archive) {
  GT_TEMPLATE_CHECK(template);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
    gt_alignment_realign_hamming(alignment,sequence_archive);
  }
  if (gt_template_get_num_blocks(template)>1) gt_template_recalculate_counters(template);
}
GT_INLINE void gt_template_realign_levenshtein(gt_template* const template,gt_sequence_archive* const sequence_archive) {
  GT_TEMPLATE_CHECK(template);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
    gt_alignment_realign_levenshtein(alignment,sequence_archive);
  }
  if (gt_template_get_num_blocks(template)>1) gt_template_recalculate_counters(template);
}
GT_INLINE void gt_template_realign_weighted(
    gt_template* const template,gt_sequence_archive* const sequence_archive,int32_t (*gt_weigh_fx)(char*,char*)) {
  GT_TEMPLATE_CHECK(template);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  GT_NULL_CHECK(gt_weigh_fx);
  GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
    gt_alignment_realign_weighted(alignment,sequence_archive,gt_weigh_fx);
  }
  if (gt_template_get_num_blocks(template)>1) gt_template_recalculate_counters(template);
}

/*
 * Template trimming
 */
GT_INLINE void gt_template_hard_trim(gt_template* const template,const uint64_t left,const uint64_t right) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
    gt_alignment_hard_trim(alignment,left,right);
  }
  if (gt_template_is_paired_end(template)) gt_template_recalculate_counters(template);
}
GT_INLINE void gt_template_restore_trim(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
    gt_alignment_restore_trim(alignment);
  }
  if (gt_template_is_paired_end(template)) gt_template_recalculate_counters(template);
}

///*
// * Template trimming
// */
//GT_INLINE void gt_template_trim(gt_template* const template,uint64_t const left,uint64_t const right,uint64_t const min_length,const bool set_extra) {
//  GT_TEMPLATE_CHECK(template);
//  GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
//    gt_alignment_hard_trim(alignment,left,right,min_length,set_extra);
//  }
//  if (set_extra) { // TODO: Why?
//    gt_alignment* first = gt_template_get_block(template,0);
//    // Transfer extra attribute
//    gt_string* alignment_extra = gt_attributes_get(first->attributes,GT_ATTR_ID_TAG_EXTRA);
//    if (gt_attributes_is_contained(template->attributes,GT_ATTR_ID_TAG_EXTRA)) {
//      gt_string* old = gt_attributes_get(template->attributes,GT_ATTR_ID_TAG_EXTRA);
//      gt_string_copy(old,alignment_extra);
//    } else {
//      gt_string* const extra = gt_string_new(8 + (2*left) + (2*right));
//      gt_string_copy(extra,alignment_extra);
//      gt_attributes_add_string(template->attributes,GT_ATTR_ID_TAG_EXTRA,extra);
//    }
//  }
//}

/*
 * Template/Alignment Placeholder (mmap/map)
 */
GT_INLINE gt_map_placeholder* gt_map_placeholder_new() {
  gt_map_placeholder* const map_placeholder = gt_alloc(gt_map_placeholder);
  gt_map_placeholder_clear(map_placeholder);
  return map_placeholder;
}
GT_INLINE void gt_map_placeholder_delete(gt_map_placeholder* const map_placeholder) {
  GT_NULL_CHECK(map_placeholder);
  gt_free(map_placeholder);
}
GT_INLINE void gt_map_placeholder_clear(gt_map_placeholder* const map_placeholder) {
  GT_NULL_CHECK(map_placeholder);
  memset(map_placeholder,0,sizeof(gt_map_placeholder));
}
GT_INLINE void gt_map_placeholder_set_sam_fields(gt_map_placeholder* const map_placeholder,
    const bool not_passing_QC,const bool PCR_duplicate,const uint32_t hard_trim_left,const uint32_t hard_trim_right) {
  GT_NULL_CHECK(map_placeholder);
  map_placeholder->not_passing_QC = not_passing_QC;
  map_placeholder->PCR_duplicate = PCR_duplicate;
  map_placeholder->hard_trim_left = hard_trim_left;
  map_placeholder->hard_trim_right = hard_trim_right;
}

/*
 * Compare Template/Alignment Placeholders _> (int (*)(const void *,const void *))
 * i.e. (Sorting from smaller to bigger, <) := a-b
 *   cmp(a,b) := -n if (a<b)
 *                n if (a>b)
 *                0 if (a==b)
 */
#define GT_MAP_PLACEHOLDER_CMP_CHECK_BOUNDARY_CONDITIONS(ph_map_a,ph_map_b) \
  /* Boundary conditions (Different type, same map, unmapped) */ \
  if (ph_map_a->type!=ph_map_b->type) return ((int)ph_map_a->type - (int)ph_map_b->type); /* Check different type */ \
  if (ph_map_a->map==ph_map_b->map) return 0; /* Check same map */ \
  if (ph_map_a->map==NULL || ph_map_b->map==NULL) return (ph_map_a->map==NULL) ? 1 : -1; /* Check unmapped */
int gt_map_placeholder_cmp_coordinates(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b) {
  // NOTE: Doesn't check chromosome
  // Boundary conditions (same map, unmapped)
  if (ph_map_a->map==ph_map_b->map) return 0; /* Check same map */ \
  if (ph_map_a->map==NULL || ph_map_b->map==NULL) return (ph_map_a->map==NULL) ? 1 : -1; /* Check unmapped */
  // Here, different map with same type and both mapped
  const int64_t coordinate_a = gt_map_get_global_coordinate(ph_map_a->map);
  const int64_t coordinate_b = gt_map_get_global_coordinate(ph_map_b->map);
  return coordinate_a-coordinate_b;
}

int gt_map_placeholder_cmp_map_distance(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b) {
  // TODO
  return 0;
}
int gt_map_placeholder_cmp_map_phred_scores(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b) {
  /*
   * Compare phred scores ( -10*log(Pr{mapping position is wrong},10) )
   *   - 0, 255 values are not allowed
   *   - 1..254 = bad_score..good_score
   */
  GT_MAP_PLACEHOLDER_CMP_CHECK_BOUNDARY_CONDITIONS(ph_map_a,ph_map_b);
  // Here, different map with same type and both mapped
  const int score_a = ph_map_a->map->phred_score;
  const int score_b = ph_map_b->map->phred_score;
  // No scored
  if (score_a==GT_MAP_NO_PHRED_SCORE || score_b==GT_MAP_NO_PHRED_SCORE) return score_a-score_b;
  // Regular compare (including case (score_a==0 || score_b==0))
  return score_b-score_a;
}
int gt_map_placeholder_cmp_map_gt_scores_ascending(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b) {
  GT_MAP_PLACEHOLDER_CMP_CHECK_BOUNDARY_CONDITIONS(ph_map_a,ph_map_b);
  // Here, different map with same type and both mapped
  const uint64_t score_a = ph_map_a->map->gt_score;
  const uint64_t score_b = ph_map_b->map->gt_score;
  return (score_a<score_b) ? -1 : ( (score_a>score_b) ? 1 : 0);
}
int gt_map_placeholder_cmp_map_gt_scores_descending(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b) {
  GT_MAP_PLACEHOLDER_CMP_CHECK_BOUNDARY_CONDITIONS(ph_map_a,ph_map_b);
  // Here, different map with same type and both mapped
  const uint64_t score_a = ph_map_a->map->gt_score;
  const uint64_t score_b = ph_map_b->map->gt_score;
  return (score_a<score_b) ? 1 : ( (score_a>score_b) ? -1 : 0);
}

int gt_map_placeholder_cmp_mmap_phred_scores(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b) {
  // TODO
  return 0;
}
int gt_map_placeholder_cmp_mmap_gt_scores_ascending(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b) {
  // TODO
  return 0;
}
int gt_map_placeholder_cmp_mmap_gt_scores_descending(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b) {
  // TODO
  return 0;
}

/*
 * Template/Alignment Placeholder (mmap/map)
 *
 * Fills the placeholder's vector @mmap_placeholder with all the mmaps from template
 *   - If @include_mate_placeholder is set, inserts a placeholder for the mate
 *       ({ph.paired_end.paired_end_position==1}). i.e. As for SAM output like operations
 *   - If @split_segments is set, individual segments are inserted into @mmap_placeholder (taking into account quimeras/segments)
 *   - If @primary_mmap_end1_pos is not null, sets it to the position in the array of the primary/best mmap(end/1) placeholder
 *       Sorts w.r.t @gt_ph_cmp_fx sorting. If @gt_ph_cmp_fx is null, then is set to the first one added.
 *   - If @primary_mmap_end2_pos is not null, sets it to the position in the array of the primary/best mmap(end/2) placeholder
 *       Sorts w.r.t @gt_ph_cmp_fx sorting. If @gt_ph_cmp_fx is null, then is set to the first one added.
 *       If @include_mate_placeholder is not set, then @primary_mmap_end1_pos==@primary_mmap_end2_pos
 *   - If @placeholder_template is not null, then its values are used as defaults (not_passing_QC,PCR_duplicate,...)
 */
#define GT_MAP_PLACEHOLDER_SET_BEST_PH(best_mmap_ph,best_mmap_ph_position,mmap_ph,mmap_ph_position) \
  best_mmap_ph = mmap_ph; \
  best_mmap_ph_position = mmap_ph_position
#define GT_MAP_PLACEHOLDER_CMP_BEST_PH(best_mmap_ph,best_mmap_ph_position,mmap_ph,mmap_ph_position) { \
  if (gt_expect_false(best_mmap_ph_position==UINT64_MAX)) { \
    GT_MAP_PLACEHOLDER_SET_BEST_PH(best_mmap_ph,best_mmap_ph_position,mmap_ph,mmap_ph_position); \
  } else if (gt_ph_cmp_fx!=NULL) { \
    if (gt_ph_cmp_fx(&best_mmap_ph,&mmap_ph)>0) { \
      GT_MAP_PLACEHOLDER_SET_BEST_PH(best_mmap_ph,best_mmap_ph_position,mmap_ph,mmap_ph_position); \
    } \
  } \
}
GT_INLINE void gt_map_placeholder_add_mmap(
    gt_map* const map_endA,gt_map* const map_endB,gt_string* const read_endA,const uint64_t paired_end_position,
    gt_vector* const mmap_placeholder,const bool split_segments,
    int (*gt_ph_cmp_fx)(gt_map_placeholder* const,gt_map_placeholder* const),const bool cmp_with_best,
    gt_map_placeholder* const best_mmap_ph,uint64_t* const best_mmap_ph_position,
    gt_map_placeholder* const mmap_ph) {
  uint64_t num_placeholders = gt_vector_get_used(mmap_placeholder);
  // Note that (map_endA!=NULL) || (map_endB!=NULL) must hold
  mmap_ph->paired_end.paired_end_position = paired_end_position;
  if (map_endA==NULL) {
    /*
     * End/1 Unmapped
     *   MAP => Don't unfold quimeras (never needed)
     *   SAM => Don't unfold quimeras (Is unmapped, so unfold is not required)
     */
    mmap_ph->type = GT_MMAP_PLACEHOLDER_UNPAIRED;
    mmap_ph->map = NULL;
    mmap_ph->paired_end.mate = map_endB; // SAM doen't care of this (as it's unmapped). But MAP does.
    mmap_ph->hard_trim_left = 0;
    mmap_ph->hard_trim_right = 0;
    gt_vector_insert(mmap_placeholder,*mmap_ph,gt_map_placeholder);
    // Pick primary alignment
    if (cmp_with_best) GT_MAP_PLACEHOLDER_CMP_BEST_PH(*best_mmap_ph,*best_mmap_ph_position,*mmap_ph,num_placeholders);
  } else {
    /*
     * End/1 Mapped
     */
    if (map_endB==NULL) {
      /*
       * End/2 Unmapped
       */
      mmap_ph->type = GT_MMAP_PLACEHOLDER_UNPAIRED;
      GT_MAP_SEGMENT_ITERATOR(map_endA,map_segment_iterator_end1) {
        mmap_ph->map = gt_map_segment_iterator_get_map(&map_segment_iterator_end1);
        mmap_ph->paired_end.mate = NULL;
        mmap_ph->hard_trim_left = (split_segments) ? gt_map_segment_iterator_get_accumulated_offset(&map_segment_iterator_end1) : 0;
        mmap_ph->hard_trim_right = (split_segments) ? gt_map_segment_iterator_get_remaining_bases(&map_segment_iterator_end1,read_endA) : 0;
        gt_vector_insert(mmap_placeholder,*mmap_ph,gt_map_placeholder);
        // Pick primary alignment
        if (cmp_with_best) GT_MAP_PLACEHOLDER_CMP_BEST_PH(*best_mmap_ph,*best_mmap_ph_position,*mmap_ph,num_placeholders);
        ++num_placeholders;
        if (!split_segments) break; // Break if quimeras are not supposed to be unfolded
      }
    } else {
      /*
       * End/2 Mapped
       */
      mmap_ph->type = GT_MMAP_PLACEHOLDER_PAIRED;
      GT_MAP_SEGMENT_ITERATOR(map_endA,map_segment_iterator_end1) {
        GT_MAP_SEGMENT_ITERATOR(map_endB,map_segment_iterator_end2) {
          mmap_ph->map = gt_map_segment_iterator_get_map(&map_segment_iterator_end1);
          mmap_ph->paired_end.mate = gt_map_segment_iterator_get_map(&map_segment_iterator_end2);
          mmap_ph->hard_trim_left = (split_segments) ? gt_map_segment_iterator_get_accumulated_offset(&map_segment_iterator_end1) : 0;
          mmap_ph->hard_trim_right = (split_segments) ? gt_map_segment_iterator_get_remaining_bases(&map_segment_iterator_end1,read_endA) : 0;
          gt_vector_insert(mmap_placeholder,*mmap_ph,gt_map_placeholder);
          // Pick primary alignment
          if (cmp_with_best) GT_MAP_PLACEHOLDER_CMP_BEST_PH(*best_mmap_ph,*best_mmap_ph_position,*mmap_ph,num_placeholders);
          ++num_placeholders;
          if (!split_segments) break; // Break if quimeras are not supposed to be unfolded
        }
        if (!split_segments) break; // Break if quimeras are not supposed to be unfolded
      }
    }
  }
}
GT_INLINE void gt_map_placeholder_build_from_template(
    gt_template* const template,gt_vector* const mmap_placeholder,
    const bool include_mate_placeholder,const bool split_segments,const uint64_t max_num_maps,
    int (*gt_ph_cmp_fx)(gt_map_placeholder* const,gt_map_placeholder* const),
    uint64_t* const primary_mmap_end1_pos,uint64_t* const primary_mmap_end2_pos,
    gt_map_placeholder* const placeholder_template) {
  GT_TEMPLATE_CHECK(template);
  GT_VECTOR_CHECK(mmap_placeholder);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    placeholder_template->single_end.template = template;
    gt_map_placeholder_build_from_alignment(alignment,mmap_placeholder,split_segments,max_num_maps,
        gt_ph_cmp_fx,primary_mmap_end1_pos,placeholder_template);
    if (primary_mmap_end1_pos!=NULL && primary_mmap_end2_pos!=NULL) *primary_mmap_end2_pos=*primary_mmap_end1_pos;
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Prepare placeholder template
  gt_map_placeholder mmap_ph;
  if (placeholder_template!=NULL) mmap_ph=*placeholder_template;
  mmap_ph.paired_end.template = template;
  mmap_ph.secondary_alignment = true;
  // Gather initial data
  const uint64_t num_initial_placeholders = gt_vector_get_used(mmap_placeholder);
  gt_string* const read[2] = {gt_template_get_end1(template)->read,gt_template_get_end2(template)->read};
  // Add mmaps
  uint64_t num_mmaps_added = 0;
  gt_map_placeholder best_mmap_ph_end1, best_mmap_ph_end2;
  uint64_t best_mmap_ph_end1_pos = UINT64_MAX, best_mmap_ph_end2_pos = UINT64_MAX;
  if (gt_template_get_num_mmaps(template)==0) { // Unmmapped
    // Include unampped mmap
    mmap_ph.paired_end.mmap_attributes = NULL;
    gt_map_placeholder_add_mmap(NULL,NULL,read[0],0,
        mmap_placeholder,split_segments,gt_ph_cmp_fx,primary_mmap_end1_pos!=NULL,
        &best_mmap_ph_end1,&best_mmap_ph_end1_pos,&mmap_ph); // End/1
    if (include_mate_placeholder) {
      gt_map_placeholder_add_mmap(NULL,NULL,read[0],1,
          mmap_placeholder,split_segments,gt_ph_cmp_fx,primary_mmap_end2_pos!=NULL,
          &best_mmap_ph_end2,&best_mmap_ph_end2_pos,&mmap_ph); // End/2
    }
  } else {
    // Include mmaps
    bool end1_record_included = false, end2_record_included = false;
    GT_TEMPLATE_ITERATE_MMAP__ATTR_(template,mmap,mmap_attr) {
      mmap_ph.paired_end.mmap_attributes = mmap_attr;
      if (include_mate_placeholder) {
        if (mmap[0]!=NULL) {
          end1_record_included = true;
          gt_map_placeholder_add_mmap(mmap[0],mmap[1],read[0],0,
              mmap_placeholder,split_segments,gt_ph_cmp_fx,primary_mmap_end1_pos!=NULL,
              &best_mmap_ph_end1,&best_mmap_ph_end1_pos,&mmap_ph); // End/1
        }
        if (mmap[1]!=NULL) {
          end2_record_included = true;
          gt_map_placeholder_add_mmap(mmap[1],mmap[0],read[1],1,
              mmap_placeholder,split_segments,gt_ph_cmp_fx,primary_mmap_end2_pos!=NULL,
              &best_mmap_ph_end2,&best_mmap_ph_end2_pos,&mmap_ph); // End/2
        }
      } else {
        gt_map_placeholder_add_mmap(mmap[0],mmap[1],read[0],0,
            mmap_placeholder,split_segments,gt_ph_cmp_fx,primary_mmap_end1_pos!=NULL,
            &best_mmap_ph_end1,&best_mmap_ph_end1_pos,&mmap_ph); // End/1
      }
      if (++num_mmaps_added > max_num_maps) break;
    }
    // Check if at least we have a record per each end
    if (include_mate_placeholder) {
      if (!end1_record_included) {
        gt_map_placeholder_add_mmap(NULL,NULL,read[0],0,
            mmap_placeholder,split_segments,gt_ph_cmp_fx,primary_mmap_end1_pos!=NULL,
            &best_mmap_ph_end1,&best_mmap_ph_end1_pos,&mmap_ph); // End/1
      }
      if (!end2_record_included) {
        gt_map_placeholder_add_mmap(NULL,NULL,read[1],1,
            mmap_placeholder,split_segments,gt_ph_cmp_fx,primary_mmap_end2_pos!=NULL,
            &best_mmap_ph_end2,&best_mmap_ph_end2_pos,&mmap_ph); // End/2
      }
    }
  }
  // Set primary alignment
  if (primary_mmap_end1_pos!=NULL) {
    *primary_mmap_end1_pos = best_mmap_ph_end1_pos;
    gt_vector_get_elm(mmap_placeholder,best_mmap_ph_end1_pos,gt_map_placeholder)->secondary_alignment = false;
    GT_SWAP(*(gt_vector_get_elm(mmap_placeholder,best_mmap_ph_end1_pos,gt_map_placeholder)),
            *(gt_vector_get_elm(mmap_placeholder,num_initial_placeholders,gt_map_placeholder)));
  }
  if (primary_mmap_end2_pos!=NULL) {
    if (include_mate_placeholder) {
      *primary_mmap_end2_pos = best_mmap_ph_end2_pos;
      gt_vector_get_elm(mmap_placeholder,best_mmap_ph_end2_pos,gt_map_placeholder)->secondary_alignment = false;
      GT_SWAP(*(gt_vector_get_elm(mmap_placeholder,best_mmap_ph_end2_pos,gt_map_placeholder)),
              *(gt_vector_get_elm(mmap_placeholder,num_initial_placeholders+1,gt_map_placeholder))); // Courtesy
    } else {
      if (primary_mmap_end1_pos!=NULL) *primary_mmap_end2_pos = *primary_mmap_end1_pos;
    }
  }
}
/*
 * Fills the placeholder's vector @mmap_placeholder with all the maps from alignment
 *   - If @split_segments is set, individual segments are inserted into @mmap_placeholder (taking into account quimeras/segments)
 *   - If @primary_map_position is not null, sets it to the position in the array of the primary/best alignment placeholder
 *       Sorts w.r.t @gt_ph_cmp_fx sorting. If @gt_ph_cmp_fx is null, then is set to the first one added.
 *       If @include_mate_placeholder is not set, then @primary_mmap_end1_pos==@primary_mmap_end2_pos
 *   - If @placeholder_template is not null, then its values are used as defaults (not_passing_QC,PCR_duplicate,...)
 */
GT_INLINE void gt_map_placeholder_add_map(
    gt_map* const map,gt_string* const read,
    gt_vector* const mmap_placeholder,const bool split_segments,
    int (*gt_ph_cmp_fx)(gt_map_placeholder* const,gt_map_placeholder* const),const bool cmp_with_best,
    gt_map_placeholder* const best_mmap_ph,uint64_t* const best_mmap_ph_position,
    gt_map_placeholder* const mmap_ph) {
  uint64_t num_placeholders = gt_vector_get_used(mmap_placeholder);
  mmap_ph->type = GT_MAP_PLACEHOLDER;
  if (map==NULL) { // Unmapped
    mmap_ph->map = NULL;
    mmap_ph->hard_trim_left = 0;
    mmap_ph->hard_trim_right = 0;
    gt_vector_insert(mmap_placeholder,*mmap_ph,gt_map_placeholder);
    // Pick primary alignment
    if (cmp_with_best) GT_MAP_PLACEHOLDER_CMP_BEST_PH(*best_mmap_ph,*best_mmap_ph_position,*mmap_ph,num_placeholders);
  } else {
    GT_MAP_SEGMENT_ITERATOR(map,map_segment_iterator) {
      mmap_ph->map = gt_map_segment_iterator_get_map(&map_segment_iterator);
      mmap_ph->hard_trim_left = (split_segments) ? gt_map_segment_iterator_get_accumulated_offset(&map_segment_iterator) : 0;
      mmap_ph->hard_trim_right = (split_segments) ? gt_map_segment_iterator_get_remaining_bases(&map_segment_iterator,read) : 0;
      gt_vector_insert(mmap_placeholder,*mmap_ph,gt_map_placeholder);
      // Pick primary alignment
      if (best_mmap_ph_position!=NULL) {
        GT_MAP_PLACEHOLDER_CMP_BEST_PH(*best_mmap_ph,*best_mmap_ph_position,*mmap_ph,num_placeholders);
      }
      ++num_placeholders;
      // Break if quimeras are not supposed to be unfolded
      if (!split_segments) break;
    }
  }
}
GT_INLINE void gt_map_placeholder_build_from_alignment(
    gt_alignment* const alignment,gt_vector* const mmap_placeholder,const bool split_segments,const uint64_t max_num_maps,
    int (*gt_ph_cmp_fx)(gt_map_placeholder* const,gt_map_placeholder* const),uint64_t* const primary_map_position,
    gt_map_placeholder* const placeholder_template) {
  GT_ALIGNMENT_CHECK(alignment);
  GT_VECTOR_CHECK(mmap_placeholder);
  // Prepare placeholder template
  gt_map_placeholder mmap_ph;
  if (placeholder_template!=NULL) mmap_ph=*placeholder_template;
  const uint64_t num_maps = gt_alignment_get_num_maps(alignment);
  const uint64_t num_initial_placeholders = gt_vector_get_used(mmap_placeholder);
  mmap_ph.single_end.alignment = alignment;
  mmap_ph.secondary_alignment = true;
  // Add maps
  uint64_t num_mmaps_added = 0;
  gt_map_placeholder best_mmap_ph;
  uint64_t best_mmap_ph_pos = UINT64_MAX;
  if (num_maps==0) {
    // Include unmapped
    gt_map_placeholder_add_map(NULL,alignment->read,mmap_placeholder,split_segments,
        gt_ph_cmp_fx,primary_map_position!=NULL,&best_mmap_ph,&best_mmap_ph_pos,&mmap_ph);
  } else {
    // Include maps
    GT_ALIGNMENT_ITERATE(alignment,map) {
      gt_map_placeholder_add_map(map,alignment->read,mmap_placeholder,split_segments,
          gt_ph_cmp_fx,primary_map_position!=NULL,&best_mmap_ph,&best_mmap_ph_pos,&mmap_ph);
      if (++num_mmaps_added > max_num_maps) break;
    }
  }
  // Set primary alignment
  if (primary_map_position!=NULL) {
    *primary_map_position = best_mmap_ph_pos;
    gt_vector_get_elm(mmap_placeholder,best_mmap_ph_pos,gt_map_placeholder)->secondary_alignment = false;
    GT_SWAP(*(gt_vector_get_elm(mmap_placeholder,best_mmap_ph_pos,gt_map_placeholder)),
            *(gt_vector_get_elm(mmap_placeholder,num_initial_placeholders,gt_map_placeholder))); // Courtesy
  }
}

