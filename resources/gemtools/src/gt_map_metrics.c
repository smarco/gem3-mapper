/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_metrics.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_map_metrics.h"

/*
 * Metrics :: Local(over first block) / Global(over all blocks)
 */
// Length of the mapping into the reference sequence (i.e. size over the chromosome)
GT_INLINE uint64_t gt_map_get_length(gt_map* const map) {
  GT_MAP_CHECK(map);
  // Indels are taken into account to calculate the length
  int64_t length = map->base_length;
  GT_MISMS_ITERATE(map,misms_it) {
    switch (misms_it->misms_type) {
      case MISMS: break;
      case INS:
        length += gt_misms_get_size(misms_it);
        break;
      case DEL:
        length -= gt_misms_get_size(misms_it);
        break;
    }
  }
  gt_cond_fatal_error(length<0,MAP_NEG_LENGTH);
  return (uint64_t)length;
}
GT_INLINE uint64_t gt_map_get_segment_length(gt_map* const map) {
  GT_MAP_CHECK(map);
  uint64_t length = 0;
  GT_MAP_SEGMENT_ITERATE(map,map_it) {
    length += gt_map_get_length(map_it) + gt_map_get_junction_size(map_it);
  }
  return length;
}
GT_INLINE uint64_t gt_map_get_global_length(gt_map* const map) {
  GT_MAP_CHECK(map);
  uint64_t length = 0;
  GT_MAP_ITERATE(map,map_it) {
    length += gt_map_get_length(map_it) + gt_map_get_junction_size(map_it);
  }
  return length;
}
GT_INLINE uint64_t gt_map_get_global_base_length(gt_map* const map) {
  GT_MAP_CHECK(map);
  uint64_t base_length = 0;
  GT_MAP_ITERATE(map,map_it) {
    base_length += gt_map_get_base_length(map_it);
  }
  return base_length;
}
// GEM Distance (Number of Mismatches/Insert/Delete/Split operations)
GT_INLINE uint64_t gt_map_get_distance(gt_map* const map) {
  GT_MAP_CHECK(map);
  return gt_vector_get_used(map->mismatches);
}
GT_INLINE uint64_t gt_map_get_segment_distance(gt_map* const map) {
  GT_MAP_CHECK(map);
  uint64_t distance = 0;
  GT_MAP_SEGMENT_ITERATE(map,map_it) {
    distance += gt_map_get_distance(map_it) + (gt_map_has_next_block(map_it)?1:0);
  }
  return distance;
}
GT_INLINE uint64_t gt_map_get_global_distance(gt_map* const map) {
  GT_MAP_CHECK(map);
  uint64_t distance = 0;
  GT_MAP_ITERATE(map,map_it) {
    distance += gt_map_get_distance(map_it) + (gt_map_has_next_block(map_it)?1:0);
  }
  return distance;
}
GT_INLINE uint64_t gt_map_get_no_split_distance(gt_map* const map){
  return gt_map_get_global_distance(map) - (gt_map_get_num_blocks(map)-1);
}
// SWG-Distance
GT_INLINE int64_t gt_map_get_swg_score(gt_map* const map) {
  GT_MAP_CHECK(map);
  // Indels are taken into account to calculate the length
  int64_t swg_score = map->base_length;
  GT_MISMS_ITERATE(map,misms_it) {
    switch (misms_it->misms_type) {
      case MISMS:
        swg_score -= 4;
        break;
      case INS:
      case DEL:
        swg_score -= 6 + gt_misms_get_size(misms_it);
        break;
    }
  }
  return swg_score;
}
// Bases aligned (Mismatches not included)
GT_INLINE uint64_t gt_map_get_bases_aligned(gt_map* const map) {
  GT_MAP_CHECK(map);
  int64_t bases_aligned = map->base_length;
  GT_MISMS_ITERATE(map,misms_it) {
    switch (misms_it->misms_type) {
      case INS:
        break;
      case MISMS:
        --bases_aligned;
        break;
      case DEL:
        bases_aligned -= gt_misms_get_size(misms_it);
        break;
    }
  }
  gt_cond_fatal_error(bases_aligned<0,MAP_NEG_MAPPED_BASES);
  return (uint64_t)bases_aligned;
}
GT_INLINE uint64_t gt_map_get_segment_bases_aligned(gt_map* const map) {
  GT_MAP_CHECK(map);
  uint64_t bases_aligned = 0;
  GT_MAP_SEGMENT_ITERATE(map,map_it) {
    bases_aligned += gt_map_get_bases_aligned(map_it);
  }
  return bases_aligned;
}
GT_INLINE uint64_t gt_map_get_global_bases_aligned(gt_map* const map) {
  GT_MAP_CHECK(map);
  uint64_t bases_aligned = 0;
  GT_MAP_ITERATE(map,map_it) {
    bases_aligned += gt_map_get_bases_aligned(map_it);
  }
  return bases_aligned;
}
/*
 * Distance procedures
 */
GT_INLINE uint64_t gt_map_get_levenshtein_distance(gt_map* const map) {
  GT_MAP_CHECK(map);
  uint64_t lev_distance = 0;
  GT_MISMS_ITERATE(map,misms_it) {
    switch (misms_it->misms_type) {
      case MISMS:
        ++lev_distance;
        break;
      case INS:
      case DEL:
        lev_distance += gt_misms_get_size(misms_it);
        break;
    }
  }
  return lev_distance;
}
GT_INLINE uint64_t gt_map_get_segment_levenshtein_distance(gt_map* const map) {
  GT_MAP_CHECK(map);
  uint64_t lev_distance = 0;
  GT_MAP_SEGMENT_ITERATE(map,map_it) {
    lev_distance += gt_map_get_levenshtein_distance(map_it);
  }
  return lev_distance;
}
GT_INLINE uint64_t gt_map_get_global_levenshtein_distance(gt_map* const map) {
  GT_MAP_CHECK(map);
  uint64_t lev_distance = 0;
  GT_MAP_ITERATE(map,map_it) {
    lev_distance += gt_map_get_levenshtein_distance(map_it);
  }
  return lev_distance;
}
/*
 * MMap/Vector Metrics
 */
// MMap based (paired mapping, ...)
GT_INLINE uint64_t gt_mmap_get_global_distance(gt_map** const mmap,const uint64_t num_blocks) {
  GT_ZERO_CHECK(num_blocks);
  uint64_t i, distance=0;
  for (i=0;i<num_blocks;++i) {
    distance += (mmap[i]!=NULL) ? gt_map_get_global_distance(mmap[i]) : 0;
  }
  return distance;
}
GT_INLINE uint64_t gt_mmap_get_global_levenshtein_distance(gt_map** const mmap,const uint64_t num_blocks) {
  GT_ZERO_CHECK(num_blocks);
  uint64_t i, distance=0;
  for (i=0;i<num_blocks;++i) {
    distance += (mmap[i]!=NULL) ? gt_map_get_global_levenshtein_distance(mmap[i]) : 0;
  }
  return distance;
}
// Vector based ( Metrics out of a set of maps )
GT_INLINE uint64_t gt_map_vector_get_length(gt_vector* const maps) {
  GT_NULL_CHECK(maps);
  uint64_t length = 0;
  GT_VECTOR_ITERATE(maps,map,map_pos,gt_map*) {
    length += (maps!=NULL) ? gt_map_get_global_length(*map) : 0;
  }
  return length;
}
GT_INLINE uint64_t gt_map_vector_get_distance(gt_vector* const maps) {
  GT_NULL_CHECK(maps);
  uint64_t distance = 0;
  GT_VECTOR_ITERATE(maps,map,map_pos,gt_map*) {
    distance += (maps!=NULL) ? gt_map_get_global_distance(*map) : 0;
  }
  return distance;
}
// Insert/Template Size
GT_INLINE void gt_map_get_last_block__acc_base_length(gt_map* const map,gt_map** const last_block,uint64_t* const acc_base_length) {
  uint64_t _acc_base_length = 0;
  gt_map* _last_block = map;
  GT_MAP_SEGMENT_ITERATE(map,map_block) {
    _last_block=map_block;
    _acc_base_length+=gt_map_get_base_length(map_block);
  }
  *last_block = _last_block;
  *acc_base_length = _acc_base_length;
}
GT_INLINE int64_t gt_map_get_observed_template_size(gt_map* const map_a,gt_map* const map_b) {
  GT_MAP_CHECK(map_a);
  GT_MAP_CHECK(map_b);
  if (gt_expect_false(!gt_string_equals(map_a->seq_name,map_b->seq_name))) return 0;
  gt_map *right_block_a, *right_block_b;
  gt_map *left_block_a,  *left_block_b;
  uint64_t map_length_a, map_length_b;
  // Get rightmost/leftmost coordinates
  if (gt_map_get_strand(map_a)==FORWARD) {
    right_block_a = map_a;
    gt_map_get_last_block__acc_base_length(map_a,&left_block_a,&map_length_a);
  } else {
    left_block_a = map_a;
    gt_map_get_last_block__acc_base_length(map_a,&right_block_a,&map_length_a);
  }
  if (gt_map_get_strand(map_b)==FORWARD) {
    right_block_b = map_b;
    gt_map_get_last_block__acc_base_length(map_b,&left_block_b,&map_length_b);
  } else {
    left_block_b = map_b;
    gt_map_get_last_block__acc_base_length(map_b,&right_block_b,&map_length_b);
  }
  const uint64_t rightmost_a = gt_map_get_position(right_block_a)+gt_map_get_length(right_block_a);
  const uint64_t rightmost_b = gt_map_get_position(right_block_b)+gt_map_get_length(right_block_b);
  const uint64_t leftmost_a = gt_map_get_position(left_block_a);
  const uint64_t leftmost_b = gt_map_get_position(left_block_b);
  const uint64_t rightmost = GT_MAX(rightmost_a,rightmost_b);
  if (leftmost_a < leftmost_b) { // Determine the leftmost map (thus the sign of the return value)
    return (rightmost-leftmost_a);
  } else {
    return -(rightmost-leftmost_b);
  }
}
/*
 * Strict Map compare functions
 *   1.- Same sequence name
 *   2.- Same strand
 *   3.- Same number of blocks
 *   4.- for mapBlock in map {
 *         4.1 Same begin position
 *         4.2 Same end position
 *       }
 */
GT_INLINE int64_t gt_map_sm_cmp(gt_map* const map_1,gt_map* const map_2) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2);
  const int64_t begin_distance = ((int64_t)gt_map_get_begin_mapping_position(map_1)) - ((int64_t)gt_map_get_begin_mapping_position(map_2));
  if (begin_distance != 0) return 1;
  const int64_t end_distance = ((int64_t)gt_map_get_end_mapping_position(map_1)) - (int64_t)(gt_map_get_end_mapping_position(map_2));
  if (end_distance != 0) return 1;
  if (map_1->next_block.map==NULL && map_2->next_block.map==NULL) return 0;
  if (map_1->next_block.map!=NULL && map_2->next_block.map!=NULL) {
    return gt_map_sm_cmp(map_1->next_block.map,map_2->next_block.map);
  }
  return 1;
}
GT_INLINE int64_t gt_map_cmp(gt_map* const map_1,gt_map* const map_2) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2);
  if (gt_string_cmp(map_1->seq_name,map_2->seq_name)!=0) {
    return 1;
  } else {
    if (map_1->strand==map_2->strand) {
      return gt_map_sm_cmp(map_1,map_2);
    } else {
      return 1;
    }
  }
}
/*
 * General Map compare function (differs from the standard one at @gt_map_cmp)
 *   1.- Same sequence name
 *   2.- Same strand
 *   3.- Same number of blocks
 *   4.- If (numBlocks == 1) {
 *           4.1 |begin_a-begin_b|<=range_tolerated || |end_a-end_b|<=range_tolerated
 *       } else {
 *         for mapBlock in map {
 *           5.1 |begin_a-begin_b|<=range_tolerated
 *           5.2 |end_a-end_b|<=range_tolerated
 *         }
 *       }
 */
GT_INLINE int64_t gt_map_range_sm_cmp(gt_map* const map_1,gt_map* const map_2,const uint64_t range_tolerated,const uint64_t num_maps_left) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2);
  const int64_t begin_distance = ((int64_t)gt_map_get_begin_mapping_position(map_1)) - ((int64_t)gt_map_get_begin_mapping_position(map_2));
  const int64_t end_distance = ((int64_t)gt_map_get_end_mapping_position(map_1)) - (int64_t)(gt_map_get_end_mapping_position(map_2));
  if (GT_ABS(begin_distance)<=range_tolerated && GT_ABS(end_distance)<=range_tolerated) {
    return (num_maps_left==1) ? 0 : gt_map_range_sm_cmp(map_1->next_block.map,map_2->next_block.map,range_tolerated,num_maps_left-1);
  } else {
    return GT_ABS(begin_distance)+GT_ABS(end_distance);
  }
}
GT_INLINE int64_t gt_map_range_cmp(gt_map* const map_1,gt_map* const map_2,const uint64_t range_tolerated) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2);
  int64_t cmp_tags = gt_string_cmp(map_1->seq_name,map_2->seq_name);
  if (cmp_tags!=0) {
    return cmp_tags;
  } else {
    if (map_1->strand==map_2->strand) {
      const uint64_t num_blocks_map_1 = gt_map_get_num_blocks(map_1);
      const uint64_t num_blocks_map_2 = gt_map_get_num_blocks(map_2);
      if (num_blocks_map_1==num_blocks_map_2) {
        if (num_blocks_map_1==1) { // Standard Mapping
          const int64_t begin_distance = ((int64_t)gt_map_get_begin_mapping_position(map_1)) - ((int64_t)gt_map_get_begin_mapping_position(map_2));
          const int64_t end_distance = ((int64_t)gt_map_get_end_mapping_position(map_1)) - (int64_t)(gt_map_get_end_mapping_position(map_2));
          if (GT_ABS(begin_distance)<=range_tolerated || GT_ABS(end_distance)<=range_tolerated) return 0;
          return GT_ABS(begin_distance)+GT_ABS(end_distance);
        } else { // Split Maps Involved
          return gt_map_range_sm_cmp(map_1,map_2,range_tolerated,num_blocks_map_1);
        }
      } else { // Different splits
        return (num_blocks_map_1-num_blocks_map_2);
      }
    } else {
      return map_1->strand==FORWARD ? 1 : -1;
    }
  }
}
/*
 * CIGAR Map compare function
 */
GT_INLINE int64_t gt_map_cmp_cigar(gt_map* const map_1,gt_map* const map_2) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2);
  const uint64_t num_misms_1 = gt_vector_get_used(map_1->mismatches);
  const uint64_t num_misms_2 = gt_vector_get_used(map_2->mismatches);
  if (num_misms_1 != num_misms_2) {
    return 1;
  } else {
    gt_misms* const misms_1 = gt_vector_get_mem(map_1->mismatches,gt_misms);
    gt_misms* const misms_2 = gt_vector_get_mem(map_1->mismatches,gt_misms);
    uint64_t i;
    for (i=0;i<num_misms_1;++i) {
      if (misms_1[i].misms_type != misms_2[i].misms_type ||
          misms_1[i].position != misms_2[i].position) return 1;
      if (misms_1[i].misms_type == MISMS &&
          misms_1[i].base != misms_2[i].base) return 1;
      if (misms_1[i].size != misms_2[i].size) return 1;
    }
  }
  return 0;
}
/*
 * MMap compare functions (based on Map compare functions)
 */
GT_INLINE int64_t gt_mmap_cmp(gt_map** const map_1,gt_map** const map_2,const uint64_t num_maps) {
  GT_NULL_CHECK(map_1); GT_NULL_CHECK(map_2);
  uint64_t i;
  for (i=0;i<num_maps;++i) {
    const int64_t map_cmp = gt_map_cmp(map_1[i],map_2[i]);
    if (map_cmp!=0) return map_cmp;
  }
  return 0;
}
GT_INLINE int64_t gt_mmap_range_cmp(
    gt_map** const map_1,gt_map** const map_2,const uint64_t num_maps,const uint64_t range_tolerated) {
  GT_NULL_CHECK(map_1); GT_NULL_CHECK(map_2);
  uint64_t i;
  for (i=0;i<num_maps;++i) {
    const int64_t map_cmp = gt_map_range_cmp(map_1[i],map_2[i],range_tolerated);
    if (map_cmp!=0) return map_cmp;
  }
  return 0;
}
GT_INLINE bool gt_map_less_than(gt_map* const map_1,gt_map* const map_2) {
  GT_MAP_CHECK(map_1); GT_MAP_CHECK(map_2);
  if (gt_map_get_global_bases_aligned(map_1) < gt_map_get_global_bases_aligned(map_2)) return true;
  if (gt_map_get_num_indels(map_1) < gt_map_get_num_indels(map_2)) return true;
  if (gt_map_get_global_levenshtein_distance(map_1) < gt_map_get_global_levenshtein_distance(map_2)) return true;
  return false;
}
