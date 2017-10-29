/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */


#include "gt_map.h"

#define GT_MAP_NUM_INITIAL_MISMS 4
#define GT_MAP_INITIAL_SEQ_NAME_SIZE 10

/*
 * Setup
 */
GT_INLINE gt_map* gt_map_new() {
  gt_map* map = gt_alloc(gt_map);
  map->seq_name = gt_string_new(GT_MAP_INITIAL_SEQ_NAME_SIZE);
  map->position = 0;
  map->base_length = 0;
  map->gt_score = GT_MAP_NO_GT_SCORE;
  map->phred_score = GT_MAP_NO_PHRED_SCORE;
  map->mismatches = gt_vector_new(GT_MAP_NUM_INITIAL_MISMS,sizeof(gt_misms));
  map->next_block.map = NULL;
  map->attributes = NULL;
  return map;
}
GT_INLINE void gt_map_clear(gt_map* const map) {
  GT_MAP_CHECK(map);
  gt_string_clear(map->seq_name);
  map->position = 0;
  map->base_length = 0;
  map->gt_score = GT_MAP_NO_GT_SCORE;
  map->phred_score = GT_MAP_NO_PHRED_SCORE;
  gt_map_clear_misms(map);
  map->next_block.map = NULL;
  if (map->attributes!=NULL) gt_attributes_clear(map->attributes);
}
GT_INLINE void gt_map_block_delete(gt_map* const map) {
  GT_MAP_CHECK(map);
  gt_string_delete(map->seq_name);
  gt_vector_delete(map->mismatches);
  if (map->attributes!=NULL) gt_attributes_delete(map->attributes);
  gt_free(map);
}
GT_INLINE void gt_map_delete(gt_map* const map) {
  GT_MAP_CHECK(map);
  if (map->next_block.map != NULL) gt_map_delete(map->next_block.map);
  gt_map_block_delete(map);
}
/*
 * MapBlock Accessors
 */
GT_INLINE char* gt_map_get_seq_name(gt_map* const map) {
  GT_MAP_CHECK(map);
  return gt_string_get_string(map->seq_name);
}
GT_INLINE uint64_t gt_map_get_seq_name_length(gt_map* const map) {
  GT_MAP_CHECK(map);
  return gt_string_get_length(map->seq_name);
}
GT_INLINE gt_string* gt_map_get_string_seq_name(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->seq_name;
}
GT_INLINE void gt_map_set_seq_name(gt_map* const map,const char* const seq_name,const uint64_t length) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(seq_name);
  gt_string_set_nstring_static(map->seq_name,seq_name,length);
}
GT_INLINE void gt_map_set_string_seq_name(gt_map* const map,gt_string* const seq_name) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(seq_name);
  gt_string_set_nstring(map->seq_name,gt_string_get_string(seq_name),gt_string_get_length(seq_name));
}
GT_INLINE gt_strand gt_map_get_strand(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->strand;
}
GT_INLINE void gt_map_set_strand(gt_map* const map,const gt_strand strand) {
  GT_MAP_CHECK(map);
  map->strand = strand;
}
GT_INLINE uint64_t gt_map_get_base_length(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->base_length;
}
GT_INLINE void gt_map_set_base_length(gt_map* const map,const uint64_t length) {
  GT_MAP_CHECK(map);
  map->base_length = length;
}
// Most right-hand position of the current map block (NOT GLOBAL COORDINATE --> gt_map_get_global_coordinate())
GT_INLINE uint64_t gt_map_get_position(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->position;
}
GT_INLINE void gt_map_set_position(gt_map* const map,const uint64_t position) {
  GT_MAP_CHECK(map);
  map->position = position;
}

//Utils
GT_INLINE uint64_t gt_map_get_min_intron_length(gt_map* const map) {
  uint64_t min_intron_length = UINT64_MAX;
  GT_MAP_ITERATE(map, map_block){
    if(gt_map_has_next_block(map_block)){
      int64_t l = gt_map_get_junction_size(map_block);
      if(l>=0 && min_intron_length > l){
        min_intron_length = l;
      }
    }
  }
  return min_intron_length;
}
GT_INLINE uint64_t gt_map_get_min_block_length(gt_map* const map) {
  uint64_t min_block_length = UINT64_MAX;
  GT_MAP_ITERATE(map, map_block){
    int64_t l = gt_map_get_base_length(map_block);
    if(l>=0 && min_block_length > l){
      min_block_length = l;
    }
  }
  return min_block_length;
}

// Trim helpers
GT_INLINE uint64_t gt_map_get_left_trim_length(gt_map* const map) {
  if (gt_map_get_num_misms(map)>0) {
    gt_misms* const first_misms = gt_map_get_misms(map,0);
    if (first_misms->position==0 && first_misms->misms_type==DEL) return first_misms->size;
  }
  return 0;
}
GT_INLINE uint64_t gt_map_get_right_trim_length(gt_map* const map) {
  const uint64_t num_misms = gt_map_get_num_misms(map);
  if (num_misms>0) {
    gt_misms* const last_misms = gt_map_get_misms(map,num_misms-1);
    if (last_misms->position+last_misms->size==map->base_length && last_misms->misms_type==DEL) return last_misms->size;
  }
  return 0;
}
// Begin/End Mapping Position
GT_INLINE uint64_t gt_map_get_begin_mapping_position(gt_map* const map) {
  return map->position; //-gt_map_get_left_trim_length(map); // FIXME
}
GT_INLINE uint64_t gt_map_get_end_mapping_position(gt_map* const map) {
  return map->position+gt_map_get_length(map) - 1; // FIXME No del counted
}
/*
 * Mismatch Handlers
 */
GT_INLINE void gt_map_add_misms(gt_map* const map,gt_misms* misms) {
  GT_MAP_CHECK(map);
  gt_vector_insert(map->mismatches,*misms,gt_misms);
}
GT_INLINE void gt_map_clear_misms(gt_map* const map) {
  GT_MAP_CHECK(map);
  gt_vector_clear(map->mismatches);
}
GT_INLINE gt_misms* gt_map_get_misms(gt_map* const map,const uint64_t offset) {
  GT_MAP_CHECK(map);
  return gt_vector_get_elm(map->mismatches,offset,gt_misms);
}
GT_INLINE void gt_map_set_misms(gt_map* const map,gt_misms* misms,const uint64_t offset) {
  GT_MAP_CHECK(map);
  gt_vector_set_elm(map->mismatches,offset,gt_misms,*misms);
}
GT_INLINE uint64_t gt_map_get_num_misms(gt_map* const map) {
  GT_MAP_CHECK(map);
  return gt_vector_get_used(map->mismatches);
}
GT_INLINE void gt_map_set_num_misms(gt_map* const map,const uint64_t num_misms) {
  GT_MAP_CHECK(map);
  gt_vector_set_used(map->mismatches,num_misms);
}
// Trim
GT_INLINE void gt_map_left_trim(gt_map* const map,const uint64_t length) {
  GT_MAP_CHECK(map);
  // Find trimmed mismatches
  uint64_t removed_misms = 0;
  GT_MISMS_ITERATE(map,misms) {
    if (misms->position>length) break;
    if (misms->misms_type==DEL && misms->position+misms->size>length) {
      misms->size = misms->position+misms->size-length;
      break;
    } else {
      ++removed_misms;
    }
  }
  // Remove trimmed mismatches
  if (removed_misms > 0) {
    gt_misms* const misms_vector = gt_vector_get_mem(map->mismatches,gt_misms);
    const uint64_t remaining_misms = gt_vector_get_used(map->mismatches)-removed_misms;
    uint64_t i;
    for (i=0;i<remaining_misms;++i) {
      misms_vector[i] = misms_vector[i+removed_misms];
      misms_vector[i].position -= length; // Adjust Position
    }
    gt_vector_set_used(map->mismatches,remaining_misms);
  } else {
    GT_MISMS_ITERATE(map,misms) {
      misms->position = (misms->position>length) ? misms->position-length : 0; // Adjust Position
    }
  }
  // Reduce base length
  map->base_length -= length;
  map->position += length;
}
GT_INLINE void gt_map_right_trim(gt_map* const map,const uint64_t length) {
  GT_MAP_CHECK(map);
  // Find trimmed mismatches
  const uint64_t misms_used = gt_vector_get_used(map->mismatches);
  gt_misms* misms_vector = gt_vector_get_mem(map->mismatches,gt_misms);
  const uint64_t new_eff_length = gt_map_get_base_length(map)-length;
  uint64_t misms_removed=0;
  int64_t i;
  for (i=misms_used-1;i>=0;--i) {
    if (misms_vector[i].position < new_eff_length) {
      if (misms_vector[i].misms_type==DEL && misms_vector[i].position+misms_vector[i].size > new_eff_length) {
        misms_vector[i].size = new_eff_length-misms_vector[i].position;
      }
      break;
    }
    ++misms_removed;
  }
  // Remove trimmed mismatches
  gt_vector_set_used(map->mismatches,misms_used-misms_removed);
  // Reduce base length
  map->base_length -= length;
}
GT_INLINE void gt_map_restore_left_trim(gt_map* const map,const uint64_t length) {
  GT_MAP_CHECK(map);
  // Check first misms
  bool inserted_trim = false;
  if (gt_vector_get_used(map->mismatches) > 0) {
    gt_misms* const first_misms = gt_map_get_misms(map,0);
    if (first_misms->position==0 && first_misms->misms_type==DEL) {
      // Add trim to the alredy existing deletion at the beginning
      first_misms->size += length;
      inserted_trim = true;
    }
  }
  if (!inserted_trim) {
    // Allocate one more misms
    gt_vector_reserve_additional(map->mismatches,1);
    // Shift all misms 1-right
    const uint64_t num_misms = gt_vector_get_used(map->mismatches);
    gt_misms* const misms_vector = gt_vector_get_mem(map->mismatches,gt_misms);
    int64_t i;
    for (i=num_misms-1;i>=0;--i) misms_vector[i+1] = misms_vector[i];
    gt_vector_add_used(map->mismatches,1);
    // Add the trim
    misms_vector[0].misms_type = DEL;
    misms_vector[0].position = 0;
    misms_vector[0].size = length;
  }
  // Adjust the rest of positions
  GT_VECTOR_ITERATE_OFFSET(map->mismatches,misms,misms_pos,1,gt_misms) {
    misms->position += length;
  }
  // Adjust base length
  map->base_length += length;
  // map->position -= length; // BUG
}
GT_INLINE void gt_map_restore_right_trim(gt_map* const map,const uint64_t length) {
  GT_MAP_CHECK(map);
  // Check last misms
  const uint64_t num_misms = gt_vector_get_used(map->mismatches);
  if (num_misms > 0) {
    gt_misms* const last_misms = gt_map_get_misms(map,num_misms-1);
    if (last_misms->position+last_misms->size==map->base_length && last_misms->misms_type==DEL) {
      last_misms->size += length;
      map->base_length += length;
      return;
    }
  }
  // Add the trim
  gt_misms misms;
  misms.misms_type = DEL;
  misms.position = map->base_length;
  misms.size = length;
  gt_map_add_misms(map,&misms);
  // Adjust base length
  map->base_length += length;
}
// Counting
GT_INLINE uint64_t gt_map_get_num_mismatch_bases(gt_map* const map) {
  GT_MAP_CHECK(map);
  uint64_t count = 0;
  GT_MISMS_ITERATE(map,misms_it) {
    switch (misms_it->misms_type) {
      case MISMS:
        ++count;
        break;
      case INS:
      case DEL:
        break;
    }
  }
  return count;
}
GT_INLINE uint64_t gt_map_get_num_indels(gt_map* const map) {
  GT_MAP_CHECK(map);
  uint64_t count = 0;
  GT_MISMS_ITERATE(map,misms_it) {
    switch (misms_it->misms_type) {
      case MISMS: break;
      case INS:
      case DEL:
        ++count;
        break;
    }
  }
  return count;
}
GT_INLINE uint64_t gt_map_get_num_insertions(gt_map* const map) {
  GT_MAP_CHECK(map);
  uint64_t count = 0;
  GT_MISMS_ITERATE(map,misms_it) {
    switch (misms_it->misms_type) {
      case MISMS: break;
      case INS:
        ++count;
        break;
      case DEL:   break;
    }
  }
  return count;
}
GT_INLINE uint64_t gt_map_get_num_deletions(gt_map* const map) {
  GT_MAP_CHECK(map);
  uint64_t count = 0;
  GT_MISMS_ITERATE(map,misms_it) {
    switch (misms_it->misms_type) {
      case MISMS: break;
      case INS:   break;
      case DEL:
        ++count;
        break;
    }
  }
  return count;
}
/*
 * Attributes (lazy allocation of the attributes field)
 */
GT_INLINE void* gt_map_attribute_get(gt_map* const map,char* const attribute_id) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(attribute_id);
  return (map->attributes==NULL) ? NULL : gt_attributes_get(map->attributes,attribute_id);
}
GT_INLINE void gt_map_attribute_set_(gt_map* const map,char* const attribute_id,void* const attribute,const size_t element_size) {
  GT_MAP_CHECK(map);
  if (map->attributes==NULL) {
    map->attributes = gt_attributes_new();
  }
  gt_attributes_add_primitive(map->attributes,attribute_id,attribute,element_size);
}
/*
 * Multiple BlockMap Handlers
 *   One MapBlock is a continuous mapping of a read (or subread) to
 *     a reference sequence (i.e. chromosome) on a specific position and strand.
 *   A Map is a list of MapBlock connected
 *   Map = (map.block1,map.block2,map.block3)
 *     MapBlock_1 := map.block1
 *     MapBlock_2 := map.block2
 *     MapBlock_3 := map.block3
 *   This concept is essential as to handle SPLIT-MAPS (i.e. RNA mappings)
 */
GT_INLINE uint64_t gt_map_get_num_blocks(gt_map* const map) {
  GT_MAP_CHECK(map);
  uint64_t num_blocks = 1;
  gt_map* aux_map = map;
  while (gt_map_has_next_block(aux_map)) {
    aux_map = gt_map_get_next_block(aux_map);
    ++num_blocks;
  }
  return num_blocks;
}
GT_INLINE bool gt_map_has_next_block(gt_map* const map) {
  GT_MAP_CHECK(map);
  return (map->next_block.map!=NULL);
}
GT_INLINE gt_map* gt_map_get_next_block(gt_map* const map) {
  GT_MAP_CHECK(map);
  return map->next_block.map;
}
GT_INLINE gt_map* gt_map_get_last_block(gt_map* const map) {
  GT_MAP_CHECK(map);
  gt_map* aux_map = map;
  while (gt_map_has_next_block(aux_map)) {
    aux_map = gt_map_get_next_block(aux_map);
  }
  return aux_map;
}
GT_INLINE void gt_map_set_next_block(
    gt_map* const map,gt_map* const next_map,const gt_junction_t junction,const int64_t junction_size) {
  GT_MAP_CHECK(map);
  map->next_block.junction = (next_map!=NULL) ? junction : NO_JUNCTION;
  map->next_block.junction_size = junction_size;
  map->next_block.map = next_map;
}
GT_INLINE void gt_map_insert_next_block(
    gt_map* const map,gt_map* const next_map,const gt_junction_t junction,const int64_t junction_size) {
  GT_MAP_CHECK(map);
  GT_MAP_CHECK(next_map);
  gt_map_junction aux_next_block = map->next_block;
  map->next_block.junction = junction;
  map->next_block.junction_size = junction_size;
  map->next_block.map = next_map;
  next_map->next_block = aux_next_block;
}
GT_INLINE void gt_map_append_block(
    gt_map* const map,gt_map* const next_map,const gt_junction_t junction,const int64_t junction_size) {
  GT_MAP_CHECK(map);
  GT_MAP_CHECK(next_map);
  gt_map* aux_map = gt_map_get_last_block(map);
  gt_map_set_next_block(aux_map,next_map,junction,junction_size);
}
// Junctions
GT_INLINE gt_junction_t gt_map_get_junction(gt_map* const map) {
  GT_MAP_CHECK(map);
  GT_MAP_NEXT_BLOCK_CHECK(map);
  return map->next_block.junction;
}
GT_INLINE int64_t gt_map_get_junction_size(gt_map* const map) {
  GT_MAP_CHECK(map);
  return (gt_map_has_next_block(map)) ? map->next_block.junction_size : 0;
}
// Reverse Blocks / Mismatches-Within-Block
GT_INLINE uint64_t gt_map_reverse_blocks_positions_(gt_map* const map,const uint64_t start_position) {
  GT_MAP_CHECK(map);
  if (gt_map_has_next_block(map)) {
    const uint64_t position =
        gt_map_reverse_blocks_positions_(gt_map_get_next_block(map),start_position) + gt_map_get_junction_size(map);
    gt_map_set_position(map,position);
    return position+gt_map_get_length(map);
  } else {
    gt_map_set_position(map,start_position);
    return start_position+gt_map_get_length(map);
  }
}
GT_INLINE void gt_map_reverse_blocks_positions(gt_map* const head_map,const uint64_t start_position) {
  GT_MAP_CHECK(head_map);
  gt_map_reverse_blocks_positions_(head_map,start_position);
}
#define GT_MAP_REVERSE_MISMS_ADJUST_POS(misms,base_length) \
  switch (misms->misms_type) { \
    case MISMS: misms->position = base_length - (misms->position + 1); misms->base = gt_get_complement(misms->base); break; \
    case DEL: misms->position = base_length - (misms->position + misms->size); break; \
    case INS: misms->position = base_length - misms->position; break; \
    default: break; \
  }
GT_INLINE void gt_map_reverse_misms(gt_map* const map) {
  GT_MAP_CHECK(map);
  // Flip all mismatches
  uint64_t z;
  const uint64_t base_length = gt_map_get_base_length(map);
  const uint64_t num_misms = gt_map_get_num_misms(map);
  const uint64_t mid_point = num_misms/2;
  for (z=0;z<mid_point;++z) {
    gt_misms* misms_a = gt_map_get_misms(map,z);
    gt_misms* misms_b = gt_map_get_misms(map,num_misms-1-z);
    // Correct position
    GT_MAP_REVERSE_MISMS_ADJUST_POS(misms_a,base_length);
    GT_MAP_REVERSE_MISMS_ADJUST_POS(misms_b,base_length);
    // Flip mismatches
    gt_misms misms = *misms_a;
    gt_map_set_misms(map,misms_b,z);
    gt_map_set_misms(map,&misms,num_misms-1-z);
  }
  if (num_misms%2) {
    gt_misms* misms = gt_map_get_misms(map,mid_point);
    GT_MAP_REVERSE_MISMS_ADJUST_POS(misms,base_length);
  }
}
/*
 * Map Segment Handlers
 *  One map segment is said to be a list of connected maps within the same
 *  reference sequence (i.e. chromosome) mapping the same strand.
 *     Map = (map.block1,map.block2,map.block3)
 *       s.t. map.block1.seq_name == map.block2.seq_name && map.block1.strand == map.block2.strand
 *            map.block2.seq_name != map.block3.seq_name && map.block2.strand != map.block3.strand
 *     MapSegment_1 := {map.block1,map.block2}
 *     MapSegment_2 := {map.block3}
 *   This concept is essential as to handle properly QUIMERAS
 */
GT_INLINE uint64_t gt_map_segment_get_num_segments(gt_map* const map) {
  GT_MAP_CHECK(map);
  gt_map* current_segment = map;
  uint64_t num_segments = 1;
  while ((current_segment=gt_map_segment_get_next_segment(current_segment))!=NULL) {
    ++num_segments;
  }
  return num_segments;
}
GT_INLINE gt_map* gt_map_segment_get_next_block(gt_map* const map) {
  GT_MAP_CHECK(map);
  gt_map* const next_map = gt_map_get_next_block(map);
  return (next_map==NULL || !GT_MAP_IS_SAME_SEGMENT(map,next_map)) ? NULL : next_map;
}
GT_INLINE gt_map* gt_map_segment_get_last_block(gt_map* const map) {
  GT_MAP_CHECK(map);
  gt_map* current_map = map;
  do {
    gt_map* const next_map = gt_map_get_next_block(current_map);
    if (next_map==NULL || !GT_MAP_IS_SAME_SEGMENT(current_map,next_map)) break;
    current_map = next_map;
  } while (true);
  return current_map;
}
GT_INLINE gt_map* gt_map_segment_get_next_segment(gt_map* const map) {
  GT_MAP_CHECK(map);
  return gt_map_get_next_block(gt_map_segment_get_last_block(map));
}
// Most right-hand position of the first/last segment if forward/reverse strand
GT_INLINE uint64_t gt_map_get_global_coordinate(gt_map* const map) {
  GT_MAP_CHECK(map);
  return (map->strand==FORWARD) ? gt_map_get_position(map): gt_map_get_position(gt_map_segment_get_last_block(map));
}

/*
 * Map's Mismatches iterator
 *   Iterate over the mismatches(M/I/D) of a map
 */
GT_INLINE void gt_map_new_misms_iterator(gt_map* const map,gt_map_mism_iterator* const map_mism_iterator) {
  GT_NULL_CHECK(map_mism_iterator);
  GT_MAP_CHECK(map);
  map_mism_iterator->map = map;
  map_mism_iterator->next_pos = 0;
  map_mism_iterator->total_pos = gt_vector_get_used(map->mismatches);
  map_mism_iterator->next_misms = gt_vector_get_mem(map->mismatches,gt_misms);
}
GT_INLINE gt_misms* gt_map_next_misms(gt_map_mism_iterator* const map_mism_iterator) {
  GT_NULL_CHECK(map_mism_iterator);
  GT_MAP_CHECK(map_mism_iterator->map);
  gt_map* const map = map_mism_iterator->map;
  if (gt_expect_true(map_mism_iterator->next_pos<gt_vector_get_used(map->mismatches))) {
    gt_misms* const misms = gt_vector_get_elm(map->mismatches,map_mism_iterator->next_pos,gt_misms);
    ++map_mism_iterator->next_pos;
    return misms;
  } else {
    return NULL;
  }
}
GT_INLINE uint64_t gt_map_next_misms_pos(gt_map_mism_iterator* const map_mism_iterator) {
  GT_NULL_CHECK(map_mism_iterator);
  GT_MAP_CHECK(map_mism_iterator->map);
  return map_mism_iterator->next_pos;
}
/*
 * Map's Blocks iterator
 *   Iterate over the map blocks of a map
 *     Map = (map.block1,map.block2)
 *     GT_MAP_BLOCKS_ITERATE := {map.block1,map.block2}
 */
GT_INLINE void gt_map_new_block_iterator(gt_map* const map,gt_map_block_iterator* const map_block_iterator) {
  GT_NULL_CHECK(map_block_iterator);
  GT_MAP_CHECK(map);
  map_block_iterator->map = map;
  map_block_iterator->next_map = map;
}
GT_INLINE gt_map* gt_map_next_block(gt_map_block_iterator* const map_block_iterator) {
  GT_NULL_CHECK(map_block_iterator);
  gt_map* returned_map = map_block_iterator->next_map;
  map_block_iterator->next_map = (returned_map!=NULL) ? returned_map->next_block.map : NULL;
  return returned_map;
}
/*
 * Map's Segments iterator
 *   Iterate over the connected segments of a map
 *     Map = (map.block1,map.block2,map.block3)
 *       s.t. map.block1.seq_name == map.block2.seq_name && map.block1.strand == map.block2.strand
 *            map.block2.seq_name != map.block3.seq_name && map.block2.strand != map.block3.strand
 *     GT_MAP_BLOCKS_ITERATE := {map.block1,map.block3}
 */
GT_INLINE void gt_map_new_segment_iterator(gt_map* const map,gt_map_segment_iterator* const map_segment_iterator) {
  GT_NULL_CHECK(map_segment_iterator);
  GT_MAP_CHECK(map);
  map_segment_iterator->next_map = map;
  map_segment_iterator->accumulated_offset = 0;
  map_segment_iterator->total_base_length = 0;
  map_segment_iterator->begin_map_segment = NULL;
  map_segment_iterator->end_map_segment = NULL;
  map_segment_iterator->total_map_blocks = 0;
}
GT_INLINE bool gt_map_next_segment(gt_map_segment_iterator* const map_segment_iterator) {
  GT_NULL_CHECK(map_segment_iterator);
  if (map_segment_iterator->next_map == NULL) return false;
  // Accumulate offset
  map_segment_iterator->accumulated_offset += map_segment_iterator->total_base_length;
  // Init
  gt_map* const base_map = map_segment_iterator->next_map;
  map_segment_iterator->begin_map_segment = base_map;
  map_segment_iterator->end_map_segment = base_map;
  map_segment_iterator->total_base_length = gt_map_get_base_length(base_map);
  map_segment_iterator->total_map_blocks = 1;
  // Traverse the connected segment
  gt_map* map = base_map;
  gt_map* next_map = NULL;
  do {
    next_map = gt_map_get_next_block(map);
    // Check not connected or end
    if (next_map==NULL || !GT_MAP_IS_SAME_SEGMENT(base_map,next_map)) break;
    // Add Block
    map_segment_iterator->total_base_length += gt_map_get_base_length(next_map);
    ++(map_segment_iterator->total_map_blocks);
    map_segment_iterator->end_map_segment = next_map;
    // Next Block
    map = next_map;
  } while (true);
  // Set next & return
  map_segment_iterator->next_map = next_map;
  return true;
}
GT_INLINE uint64_t gt_map_segment_iterator_get_accumulated_offset(gt_map_segment_iterator* const map_segment_iterator) {
  GT_MAP_SEGMENT_ITERATOR_CHECK(map_segment_iterator);
  return map_segment_iterator->accumulated_offset;
}
GT_INLINE uint64_t gt_map_segment_iterator_get_remaining_bases(gt_map_segment_iterator* const map_segment_iterator,gt_string* const read) {
  GT_MAP_SEGMENT_ITERATOR_CHECK(map_segment_iterator);
  if (read!=NULL) {
    const int64_t remaining_bases = (int64_t)gt_string_get_length(read)-
        (int64_t)(map_segment_iterator->accumulated_offset+map_segment_iterator->total_base_length);
    return (remaining_bases>0) ? remaining_bases : 0;
  } else {
    return 0;
  }
}
GT_INLINE uint64_t gt_map_segment_iterator_get_position(gt_map_segment_iterator* const map_segment_iterator) {
  return (map_segment_iterator->begin_map_segment->strand==FORWARD) ?
      gt_map_get_position(map_segment_iterator->begin_map_segment) : gt_map_get_position(map_segment_iterator->end_map_segment);
}
GT_INLINE uint64_t gt_map_segment_iterator_get_total_map_blocks(gt_map_segment_iterator* const map_segment_iterator) {
  GT_MAP_SEGMENT_ITERATOR_CHECK(map_segment_iterator);
  return map_segment_iterator->total_map_blocks;
}
GT_INLINE uint64_t gt_map_segment_iterator_get_total_base_length(gt_map_segment_iterator* const map_segment_iterator) {
  GT_MAP_SEGMENT_ITERATOR_CHECK(map_segment_iterator);
  return map_segment_iterator->total_base_length;
}
GT_INLINE gt_map* gt_map_segment_iterator_get_map(gt_map_segment_iterator* const map_segment_iterator) {
  GT_MAP_SEGMENT_ITERATOR_CHECK(map_segment_iterator);
  return map_segment_iterator->begin_map_segment;
}
/*
 * Miscellaneous
 */
GT_INLINE gt_map* gt_map_copy(gt_map* const map) {
  GT_MAP_CHECK(map);
  gt_map* map_cpy = gt_map_new();
  gt_string_copy(map_cpy->seq_name,map->seq_name);
  map_cpy->position = map->position;
  map_cpy->base_length = map->base_length;
  map_cpy->strand = map->strand;
  map_cpy->phred_score = map->phred_score;
  map_cpy->gt_score = map->gt_score;
  gt_vector_copy(map_cpy->mismatches,map->mismatches);
  // Copy next blocks
  map_cpy->next_block = map->next_block;
  map_cpy->next_block.junction = map->next_block.junction;
  map_cpy->next_block.junction_size = map->next_block.junction_size;
  if (map->next_block.map!=NULL) map_cpy->next_block.map = gt_map_copy(map->next_block.map);
  return map_cpy;
}
GT_INLINE gt_map** gt_mmap_array_copy(gt_map** mmap,const uint64_t num_blocks) {
  GT_ZERO_CHECK(num_blocks);
  gt_map** mmap_copy = gt_calloc(num_blocks,gt_map*,false);
  uint64_t i;
  for (i=0;i<num_blocks;++i) {
    if (mmap[i] != NULL) {
      mmap_copy[i] = gt_map_copy(mmap[i]);
    } else {
      mmap_copy[i] = NULL;
    }
  }
  return mmap_copy;
}
