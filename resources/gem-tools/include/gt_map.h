/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_MAP_H_
#define GT_MAP_H_

#include "gt_essentials.h"
#include "gt_attributes.h"

#include "gt_misms.h"
#include "gt_dna_string.h"

/*
 * Constants
 */
#define GT_MAP_NO_GT_SCORE UINT64_MAX
#define GT_MAP_NO_PHRED_SCORE 255

/*
 * Junction (gt_map_junction)
 */
// Types of junctions between map blocks
typedef enum { NO_JUNCTION, SPLICE, POSITIVE_SKIP, NEGATIVE_SKIP, QUIMERA } gt_junction_t;
typedef struct _gt_map gt_map; // Forward declaration of gt_map
typedef struct {
  gt_map* map;
  gt_junction_t junction;
  int64_t junction_size;
} gt_map_junction;

/*
 * Map (gt_map)
 */
struct _gt_map {
  /* Sequence-name(Chromosome/Contig/...), position and strand */
  gt_string* seq_name;
  uint64_t position;
  uint64_t base_length; // Length not including indels
  gt_strand strand;
  /* Metrics */
  uint64_t gt_score;
  uint8_t phred_score;
  /* Mismatches/Indels */
  gt_vector* mismatches; /* (misms_t) */
  /* Multiple Block Map (splice-map, local alignments, ...) */
  gt_map_junction next_block;
  /* Attributes */
  gt_attributes* attributes;
};

// Iterators
typedef struct {
  gt_map* map;
  uint64_t next_pos;
  uint64_t total_pos;
  gt_misms* next_misms;
} gt_map_mism_iterator;
typedef struct {
  gt_map* map;
  gt_map* next_map;
} gt_map_block_iterator;
typedef struct {
  uint64_t accumulated_offset;
  gt_map* begin_map_segment;
  gt_map* end_map_segment;
  uint64_t total_map_blocks;
  uint64_t total_base_length;
  gt_map* next_map;
} gt_map_segment_iterator;

/*
 * Checkers
 */
#define GT_MAP_CHECK(map) \
  GT_NULL_CHECK(map); \
  GT_NULL_CHECK((map)->mismatches)
#define GT_MAP_NEXT_BLOCK_CHECK(map) \
  GT_MAP_CHECK(map->next_block.map)
#define GT_MAP_SEGMENT_ITERATOR_CHECK(map_segment_iterator) \
  GT_NULL_CHECK(map_segment_iterator); \
  GT_MAP_CHECK((map_segment_iterator)->begin_map_segment); \
  GT_ZERO_CHECK((map_segment_iterator)->total_base_length); \
  GT_ZERO_CHECK((map_segment_iterator)->total_map_blocks)

/*
 * Setup
 */
GT_INLINE gt_map* gt_map_new(void);
GT_INLINE void gt_map_clear(gt_map* const map);
GT_INLINE void gt_map_delete(gt_map* const map);

/*
 * MapBlock Accessors
 */
GT_INLINE char* gt_map_get_seq_name(gt_map* const map);
GT_INLINE uint64_t gt_map_get_seq_name_length(gt_map* const map);
GT_INLINE gt_string* gt_map_get_string_seq_name(gt_map* const map);
GT_INLINE void gt_map_set_seq_name(gt_map* const map,const char* const seq_name,const uint64_t length);
GT_INLINE void gt_map_set_string_seq_name(gt_map* const map,gt_string* const seq_name);
GT_INLINE gt_strand gt_map_get_strand(gt_map* const map);
GT_INLINE void gt_map_set_strand(gt_map* const map,const gt_strand strand);
// Length of the base read (no indels)
GT_INLINE uint64_t gt_map_get_base_length(gt_map* const map);
GT_INLINE void gt_map_set_base_length(gt_map* const map,const uint64_t length);
// Most right-hand position of the current map block (NOT GLOBAL COORDINATE --> gt_map_get_global_coordinate())
GT_INLINE uint64_t gt_map_get_position(gt_map* const map);
GT_INLINE void gt_map_set_position(gt_map* const map,const uint64_t position);
// Trim helpers
GT_INLINE uint64_t gt_map_get_left_trim_length(gt_map* const map);
GT_INLINE uint64_t gt_map_get_right_trim_length(gt_map* const map);
// Begin/End Mapping Position
GT_INLINE uint64_t gt_map_get_begin_mapping_position(gt_map* const map);
GT_INLINE uint64_t gt_map_get_end_mapping_position(gt_map* const map);

/*
 * Mismatch Handlers
 */
GT_INLINE void gt_map_add_misms(gt_map* const map,gt_misms* misms);
GT_INLINE void gt_map_clear_misms(gt_map* const map);
GT_INLINE gt_misms* gt_map_get_misms(gt_map* const map,const uint64_t offset);
GT_INLINE void gt_map_set_misms(gt_map* const map,gt_misms* misms,const uint64_t offset);
GT_INLINE uint64_t gt_map_get_num_misms(gt_map* const map);
GT_INLINE void gt_map_set_num_misms(gt_map* const map,const uint64_t num_misms);
// Trim
GT_INLINE void gt_map_left_trim(gt_map* const map,const uint64_t length);
GT_INLINE void gt_map_right_trim(gt_map* const map,const uint64_t length);
GT_INLINE void gt_map_restore_left_trim(gt_map* const map,const uint64_t length);
GT_INLINE void gt_map_restore_right_trim(gt_map* const map,const uint64_t length);
// Counting
GT_INLINE uint64_t gt_map_get_num_mismatch_bases(gt_map* const map);
GT_INLINE uint64_t gt_map_get_num_indels(gt_map* const map);
GT_INLINE uint64_t gt_map_get_num_insertions(gt_map* const map);
GT_INLINE uint64_t gt_map_get_num_deletions(gt_map* const map);

/*
 * Attributes (lazy allocation of the attributes field)
 */
GT_INLINE void* gt_map_attribute_get(gt_map* const map,char* const attribute_id);
GT_INLINE void gt_map_attribute_set_(gt_map* const map,char* const attribute_id,void* const attribute,const size_t element_size);
#define gt_map_attribute_set(map,attribute_id,attribute,element_type) \
    gt_map_attribute_set_(map,attribute_id,(void*)attribute,sizeof(element_type))

// Utils
GT_INLINE uint64_t gt_map_get_min_intron_length(gt_map* const map);
GT_INLINE uint64_t gt_map_get_min_block_length(gt_map* const map);

/*
 * Multiple BlockMap Handlers
 *   One MapBlock is a continuous mapping of a read (or subread) to
 *     a reference sequence (i.e. chromosome) on a specific position and strand.
 *   A Map is a list of MapBlock connected
 *   Map = (map.block1,map.block2,map.block3)
 *     MapBlock_1 := map.block1
 *     MapBlock_2 := map.block2
 *     MapBlock_3 := map.block3
 *   This concept is essential as to handle SPLIT-MAPS (i.e. RNA mappings) and QUIMERAS
 */
GT_INLINE uint64_t gt_map_get_num_blocks(gt_map* const map);
GT_INLINE bool gt_map_has_next_block(gt_map* const map);
GT_INLINE gt_map* gt_map_get_next_block(gt_map* const map);
GT_INLINE gt_map* gt_map_get_last_block(gt_map* const map);
GT_INLINE void gt_map_set_next_block(
    gt_map* const map,gt_map* const next_map,const gt_junction_t junction,const int64_t junction_size);
GT_INLINE void gt_map_insert_next_block(
    gt_map* const map,gt_map* const next_map,const gt_junction_t junction,const int64_t junction_size);
GT_INLINE void gt_map_append_block(
    gt_map* const map,gt_map* const next_map,const gt_junction_t junction,const int64_t junction_size);
// Junctions
GT_INLINE gt_junction_t gt_map_get_junction(gt_map* const map);
GT_INLINE int64_t gt_map_get_junction_size(gt_map* const map);
// Reverse Blocks / Mismatches-Within-Block
GT_INLINE void gt_map_reverse_blocks_positions(gt_map* const head_map,const uint64_t start_position);
GT_INLINE void gt_map_reverse_misms(gt_map* const map);

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
#define GT_MAP_IS_SAME_SEGMENT(map_1,map_2) \
  (gt_string_equals(gt_map_get_string_seq_name(map_1),gt_map_get_string_seq_name(map_2)) && \
   gt_map_get_strand(map_1)==gt_map_get_strand(map_2))
GT_INLINE uint64_t gt_map_segment_get_num_segments(gt_map* const map);
GT_INLINE gt_map* gt_map_segment_get_next_block(gt_map* const map);
GT_INLINE gt_map* gt_map_segment_get_last_block(gt_map* const map);
GT_INLINE gt_map* gt_map_segment_get_next_segment(gt_map* const map);
// Most right-hand position of the first/last segment if forward/reverse strand
//   NOTE: In case of doubt, use this one to get the coordinate of a map (gives MAP,SAM,... complaint position)
GT_INLINE uint64_t gt_map_get_global_coordinate(gt_map* const map);

/*
 * Map's Mismatches iterator
 *   Iterate over the mismatches(M/I/D) of a map
 *   {
 *     gt_map_mism_iterator misms_iterator;
 *     gt_misms* misms;
 *     gt_map_new_misms_iterator(map_block,&misms_iterator);
 *     while ((misms=gt_map_next_misms(&misms_iterator))) {
 *       // Map Block. Iterate over all mismatches
 *       ..code(misms)..
 *     }
 *   }
 */
GT_INLINE void gt_map_new_misms_iterator(gt_map* const map,gt_map_mism_iterator* const map_mism_iterator);
GT_INLINE gt_misms* gt_map_next_misms(gt_map_mism_iterator* const map_mism_iterator);
GT_INLINE uint64_t gt_map_next_misms_pos(gt_map_mism_iterator* const map_mism_iterator);
#define GT_MISMS_ITERATE(map_block,misms) \
  GT_VECTOR_ITERATE(map_block->mismatches,misms,misms_##vector_position,gt_misms)

/*
 *  Iterate over the map blocks of a map
 *     Map = (map.block1,map.block2)
 *     GT_MAP_BLOCKS_ITERATE := {map.block1,map.block2}
 *  i.e.
 *    gt_map_block_iterator map_block_iterator;
 *    gt_map* map_block;
 *    gt_map_new_block_iterator(map,&(map_block_iterator));
 *    while ((map_block=gt_map_next_block(&(map_block_iterator)))) {
 *
 *    }
 *  Otherwise, use the MACRO Iterator
 *    GT_MAP_ITERATE(map,map_block) {
 *      ..code(map_block)..
 *    }
 */
GT_INLINE void gt_map_new_block_iterator(gt_map* const map,gt_map_block_iterator* const map_block_iterator);
GT_INLINE gt_map* gt_map_next_block(gt_map_block_iterator* const map_block_iterator);
#define GT_MAP_ITERATE(map,map_block) \
  gt_map* map_block; \
  for (map_block=map;map_block!=NULL;map_block=gt_map_get_next_block(map_block))
#define GT_MAP_ITERATE_MAP_BLOCK(map,map_block) GT_MAP_ITERATE(map,map_block) /* Just an alias */

/*
 * Map's Segments iterator
 *   Iterate over the connected segments of a map
 *     Map = (map.block1,map.block2,map.block3)
 *       s.t. map.block1.seq_name == map.block2.seq_name && map.block1.strand == map.block2.strand
 *            map.block2.seq_name != map.block3.seq_name && map.block2.strand != map.block3.strand
 *     GT_MAP_BLOCKS_ITERATE := {map.block1,map.block3}
 *   Example of usage of the iterator
 *     GT_MAP_SEGMENT_ITERATOR(map,map_segment_iterator) {
 *       ..code(map_segment_iterator.begin_map_segment,map_segment_iterator.end_map_segment)..
 *     }
 *   Example of simple iteration over the segments of a map
 *     GT_MAP_ITERATE_SEGMENT(map,map_segment) {
 *       ..code(map_segment)..
 *     }
 */
GT_INLINE void gt_map_new_segment_iterator(gt_map* const map,gt_map_segment_iterator* const map_segment_iterator);
GT_INLINE bool gt_map_next_segment(gt_map_segment_iterator* const map_segment_iterator);
GT_INLINE uint64_t gt_map_segment_iterator_get_accumulated_offset(gt_map_segment_iterator* const map_segment_iterator);
GT_INLINE uint64_t gt_map_segment_iterator_get_remaining_bases(gt_map_segment_iterator* const map_segment_iterator,gt_string* const read);
GT_INLINE uint64_t gt_map_segment_iterator_get_position(gt_map_segment_iterator* const map_segment_iterator);
GT_INLINE uint64_t gt_map_segment_iterator_get_total_map_blocks(gt_map_segment_iterator* const map_segment_iterator);
GT_INLINE uint64_t gt_map_segment_iterator_get_total_base_length(gt_map_segment_iterator* const map_segment_iterator);
GT_INLINE gt_map* gt_map_segment_iterator_get_map(gt_map_segment_iterator* const map_segment_iterator);
#define GT_MAP_SEGMENT_ITERATOR(map,map_segment_iterator) \
  gt_map_segment_iterator map_segment_iterator; \
  gt_map_new_segment_iterator(map,&map_segment_iterator); \
  while (gt_map_next_segment(&map_segment_iterator))
#define GT_MAP_ITERATE_SEGMENT(map,map_segment) \
  gt_map* map_segment; \
  for (map_segment=map;map_segment!=NULL;map_segment=gt_map_segment_get_next_segment(map_segment))

/*
 * Iterator over the MapBlocks of a Map Segment
 * Example-1 (iterating over the MapBlocks of a segment)
 *   GT_MAP_SEGMENT_ITERATE(map_segment,map_block) {
 *     ..code(map_block)..
 *   }
 * Example-2 (Nested iteration.
 *   First over the segments of a Map.
 *   Second over the MapBlocks of the segment)
 *   GT_MAP_ITERATE_SEGMENT(map,map_segment) {
 *     GT_MAP_SEGMENT_ITERATE(map_segment,map_block) {
 *       ..code(map_block)..
 *     }
 *   }
 */
#define GT_MAP_SEGMENT_ITERATE(map_segment,map_block) \
  gt_map* map_block; \
  for (map_block=map_segment;map_block!=NULL;map_block=gt_map_segment_get_next_block(map_block))

/*
 * Miscellaneous
 */
GT_INLINE gt_map* gt_map_copy(gt_map* map);
GT_INLINE gt_map** gt_mmap_array_copy(gt_map** mmap,const uint64_t num_blocks);

/*
 * Map check/recover operators (handy when traversing mismatches)
 */
#define GT_MAP_CHECK__RELOAD_MISMS_PTR(map,misms_offset,misms_ptr,num_misms) { \
  if (misms_offset<num_misms) { \
    misms_ptr=gt_map_get_misms(map,misms_offset); \
  } else { \
    misms_ptr=NULL; \
  } \
}
#define GT_MAP_CHECK__RELOAD_MISMS(map,misms_offset,misms,num_misms) { \
  if (misms_offset<num_misms) { \
    misms=*gt_map_get_misms(map,misms_offset); \
  } else { \
    misms.position=UINT64_MAX; \
  } \
}

/*
 * Error Messages
 */
#define GT_ERROR_MAP_MISMS_NOT_PARSED "Map's mismatches not parsed yet"
#define GT_ERROR_MAP_NEG_LENGTH "Negative Map total length"
#define GT_ERROR_MAP_NEG_MAPPED_BASES "Negative number of bases mapped"

/*
 * Map Metrics
 */
#include "gt_map_metrics.h"

/*
 * Map Utils
 */
#include "gt_map_utils.h"

#endif /* GT_MAP_H_ */
