/*
 * PROJECT: GEM-Tools library
 * FILE: gt_alignment.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_ALIGNMENT_H_
#define GT_ALIGNMENT_H_

#include "gt_essentials.h"
#include "gt_attributes.h"

#include "gt_map.h"

#include "gt_input_parser.h"

// Alignment itself
typedef struct _gt_alignment_dictionary gt_alignment_dictionary; // Forward declaration of gt_alignment_dictionary
typedef struct {
  /* IDs */
  uint32_t alignment_id;
  uint32_t in_block_id;
  /* Basic components */
  gt_string* tag;
  gt_string* read;
  gt_string* qualities;
  gt_vector* counters;
  /* Maps structures */
  gt_vector* maps; /* (gt_map*) */
  /* Attibutes */
  gt_attributes* attributes;
  /* Hashed Dictionary */
  gt_alignment_dictionary* alg_dictionary;
} gt_alignment;

// Iterator
typedef struct {
  gt_alignment* alignment;
  uint64_t next_pos;
} gt_alignment_map_iterator;

// Map Dictionary (For Fast Indexing)
typedef struct {
  gt_ihash* begin_position; /* (uint64_t) */
  gt_ihash* end_position;   /* (uint64_t) */
} gt_alignment_dictionary_element;
struct _gt_alignment_dictionary {
  gt_shash* maps_dictionary; /* (gt_alignment_dictionary_element*) */
  gt_alignment* alignment;
};

/*
 * Checkers
 */
#define GT_ALIGNMENT_CHECK(alignment) \
  GT_NULL_CHECK(alignment); \
  GT_STRING_CHECK(alignment->tag); \
  GT_STRING_CHECK(alignment->read); \
  GT_STRING_CHECK(alignment->qualities); \
  GT_VECTOR_CHECK(alignment->counters); \
  GT_VECTOR_CHECK(alignment->maps); \
  GT_ATTRIBUTES_CHECK(alignment->attributes)
#define GT_ALIGNMENT_DICTIONARY_CHECK(alignment_dictionary) \
  GT_NULL_CHECK(alignment_dictionary); \
  GT_HASH_CHECK(alignment_dictionary->maps_dictionary)

/*
 * Setup
 */
GT_INLINE gt_alignment* gt_alignment_new(void);
GT_INLINE void gt_alignment_clear(gt_alignment* const alignment);
GT_INLINE void gt_alignment_delete(gt_alignment* const alignment);

/*
 * Accessors
 */
GT_INLINE char* gt_alignment_get_tag(gt_alignment* const alignment);
GT_INLINE gt_string* gt_alignment_get_string_tag(gt_alignment* const alignment);
GT_INLINE void gt_alignment_set_tag(gt_alignment* const alignment,char* const tag,const uint64_t length);
GT_INLINE uint64_t gt_alignment_get_tag_length(gt_alignment* const alignment);

GT_INLINE char* gt_alignment_get_read(gt_alignment* const alignment);
GT_INLINE void gt_alignment_set_read(gt_alignment* const alignment,char* const read,const uint64_t length);
GT_INLINE uint64_t gt_alignment_get_read_length(gt_alignment* const alignment);

GT_INLINE char* gt_alignment_get_qualities(gt_alignment* const alignment);
GT_INLINE void gt_alignment_set_qualities(gt_alignment* const alignment,char* const qualities,const uint64_t length);
GT_INLINE bool gt_alignment_has_qualities(gt_alignment* const alignment);

GT_INLINE gt_vector* gt_alignment_get_counters_vector(gt_alignment* const alignment);
GT_INLINE void gt_alignment_set_counters_vector(gt_alignment* const alignment,gt_vector* const counters);
GT_INLINE uint64_t gt_alignment_get_num_counters(gt_alignment* const alignment);
GT_INLINE uint64_t gt_alignment_get_counter(gt_alignment* const alignment,const uint64_t stratum);
GT_INLINE void gt_alignment_set_counter(gt_alignment* const alignment,const uint64_t stratum,const uint64_t value);
GT_INLINE void gt_alignment_dec_counter(gt_alignment* const alignment,const uint64_t stratum);
GT_INLINE void gt_alignment_inc_counter(gt_alignment* const alignment,const uint64_t stratum);

/*
 * Attribute accessors
 */
GT_INLINE uint64_t gt_alignment_get_mcs(gt_alignment* const alignment);
GT_INLINE void gt_alignment_set_mcs(gt_alignment* const alignment,uint64_t max_complete_strata);
GT_INLINE void gt_alignment_set_not_unique_flag(gt_alignment* const alignment,bool is_not_unique);
GT_INLINE bool gt_alignment_get_not_unique_flag(gt_alignment* const alignment);
GT_INLINE int64_t gt_alignment_get_pair(gt_alignment* const alignment);
GT_INLINE void gt_alignment_set_map_primary(gt_alignment* const alignment,gt_map* const map);
GT_INLINE gt_map* gt_alignment_get_map_primary(gt_alignment* const alignment);

/*
 * Maps Handlers
 */
GT_INLINE uint64_t gt_alignment_get_num_maps(gt_alignment* const alignment);
GT_INLINE void gt_alignment_add_map(gt_alignment* const alignment,gt_map* const map);
GT_INLINE void gt_alignment_add_map_gt_vector(gt_alignment* const alignment,gt_vector* const map_vector);
GT_INLINE gt_map* gt_alignment_get_map(gt_alignment* const alignment,const uint64_t position);
GT_INLINE void gt_alignment_set_map(gt_alignment* const alignment,gt_map* const map,const uint64_t position);
GT_INLINE void gt_alignment_clear_maps(gt_alignment* const alignment);

GT_INLINE bool gt_alignment_locate_map_reference(gt_alignment* const alignment,gt_map* const map,uint64_t* const position);

/*
 * Miscellaneous
 */
GT_INLINE gt_alignment* gt_alignment_copy(gt_alignment* const alignment,const bool copy_maps);

// Alignment's Maps iterator
GT_INLINE void gt_alignment_new_map_iterator(gt_alignment* const alignment,gt_alignment_map_iterator* const alignment_map_iterator);
GT_INLINE gt_map* gt_alignment_next_map(gt_alignment_map_iterator* const alignment_map_iterator);
GT_INLINE uint64_t gt_alignment_next_map_pos(gt_alignment_map_iterator* const alignment_map_iterator);

/*
 * Iterate over the map of an alignment
 *   Alignment = {(map1),(map2.block1,map2.block2),(map3),(map4)}
 *   GT_ALIGNMENT_ITERATE := {(map1),(map2.block1,map2.block2),(map3),(map4)}
 */
#define GT_ALIGNMENT_ITERATE(alignment,map) \
  /* Alignment. Iterate over all maps */ \
  gt_alignment_map_iterator __##map##_iterator; \
  gt_map* map; \
  gt_alignment_new_map_iterator(alignment,&(__##map##_iterator)); \
  while ((map=gt_alignment_next_map(&(__##map##_iterator))))

/*
 * Map Dictionary (For Fast Indexing)
 */
GT_INLINE gt_alignment_dictionary* gt_alignment_dictionary_new(gt_alignment* const alignment);
GT_INLINE void gt_alignment_dictionary_delete(gt_alignment_dictionary* const alignment_dictionary);
GT_INLINE bool gt_alignment_dictionary_try_add(
    gt_alignment_dictionary* const alignment_dictionary,gt_map* const map,
    const uint64_t begin_position,const uint64_t end_position,
    gt_alignment_dictionary_element** alg_dicc_elem,
    gt_ihash_element** ihash_element_b,gt_ihash_element** ihash_element_e,uint64_t const vector_position);
GT_INLINE void gt_alignment_dictionary_record_position(
    gt_alignment_dictionary* const alignment_dictionary,
    const uint64_t begin_position,const uint64_t end_position,
    gt_alignment_dictionary_element* const alg_dicc_elem,
    gt_ihash_element* const ihash_element_b,gt_ihash_element* const ihash_element_e,const uint64_t vector_position);

#endif /* GT_ALIGNMENT_H_ */
