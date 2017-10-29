/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Data structure modeling sequences' templates.
 *   That is, set of alignments and relationships between their maps.
 */

#ifndef GT_TEMPLATE_H_
#define GT_TEMPLATE_H_

#include "gt_essentials.h"
#include "gt_alignment.h"
#include "gt_input_parser.h"

// Codes gt_status
#define GT_TEMPLATE_OK 1
#define GT_TEMPLATE_FAIL 0

// General types of alignment/template
typedef enum { GT_PE_MAPPED_PAIRED, GT_PE_MAPPED_UNPAIRED, GT_PE_UNMAPPED,
               GT_SE_MAPPED, GT_SE_UNMAPPED } gt_template_t;

typedef struct _gt_template_dictionary gt_template_dictionary; // Forward declaration of gt_alignment_dictionary

typedef struct {
  uint32_t template_id;
  uint32_t in_block_id;
  gt_string* tag;
  gt_alignment* alignment_end1;
  gt_alignment* alignment_end2;
  gt_vector* counters; /* (uint64_t) */
  gt_vector* mmaps; /* (gt_mmap) */
  gt_attributes* attributes;
  /* Hashed Dictionary */
  gt_template_dictionary* alg_dictionary;
} gt_template;
typedef struct {
  uint64_t distance;
  uint64_t gt_score;
  uint8_t phred_score;
  // gt_attributes* attributes; // Not needed for now
} gt_mmap_attributes;
typedef struct {
  gt_map* mmap[2];
  gt_mmap_attributes attributes;
} gt_mmap;

// Map Dictionary (For Fast Indexing)
typedef struct {
  gt_map** mmap; /* the map*/
  gt_mmap_attributes* mmap_attrs; /* the map attributes*/
  uint64_t pos; /* the map position in the alignment map list */
} gt_template_dictionary_map_element;
struct _gt_template_dictionary {
  gt_shash* refs_dictionary; /* (gt_template_dictionary_map_element*) */
  gt_template* template;
};


// Iterators
typedef struct {
  gt_template* template;
  uint64_t next_alignment_end;
} gt_template_alignment_iterator;
typedef struct {
  gt_template* template;
  uint64_t next_mmap_position;
} gt_template_maps_iterator;

/*
 * Checkers & Errors
 */
#define GT_TEMPLATE_CHECK(template) \
  GT_NULL_CHECK(template); \
  GT_STRING_CHECK(template->tag); \
  GT_VECTOR_CHECK(template->counters); \
  GT_VECTOR_CHECK(template->mmaps); \
  GT_ATTRIBUTES_CHECK(template->attributes)
#define GT_TEMPLATE_COMMON_CONSISTENCY_ERROR(template_A,template_B) \
  gt_cond_fatal_error(gt_template_get_num_blocks(template_A)!=gt_template_get_num_blocks(template_B),TEMPLATE_INCONSISTENT_NUM_BLOCKS)
#define GT_MMAP_ARRAY_CHECK(mmap) \
  gt_cond_fatal_error(mmap[0]==NULL && mmap[1]==NULL,TEMPLATE_MMAP_NULL);
#define GT_MMAP_CHECK(mmap) \
  GT_MMAP_ARRAY_CHECK(mmap->mmap)
#define GT_TEMPLATE_DICTIONARY_CHECK(template_dictionary) \
  GT_NULL_CHECK(template_dictionary); \
  GT_HASH_CHECK(template_dictionary->refs_dictionary)

/*
 * Reduction to single alignment
 */
#define GT_TEMPLATE_REDUCTION(template,reduced_alignment) \
  gt_alignment* const reduced_alignment = gt_template_get_block((template),0)
#define GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,reduced_alignment) \
  if (gt_expect_false(gt_template_get_num_blocks((template))==1)) { \
    GT_TEMPLATE_REDUCTION(template,reduced_alignment);
#define GT_TEMPLATE_IF_SE_ALINGMENT(template) \
  if (gt_expect_false(gt_template_get_num_blocks((template))==1))
#define GT_TEMPLATE_END_REDUCTION }
#define GT_TEMPLATE_END_REDUCTION__RETURN return;}
#define GT_TEMPLATE_REDUCTION_ELSE } else

#define GT_TEMPLATE_REDUCE_BOTH_ENDS(template,alignment_end1,alignment_end2) \
    gt_alignment* const alignment_end1 = gt_template_get_block((template),0); \
    gt_alignment* const alignment_end2 = gt_template_get_block((template),1)

/*
 * Setup
 */
GT_INLINE gt_template* gt_template_new(void);
GT_INLINE void gt_template_clear(gt_template* const template,const bool delete_alignments);
GT_INLINE void gt_template_delete(gt_template* const template);

/*
 * Accessors
 */
GT_INLINE char* gt_template_get_tag(gt_template* const template);
GT_INLINE gt_string* gt_template_get_string_tag(gt_template* const template);
GT_INLINE void gt_template_set_tag(gt_template* const template,char* const tag,const uint64_t length);
GT_INLINE uint64_t gt_template_get_total_length(gt_template* const template);
GT_INLINE int64_t gt_template_get_pair(gt_template* const template);

/* Blocks (template's alignments) */
GT_INLINE uint64_t gt_template_get_num_blocks(gt_template* const template);
GT_INLINE void gt_template_add_block(gt_template* const template,gt_alignment* const alignment);
GT_INLINE gt_alignment* gt_template_get_block(gt_template* const template,const uint64_t position);
GT_INLINE gt_alignment* gt_template_get_block_dyn(gt_template* const template,const uint64_t position);
GT_INLINE void gt_template_delete_blocks(gt_template* const template);
/* SingleEnd/PairedEnd API (just adaptors) */
#define gt_template_get_end1(template) ((gt_alignment*)gt_template_get_block(template,0))
#define gt_template_get_end2(template) ((gt_alignment*)gt_template_get_block(template,1))
GT_INLINE void gt_template_set_end1(gt_template* const template,gt_alignment* const alignment);
GT_INLINE void gt_template_set_end2(gt_template* const template,gt_alignment* const alignment);
#define gt_template_is_single_end(template) (gt_template_get_num_blocks(template)==1)
#define gt_template_is_paired_end(template) (gt_template_get_num_blocks(template)==2)
/* Counters */
GT_INLINE gt_vector* gt_template_get_counters_vector(gt_template* const template);
GT_INLINE void gt_template_set_counters_vector(gt_template* const template,gt_vector* const counters);
GT_INLINE uint64_t gt_template_get_num_counters(gt_template* const template);
GT_INLINE uint64_t gt_template_get_counter(gt_template* const template,const uint64_t stratum);
GT_INLINE void gt_template_set_counter(gt_template* const template,const uint64_t stratum,const uint64_t value);
GT_INLINE void gt_template_dec_counter(gt_template* const template,const uint64_t stratum);
GT_INLINE void gt_template_inc_counter(gt_template* const template,const uint64_t stratum);

/*
 * Attribute accessors
 */
GT_INLINE uint64_t gt_template_get_mcs(gt_template* const template);
GT_INLINE void gt_template_set_mcs(gt_template* const template,uint64_t max_complete_strata);
GT_INLINE bool gt_template_has_qualities(gt_template* const template);
GT_INLINE bool gt_template_get_not_unique_flag(gt_template* const template);
GT_INLINE void gt_template_set_not_unique_flag(gt_template* const template,bool is_not_unique);
GT_INLINE void gt_template_set_mmap_primary(gt_template* const template,gt_map** const mmap);
GT_INLINE gt_map** gt_template_get_mmap_primary(gt_template* const template);

/*
 * Multi-maps handlers (Map's relation)
 */
GT_INLINE uint64_t gt_template_get_num_mmaps(gt_template* const template);
GT_INLINE void gt_template_clear_mmaps(gt_template* const template);
/* MMap attributes */
GT_INLINE void gt_template_mmap_attributes_clear(gt_mmap_attributes* const mmap_attributes);
/* MMap record */
GT_INLINE gt_mmap* gt_template_get_mmap(gt_template* const template,const uint64_t position);
GT_INLINE void gt_template_set_mmap(gt_template* const template,const uint64_t position,gt_mmap* const mmap);
GT_INLINE void gt_template_add_mmap(gt_template* const template,gt_mmap* const mmap);
/* MMap array */
GT_INLINE gt_map** gt_template_get_mmap_array(
    gt_template* const template,const uint64_t position,gt_mmap_attributes** mmap_attributes);
GT_INLINE void gt_template_set_mmap_array(
    gt_template* const template,const uint64_t position,gt_map** const mmap,gt_mmap_attributes* const mmap_attributes);
GT_INLINE void gt_template_add_mmap_array(
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attributes);
/* MMap individual ends*/
GT_INLINE void gt_template_set_mmap_ends(
    gt_template* const template,const uint64_t position,
    gt_map* const map_end1,gt_map* const map_end2,gt_mmap_attributes* const mmap_attributes);
GT_INLINE void gt_template_add_mmap_ends(
    gt_template* const template,
    gt_map* const map_end1,gt_map* const map_end2,gt_mmap_attributes* const mmap_attributes);
/* MMap gt_vector */
GT_INLINE void gt_template_get_mmap_gtvector(
    gt_template* const template,const uint64_t position,gt_vector* const mmap,gt_mmap_attributes* const mmap_attributes);
GT_INLINE void gt_template_add_mmap_gtvector(
    gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attributes);

/**
 * indexing dictionary
 */
GT_INLINE gt_template_dictionary* gt_template_dictionary_new(gt_template* const template);
GT_INLINE void gt_template_dictionary_delete(gt_template_dictionary* const template_dictionary);
GT_INLINE void gt_template_dictionary_add_ref(gt_template_dictionary* const template_dictionary, gt_template* const template);

/*
 * Miscellaneous
 */
GT_INLINE gt_template* gt_template_dup(gt_template* const template,const bool copy_maps,const bool copy_mmaps);
GT_INLINE void gt_template_copy(gt_template* const template_dst,gt_template* const template_src,const bool copy_maps,const bool copy_mmaps);
GT_INLINE void gt_template_swap(gt_template* const template_a,gt_template* const template_b);

/*
 * Template's Alignments iterator (end1,end2, ... )
 */
GT_INLINE void gt_template_new_alignment_iterator(gt_template* const template,gt_template_alignment_iterator* const template_alignment_iterator);
GT_INLINE gt_alignment* gt_template_next_alignment(gt_template_alignment_iterator* const template_alignment_iterator);
GT_INLINE uint64_t gt_template_next_alignment_pos(gt_template_alignment_iterator* const template_alignment_iterator);
/*
 * Template's Multi-maps iterator ( (end1:map1,end2:map1) , (end1:map2,end2:map2) , ... )
 */
GT_INLINE void gt_template_new_mmap_iterator(gt_template* const template,gt_template_maps_iterator* const template_maps_iterator);
GT_INLINE gt_status gt_template_next_mmap(
    gt_template_maps_iterator* const template_maps_iterator,
    gt_map*** const mmap,gt_mmap_attributes** const mmap_attributes);
GT_INLINE uint64_t gt_template_next_mmap_pos(gt_template_maps_iterator* const template_maps_iterator);

/*
 * Iterate over the map(s) of the template
 *   (Eg. Single End => maps)
 *   (Eg. Paired Alignment => pairs of maps (map_end1,map_end2) )
 *   Template = Alignment_end1 + Alignment_end2 + {(end1.map1,end2.map1),(end1.map2,end2.map2),...}
 *   GT_TEMPLATE_ITERATE(template) := {(end1.map1,end2.map1),(end1.map2,end2.map2)}
 *   GT_TEMPLATE_ITERATE_(template,map_array) {
 *     ..code(map_array)..
 *   }
 * Also..
 *   GT_TEMPLATE_ITERATE_MMAP(template,map_array) // Just an explicit alias of GT_TEMPLATE_ITERATE
 *   GT_TEMPLATE_ITERATE_MMAP__ATTR_(template,map_array,map_array_attr) {
 *     .. to access mmap attributes {distance, score, ...}
 *   }
 */
#define GT_TEMPLATE_ITERATE_(template,mmap) \
  gt_map** mmap; \
  gt_template_maps_iterator __##template##_maps_iterator; \
  gt_template_new_mmap_iterator(template,&(__##template##_maps_iterator)); \
  while (gt_template_next_mmap(&(__##template##_maps_iterator),&mmap,NULL))
#define GT_TEMPLATE_ITERATE(template,mmap) \
  const uint64_t __##mmap##_num_blocks = gt_template_get_num_blocks(template); \
  GT_TEMPLATE_ITERATE_(template,mmap)
#define GT_TEMPLATE_ITERATE_MMAP(template,mmap) GT_TEMPLATE_ITERATE(template,mmap) /* Just an alias */
#define GT_TEMPLATE_ITERATE_MMAP_(template,mmap) GT_TEMPLATE_ITERATE_(template,mmap) /* Just an alias */
#define GT_TEMPLATE_ITERATE_MMAP__ATTR_(template,mmap,mmap_attribute) \
  gt_map** mmap; \
  gt_mmap_attributes *mmap_attribute = NULL; \
  gt_template_maps_iterator __##template##_maps_iterator; \
  gt_template_new_mmap_iterator(template,&(__##template##_maps_iterator)); \
  while (gt_template_next_mmap(&(__##template##_maps_iterator),&mmap,&mmap_attribute))
#define GT_TEMPLATE_ITERATE_MMAP__ATTR(template,mmap,map_array_attr) \
  const uint64_t __##mmap##_num_blocks = gt_template_get_num_blocks(template); \
  GT_TEMPLATE_ITERATE_MMAP__ATTR_(template,mmap,map_array_attr)

/*
 * Iterate over array of maps (mmap)
 *   map_array = (end1.map1,end2.map1)
 *   GT_MMAP_ITERATE_ENDS(map_array) := {end1.map1,end2.map1}
 *
 *   GT_MMAP_ITERATE_ENDS(mmap_array,map,end_position) {
 *     ..code(map,end_position)..
 *   }
 *
 * Used in conjunction with {GT_TEMPLATE_ITERATE,GT_TEMPLATE_ITERATE_MMAP,GT_TEMPLATE_ITERATE_MMAP__ATTR}
 * can iterate over the map ends of a paired template. i.e
 *   GT_TEMPLATE_ITERATE(template,mmap) {
 *     GT_MMAP_ITERATE(mmap,map,end_position) {
 *       ..code(map,end_position)..
 *     }
 *   }
 */
#define GT_MMAP_ITERATE_ENDS(mmap,num_blocks,map,end_position) \
  const uint64_t _num_blocks_##mmap = num_blocks; \
  uint64_t end_position; \
  gt_map* map; \
  for (end_position=0,map=mmap[0];end_position<(_num_blocks_##mmap);map=mmap[++end_position])
#define GT_MMAP_ITERATE(mmap,map,end_position) \
  GT_MMAP_ITERATE_ENDS(mmap,__##mmap##_num_blocks,map,end_position)

/*
 * Iterate over the alignment of a template (individual blocks)
 *   Template = Alignment_end1 + Alignment_end2 + {(end1.map1,end2.map1),(end1.map2,end2.map2),...}
 *   GT_TEMPLATE_ITERATE_ALIGNMENT(template) := {Alignment_end1,Alignment_end2}
 */
#define GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) \
  /* Template. Iterate over all alignments */ \
  gt_template_alignment_iterator alignment##_iterator; \
  gt_alignment* alignment; \
  gt_template_new_alignment_iterator(template,&(alignment##_iterator)); \
  while ((alignment=gt_template_next_alignment(&(alignment##_iterator))))

/*
 * Error Messages. Template/Alignment/Map/Misms errors
 */
#define GT_ERROR_MISMS_TRANSITION "Incorrect mismatch transition (Same base at both sides)"
#define GT_ERROR_MISMS_SPLICE_POS "Splicing distance must be positive (non-zero)"
#define GT_ERROR_COUNTERS_POS_STRATUM "Stratum must be strictly positive (stratum>0)"
#define GT_ERROR_ALIGNMENT_READ_QUAL_LENGTH "Read and quality length differs"
#define GT_ERROR_ALIGNMENT_MAPS_NOT_PARSED "Alignment's maps not parsed yet"
#define GT_ERROR_ALIGNMENT_INCONSISTENT_COUNTERS "Alignment inconsistency. Maps inconsistent with counters values"
#define GT_ERROR_TEMPLATE_ZERO_BLOCKS "Zero alignment blocks (num_blocks_template==0)"
#define GT_ERROR_TEMPLATE_TOO_MANY_BLOCKS "Template contains already two ends"
#define GT_ERROR_TEMPLATE_BLOCKS_EXCESS "Template blocks exceeds two ends"
#define GT_ERROR_TEMPLATE_MMAP_NULL "Template mmap is null (all maps are null)"
#define GT_ERROR_TEMPLATE_INCONSISTENT_MMAPS_ALIGNMENT "Template inconsistency. Multimaps' members must be contained by single alignments"
#define GT_ERROR_TEMPLATE_INCONSISTENT_MMAPS_ATTRB_RELATION "Template inconsistency. Incorrect number of mmaps and mmaps' attributes"
#define GT_ERROR_TEMPLATE_INCONSISTENT_NUM_MAPS_RELATION "Template inconsistency. Incorrect number of matches' elements (check num_blocks_template)"
#define GT_ERROR_TEMPLATE_INCONSISTENT_NUM_BLOCKS "Template inconsistency. Number of blocks must be the same across templates"
#define GT_ERROR_TEMPLATE_INCONSISTENT_COUNTERS "Template inconsistency. MMaps inconsistent with counters values"
#define GT_ERROR_TEMPLATE_ADD_BAD_NUM_BLOCKS "Trying to add wrong number of blocks to the template"
#define GT_ERROR_PALIGN_BAD_NUM_BLOCKS "Invalid Paired-alignment. Wrong number of alignment blocks (%"PRIu64")"
#define GT_ERROR_TEMPLATE_NOT_SCORED "Alignments have no valid score.  MAPQ can not be calculated unless the map file has been processed using score_reads"

#endif /* GT_TEMPLATE_H_ */
