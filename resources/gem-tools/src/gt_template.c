/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Data structure modeling sequences' templates.
 *   That is, set of alignments and relationships between their maps.
 */

#include "gt_template.h"
#include "gt_sam_attributes.h"

#define GT_TEMPLATE_TAG_INITIAL_LENGTH 100
#define GT_TEMPLATE_NUM_INITIAL_COUNTERS 10
#define GT_TEMPLATE_NUM_INITIAL_BLOCKS 2
#define GT_TEMPLATE_NUM_INITIAL_MMAPS 20

/*
 * Setup
 */
GT_INLINE gt_template* gt_template_new() {
  gt_template* template = gt_alloc(gt_template);
  template->template_id = UINT32_MAX;
  template->in_block_id = UINT32_MAX;
  template->tag = gt_string_new(GT_TEMPLATE_TAG_INITIAL_LENGTH);
  template->alignment_end1=NULL;
  template->alignment_end2=NULL;
  template->counters = gt_vector_new(GT_TEMPLATE_NUM_INITIAL_COUNTERS,sizeof(uint64_t));
  template->mmaps = gt_vector_new(GT_TEMPLATE_NUM_INITIAL_MMAPS,sizeof(gt_mmap));
  template->attributes = gt_attributes_new();
  template->alg_dictionary = NULL;
  return template;
}
GT_INLINE void gt_template_clear_handler(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  gt_string_clear(template->tag);
  gt_attributes_clear(template->attributes);
}
GT_INLINE void gt_template_clear(gt_template* const template,const bool delete_alignments) {
  GT_TEMPLATE_CHECK(template);
  if (delete_alignments) gt_template_delete_blocks(template);
  gt_vector_clear(template->counters);
  gt_vector_clear(template->mmaps);
  gt_template_clear_handler(template);
}
GT_INLINE void gt_template_delete(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  gt_string_delete(template->tag);
  gt_template_delete_blocks(template);
  gt_vector_delete(template->counters);
  gt_vector_delete(template->mmaps);
  gt_attributes_delete(template->attributes);
  gt_free(template);
}

/*
 * Accessors
 */
GT_INLINE char* gt_template_get_tag(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_tag(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return gt_string_get_string(template->tag);
}
GT_INLINE gt_string* gt_template_get_string_tag(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_string_tag(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return template->tag;
}
GT_INLINE void gt_template_set_tag(gt_template* const template,char* const tag,const uint64_t length) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(tag);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_set_tag(alignment,tag,length);
  } GT_TEMPLATE_END_REDUCTION;
  gt_string_set_nstring(template->tag,tag,length);
}
GT_INLINE uint64_t gt_template_get_total_length(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_read_length(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  uint64_t total_length = 0;
  GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
    total_length += gt_alignment_get_read_length(alignment);
  }
  return total_length;
}
GT_INLINE int64_t gt_template_get_pair(gt_template* const template){
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_pair(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return *((int64_t*)gt_attributes_get(template->attributes,GT_ATTR_ID_TAG_PAIR));
}
/* Blocks (single alignments) */
GT_INLINE uint64_t gt_template_get_num_blocks(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  return (gt_expect_true(template->alignment_end1!=NULL) ? (template->alignment_end2!=NULL ? 2 : 1 ) : 0 );
}
GT_INLINE void gt_template_add_block(gt_template* const template,gt_alignment* const alignment) {
  GT_TEMPLATE_CHECK(template);
  if (template->alignment_end1==NULL) {
    template->alignment_end1 = alignment;
  } else if (template->alignment_end2==NULL) {
    template->alignment_end2 = alignment;
  } else {
    gt_fatal_error(TEMPLATE_TOO_MANY_BLOCKS);
  }
}
GT_INLINE gt_alignment* gt_template_get_block(gt_template* const template,const uint64_t position) {
  GT_TEMPLATE_CHECK(template);
  gt_cond_fatal_error(position>1,TEMPLATE_BLOCKS_EXCESS);
  return (position==0 ? template->alignment_end1 : template->alignment_end2);
}
GT_INLINE gt_alignment* gt_template_get_block_dyn(gt_template* const template,const uint64_t position) {
  GT_TEMPLATE_CHECK(template);
  if (position==0) {
    if (template->alignment_end1==NULL) template->alignment_end1 = gt_alignment_new();
    return template->alignment_end1;
  } else if (position==1) {
    if (template->alignment_end1==NULL) template->alignment_end1 = gt_alignment_new();
    if (template->alignment_end2==NULL) template->alignment_end2 = gt_alignment_new();
    return template->alignment_end2;
  } else {
    gt_fatal_error(TEMPLATE_BLOCKS_EXCESS);
  }
}
GT_INLINE void gt_template_delete_blocks(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  if (template->alignment_end1!=NULL) {
    gt_alignment_delete(template->alignment_end1);
    template->alignment_end1=NULL;
  }
  if (template->alignment_end2!=NULL) {
    gt_alignment_delete(template->alignment_end2);
    template->alignment_end2=NULL;
  }
}
/* SingleEnd/PairedEnd API (just adaptors) */
GT_INLINE void gt_template_set_end1(gt_template* const template,gt_alignment* const alignment) {
  GT_TEMPLATE_CHECK(template);
  GT_ALIGNMENT_CHECK(alignment);
  template->alignment_end1=alignment;
}
GT_INLINE void gt_template_set_end2(gt_template* const template,gt_alignment* const alignment) {
  GT_TEMPLATE_CHECK(template);
  GT_ALIGNMENT_CHECK(alignment);
  template->alignment_end2=alignment;
}
/* Counters */
GT_INLINE gt_vector* gt_template_get_counters_vector(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_counters_vector(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return template->counters;
}
GT_INLINE void gt_template_set_counters_vector(gt_template* const template,gt_vector* const counters) {
  GT_TEMPLATE_CHECK(template);
  GT_VECTOR_CHECK(counters);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_set_counters_vector(alignment,counters);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  template->counters = counters;
}
GT_INLINE uint64_t gt_template_get_num_counters(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_num_counters(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return gt_vector_get_used(template->counters);
}
GT_INLINE uint64_t gt_template_get_counter(gt_template* const template,const uint64_t stratum) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_counter(alignment,stratum);
  } GT_TEMPLATE_END_REDUCTION;
  return *gt_vector_get_elm(template->counters,stratum,uint64_t);
}
GT_INLINE void gt_template_dynamically_allocate_counter(gt_template* const template,const uint64_t stratum) {
  GT_TEMPLATE_CHECK(template);
  const uint64_t used_strata = stratum+1;
  gt_vector_reserve(template->counters,used_strata,true);
  if (gt_vector_get_used(template->counters)<used_strata) {
    gt_vector_set_used(template->counters,used_strata);
  }
}
GT_INLINE void gt_template_inc_counter(gt_template* const template,const uint64_t stratum) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_inc_counter(alignment,stratum);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_template_dynamically_allocate_counter(template,stratum);
  ++(*gt_vector_get_elm(template->counters,stratum,uint64_t));
}
GT_INLINE void gt_template_dec_counter(gt_template* const template,const uint64_t stratum) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_dec_counter(alignment,stratum);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_template_dynamically_allocate_counter(template,stratum);
  --(*gt_vector_get_elm(template->counters,stratum,uint64_t));
}
GT_INLINE void gt_template_set_counter(gt_template* const template,const uint64_t stratum,const uint64_t value) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_set_counter(alignment,stratum,value);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_template_dynamically_allocate_counter(template,stratum);
  *gt_vector_get_elm(template->counters,stratum,uint64_t) = value;
}

/*
 * Predefined attributes
 */
GT_INLINE uint64_t gt_template_get_mcs(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_mcs(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  uint64_t* mcs = gt_attributes_get(template->attributes,GT_ATTR_ID_MAX_COMPLETE_STRATA);
  if (mcs == NULL) return UINT64_MAX;
  return *mcs;
}
GT_INLINE void gt_template_set_mcs(gt_template* const template,uint64_t max_complete_strata) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_set_mcs(alignment,max_complete_strata);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_attributes_add(template->attributes,GT_ATTR_ID_MAX_COMPLETE_STRATA,&max_complete_strata,uint64_t);
}
GT_INLINE bool gt_template_has_qualities(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_has_qualities(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  const uint64_t num_blocks = gt_template_get_num_blocks(template);
  uint64_t i;
  for (i=0;i<num_blocks;++i) {
    gt_alignment* alignment = gt_template_get_block(template,i);
    if (!gt_alignment_has_qualities(alignment)) return false;
  }
  return true;
}
GT_INLINE bool gt_template_get_not_unique_flag(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_not_unique_flag(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  bool* const not_unique_flag = gt_attributes_get(template->attributes,GT_ATTR_ID_NOT_UNIQUE);
  if (not_unique_flag==NULL) return false;
  return *not_unique_flag;
}
GT_INLINE void gt_template_set_not_unique_flag(gt_template* const template,bool is_not_unique) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_set_not_unique_flag(alignment,is_not_unique);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_attributes_add(template->attributes,GT_ATTR_ID_NOT_UNIQUE,&is_not_unique,bool);
}
GT_INLINE gt_map** gt_template_get_mmap_primary(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  // NOTE: No reduction performed
  return gt_attributes_get(template->attributes,GT_ATTR_ID_SAM_PRIMARY_ALIGNMENT);
}
GT_INLINE void gt_template_set_mmap_primary(gt_template* const template,gt_map** const mmap) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(mmap);
  // NOTE: No reduction performed
  gt_attributes_add(template->attributes,GT_ATTR_ID_SAM_PRIMARY_ALIGNMENT,mmap,gt_map**);
}

/*
 * Multi-maps handlers (Map's relation)
 */
GT_INLINE uint64_t gt_template_get_num_mmaps(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_alignment_get_num_maps(alignment);
  } GT_TEMPLATE_END_REDUCTION;
  return gt_vector_get_used(template->mmaps);
}
GT_INLINE void gt_template_clear_mmaps(gt_template* const template) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_clear_maps(alignment);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_vector_clear(template->mmaps);
}
/* MMap attributes */
GT_INLINE void gt_template_mmap_attributes_clear(gt_mmap_attributes* const mmap_attributes) {
  GT_NULL_CHECK(mmap_attributes);
  mmap_attributes->distance = 0;
  mmap_attributes->gt_score = GT_MAP_NO_GT_SCORE;
  mmap_attributes->phred_score = GT_MAP_NO_PHRED_SCORE;
}
/* MMap record */
GT_INLINE gt_mmap* gt_template_get_mmap(gt_template* const template,const uint64_t position) {
  GT_TEMPLATE_CHECK(template);
  return gt_vector_get_elm(template->mmaps,position,gt_mmap);
}
GT_INLINE void gt_template_set_mmap(gt_template* const template,const uint64_t position,gt_mmap* const mmap) {
  GT_TEMPLATE_CHECK(template);
  GT_MMAP_CHECK(mmap);
  gt_vector_set_elm(template->mmaps,position,gt_mmap,*mmap);
}
GT_INLINE void gt_template_add_mmap(gt_template* const template,gt_mmap* const mmap) {
  GT_TEMPLATE_CHECK(template);
  GT_MMAP_CHECK(mmap);
  gt_vector_insert(template->mmaps,*mmap,gt_mmap);
}
/* MMap array */
GT_INLINE gt_map** gt_template_get_mmap_array(
    gt_template* const template,const uint64_t position,gt_mmap_attributes** mmap_attributes) {
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_vector_get_elm(alignment->maps,position,gt_map*);
  } GT_TEMPLATE_END_REDUCTION;
  // Retrieve the mmap from the mmap vector
  gt_mmap* mmap_ph = gt_vector_get_elm(template->mmaps,position,gt_mmap);
  if (mmap_attributes) *mmap_attributes = &mmap_ph->attributes;
  return mmap_ph->mmap;
}
GT_INLINE void gt_template_set_mmap_array(
    gt_template* const template,const uint64_t position,gt_map** const mmap,gt_mmap_attributes* const mmap_attributes) {
  GT_TEMPLATE_CHECK(template);
  GT_MMAP_ARRAY_CHECK(mmap);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_set_map(alignment,mmap[0],position);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_mmap* mmap_ph = gt_vector_get_elm(template->mmaps,position,gt_mmap);
  mmap_ph->mmap[0] = mmap[0];
  mmap_ph->mmap[1] = mmap[1];
  if (mmap_attributes!=NULL) {
    mmap_ph->attributes = *mmap_attributes;
  } else {
    gt_template_mmap_attributes_clear(&mmap_ph->attributes);
  }
}
GT_INLINE void gt_template_add_mmap_array(
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attributes) {
  GT_TEMPLATE_CHECK(template);
  GT_MMAP_ARRAY_CHECK(mmap);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_add_map(alignment,mmap[0]);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Allocate new mmap element
  gt_vector_reserve_additional(template->mmaps,1);
  gt_vector_inc_used(template->mmaps);
  gt_mmap* mmap_ph = gt_vector_get_last_elm(template->mmaps,gt_mmap);
  // Store MMap array
  mmap_ph->mmap[0] = mmap[0];
  mmap_ph->mmap[1] = mmap[1];
  if (mmap_attributes!=NULL) {
    mmap_ph->attributes = *mmap_attributes;
  } else {
    gt_template_mmap_attributes_clear(&mmap_ph->attributes);
  }
}
/* MMap individual ends*/
GT_INLINE void gt_template_set_mmap_ends(
    gt_template* const template,const uint64_t position,
    gt_map* const map_end1,gt_map* const map_end2,gt_mmap_attributes* const mmap_attributes) {
  GT_TEMPLATE_CHECK(template);
  gt_cond_fatal_error(map_end1==NULL && map_end2==NULL,TEMPLATE_MMAP_NULL);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_set_map(alignment,map_end1,position);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  gt_mmap* mmap_ph = gt_vector_get_elm(template->mmaps,position,gt_mmap);
  mmap_ph->mmap[0] = map_end1;
  mmap_ph->mmap[1] = map_end2;
  if (mmap_attributes!=NULL) {
    mmap_ph->attributes = *mmap_attributes;
  } else {
    gt_template_mmap_attributes_clear(&mmap_ph->attributes);
  }
}
GT_INLINE void gt_template_add_mmap_ends(
    gt_template* const template,
    gt_map* const map_end1,gt_map* const map_end2,gt_mmap_attributes* const mmap_attributes) {
  GT_TEMPLATE_CHECK(template);
  gt_cond_fatal_error(map_end1==NULL && map_end2==NULL,TEMPLATE_MMAP_NULL);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_add_map(alignment,map_end1);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Allocate new mmap element
  gt_vector_reserve_additional(template->mmaps,1);
  gt_vector_inc_used(template->mmaps);
  gt_mmap* mmap_ph = gt_vector_get_last_elm(template->mmaps,gt_mmap);
  // Store MMap array
  mmap_ph->mmap[0] = map_end1;
  mmap_ph->mmap[1] = map_end2;
  if (mmap_attributes!=NULL) {
    mmap_ph->attributes = *mmap_attributes;
  } else {
    gt_template_mmap_attributes_clear(&mmap_ph->attributes);
  }
}
/* MMap gt_vector */
GT_INLINE void gt_template_get_mmap_gtvector(
    gt_template* const template,const uint64_t position,gt_vector* const mmap,gt_mmap_attributes* const mmap_attributes) {
  GT_TEMPLATE_CHECK(template);
  GT_VECTOR_CHECK(mmap);
  // Handle reduction to alignment
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_vector_prepare(mmap,gt_map*,1); // Reset output vector
    gt_vector_insert(mmap,gt_alignment_get_map(alignment,position),gt_map*);
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Retrieve the maps from the mmap vector
  gt_vector_prepare(mmap,gt_map*,2); // Reset output vector
  gt_mmap* mmap_ph = gt_vector_get_elm(template->mmaps,position,gt_mmap);
  gt_vector_insert(mmap,mmap_ph->mmap[0],gt_map*);
  gt_vector_insert(mmap,mmap_ph->mmap[1],gt_map*);
  if (mmap_attributes) *mmap_attributes = mmap_ph->attributes;
}
GT_INLINE void gt_template_add_mmap_gtvector(
    gt_template* const template,gt_vector* const mmap_vector,gt_mmap_attributes* const mmap_attributes) {
  GT_TEMPLATE_CHECK(template);
  GT_VECTOR_CHECK(mmap_vector);
  // Handle reduction to alignment
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    gt_alignment_add_map(alignment,*gt_vector_get_elm(mmap_vector,0,gt_map*));
  } GT_TEMPLATE_END_REDUCTION__RETURN;
  // Allocate new mmap element
  gt_vector_reserve_additional(template->mmaps,1);
  gt_vector_inc_used(template->mmaps);
  gt_mmap* mmap_ph = gt_vector_get_last_elm(template->mmaps,gt_mmap);
  // Store MMap array
  mmap_ph->mmap[0] = *gt_vector_get_elm(mmap_vector,0,gt_map*);
  mmap_ph->mmap[1] = *gt_vector_get_elm(mmap_vector,1,gt_map*);
  if (mmap_attributes!=NULL) {
    mmap_ph->attributes = *mmap_attributes;
  } else {
    gt_template_mmap_attributes_clear(&mmap_ph->attributes);
  }
}
/**
 * dictionary based indexing
 */
GT_INLINE gt_template_dictionary* gt_template_dictionary_new(gt_template* const template){
  gt_template_dictionary* template_dictionary = gt_alloc(gt_template_dictionary);
  template_dictionary->refs_dictionary = gt_shash_new();
  template_dictionary->template = template;
  return template_dictionary;
}
GT_INLINE void gt_template_dictionary_delete(gt_template_dictionary* const template_dictionary){
  GT_TEMPLATE_DICTIONARY_CHECK(template_dictionary);
  GT_SHASH_BEGIN_ELEMENT_ITERATE(template_dictionary->refs_dictionary,alg_dicc_elem,gt_vector) {
    GT_VECTOR_ITERATE(alg_dicc_elem, e, c, gt_template_dictionary_map_element*){
      gt_free(*e);
    }
    gt_vector_delete(alg_dicc_elem);
  } GT_SHASH_END_ITERATE;
  gt_shash_delete(template_dictionary->refs_dictionary,false);
  gt_free(template_dictionary);
}

GT_INLINE void gt_template_dictionary_add_ref(gt_template_dictionary* const template_dictionary, gt_template* const template){
  GT_TEMPLATE_DICTIONARY_CHECK(template_dictionary);
  uint64_t pos = 0;
  GT_TEMPLATE_ITERATE_MMAP__ATTR_(template, mmap, mmap_attr){
    char* ref = gt_map_get_seq_name(mmap[0]);
    gt_vector* dictionary_elements = NULL;
    if(!gt_shash_is_contained(template_dictionary->refs_dictionary, ref)){
      dictionary_elements = gt_vector_new(16, sizeof(gt_template_dictionary_map_element*));
      gt_shash_insert(template_dictionary->refs_dictionary, ref, dictionary_elements, gt_vector);
    }else{
      dictionary_elements = gt_shash_get(template_dictionary->refs_dictionary, ref, gt_vector);
    }
    gt_template_dictionary_map_element* e = gt_alloc(gt_template_dictionary_map_element);
    e->mmap = mmap;
    e->mmap_attrs = mmap_attr;
    e->pos = pos;
    gt_vector_insert(dictionary_elements, e, gt_template_dictionary_map_element*);
    pos++;
  }
}


/*
 * Miscellaneous
 */
GT_INLINE void gt_template_copy_handler(gt_template* template_dst,gt_template* const template_src) {
  GT_TEMPLATE_CHECK(template_src);
  template_dst->template_id = template_src->template_id;
  template_dst->in_block_id = template_src->in_block_id;
  gt_string_copy(template_dst->tag,template_src->tag);
  gt_attributes_copy(template_dst->attributes,template_src->attributes); // Copy templates' attributes
}
GT_INLINE void gt_template_copy_blocks(gt_template* template_dst,gt_template* const template_src,const bool copy_maps) {
  GT_TEMPLATE_CHECK(template_src);
  gt_template_delete_blocks(template_dst);
  const uint64_t num_blocks = gt_template_get_num_blocks(template_src);
  uint64_t i;
  for (i=0;i<num_blocks;++i) {
    gt_alignment* const alignment_copy = gt_alignment_copy(gt_template_get_block(template_src,i),copy_maps);
    gt_template_add_block(template_dst,alignment_copy);
  }
}
GT_INLINE gt_template* gt_template_dup(gt_template* const template,const bool copy_maps,const bool copy_mmaps) {
  GT_TEMPLATE_CHECK(template);
  gt_template* template_cpy = gt_template_new();
  gt_template_copy(template_cpy,template,copy_maps,copy_mmaps);
  return template_cpy;
}
GT_INLINE void gt_template_copy(gt_template* const template_dst,gt_template* const template_src,const bool copy_maps,const bool copy_mmaps) {
  GT_TEMPLATE_CHECK(template_dst);
  GT_TEMPLATE_CHECK(template_src);
  // Copy handler
  gt_template_copy_handler(template_dst,template_src);
  // Copy blocks
  gt_template_copy_blocks(template_dst,template_src,copy_maps);
  // Copy mmaps
  if (copy_maps && copy_mmaps && gt_template_get_num_blocks(template_src)>1) {
    // Copy counters
    gt_vector_copy(template_dst->counters,template_src->counters);
    // Copy mmaps & mmaps_attributes
    GT_TEMPLATE_ITERATE_MMAP__ATTR_(template_src,mmap,mmap_attributes) {
      // Locate maps
      uint64_t map_pos_end1, map_pos_end2;
      gt_cond_fatal_error(!gt_alignment_locate_map_reference(
          gt_template_get_end1(template_src),mmap[0],&map_pos_end1),TEMPLATE_INCONSISTENT_MMAPS_ALIGNMENT);
      gt_cond_fatal_error(!gt_alignment_locate_map_reference(
          gt_template_get_end2(template_src),mmap[1],&map_pos_end2),TEMPLATE_INCONSISTENT_MMAPS_ALIGNMENT);
      gt_map* const homologe_map_end1 = gt_alignment_get_map(gt_template_get_end1(template_dst),map_pos_end1);
      gt_map* const homologe_map_end2 = gt_alignment_get_map(gt_template_get_end2(template_dst),map_pos_end2);
      // Add mmap
      gt_template_add_mmap_ends(template_dst,homologe_map_end1,homologe_map_end2,mmap_attributes);
    }
  }
}
GT_INLINE void gt_template_swap(gt_template* const template_a,gt_template* const template_b) {
  GT_SWAP(template_a->template_id,template_b->template_id);
  GT_SWAP(template_a->in_block_id,template_b->in_block_id);
  GT_SWAP(template_a->tag,template_b->tag);
  GT_SWAP(template_a->alignment_end1,template_b->alignment_end1);
  GT_SWAP(template_a->alignment_end2,template_b->alignment_end2);
  GT_SWAP(template_a->counters,template_b->counters);
  GT_SWAP(template_a->mmaps,template_b->mmaps);
  GT_SWAP(template_a->attributes,template_b->attributes);
}
/*
 * Template's Alignments iterator (end1,end2, ... )
 */
GT_INLINE void gt_template_new_alignment_iterator(
    gt_template* const template,gt_template_alignment_iterator* const template_alignment_iterator) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(template_alignment_iterator);
  template_alignment_iterator->template = template;
  template_alignment_iterator->next_alignment_end = 0;
}
GT_INLINE gt_alignment* gt_template_next_alignment(gt_template_alignment_iterator* const template_alignment_iterator) {
  GT_NULL_CHECK(template_alignment_iterator);
  GT_TEMPLATE_CHECK(template_alignment_iterator->template);
  gt_template* const template = template_alignment_iterator->template;
  if (gt_expect_true(template_alignment_iterator->next_alignment_end < gt_template_get_num_blocks(template))) {
    gt_alignment* const alignment = gt_template_get_block(template,template_alignment_iterator->next_alignment_end);
    ++template_alignment_iterator->next_alignment_end;
    return alignment;
  } else {
    return NULL;
  }
}
GT_INLINE uint64_t gt_template_next_alignment_pos(gt_template_alignment_iterator* const template_alignment_iterator) {
  GT_NULL_CHECK(template_alignment_iterator);
  GT_TEMPLATE_CHECK(template_alignment_iterator->template);
  return template_alignment_iterator->next_alignment_end;
}

/*
 * Template's Maps iterator ( (end1:map1,end2:map1) , (end1:map2,end2:map2) , ... )
 */
GT_INLINE void gt_template_new_mmap_iterator(
    gt_template* const template,gt_template_maps_iterator* const template_maps_iterator) {
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(template_maps_iterator);
  template_maps_iterator->template = template;
  template_maps_iterator->next_mmap_position = 0;
}
GT_INLINE gt_status gt_template_next_mmap(
    gt_template_maps_iterator* const template_maps_iterator,
    gt_map*** const mmap,gt_mmap_attributes** const mmap_attributes) {
  GT_NULL_CHECK(template_maps_iterator);
  GT_TEMPLATE_CHECK(template_maps_iterator->template);
  GT_NULL_CHECK(mmap);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_maps_iterator->template,alignment) {
    if (gt_expect_true(template_maps_iterator->next_mmap_position < gt_alignment_get_num_maps(alignment))) {
      *mmap = gt_vector_get_elm(alignment->maps,template_maps_iterator->next_mmap_position,gt_map*);
      ++template_maps_iterator->next_mmap_position;
      return GT_TEMPLATE_OK;
    } else {
      *mmap = NULL;
      return GT_TEMPLATE_FAIL;
    }
  } GT_TEMPLATE_END_REDUCTION;
  gt_template* const template = template_maps_iterator->template;
  if (gt_expect_true(template_maps_iterator->next_mmap_position<gt_vector_get_used(template->mmaps))) {
    gt_mmap* const mmap_ph = gt_vector_get_elm(template->mmaps,template_maps_iterator->next_mmap_position,gt_mmap);
    *mmap = mmap_ph->mmap;
    if (mmap_attributes) *mmap_attributes = &mmap_ph->attributes;
    ++template_maps_iterator->next_mmap_position;
    return GT_TEMPLATE_OK;
  } else {
    *mmap = NULL;
    return GT_TEMPLATE_FAIL;
  }
}
GT_INLINE uint64_t gt_template_next_mmap_pos(gt_template_maps_iterator* const template_maps_iterator) {
  GT_NULL_CHECK(template_maps_iterator);
  GT_TEMPLATE_CHECK(template_maps_iterator->template);
  return template_maps_iterator->next_mmap_position;
}
