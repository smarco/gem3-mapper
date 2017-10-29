/*
 * PROJECT: GEM-Tools library
 * FILE: gt_sequence_archive.c
 * DATE: 3/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Data structures needed to store a dictionary of DNA-sequences(chromosomes,contigs,etc) indexed by tag
 */

#include "gt_sequence_archive.h"
#include "gt_gemIdx_loader.h"

#define GT_SEQ_ARCHIVE_NUM_BLOCKS 15000
#define GT_SEQ_ARCHIVE_BLOCK_SIZE GT_BUFFER_SIZE_256K
#define GT_SEQ_ARCHIVE_NUM_INITIAL_BED_INTERVALS 5

/*
 * SequenceARCHIVE Constructor
 */
GT_INLINE gt_sequence_archive* gt_sequence_archive_new(const gt_sequence_archive_t sequence_archive_type) {
  gt_sequence_archive* seq_archive = gt_alloc(gt_sequence_archive);
  seq_archive->sequences = gt_shash_new();
  seq_archive->sequence_archive_type = sequence_archive_type;
  if (sequence_archive_type == GT_BED_ARCHIVE) {
    seq_archive->bed_intervals = gt_shash_new();
  }
  seq_archive->mm = NULL;
  return seq_archive;
}
GT_INLINE void gt_sequence_archive_clear(gt_sequence_archive* const seq_archive) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  // Clear Sequences
  GT_SHASH_BEGIN_ELEMENT_ITERATE(seq_archive->sequences,sequence,gt_segmented_sequence) {
    gt_segmented_sequence_delete(sequence);
  } GT_SHASH_END_ITERATE;
  gt_shash_clear(seq_archive->sequences,false);
  if (seq_archive->sequence_archive_type == GT_BED_ARCHIVE) {
    GT_SHASH_BEGIN_ELEMENT_ITERATE(seq_archive->bed_intervals,interval_vector,gt_vector) {
      gt_vector_clear(interval_vector);
    } GT_SHASH_END_ITERATE;
    gt_shash_clear(seq_archive->bed_intervals,false); // Clear BED intervals
  }
  // Free MM
  if (seq_archive->mm!=NULL) {
    gt_mm_free(seq_archive->mm);
    seq_archive->mm = NULL;
  }
}
GT_INLINE void gt_sequence_archive_delete(gt_sequence_archive* const seq_archive) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  // Delete Sequences
  GT_SHASH_BEGIN_ELEMENT_ITERATE(seq_archive->sequences,sequence,gt_segmented_sequence) {
    gt_segmented_sequence_delete(sequence);
  } GT_SHASH_END_ITERATE;
  gt_shash_delete(seq_archive->sequences,false);
  if (seq_archive->sequence_archive_type == GT_BED_ARCHIVE) {
    GT_SHASH_BEGIN_ELEMENT_ITERATE(seq_archive->bed_intervals,interval_vector,gt_vector) {
      gt_vector_delete(interval_vector);
    } GT_SHASH_END_ITERATE;
    gt_shash_delete(seq_archive->bed_intervals,false); // Clear BED intervals
  }
  // Free MM
  if (seq_archive->mm!=NULL) gt_mm_free(seq_archive->mm);
  // Free handler
  gt_free(seq_archive);
}

/*
 * SequenceARCHIVE handler
 */
GT_INLINE void gt_sequence_archive_add_segmented_sequence(gt_sequence_archive* const seq_archive,gt_segmented_sequence* const sequence) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  GT_SEGMENTED_SEQ_CHECK(sequence);
  gt_shash_insert(seq_archive->sequences,gt_string_get_string(sequence->seq_name),sequence,gt_segmented_sequence);
}
GT_INLINE void gt_sequence_archive_remove_segmented_sequence(gt_sequence_archive* const seq_archive,char* const seq_id) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  gt_segmented_sequence* seg_seq = gt_sequence_archive_get_segmented_sequence(seq_archive,seq_id);
  if (seg_seq != NULL) gt_segmented_sequence_delete(seg_seq);
  gt_shash_remove(seq_archive->sequences,seq_id,false);
}
GT_INLINE gt_segmented_sequence* gt_sequence_archive_get_segmented_sequence(gt_sequence_archive* const seq_archive,char* const seq_id) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  return gt_shash_get(seq_archive->sequences,seq_id,gt_segmented_sequence);
}
/* GT_BED_ARCHIVE */
GT_INLINE void gt_sequence_archive_add_bed_sequence(gt_sequence_archive* const seq_archive,gt_segmented_sequence* const sequence) {
  GT_SEQUENCE_BED_ARCHIVE_CHECK(seq_archive);
  GT_SEGMENTED_SEQ_CHECK(sequence);
  gt_shash_insert(seq_archive->sequences,gt_string_get_string(sequence->seq_name),sequence,gt_segmented_sequence);
}
GT_INLINE void gt_sequence_archive_remove_bed_sequence(gt_sequence_archive* const seq_archive,char* const seq_id) {
  GT_SEQUENCE_BED_ARCHIVE_CHECK(seq_archive);
  gt_segmented_sequence* seg_seq = gt_shash_get(seq_archive->sequences,seq_id,gt_segmented_sequence);
  if (seg_seq != NULL) gt_segmented_sequence_delete(seg_seq);
  gt_shash_remove(seq_archive->sequences,seq_id,false);
}
GT_INLINE gt_vector* gt_sequence_archive_get_bed_intervals_vector_dyn(gt_sequence_archive* const seq_archive,char* const seq_id) {
  GT_SEQUENCE_BED_ARCHIVE_CHECK(seq_archive);
  gt_vector* interval_vector = gt_sequence_archive_get_bed_intervals_vector(seq_archive,seq_id);
  if (interval_vector==NULL) {
    interval_vector = gt_vector_new(GT_SEQ_ARCHIVE_NUM_INITIAL_BED_INTERVALS,sizeof(gem_loc_t));
    gt_shash_insert(seq_archive->bed_intervals,seq_id,interval_vector,gt_vector*);
  }
  return interval_vector;
}
GT_INLINE gt_vector* gt_sequence_archive_get_bed_intervals_vector(gt_sequence_archive* const seq_archive,char* const seq_id) {
  GT_SEQUENCE_BED_ARCHIVE_CHECK(seq_archive);
  return gt_shash_get(seq_archive->bed_intervals,seq_id,gt_vector);
}

/*
 * SequenceARCHIVE High-level Retriever
 */
GT_INLINE gt_status gt_sequence_archive_get_sequence_string(
    gt_sequence_archive* const seq_archive,char* const seq_id,const gt_strand strand,
    const uint64_t position,const uint64_t length,gt_string* const string) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  GT_NULL_CHECK(seq_id);
  GT_ZERO_CHECK(length);
  GT_STRING_CHECK_NO_STATIC(string);
  gt_status error_code;
  // Retrieve the sequence
  gt_segmented_sequence* seg_seq = NULL;
  switch (seq_archive->sequence_archive_type) {
  case GT_CDNA_ARCHIVE:
    seg_seq = gt_sequence_archive_get_segmented_sequence(seq_archive,seq_id);
    if (seg_seq==NULL) return GT_SEQUENCE_NOT_FOUND;
    // Get the actual chunk
    if ((error_code=gt_segmented_sequence_get_sequence(seg_seq,position,length,string))) return error_code;
    break;
  case GT_BED_ARCHIVE:
    if ((error_code=gt_gemIdx_get_bed_sequence_string(seq_archive,seq_id,position,length,string)) < 0) {
      if (error_code==GT_GEMIDX_SEQ_NOT_FOUND) gt_error(GEMIDX_SEQ_ARCHIVE_NOT_FOUND,seq_id);
      if (error_code==GT_GEMIDX_INTERVAL_NOT_FOUND) gt_error(GEMIDX_INTERVAL_NOT_FOUND,seq_id);
      return -1;
    }
    break;
  default:
    gt_fatal_error(NOT_IMPLEMENTED);
    break;
  }
  // RC (if needed)
  if (strand==REVERSE) gt_dna_string_reverse_complement(string);
  return 0;
}
GT_INLINE gt_status gt_sequence_archive_retrieve_sequence_chunk(
    gt_sequence_archive* const seq_archive,char* const seq_id,const gt_strand strand,
    const uint64_t position,const uint64_t length,const uint64_t extra_length,gt_string* const string) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  GT_ZERO_CHECK(position);
  GT_NULL_CHECK(seq_id);
  GT_STRING_CHECK_NO_STATIC(string);
  gt_status error_code;
  // Retrieve the sequence
  gt_segmented_sequence* seg_seq = gt_sequence_archive_get_segmented_sequence(seq_archive,seq_id);
  if (seg_seq==NULL) {
    gt_error(SEQ_ARCHIVE_NOT_FOUND,seq_id);
    return GT_SEQUENCE_NOT_FOUND;
  }
  // Calculate sequence's boundaries
  const uint64_t sequence_total_length = seg_seq->sequence_total_length;
  uint64_t init_position, total_length;
  if (position-1 >= seg_seq->sequence_total_length) { // Check position
    gt_error(SEQ_ARCHIVE_POS_OUT_OF_RANGE,position-1);
    return GT_SEQUENCE_POS_OUT_OF_RANGE;
  }
  // Adjust init_position,total_length wrt strand
  init_position = position-1;
  if (strand==REVERSE) {
    init_position = (extra_length>init_position) ? 0 : init_position-extra_length;
  }
  total_length = length+extra_length;
  if (total_length >= sequence_total_length) total_length = seg_seq->sequence_total_length-1;
  // Get the sequence string
  switch (seq_archive->sequence_archive_type) {
  case GT_CDNA_ARCHIVE:
    // Get the actual chunk
    //  gt_status error_code; /* Error checking & reporting version */
    //  if ((error_code=gt_segmented_sequence_get_sequence(seg_seq,init_position,total_length,string))) {
    //    gt_error(SEQ_ARCHIVE_CHUNK_OUT_OF_RANGE,init_position,init_position+total_length,seq_id);
    //    return GT_SEQ_ARCHIVE_CHUNK_OUT_OF_RANGE;
    //  }
    gt_segmented_sequence_get_sequence(seg_seq,init_position,total_length,string);
    break;
  case GT_BED_ARCHIVE:
    if ((error_code=gt_gemIdx_get_bed_sequence_string(seq_archive,seq_id,init_position,total_length,string)) < 0) {
      if (error_code==GT_GEMIDX_SEQ_NOT_FOUND) gt_error(GEMIDX_SEQ_ARCHIVE_NOT_FOUND,seq_id);
      if (error_code==GT_GEMIDX_INTERVAL_NOT_FOUND) gt_error(GEMIDX_INTERVAL_NOT_FOUND,seq_id);
      return -1;
    }
    break;
  default:
    gt_fatal_error(NOT_IMPLEMENTED);
    break;
  }
  // RC (if needed)
  if (strand==REVERSE) gt_dna_string_reverse_complement(string);
  return 0;
}


/*
 * SequenceARCHIVE sorting functions
 */
int gt_sequence_archive_lexicographical_sort_fx(char *a,char *b) {
  /*
   * return (int) -1 if (a < b)
   * return (int)  0 if (a == b)
   * return (int)  1 if (a > b)
   */
  return gt_strcmp(a,b);
}
int gt_sequence_archive_karyotypic_sort_fx(char *a,char *b) {
  /*
   * Karyotypic order: 1, 2, ..., 10, 11, ..., 20, 21, 22, X, Y with M either leading or trailing these contigs
   */
  const uint64_t alen = strlen(a);
  const uint64_t alast = alen>1 ? alen-1 : 0;
  const uint64_t blen = strlen(b);
  const uint64_t blast = blen>1 ? blen-1 : 0;
  const int str_cmp_ab = gt_strcmp(a,b);
  if (str_cmp_ab==0) return 0;
  // Chromosome M
  if (strncmp(a,"chrM",4)==0 || a[alast]=='M' || a[alast]=='m') return INT16_MAX;
  if (strncmp(b,"chrM",4)==0 || b[blast]=='M' || b[blast]=='m') return -INT16_MAX;
  // Chromosome Y
  if (strncmp(a,"chrY",4)==0 || a[alast]=='Y' || a[alast]=='y') return INT16_MAX-1;
  if (strncmp(b,"chrY",4)==0 || b[blast]=='Y' || b[blast]=='y') return -(INT16_MAX-1);
  // Chromosome X
  if (strncmp(a,"chrX",4)==0 || a[alast]=='X' || a[alast]=='x') return INT16_MAX-2;
  if (strncmp(b,"chrX",4)==0 || b[blast]=='X' || b[blast]=='x') return -(INT16_MAX-2);
  // Other Chromosome
  return str_cmp_ab;
}

#define gt_string_cmp_wrapper(arg1,arg2) gt_string_cmp((char*)arg1,(char*)arg2)
#define gt_sequence_archive_lexicographical_sort_fx_wrapper(arg1,arg2) gt_sequence_archive_lexicographical_sort_fx((char*)arg1,(char*)arg2)
#define gt_sequence_archive_karyotypic_sort_fx_wrapper(arg1,arg2) gt_sequence_archive_karyotypic_sort_fx((char*)arg1,(char*)arg2)
GT_INLINE void gt_sequence_archive_sort(gt_sequence_archive* const seq_archive,int (*gt_string_cmp)(char*,char*)) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  HASH_SORT(seq_archive->sequences->shash_head,gt_string_cmp_wrapper);
}
GT_INLINE void gt_sequence_archive_lexicographical_sort(gt_sequence_archive* const seq_archive) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  HASH_SORT(seq_archive->sequences->shash_head,gt_sequence_archive_lexicographical_sort_fx_wrapper);
}
GT_INLINE void gt_sequence_archive_karyotypic_sort(gt_sequence_archive* const seq_archive) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  HASH_SORT(seq_archive->sequences->shash_head,gt_sequence_archive_karyotypic_sort_fx_wrapper);
}

/*
 * SequenceARCHIVE Iterator
 */
GT_INLINE void gt_sequence_archive_new_iterator(
    gt_sequence_archive* const seq_archive,gt_sequence_archive_iterator* const seq_archive_iterator) {
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive);
  GT_NULL_CHECK(seq_archive_iterator);
  seq_archive_iterator->sequence_archive = seq_archive;
  seq_archive_iterator->shash_it = seq_archive->sequences->shash_head;
}
GT_INLINE bool gt_sequence_archive_iterator_eos(gt_sequence_archive_iterator* const seq_archive_iterator) {
  GT_SEQUENCE_ARCHIVE_ITERATOR_CHECK(seq_archive_iterator);
  return seq_archive_iterator->shash_it==NULL;
}
GT_INLINE gt_segmented_sequence* gt_sequence_archive_iterator_next(gt_sequence_archive_iterator* const seq_archive_iterator) {
  GT_SEQUENCE_ARCHIVE_ITERATOR_CHECK(seq_archive_iterator);
  if (seq_archive_iterator->shash_it) {
    gt_segmented_sequence* elm =  seq_archive_iterator->shash_it->element;
    seq_archive_iterator->shash_it = seq_archive_iterator->shash_it->hh.next;
    return elm;
  } else {
    return NULL;
  }
}
