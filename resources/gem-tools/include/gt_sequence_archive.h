/*
 * PROJECT: GEM-Tools library
 * FILE: gt_sequence_archive.h
 * DATE: 3/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Data structures needed to store a dictionary of DNA-sequences(chromosomes,contigs,etc) indexed by tag
 */

#ifndef GT_SEQUENCE_ARCHIVE_H_
#define GT_SEQUENCE_ARCHIVE_H_

#include "gt_essentials.h"
#include "gt_segmented_sequence.h"
#include "gt_dna_string.h"

/*
 * Codes gt_status
 */
#define GT_SEQUENCE_OK 0
#define GT_SEQUENCE_NOT_FOUND 10

/*
 * SequenceARCHIVE
 */
typedef enum { GT_CDNA_ARCHIVE, GT_BED_ARCHIVE } gt_sequence_archive_t; // TODO: GT_RAW_ARCHIVE
typedef struct {
  gt_sequence_archive_t sequence_archive_type;
  /* GT_CDNA_ARCHIVE */
  gt_shash* sequences; /* (gt_segmented_sequence*<gt_compact_dna_string>) */
  /* GT_BED_ARCHIVE */
  gt_shash* bed_intervals; /* (gt_vector*<gem_loc_t>) */
  uint64_t* bed; /* (GEMBitmap*) */
  /* Memory Management */
  gt_mm* mm;
} gt_sequence_archive;

typedef struct {
  gt_sequence_archive* sequence_archive;
  gt_shash_element *shash_it;
} gt_sequence_archive_iterator;

/*
 * Checkers
 */
#define GT_SEQUENCE_ARCHIVE_CHECK(seq_archive) \
    GT_NULL_CHECK(seq_archive); \
    GT_HASH_CHECK(seq_archive->sequences)
/*
 *  if (seq_archive->sequence_archive_type == GT_BED_ARCHIVE) { \
 *    GT_NULL_CHECK(seq_archive->bed); \
 *    GT_HASH_CHECK(seq_archive->bed_intervals); \
 *  }
 */
#define GT_SEQUENCE_CDNA_ARCHIVE_CHECK(seq_archive) \
    GT_NULL_CHECK(seq_archive); \
    GT_HASH_CHECK(seq_archive->sequences); \
    if (seq_archive->sequence_archive_type!=GT_CDNA_ARCHIVE) gt_fatal_error(SEQ_ARCHIVE_WRONG_TYPE)
#define GT_SEQUENCE_BED_ARCHIVE_CHECK(seq_archive) \
    GT_NULL_CHECK(seq_archive); \
    GT_HASH_CHECK(seq_archive->sequences); \
    if (seq_archive->sequence_archive_type!=GT_BED_ARCHIVE) gt_fatal_error(SEQ_ARCHIVE_WRONG_TYPE)
#define GT_SEQUENCE_ARCHIVE_ITERATOR_CHECK(seq_archive_iterator) \
  GT_SEQUENCE_ARCHIVE_CHECK(seq_archive_iterator->sequence_archive)

/*
 * SequenceARCHIVE Constructor
 */
GT_INLINE gt_sequence_archive* gt_sequence_archive_new(const gt_sequence_archive_t sequence_archive_type);
GT_INLINE void gt_sequence_archive_clear(gt_sequence_archive* const seq_archive);
GT_INLINE void gt_sequence_archive_delete(gt_sequence_archive* const seq_archive);
/*
 * SequenceARCHIVE handler
 */
/* GT_CDNA_ARCHIVE */
GT_INLINE void gt_sequence_archive_add_segmented_sequence(gt_sequence_archive* const seq_archive,gt_segmented_sequence* const sequence);
GT_INLINE void gt_sequence_archive_remove_segmented_sequence(gt_sequence_archive* const seq_archive,char* const seq_id);
GT_INLINE gt_segmented_sequence* gt_sequence_archive_get_segmented_sequence(gt_sequence_archive* const seq_archive,char* const seq_id);
/* GT_BED_ARCHIVE */
GT_INLINE void gt_sequence_archive_add_bed_sequence(gt_sequence_archive* const seq_archive,gt_segmented_sequence* const sequence);
GT_INLINE void gt_sequence_archive_remove_bed_sequence(gt_sequence_archive* const seq_archive,char* const seq_id);
GT_INLINE gt_vector* gt_sequence_archive_get_bed_intervals_vector_dyn(gt_sequence_archive* const seq_archive,char* const seq_id);
GT_INLINE gt_vector* gt_sequence_archive_get_bed_intervals_vector(gt_sequence_archive* const seq_archive,char* const seq_id);

/*
 * SequenceARCHIVE retrieve string sequences
 */
GT_INLINE gt_status gt_sequence_archive_get_sequence_string(
    gt_sequence_archive* const seq_archive,char* const seq_id,const gt_strand strand,
    const uint64_t position,const uint64_t length,gt_string* const string);
GT_INLINE gt_status gt_sequence_archive_retrieve_sequence_chunk(
    gt_sequence_archive* const seq_archive,char* const seq_id,const gt_strand strand,
    const uint64_t position,const uint64_t length,const uint64_t extra_length,gt_string* const string);

/*
 * SequenceARCHIVE sorting functions
 */
GT_INLINE void gt_sequence_archive_sort(gt_sequence_archive* const seq_archive,int (*gt_string_cmp)(char*,char*));
GT_INLINE void gt_sequence_archive_lexicographical_sort(gt_sequence_archive* const seq_archive);
GT_INLINE void gt_sequence_archive_karyotypic_sort(gt_sequence_archive* const seq_archive);
/*
 * SequenceARCHIVE Iterator
 */
GT_INLINE void gt_sequence_archive_new_iterator(
    gt_sequence_archive* const seq_archive,gt_sequence_archive_iterator* const seq_archive_iterator);
GT_INLINE bool gt_sequence_archive_iterator_eos(gt_sequence_archive_iterator* const seq_archive_iterator);
GT_INLINE gt_segmented_sequence* gt_sequence_archive_iterator_next(gt_sequence_archive_iterator* const seq_archive_iterator);
GT_INLINE gt_segmented_sequence* gt_sequence_archive_iterator_previous(gt_sequence_archive_iterator* const seq_archive_iterator);

/*
 * Error Messages. Sequence Archive
 */
#define GT_ERROR_SEQ_ARCHIVE_WRONG_TYPE "Wrong sequence archive type"
#define GT_ERROR_SEQ_ARCHIVE_NOT_FOUND "Sequence '%s' not found in reference archive"
#define GT_ERROR_SEQ_ARCHIVE_POS_OUT_OF_RANGE "Requested position '%"PRIu64"' out of sequence boundaries"
#define GT_ERROR_SEQ_ARCHIVE_CHUNK_OUT_OF_RANGE "Requested sequence string [%"PRIu64",%"PRIu64") out of sequence '%s' boundaries"

#endif /* GT_SEQUENCE_ARCHIVE_H_ */
