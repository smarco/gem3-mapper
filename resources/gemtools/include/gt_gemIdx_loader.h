/*
 * PROJECT: GEM-Tools library
 * FILE: gt_gemIdx_loader.h
 * DATE: 01/02/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_GEM_INDEX_LOADER_H_
#define GT_GEM_INDEX_LOADER_H_

#include "gt_commons.h"
#include "gt_compact_dna_string.h"
#include "gt_sequence_archive.h"

/*
 * Error Codes
 */
#define GT_GEMIDX_SEQ_NOT_FOUND      -1
#define GT_GEMIDX_INTERVAL_NOT_FOUND -2

/*
 * Auxiliary Data Structures (as to read the idx)
 */
typedef struct {
  uint64_t bot; // Global bottom location (over all text)
  uint64_t top; // Global top location (over all text)
  int64_t sequence_offset; // Offset relative to the sequence/chromosome
  int64_t tag_offset;
} gem_loc_t;

/*
 * Setup
 */
GT_INLINE void gt_gemIdx_load_archive(
    char* const index_file_name,gt_sequence_archive* const sequence_archive,const bool load_sequences);

/*
 * Retrieve sequences from GEMindex
 */
GT_INLINE int64_t gt_gemIdx_get_bed_sequence_string(
  gt_sequence_archive* const sequence_archive,char* const seq_id,
  const uint64_t position,const uint64_t length,gt_string* const string);

/*
 * Error Messages
 */
#define GT_ERROR_GEMIDX_SEQ_ARCHIVE_NOT_FOUND "GEMIdx. Sequence '%s' not found in reference archive"
#define GT_ERROR_GEMIDX_INTERVAL_NOT_FOUND "GEMIdx. Interval relative to sequence '%s' not found in reference archive"

#endif /* GT_GEM_INDEX_LOADER_H_ */
