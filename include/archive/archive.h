/*
 * PROJECT: GEMMapper
 * FILE: archivex.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: GEM-index data structures and operators
 */

#ifndef ARCHIVE_H_
#define ARCHIVE_H_

#include "utils/essentials.h"
#include "archive/locator.h"
#include "archive/archive_text.h"
#include "fm_index/fm_index.h"

/*
 * Archive Model & Version
 */
#define ARCHIVE_MODEL_NO 5007ull

/*
 * Archive
 */
typedef enum {
  archive_dna=0,
  archive_dna_bisulfite=UINT64_MAX,
} archive_type;
typedef struct {
  // Meta-information
  archive_type type;                // Archive type
  bool indexed_complement;          // RC indexed explicitly
  uint64_t ns_threshold;            // Stretches of Ns (|Ns| >= ns_threshold) are not indexed
  bool indexed_reverse_text;        // Indexed reverse text (backwards text)
  // Locator
  locator_t* locator;               // Sequence Locator
  // Text
  archive_text_t* text;             // Archive Text
  // Index
  fm_index_t* fm_index;             // FM-Index
  // MM
  mm_t* mm;                         // MM
} archive_t;

/*
 * Archive Loader
 */
archive_t* archive_read(
    char* const file_name,
    const bool read_text_only);
void archive_delete(archive_t* const archive);

/*
 * Archive Accessors
 */
uint64_t archive_get_size(const archive_t* const archive);
uint64_t archive_get_index_length(const archive_t* const archive);

/*
 * Display
 */
void archive_print(
    FILE* const stream,
    const archive_t* const archive);

#endif /* ARCHIVE_H_ */
