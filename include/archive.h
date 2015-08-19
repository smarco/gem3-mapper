/*
 * PROJECT: GEMMapper
 * FILE: archivex.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: GEM-index data structures and operators
 */

#ifndef ARCHIVE_H_
#define ARCHIVE_H_

#include "essentials.h"

#include "locator.h"
#include "archive_text.h"
#include "fm_index.h"

/*
 * Archive Model & Version
 */
#define ARCHIVE_MODEL_NO 5000ull

/*
 * Checker
 */
#define ARCHIVE_CHECK(archive) \
  GEM_CHECK_NULL(archive); \
  LOCATOR_CHECK(archive->locator); \
  ARCHIVE_TEXT_CHECK(archive->text); \
  FM_INDEX_CHECK(archive->fm_index)

#define ARCHIVE_CHECK_INDEX_POSITION(archive,index_position) \
  gem_fatal_check((int64_t)index_position < 0 || (int64_t)index_position >= archive_get_index_length(archive), \
      ARCHIVE_INDEX_OOB,index_position,archive_get_index_length(archive))

/*
 * Archive
 */
typedef enum {
  Iupac_dna=0,
  Iupac_colorspace_dna=1,
  Iupac_range=UINT64_MAX
} archive_filter_type;
typedef struct {
  // Meta-information
  archive_filter_type filter_type;  // Filter applied to the original text (MFasta)
  bool indexed_complement;          // RC indexed explicitly
  uint64_t ns_threshold;            // Stretches of Ns equal or larger than the threshold have not been indexed
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
archive_t* archive_read(char* const file_name,const bool do_tests,const bool verbose);
void archive_delete(archive_t* const archive);

/*
 * Archive Accessors
 */
uint64_t archive_get_size(const archive_t* const archive);
uint64_t archive_get_index_length(const archive_t* const archive);

/*
 * Display
 */
void archive_print(FILE* const stream,const archive_t* const archive);

/*
 * Error Messages
 */
#define GEM_ERROR_ARCHIVE_WRONG_MODEL_NO "Archive error. Wrong GEM-Index Model %"PRIu64" (Expected model %"PRIu64"). Please rebuild index (gem-indexer)"
#define GEM_ERROR_ARCHIVE_INDEX_OOB "Archive error. Index position (%"PRIu64") out-of-bounds [0,%"PRIu64")]"

#endif /* ARCHIVE_H_ */
