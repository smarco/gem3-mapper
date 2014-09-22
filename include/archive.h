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
#include "sampled_rl.h"
#include "graph_text.h"
#include "fm_index.h"

/*
 * Constants
 */
#define GEM_INDEX_VERSION 5

/*
 * Checker
 */
#define ARCHIVE_CHECK(archive) \
  GEM_CHECK_NULL(archive); \
  LOCATOR_CHECK(archive->locator); \
  DNA_TEXT_CHECK(archive->enc_text); \
  FM_INDEX_CHECK(archive->fm_index)
// FIXME GRAPH_TEXT_CHECK(archive->graph);

#define ARCHIVE_CHECK_INDEX_POSITION(archive,index_position) \
  gem_fatal_check((int64_t)index_position < 0 || (int64_t)index_position >= archive_get_index_length(archive), \
      ARCHIVE_INDEX_OOB,index_position,archive_get_index_length(archive))

/*
 * Archive
 */
typedef enum { fm_dna_classic=0, fm_dna_run_length=1, fm_dna_graph=2 } index_t;
typedef enum { Iupac_dna=0, Iupac_colorspace_dna=1 } filter_t; // TODO Pass_through
typedef struct {
  // Meta-information
  index_t index_type;         // Index/Archive type (architecture)
  filter_t filter_type;       // Filter applied to the original text (MFasta)
  bool indexed_complement;    // RC indexed explicitly
  uint64_t ns_threshold;      // Stretches of Ns equal or larger than the threshold have not been indexed
  // Locator
  locator_t* locator;         // Sequence Locator
  // Graph
  graph_text_t* graph;        // Graph (text + graph = hypertext)
  // Indexed text
  dna_text_t* enc_text;       // Index-Text
  sampled_rl_t* sampled_rl;   // Sampled RL-Index Positions
  // Index
  fm_index_t* fm_index;       // FM-Index
  // MM
  mm_t* mm;
} archive_t;

/*
 * Archive Loader
 */
GEM_INLINE archive_t* archive_read(char* const file_name,const bool do_tests,const bool verbose);
GEM_INLINE void archive_delete(archive_t* const archive);

/*
 * Archive Accessors
 */
GEM_INLINE uint64_t archive_get_size(const archive_t* const archive);
GEM_INLINE uint64_t archive_get_index_length(const archive_t* const archive);
GEM_INLINE bool archive_is_indexed_complement(const archive_t* const archive);

/*
 * Display
 */
GEM_INLINE void archive_print(FILE* const stream,const archive_t* const archive);

/*
 * Error Messages
 */
#define GEM_ERROR_ARCHIVE_INDEX_OOB "Archive. Index position (%lu) out-of-bounds [0,%lu)]"

#endif /* ARCHIVE_H_ */
