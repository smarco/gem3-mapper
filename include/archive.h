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
#include "dna_text.h"
#include "fm_index.h"

#include "sequence.h"
#include "approximate_search.h"
#include "matches.h"

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

#define ARCHIVE_SEARCH_CHECK(archive_search) \
  GEM_CHECK_NULL(archive_search); \
  ARCHIVE_CHECK(archive_search->archive)

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
 * Archive Approximate Search
 */
typedef struct {
  /* Archive */
  archive_t* archive;
  /* Sequence */
  sequence_t* sequence;    // Input
  sequence_t* rc_sequence; // Generated
  /* Approximate Search */
  approximate_search_parameters_t search_parameters; // Search parameters
  approximate_search_t* forward_search_state;        // Forward Search State
  approximate_search_t* reverse_search_state;        // Reverse Search State
  /* Archive search control (Flow control) */
  bool search_reverse;
  bool probe_strand;
  /* Matches */
  matches_t* matches;
  /* MM */
  mm_stack_t* mm_stack; // Memory stack
} archive_search_t;

/*
 * Archive Loader
 */
GEM_INLINE archive_t* archive_read(char* const file_name,const bool do_tests,const bool verbose);
GEM_INLINE void archive_delete(archive_t* const archive);

/*
 * Archive Accessors
 */
GEM_INLINE uint64_t archive_get_index_length(const archive_t* const archive);

/*
 * Archive Search Setup
 */
GEM_INLINE archive_search_t* archive_search_new(archive_t* const archive);
GEM_INLINE void archive_search_delete(archive_search_t* const archive_search);
// [Initialize]
GEM_INLINE void archive_search_prepare_sequence(archive_search_t* const archive_search,sequence_t* const sequence);
// [Accessors]
GEM_INLINE approximate_search_parameters_t* archive_search_get_search_parameters(archive_search_t* const archive_search);
GEM_INLINE matches_t* archive_search_get_matches(archive_search_t* const archive_search);


/*
 * SingleEnd Indexed Search (SE Online Approximate String Search)
 */
GEM_INLINE void archive_search_single_end(archive_search_t* const archive_search,sequence_t* const sequence);

///*
// * PE Pairing (based on SE matches)
// */
//GEM_INLINE void archive_extend__pair_matches(
//    const archive_t* const archive,multimatches* const multimatches,
//    const uint64_t max_stratum,const uint64_t id_matches_from,const uint64_t id_matches_to,
//    fmi_extend_parameters* const extend_parameters,vector_pool* const mpool);
//

/*
 * Select Matches (Retrieving & Processing matches)
 *   - 1. Expand interval-matches (compacted)
 *   - 2. Transform CIGAR of reverse matches
 *   - 3. Sort matches wrt distance
 */
GEM_INLINE void archive_search_select_matches(
    archive_search_t* const archive_search,
    const uint64_t max_decoded_matches,const uint64_t min_decoded_strata,
    const uint64_t min_reported_matches,const uint64_t max_reported_matches);

/*
 * Decoding Matches Iterator (Decodes matches on the fly)
 */
typedef struct {
  // TODO
} matches_decode_iterator_t;
GEM_INLINE void matches_decode_iterator_new(
    matches_decode_iterator_t* const matches_iterator,archive_search_t* const archive_search,
    const uint64_t max_decoded_matches_stratum_wise,const uint64_t min_decoded_strata,
    const uint64_t min_decoded_matches,const uint64_t max_decoded_matches);

// TODO
GEM_INLINE void matches_decode_iterator_get_XXX(matches_decode_iterator_t* const matches_iterator);

GEM_INLINE void matches_decode_iterator_eoi(matches_decode_iterator_t* const matches_iterator);
GEM_INLINE void matches_decode_iterator_next(matches_decode_iterator_t* const matches_iterator);

/*
 * Display
 */
GEM_INLINE void archive_print(FILE* const stream,const archive_t* const archive);

/*
 * Error Messages
 */
#define GEM_ERROR_ARCHIVE_INDEX_OOB "Archive. Index position (%lu) out-of-bounds [0,%lu)]"

#endif /* ARCHIVE_H_ */
