/*
 * PROJECT: GEMMapper
 * FILE: archive_search.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_SEARCH_H_
#define ARCHIVE_SEARCH_H_

#include "essentials.h"

#include "archive.h"
#include "sequence.h"
#include "matches.h"
#include "approximate_search.h"

//#include "fm_index.h"
//#include "locator.h"
//#include "sampled_rl.h"
//#include "graph_text.h"
//#include "dna_text.h"

/*
 * Checker
 */
#define ARCHIVE_SEARCH_CHECK(archive_search) \
  GEM_CHECK_NULL(archive_search); \
  ARCHIVE_CHECK(archive_search->archive)

/*
 * Archive Approximate-Search
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
  mm_stack_t* mm_stack;
} archive_search_t;

/*
 * Archive Search Setup
 */
GEM_INLINE archive_search_t* archive_search_new(archive_t* const archive,mm_stack_t* const mm_stack);
GEM_INLINE void archive_search_clear(archive_search_t* const archive_search);
GEM_INLINE void archive_search_delete(archive_search_t* const archive_search);
GEM_INLINE void archive_search_set_mm_stack(archive_search_t* const archive_search,mm_stack_t* const mm_stack);
// [Accessors]
GEM_INLINE approximate_search_parameters_t* archive_search_get_search_parameters(archive_search_t* const archive_search);
GEM_INLINE sequence_t* archive_search_get_sequence(const archive_search_t* const archive_search);
GEM_INLINE matches_t* archive_search_get_matches(const archive_search_t* const archive_search);
GEM_INLINE uint64_t archive_search_get_num_potential_canditates(const archive_search_t* const archive_search);

/*
 * SingleEnd Indexed Search (SE Online Approximate String Search)
 */
GEM_INLINE void archive_search_prepare_sequence(archive_search_t* const archive_search);
GEM_INLINE void archive_search_single_end(archive_search_t* const archive_search);

///*
// * PE Pairing (based on SE matches)
// */
//GEM_INLINE void archive_extend__pair_matches(
//    const archive_t* const archive,multimatches* const multimatches,
//    const uint64_t max_stratum,const uint64_t id_matches_from,const uint64_t id_matches_to,
//    fmi_extend_parameters* const extend_parameters,vector_pool* const mpool);

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
 * Error Messages
 */
//#define GEM_ERROR_ARCHIVE_SEARCH_

#endif /* ARCHIVE_SEARCH_H_ */
