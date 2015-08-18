/*
 * PROJECT: GEMMapper
 * FILE: archive_search.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef ARCHIVE_SEARCH_H_
#define ARCHIVE_SEARCH_H_

#include "essentials.h"

#include "archive.h"
#include "select_parameters.h"
#include "approximate_search.h"
#include "sequence.h"
#include "matches_classify.h"
#include "mm_search.h"

/*
 * Checker
 */
#define ARCHIVE_SEARCH_CHECK(archive_search) \
  GEM_CHECK_NULL(archive_search); \
  ARCHIVE_CHECK(archive_search->archive)

/*
 * Archive Search State
 */
typedef enum {
  archive_search_pe_begin,                     // Beginning of the search
  archive_search_pe_search_end1,               // Generate candidates for end1
  archive_search_pe_search_end2,               // Generate candidates for end2
  archive_search_pe_recovery,                  // Recover by extension when needed
  archive_search_pe_find_pairs,                // Cross-link matches from both ends
  archive_search_pe_end                        // End of the current workflow
} archive_search_pe_state_t;

/*
 * Archive Search
 */
typedef struct {
  /* Archive */
  archive_t* archive;                        // Archive
  /* Archive Paired-End-Search (Only end/1 used in PE search) */
  archive_search_pe_state_t pe_search_state; // Search State
  bool pair_searched;                        // Paired search performed
  bool pair_extended;                        // Paired estension performed
  matches_class_t end_class;
  /* Sequence */
  sequence_t sequence;                       // Input
  sequence_t rc_sequence;                    // Generated
  /* Parameters */
  bool emulate_rc_search;                    // Flow control
  bool probe_strand;                         // Flow control
  as_parameters_t as_parameters;             // Approximated-Search parameters (evaluated to read-length)
  select_parameters_t* select_parameters;    // Select parameters
  /* Approximate Search */
  approximate_search_t forward_search_state; // Forward Search State
  approximate_search_t reverse_search_state; // Reverse Search State
  /* Text-Collection */
  text_collection_t* text_collection;
  /* Stats */
  mapper_stats_t* mapper_stats;              // Mapping statistics
  /* MM */
  mm_stack_t* mm_stack;                      // MM-Stack
} archive_search_t;

/*
 * Setup
 */
GEM_INLINE archive_search_t* archive_search_new(
    archive_t* const archive,search_parameters_t* const search_parameters,
    select_parameters_t* const select_parameters);
GEM_INLINE void archive_search_configure(
    archive_search_t* const archive_search,const sequence_end_t sequence_end,
    mm_search_t* const mm_search);
GEM_INLINE void archive_search_reset(archive_search_t* const archive_search);
GEM_INLINE void archive_search_delete(archive_search_t* const archive_search);

/*
 * Accessors
 */
GEM_INLINE sequence_t* archive_search_get_sequence(const archive_search_t* const archive_search);
GEM_INLINE uint64_t archive_search_get_search_canditates(const archive_search_t* const archive_search);
GEM_INLINE uint64_t archive_search_get_search_exact_matches(const archive_search_t* const archive_search);
GEM_INLINE uint64_t archive_search_get_max_region_length(const archive_search_t* const archive_search);
GEM_INLINE uint64_t archive_search_get_num_zero_regions(const archive_search_t* const archive_search);

GEM_INLINE bool archive_search_finished(const archive_search_t* const archive_search);

/*
 * Utils
 */
GEM_INLINE void archive_search_hold_verification_candidates(archive_search_t* const archive_search);
GEM_INLINE void archive_search_release_verification_candidates(archive_search_t* const archive_search);

#endif /* ARCHIVE_SEARCH_H_ */
