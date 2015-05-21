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

#include "sequence.h"
#include "matches.h"
#include "paired_matches.h"

#include "archive.h"
#include "archive_select_parameters.h"
#include "approximate_search.h"
#include "mm_search.h"

/*
 * Checker
 */
#define ARCHIVE_SEARCH_CHECK(archive_search) \
  GEM_CHECK_NULL(archive_search); \
  ARCHIVE_CHECK(archive_search->archive)

/*
 * Archive Approximate-Search
 */
typedef enum {
  archive_search_pe_begin,                     // Beginning of the search
  archive_search_pe_search_end1,               // Generate candidates for end1
  archive_search_pe_search_end2,               // Generate candidates for end2
  archive_search_pe_extend_end1,               // Try to extend end1
  archive_search_pe_extended_end1,             // End1 extended
  archive_search_pe_extend_end2,               // Try to extend end2
  archive_search_pe_extended_end2,             // End2 extended
  archive_search_pe_paired_filtering_discard,  // Discard paired-filtering candidates
  archive_search_pe_paired_filtering_verify,   // Verify paired-filtering candidates
  archive_search_pe_paired_filtering_verified, // Paired-filtering candidates had been verified
  archive_search_pe_both_ends_verify,          // Verify candidates for both ends
  archive_search_pe_both_ends_verified,        // Candidates of both ends had been verified
  archive_search_pe_recovery,                  // Paired-end recovery by extension
  archive_search_pe_end                        // End of the current workflow
} archive_search_pe_state_t;
typedef struct {
  /* Archive */
  archive_t* archive;                        // Archive
  /* Archive Paired-End-Search (Only end/1 used in PE search) */
  archive_search_pe_state_t pe_search_state; // Search State
  bool try_extending_end;                    // Tried extension
  bool paired_filtering;                     // Used paired-filtering
  /* Sequence */
  sequence_t sequence;                       // Input
  sequence_t rc_sequence;                    // Generated
  /* Parameters */
  bool emulate_rc_search;                    // Flow control
  bool probe_strand;                         // Flow control
  search_actual_parameters_t search_actual_parameters; // Search parameters (evaluated to read-length)
  select_parameters_t* select_parameters;              // Select parameters
  /* Approximate Search */
  approximate_search_t forward_search_state; // Forward Search State
  approximate_search_t reverse_search_state; // Reverse Search State
  /* Text-Collection */
  text_collection_t* text_collection;
  /* Stats */
  mapper_stats_t* mapper_stats;            // Mapping statistics
  /* MM */
  mm_stack_t* mm_stack;                      // MM-Stack
} archive_search_t;

/*
 * Archive Search Setup
 */
GEM_INLINE archive_search_t* archive_search_new(
    archive_t* const archive,search_parameters_t* const search_parameters,
    select_parameters_t* const select_parameters);
GEM_INLINE void archive_search_configure(archive_search_t* const archive_search,mm_search_t* const mm_search);
GEM_INLINE void archive_search_pe_configure(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    mm_search_t* const mm_search);
GEM_INLINE void archive_search_reset(archive_search_t* const archive_search);
GEM_INLINE void archive_search_delete(archive_search_t* const archive_search);
// [Accessors]
GEM_INLINE sequence_t* archive_search_get_sequence(const archive_search_t* const archive_search);
GEM_INLINE uint64_t archive_search_get_search_canditates(const archive_search_t* const archive_search);

/*
 * SingleEnd Indexed Search
 */
// Step-wise
GEM_INLINE void archive_search_generate_candidates(archive_search_t* const archive_search);
GEM_INLINE void archive_search_verify_candidates(archive_search_t* const archive_search,matches_t* const matches);
GEM_INLINE void archive_search_finish_search(archive_search_t* const archive_search,matches_t* const matches);
GEM_INLINE void archive_search_copy_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer);
GEM_INLINE void archive_search_retrieve_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer,matches_t* const matches);
// SE Online Approximate String Search
GEM_INLINE void archive_search_single_end(archive_search_t* const archive_search,matches_t* const matches);

/*
 * PairedEnd Indexed Search
 */
// Step-wise
GEM_INLINE void archive_search_pe_generate_candidates(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);
GEM_INLINE void archive_search_pe_finish_search(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);
// PE Online Approximate String Search
GEM_INLINE void archive_search_paired_end(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

/*
 * Errors
 */
#define GEM_ERROR_ARCHIVE_SEARCH_INDEX_COMPLEMENT_REQUIRED "Archive Search. Explicit indexed complement required"

#endif /* ARCHIVE_SEARCH_H_ */
