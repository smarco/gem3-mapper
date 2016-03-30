/*
 * PROJECT: GEMMapper
 * FILE: archive_search.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef ARCHIVE_SEARCH_H_
#define ARCHIVE_SEARCH_H_

#include "utils/essentials.h"
#include "archive/archive.h"
#include "approximate_search/approximate_search.h"
#include "approximate_search/mm_search.h"
#include "data_structures/sequence.h"
#include "matches/matches_classify.h"

/*
 * Archive Search State
 */
typedef enum {
  archive_search_pe_begin = 0,                     // Beginning of the search
  archive_search_pe_search_end1 = 1,               // Generate candidates for end1
  archive_search_pe_search_end2 = 2,               // Generate candidates for end2
  archive_search_pe_recovery = 3,                  // Recover by extension when needed
  archive_search_pe_find_pairs = 4,                // Cross-link matches from both ends
  archive_search_pe_end = 5                        // End of the current workflow
} archive_search_pe_state_t;
extern const char* archive_search_pe_state_label[7];

/*
 * Archive Search
 */
typedef struct {
  /* Archive */
  archive_t* archive;                        // Archive
  /* Archive Paired-End-Search (Only end/1 used in PE search) */
  archive_search_pe_state_t pe_search_state; // Search State
  bool pair_searched;                        // Paired search performed
  bool pair_extended;                        // Paired extension performed
  bool pair_extended_shortcut;               // Paired extension performed (to shortcut)
  /* Sequence */
  sequence_t sequence;                       // Input
  sequence_t rc_sequence;                    // Generated
  /* Parameters */
  search_parameters_t search_parameters;     // Search parameters
  bool emulate_rc_search;                    // Flow control
  bool probe_strand;                         // Flow control
  bool buffered_search;                      // Buffered Search
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
void archive_search_init(
    archive_search_t* const archive_search,
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    const bool buffered_search,
    mm_stack_t* const mm_stack);
void archive_search_destroy(archive_search_t* const archive_search);

void archive_search_se_new(
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    const bool buffered_search,
    mm_stack_t* const mm_stack,
    archive_search_t** const archive_search);
void archive_search_pe_new(
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    const bool buffered_search,
    mm_stack_t* const mm_stack,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2);
void archive_search_reset(archive_search_t* const archive_search);
void archive_search_delete(archive_search_t* const archive_search);
/*
 * Memory Injection (Support Data Structures)
 */
void archive_search_inject_mm_stack(
    archive_search_t* const archive_search,
    mm_stack_t* const mm_stack);
void archive_search_inject_mapper_stats(
    archive_search_t* const archive_search,
    mapper_stats_t* mapper_stats);
void archive_search_inject_interval_set(
    archive_search_t* const archive_search,
    interval_set_t* const interval_set);
void archive_search_inject_text_collection(
    archive_search_t* const archive_search,
    text_collection_t* const text_collection);
void archive_search_inject_filtering_candidates(
    archive_search_t* const archive_search,
    filtering_candidates_t* const filtering_candidates_forward,
    filtering_candidates_t* const filtering_candidates_reverse,
    text_collection_t* const text_collection,
    mm_stack_t* const mm_stack);

/*
 * Accessors
 */
sequence_t* archive_search_get_sequence(const archive_search_t* const archive_search);
bool archive_search_finished(const archive_search_t* const archive_search);

uint64_t archive_search_get_max_region_length(const archive_search_t* const archive_search);
uint64_t archive_search_get_num_zero_regions(const archive_search_t* const archive_search);

uint64_t archive_search_get_num_regions_profile(const archive_search_t* const archive_search);
uint64_t archive_search_get_num_decode_candidates(const archive_search_t* const archive_search);
uint64_t archive_search_get_num_verify_candidates(const archive_search_t* const archive_search);

/*
 * Errors
 */
#define GEM_ERROR_ARCHIVE_SEARCH_INDEX_COMPLEMENT_REQUIRED "Archive Search. Explicit indexed complement required"

#endif /* ARCHIVE_SEARCH_H_ */
