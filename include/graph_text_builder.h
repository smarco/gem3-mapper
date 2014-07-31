/*
 * PROJECT: GEMMapper
 * FILE: graph_text_builder.h
 * DATE: 01/02/2014
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#ifndef GRAPH_TEXT_BUILDER_H_
#define GRAPH_TEXT_BUILDER_H_

#include "essentials.h"
#include "locator.h"
#include "graph_text.h"

/*
 * Constants
 */
#define GRAPH_TEXT_BUILDER_POSITION_UNKNOWN UINT64_MAX

/*
 * Checkers
 */
#define GRAPH_TEXT_BUILDER_CHECK(graph_text_builder) GEM_CHECK_NULL(graph_text_builder)

/*
 * Graph Text Builder
 */
typedef struct _graph_link_t graph_link_t; // Forward declaration of graph_link_t
struct _graph_link_t {
  edge_attributes_t attributes; // Type
  int64_t tag_id;               // Tag ID (Translated from sequence string) [Negative if Reverse/Complement]
  uint64_t position_text;       // Position in the text
  uint64_t position_index;      // Position in the index
  graph_link_t* link;           // Link {to,from}
};
typedef struct {
  uint64_t offset;
  uint64_t num_links;
} graph_links_index_t;
typedef struct {
  /* Graph links */
  svector_t* graph_links;                       // (graph_link_t)
  svector_iterator_t graph_links_iterator;      // Graph-links Writing Iterator
  /* Sorted Graph links */
  graph_link_t** link_table_sorted_mem;         // pointer to Memory Allocated
  graph_link_t** link_table_sorted;             // Sorted link-table wrt (tagID,positionText)
  graph_link_t link_table_sentinel;             // Dummy link at the end/beginning of the table
  uint64_t link_table_length;                   // Total number of links after processing
  graph_links_index_t* graph_links_idx_forward; // Index starting position for a positive/forward Tag-ID
  graph_links_index_t* graph_links_idx_reverse; // Index starting position for a negative/reverse Tag-ID
  uint64_t graph_links_idx_length;              // Total length of each graph-link-idx table
  /* MM */
  mm_slab_t* mm_slab;
} graph_text_builder_t;
typedef struct {
  /* Iterator */
  strand_t strand;                 // Traversal direction {forward,reverse}
  int64_t tag_id;                  // Tag-ID of the links
  graph_link_t** next_link;        // Pointer to next link (links-sorted table by text-position)
  bool eoi;                        // No more link from this sequence
  /* Link-table chunk */
  graph_link_t** link_table_chunk; // Link-table chunk (sharing same tagID & text-position)
  uint64_t num_links;              // Number of links (sharing same tagID & text-position)
  bool has_overlapping_jump;       // Link-table has an overlapping jump (overlapping SNV)
  bool has_non_overlapping_jump;   // Link-table has a general jump (non-overlapping SNV or general jump)
} graph_sorted_link_locator_t;

/*
 * Setup
 */
GEM_INLINE graph_text_builder_t* graph_text_builder_new(mm_slab_t* const mm_slab);
GEM_INLINE void graph_text_builder_delete(graph_text_builder_t* const graph_builder);

/*
 * Step-0:: Adding links {Jumps/SNVs}
 */
GEM_INLINE uint64_t graph_text_builder_get_num_links(graph_text_builder_t* const graph_builder);

GEM_INLINE void graph_text_builder_add_general_link(
    graph_text_builder_t* const graph_builder,
    graph_link_t* const graph_link_from,graph_link_t* const graph_link_to);
GEM_INLINE void graph_text_builder_add_snv_link(
    graph_text_builder_t* const graph_builder,
    graph_link_t* const graph_link,const char character,const bool overlaps_reference);

/*
 * Step-1:: Sort all links (as to solve pending text-positions => index-positions)
 */
GEM_INLINE void graph_text_builder_link_table_sort(graph_text_builder_t* const graph_builder,const uint64_t total_tags);
GEM_INLINE uint64_t graph_text_builder_link_table_get_length(graph_text_builder_t* const graph_builder,const int64_t tag_id);

/*
 * Step-2:: Locate links (as to solve pending index-positions)
 */
GEM_INLINE void graph_text_builder_link_locator_iterate_forward(
    graph_text_builder_t* const graph_builder,
    graph_sorted_link_locator_t* const link_locator,const int64_t tag_id);
GEM_INLINE void graph_text_builder_link_locator_iterate_backward(
    graph_text_builder_t* const graph_builder,
    graph_sorted_link_locator_t* const link_locator,const int64_t tag_id);
GEM_INLINE bool graph_text_builder_link_locator_find(
    graph_text_builder_t* const graph_builder,
    graph_sorted_link_locator_t* const link_locator,const int64_t tag_id,const uint64_t text_position);
GEM_INLINE bool graph_text_builder_link_locator_cmp(graph_sorted_link_locator_t* const link_locator,const uint64_t text_position);

GEM_INLINE bool graph_text_builder_link_locator_has_jump_non_overlapping(graph_sorted_link_locator_t* const link_locator);
GEM_INLINE bool graph_text_builder_link_locator_has_jump_overlapping(graph_sorted_link_locator_t* const link_locator);
GEM_INLINE char graph_text_builder_link_locator_jump_overlapping_get_reference_char(graph_sorted_link_locator_t* const link_locator);
GEM_INLINE void graph_text_builder_link_locator_solve_jump_non_overlapping(
    graph_sorted_link_locator_t* const link_locator,const uint64_t index_position);
GEM_INLINE void graph_text_builder_link_locator_solve_jump_overlapping(
    graph_sorted_link_locator_t* const link_locator,const uint64_t index_position,const char reference_char);

/*
 * Step-3:: Write
 */
GEM_INLINE void graph_text_builder_write(
    fm_t* const file_manager,graph_text_builder_t* const graph_builder,locator_builder_t* const locator);

/*
 * Display
 */
GEM_INLINE void graph_text_builder_print(
    FILE* const stream,graph_text_builder_t* const graph_builder,
    locator_builder_t* const locator,const bool display_links);
GEM_INLINE void graph_text_builder_link_table_print(
    FILE* const stream,graph_text_builder_t* const graph_builder,
    locator_builder_t* const locator,const bool display_links);

/*
 * Errors
 */
#define GEM_ERROR_GRAPH_TEXT_BUILDER_POSITION_UNSOLVED "Building Graph-Text. Position unsolved (%s:%c:%lu)"
#define GEM_ERROR_GRAPH_TEXT_BUILDER_SNV_UNKNOWN_REFERENCE_CHAR "Building Graph-Text. Unknown reference char"

#endif /* GRAPH_TEXT_BUILDER_H_ */
