/*
 * PROJECT: GEMMapper
 * FILE: graph_text.h
 * DATE: 01/02/2014
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#ifndef GRAPH_TEXT_H_
#define GRAPH_TEXT_H_

#include "essentials.h"

/*
 * Checkers
 */
#define GRAPH_TEXT_CHECK(graph_text) \
  GEM_CHECK_NULL(graph_text)

/*
 * Graph Text
 */
typedef uint64_t edge_attributes_t; // Edge Attributes (Overlapping/Non-Overlapping;SNV/Source/Destiny;...)
typedef struct {
  uint64_t position_src;        // Source Position
  uint64_t position_dst;        // Destiny Position
  edge_attributes_t attributes; // Attributes
} edge_t;
typedef struct {
  /* ID & type */
  uint64_t non_overlapping_id;       // ID (Only accounting non-overlapping vertices)
  edge_attributes_t edge_attributes; // Aggregated attributes of all IN/OUT edges
  /* Location & edges */
  uint64_t position;                 // Index position of the vertex
  uint64_t edge_begin_idx;           // Begin index of the IN/OUT edges in the edges-table
  uint64_t edge_end_idx;             // End index of the IN/OUT edges in the edges-table (not included)
} vertex_t;
typedef struct {
  /* Edges */
  uint64_t num_edges;     // Total number of edges
  uint64_t num_vertices;  // Total number of vertices
  vertex_t* vertices;     // Vertices
  edge_t* edges;          // Edges
  /* Fast-Lookup */
  ihash_t* vertex_hash;   // Vertices indexed by position_src (vertex_t*)
  /* MM */
  mm_t* mm;
} graph_text_t;

/*
 * Setup
 */
graph_text_t* graph_text_read(fm_t* const file_manager);
graph_text_t* graph_text_read_mem(mm_t* const memory_manager);
void graph_text_delete(graph_text_t* const graph);

/*
 * Accessors
 */
uint64_t graph_text_get_size(graph_text_t* const graph);
vertex_t* graph_text_get_vertex(graph_text_t* const graph,const uint64_t vertex_index);
edge_t* graph_text_get_edge(graph_text_t* const graph,const uint64_t edge_index);

/*
 * Lookup vertex
 */
uint64_t graph_text_lookup_vertex_index(graph_text_t* const graph,const uint64_t index_position);
uint64_t graph_text_lookup_previous_vertex_index(
    graph_text_t* const graph,const uint64_t index_position,const uint64_t max_distance);
uint64_t graph_text_lookup_next_vertex_index(
    graph_text_t* const graph,const uint64_t index_position,const uint64_t max_distance);

/*
 * Edge Attributes
 *   An edge can be overlapping (necessarily an SNV) or non-overlapping (SNV or jump)
 *   Classification
 *   - TYPE I (SNV Overlapping)
 *       Overlaps on character of the reference {ACGTE} over {ACGTN}
 *   - TYPE II (SNV Non-Overlapping)
 *       Creates a one-lasso at the position. {ACGTE} over {E}
 *   - TYPE III (General edge)
 *       General edge to another position
 */
#define EDGE_ATTRIBUTES_NULL       (0ull)

#define EDGE_SVN_A                 (1ull)
#define EDGE_SVN_C                 (1ull<<1)
#define EDGE_SVN_G                 (1ull<<2)
#define EDGE_SVN_T                 (1ull<<3)
#define EDGE_SVN_N                 (1ull<<4)
#define EDGE_SVN_DEL               (1ull<<5)

#define EDGE_REFERENCE_A           (1ull<<6)
#define EDGE_REFERENCE_C           (1ull<<7)
#define EDGE_REFERENCE_G           (1ull<<8)
#define EDGE_REFERENCE_T           (1ull<<9)
#define EDGE_REFERENCE_N           (1ull<<10)

#define EDGE_OVERLAPPING           (1ull<<11)
#define EDGE_SOURCE                (1ull<<12)
#define EDGE_DESTINY               (1ull<<13)

edge_attributes_t edge_attributes_snv_new(const char snv_character);
edge_attributes_t edge_attributes_snv_add(const edge_attributes_t attributes,const char snv_character);

edge_attributes_t edge_attributes_set_overlapping(const edge_attributes_t attributes);
edge_attributes_t edge_attributes_snv_set_reference(const edge_attributes_t attributes,const char reference_character);
char edge_attributes_snv_get_reference(const edge_attributes_t attributes);

bool edge_attributes_is_jump(const edge_attributes_t attributes);
bool edge_attributes_is_jump_source(const edge_attributes_t attributes);
bool edge_attributes_is_jump_destiny(const edge_attributes_t attributes);
bool edge_attributes_is_svn(const edge_attributes_t attributes);
bool edge_attributes_is_overlapping(const edge_attributes_t attributes);

/*
 * Display
 */
void graph_text_print(FILE* const stream,graph_text_t* const graph,const bool full_dump);

void edge_print(FILE* const stream,const edge_t* const edge);
void edge_attributes_print(FILE* const stream,const edge_attributes_t attributes);

/*
 * Error Messages
 */
#define GEM_ERROR_GRAPH_TEXT_EDGE_GROUP_OOB "Graph-Text. Requested edge-group index (%"PRIu64") out-of-bounds [0,%"PRIu64")]"
#define GEM_ERROR_GRAPH_TEXT_EDGE_OOB "Graph-Text. Requested edge index (%"PRIu64") out-of-bounds [0,%"PRIu64")]"

#endif /* GRAPH_TEXT_H_ */
