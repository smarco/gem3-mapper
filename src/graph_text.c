/*
 * PROJECT: GEMMapper
 * FILE: graph_edges.c
 * DATE: 01/02/2014
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#include "graph_text.h"
#include "dna_string.h"

/*
 * Setup
 */
GEM_INLINE void graph_text_initialize(graph_text_t* const graph) {
  GRAPH_TEXT_CHECK(graph);
  // Initialize Hash Vertex Positions
  graph->vertex_hash = ihash_new();
  uint64_t i;
  for (i=0;i<graph->num_vertices;++i) {
    ihash_insert(graph->vertex_hash,graph->vertices[i].position,(graph->vertices+i));
  }
}
GEM_INLINE graph_text_t* graph_text_read(fm_t* const file_manager) {
  // Allocate
  graph_text_t* const graph = mm_alloc(graph_text_t);
  // Read Meta-Data
  graph->num_vertices = fm_read_uint64(file_manager);
  graph->num_edges = fm_read_uint64(file_manager);
  // Read Vertices+Edges
  const uint64_t vertices_size = graph->num_vertices*sizeof(vertex_t);
  const uint64_t edge_size = graph->num_edges*sizeof(edge_t);
  graph->mm = fm_load_mem(file_manager,vertices_size+edge_size);
  graph->vertices = mm_read_mem(graph->mm,vertices_size);
  graph->edges = mm_read_mem(graph->mm,edge_size);
  // Initialize
  graph_text_initialize(graph);
  // Return
  return graph;
}
GEM_INLINE graph_text_t* graph_text_read_mem(mm_t* const memory_manager) {
  // Allocate
  graph_text_t* const graph = mm_alloc(graph_text_t);
  // Read Meta-Data
  graph->num_vertices = mm_read_uint64(memory_manager);
  graph->num_edges = mm_read_uint64(memory_manager);
  // Read Edges
  const uint64_t vertices_size = graph->num_vertices*sizeof(vertex_t);
  const uint64_t edge_size = graph->num_edges*sizeof(edge_t);
  graph->mm = NULL;
  graph->vertices = mm_read_mem(memory_manager,vertices_size);
  graph->edges = mm_read_mem(memory_manager,edge_size);
  // Initialize
  graph_text_initialize(graph);
  // Return
  return graph;
}
GEM_INLINE void graph_text_delete(graph_text_t* const graph) {
  GRAPH_TEXT_CHECK(graph);
  // Delete hash
  ihash_delete(graph->vertex_hash);
  // Free Mem & Handler
  if (graph->mm) mm_bulk_free(graph->mm);
  mm_free(graph);
}
/*
 * Accessors
 */
GEM_INLINE uint64_t graph_text_get_size(graph_text_t* const graph) {
  GRAPH_TEXT_CHECK(graph);
  return graph->num_edges*sizeof(edge_t) + ihash_get_size(graph->vertex_hash);
}
GEM_INLINE vertex_t* graph_text_get_vertex(graph_text_t* const graph,const uint64_t vertex_index) {
  GRAPH_TEXT_CHECK(graph);
  gem_fatal_check(vertex_index>=graph->num_vertices,GRAPH_TEXT_EDGE_GROUP_OOB,vertex_index,graph->num_vertices);
  return graph->vertices+vertex_index;
}
GEM_INLINE edge_t* graph_text_get_edge(graph_text_t* const graph,const uint64_t edge_index) {
  GRAPH_TEXT_CHECK(graph);
  gem_fatal_check(edge_index>=graph->num_edges,GRAPH_TEXT_EDGE_OOB,edge_index,graph->num_edges);
  return graph->edges+edge_index;
}
/*
 * Lookup vertex
 */
GEM_INLINE uint64_t graph_text_lookup_vertex_index(graph_text_t* const graph,const uint64_t index_position) {
  vertex_t* const vertex = ihash_get(graph->vertex_hash,index_position,vertex_t);
  if (gem_expect_false(vertex==NULL)) return UINT64_MAX;
  return vertex-graph->vertices;
}
GEM_INLINE uint64_t graph_text_lookup_vertex(graph_text_t* const graph,const uint64_t index_position) {
  GRAPH_TEXT_CHECK(graph);
  /*
   * Binary search of the edge. Returns the index of the group which
   *   position is closest to @index_position from the left
   */
  const vertex_t* const vertices = graph->vertices;
  uint64_t lo = 0;
  uint64_t hi = graph->num_edges-1;
  do {
    const uint64_t half = (hi+lo)/2;
    if (index_position > vertices[half].position) {
      lo = half+1;
    } else {
      hi = half;
    }
  } while (hi > lo);
  return lo;
}
GEM_INLINE uint64_t graph_text_lookup_previous_vertex_index(
    graph_text_t* const graph,const uint64_t index_position,const uint64_t max_distance) {
  GRAPH_TEXT_CHECK(graph);
  // Binary search of the edge
  int64_t sentinel = graph_text_lookup_vertex(graph,index_position);
  // Find Nearest Jump Edge (to the right)
  const vertex_t* const vertices = graph->vertices;
  while (sentinel >= 0) {
    if ((int64_t)index_position-(int64_t)vertices[sentinel].position > (int64_t)max_distance) return UINT64_MAX; // No Edge-Group
    if (edge_attributes_is_jump_source(vertices[sentinel].edge_attributes)) return sentinel;
    --(sentinel);
  }
  return UINT64_MAX; // No Edge-Group
}
GEM_INLINE uint64_t graph_text_lookup_next_vertex_index(
    graph_text_t* const graph,const uint64_t index_position,const uint64_t max_distance) {
  GRAPH_TEXT_CHECK(graph);
  // Binary search of the edge
  int64_t sentinel = graph_text_lookup_vertex(graph,index_position);
  // Find Nearest Jump Edge (to the right)
  const vertex_t* const vertices = graph->vertices;
  while (sentinel < graph->num_vertices) {
    if ((int64_t)vertices[sentinel].position-(int64_t)index_position > (int64_t)max_distance) return UINT64_MAX; // No Edge-Group
    if (edge_attributes_is_jump_source(vertices[sentinel].edge_attributes)) return sentinel;
    ++(sentinel);
  }
  return UINT64_MAX; // No Edge-Group
}
/*
 * Edge Attributes
 */
#define edge_attributes_snv_encode(character) (edge_attributes_snv_encode_table[(int)(character)])
const edge_attributes_t edge_attributes_snv_encode_table[256] =
{
    [0 ... 255] = 0,
    [0] = EDGE_SVN_DEL,
    [DNA_CHAR_A] = EDGE_SVN_A, [DNA_CHAR_C] = EDGE_SVN_C,
    [DNA_CHAR_G] = EDGE_SVN_G, [DNA_CHAR_T] = EDGE_SVN_T,
    [DNA_CHAR_N] = EDGE_SVN_N
};
#define edge_attributes_reference_encode(character) (edge_attributes_reference_encode_table[(int)(character)])
const edge_attributes_t edge_attributes_reference_encode_table[256] =
{
    [0 ... 255] = 0,
    [DNA_CHAR_A] = EDGE_REFERENCE_A, [DNA_CHAR_C] = EDGE_REFERENCE_C,
    [DNA_CHAR_G] = EDGE_REFERENCE_G, [DNA_CHAR_T] = EDGE_REFERENCE_T,
    [DNA_CHAR_N] = EDGE_REFERENCE_N
};
GEM_INLINE edge_attributes_t edge_attributes_snv_new(const char snv_character) {
  return edge_attributes_snv_encode(snv_character);
}
GEM_INLINE edge_attributes_t edge_attributes_snv_add(const edge_attributes_t attributes,const char snv_character) {
  return attributes | edge_attributes_snv_encode(snv_character);
}
GEM_INLINE edge_attributes_t edge_attributes_set_overlapping(const edge_attributes_t attributes) {
  return attributes | EDGE_OVERLAPPING;
}
GEM_INLINE edge_attributes_t edge_attributes_snv_set_reference(const edge_attributes_t attributes,const char reference_character) {
  return attributes | edge_attributes_reference_encode(reference_character);
}
GEM_INLINE char edge_attributes_snv_get_reference(const edge_attributes_t attributes) {
  if (attributes & (EDGE_REFERENCE_A | EDGE_REFERENCE_C |
                   EDGE_REFERENCE_G | EDGE_REFERENCE_T |
                   EDGE_REFERENCE_N)) {
    if (attributes & EDGE_REFERENCE_A) return DNA_CHAR_A;
    if (attributes & EDGE_REFERENCE_C) return DNA_CHAR_C;
    if (attributes & EDGE_REFERENCE_G) return DNA_CHAR_G;
    if (attributes & EDGE_REFERENCE_T) return DNA_CHAR_T;
    if (attributes & EDGE_REFERENCE_N) return DNA_CHAR_N;
  }
  return 0;
}
GEM_INLINE bool edge_attributes_is_jump(const edge_attributes_t attributes) {
  return edge_attributes_is_jump_source(attributes) || edge_attributes_is_jump_destiny(attributes);
}
GEM_INLINE bool edge_attributes_is_jump_source(const edge_attributes_t attributes) {
  return attributes == EDGE_SOURCE;
}
GEM_INLINE bool edge_attributes_is_jump_destiny(const edge_attributes_t attributes) {
  return attributes == EDGE_DESTINY;
}
GEM_INLINE bool edge_attributes_is_svn(const edge_attributes_t attributes) {
  return attributes &
      (EDGE_SVN_A | EDGE_SVN_C | EDGE_SVN_G |
       EDGE_SVN_T | EDGE_SVN_N | EDGE_SVN_DEL |
       EDGE_REFERENCE_A | EDGE_REFERENCE_C |
       EDGE_REFERENCE_G | EDGE_REFERENCE_T |
       EDGE_REFERENCE_N | EDGE_OVERLAPPING);
}
GEM_INLINE bool edge_attributes_is_overlapping(const edge_attributes_t attributes) {
  return attributes & EDGE_OVERLAPPING;
}
/*
 * Display
 */
GEM_INLINE void graph_text_print(FILE* const stream,graph_text_t* const graph,const bool full_dump) {
  // Calculate magnitudes
  const uint64_t edge_table_size = graph->num_edges*sizeof(edge_t);
  const uint64_t edge_hash_size = ihash_get_size(graph->vertex_hash);
  const uint64_t total_size = edge_table_size + edge_hash_size;
  // Display graph info
  tab_fprintf(stream,"[GEM]>Graph.Text\n");
  tab_fprintf(stream,"  => Edges %lu\n",graph->num_edges);
  tab_fprintf(stream,"    => Edges.Size %lu MB (%2.3f%%)\n",CONVERT_B_TO_MB(edge_table_size),PERCENTAGE(edge_table_size,total_size));
  tab_fprintf(stream,"    => Edges.Hash.Size %lu MB (%2.3f%%)\n",CONVERT_B_TO_MB(edge_hash_size),PERCENTAGE(edge_hash_size,total_size));
  // Dump
  if (full_dump) {
    tab_fprintf(stream,"  => Edges.Groups\n");
    uint64_t i;
    for (i=0;i<graph->num_vertices;++i) {
      tab_fprintf(stream,"    #%ld [%lu,%lu)\t{",
          graph->vertices[i].non_overlapping_id,
          graph->vertices[i].edge_begin_idx,graph->vertices[i].edge_end_idx);
      edge_attributes_print(stream,graph->vertices[i].edge_attributes);
      fprintf(stream,"}\n");
    }
  }
  if (full_dump) {
    tab_fprintf(stream,"  => Edges.List\n");
    uint64_t i;
    for (i=0;i<graph->num_edges;++i) {
      tab_fprintf(stream,"    [%lu]\t",i);
      edge_print(stream,graph->edges+i);
      fprintf(stream,"\n");
    }
  }
  // Flush
  fflush(stream);
}
GEM_INLINE void edge_print(FILE* const stream,const edge_t* const edge) {
  if (edge_attributes_is_svn(edge->attributes)) { // SNVs
    fprintf(stream,"{"); edge_attributes_print(stream,edge->attributes); fprintf(stream,"}");
    fprintf(stream,"(idx:%lu)",edge->position_src);
  } else { // General Link
    if (edge->attributes==EDGE_SOURCE) {
      fprintf(stream,"{Source}\t(idx:%lu)->(idx:%lu)",edge->position_src,edge->position_dst);
    } else {
      fprintf(stream,"{Destiny}\t(idx:%lu)<-(idx:%lu)",edge->position_src,edge->position_dst);
    }
  }
}
GEM_INLINE void edge_attributes_print(FILE* const stream,const edge_attributes_t attributes) {
  if (attributes==EDGE_ATTRIBUTES_NULL) {
    fprintf(stream,"Null");
    fflush(stream);
    return;
  }
  if (attributes==EDGE_SOURCE) {
    fprintf(stream,"Source");
    fflush(stream);
  }
  if (attributes==EDGE_DESTINY) {
    fprintf(stream,"Destiny");
    fflush(stream);
  }
  fprintf(stream,"SNV[");
  if (attributes & EDGE_SVN_A)   fprintf(stream,"a");
  if (attributes & EDGE_SVN_C)   fprintf(stream,"c");
  if (attributes & EDGE_SVN_G)   fprintf(stream,"g");
  if (attributes & EDGE_SVN_T)   fprintf(stream,"t");
  if (attributes & EDGE_SVN_N)   fprintf(stream,"n");
  if (attributes & EDGE_SVN_DEL) fprintf(stream,"E");
  if (attributes & EDGE_REFERENCE_A) fprintf(stream,"A");
  if (attributes & EDGE_REFERENCE_C) fprintf(stream,"C");
  if (attributes & EDGE_REFERENCE_G) fprintf(stream,"G");
  if (attributes & EDGE_REFERENCE_T) fprintf(stream,"T");
  if (attributes & EDGE_REFERENCE_N) fprintf(stream,"N");
  if (attributes & EDGE_OVERLAPPING) {
    fprintf(stream,"]+");
  } else {
    fprintf(stream,"]");
  }
  // Flush
  fflush(stream);
}
