/*
 * PROJECT: GEMMapper
 * FILE: graph_text_builder.c
 * DATE: 01/02/2014
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#include "graph_text_builder.h"

/*
 * Builder
 */
GEM_INLINE graph_text_builder_t* graph_text_builder_new(mm_slab_t* const mm_slab) {
  MM_SLAB_CHECK(mm_slab);
  // Allocate handler
  graph_text_builder_t* const graph_builder = mm_alloc(graph_text_builder_t);
  // Graph links
  graph_builder->graph_links = svector_new(mm_slab,graph_link_t);
  svector_iterator_new(&graph_builder->graph_links_iterator,
      graph_builder->graph_links,SVECTOR_WRITE_ITERATOR,0);
  // Sorted Graph links
  graph_builder->link_table_sorted = NULL; // Later
  graph_builder->link_table_length = 0;
  graph_builder->graph_links_idx_forward = NULL; // Later
  graph_builder->graph_links_idx_reverse = NULL; // Later
  // Return
  return graph_builder;
}
GEM_INLINE void graph_text_builder_delete(graph_text_builder_t* const graph_builder) {
  GRAPH_TEXT_BUILDER_CHECK(graph_builder);
  svector_delete(graph_builder->graph_links);
  if (graph_builder->link_table_sorted_mem!=NULL) mm_free(graph_builder->link_table_sorted_mem);
  if (graph_builder->graph_links_idx_forward!=NULL) mm_free(graph_builder->graph_links_idx_forward);
  if (graph_builder->graph_links_idx_reverse!=NULL) mm_free(graph_builder->graph_links_idx_reverse);
  mm_free(graph_builder);
}
/*
 * Step-0:: Adding links {Jumps/SNVs}
 */
GEM_INLINE uint64_t graph_text_builder_get_num_links(graph_text_builder_t* const graph_builder) {
  GRAPH_TEXT_BUILDER_CHECK(graph_builder);
  return svector_get_used(graph_builder->graph_links);
}
GEM_INLINE void graph_text_builder_add_general_link(
    graph_text_builder_t* const graph_builder,
    graph_link_t* const graph_link_from,graph_link_t* const graph_link_to) {
  GRAPH_TEXT_BUILDER_CHECK(graph_builder);
  // Add Forward link
  graph_link_t* const link_from = svector_iterator_get_element(&graph_builder->graph_links_iterator,graph_link_t);
  *link_from = *graph_link_from; // Copy
  link_from->attributes = EDGE_SOURCE;
  svector_write_iterator_next(&graph_builder->graph_links_iterator); // Next
  // Add Backward link
  graph_link_t* const link_to = svector_iterator_get_element(&graph_builder->graph_links_iterator,graph_link_t);
  *link_to = *graph_link_to; // Copy
  link_to->attributes = EDGE_DESTINY;
  svector_write_iterator_next(&graph_builder->graph_links_iterator); // Next
  // Link each other
  link_from->link = link_to;
  link_to->link = link_from;
}
GEM_INLINE void graph_text_builder_add_snv_link(
    graph_text_builder_t* const graph_builder,
    graph_link_t* const graph_link,const char character,const bool overlaps_reference) {
  GRAPH_TEXT_BUILDER_CHECK(graph_builder);
  // Add SNV link
  graph_link_t* const snv = svector_iterator_get_element(&graph_builder->graph_links_iterator,graph_link_t);
  *snv = *graph_link; // Copy
  // Set SNV type
  snv->attributes = edge_attributes_snv_new(character);
  if (overlaps_reference) {
    snv->attributes = edge_attributes_set_overlapping(snv->attributes);
  } else {
    snv->attributes = edge_attributes_snv_add(snv->attributes,0);
  }
  // Close link
  snv->link = NULL;
  // Next
  svector_write_iterator_next(&graph_builder->graph_links_iterator);
}
/*
 * Step-1:: Sort all links (as to solve pending text-positions => index-positions)
 */
int graph_text_builder_link_table_cmp_text_position(const graph_link_t** const ptr_link_a,const graph_link_t** const ptr_link_b) {
  const graph_link_t* const link_a = *ptr_link_a;
  const graph_link_t* const link_b = *ptr_link_b;
  // Compare Tags
  if (link_a->tag_id < link_b->tag_id) return -1;
  if (link_a->tag_id > link_b->tag_id) return  1;
  // Same tag
  const int64_t diff = ((int64_t)link_a->position_text - (int64_t)link_b->position_text);
  if (diff != 0) return diff;
  // Same position (SNVs first)
  if (link_a->link!=NULL && link_b->link==NULL) return -1;
  if (link_a->link==NULL && link_b->link!=NULL) return  1;
  // Same position and both SNVs (non-overlapping first)
  bool is_link_a_snv_overlapping = edge_attributes_is_overlapping(link_a->attributes);
  bool is_link_b_snv_overlapping = edge_attributes_is_overlapping(link_b->attributes);
  if (!(is_link_a_snv_overlapping ^ is_link_b_snv_overlapping)) return 0; // Same link
  if (!is_link_a_snv_overlapping &&  is_link_b_snv_overlapping) {
    return -1;
  } else {
    return 1;
  }
}
int graph_text_builder_link_table_cmp_index_position(const graph_link_t** const ptr_link_a,const graph_link_t** const ptr_link_b) {
  const graph_link_t* const link_a = *ptr_link_a;
  const graph_link_t* const link_b = *ptr_link_b;
  const int64_t diff = ((int64_t)link_a->position_index - (int64_t)link_b->position_index);
  if (diff != 0) return diff;
  // Same position (SNVs first)
  if (link_a->link!=NULL && link_b->link==NULL) return -1;
  if (link_a->link==NULL && link_b->link!=NULL) return  1;
  // Same position and both SNVs (non-overlapping first)
  bool is_link_a_snv_overlapping = edge_attributes_is_overlapping(link_a->attributes);
  bool is_link_b_snv_overlapping = edge_attributes_is_overlapping(link_b->attributes);
  if (!(is_link_a_snv_overlapping ^ is_link_b_snv_overlapping)) return 0; // Same link
  if (!is_link_a_snv_overlapping &&  is_link_b_snv_overlapping) {
    return -1;
  } else {
    return 1;
  }
}
GEM_INLINE void graph_text_builder_link_table_sort(graph_text_builder_t* const graph_builder,const uint64_t maximum_tag_id) {
  GRAPH_TEXT_BUILDER_CHECK(graph_builder);
  // Allocate sorted graph-links tables
  const uint64_t link_table_length = svector_get_used(graph_builder->graph_links);
  graph_builder->link_table_sorted_mem =  mm_calloc(link_table_length+2,graph_link_t*,false);
  graph_builder->graph_links_idx_forward = mm_calloc(maximum_tag_id+1,graph_links_index_t,true);
  graph_builder->graph_links_idx_reverse = mm_calloc(maximum_tag_id+1,graph_links_index_t,true);
  graph_builder->graph_links_idx_length = maximum_tag_id;
  graph_builder->link_table_length = link_table_length;
  // Fill graph-links table
  graph_builder->link_table_sentinel.tag_id = UINT64_MAX; // Setup sentinel
  graph_builder->link_table_sorted_mem[0] = &graph_builder->link_table_sentinel; // Set first sentinel link (dummy)
  graph_builder->link_table_sorted = graph_builder->link_table_sorted_mem+1; // Skip first sentinel
  uint64_t i = 0;
  svector_iterator_new(&graph_builder->graph_links_iterator,graph_builder->graph_links,SVECTOR_READ_ITERATOR,0);
  while (!svector_read_iterator_eoi(&graph_builder->graph_links_iterator)) {
    graph_builder->link_table_sorted[i++] = svector_iterator_get_element(&graph_builder->graph_links_iterator,graph_link_t);
    svector_read_iterator_next(&graph_builder->graph_links_iterator);
  }
  graph_builder->link_table_sorted[i] = &graph_builder->link_table_sentinel; // Set last sentinel link (dummy)
  // Sort all links
  qsort(graph_builder->link_table_sorted,link_table_length,sizeof(graph_link_t*),
      (int (*)(const void *,const void *))graph_text_builder_link_table_cmp_text_position);
  // Build index table
  int64_t last_tag_id = INT64_MAX;
  for (i=0;i<link_table_length;++i) {
    const int64_t tag_id = graph_builder->link_table_sorted[i]->tag_id;
    if (tag_id != last_tag_id) {
      if (tag_id >= 0) { // Initialize
        graph_builder->graph_links_idx_forward[+tag_id].offset = i;
        graph_builder->graph_links_idx_forward[+tag_id].num_links = 1;
      } else {
        graph_builder->graph_links_idx_reverse[-tag_id].offset = i;
        graph_builder->graph_links_idx_reverse[-tag_id].num_links = 1;
      }
      last_tag_id = tag_id; // Set new last seen tag
    } else {
      if (tag_id >= 0) { // Increment count
        ++(graph_builder->graph_links_idx_forward[+tag_id].num_links);
      } else {
        ++(graph_builder->graph_links_idx_reverse[-tag_id].num_links);
      }
    }
  }
}
GEM_INLINE uint64_t graph_text_builder_link_table_get_length(graph_text_builder_t* const graph_builder,const int64_t tag_id) {
  return (tag_id > 0) ?
      graph_builder->graph_links_idx_forward[+tag_id].num_links :
      graph_builder->graph_links_idx_reverse[-tag_id].num_links;
}
/*
 * Step-2:: Query links (as to solve pending index-positions)
 */
GEM_INLINE void graph_text_builder_link_locator_iterate_forward(
    graph_text_builder_t* const graph_builder,
    graph_sorted_link_locator_t* const link_locator,const int64_t tag_id) {
  link_locator->strand = Forward;
  link_locator->tag_id = tag_id;
  // Check tagId in table
  if (ABS(tag_id) > graph_builder->graph_links_idx_length) {
    link_locator->next_link = NULL;
    link_locator->eoi = true;
  } else {
    // Check number of links
    graph_links_index_t* const link_table_index = (tag_id > 0) ?
        graph_builder->graph_links_idx_forward + ( tag_id) :
        graph_builder->graph_links_idx_reverse + (-tag_id) ;
    if (link_table_index->num_links > 0) {
      link_locator->next_link = graph_builder->link_table_sorted + link_table_index->offset;
      link_locator->eoi = false;
    } else {
      link_locator->next_link = NULL;
      link_locator->eoi = true;
    }
  }
}
GEM_INLINE void graph_text_builder_link_locator_iterate_backward(
    graph_text_builder_t* const graph_builder,
    graph_sorted_link_locator_t* const link_locator,const int64_t tag_id) {
  link_locator->strand = Reverse;
  link_locator->tag_id = tag_id;
  // Check tagId in table
  if (ABS(tag_id) > graph_builder->graph_links_idx_length) {
    link_locator->next_link = NULL;
    link_locator->eoi = true;
  } else {
    // Check number of links
    graph_links_index_t* const link_table_index = (tag_id > 0) ?
        graph_builder->graph_links_idx_forward + ( tag_id) :
        graph_builder->graph_links_idx_reverse + (-tag_id) ;
    if (link_table_index->num_links > 0) {
      link_locator->next_link = graph_builder->link_table_sorted + (link_table_index->offset + link_table_index->num_links - 1);
      link_locator->eoi = false;
    } else {
      link_locator->next_link = NULL;
      link_locator->eoi = true;
    }
  }
}
GEM_INLINE void graph_text_builder_link_locator_analyze_chunk_forward(graph_sorted_link_locator_t* const link_locator,const uint64_t text_position) {
  // Analyze all links related to text-position
  link_locator->link_table_chunk = link_locator->next_link;
  link_locator->num_links = 0;
  link_locator->has_overlapping_jump = false;
  link_locator->has_non_overlapping_jump = false;
  // Traverse all link related to text-position
  graph_link_t* iterate_link = *(link_locator->next_link);
  do {
    // Test link
    if (edge_attributes_is_overlapping(iterate_link->attributes)) {
      link_locator->has_overlapping_jump = true;
    } else {
      link_locator->has_non_overlapping_jump = true;
    }
    // Inc & load link
    ++(link_locator->num_links);
    ++(link_locator->next_link); // Forward
    // Check EOI
    iterate_link = *(link_locator->next_link); // Next
    if (iterate_link->tag_id!=link_locator->tag_id) {
      link_locator->eoi = true;
      break;
    }
  } while (iterate_link->position_text==text_position);
}
GEM_INLINE void graph_text_builder_link_locator_analyze_chunk_reverse(graph_sorted_link_locator_t* const link_locator,const uint64_t text_position) {
  // Analyze all links related to text-position
  link_locator->num_links = 0;
  link_locator->has_overlapping_jump = false;
  link_locator->has_non_overlapping_jump = false;
  // Traverse all link related to text-position
  graph_link_t* iterate_link = *(link_locator->next_link);
  do {
    // Test link
    if (edge_attributes_is_overlapping(iterate_link->attributes)) {
      link_locator->has_overlapping_jump = true;
    } else {
      link_locator->has_non_overlapping_jump = true;
    }
    // Inc & load link
    ++(link_locator->num_links);
    --(link_locator->next_link); // Backward
    // Check EOI
    iterate_link = *(link_locator->next_link); // Next
    if (iterate_link->tag_id!=link_locator->tag_id) {
      link_locator->eoi = true;
      break;
    }
  } while (iterate_link->position_text==text_position);
  link_locator->link_table_chunk = link_locator->next_link+1;
}
GEM_INLINE bool graph_text_builder_link_locator_find(
    graph_text_builder_t* const graph_builder,
    graph_sorted_link_locator_t* const link_locator,const int64_t tag_id,const uint64_t text_position) {
  link_locator->tag_id = tag_id;
  graph_links_index_t* const link_table_index = (tag_id > 0) ?
      graph_builder->graph_links_idx_forward + ( tag_id) :
      graph_builder->graph_links_idx_reverse + (-tag_id) ;
  if (link_table_index->offset != UINT64_MAX) {
    graph_link_t** const link_table_chunk = link_locator->link_table_chunk+link_table_index->offset;
    // Binary Search looking for text_position
    uint64_t inf=0, sup=link_table_index->num_links-1;
    while (sup > inf) {
      const uint64_t half=(sup+inf)/2;
      graph_link_t* const graph_link_half = *(link_table_chunk+half);
      if (graph_link_half->position_text < text_position) {
        inf = half+1;
      } else { // link_table_chunk_half->position_text >= text_position
        sup = half;
      }
    }
    // Check text_position
    graph_link_t* const graph_link = *(link_table_chunk+inf);
    if (graph_link->position_text != text_position) {
      // Found position
      link_locator->next_link = link_table_chunk+inf;
      link_locator->eoi = false;
      // Analyze all links related to text-position
      graph_text_builder_link_locator_analyze_chunk_forward(link_locator,text_position);
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}
GEM_INLINE bool graph_text_builder_link_locator_cmp(graph_sorted_link_locator_t* const link_locator,const uint64_t text_position) {
  // Check EOI & text_position
  if (gem_expect_true(link_locator->eoi || (*(link_locator->next_link))->position_text!=text_position)) return false;
  // Analyze all links related to text-position
  if (link_locator->strand == Forward) {
    graph_text_builder_link_locator_analyze_chunk_forward(link_locator,text_position);
  } else {
    graph_text_builder_link_locator_analyze_chunk_reverse(link_locator,text_position);
  }
  // Return cmp true
  return true;
}
GEM_INLINE bool graph_text_builder_link_locator_has_jump_non_overlapping(graph_sorted_link_locator_t* const link_locator) {
  return link_locator->has_non_overlapping_jump;
}
GEM_INLINE bool graph_text_builder_link_locator_has_jump_overlapping(graph_sorted_link_locator_t* const link_locator) {
  return link_locator->has_overlapping_jump;
}
GEM_INLINE char graph_text_builder_link_locator_jump_overlapping_get_reference_char(
    graph_sorted_link_locator_t* const link_locator) {
  uint64_t i;
  for (i=0;i<link_locator->num_links;++i) {
    graph_link_t* const iterate_link = *(link_locator->link_table_chunk+i);
    if (edge_attributes_is_overlapping(iterate_link->attributes)) {
      const char reference_char = edge_attributes_snv_get_reference(iterate_link->attributes);
      if (reference_char!=0) return reference_char;
    }
  }
  gem_fatal_error(GRAPH_TEXT_BUILDER_SNV_UNKNOWN_REFERENCE_CHAR);
  return 0;
}
GEM_INLINE void graph_text_builder_link_locator_solve_jump_non_overlapping(
    graph_sorted_link_locator_t* const link_locator,const uint64_t index_position) {
  uint64_t i;
  for (i=0;i<link_locator->num_links;++i) {
    graph_link_t* const iterate_link = *(link_locator->link_table_chunk+i);
    if (!edge_attributes_is_overlapping(iterate_link->attributes)) {
      iterate_link->position_index = index_position;
    }
  }
}
GEM_INLINE void graph_text_builder_link_locator_solve_jump_overlapping(
    graph_sorted_link_locator_t* const link_locator,const uint64_t index_position,const char reference_char) {
  uint64_t i;
  graph_link_t* iterate_link = NULL;
  for (i=0;i<link_locator->num_links;++i) {
    iterate_link = *(link_locator->link_table_chunk+i);
    if (edge_attributes_is_overlapping(iterate_link->attributes)) {
      iterate_link->position_index = index_position;
    }
  }
  // Overload last SNV-link with the reference character
  if (iterate_link) iterate_link->attributes = edge_attributes_snv_set_reference(iterate_link->attributes,reference_char);
}
/*
 * Step-3:: Write
 */
GEM_INLINE void graph_text_builder_write(
    fm_t* const file_manager,graph_text_builder_t* const graph_builder,locator_builder_t* const locator) {
  FM_CHECK(file_manager);
  GRAPH_TEXT_BUILDER_CHECK(graph_builder);
  /*
   * Preprocess
   */
  // Sort all links wrt index-position
  const uint64_t graph_links_table_size = svector_get_used(graph_builder->graph_links);
  qsort(graph_builder->link_table_sorted,graph_links_table_size,sizeof(graph_link_t*),
      (int (*)(const void *,const void *))graph_text_builder_link_table_cmp_index_position);
  /*
   * Process all graph-links
   *   - Check all links
   *   - Collapse all SNVs at the same position
   *   - Count link-groups
   */
  // Traverse all graph-links
  uint64_t num_edges=0, num_vertices=0;
  uint64_t i, last_position=UINT64_MAX;
  graph_link_t* last_graph_link = NULL;
  for (i=0;i<graph_links_table_size;++i) {
    // Check link (Force it to be solved (general graph)) // TODO
    graph_link_t* const graph_link = graph_builder->link_table_sorted[i];
    gem_cond_fatal_error(graph_link->position_index==GRAPH_TEXT_BUILDER_POSITION_UNKNOWN,
        GRAPH_TEXT_BUILDER_POSITION_UNSOLVED,
        locator_builder_get_tag(locator,graph_link->tag_id)->tag,
        (graph_link->tag_id > 0) ? '+' : '-',graph_link->position_text);
    /*
     * Collapse SNVs at the same position
     *   - Overlapping SNV
     *   - Non-Overlapping SNV
     */
    if (last_graph_link!=NULL &&
        graph_link->position_index==last_graph_link->position_index && /* Same Index Position */
        last_graph_link->link==NULL && graph_link->link==NULL && /* Both are SNV */
        !(edge_attributes_is_overlapping(last_graph_link->attributes) ^
          edge_attributes_is_overlapping(graph_link->attributes)) /* Both SNV type */ ) {
      // SNVs at the same position (Collapse)
      last_graph_link->attributes = last_graph_link->attributes | graph_link->attributes;
      // Nullify the current link
      graph_link->attributes = EDGE_ATTRIBUTES_NULL;
    } else {
      // Check vertex (different position => different vertex)
      if (graph_link->position_index!=last_position) {
        last_position = graph_link->position_index;
        ++num_vertices;
      }
      ++num_edges;
      last_graph_link = graph_link;
    }
  }
  graph_builder->link_table_length = num_edges;
  /*
   * Write Graph-Text
   */
  // Write Meta-Data
  fm_write_uint64(file_manager,num_vertices);
  fm_write_uint64(file_manager,num_edges);
  // Write Vertices table
  uint64_t non_overlapping_id=0;
  vertex_t vertex;
  num_edges = 0;
  last_position = UINT64_MAX;
  for (i=0;i<graph_links_table_size;++i) {
    graph_link_t* const graph_link = graph_builder->link_table_sorted[i];
    if (graph_link->attributes==EDGE_ATTRIBUTES_NULL) continue; // Skip null links
    // Check new vertex
    if (graph_link->position_index!=last_position) {
      // Dump vertex
      if (last_position!=UINT64_MAX) {
        vertex.non_overlapping_id = (edge_attributes_is_overlapping(vertex.edge_attributes)) ? UINT64_MAX : non_overlapping_id++;
        vertex.edge_end_idx = num_edges;
        fm_write_mem(file_manager,&vertex,sizeof(vertex_t)); // Write Vertex
      }
      // New vertex
      vertex.edge_begin_idx = num_edges;
      vertex.position = graph_link->position_index;
      vertex.edge_attributes = EDGE_ATTRIBUTES_NULL;
      last_position = graph_link->position_index;
    }
    // Compile group properties
    vertex.edge_attributes |= graph_link->attributes;
    ++num_edges;
  }
  if (last_position!=UINT64_MAX) { // Dump last vertex
    vertex.non_overlapping_id = (edge_attributes_is_overlapping(vertex.edge_attributes)) ? UINT64_MAX : non_overlapping_id++;
    vertex.edge_end_idx = num_edges;
    fm_write_mem(file_manager,&vertex,sizeof(vertex_t)); // Write Vertex
  }
  // Write Edges table
  for (i=0;i<graph_links_table_size;++i) {
    graph_link_t* const graph_link = graph_builder->link_table_sorted[i];
    if (graph_link->attributes==EDGE_ATTRIBUTES_NULL) continue; // Skip null links
    // Create edge
    edge_t edge;
    edge.position_src = graph_link->position_index;
    edge.position_dst = (graph_link->link == NULL) ? UINT64_MAX : graph_link->link->position_index;
    edge.attributes = graph_link->attributes;
    fm_write_mem(file_manager,&edge,sizeof(edge_t)); // Write Edge
  }
}
/*
 * Display
 */
GEM_INLINE void graph_text_builder_text_link_print(
    FILE* const stream,graph_link_t* const graph_link,locator_builder_t* const locator) {
  if (graph_link->link==NULL) { // SNVs
    fprintf(stream,"{"); edge_attributes_print(stream,graph_link->attributes); fprintf(stream,"}");
    fprintf(stream,"(%s:%c:%lu)",
        locator_builder_get_tag(locator,graph_link->tag_id)->tag,
        (graph_link->tag_id > 0) ? '+' : '-',graph_link->position_text);
  } else { // General Link
    graph_link_t* link = graph_link;
    fprintf(stream,"{"); edge_attributes_print(stream,link->attributes); fprintf(stream,"}");
    fprintf(stream,"(%s:%c:%lu)",
        locator_builder_get_tag(locator,link->tag_id)->tag,
        (link->tag_id > 0) ? '+' : '-',link->position_text);
    link = graph_link->link;
    fprintf(stream,"\t{"); edge_attributes_print(stream,link->attributes); fprintf(stream,"}");
    fprintf(stream,"(%s:%c:%lu)",
        locator_builder_get_tag(locator,link->tag_id)->tag,
        (link->tag_id > 0) ? '+' : '-',link->position_text);
  }
}
GEM_INLINE void graph_text_builder_index_link_print(FILE* const stream,graph_link_t* const graph_link) {
  if (graph_link->link==NULL) { // SNVs
    fprintf(stream,"{"); edge_attributes_print(stream,graph_link->attributes); fprintf(stream,"}");
    fprintf(stream,"(idx:%lu)",graph_link->position_index);
  } else { // General Link
    graph_link_t* link = graph_link;
    fprintf(stream,"{"); edge_attributes_print(stream,link->attributes); fprintf(stream,"}");
    fprintf(stream,"(idx:%lu)",link->position_index);
    link = graph_link->link;
    fprintf(stream,"\t{"); edge_attributes_print(stream,link->attributes); fprintf(stream,"}");
    fprintf(stream,"(idx:%lu)",link->position_index);
  }
}
GEM_INLINE void graph_text_builder_print(
    FILE* const stream,graph_text_builder_t* const graph_builder,
    locator_builder_t* const locator,const bool display_links) {
  GEM_CHECK_NULL(stream);
  GRAPH_TEXT_BUILDER_CHECK(graph_builder);
  // Print graph info
  tab_fprintf(stream,"[GEM]>Graph.Builder\n");
  tab_fprintf(stream,"  => Num.Links  %lu\n",svector_get_used(graph_builder->graph_links));
  tab_fprintf(stream,"  => Size.Links %lu MB\n",CONVERT_B_TO_MB(graph_builder->link_table_length*sizeof(graph_link_t)));
  tab_fprintf(stream,"  => Size.Edges %lu MB\n",CONVERT_B_TO_MB(graph_builder->link_table_length*sizeof(edge_t)));
  // Display all links
  if (display_links) {
    tab_fprintf(stream,"  => Graph.Builder.links\n");
    uint64_t i = 0;
    svector_iterator_t graph_links_iterator;
    svector_iterator_new(&graph_links_iterator,graph_builder->graph_links,SVECTOR_READ_ITERATOR,0);
    while (!svector_read_iterator_eoi(&graph_links_iterator)) {
      graph_link_t* const graph_link = svector_iterator_get_element(&graph_links_iterator,graph_link_t);
      if (graph_link->attributes==EDGE_ATTRIBUTES_NULL) continue;
      tab_fprintf(stream,"       [%03lu]\t",i);
      graph_text_builder_text_link_print(stream,graph_link,locator);
      fprintf(stream,"\n");
      // Next
      svector_read_iterator_next(&graph_links_iterator);
      i++;
    }
  }
  // Flush
  fflush(stream);
}
GEM_INLINE void graph_text_builder_link_table_print(
    FILE* const stream,graph_text_builder_t* const graph_builder,
    locator_builder_t* const locator,const bool display_links) {
  GEM_CHECK_NULL(stream);
  GRAPH_TEXT_BUILDER_CHECK(graph_builder);
  // Print graph info
  const uint64_t graph_links_table_size = svector_get_used(graph_builder->graph_links);
  tab_fprintf(stream,"[GEM]>Graph.Builder {Sorted}\n");
  tab_fprintf(stream,"  => Num.Links  %lu\n",graph_links_table_size);
  tab_fprintf(stream,"  => Size.Links %lu MB\n",CONVERT_B_TO_MB(graph_builder->link_table_length*sizeof(graph_link_t)));
  tab_fprintf(stream,"  => Size.Edges %lu MB\n",CONVERT_B_TO_MB(graph_builder->link_table_length*sizeof(edge_t)));
  // Display all links
  if (display_links) {
    tab_fprintf(stream,"  => Graph.Builder.links\n");
    uint64_t i;
    for (i=0;i<graph_links_table_size;++i) {
      graph_link_t* const graph_link = graph_builder->link_table_sorted[i];
      if (graph_link->attributes==EDGE_ATTRIBUTES_NULL) continue;
      tab_fprintf(stream,"       [%03lu]\t",i);
      graph_text_builder_index_link_print(stream,graph_link);
      fprintf(stream,"\n");
    }
  }
  // Flush
  fflush(stream);
}
