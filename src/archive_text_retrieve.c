/*
 * PROJECT: GEMMapper
 * FILE: archive_text_retrieve.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "archive_text_retrieve.h"

///*
// * Position Locator
// */
//GEM_INLINE void archive_text_append_starting_position_at_interval(
//    const locator_interval_t* const loc_interval,const uint64_t position_no_offset,
//    svector_iterator_t* const positions_iterator) {
//  const int64_t starting_position = MAX(loc_interval->begin_position,position_no_offset);
//  *svector_iterator_get_element(positions_iterator,uint64_t) = starting_position;
//  svector_write_iterator_next(positions_iterator); // Next
//}
//GEM_INLINE void archive_text_append_starting_position_rec(
//    archive_t* const archive,
//    const locator_interval_t* const loc_interval,const int64_t current_vertex_index,
//    const int64_t current_index_position,const int64_t current_right_offset,
//    svector_iterator_t* const positions_iterator) {
//  GEM_CHECK_NOT_NEGATIVE(current_index_position);
//  GEM_CHECK_NOT_NEGATIVE(current_right_offset);
//  // (1) No edge-group
//  const int64_t position_no_offset =
//      (current_index_position > current_right_offset) ? current_index_position-current_right_offset : 0;
//  if (gem_expect_false(current_vertex_index < 0)) {
//    archive_text_append_starting_position_at_interval(loc_interval,position_no_offset,positions_iterator);
//    return;
//  }
//  // (2) Vertex not contained in the interval
//  const vertex_t* vertex = graph_text_get_vertex(archive->graph,current_vertex_index);
//  if (vertex->position < loc_interval->begin_position) {
//    // Vertex not contained in the interval
//    archive_text_append_starting_position_at_interval(loc_interval,position_no_offset,positions_iterator);
//    return;
//  }
//  // (3) Vertex not contained in the text interval
//  if (vertex->position < position_no_offset) {
//    // Vertex not contained in the text interval
//    *svector_iterator_get_element(positions_iterator,uint64_t) = position_no_offset;
//    svector_write_iterator_next(positions_iterator); // Next
//    return;
//  }
//  // (4) (Check the properties of edge-group in the text interval)
//  const int64_t distance_to_edge = current_index_position - vertex->position;
//  const int64_t new_right_offset = current_right_offset - distance_to_edge;
//  GEM_CHECK_NOT_NEGATIVE(new_right_offset);
//  // (4.1) Destiny-jump (follow jump connections backwards)
//  if (edge_attributes_is_jump_destiny(vertex->edge_attributes)) {
//    uint64_t edge_idx;
//    for (edge_idx=vertex->edge_begin_idx;edge_idx<vertex->edge_end_idx;++edge_idx) {
//      const edge_t* const edge = graph_text_get_edge(archive->graph,edge_idx);
//      const uint64_t vertex_dst = graph_text_lookup_vertex_index(archive->graph,edge->position_dst);
//      const uint64_t locator_interval_dst_index = locator_lookup_interval_index(archive->locator,edge->position_dst);
//      const locator_interval_t* const loc_interval_dst = locator_get_interval(archive->locator,locator_interval_dst_index);
//      archive_text_append_starting_position_rec(
//          archive,loc_interval_dst,vertex_dst,
//          edge->position_dst,new_right_offset,positions_iterator);
//    }
//  }
//  // (4.2) Overlapping SNV (Follow the reference)
//  if (edge_attributes_is_overlapping(vertex->edge_attributes)) {
//    archive_text_append_starting_position_rec(
//        archive,loc_interval,current_vertex_index-1,
//        vertex->position,new_right_offset,positions_iterator);
//    return;
//  }
//  // (4.3) Source-Jump or non-overlapping SNV (Ghost character; Follow the reference)
//  archive_text_append_starting_position_rec(
//      archive,loc_interval,current_vertex_index-1,
//      vertex->position,new_right_offset+1,positions_iterator);
//}
//GEM_INLINE void archive_text_append_starting_position(
//    archive_t* const archive,const uint64_t index_position,const uint64_t right_offset,
//    svector_iterator_t* const positions_iterator) {
//  ARCHIVE_CHECK(archive);
//  SEGMENTED_VECTOR_ITERATOR_CHECK(positions_iterator);
//  // Get locator interval boundaries
//  const uint64_t locator_interval_index = locator_lookup_interval_index(archive->locator,index_position);
//  const locator_interval_t* const loc_interval = locator_get_interval(archive->locator,locator_interval_index);
//  // Check graph
//  if (archive->graph == NULL) {
//    archive_text_append_starting_position_rec(archive,loc_interval,-1,index_position,right_offset,positions_iterator);
//  } else {
//    // Get closest previous jump
//    uint64_t vertex_index = graph_text_lookup_previous_vertex_index(archive->graph,index_position,right_offset);
//    archive_text_append_starting_position_rec(archive,loc_interval,vertex_index,index_position,right_offset,positions_iterator);
//  }
//}
///*
// * Text Retriever
// */
//GEM_INLINE void archive_text_retrieve_continuous_text_within_interval(
//    const locator_interval_t* const loc_interval,const uint64_t current_index_position,
//    const uint64_t remaining_length,text_block_t* const current_text_block) {
//  // TODO
//}
//GEM_INLINE void archive_text_retrieve_rec(
//    archive_t* const archive,
//    const locator_interval_t* const loc_interval,const int64_t current_vertex_index,
//    const uint64_t current_index_position,const int64_t remaining_length,
//    text_collection_t* const text_collection,text_block_t* current_text_block,const uint64_t num_trace_blocks) {
//  ARCHIVE_CHECK_INDEX_POSITION(archive,current_index_position);
//  GEM_CHECK_NOT_NEGATIVE(remaining_length);
//  // Request new text-block if we come from a jump
//  if (current_text_block==NULL) current_text_block = text_collection_new_block(text_collection);
//  // (1) No Vertex not contained in the interval
//  if (gem_expect_false(current_vertex_index < 0)) {
//    // Retrieve continuous text-block
//    archive_text_retrieve_continuous_text_within_interval(loc_interval,current_index_position,remaining_length,current_text_block);
//    // Request new trace, set the last trace-block & return
//    text_trace_t* const text_trace = text_collection_new_trace(text_collection,num_trace_blocks+1);
//    *text_collection_get_trace_block(text_collection,text_trace,num_trace_blocks) = current_text_block;
//    return;
//  }
//  const int64_t end_position = current_index_position+remaining_length-1;
//  const vertex_t* vertex = graph_text_get_vertex(archive->graph,current_vertex_index);
//  // (2) Vertex not contained in the locator-interval
//  // (3) Vertex beyond the text scope
//  if (loc_interval->end_position <= vertex->position || end_position <= vertex->position) {
//    // Retrieve continuous text-block
//    archive_text_retrieve_continuous_text_within_interval(loc_interval,current_index_position,remaining_length,current_text_block);
//    // Request new trace, set the last trace-block & return
//    text_trace_t* const text_trace = text_collection_new_trace(text_collection,num_trace_blocks+1);
//    *text_collection_get_trace_block(text_collection,text_trace,num_trace_blocks) = current_text_block;
//    return;
//  }
//  // (4) Vertex in between (Check the properties of edge-group in the text interval)
//  const int64_t distance_to_edge = current_index_position-vertex->position;
//  const int64_t new_remaining_length = remaining_length - distance_to_edge;
//  GEM_CHECK_NOT_NEGATIVE(new_remaining_length);
//  // Retrieve continuous text-block
//  const uint64_t initial_num_traces = text_collection_get_num_traces(text_collection);
//
//  // TODO [current_index_position,current_index_position+distance_to_edge)
//
//  // (4.1) Destiny-jump (follow jump connections backwards)
//  if (edge_attributes_is_jump_source(vertex->edge_attributes)) {
//    uint64_t edge_idx;
//    for (edge_idx=vertex->edge_begin_idx;edge_idx<vertex->edge_end_idx;++edge_idx) {
//      const edge_t* const edge = graph_text_get_edge(archive->graph,edge_idx);
//      const uint64_t vertex_dst = graph_text_lookup_vertex_index(archive->graph,edge->position_dst);
//      const uint64_t locator_interval_dst_index = locator_lookup_interval_index(archive->locator,edge->position_dst);
//      const locator_interval_t* const loc_interval_dst = locator_get_interval(archive->locator,locator_interval_dst_index);
//      archive_text_retrieve_rec(archive,loc_interval_dst,vertex_dst,
//          edge->position_dst,new_remaining_length,text_collection,current_text_block,num_trace_blocks+1);
//    }
//  }
//  if (edge_attributes_is_overlapping(vertex->edge_attributes)) {
//    // (4.2) Overlapping SNV (Follow the reference)
//    archive_text_retrieve_rec(archive,loc_interval,current_vertex_index+1,
//        vertex->position,new_remaining_length,text_collection,current_text_block,num_trace_blocks);
//  } else {
//    // (4.3) Source-Jump or non-overlapping SNV (Ghost character; Follow the reference)
//    archive_text_retrieve_rec(archive,loc_interval,current_vertex_index+1,
//        vertex->position,new_remaining_length+1,text_collection,current_text_block,num_trace_blocks);
//  }
//  // Add text-block to all text-traces generated
//
//  // FIXME Not all instances add text-block
//
//  const uint64_t total_num_traces = text_collection_get_num_traces(text_collection);
//  uint64_t i;
//  for (i=initial_num_traces;i<total_num_traces;++i) {
//    text_trace_t* const text_trace = text_collection_get_trace(text_collection,i);
//    *text_collection_get_trace_block(text_collection,text_trace,num_trace_blocks) = current_text_block;
//  }
//
//
//
//
//}
GEM_INLINE uint64_t archive_text_retrieve(
    const locator_t* const locator,graph_text_t* const graph,const dna_text_t* const enc_text,
    text_collection_t* const text_collection,const uint64_t text_position,const uint64_t length,
    uint64_t* const text_trace_first_offset,mm_stack_t* const mm_stack) {
  // Allocate text-trace
  const uint64_t text_trace_offset = text_collection_new_trace(text_collection);
  *text_trace_first_offset = text_trace_offset;
  // Retrieve sequence
  text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
  text_trace->length = length;
  text_trace->text = dna_text_retrieve_sequence(enc_text,text_position,text_trace->length,mm_stack);
  // [[ TODO GRAPH COMPILANT]]
  // Return
  return 1; // Number of traces retrieved
}






