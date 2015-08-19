/*
 * PROJECT: GEMMapper
 * FILE: archive_text.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "archive_text.h"

/*
 * Archive-Text Model & Version
 */
#define ARCHIVE_TEXT_MODEL_NO  1002ull

/*
 * Builder
 */
GEM_INLINE void archive_text_write(
    fm_t* const file_manager,dna_text_t* const enc_text,
    const bool explicit_complement,const uint64_t forward_text_length,
    sampled_rl_t* const sampled_rl,graph_text_builder_t* const graph,const bool verbose) {
  // Write Header (Meta-data)
  fm_write_uint64(file_manager,ARCHIVE_TEXT_MODEL_NO);
  fm_write_uint64(file_manager,(graph!=NULL) ? 1ul : 0ul); // Hypertext
  fm_write_uint64(file_manager,(sampled_rl!=NULL) ? 1ul : 0ul); // RL-Text
  fm_write_uint64(file_manager,(explicit_complement) ? 1ul : 0ul); // Explicit RC-text
  fm_write_uint64(file_manager,forward_text_length); // Total length of the forward text
  // Graph
  // TODO
  // Text
  if (explicit_complement) {
    // Save all (except extra separator)
    dna_text_write_chunk(file_manager,enc_text,dna_text_get_length(enc_text)-1);
    if (verbose) dna_text_print(gem_info_get_stream(),enc_text,dna_text_get_length(enc_text)-1);
  } else {
    // Save just forward text
    dna_text_write_chunk(file_manager,enc_text,forward_text_length);
    if (verbose) dna_text_print(gem_info_get_stream(),enc_text,forward_text_length);
  }
  // Sampled RL-Index Positions
  if (sampled_rl!=NULL) sampled_rl_write(file_manager,sampled_rl);
}
/*
 * Loader
 */
GEM_INLINE archive_text_t* archive_text_read_mem(mm_t* const memory_manager) {
  // Allocate
  archive_text_t* const archive_text = mm_alloc(archive_text_t);
  // Read Header
  const uint64_t archive_text_model_no = mm_read_uint64(memory_manager);
  gem_cond_fatal_error(archive_text_model_no!=ARCHIVE_TEXT_MODEL_NO,
      ARCHIVE_TEXT_WRONG_MODEL_NO,archive_text_model_no,ARCHIVE_TEXT_MODEL_NO);
  archive_text->hypertext = (mm_read_uint64(memory_manager)==1); // Hypertext
  archive_text->run_length = (mm_read_uint64(memory_manager)==1); // RL-Text
  archive_text->explicit_complement = (mm_read_uint64(memory_manager)==1); // Explicit RC-text
  archive_text->forward_text_length = mm_read_uint64(memory_manager); // Total length of the forward text
  // Graph
  // TODO archive->graph = (archive->index_type == fm_dna_graph) ? graph_text_read_mem(archive->mm) : NULL;
  // Text
  archive_text->enc_text = dna_text_read_mem(memory_manager);
  // Load Sampled RL-Index Positions
  if (archive_text->run_length) archive_text->sampled_rl = sampled_rl_read_mem(memory_manager);
  // Return
  return archive_text;
}
GEM_INLINE void archive_text_delete(archive_text_t* const archive_text) {
  // Delete Graph
  // TODO if (archive->index_type == fm_dna_graph) graph_text_delete(archive->graph);
  // Delete Text
  dna_text_delete(archive_text->enc_text);
  // Delete Sampled RL-Index Positions
  if (archive_text->run_length) sampled_rl_delete(archive_text->sampled_rl);
  mm_free(archive_text);
}
/*
 * Accessors
 */
GEM_INLINE uint64_t archive_text_get_size(archive_text_t* const archive_text) {
  const uint64_t graph_size = archive_text->hypertext ? graph_text_get_size(archive_text->graph) : 0;
  const uint64_t text_size = dna_text_get_size(archive_text->enc_text);
  const uint64_t sampled_rl_size = 0; // TODO
  return graph_size+text_size+sampled_rl_size;
}
GEM_INLINE strand_t archive_text_get_position_strand(
    archive_text_t* const archive_text,const uint64_t index_position) {
  return index_position < archive_text->forward_text_length ? Forward : Reverse;
}
GEM_INLINE uint64_t archive_text_get_unitary_projection(
    archive_text_t* const archive_text,const uint64_t index_position) {
  return 2*archive_text->forward_text_length - index_position - 2;
}
GEM_INLINE uint64_t archive_text_get_projection(
    archive_text_t* const archive_text,const uint64_t index_position,const uint64_t length) {
  return 2*archive_text->forward_text_length - index_position - length - 1;
}
/*
 * Archive Text Retriever
 */
GEM_INLINE uint64_t archive_text_retrieve(
    archive_text_t* const archive_text,const text_collection_t* const text_collection,
    const uint64_t index_position,const uint64_t length,
    const bool reverse_complement_text,mm_stack_t* const mm_stack) {
  // Allocate text-trace
  const uint64_t text_trace_offset = text_collection_new_trace(text_collection);
  // Retrieve sequence
  text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
  text_trace->length = length;
  if (index_position < archive_text->forward_text_length || archive_text->explicit_complement) {
    if (reverse_complement_text) {
      if (archive_text->explicit_complement) {
        const uint64_t position_fprojection = archive_text_get_projection(archive_text,index_position,length);
        text_trace->text = dna_text_retrieve_sequence(archive_text->enc_text,position_fprojection,length,mm_stack);
      } else {
        // Reverse-Complement the text
        uint8_t* const text = dna_text_retrieve_sequence(archive_text->enc_text,index_position,length,mm_stack);
        text_trace->text = mm_stack_calloc(mm_stack,length,uint8_t,false);
        uint64_t i_forward, i_backward;
        for (i_forward=0,i_backward=length-1;i_forward<length;++i_forward,--i_backward) {
          text_trace->text[i_forward] = dna_encoded_complement(text[i_backward]);
        }
      }
    } else {
      text_trace->text = dna_text_retrieve_sequence(archive_text->enc_text,index_position,length,mm_stack);
    }
  } else {
    // Forward projection
    const uint64_t position_fprojection = archive_text_get_projection(archive_text,index_position,length);
    uint8_t* const text = dna_text_retrieve_sequence(archive_text->enc_text,position_fprojection,length,mm_stack);
    if (reverse_complement_text) {
      text_trace->text = text;
    } else {
      // Reverse-Complement the text
      text_trace->text = mm_stack_calloc(mm_stack,length,uint8_t,false);
      uint64_t i_forward, i_backward;
      for (i_forward=0,i_backward=length-1;i_forward<length;++i_forward,--i_backward) {
        text_trace->text[i_forward] = dna_encoded_complement(text[i_backward]);
      }
    }
  }
  // Return
  return text_trace_offset;
}
/*
 * Display
 */
GEM_INLINE void archive_text_print(FILE* const stream,const archive_text_t* const archive_text) {
  const uint64_t graph_size = archive_text->hypertext ? graph_text_get_size(archive_text->graph) : 0;
  const uint64_t text_size = dna_text_get_size(archive_text->enc_text);
  const uint64_t sampled_rl_size = 0; // TODO
  const uint64_t archive_text_size = graph_size + text_size + sampled_rl_size;
  // Display
  tab_fprintf(stream,"[GEM]>Archive.Text\n");
  tab_fprintf(stream,"  => Archive.HyperText           %s\n",(archive_text->hypertext)?"Yes":"No");
  tab_fprintf(stream,"  => Archive.RL.Text             %s\n",(archive_text->run_length)?"Yes":"No");
  tab_fprintf(stream,"  => Archive.ExplicitComplement  %s\n",(archive_text->explicit_complement)?"Yes":"No");
  tab_fprintf(stream,"  => Archive.Length              %"PRIu64"\n",dna_text_get_length(archive_text->enc_text));
  tab_fprintf(stream,"    => Archive.Forward.Length    %"PRIu64"\n",archive_text->forward_text_length);
  tab_fprintf(stream,"  => Archive.Text.Size %"PRIu64" MB (100%%)\n",CONVERT_B_TO_MB(archive_text_size));
  tab_fprintf(stream,"    => Graph.Size      %"PRIu64" MB (%2.3f%%)\n",CONVERT_B_TO_MB(graph_size),PERCENTAGE(graph_size,archive_text_size));
  tab_fprintf(stream,"    => Text.Size       %"PRIu64" MB (%2.3f%%)\n",CONVERT_B_TO_MB(text_size),PERCENTAGE(text_size,archive_text_size));
  tab_fprintf(stream,"    => SampledRL.Size  %"PRIu64" MB (%2.3f%%)\n",CONVERT_B_TO_MB(sampled_rl_size),PERCENTAGE(sampled_rl_size,archive_text_size));
  /*
   * Components Display
   */
  // Graph
  // TODO
  // Archive Text
  tab_global_inc();
  dna_text_print(stream,archive_text->enc_text,dna_text_get_length(archive_text->enc_text));
  tab_global_dec();
  // Sampled RL-Text
  tab_global_inc();
  sampled_rl_print(stream,archive_text->sampled_rl);
  tab_global_dec();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

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


