/*
 * PROJECT: GEMMapper
 * FILE: archive_text.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "archive/archive_text.h"
#include "archive/archive_text_rl.h"

/*
 * Archive-Text Model & Version
 */
#define ARCHIVE_TEXT_MODEL_NO  1006ull

/*
 * Builder
 */
void archive_text_write(
    fm_t* const file_manager,
    dna_text_t* const enc_text,
    const bool explicit_complement,
    const uint64_t forward_text_length,
    sampled_rl_t* const sampled_rl,
    const bool verbose) {
  // Write Header (Meta-data)
  fm_write_uint64(file_manager,ARCHIVE_TEXT_MODEL_NO);
  fm_write_uint64(file_manager,(sampled_rl!=NULL) ? 1ul : 0ul); // RL-Text
  fm_write_uint64(file_manager,(explicit_complement) ? 1ul : 0ul); // Explicit RC-text
  fm_write_uint64(file_manager,forward_text_length); // Total length of the forward text
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
archive_text_t* archive_text_read_mem(mm_t* const memory_manager) {
  // Allocate
  archive_text_t* const archive_text = mm_alloc(archive_text_t);
  // Read Header
  const uint64_t archive_text_model_no = mm_read_uint64(memory_manager);
  gem_cond_fatal_error(archive_text_model_no!=ARCHIVE_TEXT_MODEL_NO,
      ARCHIVE_TEXT_WRONG_MODEL_NO,archive_text_model_no,(uint64_t)ARCHIVE_TEXT_MODEL_NO);
  archive_text->run_length = (mm_read_uint64(memory_manager)==1); // RL-Text
  archive_text->explicit_complement = (mm_read_uint64(memory_manager)==1); // Explicit RC-text
  archive_text->forward_text_length = mm_read_uint64(memory_manager); // Total length of the forward text
  // Text
  archive_text->enc_text = dna_text_read_mem(memory_manager);
  // Load Sampled RL-Index Positions
  if (archive_text->run_length) archive_text->sampled_rl = sampled_rl_read_mem(memory_manager);
  // Return
  return archive_text;
}
void archive_text_delete(archive_text_t* const archive_text) {
  // Delete Text
  dna_text_delete(archive_text->enc_text);
  // Delete Sampled RL-Index Positions
  if (archive_text->run_length) sampled_rl_delete(archive_text->sampled_rl);
  mm_free(archive_text);
}
/*
 * Accessors
 */
uint64_t archive_text_get_size(archive_text_t* const archive_text) {
  const uint64_t graph_size = 0;
  const uint64_t text_size = dna_text_get_size(archive_text->enc_text);
  const uint64_t sampled_rl_size = 0; // TODO
  return graph_size+text_size+sampled_rl_size;
}
strand_t archive_text_get_position_strand(
    archive_text_t* const archive_text,
    const uint64_t index_position) {
  return index_position < archive_text->forward_text_length ? Forward : Reverse;
}
uint64_t archive_text_get_unitary_projection(
    archive_text_t* const archive_text,
    const uint64_t index_position) {
  return 2*archive_text->forward_text_length - index_position - 2;
}
uint64_t archive_text_get_projection(
    archive_text_t* const archive_text,
    const uint64_t index_position,
    const uint64_t length) {
  return 2*archive_text->forward_text_length - index_position - length - 1;
}
/*
 * Archive Text Retriever
 */
void archive_text_retrieve(
    archive_text_t* const archive_text,
    const uint64_t text_position,
    const uint64_t text_length,
    const bool reverse_complement_text,
    const bool run_length_text,
    text_trace_t* const text_trace,
    mm_stack_t* const mm_stack) {
  // Retrieve text
  text_trace->text_length = text_length;
  if (text_position < archive_text->forward_text_length || archive_text->explicit_complement) {
    if (reverse_complement_text) {
      if (archive_text->explicit_complement) {
        const uint64_t position_fprojection = archive_text_get_projection(archive_text,text_position,text_length);
        text_trace->text = dna_text_retrieve_sequence(
            archive_text->enc_text,position_fprojection,text_length,mm_stack);
      } else {
        // Reverse-Complement the text
        uint8_t* const text = dna_text_retrieve_sequence(archive_text->enc_text,text_position,text_length,mm_stack);
        text_trace->text = mm_stack_calloc(mm_stack,text_length,uint8_t,false);
        uint64_t i_forward, i_backward;
        for (i_forward=0,i_backward=text_length-1;i_forward<text_length;++i_forward,--i_backward) {
          text_trace->text[i_forward] = dna_encoded_complement(text[i_backward]);
        }
      }
    } else {
      text_trace->text = dna_text_retrieve_sequence(archive_text->enc_text,text_position,text_length,mm_stack);
    }
  } else {
    // Forward projection
    const uint64_t position_fprojection = archive_text_get_projection(archive_text,text_position,text_length);
    uint8_t* const text = dna_text_retrieve_sequence(archive_text->enc_text,position_fprojection,text_length,mm_stack);
    if (reverse_complement_text) {
      text_trace->text = text;
    } else {
      // Reverse-Complement the text
      text_trace->text = mm_stack_calloc(mm_stack,text_length,uint8_t,false);
      uint64_t i_forward, i_backward;
      for (i_forward=0,i_backward=text_length-1;i_forward<text_length;++i_forward,--i_backward) {
        text_trace->text[i_forward] = dna_encoded_complement(text[i_backward]);
      }
    }
  }
  // Compute RL-text
  if (run_length_text) {
    // Allocate RL-Encoded Text
    text_trace->rl_text = mm_stack_calloc(mm_stack,text_length,uint8_t,false);
    text_trace->rl_runs_acc = mm_stack_calloc(mm_stack,text_length,uint32_t,false);
    // RL encode
    archive_text_rl_encode(
        text_trace->text,text_length,text_trace->rl_text,
        &text_trace->rl_text_length,text_trace->rl_runs_acc);
  }
}
uint64_t archive_text_retrieve_collection(
    archive_text_t* const archive_text,
    const text_collection_t* const text_collection,
    const uint64_t text_position,
    const uint64_t text_length,
    const bool reverse_complement_text,
    const bool run_length_text,
    mm_stack_t* const mm_stack) {
  // Allocate text-trace
  const uint64_t text_trace_offset = text_collection_new_trace(text_collection);
  text_trace_t* const text_trace = text_collection_get_trace(text_collection,text_trace_offset);
  // Retrieve sequence
  archive_text_retrieve(archive_text,text_position,text_length,
      reverse_complement_text,run_length_text,text_trace,mm_stack);
  // Return
  return text_trace_offset;
}
/*
 * Display
 */
void archive_text_print(
    FILE* const stream,
    const archive_text_t* const archive_text) {
  const uint64_t text_size = dna_text_get_size(archive_text->enc_text);
  const uint64_t sampled_rl_size = 0; // TODO
  const uint64_t archive_text_size = text_size + sampled_rl_size;
  // Display
  tab_fprintf(stream,"[GEM]>Archive.Text\n");
  tab_fprintf(stream,"  => Archive.RL.Text             %s\n",(archive_text->run_length)?"Yes":"No");
  tab_fprintf(stream,"  => Archive.ExplicitComplement  %s\n",(archive_text->explicit_complement)?"Yes":"No");
  tab_fprintf(stream,"  => Archive.Length              %"PRIu64"\n",dna_text_get_length(archive_text->enc_text));
  tab_fprintf(stream,"    => Archive.Forward.Length    %"PRIu64"\n",archive_text->forward_text_length);
  tab_fprintf(stream,"  => Archive.Text.Size %"PRIu64" MB (100%%)\n",CONVERT_B_TO_MB(archive_text_size));
  tab_fprintf(stream,"    => Text.Size       %"PRIu64" MB (%2.3f%%)\n",CONVERT_B_TO_MB(text_size),PERCENTAGE(text_size,archive_text_size));
  tab_fprintf(stream,"    => SampledRL.Size  %"PRIu64" MB (%2.3f%%)\n",CONVERT_B_TO_MB(sampled_rl_size),PERCENTAGE(sampled_rl_size,archive_text_size));
  /*
   * Components Display
   */
  // Archive Text
  tab_global_inc();
  dna_text_print(stream,archive_text->enc_text,dna_text_get_length(archive_text->enc_text));
  tab_global_dec();
  // Sampled RL-Text
  tab_global_inc();
  sampled_rl_print(stream,archive_text->sampled_rl);
  tab_global_dec();
}
