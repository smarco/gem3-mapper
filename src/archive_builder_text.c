/*
 * PROJECT: GEMMapper
 * FILE: archive_builder_text.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_builder_text.h"
#include "archive_builder_text_parser.h"

/*
 * Generate DNA-Text
 */
GEM_INLINE void archive_builder_text_process(
    archive_builder_t* const archive_builder,input_file_t* const input_multifasta,
    const bool dump_locator_intervals,const bool dump_indexed_text,const bool verbose) {
  /*
   * Text
   */
  // Inspect Text
  archive_builder_inspect_text(archive_builder,input_multifasta,verbose);
  // Generate Text (Forward)
  archive_builder_generate_forward_text(archive_builder,input_multifasta,verbose);
  // Generate C2T & G2A Texts
  if (archive_builder->type==archive_dna_bisulfite) {
    archive_builder_generate_bisulfite_text(archive_builder,verbose);
  }
  // Set forward-text length
  archive_builder->forward_text_length = archive_builder->parsing_state.index_position;
  // Generate RC-Text
  if (archive_builder->indexed_complement == index_complement_yes) {
    archive_builder_generate_rc_text(archive_builder,verbose);
  }
  archive_builder_generate_text_add_separator(archive_builder); // Add extra separator (Close full-text)
  // Set full-text length
  dna_text_set_length(archive_builder->enc_text,archive_builder->parsing_state.index_position);
  /*
   * DEBUG
   */
  // DEBUG locator
  locator_builder_print(gem_info_get_stream(),archive_builder->locator,dump_locator_intervals); // Locator
  // DEBUG index_text
  if (dump_indexed_text) archive_builder_text_dump(archive_builder,".text");
  /*
   * Write Header & Locator
   */
  archive_builder_write_header(archive_builder);
  archive_builder_write_locator(archive_builder);
  locator_builder_delete(archive_builder->locator); // Free Locator
}
/*
 * Run-length Text (Apply RL to the text)
 */
GEM_INLINE void archive_builder_text_apply_run_length(
    archive_builder_t* const archive_builder,const bool dump_run_length_text,const bool verbose) {
  // Allocate RL-text (Circular BWT extra)
  const uint64_t enc_text_length = dna_text_get_length(archive_builder->enc_text);
  archive_builder->enc_rl_text = dna_text_padded_new(enc_text_length,2,SA_BWT_PADDED_LENGTH);
  const uint64_t max_num_samples = DIV_CEIL(enc_text_length,SAMPLED_RL_SAMPLING_RATE);
  archive_builder->sampled_rl = sampled_rl_new(SAMPLED_RL_SAMPLING_RATE,max_num_samples,enc_text_length);
  // Compact the text into the RL-text
  const uint8_t* const enc_text = dna_text_get_text(archive_builder->enc_text);
  uint8_t* const enc_rl_text = dna_text_get_text(archive_builder->enc_rl_text);
  uint64_t text_position, rl_text_position=1;
  uint64_t run_length=1, num_rl_samples=1;
  enc_rl_text[0] = enc_text[0]; // Add character
  sampled_rl_sample(archive_builder->sampled_rl,0,0);
  for (text_position=1;text_position<enc_text_length;++text_position) {
    if (enc_text[text_position] == enc_text[text_position-1] && run_length < SAMPLED_RL_MAX_RUN_LENGTH) {
      ++run_length; // Add RL-counter
    } else {
      enc_rl_text[rl_text_position] = enc_text[text_position]; // Add character
      run_length = 1; // Reset RL-counter
      // Store RL samples
      if (rl_text_position%SAMPLED_RL_SAMPLING_RATE==0) {
        sampled_rl_sample(archive_builder->sampled_rl,num_rl_samples++,text_position);
      }
      ++rl_text_position;
    }
  }
  // Print last run
  dna_text_set_length(archive_builder->enc_rl_text,rl_text_position);
  // Swap text with rl-text (as to SA-index the rl-text)
  SWAP(archive_builder->enc_text,archive_builder->enc_rl_text);
}
/*
 * Display
 */
GEM_INLINE void archive_builder_text_dump(archive_builder_t* const archive_builder,const char* const extension) {
  // Open file
  char* const indexed_text_file_name = gem_strcat(archive_builder->output_file_name_prefix,extension);
  FILE* const indexed_text_file = fopen(indexed_text_file_name,"w");
  dna_text_pretty_print_content(indexed_text_file,archive_builder->enc_text,80);
  // Close & release
  fclose(indexed_text_file);
  free(indexed_text_file_name);
}
