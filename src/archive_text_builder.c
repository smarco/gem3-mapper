/*
 * PROJECT: GEMMapper
 * FILE: archive_builder.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_builder.h"
#include "archive_text_builder.h"

/*
 * Constants
 */
#define ARCHIVE_BUILDER_RL_MAX_RUN_LENGTH 5
#define ARCHIVE_BUILDER_RL_SAMPLING_RATE  (1<<6) /* 2^6=64 (HG => 460 MB)*/

/*
 * DEBUG
 */
GEM_INLINE void archive_builder_dump_index_text(archive_builder_t* const archive_builder) {
  // Open file
  char* const indexed_text_file_name = gem_strcat(archive_builder->output_file_name_prefix,".text");
  FILE* const indexed_text_file = fopen(indexed_text_file_name,"w");
  dna_text_builder_pretty_print_content(indexed_text_file,archive_builder->enc_text,80);
  // Close & release
  fclose(indexed_text_file);
  free(indexed_text_file_name);
}
/*
 * Archive Builder. Inspect Text
 */
GEM_INLINE void archive_builder_inspect_text(
    archive_builder_t* const archive_builder,input_file_t* const input_multifasta,const bool verbose) {
  // Prepare ticket
  ticker_t ticker;
  ticker_percentage_reset(&ticker,verbose,"Inspecting MultiFASTA",0,0,true);
  // MultiFASTA Read cycle
  vector_t* const line_buffer = vector_new(200,char);
  uint64_t enc_text_length = 0; // input_file_get_size(input_multifasta);
  while (!input_file_eof(input_multifasta)) {
    // Get line
    input_file_get_lines(input_multifasta,line_buffer,1);
    const char line_begin = *vector_get_mem(line_buffer,char);
    if (gem_expect_true(line_begin!=FASTA_TAG_BEGIN)) {
      enc_text_length += vector_get_used(line_buffer)-1;
    } else {
      ++enc_text_length; // Separator
    }
  }
  ++enc_text_length; // Separator
  vector_delete(line_buffer);
  ticker_finish(&ticker);
  // Configure RC generation
  if (archive_builder->indexed_complement == index_complement_auto) {
    archive_builder->indexed_complement =
        (enc_text_length <= archive_builder->complement_size_threshold) ? index_complement_yes : index_complement_no;
  }
  if (archive_builder->indexed_complement == index_complement_yes) {
    enc_text_length = 2*enc_text_length + 1; // Add complement length
  } else {
    ++enc_text_length; // Add extra separator (Close text)
  }
  // Rewind input MULTIFASTA
  input_file_rewind(input_multifasta);
  input_multifasta_state_clear(&(archive_builder->parsing_state));
  // Log
  gem_info("Inspected text %lu characters (%s). Requesting %lu MB (enc_text)",enc_text_length,
      (archive_builder->indexed_complement==index_complement_yes) ? "index_complement=yes" : "index_complement=no",
      CONVERT_B_TO_MB(enc_text_length));
  // Allocate Text (Circular BWT extra)
  archive_builder->enc_text = dna_text_builder_padded_new(enc_text_length,2,SA_BWT_PADDED_LENGTH);
}
/*
 * Archive Builder. Text Generation Low-Level Building-Blocks
 */
GEM_INLINE void archive_builder_generate_text_add_character(archive_builder_t* const archive_builder,const uint8_t char_enc) {
  // Update parsing-location
  ++(archive_builder->parsing_state.text_interval_length);
  ++(archive_builder->parsing_state.index_interval_length);
  // Add to Index-Text
  dna_text_builder_set_char(archive_builder->enc_text,(archive_builder->parsing_state.index_position)++,char_enc);
}
GEM_INLINE void archive_builder_generate_text_filter__add_character(archive_builder_t* const archive_builder,const uint8_t char_enc) {
  uint8_t filtered_char_enc = char_enc;
  // Check colorspace
  if (gem_expect_false(archive_builder->filter_type == Iupac_colorspace_dna)) {
    filtered_char_enc = dna_encoded_colorspace(archive_builder->parsing_state.last_char,char_enc);
    archive_builder->parsing_state.last_char = char_enc;
  }
  archive_builder_generate_text_add_character(archive_builder,filtered_char_enc);
}
GEM_INLINE void archive_builder_generate_text_add_separator(archive_builder_t* const archive_builder) {
  // Add to Index-Text
  archive_builder->parsing_state.last_char = ENC_DNA_CHAR_SEP;
  dna_text_builder_set_char(archive_builder->enc_text,(archive_builder->parsing_state.index_position)++,ENC_DNA_CHAR_SEP);
  locator_builder_skip_index(archive_builder->locator,1); // Skip Separator
}
GEM_INLINE void archive_builder_generate_text_add_Ns(archive_builder_t* const archive_builder) {
  // Below threshold, restore all Ns
  input_multifasta_state_t* const parsing_state = &(archive_builder->parsing_state);
  uint64_t i;
  for (i=0;i<parsing_state->ns_pending;++i) {
    archive_builder_generate_text_filter__add_character(archive_builder,ENC_DNA_CHAR_N);
  }
  parsing_state->ns_pending = 0;
}
/*
 * Archive Builder. Text Generation High-Level Building-Blocks
 */
GEM_INLINE void archive_builder_generate_text_close_sequence(archive_builder_t* const archive_builder) {
  locator_builder_t* const locator = archive_builder->locator;
  input_multifasta_state_t* const parsing_state = &(archive_builder->parsing_state);
  if (parsing_state->index_interval_length > 0) {
    // Close interval
    locator_builder_close_interval(locator,
        parsing_state->text_interval_length,
        parsing_state->index_interval_length,locator_interval_regular);
    // Close sequence (Add 1 separator)
    archive_builder_generate_text_add_separator(archive_builder);
    // Open new interval
    locator_builder_open_interval(locator,parsing_state->tag_id);
  }
  // Reset length
  input_multifasta_state_reset_interval(parsing_state);
}
GEM_INLINE void archive_builder_generate_text_process_unknowns(archive_builder_t* const archive_builder) {
  input_multifasta_state_t* const parsing_state = &(archive_builder->parsing_state);
  // Check Ns
  const uint64_t ns_pending = parsing_state->ns_pending;
  if (ns_pending > 0) {
    if (ns_pending < archive_builder->ns_threshold) {
      // Add Explicit Ns
      archive_builder_generate_text_add_Ns(archive_builder);
    } else {
      if (parsing_state->index_interval_length > 0) {
        // Close interval
        locator_builder_close_interval(archive_builder->locator,
            parsing_state->text_interval_length,
            parsing_state->index_interval_length,locator_interval_regular);
        // Close sequence (Add 1 separator)
        archive_builder_generate_text_add_separator(archive_builder);
        // Open new interval
        locator_builder_open_interval(archive_builder->locator,parsing_state->tag_id);
      }
      // Add Ns Interval
      locator_builder_close_interval(archive_builder->locator,ns_pending,0,locator_interval_unknown);
      // Open new interval
      locator_builder_open_interval(archive_builder->locator,parsing_state->tag_id);
      // Reset
      parsing_state->ns_pending = 0;
      input_multifasta_state_reset_interval(parsing_state);
    }
  }
}
GEM_INLINE void archive_builder_generate_text_add_sequence(
    archive_builder_t* const archive_builder,
    input_file_t* const input_multifasta,vector_t* const tag) {
  // Close sequence
  if (archive_builder->parsing_state.multifasta_read_state==Reading_sequence) {
    gem_cond_fatal_error(input_multifasta_get_text_sequence_length(&archive_builder->parsing_state)==0,
        MULTIFASTA_SEQ_EMPTY,PRI_input_file_content(input_multifasta));
    archive_builder_generate_text_process_unknowns(archive_builder); // Check Ns
    archive_builder_generate_text_close_sequence(archive_builder);   // Close last sequence
  }
  // Parse TAG (Skip separators)
  const uint64_t tag_buffer_length = vector_get_used(tag);
  char* const tag_buffer = vector_get_mem(tag,char)+1;
  uint64_t tag_length;
  for (tag_length=0;tag_length<tag_buffer_length;++tag_length) {
    if (MFASTA_IS_ANY_TAG_SEPARATOR(tag_buffer[tag_length])) break;
  }
  tag_buffer[tag_length] = EOS;
  gem_cond_fatal_error(tag_length==0,MULTIFASTA_TAG_EMPTY,PRI_input_file_content(input_multifasta));
  // Add to locator
  const int64_t tag_id = locator_builder_add_sequence(archive_builder->locator,tag_buffer,tag_length);
  // Open interval
  locator_builder_open_interval(archive_builder->locator,tag_id);
  // Begin new text-sequence (Expect sequence after TAG)
  input_multifasta_state_begin_sequence(&archive_builder->parsing_state);
}
GEM_INLINE void archive_builder_generate_text_process_character(
    archive_builder_t* const archive_builder,input_file_t* const input_multifasta,const char current_char) {
  // Check Character
  if (current_char==DNA_CHAR_N || !is_extended_dna(current_char)) { // Handle Ns
    gem_cond_fatal_error(!is_iupac_code(current_char),MULTIFASTA_INVALID_CHAR,PRI_input_file_content(input_multifasta),current_char);
    ++(archive_builder->parsing_state.ns_pending);
  } else { // Other character
    // Handle pending Ns
    archive_builder_generate_text_process_unknowns(archive_builder);
    // Add character
    archive_builder_generate_text_filter__add_character(archive_builder,dna_encode(current_char));
  }
}
/*
 * Archive Builder. Generate Text & RC-Text
 */
GEM_INLINE uint64_t archive_builder_generate_text(
    archive_builder_t* const archive_builder,input_file_t* const input_multifasta,const bool verbose) {
  // Check MultiFASTA
  input_file_check_buffer(input_multifasta);
  gem_cond_fatal_error((char)input_file_get_current_char(input_multifasta)!=FASTA_TAG_BEGIN,
      MULTIFASTA_BEGINNING_TAG,PRI_input_file_content(input_multifasta));
  // Prepare ticket
  ticker_t ticker;
  ticker_count_reset(&ticker,verbose,"Reading MultiFASTA",0,10000000,true);
  ticker_add_process_label(&ticker,"","bases parsed");
  ticker_add_finish_label(&ticker,"Total","bases parsed");
  // MultiFASTA Read cycle
  vector_t* const line_buffer = vector_new(200,char);
  input_multifasta_state_t* const parsing_state = &(archive_builder->parsing_state);
  while (!input_file_eof(input_multifasta)) {
    // Get line
    input_file_get_lines(input_multifasta,line_buffer,1);
    const uint64_t line_length = vector_get_used(line_buffer)-1;
    const char* const line = vector_get_mem(line_buffer,char);
    // Parse line
    if (line[0]==FASTA_TAG_BEGIN) { // Archive builder parse Tag
      archive_builder_generate_text_add_sequence(archive_builder,input_multifasta,line_buffer);
    } else { // Archive builder parse content
      parsing_state->multifasta_read_state = Reading_sequence;
      // Process characters
      uint64_t i;
      for (i=0;i<line_length;++i) {
        archive_builder_generate_text_process_character(archive_builder,input_multifasta,line[i]);
        ++(parsing_state->text_position); // Inc Text Position
      }
      ticker_update(&ticker,line_length);
    }
  }
  // Close sequence
  archive_builder_generate_text_close_sequence(archive_builder);
  gem_cond_fatal_error(input_multifasta_get_text_sequence_length(&archive_builder->parsing_state)==0,
      MULTIFASTA_SEQ_EMPTY,PRI_input_file_content(input_multifasta)); // Check sequence not null // FIXME
  // Free
  vector_delete(line_buffer);
  // Ticker banner
  ticker_finish(&ticker);
  // Return the length of the text
  return archive_builder->parsing_state.index_position;
}
GEM_INLINE void archive_builder_generate_rc_text(archive_builder_t* const archive_builder,const bool verbose) {
  // Prepare ticker
  ticker_t ticker_rc;
  ticker_percentage_reset(&ticker_rc,verbose,"Generating Text (explicit Reverse-Complement)",0,100,true);
  // Traverse all reference intervals (num_base_intervals are those from the graph file; not RC)
  locator_builder_t* const locator = archive_builder->locator;
  const uint64_t num_intervals = locator_builder_get_num_intervals(archive_builder->locator);
  uint64_t i;
  for (i=0;i<num_intervals;++i) {
    // Retrieve interval
    locator_interval_t* const locator_interval = locator_builder_get_interval(archive_builder->locator,i);
    const uint64_t interval_length = locator_interval_get_index_length(locator_interval);
    // Add RC interval to locator
    locator_builder_add_rc_interval(locator,locator_interval);
    // Generate RC-text
    uint64_t i = 0, text_position = locator_interval->end_position-1;
    for (i=0;i<interval_length;++i) {
      // Add character (k-mer counting) [Check Colorspace]
      const uint8_t filtered_char_enc =
          (gem_expect_false(archive_builder->filter_type == Iupac_colorspace_dna)) ?
              dna_text_builder_get_char(archive_builder->enc_text,text_position) :
              dna_encoded_complement(dna_text_builder_get_char(archive_builder->enc_text,text_position));
      archive_builder_generate_text_add_character(archive_builder,filtered_char_enc);
      --text_position; // Next
    }
    // Add Separator
    archive_builder_generate_text_add_separator(archive_builder);
  }
  // Add extra separator (Close text)
  archive_builder_generate_text_add_separator(archive_builder);
  ticker_finish(&ticker_rc);
}
/*
 * STEP1 Archive Build :: Process MultiFASTA file
 *   1. MultiFASTA Read cycle
 *     1.1 Filter UIPAC-DNA bases
 *     1.2 Strip Ns
 *   2. Generate Locator
 *   3. Generate Index-Text
 *   4. Write (1/4) :: Header & Locator
 */
GEM_INLINE void archive_builder_process_multifasta(
    archive_builder_t* const archive_builder,input_file_t* const input_multifasta,
    const bool dump_locator_intervals,const bool dump_indexed_text,const bool verbose) {
  /*
   * Inspect Text
   */
  archive_builder_inspect_text(archive_builder,input_multifasta,verbose);
  /*
   * Generate Text
   */
  // Generate Text (Forward)
  archive_builder_generate_text(archive_builder,input_multifasta,verbose);
  // Generate RC-Text (Reverse)
  if (archive_builder->indexed_complement == index_complement_yes) {
    archive_builder_generate_rc_text(archive_builder,verbose);
  } else {
    archive_builder_generate_text_add_separator(archive_builder);
  }
  // Set the precise text-length
  dna_text_builder_set_length(archive_builder->enc_text,archive_builder->parsing_state.index_position);
  /*
   * DEBUG
   */
  // DEBUG locator
  locator_builder_print(gem_info_get_stream(),archive_builder->locator,dump_locator_intervals); // Locator
  // DEBUG index_text
  if (dump_indexed_text) archive_builder_dump_index_text(archive_builder);
  /*
   * Write (1/3) :: Header & Locator
   */
  archive_builder_write_header(archive_builder);
  archive_builder_write_locator(archive_builder);
  locator_builder_delete(archive_builder->locator); // Free Locator
}
/*
 * STEP1 Archive Build :: Process MultiFASTA file
 *  ...
 *   6. Generate the RL-text
 */
GEM_INLINE void archive_builder_process_run_length_text(
    archive_builder_t* const archive_builder,const bool dump_run_length_text,const bool verbose) {
  // Allocate RL-text (Circular BWT extra)
  const uint64_t enc_text_length = dna_text_builder_get_length(archive_builder->enc_text);
  archive_builder->enc_rl_text = dna_text_builder_padded_new(enc_text_length,2,SA_BWT_PADDED_LENGTH);
  const uint64_t max_num_samples = DIV_CEIL(enc_text_length,ARCHIVE_BUILDER_RL_SAMPLING_RATE);
  uint64_t* const sampled_rl_text = mm_calloc(max_num_samples,uint64_t,true);
  archive_builder->sampled_rl_text = sampled_rl_text;
  // Compact the text into the RL-text
  const uint8_t* const enc_text = dna_text_builder_get_text(archive_builder->enc_text);
  uint8_t* const enc_rl_text = dna_text_builder_get_text(archive_builder->enc_rl_text);
  uint64_t text_position, rl_text_position;
  uint64_t run_length=1, num_rl_samples=0;
  uint8_t last_char_enc = enc_text[0]; // Get first character
  for (text_position=1,rl_text_position=0;text_position<enc_text_length;++text_position) {
    const uint8_t char_enc = enc_text[text_position]; // Get character
    if (char_enc == last_char_enc) {
      if (gem_expect_false(run_length == ARCHIVE_BUILDER_RL_MAX_RUN_LENGTH)) {
        // Store RL samples
        if (rl_text_position%ARCHIVE_BUILDER_RL_SAMPLING_RATE==0) {
          sampled_rl_text[num_rl_samples++] = text_position;
        }
        // Add character
        enc_rl_text[rl_text_position++] = last_char_enc;
        // Reset RL-counter
        run_length = 1;
      } else {
        ++run_length;
      }
    } else {
      // Store RL samples
      if (rl_text_position%ARCHIVE_BUILDER_RL_SAMPLING_RATE==0) {
        sampled_rl_text[num_rl_samples++] = text_position;
      }
      // Add character
      enc_rl_text[rl_text_position++] = last_char_enc;
      // Reset RL-counter
      run_length = 1;
      last_char_enc = char_enc;
    }
  }
  // Print last run
  enc_rl_text[rl_text_position++] = last_char_enc;
  dna_text_builder_set_length(archive_builder->enc_rl_text,rl_text_position);
  // Swap text with rl-text (as to SA-index the rl-text)
  SWAP(archive_builder->enc_text,archive_builder->enc_rl_text);
}
