/*
 * PROJECT: GEMMapper
 * FILE: archive_builder.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_builder.h"
#include "archive_text_builder.h"

/*
 * DEBUG
 */
GEM_INLINE void archive_builder_dump_index_text(archive_builder_t* const archive_builder) {
  // Open file
  char* const indexed_text_file_name = gem_strcat(archive_builder->output_file_name_prefix,".text");
  FILE* const indexed_text_file = fopen(indexed_text_file_name,"w");
  dna_text_pretty_print_content(indexed_text_file,archive_builder->enc_text,80);
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
    input_file_get_line(input_multifasta,line_buffer);
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
  gem_log("Inspected text %lu characters (%s). Requesting %lu MB (enc_text)",enc_text_length,
      (archive_builder->indexed_complement==index_complement_yes) ? "index_complement=yes" : "index_complement=no",
      CONVERT_B_TO_MB(enc_text_length));
  // Allocate Text (Circular BWT extra)
  archive_builder->enc_text = dna_text_padded_new(enc_text_length,2,SA_BWT_PADDED_LENGTH);
}
/*
 * Archive Builder. Text Generation Low-Level Building-Blocks
 */
GEM_INLINE void archive_builder_generate_text_add_character(archive_builder_t* const archive_builder,const uint8_t char_enc) {
  // Update parsing-location
  ++(archive_builder->parsing_state.text_interval_length);
  ++(archive_builder->parsing_state.index_interval_length);
  // Add to Index-Text
  dna_text_set_char(archive_builder->enc_text,(archive_builder->parsing_state.index_position)++,char_enc);
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
  dna_text_set_char(archive_builder->enc_text,(archive_builder->parsing_state.index_position)++,ENC_DNA_CHAR_SEP);
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
    locator_builder_skip_index(locator,1);
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
        locator_builder_skip_index(archive_builder->locator,1);
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
    input_file_get_line(input_multifasta,line_buffer);
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
  const uint64_t num_intervals = locator_builder_get_num_intervals(archive_builder->locator);
  uint64_t i;
  for (i=0;i<num_intervals;++i) {
    // Retrieve interval
    locator_interval_t* const locator_interval = locator_builder_get_interval(archive_builder->locator,i);
    const uint64_t interval_length = locator_interval_get_index_length(locator_interval);
    // Add RC interval to locator
    locator_builder_add_rc_interval(archive_builder->locator,locator_interval);
    // Generate RC-text
    uint64_t i = 0, text_position = locator_interval->end_position-1;
    for (i=0;i<interval_length;++i) {
      // Add character (k-mer counting) [Check Colorspace]
      const uint8_t filtered_char_enc =
          (gem_expect_false(archive_builder->filter_type == Iupac_colorspace_dna)) ?
              dna_text_get_char(archive_builder->enc_text,text_position) :
              dna_encoded_complement(dna_text_get_char(archive_builder->enc_text,text_position));
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
  dna_text_set_length(archive_builder->enc_text,archive_builder->parsing_state.index_position);
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
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
//  // Allocate RL-text
//  archive_builder->rl_text = cdna_text_new(mm_pool_get_slab(mm_pool_32MB));
//  // Init SA-Builder
//  sa_builder_count_begin(archive_builder->sa_builder);
//  // Init iterator & run-counters
//  cdna_text_iterator_t iterator;
//  cdna_text_iterator_init(&iterator,archive_builder->rl_text,0);
//  gem_cond_fatal_error(cdna_text_iterator_eoi(&iterator),
//      ARCHIVE_BUILDER_SEQUENCE_MIN_LENGTH,0ul,archive_builder->sa_builder->kmer_length);
//  uint8_t last_char_enc = dna_decode(cdna_text_iterator_get_char_encoded(&iterator)); // Get character
//  uint64_t run_length = 1, text_position = 0;
//  // Run-Length Compact the text
//  while (!cdna_text_iterator_eoi(&iterator)) {
//    ++text_position;
//    const uint8_t char_enc = dna_decode(cdna_text_iterator_get_char_encoded(&iterator)); // Get character
//    if (char_enc == last_char_enc) {
//      ++run_length;
//    } else {
//      // Add character
//      sa_builder_count_suffix(archive_builder->sa_builder,last_char_enc); // Add to SA-Builder
//      cdna_text_add_char(archive_builder->rl_text,last_char_enc); // Add to Index-Text
//      cdna_text_iterator_next_char(&iterator); // Next
//      // Reset RL-counter
//      last_char_enc = char_enc;
//      run_length = 1;
//    }
//
//    // TODO
//    // Store RL samples
//
//  }
//  // Print last run
//  sa_builder_count_suffix(archive_builder->sa_builder,last_char_enc); // Add to SA-Builder
//  cdna_text_add_char(archive_builder->rl_text,last_char_enc); // Add to Index-Text
//  cdna_text_iterator_next_char(&iterator); // Next
//  // Finish counting k-mers & Close text
//  sa_builder_count_end(archive_builder->sa_builder);
//  cdna_text_close(archive_builder->rl_text);
//  // Check Minimum bases processed
//  const uint64_t kmer_length = archive_builder->sa_builder->kmer_length;
//  const uint64_t rl_index_length = cdna_text_get_length(archive_builder->rl_text);
//  gem_cond_fatal_error(rl_index_length < kmer_length,ARCHIVE_BUILDER_SEQUENCE_MIN_LENGTH,rl_index_length,kmer_length);
//  // Swap text with rl-text (as to SA-index the rl-text)
//  SWAP(archive_builder->text,archive_builder->rl_text);
}
