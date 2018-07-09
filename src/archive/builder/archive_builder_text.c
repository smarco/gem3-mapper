/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Archive builder module to parse a MFASTA (multi-FASTA) input file
 *   and generate the text-index (raw text index) and meta-info
 *     1. MultiFASTA Read cycle
 *       1.1 Filter UIPAC-DNA bases
 *       1.2 Strip Ns
 *     2. Generate Locator
 *     3. Generate Index-Text
 *     4. Write Header & Locator
 */

#include "text/dna_text.h"
#include "archive/builder/archive_builder_text.h"
#include "archive/archive_text_rl.h"

/*
 * Archive Builder. Text Generation Low-Level Building-Blocks
 */
void archive_builder_text_add_character(
    archive_builder_t* const archive_builder,
    const uint8_t char_enc) {
  // Parameters
  input_multifasta_file_t* const input_multifasta_file = &archive_builder->input_multifasta_file;
  input_multifasta_state_t* const parsing_state = &input_multifasta_file->parsing_state;
  // Update parsing-location
  ++(parsing_state->text_interval_length);
  ++(parsing_state->index_interval_length);
  // Add to Index-Text
  dna_text_set_char(archive_builder->enc_text,(parsing_state->index_position)++,char_enc);
}
void archive_builder_text_filter__add_character(
    archive_builder_t* const archive_builder,
    const uint8_t char_enc) {
  uint8_t filtered_char_enc = char_enc;
  archive_builder_text_add_character(archive_builder,filtered_char_enc);
}
void archive_builder_text_add_separator(archive_builder_t* const archive_builder) {
  // Parameters
  input_multifasta_file_t* const input_multifasta_file = &archive_builder->input_multifasta_file;
  input_multifasta_state_t* const parsing_state = &input_multifasta_file->parsing_state;
  // Add to Index-Text
  parsing_state->last_char = ENC_DNA_CHAR_SEP;
  dna_text_set_char(archive_builder->enc_text,(parsing_state->index_position)++,ENC_DNA_CHAR_SEP);
  locator_builder_skip_index(archive_builder->locator,1); // Skip Separator
}
void archive_builder_text_add_Ns(
    archive_builder_t* const archive_builder) {
  // Parameters
  input_multifasta_file_t* const input_multifasta_file = &archive_builder->input_multifasta_file;
  input_multifasta_state_t* const parsing_state = &input_multifasta_file->parsing_state;
  // Below threshold, restore all Ns
  uint64_t i;
  for (i=0;i<parsing_state->ns_pending;++i) {
    archive_builder_text_filter__add_character(archive_builder,ENC_DNA_CHAR_N);
  }
  parsing_state->ns_pending = 0;
}
/*
 * Text Generation. High-Level Building-Blocks
 */
void archive_builder_text_close_sequence(
    archive_builder_t* const archive_builder) {
  // Parameters
  input_multifasta_file_t* const input_multifasta_file = &archive_builder->input_multifasta_file;
  input_multifasta_state_t* const parsing_state = &input_multifasta_file->parsing_state;
  locator_builder_t* const locator = archive_builder->locator;
  if (parsing_state->index_interval_length > 0) {
    // Close interval
    locator_builder_close_interval(locator,parsing_state->text_interval_length,
        parsing_state->index_interval_length,parsing_state->interval_type,
        parsing_state->strand,parsing_state->bs_strand);
    // Close sequence (Add 1 separator)
    archive_builder_text_add_separator(archive_builder);
    // Open new interval
    locator_builder_open_interval(locator,parsing_state->tag_id);
  }
  // Reset length
  input_multifasta_state_reset_interval(parsing_state);
}
void archive_builder_text_process_unknowns(
    archive_builder_t* const archive_builder) {
  // Parameters
  input_multifasta_file_t* const input_multifasta_file = &archive_builder->input_multifasta_file;
  input_multifasta_state_t* const parsing_state = &input_multifasta_file->parsing_state;
  // Check Ns
  const uint64_t ns_pending = parsing_state->ns_pending;
  if (ns_pending > 0) {
    if (ns_pending < archive_builder->ns_threshold) {
      // Add Explicit Ns
      archive_builder_text_add_Ns(archive_builder);
    } else {
      if (parsing_state->index_interval_length > 0) {
        // Close interval
        locator_builder_close_interval(
            archive_builder->locator,parsing_state->text_interval_length,
            parsing_state->index_interval_length,parsing_state->interval_type,
            parsing_state->strand,parsing_state->bs_strand);
        // Close sequence (Add 1 separator)
        archive_builder_text_add_separator(archive_builder);
        // Open new interval
        locator_builder_open_interval(archive_builder->locator,parsing_state->tag_id);
      }
      // Add Ns Interval
      locator_builder_close_interval(archive_builder->locator,
          ns_pending,0,locator_interval_uncalled,Forward,bs_strand_none);
      // Open new interval
      locator_builder_open_interval(archive_builder->locator,parsing_state->tag_id);
      // Reset
      parsing_state->ns_pending = 0;
      input_multifasta_state_reset_interval(parsing_state);
    }
  }
}
void archive_builder_text_add_sequence(
    archive_builder_t* const archive_builder,
    char* line_buffer,
    int line_length) {
  // Parameters
  input_multifasta_file_t* const input_multifasta_file = &archive_builder->input_multifasta_file;
  input_multifasta_state_t* const parsing_state = &input_multifasta_file->parsing_state;
  // Close sequence
  if (parsing_state->multifasta_read_state==Reading_sequence) {
    gem_cond_error(
        input_multifasta_get_text_sequence_length(parsing_state)==0,MULTIFASTA_SEQ_EMPTY,
        input_multifasta_file->file_name,input_multifasta_file->line_no);
    archive_builder_text_process_unknowns(archive_builder); // Check Ns
    archive_builder_text_close_sequence(archive_builder);   // Close last sequence
  }
  // Parse TAG (Skip separators & check characters)
  while ((line_length > 0) && (line_buffer[0]==SPACE || line_buffer[0]==TAB)) {
    ++line_buffer; --line_length;
  }
  uint64_t tag_length = 0;
  while ((line_length > 0) && MFASTA_VALID_CHARACTER(line_buffer[tag_length])) {
    ++tag_length; --line_length;
  }
  line_buffer[tag_length] = EOS;
  gem_cond_error(
      tag_length==0,MULTIFASTA_TAG_EMPTY,
      input_multifasta_file->file_name,
      input_multifasta_file->line_no);
  // Add to locator
  parsing_state->tag_id = locator_builder_add_sequence(archive_builder->locator,line_buffer,tag_length);
  // Open interval
  locator_builder_open_interval(archive_builder->locator,parsing_state->tag_id);
  // Begin new text-sequence (Expect sequence after TAG)
  input_multifasta_state_begin_sequence(parsing_state);
}
void archive_builder_text_process_character(
    archive_builder_t* const archive_builder,
    const char current_char) {
  // Parameters
  input_multifasta_file_t* const input_multifasta_file = &archive_builder->input_multifasta_file;
  input_multifasta_state_t* const parsing_state = &input_multifasta_file->parsing_state;
  // Check Character
  if (current_char==DNA_CHAR_N || !is_extended_dna(current_char)) { // Handle Ns
    gem_cond_error(!is_iupac_code(current_char),MULTIFASTA_INVALID_CHAR,
        input_multifasta_file->file_name,
        input_multifasta_file->line_no,current_char);
    ++(parsing_state->ns_pending);
  } else { // Other character
    // Handle pending Ns
    archive_builder_text_process_unknowns(archive_builder);
    // Add character
    archive_builder_text_filter__add_character(archive_builder,dna_encode(current_char));
  }
}
/*
 * Generate Text
 */
void archive_builder_text_generate_forward(
    archive_builder_t* const archive_builder,
    const bool verbose) {
  // Parameters
  input_multifasta_file_t* const input_multifasta_file = &archive_builder->input_multifasta_file;
  char *line_buffer = NULL;
  int line_length = 0;
  size_t line_allocated = 0;
  // Check MultiFASTA
  gem_cond_error(fm_eof(input_multifasta_file->file_manager),
      MULTIFASTA_EMPTY,input_multifasta_file->file_name);
  // Prepare ticket
  ticker_t ticker;
  ticker_count_reset(&ticker,verbose,"Reading MultiFASTA",0,100000000,true);
  ticker_add_process_label(&ticker,"","bases parsed");
  ticker_add_finish_label(&ticker,"Total","bases parsed");
  // MultiFASTA Read cycle
  input_multifasta_state_t* const parsing_state = &input_multifasta_file->parsing_state;
  input_multifasta_state_clear(parsing_state); // Init parsing state
  while (true) {
    // Get line
    line_length = fm_getline(&line_buffer,&line_allocated,input_multifasta_file->file_manager);
    ++(input_multifasta_file->line_no);
    if (line_length==-1) {
      gem_cond_error(
          !fm_eof(input_multifasta_file->file_manager),MULTIFASTA_READING,
          input_multifasta_file->file_name,input_multifasta_file->line_no);
      break;
    }
    // Check error
    gem_cond_error(
        (input_multifasta_file->line_no==1 && line_buffer[0]!=FASTA_TAG_BEGIN),
        MULTIFASTA_BEGINNING_TAG,input_multifasta_file->file_name,
        input_multifasta_file->line_no);
    // Parse line
    if (line_buffer[0]==FASTA_TAG_BEGIN) {
      // Parse Tag
      archive_builder_text_add_sequence(archive_builder,line_buffer+1,line_length-2);
    } else {
      // Process characters
      parsing_state->multifasta_read_state = Reading_sequence;
      uint64_t i;
      for (i=0;i<line_length-1;++i) {
        archive_builder_text_process_character(archive_builder,line_buffer[i]);
        ++(parsing_state->text_position); // Inc Text Position
      }
      ticker_update(&ticker,line_length);
    }
  }
  // Close sequence
  gem_cond_error(
      input_multifasta_get_text_sequence_length(&input_multifasta_file->parsing_state)==0,
      MULTIFASTA_SEQ_EMPTY,input_multifasta_file->file_name,
			input_multifasta_file->line_no); // Check sequence not null
  archive_builder_text_process_unknowns(archive_builder); // Check Ns
  archive_builder_text_close_sequence(archive_builder);
  // Free
  if (line_buffer != NULL) free(line_buffer);
  ticker_finish(&ticker);
}
/*
 * Generate RC-Text
 */
void archive_builder_text_reverse_complement(
    archive_builder_t* const archive_builder,
    const bool verbose) {
  // Prepare ticker
  ticker_t ticker_rc;
  ticker_percentage_reset(&ticker_rc,verbose,"Generating Text (explicit Reverse-Complement)",0,100,true);
  // Traverse all reference intervals (num_base_intervals are those from the graph file; not RC)
  locator_builder_t* const locator = archive_builder->locator;
  const uint64_t num_intervals = locator_builder_get_num_intervals(archive_builder->locator);
  uint64_t i;
  for (i=num_intervals;i-->0;) {
    // Retrieve interval
    locator_interval_t* const locator_interval = locator_builder_get_interval(archive_builder->locator,i);
    const uint64_t interval_length = locator_interval_get_index_length(locator_interval);
    // Add RC interval to locator
    locator_builder_add_rc_interval(locator,locator_interval);
    // Generate RC-text
    uint64_t j = 0, text_position = locator_interval->end_position-1;
    for (j=0;j<interval_length;++j) {
      const uint8_t char_enc = dna_text_get_char(archive_builder->enc_text,text_position);
      const uint8_t filtered_char_enc = dna_encoded_complement(char_enc);
      archive_builder_text_add_character(archive_builder,filtered_char_enc); // Add character
      --text_position; // Next
    }
    // Add Separator
    if (locator_interval->end_position-locator_interval->begin_position > 0) {
      archive_builder_text_add_separator(archive_builder);
    }
  }
  ticker_finish(&ticker_rc);
}
/*
 * Generate C2T & G2A Texts
 */
void archive_builder_text_generate_bisulfite(
    archive_builder_t* const archive_builder,const bool verbose) {
  // Prepare ticker
  ticker_t ticker_rc;
  ticker_percentage_reset(&ticker_rc,verbose,"Generating Bisulfite-Text (C2T & G2A)",0,100,true);
  // Traverse all reference intervals (num_base_intervals are those from the graph file; not RC)
  locator_builder_t* const locator = archive_builder->locator;
  const uint64_t num_intervals = locator_builder_get_num_intervals(archive_builder->locator);
  uint64_t i;
  for (i=0;i<num_intervals;++i) {
    // Retrieve interval & set C2T-stranded
    locator_interval_t* const locator_interval = locator_builder_get_interval(archive_builder->locator,i);
    locator_interval->bs_strand = bs_strand_C2T;
    // Add G2A-text interval to locator
    locator_builder_add_interval(locator,locator_interval->tag_id,
        locator_interval->sequence_offset,locator_interval->sequence_length,
        locator_interval->end_position-locator_interval->begin_position,
        locator_interval->type,locator_interval->strand,bs_strand_G2A);
    // Add G2A-text & convert C2T initial text
    const uint64_t interval_length = locator_interval_get_index_length(locator_interval);
    uint64_t j, text_position = locator_interval->begin_position;
    for (j=0;j<interval_length;++j) {
      // Add G2A-char
      const uint8_t char_enc = dna_text_get_char(archive_builder->enc_text,text_position);
      const uint8_t char_enc_G2A = dna_encoded_bisulfite_G2A(char_enc);
      archive_builder_text_add_character(archive_builder,char_enc_G2A); // Add character
      // Convert C2T-char
      const uint8_t char_enc_C2T = dna_encoded_bisulfite_C2T(char_enc);
      dna_text_set_char(archive_builder->enc_text,text_position,char_enc_C2T);
      // Next
      ++text_position;
    }
    // Add Separator
    if (locator_interval->end_position-locator_interval->begin_position > 0) {
      archive_builder_text_add_separator(archive_builder);
    }
  }
  ticker_finish(&ticker_rc);
}
/*
 * Inspect Text
 */
void archive_builder_text_inspect(
    archive_builder_t* const archive_builder,
    const bool verbose) {
  // Parameters
  input_multifasta_file_t* const input_multifasta_file = &archive_builder->input_multifasta_file;
  char *line_buffer = NULL;
  int line_length = 0;
  size_t line_allocated = 0;
  ticker_t ticker;
  // MultiFASTA Read cycle
  ticker_percentage_reset(&ticker,verbose,"Inspecting MultiFASTA",0,0,true); // Prepare ticket
  uint64_t enc_text_length = 0;
  while (true) {
    // Get line
    line_length = fm_getline(&line_buffer,&line_allocated,input_multifasta_file->file_manager);
    if (line_length == -1) break;
    // Account the line length
    if (line_buffer[0] != FASTA_TAG_BEGIN) {
      enc_text_length += line_length-1;
    } else {
      ++enc_text_length; // Separator
    }
  }
  ++enc_text_length; // Separator
  ticker_finish(&ticker);
  if (line_buffer != NULL) free(line_buffer); // Free
  // Configure RC generation
  if (archive_builder->type!=archive_dna_forward) {
    enc_text_length = 2*enc_text_length; // Add complement length
  }
  // Configure Bisulfite generation
  if (archive_builder->type==archive_dna_bisulfite) {
    enc_text_length = 2*enc_text_length;
  }
  ++enc_text_length; // Add extra separator (Close text)
  // Rewind input MULTIFASTA
  fm_seek(input_multifasta_file->file_manager,0);
  input_multifasta_file->line_no = 0;
  // Log
  tfprintf(gem_log_get_stream(),
      "Inspected text %"PRIu64" characters (index_complement=%s). Requesting %"PRIu64" MB (encoded text)\n",
      enc_text_length,
      (archive_builder->type==archive_dna_forward) ? "no" : "yes",
      CONVERT_B_TO_MB(enc_text_length));
  // Allocate Text (Circular BWT extra)
  archive_builder->enc_text = dna_text_padded_new(enc_text_length,2,SA_BWT_PADDED_LENGTH);
}
/*
 * Generate DNA-Text
 */
void archive_builder_text_process(
    archive_builder_t* const archive_builder,
    const bool verbose) {
  // Inspect Text
  archive_builder_text_inspect(archive_builder,verbose);
  // Generate Text (Forward)
  archive_builder_text_generate_forward(archive_builder,verbose);
  // Generate C2T & G2A Texts
  if (archive_builder->type==archive_dna_bisulfite) {
    archive_builder_text_generate_bisulfite(archive_builder,verbose);
  }
  // Set forward-text length
  input_multifasta_state_t* const parsing_state = &archive_builder->input_multifasta_file.parsing_state;
  archive_builder->forward_text_length = parsing_state->index_position;
  // Generate RC-Text
  if (archive_builder->type!=archive_dna_forward) {
    archive_builder_text_reverse_complement(archive_builder,verbose);
  }
  archive_builder_text_add_separator(archive_builder); // Add extra separator (Close full-text)
  // Set full-text length
  dna_text_set_length(archive_builder->enc_text,parsing_state->index_position);
}
/*
 * Run-length Text (Apply RL to the text)
 */
void archive_builder_text_generate_run_length(
    archive_builder_t* const archive_builder,
    const bool verbose) {
  // Allocate RL-text (Circular BWT extra)
  const uint64_t enc_text_length = dna_text_get_length(archive_builder->enc_text);
  archive_builder->enc_rl_text = dna_text_padded_new(enc_text_length,2,SA_BWT_PADDED_LENGTH);
  const uint64_t max_num_samples = DIV_CEIL(enc_text_length,SAMPLED_RL_SAMPLING_RATE);
  archive_builder->sampled_rl = sampled_rl_new(SAMPLED_RL_SAMPLING_RATE,max_num_samples,enc_text_length);
  // Locator Iterator
  uint64_t rl_begin_position = 0;
  svector_iterator_t intervals_iterator;
  svector_iterator_new(&intervals_iterator,archive_builder->locator->intervals,SVECTOR_READ_ITERATOR,0);
  // Compact the text into the RL-text
  const uint8_t* const enc_text = dna_text_get_text(archive_builder->enc_text);
  uint8_t* const enc_rl_text = dna_text_get_text(archive_builder->enc_rl_text);
  uint64_t text_position, rl_text_position=1;
  uint64_t run_length=1, num_rl_samples=1;
  enc_rl_text[0] = enc_text[0]; // Add character
  sampled_rl_sample(archive_builder->sampled_rl,0,0);
  for (text_position=1;text_position<enc_text_length;++text_position) {
    if (enc_text[text_position] == enc_text[text_position-1] && run_length < TEXT_RL_MAX_RUN_LENGTH) {
      ++run_length; // Add RL-counter
    } else {
      // Check Separator
      if (enc_text[text_position]==ENC_DNA_CHAR_SEP) {
        locator_interval_t* locator_interval = svector_iterator_get_element(&intervals_iterator,locator_interval_t);
        while (locator_interval->type==locator_interval_uncalled) { // Skip uncalled
          locator_interval->rl_begin_position = rl_begin_position;
          locator_interval->rl_end_position = rl_begin_position;
          svector_read_iterator_next(&intervals_iterator);
          locator_interval = svector_iterator_get_element(&intervals_iterator,locator_interval_t);
        }
        locator_interval->rl_begin_position = rl_begin_position;
        locator_interval->rl_end_position = rl_text_position;
        rl_begin_position = rl_text_position+1;
        svector_read_iterator_next(&intervals_iterator);
      }
      // Add character & Reset RL-counter
      enc_rl_text[rl_text_position] = enc_text[text_position];
      run_length = 1;
      // Store RL samples
      if ((rl_text_position % SAMPLED_RL_SAMPLING_RATE)==0) {
        sampled_rl_sample(archive_builder->sampled_rl,num_rl_samples++,text_position);
      }
      ++rl_text_position;
    }
  }
  // Print last run
  dna_text_set_length(archive_builder->enc_rl_text,rl_text_position);
  // Close remaining intervals
  while (!svector_read_iterator_eoi(&intervals_iterator)) {
    locator_interval_t* locator_interval = svector_iterator_get_element(&intervals_iterator,locator_interval_t);
    locator_interval->rl_begin_position = rl_begin_position;
    locator_interval->rl_end_position = rl_begin_position;
    svector_read_iterator_next(&intervals_iterator);
  }
  // Swap text with rl-text (as to SA-index the rl-text)
  SWAP(archive_builder->enc_text,archive_builder->enc_rl_text);
}
/*
 * Display
 */
void archive_builder_text_dump(
    archive_builder_t* const archive_builder,
    const char* const extension) {
  // Open file
  char* const indexed_text_file_name = gem_strcat(archive_builder->output_file_name_prefix,extension);
  FILE* const indexed_text_file = fopen(indexed_text_file_name,"w");
  dna_text_pretty_print_content(indexed_text_file,archive_builder->enc_text,80);
  // Close & release
  fclose(indexed_text_file);
  free(indexed_text_file_name);
}
