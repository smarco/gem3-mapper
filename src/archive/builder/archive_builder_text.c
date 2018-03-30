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
#include "archive/builder/archive_builder_text_parser.h"
#include "archive/archive_text_rl.h"

/*
 * Generate DNA-Text
 */
void archive_builder_text_process(
    archive_builder_t* const archive_builder,
    const bool verbose) {
  // Inspect Text
  archive_builder_inspect_text(archive_builder,verbose);
  // Generate Text (Forward)
  archive_builder_generate_forward_text(archive_builder,verbose);
  // Generate C2T & G2A Texts
  if (archive_builder->type==archive_dna_bisulfite) {
    archive_builder_generate_bisulfite_text(archive_builder,verbose);
  }
  // Set forward-text length
  archive_builder->forward_text_length = archive_builder->parsing_state.index_position;
  // Generate RC-Text
  archive_builder_generate_rc_text(archive_builder,verbose);
  archive_builder_generate_text_add_separator(archive_builder); // Add extra separator (Close full-text)
  // Set full-text length
  dna_text_set_length(archive_builder->enc_text,archive_builder->parsing_state.index_position);
}
/*
 * Run-length Text (Apply RL to the text)
 */
void archive_builder_text_apply_run_length(
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
