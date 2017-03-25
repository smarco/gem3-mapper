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
 *   Smith-Waterman-Gotoh (SWG) alignment module
 */

#include "align/align_swg_banded.h"
#include "align/align_swg.h"
#include "matches/matches_cigar.h"

/*
 * SWG Banded
 */
void align_swg_banded_compute_column(
    const uint8_t* const key,
    const uint8_t* const text,
    swg_cell_t** const dp,
    const uint64_t column,
    const uint64_t num_rows,
    const uint64_t column_start_band,
    uint64_t* const band_high_offset,
    uint64_t* const band_low_offset,
    const swg_matching_score_t* const matching_score,
    const uint64_t single_gap,
    const uint64_t gap_extension) {
  // Initialize band boundaries & update band limits
  if (*band_low_offset < num_rows) { //  Below band
    dp[column][*band_low_offset].I = SWG_SCORE_MIN;
    dp[column-1][*band_low_offset].I = SWG_SCORE_MIN;
    ++(*band_low_offset); // Swift band
  }
  // Initialize band boundaries & update band limits
  if (column >= column_start_band) {
    ++(*band_high_offset); // Swift band
    dp[column][*band_high_offset+1].D = SWG_SCORE_MIN;
  }
  // Locate the cursor at the proper cell & calculate DP
  uint64_t row;
  for (row=*band_high_offset+1;row<*band_low_offset;++row) {
    // Update DP.M
    const uint8_t enc_text = text[column-1];
    const uint8_t enc_key = key[row-1];
    const int32_t d_value = dp[column][row].D;
    const int32_t i_value = dp[column][row].I;
    const int32_t match = dp[column-1][row-1].M + (*matching_score)[enc_text][enc_key];
    const int32_t m_value = MAX(match,MAX(d_value,i_value));
    dp[column][row].M = m_value;
    // Update DP.D
    const int32_t gap_new = m_value + single_gap;
    const int32_t del_ext = d_value + gap_extension;
    dp[column][row+1].D = MAX(gap_new,del_ext);
    // Update DP.I
    const int32_t ins_ext = i_value + gap_extension;
    dp[column+1][row].I = MAX(gap_new,ins_ext);
  }
}
void align_swg_banded(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    const bool begin_free,
    const bool end_free,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_allocator_t* const mm_allocator) {
  // Parameters
  const uint8_t* const key = align_input->key;
  const uint64_t key_length = align_input->key_length;
  uint8_t* const text = align_input->text;
  uint64_t text_length = align_input->text_length;
  const swg_penalties_t* swg_penalties = align_parameters->swg_penalties;
  const uint64_t max_bandwidth = align_parameters->max_bandwidth;
  // Initialize band-limits
  if (text_length > key_length + max_bandwidth) { // Text too long for band
    if (!begin_free && !end_free) {
      match_alignment->score = SWG_SCORE_MIN; return;
    }
    if (!begin_free) {
      text_length = key_length + max_bandwidth;
    }
  }
  if (text_length + max_bandwidth <= key_length) { // Text too short for band
    match_alignment->score = SWG_SCORE_MIN; return;
  }
  PROF_START(GP_SWG_ALIGN_BANDED);
  // Allocate memory
  mm_allocator_push_state(mm_allocator); // Save allocator state
  const uint64_t num_rows = (key_length+1);
  const uint64_t num_rows_1 = num_rows-1;
  const uint64_t num_columns = (text_length+1);
  PROF_ADD_COUNTER(GP_SWG_ALIGN_BANDED_CELLS,num_columns*max_bandwidth);
  swg_cell_t** const dp = align_swg_allocate_table(num_columns+1,num_rows+1,mm_allocator);
  // Initialize DP-matrix
  const swg_matching_score_t* const matching_score = &swg_penalties->matching_score;
  const int32_t gap_extension = swg_penalties->gap_extension_score;
  const int32_t gap_open = swg_penalties->gap_open_score;
  const int32_t single_gap = gap_open + gap_extension; // g(1)
  // Initialize band
  int64_t column_start_band = max_bandwidth + 2;
  if (begin_free || column_start_band > text_length) column_start_band = text_length + 2;
  uint64_t band_high_offset = 0;
  uint64_t band_low_offset = (max_bandwidth+1 < num_rows) ? max_bandwidth+2 : num_rows;
  // Initialize first/second column
  align_swg_init_table_banded_opt(
      dp,num_columns,num_rows,column_start_band,
      band_low_offset,begin_free,single_gap,gap_extension);
  /*
   * Compute DP-matrix
   */
  // Init as full deletion
  int32_t max_score = gap_open + key_length*gap_extension;
  uint64_t column, max_score_column = 0;
  for (column=1;column<num_columns;++column) {
    // Compute column
    align_swg_banded_compute_column(
        key,text,dp,column,num_rows,column_start_band,
        &band_high_offset,&band_low_offset,matching_score,
        single_gap,gap_extension);
    // Check score
    if (end_free && band_low_offset==num_rows) {
      if (dp[column][num_rows_1].M > max_score) {
        max_score = dp[column][num_rows_1].M;
        max_score_column = column;
      }
    }
  }
  // Set alignment score/column
  if (!end_free) {
    max_score = dp[num_columns-1][num_rows_1].M;
    max_score_column = num_columns-1;
  }
  // Retrieve the alignment. Store the match (Backtrace and generate CIGAR)
  align_swg_traceback(
      align_input,dp,max_score,max_score_column,
      single_gap,begin_free,match_alignment,cigar_vector);
  // Clean-up
  mm_allocator_pop_state(mm_allocator); // Free
  PROF_STOP(GP_SWG_ALIGN_BANDED);
}
