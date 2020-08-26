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
 * SWG Compute Column Banded
 */
void align_swg_banded_compute_column(
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length,
    const bool reverse_strings,
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
    const uint8_t enc_text = text[(reverse_strings) ? text_length - column : column-1];
    const uint8_t enc_key = key[(reverse_strings) ? key_length - row : row-1];
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
void align_swg_banded_compute_column__max_score(
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length,
    const bool reverse_strings,
    swg_cell_t** const dp,
    const uint64_t column,
    const uint64_t num_rows,
    const uint64_t column_start_band,
    uint64_t* const band_high_offset,
    uint64_t* const band_low_offset,
    const swg_matching_score_t* const matching_score,
    const uint64_t single_gap,
    const uint64_t gap_extension,
    int32_t* const max_score,
    uint64_t* const max_score_row) {
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
  *max_score = SWG_SCORE_MIN;
  for (row=*band_high_offset+1;row<*band_low_offset;++row) {
    // Update DP.M
    const uint8_t enc_text = text[(reverse_strings) ? text_length - column : column-1];
    const uint8_t enc_key = key[(reverse_strings) ? key_length - row : row-1];
    const int32_t d_value = dp[column][row].D;
    const int32_t i_value = dp[column][row].I;
    const int32_t match = dp[column-1][row-1].M + (*matching_score)[enc_text][enc_key];
    const int32_t m_value = MAX(match,MAX(d_value,i_value));
    dp[column][row].M = m_value;
    // Update Max Score
    if (*max_score < m_value) {
      *max_score = m_value;
      *max_score_row = row;
    }
    // Update DP.D
    const int32_t gap_new = m_value + single_gap;
    const int32_t del_ext = d_value + gap_extension;
    dp[column][row+1].D = MAX(gap_new,del_ext);
    // Update DP.I
    const int32_t ins_ext = i_value + gap_extension;
    dp[column+1][row].I = MAX(gap_new,ins_ext);
  }
}
/*
 * SWG Banded
 */
void align_swg_banded(
    const uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    uint64_t text_length,
    const swg_penalties_t* const swg_penalties,
    const uint64_t max_bandwidth,
    const bool begin_free,
    const bool end_free,
    const bool left_gap_alignment,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_allocator_t* const mm_allocator) {
  // Check band-limits
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
  mm_allocator_push_state(mm_allocator); // Save allocator state
  // Parameters
  const swg_matching_score_t* const matching_score = &swg_penalties->matching_score;
  const int32_t gap_extension = swg_penalties->gap_extension_score;
  const int32_t gap_open = swg_penalties->gap_open_score;
  const int32_t single_gap = gap_open + gap_extension; // g(1)
  const uint64_t num_rows = (key_length+1);
  const uint64_t num_rows_1 = num_rows-1;
  const uint64_t num_columns = (text_length+1);
  PROF_ADD_COUNTER(GP_SWG_ALIGN_BANDED_CELLS,num_columns*max_bandwidth);
  // Initialize band
  int64_t column_start_band = max_bandwidth + 2;
  if (begin_free || column_start_band > text_length) column_start_band = text_length + 2;
  uint64_t band_high_offset = 0;
  uint64_t band_low_offset = (max_bandwidth+1 < num_rows) ? max_bandwidth+2 : num_rows;
  // Initialize DP-matrix
  swg_cell_t** const dp = align_swg_allocate_table(num_columns+1,num_rows+1,mm_allocator);
  align_swg_init_table_banded(
      dp,num_columns,num_rows,column_start_band,
      band_low_offset,begin_free,single_gap,gap_extension); // Initialize first/second column
  // Compute DP-matrix
  int32_t max_score = gap_open + key_length*gap_extension; // Full deletion
  uint64_t column, max_score_column = 0;
  for (column=1;column<num_columns;++column) {
    // Compute column
    align_swg_banded_compute_column(
        key,key_length,text,text_length,false,
        dp,column,num_rows,column_start_band,
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
      key,key_length,text,text_length,false,
      dp,max_score,key_length,max_score_column,
      single_gap,matching_score,
      begin_free,left_gap_alignment,
      match_alignment,cigar_vector);
  // Clean-up
  mm_allocator_pop_state(mm_allocator); // Free
  PROF_STOP(GP_SWG_ALIGN_BANDED);
}
void align_swg_banded_extend(
    const uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    uint64_t text_length,
    const bool reverse_extension,
    const bool local_extension,
    const int32_t base_score,
    const swg_penalties_t* const swg_penalties,
    const uint64_t max_bandwidth,
    const bool left_gap_alignment,
    match_alignment_t* const match_alignment,
    uint64_t* const max_key_aligned,
    uint64_t* const max_text_aligned,
    vector_t* const cigar_vector,
    mm_allocator_t* const mm_allocator) {
  // Check trivial cases
  if (key_length == 0) {
    match_alignment->score = 0;
    *max_key_aligned = 0;
    *max_text_aligned = 0;
    return;
  }
  if (text_length == 0) {
    matches_cigar_vector_append_deletion(
        cigar_vector,&match_alignment->cigar_length,
        key_length,cigar_attr_none);
    match_alignment->score =
        swg_penalties->gap_open_score +
        swg_penalties->gap_extension_score * key_length;
    *max_key_aligned = key_length;
    *max_text_aligned = 0;
    return;
  }
  // Parameters
  const swg_matching_score_t* const matching_score = &swg_penalties->matching_score;
  const int32_t gap_extension = swg_penalties->gap_extension_score;
  const int32_t gap_open = swg_penalties->gap_open_score;
  const int32_t single_gap = gap_open + gap_extension; // g(1)
  const int32_t generic_match_score = swg_penalties->generic_match_score;
  const uint64_t num_rows = (key_length+1);
  const uint64_t num_columns = (text_length+1);
  PROF_START(GP_SWG_ALIGN_BANDED);
  PROF_ADD_COUNTER(GP_SWG_ALIGN_BANDED_CELLS,num_columns*max_bandwidth);
  // Initialize band
  int64_t column_start_band = max_bandwidth + 2;
  if (column_start_band > text_length) column_start_band = text_length + 2;
  uint64_t band_high_offset = 0;
  uint64_t band_low_offset = (max_bandwidth+1 < num_rows) ? max_bandwidth+2 : num_rows;
  // Initialize DP-matrix (lazy allocation)
  mm_allocator_push_state(mm_allocator); // Save allocator state
  swg_cell_t** const dp = mm_allocator_malloc(mm_allocator,(num_columns+1)*sizeof(swg_cell_t*));
  align_swg_allocate__init_column_banded(
      dp,0,column_start_band,
      band_low_offset,band_high_offset,
      single_gap,gap_extension,mm_allocator);
  // Compute DP-matrix
  int32_t max_score = 0, max_local_score = 0;
  uint64_t max_score_column = 0, max_score_row = 0;
  uint64_t column, max_local_score_row = 0;
  if (!local_extension) {
    max_score = gap_open + key_length*gap_extension; // Init as full deletion
    max_score_column = 0;
    max_score_row = key_length;
  }
  for (column=1;column<num_columns;++column) {
    // Initialize column
    align_swg_allocate__init_column_banded(
        dp,column,column_start_band,
        band_low_offset,band_high_offset,
        single_gap,gap_extension,mm_allocator);
    // Compute column
    align_swg_banded_compute_column__max_score(
        key,key_length,text,text_length,reverse_extension,
        dp,column,num_rows,column_start_band,
        &band_high_offset,&band_low_offset,
        matching_score,single_gap,gap_extension,
        &max_local_score,&max_local_score_row);
    // Check score
    if (local_extension) {
      if (base_score + max_local_score <= 0) {
        ++column;
        break; // Drop-off
      }
    } else if (band_low_offset==num_rows) {
      if (dp[column][band_low_offset-1].M > max_score) {
        max_score_column = column;
        max_score_row = num_rows-1;
        max_score = dp[column][num_rows-1].M;
      } else {
        const int64_t key_left = num_rows - max_local_score_row;
        const int64_t max_expected_score = max_local_score + key_left*generic_match_score;
        if (max_expected_score <= max_score) break; // Drop-off
      }
    }
  }
  if (local_extension) {
    max_score_column = column-1;
    max_score_row = band_low_offset-1;
    max_score = dp[max_score_column][max_score_row].M;
  }
  *max_key_aligned = max_score_row;
  *max_text_aligned = max_score_column;
  // Retrieve the alignment. Store the match (Backtrace and generate CIGAR)
  align_swg_traceback(
      key,key_length,text,text_length,reverse_extension,
      dp,max_score,max_score_row,max_score_column,
      single_gap,matching_score,false,
      left_gap_alignment!=reverse_extension, // Account for reversed input
      match_alignment,cigar_vector);
  // Clean-up
  mm_allocator_pop_state(mm_allocator); // Free
  PROF_STOP(GP_SWG_ALIGN_BANDED);
}




