/*
 * PROJECT: GEMMapper
 * FILE: align_swg.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */


#include "align/align_swg.h"
#include "matches/matches_cigar.h"

/*
 * SWG - BackTrace
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text
 *   @align_input->text_length
 *   @match_alignment->match_position (Correction)
 *   @match_alignment->cigar_length (Cumulative)
 */
void align_swg_traceback(
    match_align_input_t* const align_input,
    swg_cell_t** const dp,
    const int32_t max_score,
    const uint64_t max_score_column,
    const int32_t single_gap,
    const int32_t gap_extension,
    const bool begin_free,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector) {
  // Parameters
  const uint8_t* const key = align_input->key;
  const uint64_t key_length = align_input->key_length;
  uint8_t* const text = align_input->text;
  const uint64_t text_length = align_input->text_length;
  // Allocate CIGAR string memory (worst case)
  vector_reserve_additional(cigar_vector,key_length+text_length); // Reserve
  cigar_element_t* cigar_buffer_sentinel = vector_get_free_elm(cigar_vector,cigar_element_t); // Sentinel
  cigar_element_t* const cigar_buffer_base = cigar_buffer_sentinel;
  cigar_buffer_sentinel->type = cigar_null; // Trick
  // Start Backtrace
  int64_t match_effective_length = key_length;
  int32_t match_alignment_score = max_score;
  int64_t h = max_score_column;
  int64_t v = key_length;
  cigar_t traceback_matrix = cigar_match;
  while (v > 0 && h > 0) {
    switch (traceback_matrix) {
      case cigar_del:
        // Traceback D-matrix
        matches_cigar_buffer_add_cigar_element(&cigar_buffer_sentinel,cigar_del,1); // Deletion <-1>@v
        if (dp[h][v].D == dp[h][v-1].M + single_gap) traceback_matrix = cigar_match;
        --v; --match_effective_length;
        break;
      case cigar_ins:
        // Traceback I-matrix
        matches_cigar_buffer_add_cigar_element(&cigar_buffer_sentinel,cigar_ins,1); // Insertion <+1>@v
        if (dp[h][v].I == dp[h-1][v].M + single_gap) traceback_matrix = cigar_match;
        --h; ++match_effective_length;
        break;
      default:
        // Traceback M-matrix
        if (dp[h][v].M == dp[h][v].D) {
          traceback_matrix = cigar_del;
        } else if (dp[h][v].M == dp[h][v].I) {
          traceback_matrix = cigar_ins;
        } else {
          if (key[v-1] != text[h-1]) {
            // Mismatch
            matches_cigar_buffer_add_mismatch(&cigar_buffer_sentinel,text[h-1]);
            --h; --v;
          } else {
            matches_cigar_buffer_add_cigar_element(&cigar_buffer_sentinel,cigar_match,1); // Match
            --h; --v;
          }
        }
        break;
    }
  }
  if (v > 0) {
    matches_cigar_buffer_add_cigar_element(&cigar_buffer_sentinel,cigar_del,v); // <-(@v+1)>@v
    match_effective_length -= v;
  }
  if (h > 0) {
    if (begin_free) {
      match_alignment->match_position += h; // We need to correct the matching_position
    } else {
      matches_cigar_buffer_add_cigar_element(&cigar_buffer_sentinel,cigar_ins,h); // <-(@h+1)>@h
      match_effective_length += h;
    }
  }
  // Set CIGAR buffer used
  if (cigar_buffer_sentinel->type!=cigar_null) ++(cigar_buffer_sentinel);
  const uint64_t num_cigar_elements = cigar_buffer_sentinel - cigar_buffer_base;
  vector_add_used(cigar_vector,num_cigar_elements);
  // Reverse CIGAR Elements
  if (num_cigar_elements > 0) {
    const uint64_t middle_point = num_cigar_elements/2;
    uint64_t i;
    for (i=0;i<middle_point;++i) {
      SWAP(cigar_buffer_base[i],cigar_buffer_base[num_cigar_elements-i-1]);
    }
  }
  // Update CIGAR/Score values
  match_alignment->cigar_length += num_cigar_elements; // Update CIGAR length
  match_alignment->effective_length = match_effective_length; // Update effective length
  match_alignment->score = match_alignment_score; // Update alignment-score
}
/*
 * SWG - Init
 */
swg_cell_t** align_swg_allocate_table(
    const uint64_t num_columns,
    const uint64_t num_rows,
    mm_stack_t* const mm_stack) {
  swg_cell_t** const dp = mm_stack_malloc(mm_stack,num_columns*sizeof(swg_cell_t*));
  const uint64_t row_size = num_rows*sizeof(swg_cell_t);
  uint64_t column;
  for (column=0;column<num_columns;++column) {
    dp[column] = mm_stack_malloc(mm_stack,row_size);
  }
  return dp;
}
void align_swg_init_table(
    swg_cell_t** const dp,
    const uint64_t num_columns,
    const uint64_t num_rows,
    const int32_t single_gap,
    const int32_t gap_extension) {
  uint64_t column, row;
  for (column=0;column<num_columns;++column) {
    dp[column][0].D = SWG_SCORE_MIN;
    dp[column][0].M = 0;
  }
  dp[0][0].I = SWG_SCORE_MIN;
  dp[0][1].I = SWG_SCORE_MIN;
  dp[0][1].M = single_gap; // g(1)
  for (row=2;row<num_rows;++row) {
    dp[0][row].I = SWG_SCORE_MIN;
    dp[0][row].M = dp[0][row-1].M + gap_extension; // g(row)
  }
}
void align_swg_init_table_banded(
    swg_cell_t** const dp,
    const uint64_t num_columns,
    const uint64_t num_rows,
    const uint64_t column_start_band,
    const uint64_t band_low_offset,
    const bool begin_free,
    const int32_t single_gap,
    const int32_t gap_extension) {
  uint64_t column, row;
  // Initialize first column
  dp[0][0].D = SWG_SCORE_MIN; // Not used
  dp[0][0].I = SWG_SCORE_MIN; // Not used
  dp[0][0].M = 0;
  dp[0][1].M = single_gap; // g(1)
  dp[0][1].I = SWG_SCORE_MIN;
  for (row=2;row<band_low_offset;++row) {
    dp[0][row].M = dp[0][row-1].M + gap_extension; // g(row)
    dp[0][row].I = SWG_SCORE_MIN;
  }
  // Initialize first row
  if (begin_free) {
    for (column=1;column<num_columns;++column) {
      dp[column][0].M = 0;
      dp[column][0].D = SWG_SCORE_MIN;
    }
  } else {
    dp[1][0].M = single_gap;
    dp[1][0].D = SWG_SCORE_MIN;
    for (column=2;column<column_start_band;++column) {
      dp[column][0].M = dp[column-1][0].M + gap_extension;
      dp[column][0].D = SWG_SCORE_MIN;
    }
  }
}
void align_swg_init_table_banded_opt(
    swg_cell_t** const dp,
    const uint64_t num_columns,
    const uint64_t num_rows,
    const uint64_t column_start_band,
    const uint64_t band_low_offset,
    const bool begin_free,
    const int32_t single_gap,
    const int32_t gap_extension) {
  dp[0][0].D = SWG_SCORE_MIN; // Not used
  dp[0][0].I = SWG_SCORE_MIN; // Not used
  dp[0][0].M = 0;
  dp[0][1].M = single_gap; // g(1)
  dp[0][1].I = SWG_SCORE_MIN;
  dp[1][1].I = SWG_SCORE_MIN;
  uint64_t column, row;
  // Initialize first/second column
  for (row=2;row<band_low_offset;++row) {
    dp[0][row].M = dp[0][row-1].M + gap_extension; // g(row)
    dp[0][row].I = SWG_SCORE_MIN;
    dp[1][row].I = dp[0][row].M + single_gap;
  }
  // Initialize first/second row
  if (begin_free) {
    for (column=1;column<num_columns;++column) {
      dp[column][0].M = 0;
      dp[column][0].D = SWG_SCORE_MIN;
      dp[column][1].D = dp[column][0].M + single_gap;
    }
  } else {
    dp[1][0].M = single_gap;
    dp[1][0].D = SWG_SCORE_MIN;
    dp[1][1].D = single_gap + single_gap;
    for (column=2;column<column_start_band;++column) {
      dp[column][0].M = dp[column-1][0].M + gap_extension;
      dp[column][0].D = SWG_SCORE_MIN;
      dp[column][1].D = dp[column][0].M + single_gap;
    }
  }
}
/*
 * Smith-waterman-gotoh Base (ie. no-optimizations)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text
 *   @align_input->text_length
 *   @align_parameters->swg_penalties
 *   @match_alignment->match_position (Adjusted)
 */
void align_swg_base(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_stack_t* const mm_stack) {
  // Parameters
  const uint8_t* const key = align_input->key;
  const uint64_t key_length = align_input->key_length;
  uint8_t* const text = align_input->text;
  const uint64_t text_length = align_input->text_length;
  const swg_penalties_t* swg_penalties = align_parameters->swg_penalties;
  // Initialize
  mm_stack_push_state(mm_stack); // Save stack state
  const uint64_t num_rows = (key_length+1);
  const uint64_t num_columns = (text_length+1);
  uint64_t column, row;
  swg_cell_t** const dp = align_swg_allocate_table(num_columns,num_rows,mm_stack);
  // Initialize first row
  const swg_matching_score_t* const matching_score = &swg_penalties->matching_score;
  const int32_t single_gap = swg_penalties->gap_open_score + swg_penalties->gap_extension_score; // g(1)
  const int32_t gap_extension = swg_penalties->gap_extension_score;
  align_swg_init_table(dp,num_columns,num_rows,single_gap,gap_extension);
  // Compute DP-matrix
  int32_t max_score = SWG_SCORE_MIN;
  uint64_t max_score_column = UINT64_MAX;
  for (column=1;column<num_columns;++column) {
    for (row=1;row<num_rows;++row) {
      // Update DP.D
      const int32_t del_new = dp[column][row-1].M + single_gap;
      const int32_t del_ext = dp[column][row-1].D + gap_extension;
      const int32_t del = MAX(del_new,del_ext);
      dp[column][row].D = del;
      // Update DP.I
      const int32_t ins_new = dp[column-1][row].M + single_gap;
      const int32_t ins_ext = dp[column-1][row].I + gap_extension;
      const int32_t ins = MAX(ins_new,ins_ext);
      dp[column][row].I = ins;
      // Update DP.M
      const uint8_t enc_text = text[column-1];
      const uint8_t enc_key = key[row-1];
      const int32_t m_match = dp[column-1][row-1].M + (*matching_score)[enc_text][enc_key];
      dp[column][row].M = MAX(m_match,MAX(ins,del));
    }
    // Check score
    if (dp[column][num_rows-1].M > max_score) {
      max_score = dp[column][num_rows-1].M;
      max_score_column = column;
    }
  }
  // DEBUG
    // align_swg_match_table_print(dp,MIN(10,key_length),MIN(10,key_length));
  // Init Match
  match_alignment->cigar_length = 0;
  match_alignment->effective_length = 0;
  match_alignment->score = 0;
  // Retrieve the alignment. Store the match (Backtrace and generate CIGAR)
  align_swg_traceback(align_input,dp,max_score,max_score_column,
      single_gap,gap_extension,true,match_alignment,cigar_vector);
  // Clean-up
  mm_stack_pop_state(mm_stack); // Free
}
/*
 * SWG Full (Computes full DP-matrix)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text
 *   @align_input->text_length
 *   @align_parameters->swg_penalties
 *   @match_alignment->match_position (Correction)
 *   @match_alignment->cigar_length (Cumulative)
 */
void align_swg_full(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    const bool begin_free,
    const bool end_free,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_stack_t* const mm_stack) {
  PROF_START(GP_SWG_ALIGN_FULL);
  // Parameters
  const uint8_t* const key = align_input->key;
  const uint64_t key_length = align_input->key_length;
  uint8_t* const text = align_input->text;
  const uint64_t text_length = align_input->text_length;
  const swg_penalties_t* swg_penalties = align_parameters->swg_penalties;
  // Allocate memory
  mm_stack_push_state(mm_stack); // Save stack state
  const uint64_t num_rows = (key_length+1);
  const uint64_t num_rows_1 = num_rows-1;
  const uint64_t num_columns = (text_length+1);
  swg_cell_t** const dp = align_swg_allocate_table(num_columns+1,num_rows+1,mm_stack);
  // Initialize DP-matrix
  const swg_matching_score_t* const matching_score = &swg_penalties->matching_score;
  const int32_t gap_extension = swg_penalties->gap_extension_score;
  const int32_t gap_open = swg_penalties->gap_open_score;
  const int32_t single_gap = gap_open + gap_extension; // g(1)
  align_swg_init_table_banded_opt(dp,num_columns,num_rows,
      num_columns,num_rows,false,single_gap,gap_extension); // Initialize first/second column
  /*
   * Compute DP-matrix
   */
  // Init as full deletion
  int32_t max_score = gap_open + key_length*gap_extension;
  uint64_t column, row, max_score_column = 0;
  for (column=1;column<num_columns;++column) {
    for (row=1;row<num_rows;++row) {
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
    // Check score
    if (end_free && row==num_rows) {
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
  align_swg_traceback(align_input,dp,max_score,max_score_column,
      single_gap,gap_extension,begin_free,match_alignment,cigar_vector);
  // Clean-up
  mm_stack_pop_state(mm_stack); // Free
  PROF_STOP(GP_SWG_ALIGN_FULL);
}
/*
 * SWG Banded
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text
 *   @align_input->text_length
 *   @align_parameters->swg_penalties
 *   @align_parameters->max_bandwidth
 *   @match_alignment->match_position (Correction)
 *   @match_alignment->cigar_length (Cumulative)
 */
void align_swg_banded(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    const bool begin_free,
    const bool end_free,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_stack_t* const mm_stack) {
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
  mm_stack_push_state(mm_stack); // Save stack state
  const uint64_t num_rows = (key_length+1);
  const uint64_t num_rows_1 = num_rows-1;
  const uint64_t num_columns = (text_length+1);
  PROF_ADD_COUNTER(GP_SWG_ALIGN_BANDED_CELLS,num_columns*max_bandwidth);
  swg_cell_t** const dp = align_swg_allocate_table(num_columns+1,num_rows+1,mm_stack);
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
  uint64_t column, row, max_score_column = 0;
  for (column=1;column<num_columns;++column) {
    // Initialize band boundaries & update band limits
    if (band_low_offset < num_rows) { //  Below band
      dp[column][band_low_offset].I = SWG_SCORE_MIN;
      dp[column-1][band_low_offset].I = SWG_SCORE_MIN;
      ++band_low_offset; // Swift band
    }
    // Initialize band boundaries & update band limits
    if (column >= column_start_band) {
      ++band_high_offset; // Swift band
      dp[column][band_high_offset+1].D = SWG_SCORE_MIN;
    }
    // Locate the cursor at the proper cell & calculate DP
    for (row=band_high_offset+1;row<band_low_offset;++row) {
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
  align_swg_traceback(align_input,dp,max_score,max_score_column,
      single_gap,gap_extension,begin_free,match_alignment,cigar_vector);
  // Clean-up
  mm_stack_pop_state(mm_stack); // Free
  PROF_STOP(GP_SWG_ALIGN_BANDED);
}
/*
 * Smith-Waterman-Gotoh - Main procedure (Dispatcher)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text
 *   @align_input->text_length
 *   @align_parameters->swg_penalties
 *   @align_parameters->allowed_enc
 *   @align_parameters->max_bandwidth
 *   @match_alignment->match_position (Adjusted)
 *   @match_alignment->cigar_length (Cumulative)
 */
void align_swg(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    const bool begin_free,
    const bool end_free,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_stack_t* const mm_stack) {
  // Parameters
  const uint8_t* const key = align_input->key;
  const uint64_t key_length = align_input->key_length;
  uint8_t* const text = align_input->text;
  const uint64_t text_length = align_input->text_length;
  const swg_penalties_t* swg_penalties = align_parameters->swg_penalties;
  const bool* const allowed_enc = align_parameters->allowed_enc;
  // Check lengths
  if (key_length == 0 && text_length == 0) {
    match_alignment->score = 0;
    match_alignment->effective_length = 0;
  } else if (key_length == 0 && text_length > 0) {
    if (begin_free) {
      // Adjust position
      match_alignment->match_position += text_length;
      match_alignment->score = 0;
      match_alignment->effective_length = 0;
    } else if (!end_free) {
      // Insertion <+@text_length>
      matches_cigar_vector_append_insertion(cigar_vector,&match_alignment->cigar_length,text_length,cigar_attr_none);
      match_alignment->score = align_swg_score_insertion(swg_penalties,text_length);
      match_alignment->effective_length = text_length;
    } else {
      match_alignment->score = 0;
      match_alignment->effective_length = 0;
    }
  } else if (key_length > 0 && text_length == 0) {
    // Deletion <-@key_length>
    matches_cigar_vector_append_deletion(cigar_vector,&match_alignment->cigar_length,key_length,cigar_attr_none);
    match_alignment->score = align_swg_score_deletion(swg_penalties,key_length);
  } else if (key_length==1 && text_length==1) {
    // Mismatch/Match
    const uint8_t key_enc = key[0];
    const uint8_t text_enc = text[0];
    if (!allowed_enc[text_enc] || text_enc != key_enc) {
      matches_cigar_vector_append_mismatch(cigar_vector,&match_alignment->cigar_length,text_enc,cigar_attr_none);
      match_alignment->score = align_swg_score_mismatch(swg_penalties);
    } else {
      matches_cigar_vector_append_match(cigar_vector,&match_alignment->cigar_length,1,cigar_attr_none);
      match_alignment->score = align_swg_score_match(swg_penalties,1);
    }
    match_alignment->effective_length = 1;
  } else {
    PROF_ADD_COUNTER(GP_SWG_ALIGN_BANDED_LENGTH,text_length);
    align_swg_banded(align_input,align_parameters,
        begin_free,end_free,match_alignment,cigar_vector,mm_stack);
  }
}
/*
 * Display
 */
void align_swg_print_table(
    swg_cell_t** const dp,
    const uint64_t num_columns,
    const uint64_t num_rows) {
  uint64_t i, j;
  for (i=0;i<=num_rows;++i) {
    for (j=0;j<=num_columns;++j) {
      printf("%+4d ",dp[j][i].M);
    }
    printf("\n");
  }
}
void align_swg_print_input(
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length) {
  uint64_t i;
  printf("KEY> ");
  for (i=0;i<key_length;++i) printf("%c",dna_decode(key[i]));
  printf("\n");
  printf("TXT> ");
  for (i=0;i<text_length;++i) printf("%c",dna_decode(text[i]));
  printf("\n");
}
