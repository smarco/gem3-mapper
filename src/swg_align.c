/*
 * PROJECT: GEMMapper
 * FILE: swg_align.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "swg_align.h"
#include "matches.h"

/*
 * Levenshtein Distance  [Dynamic Programming]
 */
#define GT_DP(i,j) dp_array[(i)*pattern_len+(j)]
GEM_INLINE void align_levenshtein_dp_matrix_print(
    uint64_t* const dp_array,const uint64_t pattern_len,const uint64_t sequence_len,
    const uint64_t pattern_limit,const uint64_t sequence_limit) {
  uint64_t i, j;
  for (j=0;j<pattern_limit;++j) {
    for (i=0;i<sequence_limit;++i) {
      fprintf(stdout,"%02"PRIu64" ",GT_DP(i,j));
    }
    fprintf(stdout,"\n");
  }
  fprintf(stdout,"\n");
}
GEM_INLINE int64_t align_levenshtein_get_distance(
    const char* const pattern,const uint64_t pattern_length,
    const char* const sequence,const uint64_t sequence_length,
    const bool ends_free,uint64_t* const position) {
  GEM_CHECK_NULL(pattern); GEM_CHECK_ZERO(pattern_length);
  GEM_CHECK_NULL(sequence); GEM_CHECK_ZERO(sequence_length);
  // Allocate DP-matrix
  const uint64_t pattern_len = pattern_length+1;
  const uint64_t sequence_len = sequence_length+1;
  uint64_t* dp_array[2];
  dp_array[0] = mm_calloc(2*pattern_len,uint64_t,false);
  dp_array[1] = dp_array[0] + pattern_len;
  // Init DP-Matrix
  uint64_t min_val = UINT64_MAX, i_pos = UINT64_MAX;
  uint64_t i, j, idx_a=0, idx_b=0;
  for (j=0;j<pattern_len;++j) dp_array[0][j]=j;
  // Calculate DP-Matrix
  for (i=1;i<sequence_len;++i) {
    // Fix indexes
    idx_a = idx_b;
    idx_b = i % 2;
    // Fix first cell
    dp_array[idx_b][0] = (ends_free) ? 0 : dp_array[idx_a][0]+1;
    // Develop row
    for (j=1;j<pattern_len;++j) {
      const uint64_t ins = dp_array[idx_a][j]   + 1;
      const uint64_t del = dp_array[idx_b][j-1] + 1;
      const uint64_t sub = dp_array[idx_a][j-1] + ((sequence[i-1]==pattern[j-1]) ? 0 : 1);
      dp_array[idx_b][j] = MIN(sub,MIN(ins,del));
    }
    // Check last cell value
    if (ends_free && dp_array[idx_b][pattern_length] < min_val) {
      min_val = dp_array[idx_b][pattern_length];
      i_pos = i;
    }
  }
  // Return results & Free
  int64_t distance = INT64_MAX;
  if (ends_free) {
    *position = i_pos-1;
    distance = min_val;
  } else {
    *position = pattern_length;
    distance = dp_array[idx_b][pattern_length];
  }
  mm_free(dp_array[0]);
  return distance;
}
GEM_INLINE void swg_align_dp_matrix_print(swg_cell_t** const dp,const uint64_t num_rows,const uint64_t num_columns) {
  uint64_t i, j;
  for (j=0;j<num_rows;++j) {
    for (i=0;i<num_columns;++i) {
      // fprintf(stdout,"%02s/%02s/%02s ",dp[i][j].M,dp[i][j].I,dp[i][j].D);
      fprintf(stdout,"%+04d ",dp[i][j].M);
    }
    fprintf(stdout,"\n");
  }
  fprintf(stdout,"\n");
}
GEM_INLINE void swg_align_match(
    const uint8_t* const key,const uint64_t key_length,const swg_penalties_t* swg_penalties,
    uint64_t* const match_position,const uint8_t* const sequence,const uint64_t sequence_length,
    const uint64_t max_distance,vector_t* const cigar_vector,uint64_t* const cigar_vector_offset,
    uint64_t* const cigar_length,uint64_t* const alignment_score,uint64_t* const effective_length,
    mm_stack_t* const mm_stack) {
  /*
   * Initialize
   */
  uint64_t column, row;
  mm_stack_push_state(mm_stack); // Save stack state
  const uint64_t row_size = (key_length+1)*sizeof(swg_cell_t);
  const uint64_t num_rows = (key_length+1);
  const uint64_t num_columns = (sequence_length+1);
  swg_cell_t** const dp = mm_stack_malloc(mm_stack,num_columns*sizeof(swg_cell_t*));
  swg_cell_t* dp_matrix = mm_stack_malloc(mm_stack,num_columns*2*row_size);
  for (column=0;column<num_columns;++column) {
    dp[column] = dp_matrix;
    dp_matrix += row_size;
  }
  // Initialize first row
  const swg_int single_gap = - (swg_penalties->gap_open_penalty + swg_penalties->gap_extension_penalty); // g(1)
  const swg_int gap_extension = - (swg_penalties->gap_extension_penalty);
  for (column=0;column<num_columns;++column) {
    dp[column][0].D = SWG_INT_MIN;
    dp[column][0].M = 0;
  }
  dp[0][0].I = SWG_INT_MIN;
  dp[0][1].I = SWG_INT_MIN;
  dp[0][1].M = single_gap; // g(1)
  for (row=2;row<num_rows;++row) {
    dp[0][row].I = SWG_INT_MIN;
    dp[0][row].M = dp[0][row-1].M + gap_extension; // g(row)
  }
  /*
   * Compute DP-matrix
   */
  swg_int max_score = SWG_INT_MIN;
  uint32_t max_score_column = SWG_INT_MIN;
  const swg_int match = swg_penalties->matching_score;
  const swg_int mismatch = - (swg_penalties->mismatch_penalty);
  for (column=1;column<num_columns;++column) {
    for (row=1;row<num_rows;++row) {
      // Update DP.D
      const swg_int del_new = dp[column][row-1].M + single_gap;
      const swg_int del_ext = dp[column][row-1].D + gap_extension;
      const swg_int del = MAX(del_new,del_ext);
      dp[column][row].D = del;
      // Update DP.I
      const swg_int ins_new = dp[column-1][row].M + single_gap;
      const swg_int ins_ext = dp[column-1][row].I + gap_extension;
      const swg_int ins = MAX(ins_new,ins_ext);
      dp[column][row].I = ins;
      // Update DP.M
      const swg_int m_match = dp[column-1][row-1].M + ((key[row-1]==sequence[column-1]) ? match : mismatch);
      dp[column][row].M = MAX(m_match,MAX(ins,del));
    }
    // Check score
    if (dp[column][num_rows-1].M > max_score) {
      max_score = dp[column][num_rows-1].M;
      max_score_column = column;
    }
  }
  //swg_align_dp_matrix_print(dp,num_rows,20);
  /*
   * Retrieve the alignment. Store the match (Backtrace and generate CIGAR)
   */
  // Set alignment-score
  *alignment_score = max_distance; // *alignment_score = max_score;
  // Allocate CIGAR string memory (worst case)
  *cigar_vector_offset = vector_get_used(cigar_vector); // Set CIGAR offset
  vector_reserve_additional(cigar_vector,key_length); // Reserve
  cigar_element_t* cigar_buffer = vector_get_free_elm(cigar_vector,cigar_element_t); // Sentinel
  cigar_element_t* const cigar_buffer_base = cigar_buffer;
  cigar_buffer->type = cigar_null; // Trick
  // Start Backtrace
  uint64_t match_effective_length = key_length;
  uint64_t h = max_score_column;
  uint64_t v = key_length;
  cigar_t traceback_matrix = cigar_match;
  while (v > 0 && h > 0) {
    switch (traceback_matrix) {
      case cigar_del:
        // Traceback D-matrix
        ALIGN_ADD_CIGAR_ELEMENT(cigar_buffer,cigar_del,1); // Deletion <-1>@v
        if (dp[h][v].D == dp[h][v-1].M + single_gap) traceback_matrix = cigar_match;
        --v; --match_effective_length;
        break;
      case cigar_ins:
        // Traceback I-matrix
        ALIGN_ADD_CIGAR_ELEMENT(cigar_buffer,cigar_ins,1); // Insertion <+1>@v
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
          if (sequence[h-1] != key[v-1]) {
            // Mismatch
            if (cigar_buffer->type!=cigar_null) ++(cigar_buffer);
            cigar_buffer->type = cigar_mismatch;
            cigar_buffer->mismatch = sequence[h-1];
            --h; --v;
          } else {
            ALIGN_ADD_CIGAR_ELEMENT(cigar_buffer,cigar_match,1); // Match
            --h; --v;
          }
        }
        break;
    }
  }
  if (v > 0) {
    ALIGN_ADD_CIGAR_ELEMENT(cigar_buffer,cigar_del,v); // <-(@v+1)>@v
    match_effective_length -= v;
  }
  if (h > 0) {
    *match_position += h; // We need to correct the matching_position
  }
  // Set effective length
  *effective_length = match_effective_length;
  // Set CIGAR buffer used
  if (cigar_buffer->type!=cigar_null) ++(cigar_buffer);
  const uint64_t num_cigar_elements = cigar_buffer - cigar_buffer_base;
  vector_add_used(cigar_vector,num_cigar_elements);
  *cigar_length = num_cigar_elements; // Set CIGAR length
  // Reverse CIGAR Elements
  if (num_cigar_elements > 0) {
    const uint64_t middle_point = num_cigar_elements/2;
    uint64_t i;
    for (i=0;i<middle_point;++i) {
      SWAP(cigar_buffer_base[i],cigar_buffer_base[num_cigar_elements-i-1]);
    }
  }
  // Clean-up
  mm_stack_pop_state(mm_stack,false); // Free
}
