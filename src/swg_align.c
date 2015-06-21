/*
 * PROJECT: GEMMapper
 * FILE: swg_align.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */


#include "swg_align.h"
#include "bpm_align.h"
#include "match_elements.h"
#include "match_scaffold.h"
#include "match_align_dto.h"
#include "matches.h"

/*
 * Constants
 */
#define SWG_SCORE_INT32_MIN (INT16_MIN)

/*
 * Check CIGAR string
 */
GEM_INLINE bool align_check_match(
    FILE* const stream,const uint8_t* const key,const uint64_t key_length,const uint8_t* const text,
    const uint64_t text_length,vector_t* const cigar_vector,uint64_t const cigar_offset,
    uint64_t const cigar_length,const bool verbose) {
  // Traverse CIGAR
  cigar_element_t* const cigar_base = vector_get_elm(cigar_vector,cigar_offset,cigar_element_t);
  uint64_t read_pos=0, text_pos=0;
  uint64_t i;
  for (i=0;i<cigar_length;++i) {
    cigar_element_t* const cigar_element = cigar_base + i;
    switch (cigar_element->type) {
      case cigar_match: {
        // Check all matching characters
        uint64_t j;
        for (j=0;j<cigar_element->length;++j) {
          if (key[read_pos] != text[text_pos]) {
            if (verbose) {
              fprintf(stream,"Align Check. Alignment not matching (key[%lu]=%c != text[%lu]=%c)\n",
                  read_pos,dna_decode(key[read_pos]),text_pos,dna_decode(text[text_pos]));
            }
            return false;
          }
          ++read_pos;
          ++text_pos;
        }
        break;
      }
      case cigar_mismatch:
        // Check mismatch
        if (key[read_pos] == text[text_pos]) {
          if (verbose) {
            fprintf(stream,"Align Check. Alignment not mismatching (key[%lu]=%c == text[%lu]=%c, CIGAR=%c)\n",
              read_pos,dna_decode(key[read_pos]),text_pos,dna_decode(text[text_pos]),dna_decode(cigar_element->mismatch));
          }
          return false;
        } else if (cigar_element->mismatch != text[text_pos]) {
          if (verbose) {
            fprintf(stream,"Align Check. Alignment not mismatching as CIGAR states (key[%lu]=%c == text[%lu]=%c, CIGAR=%c)\n",
              read_pos,dna_decode(key[read_pos]),text_pos,dna_decode(text[text_pos]),dna_decode(cigar_element->mismatch));
          }
          return false;
        }
        ++read_pos;
        ++text_pos;
        break;
      case cigar_ins:
        text_pos += cigar_element->length;
        break;
      case cigar_del:
        read_pos += cigar_element->length;
        break;
      case cigar_null:
        gem_cond_error_msg(verbose,"Align Check. CIGAR Null");
        return false;
        break;
      default:
        break;
    }
  }
  // Check alignment length
  if (read_pos != key_length) {
    if (verbose) {
      fprintf(stream,"Align Check. Alignment incorrect length (key-aligned=%lu,key-length=%lu)\n",read_pos,key_length);
    }
    return false;
  }
  if (text_pos != text_length) {
    if (verbose) {
      fprintf(stream,"Align Check. Alignment incorrect length (text-aligned=%lu,text-length=%lu)\n",text_pos,text_length);
    }
    return false;
  }
  return true;
}
/*
 * Levenshtein Align  [Dynamic Programming]
 */
#define GT_DP(i,j) dp_array[(i)*key_len+(j)]
GEM_INLINE void align_levenshtein_dp_matrix_print(
    uint64_t* const dp_array,const uint64_t key_len,const uint64_t text_len,
    const uint64_t key_limit,const uint64_t text_limit) {
  uint64_t i, j;
  for (j=0;j<key_limit;++j) {
    for (i=0;i<text_limit;++i) {
      fprintf(stdout,"%02"PRIu64" ",GT_DP(i,j));
    }
    fprintf(stdout,"\n");
  }
  fprintf(stdout,"\n");
}
GEM_INLINE int64_t align_levenshtein_get_distance(
    const char* const key,const uint64_t key_length,
    const char* const text,const uint64_t text_length,
    const bool ends_free,uint64_t* const position) {
  GEM_CHECK_NULL(key); GEM_CHECK_ZERO(key_length);
  GEM_CHECK_NULL(text); GEM_CHECK_ZERO(text_length);
  // Allocate DP-matrix
  const uint64_t key_len = key_length+1;
  const uint64_t text_len = text_length+1;
  uint64_t* dp_array[2];
  dp_array[0] = mm_calloc(2*key_len,uint64_t,false);
  dp_array[1] = dp_array[0] + key_len;
  // Init DP-Matrix
  uint64_t min_val = ALIGN_DISTANCE_INF, i_pos = ALIGN_COLUMN_INF;
  uint64_t i, j, idx_a=0, idx_b=0;
  for (j=0;j<key_len;++j) dp_array[0][j]=j;
  // Calculate DP-Matrix
  for (i=1;i<text_len;++i) {
    // Fix indexes
    idx_a = idx_b;
    idx_b = i % 2;
    // Fix first cell
    dp_array[idx_b][0] = (ends_free) ? 0 : dp_array[idx_a][0]+1;
    // Develop row
    for (j=1;j<key_len;++j) {
      const uint64_t ins = dp_array[idx_a][j]   + 1;
      const uint64_t del = dp_array[idx_b][j-1] + 1;
      const uint64_t sub = dp_array[idx_a][j-1] + ((text[i-1]==key[j-1]) ? 0 : 1);
      dp_array[idx_b][j] = MIN(sub,MIN(ins,del));
    }
    // Check last cell value
    if (ends_free && dp_array[idx_b][key_length] < min_val) {
      min_val = dp_array[idx_b][key_length];
      i_pos = i;
    }
  }
  // Return results & Free
  int64_t distance = INT64_MAX;
  if (ends_free) {
    *position = i_pos-1;
    distance = min_val;
  } else {
    *position = key_length;
    distance = dp_array[idx_b][key_length];
  }
  mm_free(dp_array[0]);
  return distance;
}
/*
 * Init SWG Query Profile
 */
GEM_INLINE void swg_allocate_query_profile_uint8(
    swg_query_profile_t* const swg_query_profile,const swg_penalties_t* swg_penalties,
    const uint64_t max_expected_key_length,mm_stack_t* const mm_stack) {
  // Compute sizes 8-bits cell
  const uint64_t num_segments_uint8 = UINT128_SIZE/UINT8_SIZE;
  const uint64_t segment_length_uint8 = DIV_CEIL(max_expected_key_length,num_segments_uint8);
  const uint64_t key_effective_length_uint8 = segment_length_uint8*num_segments_uint8;
  // Allocate 8-bits cell
  mm_stack_skip_align(mm_stack,UINT128_SIZE); // Align mm-stack & allocate
  uint8_t enc;
  for (enc=0;enc<DNA__N_RANGE;++enc) {
    swg_query_profile->query_profile_uint8[enc] = mm_stack_calloc(mm_stack,key_effective_length_uint8,uint8_t,true);
  }
}
GEM_INLINE void swg_allocate_query_profile_int16(
    swg_query_profile_t* const swg_query_profile,const swg_penalties_t* swg_penalties,
    const uint64_t max_expected_key_length,mm_stack_t* const mm_stack) {
  // Compute sizes 16-bits cell
  const uint64_t num_segments_int16 = UINT128_SIZE/UINT16_SIZE;
  const uint64_t segment_length_int16 = DIV_CEIL(max_expected_key_length,num_segments_int16);
  const uint64_t key_effective_length_uint16 = segment_length_int16*num_segments_int16;
  // Allocate 16-bits cell
  mm_stack_skip_align(mm_stack,UINT128_SIZE); // Align mm-stack & allocate
  uint8_t enc;
  for (enc=0;enc<DNA__N_RANGE;++enc) {
    swg_query_profile->query_profile_int16[enc] = (int16_t*) mm_stack_calloc(mm_stack,key_effective_length_uint16,uint8_t,true);
  }
}
GEM_INLINE void swg_init_query_profile(
    swg_query_profile_t* const swg_query_profile,const swg_penalties_t* swg_penalties,
    const uint64_t max_expected_key_length,mm_stack_t* const mm_stack) {
  // Penalties
  const int32_t generic_match_score = swg_penalties->generic_match_score;
  const int32_t generic_mismatch_score = swg_penalties->generic_mismatch_score;
  // Compute the max/min scores
  const int32_t min_profile_score = MIN(generic_match_score,generic_mismatch_score);
  const int32_t match_bias_uint8 = -min_profile_score;
  swg_query_profile->match_bias_uint8 = match_bias_uint8;
  // Allocate
  swg_allocate_query_profile_uint8(swg_query_profile,swg_penalties,max_expected_key_length,mm_stack);
  swg_allocate_query_profile_int16(swg_query_profile,swg_penalties,max_expected_key_length,mm_stack);
}
/*
 * Compile SWG Query Profile
 */
GEM_INLINE bool swg_compile_query_profile_uint8(
    swg_query_profile_t* const swg_query_profile,const swg_penalties_t* swg_penalties,
    const uint8_t* const key,const uint64_t key_length,mm_stack_t* const mm_stack) {
  // (Re)Compute sizes 8-bits cell
  const uint64_t num_segments_uint8 = UINT128_SIZE/UINT8_SIZE;
  const uint64_t segment_length_uint8 = DIV_CEIL(key_length,num_segments_uint8);
  const uint64_t key_effective_length_uint8 = segment_length_uint8*num_segments_uint8;
  swg_query_profile->segment_length_uint8 = segment_length_uint8;
  swg_query_profile->key_effective_length_uint8 = key_effective_length_uint8;
  // Compute 8b-profile bounds
  const int32_t gap_open_score = swg_penalties->gap_open_score;
  const int32_t gap_extension_score = swg_penalties->gap_extension_score;
  const int32_t match_bias_uint8 = swg_query_profile->match_bias_uint8;
  const int32_t matrix_bias_uint8 = -(gap_open_score + (int32_t)key_length*gap_extension_score);
  swg_query_profile->overflow_uint8 = (matrix_bias_uint8 + match_bias_uint8 > 255);
  if (swg_query_profile->overflow_uint8) return false;
  swg_query_profile->matrix_bias_uint8 = matrix_bias_uint8; // Lower-bound (full key deletion)
  // Compute the 8b-profile
  const uint8_t bias = swg_query_profile->match_bias_uint8;
  uint8_t enc;
  for (enc=0;enc<DNA__N_RANGE;++enc) {
    uint64_t mod, q, p;
    for (mod=0,q=0;mod<segment_length_uint8;++mod) {
      p = mod;
      while (p < key_length) {
        swg_query_profile->query_profile_uint8[enc][q++] = swg_penalties->matching_score[enc][key[p]] + bias;
        p += segment_length_uint8;
      }
      while (p < key_effective_length_uint8) {
        swg_query_profile->query_profile_uint8[enc][q++] = 0;
        p += segment_length_uint8;
      }
    }
  }
  // Return OK
  return true;
}
GEM_INLINE bool swg_compile_query_profile_int16(
    swg_query_profile_t* const swg_query_profile,const swg_penalties_t* swg_penalties,
    const uint8_t* const key,const uint64_t key_length,mm_stack_t* const mm_stack) {
  // (Re)Compute sizes 16-bits cell
  const uint64_t num_segments_int16 = UINT128_SIZE/UINT16_SIZE;
  const uint64_t segment_length_int16 = DIV_CEIL(key_length,num_segments_int16);
  const uint64_t key_effective_length_uint16 = segment_length_int16*num_segments_int16;
  swg_query_profile->segment_length_int16 = segment_length_int16;
  swg_query_profile->key_effective_length_int16 = key_effective_length_uint16;
  // Compute the 16b-profile
  uint8_t enc;
  for (enc=0;enc<DNA__N_RANGE;++enc) {
    uint64_t mod, q, p;
    for (mod=0,q=0;mod<segment_length_int16;++mod) {
      p = mod;
      while (p < key_length) {
        swg_query_profile->query_profile_int16[enc][q++] = swg_penalties->matching_score[enc][key[p]];
        p += segment_length_int16;
      }
      while (p < key_effective_length_uint16) {
        swg_query_profile->query_profile_int16[enc][q++] = 0;
        p += segment_length_int16;
      }
    }
  }
  // Return OK
  return true;
}
/*
 * SWG Score
 */
GEM_INLINE int32_t swg_score_deletion(const swg_penalties_t* const swg_penalties,const int32_t length) {
  const int32_t gap_open_score = swg_penalties->gap_open_score;
  const int32_t gap_extension = swg_penalties->gap_extension_score;
  return gap_open_score + gap_extension*length;
}
GEM_INLINE int32_t swg_score_insertion(const swg_penalties_t* const swg_penalties,const int32_t length) {
  const int32_t gap_open_score = swg_penalties->gap_open_score;
  const int32_t gap_extension = swg_penalties->gap_extension_score;
  return gap_open_score + gap_extension*length;
}
GEM_INLINE int32_t swg_score_mismatch(const swg_penalties_t* const swg_penalties) {
  return swg_penalties->generic_mismatch_score;
}
GEM_INLINE int32_t swg_score_match(const swg_penalties_t* const swg_penalties,const int32_t match_length) {
  return swg_penalties->generic_match_score * match_length;
}
GEM_INLINE int32_t swg_score_cigar_element(
    const swg_penalties_t* const swg_penalties,const cigar_element_t* const cigar_element) {
  switch (cigar_element->type) {
    case cigar_match:
      return swg_score_match(swg_penalties,cigar_element->length);
      break;
    case cigar_mismatch:
      return swg_score_mismatch(swg_penalties);
      break;
    case cigar_ins:
      return swg_score_insertion(swg_penalties,cigar_element->length);
      break;
    case cigar_del:
      return swg_score_deletion(swg_penalties,cigar_element->length);
      break;
    case cigar_null:
    default:
      GEM_INVALID_CASE();
      break;
  }
  return 0;
}
GEM_INLINE int32_t swg_score_cigar(
    const swg_penalties_t* const swg_penalties,vector_t* const cigar_vector,
    const uint64_t cigar_offset,const uint64_t cigar_length) {
  int32_t score = 0;
  // Traverse all CIGAR elements
  const cigar_element_t* const cigar_buffer = vector_get_elm(cigar_vector,cigar_offset,cigar_element_t);
  uint64_t i;
  for (i=0;i<cigar_length;++i) {
    score += swg_score_cigar_element(swg_penalties,cigar_buffer+i);
  }
  // Return score
  return score;
}
/*
 * SWG - Debug
 */
GEM_INLINE void swg_align_match_table_print(swg_cell_t** const dp,const uint64_t num_columns,const uint64_t num_rows) {
  uint64_t i, j;
  for (i=0;i<=num_rows;++i) {
    for (j=0;j<=num_columns;++j) {
      printf("%+4d ",dp[j][i].M);
    }
    printf("\n");
  }
}
GEM_INLINE void swg_align_match_input_print(
    const uint8_t* const key,const uint64_t key_length,
    const uint8_t* const text,const uint64_t text_length) {
  uint64_t i;
  printf("KEY> ");
  for (i=0;i<key_length;++i) printf("%c",dna_decode(key[i]));
  printf("\n");
  printf("TXT> ");
  for (i=0;i<text_length;++i) printf("%c",dna_decode(text[i]));
  printf("\n");
}
/*
 * SWG - BackTrace
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text
 *   @align_input->text_length
 *   @match_alignment->match_position (Correction)
 *   @match_alignment->cigar_length (Cumulative)
 */
GEM_INLINE void swg_align_match_traceback(
    match_align_input_t* const align_input,swg_cell_t** const dp,
    const int32_t max_score,const uint64_t max_score_column,
    const int32_t single_gap,const int32_t gap_extension,const bool begin_free,
    match_alignment_t* const match_alignment,vector_t* const cigar_vector) {
  // Parameters
  const uint8_t* const key = align_input->key;
  const uint64_t key_length = align_input->key_length;
  uint8_t* const text = align_input->text;
  // Allocate CIGAR string memory (worst case)
  vector_reserve_additional(cigar_vector,key_length); // Reserve
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
GEM_INLINE swg_cell_t** swg_align_match_allocate_table(
    const uint64_t num_columns,const uint64_t num_rows,mm_stack_t* const mm_stack) {
  swg_cell_t** const dp = mm_stack_malloc(mm_stack,num_columns*sizeof(swg_cell_t*));
  const uint64_t row_size = num_rows*sizeof(swg_cell_t);
  uint64_t column;
  for (column=0;column<num_columns;++column) {
    dp[column] = mm_stack_malloc(mm_stack,row_size);
  }
  return dp;
}
GEM_INLINE void swg_align_match_init_table(
    swg_cell_t** const dp,const uint64_t num_columns,const uint64_t num_rows,
    const int32_t single_gap,const int32_t gap_extension) {
  uint64_t column, row;
  for (column=0;column<num_columns;++column) {
    dp[column][0].D = SWG_SCORE_INT32_MIN;
    dp[column][0].M = 0;
  }
  dp[0][0].I = SWG_SCORE_INT32_MIN;
  dp[0][1].I = SWG_SCORE_INT32_MIN;
  dp[0][1].M = single_gap; // g(1)
  for (row=2;row<num_rows;++row) {
    dp[0][row].I = SWG_SCORE_INT32_MIN;
    dp[0][row].M = dp[0][row-1].M + gap_extension; // g(row)
  }
}
GEM_INLINE void swg_align_match_init_table_banded(
    swg_cell_t** const dp,const uint64_t num_columns,const uint64_t num_rows,
    const uint64_t column_start_band,const uint64_t band_low_offset,
    const bool begin_free,const int32_t single_gap,const int32_t gap_extension) {
  uint64_t column, row;
  // Initialize first column
  dp[0][0].D = SWG_SCORE_INT32_MIN; // Not used
  dp[0][0].I = SWG_SCORE_INT32_MIN; // Not used
  dp[0][0].M = 0;
  dp[0][1].M = single_gap; // g(1)
  dp[0][1].I = SWG_SCORE_INT32_MIN;
  for (row=2;row<band_low_offset;++row) {
    dp[0][row].M = dp[0][row-1].M + gap_extension; // g(row)
    dp[0][row].I = SWG_SCORE_INT32_MIN;
  }
  // Initialize first row
  if (begin_free) {
    for (column=1;column<num_columns;++column) {
      dp[column][0].M = 0;
      dp[column][0].D = SWG_SCORE_INT32_MIN;
    }
  } else {
    dp[1][0].M = single_gap;
    dp[1][0].D = SWG_SCORE_INT32_MIN;
    for (column=2;column<column_start_band;++column) {
      dp[column][0].M = dp[column-1][0].M + gap_extension;
      dp[column][0].D = SWG_SCORE_INT32_MIN;
    }
  }
}
GEM_INLINE void swg_align_match_init_table_banded_opt(
    swg_cell_t** const dp,const uint64_t num_columns,const uint64_t num_rows,
    const uint64_t column_start_band,const uint64_t band_low_offset,
    const bool begin_free,const int32_t single_gap,const int32_t gap_extension) {
  dp[0][0].D = SWG_SCORE_INT32_MIN; // Not used
  dp[0][0].I = SWG_SCORE_INT32_MIN; // Not used
  dp[0][0].M = 0;
  dp[0][1].M = single_gap; // g(1)
  dp[0][1].I = SWG_SCORE_INT32_MIN;
  dp[1][1].I = SWG_SCORE_INT32_MIN;
  uint64_t column, row;
  // Initialize first/second column
  for (row=2;row<band_low_offset;++row) {
    dp[0][row].M = dp[0][row-1].M + gap_extension; // g(row)
    dp[0][row].I = SWG_SCORE_INT32_MIN;
    dp[1][row].I = dp[0][row].M + single_gap;
  }
  // Initialize first/second row
  if (begin_free) {
    for (column=1;column<num_columns;++column) {
      dp[column][0].M = 0;
      dp[column][0].D = SWG_SCORE_INT32_MIN;
      dp[column][1].D = dp[column][0].M + single_gap;
    }
  } else {
    dp[1][0].M = single_gap;
    dp[1][0].D = SWG_SCORE_INT32_MIN;
    dp[1][1].D = single_gap + single_gap;
    for (column=2;column<column_start_band;++column) {
      dp[column][0].M = dp[column-1][0].M + gap_extension;
      dp[column][0].D = SWG_SCORE_INT32_MIN;
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
GEM_INLINE void swg_align_match_base(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    match_alignment_t* const match_alignment,vector_t* const cigar_vector,mm_stack_t* const mm_stack) {
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
  swg_cell_t** const dp = swg_align_match_allocate_table(num_columns,num_rows,mm_stack);
  // Initialize first row
  const matching_score_t* const matching_score = &swg_penalties->matching_score;
  const int32_t single_gap = swg_penalties->gap_open_score + swg_penalties->gap_extension_score; // g(1)
  const int32_t gap_extension = swg_penalties->gap_extension_score;
  swg_align_match_init_table(dp,num_columns,num_rows,single_gap,gap_extension);
  // Compute DP-matrix
  int32_t max_score = SWG_SCORE_INT32_MIN;
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
    // swg_align_match_table_print(dp,MIN(10,key_length),MIN(10,key_length));
  // Init Match
  match_alignment->cigar_length = 0;
  match_alignment->effective_length = 0;
  match_alignment->score = 0;
  // Retrieve the alignment. Store the match (Backtrace and generate CIGAR)
  swg_align_match_traceback(align_input,dp,max_score,max_score_column,
      single_gap,gap_extension,true,match_alignment,cigar_vector);
  // Clean-up
  mm_stack_pop_state(mm_stack,false); // Free
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
GEM_INLINE void swg_align_match_full(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    const bool begin_free,const bool end_free,match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,mm_stack_t* const mm_stack) {
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
  swg_cell_t** const dp = swg_align_match_allocate_table(num_columns+1,num_rows+1,mm_stack);
  // Initialize DP-matrix
  const matching_score_t* const matching_score = &swg_penalties->matching_score;
  const int32_t gap_extension = swg_penalties->gap_extension_score;
  const int32_t gap_open = swg_penalties->gap_open_score;
  const int32_t single_gap = gap_open + gap_extension; // g(1)
  swg_align_match_init_table_banded_opt(dp,num_columns,num_rows,
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
  swg_align_match_traceback(align_input,dp,max_score,max_score_column,
      single_gap,gap_extension,begin_free,match_alignment,cigar_vector);
  // Clean-up
  mm_stack_pop_state(mm_stack,false); // Free
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
GEM_INLINE void swg_align_match_banded(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    const bool begin_free,const bool end_free,match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,mm_stack_t* const mm_stack) {
  // Parameters
  const uint8_t* const key = align_input->key;
  const uint64_t key_length = align_input->key_length;
  uint8_t* const text = align_input->text;
  uint64_t text_length = align_input->text_length;
  const swg_penalties_t* swg_penalties = align_parameters->swg_penalties;
  const uint64_t max_bandwidth = align_parameters->max_bandwidth;
  // Initialize band-limits
  if (text_length > key_length + max_bandwidth) { // Text too long for band
    if (!begin_free && !end_free) { match_alignment->score = ALIGN_DISTANCE_INF; return; }
    if (!begin_free) text_length = key_length + max_bandwidth;
  }
  if (text_length + max_bandwidth <= key_length) { // Text too short for band
    match_alignment->score = ALIGN_DISTANCE_INF; return;
  }
  PROF_START(GP_SWG_ALIGN_BANDED);
  // Allocate memory
  mm_stack_push_state(mm_stack); // Save stack state
  const uint64_t num_rows = (key_length+1);
  const uint64_t num_rows_1 = num_rows-1;
  const uint64_t num_columns = (text_length+1);
  swg_cell_t** const dp = swg_align_match_allocate_table(num_columns+1,num_rows+1,mm_stack);
  // Initialize DP-matrix
  const matching_score_t* const matching_score = &swg_penalties->matching_score;
  const int32_t gap_extension = swg_penalties->gap_extension_score;
  const int32_t gap_open = swg_penalties->gap_open_score;
  const int32_t single_gap = gap_open + gap_extension; // g(1)
  // Initialize band
  int64_t column_start_band = max_bandwidth + 2;
  if (begin_free || column_start_band > text_length) column_start_band = text_length + 2;
  uint64_t band_high_offset = 0;
  uint64_t band_low_offset = (max_bandwidth+1 < num_rows) ? max_bandwidth+2 : num_rows;
  // Initialize first/second column
  swg_align_match_init_table_banded_opt(
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
      dp[column][band_low_offset].I = SWG_SCORE_INT32_MIN;
      dp[column-1][band_low_offset].I = SWG_SCORE_INT32_MIN;
      ++band_low_offset; // Swift band
    }
    // Initialize band boundaries & update band limits
    if (column >= column_start_band) {
      ++band_high_offset; // Swift band
      dp[column][band_high_offset+1].D = SWG_SCORE_INT32_MIN;
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
  swg_align_match_traceback(align_input,dp,max_score,max_score_column,
      single_gap,gap_extension,begin_free,match_alignment,cigar_vector);
  // Clean-up
  mm_stack_pop_state(mm_stack,false); // Free
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
GEM_INLINE void swg_align_match(
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    const bool begin_free,const bool end_free,match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,mm_stack_t* const mm_stack) {
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
      match_alignment->score = swg_score_insertion(swg_penalties,text_length);
      match_alignment->effective_length = text_length;
    } else {
      match_alignment->score = 0;
      match_alignment->effective_length = 0;
    }
  } else if (key_length > 0 && text_length == 0) {
    // Deletion <-@key_length>
    matches_cigar_vector_append_deletion(cigar_vector,&match_alignment->cigar_length,key_length,cigar_attr_none);
    match_alignment->score = swg_score_deletion(swg_penalties,key_length);
  } else if (key_length==1 && text_length==1) {
    // Mismatch/Match
    const uint8_t key_enc = key[0];
    const uint8_t text_enc = text[0];
    if (!allowed_enc[text_enc] || text_enc != key_enc) {
      matches_cigar_vector_append_mismatch(cigar_vector,&match_alignment->cigar_length,text_enc,cigar_attr_none);
      match_alignment->score = swg_score_mismatch(swg_penalties);
    } else {
      matches_cigar_vector_append_match(cigar_vector,&match_alignment->cigar_length,1,cigar_attr_none);
      match_alignment->score = swg_score_match(swg_penalties,1);
    }
    match_alignment->effective_length = 1;
  } else {
    PROF_ADD_COUNTER(GP_SWG_ALIGN_BANDED_LENGTH,text_length);
    swg_align_match_banded(align_input,align_parameters,
        begin_free,end_free,match_alignment,cigar_vector,mm_stack);
  }
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

//GEM_INLINE void swg_align_match_table_print_int16(
//    int16_t** const dp_m,const uint64_t num_columns,
//    const uint64_t num_rows,const uint64_t segment_length) {
//  const uint64_t num_segments = UINT128_SIZE/UINT16_SIZE; // Vector length
//  uint64_t i, j;
//  for (i=0;i<=num_rows;++i) {
//    const uint64_t position = (i%segment_length)*num_segments+i/segment_length;
//    for (j=0;j<=num_columns;++j) {
//      printf("%4d ",(int)dp_m[j][position]);
//    }
//    printf("\n");
//  }
//}
//GEM_INLINE int16_t** swg_align_match_allocate_table_int16(
//    const uint64_t num_columns,const uint64_t num_rows,mm_stack_t* const mm_stack) {
//  // Allocate the pointers
//  int16_t** const dp = mm_stack_malloc(mm_stack,num_columns*sizeof(int16_t*));
//  // Allocate all the columns
//  mm_stack_skip_align(mm_stack,UINT128_SIZE); // Align Stack Memory
//  const uint64_t row_size = num_rows*sizeof(int16_t);
//  uint64_t column;
//  for (column=0;column<num_columns;++column) {
//    dp[column] = mm_stack_malloc(mm_stack,row_size);
//  }
//  return dp;
//}
//GEM_INLINE void swg_align_match_init_table_int16(
//    int16_t** const dp_M,int16_t** const dp_I,const uint64_t num_columns,
//    const uint64_t num_rows,const uint64_t key_length,const bool begin_free,
//    const uint16_t single_gap,const uint16_t gap_extension) {
//  uint64_t row;
//  // Initialize first column
//  dp_M[0][0] = -single_gap;
//  dp_I[0][0] = INT16_MIN; // Not used
//  dp_I[1][0] = dp_M[0][0] - single_gap;
//  for (row=1;row<key_length;++row) {
//    dp_M[0][row] = dp_M[0][row-1] - gap_extension; // g(row)
//    dp_I[0][row] = INT16_MIN;
//    dp_I[1][row] = dp_M[0][row] - single_gap;
//  }
//  for (;row<num_rows;++row) {
//    dp_M[0][row] = INT16_MIN;
//    dp_I[0][row] = INT16_MIN;
//    dp_I[1][row] = INT16_MIN;
//  }
//}

#define _mm_extract_max_epi16(result,v128) \
  v128 = _mm_max_epi16(v128,_mm_srli_si128(v128, 8)); \
  v128 = _mm_max_epi16(v128,_mm_srli_si128(v128, 4)); \
  v128 = _mm_max_epi16(v128,_mm_srli_si128(v128, 2)); \
  result = _mm_extract_epi16(v128, 0)

//GEM_INLINE void swg_align_match_int16_simd128(
//    const uint8_t* const key,const uint64_t key_length,swg_query_profile_t* const swg_query_profile,
//    const swg_penalties_t* swg_penalties,uint64_t* const match_position,uint8_t* const text,
//    const uint64_t text_length,const bool begin_free,const bool end_free,vector_t* const cigar_vector,
//    uint64_t* const cigar_length,int64_t* const effective_length,int32_t* const alignment_score,
//    mm_stack_t* const mm_stack) {
//  // Compile query profile
//  swg_compile_query_profile_int16(swg_query_profile,swg_penalties,key,key_length,mm_stack);
//  // Allocate memory
//  mm_stack_push_state(mm_stack); // Save stack state
//  const uint64_t num_segments_16b = UINT128_SIZE/UINT16_SIZE; // Elements in each SIMD-vectors
//  const uint64_t segment_length = swg_query_profile->segment_length_int16; // Number of SIMD-vectors
//  const uint64_t key_effective_length = swg_query_profile->key_effective_length_int16; // Number of SIMD-vectors
//  const uint64_t num_columns = (text_length+1);
//  const uint64_t num_rows = key_effective_length;
//  int16_t** const dp_M = swg_align_match_allocate_table_int16(num_columns,num_rows,mm_stack);
//  int16_t** const dp_I = swg_align_match_allocate_table_int16(num_columns+1,num_rows,mm_stack);
//  int16_t** const dp_D = swg_align_match_allocate_table_int16(num_columns,num_rows,mm_stack);
//  // Initialize DP-matrix
//  const int16_t gap_extension = (-swg_penalties->gap_extension_score);
//  const int16_t gap_open = (-swg_penalties->gap_open_score);
//  const int16_t single_gap = gap_open + gap_extension; // g(1)
//  swg_align_match_init_table_int16(dp_M,dp_I,num_columns,num_rows,key_length,begin_free,single_gap,gap_extension);
//  __m128i* M_in = (__m128i*)dp_M[0], *M_out;
//  __m128i* I_in = (__m128i*)dp_I[1], *I_out;
//  __m128i* D_out;
//  __m128i* current_profile;
//  __m128i zero_v = _mm_set1_epi32(INT16_MIN);
//  __m128i gap_extension_v = _mm_set1_epi16(gap_extension);
//  __m128i single_gap_v = _mm_set1_epi16(single_gap);
//  __m128i max_score_v = zero_v, Mv, Iv, Dv, vTmp;
//  // Traverse the text
//  uint64_t i, j, max_score_column = 0;
//  int16_t max_score_local, max_score_global = 0, Mv_0 = 0;
//  for(j=0;j<text_length;j++) {
//    // Load M-in/M-out
//    M_out = (__m128i*)dp_M[j+1]; // Set M-out memory column
//    D_out = (__m128i*)dp_D[j+1]; // Set D-out memory column
//    I_out = (__m128i*)dp_I[j+2]; // Set I-out memory column
//    Dv = zero_v; // Set D-column to zero (later on corrected)
//    int16_t* Mv_mem = (int16_t*)(M_in+segment_length-1);
//    Mv = _mm_set_epi16(Mv_mem[6], Mv_mem[5], Mv_mem[4], Mv_mem[3],
//                       Mv_mem[2], Mv_mem[1], Mv_mem[0], Mv_0);
//    if (!begin_free) Mv_0 = Mv_0 - ((j==0) ? single_gap : gap_extension);
//    // Load current profile
//    current_profile = (__m128i*)swg_query_profile->query_profile_uint8[text[j]];
//    for (i=0;i<segment_length;i++) {
//      // Add the profile to Mv
//      Mv = _mm_adds_epi16(Mv,*(current_profile+i));
//      // Update Mv
//      Iv = _mm_load_si128(I_in+i); // Load Iv
//      Mv = _mm_max_epi16(Mv,Dv);
//      Mv = _mm_max_epi16(Mv,Iv);
//      max_score_v = _mm_max_epi16(max_score_v,Mv); // Pre-calculate max-score (gaps only decrease score)
//      _mm_store_si128(M_out+i,Mv); // Store Mv_i
//      // Update next Iv
//      Iv = _mm_subs_epu16(Iv,gap_extension_v);
//      vTmp = _mm_subs_epu16(Mv,single_gap_v);
//      Iv = _mm_max_epi16(Iv,vTmp);
//      _mm_store_si128(I_out+i,Iv); // Store Iv_i
//      // Update next Dv (partially)
//      Dv = _mm_subs_epu16(Dv,gap_extension_v);
//      Dv = _mm_max_epi16(Dv,vTmp);
//      _mm_store_si128(D_out+i,Dv); // Store Dv_i
//      // Load Mv for next iteration
//      Mv = _mm_load_si128(M_in+i);
//    }
//    // Lazy-F Loop
//    for (i=0;i<num_segments_16b;++i) {
//      uint64_t k;
//      // Compute the gap extend penalty for the current cell
//      Dv = _mm_slli_si128(Dv,2);
//      for (k=0;k<segment_length;++k) {
//        // Compute the current optimal value of the cell
//        vTmp = _mm_load_si128(M_out+k);
//        Mv = _mm_max_epi16(vTmp,Dv);
//        _mm_store_si128(M_out+k,Mv);
//        Mv = _mm_subs_epu16(Mv,single_gap_v);
//        Dv = _mm_subs_epu16(Mv,gap_extension_v);
//        // Check Mv unchanged
//        if(!_mm_movemask_epi8(_mm_cmpgt_epi16(Mv,Dv))) {
//          goto exit_lazy_f_loop;
//        }
//        // Compute the scores for the next cell
//        _mm_store_si128(D_out+k,Dv);
//      }
//    }
//exit_lazy_f_loop: // Exit Lazy-F Loop
//    // Max score in the column
//    _mm_extract_max_epi16(max_score_local,max_score_v);
//    if (max_score_local >= max_score_global) {
//      max_score_global = max_score_local;
//      max_score_column = j;
//    }
//    // Swap M memory
//    M_in = M_out;
//    I_in = I_out;
//  }
//  // DEBUG
//  swg_align_match_table_print_int16(dp_M,MIN(10,key_length),MIN(10,key_length),segment_length);
//  // Clean & return
//  mm_stack_pop_state(mm_stack,false); // Free
//  *alignment_score = max_score_column; // FIXME
//  *alignment_score = max_score_global; // FIXME
//}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

//GEM_INLINE void swg_align_match_table_print_uint8(
//    uint8_t** const dp_m,const uint64_t num_columns,const uint64_t num_rows,
//    const uint64_t segment_length,const uint64_t matrix_bias) {
//  const uint64_t num_segments = UINT128_SIZE; // Vector length
//  uint64_t i, j;
//  for (i=0;i<=num_rows;++i) {
//    const uint64_t position = (i%segment_length)*num_segments+i/segment_length;
//    for (j=0;j<=num_columns;++j) {
//      printf("%+4d ",(int)dp_m[j][position]-((int)matrix_bias));
//    }
//    printf("\n");
//  }
//}
//GEM_INLINE uint8_t** swg_align_match_allocate_table_uint8(
//    const uint64_t num_columns,const uint64_t num_rows,mm_stack_t* const mm_stack) {
//  // Allocate the pointers
//  uint8_t** const dp = mm_stack_malloc(mm_stack,num_columns*sizeof(uint8_t*));
//  // Allocate all the columns
//  mm_stack_skip_align(mm_stack,UINT128_SIZE); // Align Stack Memory
//  const uint64_t row_size = num_rows*sizeof(uint8_t);
//  uint64_t column;
//  for (column=0;column<num_columns;++column) {
//    dp[column] = mm_stack_malloc(mm_stack,row_size);
//  }
//  return dp;
//}
//GEM_INLINE void swg_align_match_init_table_uint8(
//    uint8_t** const dp_M,uint8_t** const dp_I,const uint64_t num_columns,
//    const uint64_t num_rows,const uint64_t key_length,const uint8_t matrix_bias,
//    const bool begin_free,const uint8_t single_gap,const uint8_t gap_extension) {
//  uint64_t row;
//  // Initialize first column
//  dp_M[0][0] = matrix_bias - single_gap;
//  dp_I[0][0] = 0; // Not used
//  dp_I[1][0] = dp_M[0][0] > single_gap ? dp_M[0][0]-single_gap : 0;
//  for (row=1;row<key_length;++row) {
//    dp_M[0][row] = dp_M[0][row-1] - gap_extension; // g(row)
//    dp_I[0][row] = 0;
//    dp_I[1][row] = dp_M[0][row] > single_gap ? dp_M[0][row]-single_gap : 0;
//  }
//  for (;row<num_rows;++row) {
//    dp_M[0][row] = 0;
//    dp_I[0][row] = 0;
//    dp_I[1][row] = 0;
//  }
//}
#define _mm_extract_max_epu8(result,v128) \
  v128 = _mm_max_epu8(v128,_mm_srli_si128(v128,8)); \
  v128 = _mm_max_epu8(v128,_mm_srli_si128(v128,4)); \
  v128 = _mm_max_epu8(v128,_mm_srli_si128(v128,2)); \
  v128 = _mm_max_epu8(v128,_mm_srli_si128(v128,1)); \
  result = _mm_extract_epi16(v128,0) & 0x00FF
//GEM_INLINE bool swg_align_match_uint8_simd128(
//    const uint8_t* const key,const uint64_t key_length,swg_query_profile_t* const swg_query_profile,
//    const swg_penalties_t* swg_penalties,uint64_t* const match_position,uint8_t* const text,
//    const uint64_t text_length,const bool begin_free,const bool end_free,vector_t* const cigar_vector,
//    uint64_t* const cigar_length,int64_t* const effective_length,int32_t* const alignment_score,
//    mm_stack_t* const mm_stack) {
//  // Compile query profile
//  if (!swg_compile_query_profile_uint8(swg_query_profile,swg_penalties,key,key_length,mm_stack)) {
//    return false; // TODO
//  }
//  // Allocate memory
//  mm_stack_push_state(mm_stack); // Save stack state
//  const uint64_t num_segments_8b = UINT128_SIZE/UINT8_SIZE; // Elements in each SIMD-vectors
//  const uint64_t segment_length = swg_query_profile->segment_length_uint8; // Number of SIMD-vectors
//  const uint64_t key_effective_length = swg_query_profile->key_effective_length_uint8; // Number of SIMD-vectors
//  const uint64_t num_columns = (text_length+1);
//  const uint64_t num_rows = key_effective_length;
//  uint8_t** const dp_M = swg_align_match_allocate_table_uint8(num_columns,num_rows,mm_stack); // TODO key_length
//  uint8_t** const dp_I = swg_align_match_allocate_table_uint8(num_columns+1,num_rows,mm_stack);
//  uint8_t** const dp_D = swg_align_match_allocate_table_uint8(num_columns,num_rows,mm_stack);
//  // Initialize DP-matrix
//  const uint8_t gap_extension = (-swg_penalties->gap_extension_score);
//  const uint8_t gap_open = (-swg_penalties->gap_open_score);
//  const uint8_t single_gap = gap_open + gap_extension; // g(1)
//  const uint8_t match_bias = swg_query_profile->match_bias_uint8;
//  const uint8_t matrix_bias = swg_query_profile->matrix_bias_uint8;
//  swg_align_match_init_table_uint8(dp_M,dp_I,
//      num_columns,num_rows,key_length,matrix_bias,begin_free,single_gap,gap_extension);
//  __m128i* M_in = (__m128i*)dp_M[0], *M_out;
//  __m128i* I_in = (__m128i*)dp_I[1], *I_out;
//  __m128i* D_out;
//  __m128i* current_profile;
//  __m128i zero_v = _mm_set1_epi32(0);
//  __m128i gap_extension_v = _mm_set1_epi8(gap_extension);
//  __m128i single_gap_v = _mm_set1_epi8(single_gap);
//  __m128i bias_v = _mm_set1_epi8(match_bias);
//  __m128i max_score_v = zero_v, Mv, Iv, Dv, vTmp;
//  // Traverse the text
//  uint64_t i, j, max_score_column = 0;
//  uint8_t max_score_local, max_score_global = 0, Mv_0 = matrix_bias;
//  for(j=0;j<text_length;j++) {
//    // Load M-in/M-out
//    M_out = (__m128i*)dp_M[j+1]; // Set M-out memory column
//    D_out = (__m128i*)dp_D[j+1]; // Set D-out memory column
//    I_out = (__m128i*)dp_I[j+2]; // Set I-out memory column
//    Dv = zero_v; // Set D-column to zero (later on corrected)
//    uint8_t* Mv_mem = (uint8_t*)(M_in+segment_length-1);
//    Mv = _mm_set_epi8(
//        Mv_mem[14],Mv_mem[13],Mv_mem[12],Mv_mem[11],
//        Mv_mem[10],Mv_mem[9], Mv_mem[8], Mv_mem[7],
//        Mv_mem[6], Mv_mem[5], Mv_mem[4], Mv_mem[3],
//        Mv_mem[2], Mv_mem[1], Mv_mem[0], Mv_0);
//    if (!begin_free) Mv_0 = Mv_0 - ((j==0) ? single_gap : gap_extension);
//    // Load current profile
//    current_profile = (__m128i*)swg_query_profile->query_profile_uint8[text[j]];
//    for (i=0;i<segment_length;i++) {
//      // Add the profile to Mv
//      Mv = _mm_adds_epu8(Mv,*(current_profile+i));
//      Mv = _mm_subs_epu8(Mv,bias_v);
//      // Update Mv
//      Iv = _mm_load_si128(I_in+i); // Load Iv
//      Mv = _mm_max_epu8(Mv,Dv);
//      Mv = _mm_max_epu8(Mv,Iv);
//      max_score_v = _mm_max_epu8(max_score_v,Mv); // Pre-calculate max-score (gaps only decrease score)
//      _mm_store_si128(M_out+i,Mv); // Store Mv_i
//      // Update next Iv
//      Iv = _mm_subs_epu8(Iv,gap_extension_v);
//      vTmp = _mm_subs_epu8(Mv,single_gap_v);
//      Iv = _mm_max_epu8(Iv,vTmp);
//      _mm_store_si128(I_out+i,Iv); // Store Iv_i
//      // Update next Dv (partially)
//      Dv = _mm_subs_epu8(Dv,gap_extension_v);
//      Dv = _mm_max_epu8(Dv,vTmp);
//      _mm_store_si128(D_out+i,Dv); // Store Dv_i
//      // Load Mv for next iteration
//      Mv = _mm_load_si128(M_in+i);
//    }
//    // Lazy-F Loop
//    for (i=0;i<num_segments_8b;++i) {
//      uint64_t k;
//      // Compute the gap extend penalty for the current cell
//      Dv = _mm_slli_si128(Dv,1);
//      for (k=0;k<segment_length;++k) {
//        // Compute the current optimal value of the cell
//        vTmp = _mm_load_si128(M_out+k);
//        Mv = _mm_max_epu8(vTmp,Dv);
//        _mm_store_si128(M_out+k,Mv);
//        Mv = _mm_subs_epu8(Mv,single_gap_v);
//        Dv = _mm_subs_epu8(Mv,gap_extension_v);
//        // Check Mv unchanged
//        if (gem_expect_false(
//             _mm_movemask_epi8(
//               _mm_cmpeq_epi8(Dv,Mv)) == 0xFFFF)) {
//          goto exit_lazy_f_loop;
//        }
//        // Compute the scores for the next cell
//        _mm_store_si128(D_out+k,Dv);
//      }
//    }
//exit_lazy_f_loop: // Exit Lazy-F Loop
//    // Max score in the column
//    _mm_extract_max_epu8(max_score_local,max_score_v);
//    if (max_score_local >= max_score_global) {
//      max_score_global = max_score_local;
//      max_score_column = j;
//      if (max_score_global + match_bias >= 255) break; // Too bad
//    }
//    // Swap M memory
//    M_in = M_out;
//    I_in = I_out;
//  }
//  // DEBUG
//  swg_align_match_table_print_uint8(dp_M,MIN(10,key_length),MIN(10,key_length),segment_length,matrix_bias);
//  // Clean & return
//  mm_stack_pop_state(mm_stack,false); // Free
//  *alignment_score = max_score_column; // FIXME
//  *alignment_score = max_score_global; // FIXME
//  return ((int)max_score_global + (int)match_bias < 255);
//}
//
//
