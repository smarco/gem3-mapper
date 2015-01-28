/*
 * PROJECT: GEMMapper
 * FILE: swg_align.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "swg_align.h"
#include "bpm_align.h"
#include "matches.h"

/*
 * Constants
 */
#define SWG_SCORE_INT32_MIN (INT16_MIN)

//uint32_t sadd32(uint32_t a, uint32_t b)
//{
//#if defined IA32
//  __asm
//  {
//    mov eax,a
//    xor edx,edx
//    add eax,b
//    setnc dl
//    dec edx
//    or eax,edx
//  }
//#elif defined ARM
//  // ARM code
//#else
//  return (a > 0xFFFFFFFF - b) ? 0xFFFFFFFF : a + b;
//#endif
//}

/*
 * Check CIGAR string
 */
GEM_INLINE bool align_check_match(
    FILE* const stream,const uint8_t* const key,const uint64_t key_length,const uint8_t* const text,
    const uint64_t text_length,vector_t* const cigar_buffer,uint64_t const cigar_offset,
    uint64_t const cigar_length,const bool verbose) {
  // Traverse CIGAR
  cigar_element_t* const cigar_base = vector_get_elm(cigar_buffer,cigar_offset,cigar_element_t);
  uint64_t read_pos=0, text_pos=0;
  uint64_t i;
  for (i=0;i<cigar_length;++i) {
    cigar_element_t* const cigar_element = cigar_base + i;
    switch (cigar_element->type) {
      case cigar_match: {
        // Check all matching characters
        uint64_t j;
        for (j=0;j<cigar_element->match_length;++j) {
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
        text_pos += cigar_element->indel.indel_length;
        break;
      case cigar_del:
      case cigar_soft_trim:
        read_pos += cigar_element->indel.indel_length;
        break;
      case cigar_null:
        gem_cond_error_msg(verbose,"Align Check. CIGAR Null");
        return false;
        break;
      default:
        break;
    }
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
 * Compile SWG Query Profile
 */
GEM_INLINE void swg_compile_query_profile(
    swg_query_profile_t* const swg_query_profile,const swg_penalties_t* swg_penalties,
    const uint8_t* const key,const uint64_t key_length,mm_stack_t* const mm_stack) {
  // Compile the query profile
  const uint64_t num_segments_8b = UINT128_SIZE/UINT8_SIZE;
  const uint64_t num_segments_16b = UINT128_SIZE/UINT16_SIZE;
  const uint64_t segment_length_8b = DIV_CEIL(key_length,num_segments_8b);
  const uint64_t segment_length_16b = DIV_CEIL(key_length,num_segments_16b);
  uint8_t enc;
  for (enc=0;enc<DNA__N_RANGE;++enc) {
    // Allocate memory for the profile
    swg_query_profile->query_profile_uint8[enc] = mm_stack_calloc(mm_stack,key_length,uint8_t,true);
    swg_query_profile->query_profile_int16[enc] = mm_stack_calloc(mm_stack,key_length,int16_t,true); // FIXME aligned 16
    // Compile the 16b-segments
    uint64_t mod, p, q;
    for (mod=0,q=0;mod<segment_length_16b;++mod) {
      p = mod;
      while (p < key_length) {
        swg_query_profile->query_profile_int16[enc][q++] = swg_penalties->matching_score[enc][key[p]];
        p += num_segments_16b;
      }
    }
    // Compile the 8b-segments
    for (mod=0,q=0;mod<segment_length_8b;++mod) {
      p = mod;
      while (p < key_length) {
        swg_query_profile->query_profile_uint8[enc][q++] = swg_penalties->matching_score[enc][key[p]];
        p += num_segments_8b;
      }
    }
  }
}
/*
 * Smith-waterman-gotoh Alignment
 */
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
GEM_INLINE void swg_align_match_traceback_32b(
    swg_cell_t** const dp,const int32_t max_score,const uint64_t max_score_column,
    const int32_t single_gap,const int32_t gap_extension,const uint8_t* const key,
    const uint64_t key_length,uint64_t* const match_position,uint8_t* const text,
    const bool begin_free,vector_t* const cigar_buffer,uint64_t* const cigar_length,
    int64_t* const effective_length,int32_t* const alignment_score) {
  // Allocate CIGAR string memory (worst case)
  vector_reserve_additional(cigar_buffer,key_length); // Reserve
  cigar_element_t* cigar_buffer_sentinel = vector_get_free_elm(cigar_buffer,cigar_element_t); // Sentinel
  cigar_element_t* const cigar_buffer_base = cigar_buffer_sentinel;
  cigar_buffer_sentinel->type = cigar_null; // Trick
  // Start Backtrace
  int64_t match_effective_length = key_length;
  int32_t match_alignment_score = max_score;
  uint64_t h = max_score_column;
  uint64_t v = key_length;
  cigar_t traceback_matrix = cigar_match;
  while (v > 0 && h > 0) {
    switch (traceback_matrix) {
      case cigar_del:
        // Traceback D-matrix
        matches_cigar_buffer_add_cigar_element(&cigar_buffer_sentinel,cigar_del,1,NULL); // Deletion <-1>@v
        if (dp[h][v].D == dp[h][v-1].M + single_gap) traceback_matrix = cigar_match;
        --v; --match_effective_length;
        break;
      case cigar_ins:
        // Traceback I-matrix
        matches_cigar_buffer_add_cigar_element(&cigar_buffer_sentinel,cigar_ins,1,text+(h-1)); // Insertion <+1>@v
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
            if (cigar_buffer_sentinel->type!=cigar_null) ++(cigar_buffer_sentinel);
            cigar_buffer_sentinel->type = cigar_mismatch;
            cigar_buffer_sentinel->mismatch = text[h-1];
            --h; --v;
          } else {
            matches_cigar_buffer_add_cigar_element(&cigar_buffer_sentinel,cigar_match,1,NULL); // Match
            --h; --v;
          }
        }
        break;
    }
  }
  if (v > 0) {
    matches_cigar_buffer_add_cigar_element(&cigar_buffer_sentinel,cigar_del,v,NULL); // <-(@v+1)>@v
    match_effective_length -= v;
  }
  if (h > 0) {
    if (begin_free) {
      *match_position += h; // We need to correct the matching_position
    } else {
      matches_cigar_buffer_add_cigar_element(&cigar_buffer_sentinel,cigar_ins,h,text); // <-(@h+1)>@h
      match_effective_length += h;
      if (h==1) {
        match_alignment_score += single_gap;
      } else {
        match_alignment_score += single_gap + (h-1)*gap_extension;
      }
    }
  }
  // Set CIGAR buffer used
  if (cigar_buffer_sentinel->type!=cigar_null) ++(cigar_buffer_sentinel);
  const uint64_t num_cigar_elements = cigar_buffer_sentinel - cigar_buffer_base;
  vector_add_used(cigar_buffer,num_cigar_elements);
  // Reverse CIGAR Elements
  if (num_cigar_elements > 0) {
    const uint64_t middle_point = num_cigar_elements/2;
    uint64_t i;
    for (i=0;i<middle_point;++i) {
      SWAP(cigar_buffer_base[i],cigar_buffer_base[num_cigar_elements-i-1]);
    }
  }
  // Update CIGAR/Score values
  *cigar_length += num_cigar_elements; // Update CIGAR length
  *effective_length += match_effective_length; // Update effective length
  *alignment_score += match_alignment_score; // Update alignment-score
}
GEM_INLINE void swg_align_match_full_32b(
    const uint8_t* const key,const uint64_t key_length,const swg_penalties_t* swg_penalties,
    uint64_t* const match_position,uint8_t* const text,const uint64_t text_length,
    vector_t* const cigar_buffer,uint64_t* const cigar_length,int64_t* const effective_length,
    int32_t* const alignment_score,mm_stack_t* const mm_stack) {
  /*
   * Initialize
   */
  uint64_t column, row;
  mm_stack_push_state(mm_stack); // Save stack state
  const uint64_t row_size = (key_length+1)*sizeof(swg_cell_t);
  const uint64_t num_rows = (key_length+1);
  const uint64_t num_columns = (text_length+1);
  swg_cell_t** const dp = mm_stack_malloc(mm_stack,num_columns*sizeof(swg_cell_t*));
  for (column=0;column<num_columns;++column) {
    dp[column] = mm_stack_malloc(mm_stack,row_size);
  }
  // Initialize first row
  const matching_score_t* const matching_score = &swg_penalties->matching_score;
  const int32_t single_gap = - (swg_penalties->gap_open_penalty + swg_penalties->gap_extension_penalty); // g(1)
  const int32_t gap_extension = - (swg_penalties->gap_extension_penalty);
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
  /*
   * Compute DP-matrix
   */
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
  // Retrieve the alignment. Store the match (Backtrace and generate CIGAR)
  swg_align_match_traceback_32b(dp,max_score,
      max_score_column,single_gap,gap_extension,key,key_length,
      match_position,text,true,cigar_buffer,cigar_length,effective_length,alignment_score);
  // Clean-up
  mm_stack_pop_state(mm_stack,false); // Free
}
GEM_INLINE void swg_align_match_banded_32b(
    const uint8_t* const key,const uint64_t key_length,const swg_penalties_t* swg_penalties,
    uint64_t* const match_position,uint8_t* const text,const uint64_t text_length,
    uint64_t max_bandwidth,const bool begin_free,const bool end_free,
    vector_t* const cigar_buffer,uint64_t* const cigar_length,int64_t* const effective_length,
    int32_t* const alignment_score,mm_stack_t* const mm_stack) {
  // Initialize band-limits
  if (max_bandwidth > key_length) max_bandwidth = key_length;
  int64_t top_bandwidth = text_length - key_length + max_bandwidth + 1;
  if (top_bandwidth < 0) {
    *alignment_score = ALIGN_DISTANCE_INF; // Impossible alignment
    return;
  }
  if (top_bandwidth > text_length) top_bandwidth = text_length;
  // Allocate memory
  uint64_t column, row;
  const uint64_t row_size = (key_length+1)*sizeof(swg_cell_t);
  const uint64_t num_rows = (key_length+1);
  const uint64_t num_columns = (text_length+1);
  mm_stack_push_state(mm_stack); // Save stack state
  swg_cell_t** const dp = mm_stack_malloc(mm_stack,num_columns*sizeof(swg_cell_t*));
  for (column=0;column<num_columns;++column) {
    dp[column] = mm_stack_malloc(mm_stack,row_size);
  }
  // Initialize DP-matrix
  const matching_score_t* const matching_score = &swg_penalties->matching_score;
  const int32_t single_gap = - (swg_penalties->gap_open_penalty + swg_penalties->gap_extension_penalty); // g(1)
  const int32_t gap_extension = - (swg_penalties->gap_extension_penalty);
  // Initialize first column
  const uint64_t num_rows_1 = num_rows-1;
  uint64_t band_high_offset = 0;
  uint64_t band_low_offset = (max_bandwidth+1 < num_rows) ? max_bandwidth+1 : num_rows_1;
  dp[0][0].D = SWG_SCORE_INT32_MIN;
  dp[0][0].I = SWG_SCORE_INT32_MIN;
  dp[0][0].M = 0;
  dp[0][1].M = single_gap; // g(1)
  for (row=2;row<band_low_offset;++row) {
    dp[0][row].I = SWG_SCORE_INT32_MIN;
    dp[0][row].M = dp[0][row-1].M + gap_extension; // g(row)
  }
  // Initialize first row
  dp[0][1].I = SWG_SCORE_INT32_MIN;
  if (begin_free) {
    for (column=1;column<top_bandwidth;++column) {
      dp[column][0].M = 0;
    }
  } else {
    dp[1][0].M = single_gap;
    for (column=2;column<top_bandwidth;++column) {
      dp[column][0].M = dp[column-1][0].M + gap_extension;
    }
  }
  /*
   * Compute DP-matrix
   */
  int32_t max_score = SWG_SCORE_INT32_MIN;
  uint64_t max_score_column = UINT64_MAX;
  for (column=1;column<num_columns;++column) {
    // Initialize band boundaries
    if (column >= top_bandwidth) {
      dp[column][band_high_offset].M = SWG_SCORE_INT32_MIN;
    }
    dp[column][band_high_offset].D = SWG_SCORE_INT32_MIN;
    // Locate the cursor at the proper cell & calculate DP
    for (row=band_high_offset+1;row<band_low_offset;++row) {
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
    if (end_free && band_low_offset==num_rows) {
      if (dp[column][num_rows_1].M > max_score) {
        max_score = dp[column][num_rows_1].M;
        max_score_column = column;
      }
    }
    // Update band limits
    if (column >= top_bandwidth) ++band_high_offset; // Swift band
    if (row < num_rows) { //  Below band
      dp[column][row].I = SWG_SCORE_INT32_MIN;
      dp[column][row].M = SWG_SCORE_INT32_MIN;
      ++band_low_offset;
    }
  }
  // Set alignment score/column
  if (!end_free) {
    max_score = dp[num_columns-1][num_rows_1].M;
    max_score_column = num_columns-1;
  }
  // Retrieve the alignment. Store the match (Backtrace and generate CIGAR)
  swg_align_match_traceback_32b(dp,max_score,max_score_column,
      single_gap,gap_extension,key,key_length,match_position,
      text,begin_free,cigar_buffer,cigar_length,effective_length,alignment_score);
  // Clean-up
  mm_stack_pop_state(mm_stack,false); // Free
}
GEM_INLINE void swg_align_match(
    const uint8_t* const key,const uint64_t key_length,const bool* const allowed_enc,
    const swg_penalties_t* swg_penalties,uint64_t* const match_position,uint8_t* const text,
    const uint64_t text_length,uint64_t max_bandwidth,const bool begin_free,const bool end_free,
    vector_t* const cigar_buffer,uint64_t* const cigar_length,int64_t* const effective_length,
    int32_t* const alignment_score,mm_stack_t* const mm_stack) {
  // fprintf(stderr,"%lu band=%lu (%lu x %lu)\n",key_length*text_length,max_bandwidth,key_length,text_length);
  // Check lengths
  if (key_length == 0 && text_length > 0) {
    if (begin_free) {
      // Adjust position
      *match_position += text_length;
    } else if (!end_free) {
      // Insertion <+@text_length>
      matches_cigar_buffer_append_indel(cigar_buffer,cigar_length,cigar_ins,text_length,text);
      *alignment_score += -swg_penalties->gap_open_penalty - swg_penalties->gap_extension_penalty*text_length;
      *effective_length += text_length;
    }
  } else if (text_length == 0) {
    if (key_length > 0) {
      // Deletion <-@key_length>
      matches_cigar_buffer_append_indel(cigar_buffer,cigar_length,cigar_del,key_length,NULL);
      *alignment_score += -swg_penalties->gap_open_penalty - swg_penalties->gap_extension_penalty*key_length;
    }
  } else {
    // Check trivial cases
    if (key_length==1 && key_length==text_length) {
      // Mismatch/Match
      const uint8_t key_enc = key[0];
      const uint8_t text_enc = text[0];
      if (!allowed_enc[text_enc] || text_enc != key_enc) {
        matches_cigar_buffer_append_mismatch(cigar_buffer,cigar_length,cigar_mismatch,text_enc);
        *alignment_score += swg_penalties->matching_score[text_enc][key_enc];
      } else {
        matches_cigar_buffer_append_match(cigar_buffer,cigar_length,1);
        *alignment_score += swg_penalties->matching_score[text_enc][key_enc];
      }
      *effective_length += 1;
    } else {
      // if (max_reachable_score <= MAX_UINT8) {} // TODO
      // if (band is useless) unbanded // TODO
      swg_align_match_banded_32b(
          key,key_length,swg_penalties,match_position,text,text_length,max_bandwidth,begin_free,end_free,
          cigar_buffer,cigar_length,effective_length,alignment_score,mm_stack);
    }
  }
}
GEM_INLINE void swg_align_match_banded_32b_opt(
    const uint8_t* const key,const uint64_t key_length,const swg_penalties_t* swg_penalties,
    uint64_t* const match_position,uint8_t* const text,const uint64_t text_length,
    uint64_t max_bandwidth,const bool begin_free,const bool end_free,
    vector_t* const cigar_buffer,uint64_t* const cigar_length,int64_t* const effective_length,
    int32_t* const alignment_score,mm_stack_t* const mm_stack) {
  // Initialize band-limits
  if (max_bandwidth > key_length) max_bandwidth = key_length;
  int64_t top_bandwidth = text_length - key_length + max_bandwidth + 1;
  if (top_bandwidth < 0) {
    *alignment_score = ALIGN_DISTANCE_INF; // Impossible alignment
    return;
  }
  if (top_bandwidth > text_length) top_bandwidth = text_length;
  // Allocate memory
  uint64_t column, row;
  const uint64_t row_size = (key_length+2)*sizeof(swg_cell_t); // (+1) Init; (+1) next D-Value
  const uint64_t num_rows = (key_length+1);                    // (+1) Init
  const uint64_t num_columns = (text_length+1);
  mm_stack_push_state(mm_stack); // Save stack state
  swg_cell_t** const dp = mm_stack_malloc(mm_stack,(num_columns+1)*sizeof(swg_cell_t*));
  for (column=0;column<=num_columns;++column) { // (+1) next I-Value
    dp[column] = mm_stack_malloc(mm_stack,row_size);
  }
  // Initialize DP-matrix
  const matching_score_t* const matching_score = &swg_penalties->matching_score;
  const int32_t single_gap = - (swg_penalties->gap_open_penalty + swg_penalties->gap_extension_penalty); // g(1)
  const int32_t gap_extension = - (swg_penalties->gap_extension_penalty);
  // Initialize first column
  const uint64_t num_rows_1 = num_rows-1;
  uint64_t band_high_offset = 0;
  uint64_t band_low_offset = (max_bandwidth+1 < num_rows) ? max_bandwidth+1 : num_rows_1;
  dp[0][0].D = SWG_SCORE_INT32_MIN;
  dp[0][0].I = SWG_SCORE_INT32_MIN;
  dp[0][0].M = 0;
  dp[0][1].M = single_gap; // g(1)
  for (row=2;row<band_low_offset;++row) {
    dp[0][row].I = SWG_SCORE_INT32_MIN;
    dp[0][row].M = dp[0][row-1].M + gap_extension; // g(row)
  }
  // Initialize first row
  dp[0][1].I = SWG_SCORE_INT32_MIN;
  if (begin_free) {
    for (column=1;column<top_bandwidth;++column) dp[column][0].M = 0;
  } else {
    dp[1][0].M = single_gap;
    for (column=2;column<top_bandwidth;++column) dp[column][0].M = dp[column-1][0].M + gap_extension;
  }
  /*
   * Compute DP-matrix
   */
  int32_t max_score = SWG_SCORE_INT32_MIN;
  uint64_t max_score_column = UINT64_MAX;
  for (column=1;column<num_columns;++column) {
    // Initialize band boundaries
    if (column >= top_bandwidth) dp[column][band_high_offset].M = SWG_SCORE_INT32_MIN;
    dp[column][band_high_offset].D = SWG_SCORE_INT32_MIN;
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
    // Update band limits
    if (column >= top_bandwidth) ++band_high_offset; // Swift band
    if (row < num_rows) { //  Below band
      dp[column][row].I = SWG_SCORE_INT32_MIN;
      dp[column][row].M = SWG_SCORE_INT32_MIN;
      ++band_low_offset;
    }
  }
  // Set alignment score/column
  if (!end_free) {
    max_score = dp[num_columns-1][num_rows_1].M;
    max_score_column = num_columns-1;
  }
  // Retrieve the alignment. Store the match (Backtrace and generate CIGAR)
  swg_align_match_traceback_32b(dp,max_score,max_score_column,
      single_gap,gap_extension,key,key_length,match_position,
      text,begin_free,cigar_buffer,cigar_length,effective_length,alignment_score);
  // Clean-up
  mm_stack_pop_state(mm_stack,false); // Free
}
