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
 *   Efficient Smith-Waterman-Gotoh (SWG) alignment module using SIMD instructions
 */

#include "align/align_swg_simd.h"

/*
 * Init SWG Query Profile
 */
void align_swg_query_profile_allocate_uint8(
    swg_query_profile_t* const swg_query_profile,
    const uint64_t max_expected_key_length,
    mm_allocator_t* const mm_allocator) {
  // Compute sizes 8-bits cell
  const uint64_t num_segments_uint8 = UINT128_SIZE/UINT8_SIZE;
  const uint64_t segment_length_uint8 = DIV_CEIL(max_expected_key_length,num_segments_uint8);
  const uint64_t key_effective_length_uint8 = segment_length_uint8*num_segments_uint8;
  // Allocate 8-bits cell
  uint8_t enc;
  for (enc=0;enc<DNA__N_RANGE;++enc) {
    swg_query_profile->query_profile_uint8[enc] =
        mm_allocator_calloc(mm_allocator,key_effective_length_uint8,uint8_t,true);
  }
}
void align_swg_query_profile_allocate_int16(
    swg_query_profile_t* const swg_query_profile,
    const uint64_t max_expected_key_length,
    mm_allocator_t* const mm_allocator) {
  // Compute sizes 16-bits cell
  const uint64_t num_segments_int16 = UINT128_SIZE/UINT16_SIZE;
  const uint64_t segment_length_int16 = DIV_CEIL(max_expected_key_length,num_segments_int16);
  const uint64_t key_effective_length_uint16 = segment_length_int16*num_segments_int16;
  // Allocate 16-bits cell
  uint8_t enc;
  for (enc=0;enc<DNA__N_RANGE;++enc) {
    swg_query_profile->query_profile_int16[enc] =
        (int16_t*)mm_allocator_calloc(mm_allocator,key_effective_length_uint16,uint8_t,true);
  }
}
void align_swg_query_profile_init(
    swg_query_profile_t* const swg_query_profile,
    const swg_penalties_t* swg_penalties,
    const uint64_t max_expected_key_length,
    mm_allocator_t* const mm_allocator) {
  // Penalties
  const int32_t generic_match_score = swg_penalties->generic_match_score;
  const int32_t generic_mismatch_score = swg_penalties->generic_mismatch_score;
  // Compute the max/min scores
  const int32_t min_profile_score = MIN(generic_match_score,generic_mismatch_score);
  const int32_t match_bias_uint8 = -min_profile_score;
  swg_query_profile->match_bias_uint8 = match_bias_uint8;
  // Allocate
  align_swg_query_profile_allocate_uint8(swg_query_profile,max_expected_key_length,mm_allocator);
  align_swg_query_profile_allocate_int16(swg_query_profile,max_expected_key_length,mm_allocator);
}
/*
 * Compile SWG Query Profile
 */
bool align_swg_query_profile_compile_uint8(
    swg_query_profile_t* const swg_query_profile,
    const swg_penalties_t* swg_penalties,
    const uint8_t* const key,
    const uint64_t key_length) {
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
bool align_swg_query_profile_compile_int16(
    swg_query_profile_t* const swg_query_profile,
    const swg_penalties_t* swg_penalties,
    const uint8_t* const key,
    const uint64_t key_length) {
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

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

//void swg_align_match_table_print_int16(
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
//int16_t** swg_align_match_allocate_table_int16(
//    const uint64_t num_columns,const uint64_t num_rows,mm_allocator_t* const mm_allocator) {
//  // Allocate the pointers
//  int16_t** const dp = mm_allocator_malloc(mm_allocator,num_columns*sizeof(int16_t*));
//  // Allocate all the columns
//  const uint64_t row_size = num_rows*sizeof(int16_t);
//  uint64_t column;
//  for (column=0;column<num_columns;++column) {
//    dp[column] = mm_allocator_malloc(mm_allocator,row_size);
//  }
//  return dp;
//}
//void swg_align_match_init_table_int16(
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

//void swg_align_match_int16_simd128(
//    const uint8_t* const key,const uint64_t key_length,swg_query_profile_t* const swg_query_profile,
//    const swg_penalties_t* swg_penalties,uint64_t* const match_position,uint8_t* const text,
//    const uint64_t text_length,const bool begin_free,const bool end_free,vector_t* const cigar_vector,
//    uint64_t* const cigar_length,int64_t* const effective_length,int32_t* const alignment_score,
//    mm_allocator_t* const mm_allocator) {
//  // Compile query profile
//  swg_compile_query_profile_int16(swg_query_profile,swg_penalties,key,key_length,mm_allocator);
//  // Allocate memory
//  mm_allocator_push_state(mm_allocator); // Save allocator state
//  const uint64_t num_segments_16b = UINT128_SIZE/UINT16_SIZE; // Elements in each SIMD-vectors
//  const uint64_t segment_length = swg_query_profile->segment_length_int16; // Number of SIMD-vectors
//  const uint64_t key_effective_length = swg_query_profile->key_effective_length_int16; // Number of SIMD-vectors
//  const uint64_t num_columns = (text_length+1);
//  const uint64_t num_rows = key_effective_length;
//  int16_t** const dp_M = swg_align_match_allocate_table_int16(num_columns,num_rows,mm_allocator);
//  int16_t** const dp_I = swg_align_match_allocate_table_int16(num_columns+1,num_rows,mm_allocator);
//  int16_t** const dp_D = swg_align_match_allocate_table_int16(num_columns,num_rows,mm_allocator);
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
//  mm_allocator_pop_state(mm_allocator,false); // Free
//  *alignment_score = max_score_column; // FIXME
//  *alignment_score = max_score_global; // FIXME
//}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

//void swg_align_match_table_print_uint8(
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
//uint8_t** swg_align_match_allocate_table_uint8(
//    const uint64_t num_columns,const uint64_t num_rows,mm_allocator_t* const mm_allocator) {
//  // Allocate the pointers
//  uint8_t** const dp = mm_allocator_malloc(mm_allocator,num_columns*sizeof(uint8_t*));
//  // Allocate all the columns
//  mm_allocator_skip_align(mm_allocator,UINT128_SIZE); // Align allocator Memory
//  const uint64_t row_size = num_rows*sizeof(uint8_t);
//  uint64_t column;
//  for (column=0;column<num_columns;++column) {
//    dp[column] = mm_allocator_malloc(mm_allocator,row_size);
//  }
//  return dp;
//}
//void swg_align_match_init_table_uint8(
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
//bool swg_align_match_uint8_simd128(
//    const uint8_t* const key,const uint64_t key_length,swg_query_profile_t* const swg_query_profile,
//    const swg_penalties_t* swg_penalties,uint64_t* const match_position,uint8_t* const text,
//    const uint64_t text_length,const bool begin_free,const bool end_free,vector_t* const cigar_vector,
//    uint64_t* const cigar_length,int64_t* const effective_length,int32_t* const alignment_score,
//    mm_allocator_t* const mm_allocator) {
//  // Compile query profile
//  if (!swg_compile_query_profile_uint8(swg_query_profile,swg_penalties,key,key_length,mm_allocator)) {
//    return false; // TODO
//  }
//  // Allocate memory
//  mm_allocator_push_state(mm_allocator); // Save allocator state
//  const uint64_t num_segments_8b = UINT128_SIZE/UINT8_SIZE; // Elements in each SIMD-vectors
//  const uint64_t segment_length = swg_query_profile->segment_length_uint8; // Number of SIMD-vectors
//  const uint64_t key_effective_length = swg_query_profile->key_effective_length_uint8; // Number of SIMD-vectors
//  const uint64_t num_columns = (text_length+1);
//  const uint64_t num_rows = key_effective_length;
//  uint8_t** const dp_M = swg_align_match_allocate_table_uint8(num_columns,num_rows,mm_allocator); // TODO key_length
//  uint8_t** const dp_I = swg_align_match_allocate_table_uint8(num_columns+1,num_rows,mm_allocator);
//  uint8_t** const dp_D = swg_align_match_allocate_table_uint8(num_columns,num_rows,mm_allocator);
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
//  mm_allocator_pop_state(mm_allocator,false); // Free
//  *alignment_score = max_score_column; // FIXME
//  *alignment_score = max_score_global; // FIXME
//  return ((int)max_score_global + (int)match_bias < 255);
//}
//
//
