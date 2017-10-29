/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_align_bpm128.c
 * DATE: 20/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Myers' Bit-vector algorithm.
 *   Implementation using SIMD instructions (vectorial extensions)
 */

#include "gt_map_align_bpm_simd.h"

#include <immintrin.h> // ALL
//#include <mmintrin.h>    // MMX
//#include <xmmintrin.h>   // SSE
//#include <emmintrin.h>   // SSE2
//#include <pmmintrin.h>   // SSE3
//#include <tmmintrin.h>   // SSSE3
//#include <smmintrin.h>   // SSE4

//#ifndef _mm_test_all_zeros
//#define _mm_test_all_zeros(M,V) _mm_testz_si128((M),(V))
//#endif

#define GT_BMP_W64_MASK 1L<<63

/*
 * Bit-compressed (Re)alignment
 *   BMP[BitParalellMayers] - Myers' Fast Bit-Vector algorithm (Levenshtein)
 */
GT_INLINE bool gt_map_block_bpm_get_distance_simd128_unfold64(
    gt_bpm_pattern* const bpm_pattern,char* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance) {
  // Pattern variables
  const uint64_t* peq = (uint64_t*)bpm_pattern->peq;
  const uint64_t num_words = bpm_pattern->pattern_num_words;
  uint64_t* const P = (uint64_t*)bpm_pattern->P;
  uint64_t* const M = (uint64_t*)bpm_pattern->M;
  const uint64_t* const level_mask = (uint64_t*)bpm_pattern->level_mask;
  int64_t* const score = (int64_t*)bpm_pattern->score;
  const int64_t* const init_score = (int64_t*)bpm_pattern->init_score;

  // Initialize search
  uint64_t min_score = UINT64_MAX, min_score_position = UINT64_MAX;
  uint64_t i;
  P[0]=UINT64_MAX;
  P[1]=UINT64_MAX;
  M[0]=0;
  M[1]=0;
  score[0] = init_score[0];
  for (i=1;i<num_words;++i) {
    const uint64_t idx = i<<1;
    P[idx]=UINT64_MAX;
    P[idx+1]=UINT64_MAX;
    M[idx]=0;
    M[idx+1]=0;
    score[i] = score[i-1] + init_score[i];
  }

  // Advance in DP-bit_encoded matrix
  uint64_t sequence_position;
  for (sequence_position=0;sequence_position<sequence_length;++sequence_position) {
    // Fetch next character
    const uint8_t enc_char = gt_cdna_encode(sequence[sequence_position]);

    // Advance all blocks
    uint64_t PHin=0,MHin=0,PHout,MHout;
    for (i=0;i<num_words;++i) {
      uint64_t* const P_mem_addr = P+(i<<1);
      uint64_t* const M_mem_addr = M+(i<<1);
      const uint64_t* const mask_mem_addr = (level_mask+(i<<1));
      const uint64_t* const Eq_mem_addr = peq+((enc_char*num_words+i)<<1);

      /*
       * Compute Block
       */
      // Xv = Eq | Mv;
      uint64_t Xv[2];
      Xv[0] = Eq_mem_addr[0] | M_mem_addr[0];
      Xv[1] = Eq_mem_addr[1] | M_mem_addr[1];
      // _Eq = Eq | MHin;
      uint64_t _Eq[2];
      _Eq[0] = Eq_mem_addr[0] | MHin;
      _Eq[1] = Eq_mem_addr[1];
      // Xh = (((_Eq & Pv) + Pv) ^ Pv) | _Eq;
      uint64_t Xh[2];
      Xh[0] = (_Eq[0] & P_mem_addr[0]);
      Xh[1] = (_Eq[1] & P_mem_addr[1]);
      Xh[0] = Xh[0] + P_mem_addr[0];
      Xh[1] = Xh[1] + P_mem_addr[1] + (Xh[0]<P_mem_addr[0]);
      Xh[0] = (Xh[0] ^ P_mem_addr[0]) | _Eq[0];
      Xh[1] = (Xh[1] ^ P_mem_addr[1]) | _Eq[1];

      /* Calculate Hout */
      // Ph = Mv | ~(Xh | Pv);
      uint64_t Ph[2];
      Ph[0] = M_mem_addr[0] | ~(Xh[0] | P_mem_addr[0]);
      Ph[1] = M_mem_addr[1] | ~(Xh[1] | P_mem_addr[1]);
      // Mh = Pv & Xh;
      uint64_t Mh[2];
      Mh[0] = P_mem_addr[0] & Xh[0];
      Mh[1] = P_mem_addr[1] & Xh[1];

      /* Account Hout that propagates for the next block */
      // PHout = (Ph & mask)!=0;
      // MHout = (Mh & mask)!=0;
      PHout = ((Ph[0] & mask_mem_addr[0]) | (Ph[1] & mask_mem_addr[1])) !=0;
      MHout = ((Mh[0] & mask_mem_addr[0]) | (Mh[1] & mask_mem_addr[1])) !=0;

      /* Hout become the Hin of the next cell */
      // Ph <<= 1;
      // Mh <<= 1;
      /* Account Hin coming from the previous block */
      // Ph |= PHin;
      // Mh |= MHin;
      Ph[1] = (Ph[1]<<1) | ((Ph[0] & GT_BMP_W64_MASK)!=0);
      Ph[0] = (Ph[0]<<1) | PHin;
      Mh[1] = (Mh[1]<<1) | ((Mh[0] & GT_BMP_W64_MASK)!=0);
      Mh[0] = (Mh[0]<<1) | MHin;

      /* Finally, generate the Vout */
      // Pv = Mh | ~(Xv | Ph);
      P_mem_addr[0] = Mh[0] | ~(Xv[0] | Ph[0]);
      P_mem_addr[1] = Mh[1] | ~(Xv[1] | Ph[1]);
      // Mv = Ph & Xv
      M_mem_addr[0] = Ph[0] & Xv[0];
      M_mem_addr[1] = Ph[1] & Xv[1];

      /* Adjust score and swap propagate Hv */
      score[i] += PHout-MHout;
      PHin=PHout;
      MHin=MHout;
    }

    // Check match
    if (score[num_words-1] < min_score) {
      min_score_position = sequence_position;
      min_score = score[num_words-1];
    }
  }
  // Return results
  if (min_score!=UINT64_MAX) {
    *distance = min_score;
    *position = min_score_position;
    return true;
  } else {
    *distance = UINT64_MAX;
    *position = UINT64_MAX;
    return false;
  }
}
GT_INLINE bool gt_map_block_bpm_get_distance_simd128(
    gt_bpm_pattern* const bpm_pattern,char* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance) {
  // Pattern variables
  const uint64_t* peq = (uint64_t*)bpm_pattern->peq;
  const uint64_t num_words = bpm_pattern->pattern_num_words;
  uint64_t* const P = (uint64_t*)bpm_pattern->P;
  uint64_t* const M = (uint64_t*)bpm_pattern->M;
  const uint64_t* const level_mask = (uint64_t*)bpm_pattern->level_mask;
  int64_t* const score = (int64_t*)bpm_pattern->score;
  const int64_t* const init_score = (int64_t*)bpm_pattern->init_score;

  // Initialize search
  uint64_t min_score = UINT64_MAX, min_score_position = UINT64_MAX;
  uint64_t i;
  P[0]=UINT64_MAX;
  P[1]=UINT64_MAX;
  M[0]=0;
  M[1]=0;
  score[0] = init_score[0];
  for (i=1;i<num_words;++i) {
    const uint64_t idx = i<<1;
    P[idx]=UINT64_MAX;
    P[idx+1]=UINT64_MAX;
    M[idx]=0;
    M[idx+1]=0;
    score[i] = score[i-1] + init_score[i];
  }

  // Advance in DP-bit_encoded matrix
  uint64_t sequence_position;
  for (sequence_position=0;sequence_position<sequence_length;++sequence_position) {
    // Fetch next character
    const uint8_t enc_char = gt_cdna_encode(sequence[sequence_position]);

    // Advance all blocks
    uint64_t PHin=0,MHin=0,PHout,MHout;
    for (i=0;i<num_words;++i) {
      uint64_t* const P_mem_addr = P+(i<<1);
      uint64_t* const M_mem_addr = M+(i<<1);
      const uint64_t* const mask_mem_addr = (level_mask+(i<<1));
      const uint64_t* const Eq_mem_addr = peq+((enc_char*num_words+i)<<1);
      // __m128i Pv = _mm_load_si128((__m128i*)P_mem_addr);
      // __m128i Mv = _mm_load_si128((__m128i*)M_mem_addr);
      // const __m128i mask = _mm_load_si128((__m128i*)mask_mem_addr);
      // const __m128i Eq = _mm_load_si128((__m128i*)Eq_mem_addr);

      /*
       * Compute Block
       */
      // Xv = Eq | Mv;
      uint64_t Xv[2];
      __m128i SIMD_Mv = _mm_load_si128((__m128i*)M_mem_addr);
      const __m128i SIMD_Eq = _mm_load_si128((__m128i*)Eq_mem_addr);
      const __m128i SIMD_Xv = _mm_or_si128(SIMD_Eq,SIMD_Mv);
      _mm_stream_si128((__m128i*)Xv,SIMD_Xv);

      // _Eq = Eq | MHin;
      const __m128i SIMD__Eq = _mm_or_si128(SIMD_Eq,_mm_set_epi64x(0,MHin));

      // Xh = (((_Eq & Pv) + Pv) ^ Pv) | _Eq;
      uint64_t Xh[2];
      __m128i SIMD_Pv = _mm_load_si128((__m128i*)P_mem_addr);
      _mm_stream_si128((__m128i*)Xh,_mm_and_si128(SIMD__Eq,SIMD_Pv)); // (_Eq & Pv)
      Xh[0] = Xh[0] + P_mem_addr[0];
      Xh[1] = Xh[1] + P_mem_addr[1] + (Xh[0]<P_mem_addr[0]);
      const __m128i SIMD_Xh = _mm_or_si128(_mm_xor_si128(_mm_load_si128((__m128i*)Xh),SIMD_Pv),SIMD__Eq);
      _mm_stream_si128((__m128i*)Xh,SIMD_Xh);

      /* Calculate Hout */
      // Ph = Mv | ~(Xh | Pv);
      __m128i SIMD_Ph = _mm_or_si128(SIMD_Mv,
          _mm_xor_si128(_mm_or_si128(SIMD_Xh,SIMD_Pv),_mm_set1_epi64((__m64)UINT64_MAX)) );
      // Mh = Pv & Xh;
      __m128i SIMD_Mh = _mm_and_si128(SIMD_Pv,SIMD_Xh);

      uint64_t Ph[2],Mh[2];
      _mm_stream_si128((__m128i*)Ph,SIMD_Ph);
      _mm_stream_si128((__m128i*)Mh,SIMD_Mh);

      /* Account Hout that propagates for the next block */
      __m128i SIMD_mask = _mm_load_si128((__m128i*)mask_mem_addr);
      PHout = (_mm_test_all_zeros(SIMD_Ph,SIMD_mask) == 0); // PHout = (Ph & mask)!=0;
      MHout = (_mm_test_all_zeros(SIMD_Mh,SIMD_mask) == 0); // MHout = (Mh & mask)!=0;

      /* Hout become the Hin of the next cell */
      // Ph <<= 1;
      // Mh <<= 1;
      /* Account Hin coming from the previous block */
      // Ph |= PHin;
      // Mh |= MHin;
      Ph[1] = (Ph[1]<<1) | ((Ph[0] & GT_BMP_W64_MASK)!=0);
      Ph[0] = (Ph[0]<<1) | PHin;
      Mh[1] = (Mh[1]<<1) | ((Mh[0] & GT_BMP_W64_MASK)!=0);
      Mh[0] = (Mh[0]<<1) | MHin;

      /* Finally, generate the Vout */
      SIMD_Ph = _mm_load_si128((__m128i*)Ph);
      SIMD_Mh = _mm_load_si128((__m128i*)Mh);
      // Pv = Mh | ~(Xv | Ph);
      SIMD_Pv = _mm_or_si128(SIMD_Mh,_mm_xor_si128(_mm_or_si128(SIMD_Xv,SIMD_Ph),_mm_set1_epi64((__m64)UINT64_MAX)));
      _mm_stream_si128((__m128i*)P_mem_addr,SIMD_Pv);
      // Mv = Ph & Xv
      SIMD_Mv = _mm_and_si128(SIMD_Ph,SIMD_Xv);
      _mm_stream_si128((__m128i*)M_mem_addr,SIMD_Mv);

      /* Adjust score and swap propagate Hv */
      score[i] += PHout-MHout;
      PHin=PHout;
      MHin=MHout;
    }

    // Check match
    if (score[num_words-1] < min_score) {
      min_score_position = sequence_position;
      min_score = score[num_words-1];
    }
  }
  // Return results
  if (min_score!=UINT64_MAX) {
    *distance = min_score;
    *position = min_score_position;
    return true;
  } else {
    *distance = UINT64_MAX;
    *position = UINT64_MAX;
    return false;
  }
}
