/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_align_bmp_generic.c
 * DATE: 20/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Myers' Bit-vector algorithm.
 *   Optimized version from the base implementation.
 *   Generic interface as to use any length word size
 */

#include "gt_map_align_bpm_generic.h"

/*
 * NEEDED DEFINITIONS
 *  - #define gt_bmp_vector_t uint8_t
 */
#ifndef gt_bmp_vector_t
  #error "GT.Generic.Code:: Type gt_bmp_vector_t should be defined"
#else

/*
 * Advance block functions
 *   const @vector Eq,mask;
 *   return (Pv,Mv,PHout,MHout);
 */
#define GT_MAP_BLOCK_BPM_ADVANCE_BLOCK(Eq,mask,Pv,Mv,PHin,MHin,PHout,MHout) \
  /* Computes modulator vector {Xv,Xh} ( cases A&C ) */ \
  const gt_bmp_vector_t Xv = Eq | Mv; \
  const gt_bmp_vector_t _Eq = Eq | MHin; \
  const gt_bmp_vector_t Xh = (((_Eq & Pv) + Pv) ^ Pv) | _Eq; \
  /* Calculate Hout */ \
  gt_bmp_vector_t Ph = Mv | ~(Xh | Pv); \
  gt_bmp_vector_t Mh = Pv & Xh; \
  /* Account Hout that propagates for the next block */ \
  PHout = (Ph & mask)!=0; \
  MHout = (Mh & mask)!=0; \
  /* Hout become the Hin of the next cell */ \
  Ph <<= 1; \
  Mh <<= 1; \
  /* Account Hin coming from the previous block */ \
  Ph |= PHin; \
  Mh |= MHin; \
  /* Finally, generate the Vout */ \
  Pv = Mh | ~(Xv | Ph); \
  Mv = Ph & Xv
/*
 * Bit-compressed (Re)alignment
 *   BMP[BitParalellMayers] - Myers' Fast Bit-Vector algorithm (Levenshtein)
 */
GT_INLINE bool GT_BMP_GENERIC_FUNCTION_NAME(gt_map_block_bpm_get_distance)(
    gt_bpm_pattern* const bpm_pattern,char* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance) {
  // Pattern variables
  const gt_bmp_vector_t* peq = (gt_bmp_vector_t*)bpm_pattern->peq;
  const uint64_t num_words = bpm_pattern->pattern_num_words;
  gt_bmp_vector_t* const P = (gt_bmp_vector_t*)bpm_pattern->P;
  gt_bmp_vector_t* const M = (gt_bmp_vector_t*)bpm_pattern->M;
  const gt_bmp_vector_t* const level_mask = (gt_bmp_vector_t*)bpm_pattern->level_mask;
  int64_t* const score = (int64_t*)bpm_pattern->score;
  const int64_t* const init_score = (int64_t*)bpm_pattern->init_score;

  // Initialize search
  const gt_bmp_vector_t word_ones = ((gt_bmp_vector_t)0) - ((gt_bmp_vector_t)1);
  uint64_t min_score = UINT64_MAX, min_score_position = UINT64_MAX;
  uint64_t i;
  P[0]=word_ones;
  M[0]=0;
  score[0] = init_score[0];
  for (i=1;i<num_words;++i) {
    P[i]=word_ones;
    M[i]=0;
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
      gt_bmp_vector_t Pv = P[i];
      gt_bmp_vector_t Mv = M[i];
      const gt_bmp_vector_t mask = level_mask[i];
      const gt_bmp_vector_t Eq = peq[GT_BPM_PEQ_IDX(enc_char,i,num_words)];
      /* Compute Block */
      GT_MAP_BLOCK_BPM_ADVANCE_BLOCK(Eq,mask,Pv,Mv,PHin,MHin,PHout,MHout);
      /* Save Block Pv,Mv */
      P[i]=Pv;
      M[i]=Mv;
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
GT_INLINE bool GT_BMP_GENERIC_FUNCTION_NAME(gt_map_block_bpm_get_distance__cutoff)(
    gt_bpm_pattern* const bpm_pattern,char* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance,const uint64_t max_distance) {
  // Pattern variables
  const uint64_t num_words = bpm_pattern->pattern_num_words;
  const uint64_t max_distance__1 = max_distance+1;
  const gt_bmp_vector_t* peq = (gt_bmp_vector_t*)bpm_pattern->peq;
  gt_bmp_vector_t* const P = (gt_bmp_vector_t*)bpm_pattern->P;
  gt_bmp_vector_t* const M = (gt_bmp_vector_t*)bpm_pattern->M;
  const gt_bmp_vector_t* const level_mask = (gt_bmp_vector_t*)bpm_pattern->level_mask;
  int64_t* const score = (int64_t*)bpm_pattern->score;
  const int64_t* const init_score = (int64_t*)bpm_pattern->init_score;

  // Initialize search
  const gt_bmp_vector_t word_ones = ((gt_bmp_vector_t)0) - ((gt_bmp_vector_t)1);
  const gt_bmp_vector_t word_length = sizeof(gt_bmp_vector_t)*8;
  const uint64_t top = num_words-1;
  uint64_t top_level = (max_distance>0) ? (max_distance+(word_length-1))/word_length : 1;
  uint64_t min_score = UINT64_MAX, min_score_position = UINT64_MAX;
  // Reset score,P,M
  uint64_t i;
  P[0]=word_ones;
  M[0]=0;
  score[0] = init_score[0];
  for (i=1;i<top_level;++i) {
    P[i]=word_ones;
    M[i]=0;
    score[i] = score[i-1] + init_score[i];
  }

  // Advance in DP-bit_encoded matrix
  uint64_t sequence_position;
  for (sequence_position=0;sequence_position<sequence_length;++sequence_position) {
    // Fetch next character
    const uint8_t enc_char = gt_cdna_encode(sequence[sequence_position]);

    // Advance all blocks
    uint64_t i;
    uint64_t PHin=0,MHin=0,PHout,MHout;
    for (i=0;i<top_level;++i) {
      gt_bmp_vector_t Pv = P[i];
      gt_bmp_vector_t Mv = M[i];
      const gt_bmp_vector_t mask = level_mask[i];
      const gt_bmp_vector_t Peq = peq[GT_BPM_PEQ_IDX(enc_char,i,num_words)];
      /* Compute Block */
      GT_MAP_BLOCK_BPM_ADVANCE_BLOCK(Peq,mask,Pv,Mv,PHin,MHin,PHout,MHout);
      /* Save Block Pv,Mv */
      P[i]=Pv;
      M[i]=Mv;
      /* Adjust score and swap propagate Hv */
      score[i] += PHout-MHout;
      PHin=PHout;
      MHin=MHout;
    }

    // Cut-off
    const uint64_t last = top_level-1;
    if (gt_expect_false(score[last]<=max_distance__1)) {
      const uint64_t last_score = score[last]+(MHin-PHin);
      const gt_bmp_vector_t Peq=peq[GT_BPM_PEQ_IDX(enc_char,top_level,num_words)];
      if (last_score<=max_distance && last<top && (MHin || (Peq & 1))) {
        // Init block V
        gt_bmp_vector_t Pv = word_ones;
        gt_bmp_vector_t Mv = 0;
        const gt_bmp_vector_t mask = level_mask[top_level];
        /* Compute Block */
        GT_MAP_BLOCK_BPM_ADVANCE_BLOCK(Peq,mask,Pv,Mv,PHin,MHin,PHout,MHout);
        /* Save Block Pv,Mv */
        P[top_level]=Pv;
        M[top_level]=Mv;
        /* Set score & increment the top level block */
        score[top_level] = last_score + init_score[top_level] + (PHout-MHout);
        ++top_level;
      } else {
        while (score[top_level-1] > (max_distance+init_score[top_level-1])) {
          --top_level;
        }
      }
    } else {
      while (score[top_level-1] > (max_distance+init_score[top_level-1])) {
        --top_level;
      }
    }

    // Check match
    if (top_level==num_words && score[top_level-1]<=max_distance) {
      if (score[top_level-1]<min_score)  {
        min_score_position = sequence_position;
        min_score = score[top_level-1];
      }
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

/* Clean up*/
#undef gt_bmp_vector_t
#undef GT_BMP_WORD_ONES
#undef GT_BMP_WORD_LENGTH
#undef GT_MAP_BLOCK_BPM_ADVANCE_BLOCK

#endif
