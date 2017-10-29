/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_align_bpm.c
 * DATE: 20/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Myers' Bit-vector algorithm.
 *   Straightforward implementation on 64bits word-vectors
 */

#include "gt_map_align_bpm.h"

// Constants
#define GT_BMP_W64_LENGTH 64
#define GT_BMP_W64_ONES UINT64_MAX

/*
 * Reset search functions
 */
GT_INLINE void gt_map_block_bpm_reset_search(
    const uint64_t num_words,uint64_t* const P,uint64_t* const M,int64_t* const score,
    const int64_t* const init_score,const uint64_t max_distance) {
  // Reset score,P,M
  uint64_t i;
  P[0]=GT_BMP_W64_ONES;
  M[0]=0;
  score[0] = init_score[0];
  for (i=1;i<num_words;++i) {
    P[i]=GT_BMP_W64_ONES;
    M[i]=0;
    score[i] = score[i-1] + init_score[i];
  }
}
GT_INLINE void gt_map_block_bpm_reset_search__cutoff(
    uint8_t* const top_level,uint64_t* const P,uint64_t* const M,int64_t* const score,
    const int64_t* const init_score,const uint64_t max_distance) {
  // Calculate the top level (maximum bit-word for cut-off purposes)
  const uint8_t y = (max_distance>0) ? (max_distance+(GT_BMP_W64_LENGTH-1))/GT_BMP_W64_LENGTH : 1;
  *top_level = y;
  // Reset score,P,M
  uint64_t i;
  P[0]=GT_BMP_W64_ONES;
  M[0]=0;
  score[0] = init_score[0];
  for (i=1;i<y;++i) {
    P[i]=GT_BMP_W64_ONES;
    M[i]=0;
    score[i] = score[i-1] + init_score[i];
  }
}
/*
 * Advance block functions
 */
int8_t T_hout_64[2][2] = {{0,-1},{1,1}};
uint8_t P_hin_64[3] = {0, 0, 1L};
uint8_t N_hin_64[3] = {1L, 0, 0};
GT_INLINE int8_t gt_map_block_bpm_advance_block(
    uint64_t Eq,const uint64_t mask,
    uint64_t Pv,uint64_t Mv,const int8_t hin,
    uint64_t* const Pv_out,uint64_t* const Mv_out) {
  uint64_t Ph, Mh;
  uint64_t Xv, Xh;
  int8_t hout=0;

  Xv = Eq | Mv;
  Eq |= N_hin_64[hin];
  Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;

  Ph = Mv | ~(Xh | Pv);
  Mh = Pv & Xh;

  hout += T_hout_64[(Ph & mask)!=0][(Mh & mask)!=0];

  Ph <<= 1;
  Mh <<= 1;

  Mh |= N_hin_64[hin];
  Ph |= P_hin_64[hin];

  Pv = Mh | ~(Xv | Ph);
  Mv = Ph & Xv;

  *Pv_out=Pv;
  *Mv_out=Mv;

  return hout;
}
/*
 * Bit-compressed (Re)alignment
 *   BMP[BitParalellMayers] - Myers' Fast Bit-Vector algorithm (Levenshtein)
 */
GT_INLINE bool gt_map_block_bpm_get_distance(
    gt_bpm_pattern* const bpm_pattern,char* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance,const uint64_t max_distance) {
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
  gt_map_block_bpm_reset_search(num_words,P,M,score,init_score,max_distance);

  // Advance in DP-bit_encoded matrix
  uint64_t sequence_position;
  for (sequence_position=0;sequence_position<sequence_length;++sequence_position) {
    // Fetch next character
    const uint8_t enc_char = gt_cdna_encode(sequence[sequence_position]);

    // Advance all blocks
    int8_t carry;
    uint64_t i;
    for (i=0,carry=0;i<num_words;++i) {
      uint64_t* const Py = P+i;
      uint64_t* const My = M+i;
      carry = gt_map_block_bpm_advance_block(
          peq[GT_BPM_PEQ_IDX(enc_char,i,num_words)],level_mask[i],*Py,*My,carry+1,Py,My);
      score[i] += carry;
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
GT_INLINE bool gt_map_block_bpm_get_distance__cutoff(
    gt_bpm_pattern* const bpm_pattern,char* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance,const uint64_t max_distance) {
  // Pattern variables
  const uint64_t* peq = (uint64_t*)bpm_pattern->peq;
  const uint64_t num_words = bpm_pattern->pattern_num_words;
  uint64_t* const P = (uint64_t*)bpm_pattern->P;
  uint64_t* const M = (uint64_t*)bpm_pattern->M;
  const uint64_t* const level_mask = (uint64_t*)bpm_pattern->level_mask;
  int64_t* const score = (int64_t*)bpm_pattern->score;
  const int64_t* const init_score = (int64_t*)bpm_pattern->init_score;

  // Initialize search
  const uint8_t top = num_words-1;
  uint8_t top_level;
  uint64_t min_score = UINT64_MAX, min_score_position = UINT64_MAX;
  gt_map_block_bpm_reset_search__cutoff(&top_level,P,M,score,init_score,max_distance);

  // Advance in DP-bit_encoded matrix
  uint64_t sequence_position;
  for (sequence_position=0;sequence_position<sequence_length;++sequence_position) {
    // Fetch next character
    const uint8_t enc_char = gt_cdna_encode(sequence[sequence_position]);

    // Advance all blocks
    int8_t carry;
    uint64_t i;
    for (i=0,carry=0;i<top_level;++i) {
      uint64_t* const Py = P+i;
      uint64_t* const My = M+i;
      carry = gt_map_block_bpm_advance_block(
          peq[GT_BPM_PEQ_IDX(enc_char,i,num_words)],level_mask[i],*Py,*My,carry+1,Py,My);
      score[i] += carry;
    }

    // Cut-off
    const uint8_t last = top_level-1;
    if ((score[last]-carry)<=max_distance && last<top &&
        ( (peq[GT_BPM_PEQ_IDX(enc_char,top_level,num_words)] & 1) || (carry<0) )  ) {
      // Init block V
      P[top_level]=GT_BMP_W64_ONES;
      M[top_level]=0;

      uint64_t* const Py = P+top_level;
      uint64_t* const My = M+top_level;
      score[top_level] = score[top_level-1] + init_score[top_level] - carry +
          gt_map_block_bpm_advance_block(peq[GT_BPM_PEQ_IDX(enc_char,top_level,num_words)],
              level_mask[top_level],*Py,*My,carry+1,Py,My);
      ++top_level;
    } else {
      while (score[top_level-1]>max_distance+init_score[top_level-1]) {
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
