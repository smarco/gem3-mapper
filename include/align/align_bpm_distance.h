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
 *   Alignment module using BPM-algorithm to compute levenshtein distance
 *   (Myers' Fast Bit-Vector algorithm to compute levenshtein distance)
 */

#ifndef ALIGN_BPM_DISTANCE_H_
#define ALIGN_BPM_DISTANCE_H_

#include "utils/essentials.h"
#include "align/align_bpm_pattern.h"
#include "filtering/candidates/filtering_candidates.h"

/*
 * Constants
 */
#define BMP_W64_LENGTH UINT64_LENGTH
#define BMP_W64_ONES   UINT64_MAX
#define BMP_W64_MASK   (1ull<<63)

/*
 * Advance block functions (Improved)
 *   const @vector Eq,mask;
 *   return (Pv,Mv,PHout,MHout);
 */
#define BPM_ADVANCE_BLOCK(Eq,mask,Pv,Mv,PHin,MHin,PHout,MHout) \
  /* Computes modulator vector {Xv,Xh} ( cases A&C ) */ \
  const uint64_t Xv = Eq | Mv; \
  const uint64_t _Eq = Eq | MHin; \
  const uint64_t Xh = (((_Eq & Pv) + Pv) ^ Pv) | _Eq; \
  /* Calculate Hout */ \
  uint64_t Ph = Mv | ~(Xh | Pv); \
  uint64_t Mh = Pv & Xh; \
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
#define BPM_ADVANCE_COLUMN(P,M,score,PEQ,enc_char,top_level,PHin,MHin,PHout,MHout) { \
  uint64_t i; \
  for (i=0;i<top_level;++i) { \
    uint64_t Pv = P[i]; \
    uint64_t Mv = M[i]; \
    const uint64_t mask = level_mask[i]; \
    const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i,enc_char)]; \
    /* Compute Block */ \
    BPM_ADVANCE_BLOCK(Eq,mask,Pv,Mv,PHin,MHin,PHout,MHout); \
    /* Save Block Pv,Mv */ \
    P[i]=Pv; \
    M[i]=Mv; \
    /* Adjust score and swap propagate Hv */ \
    score[i] += PHout-MHout; \
    PHin=PHout; \
    MHin=MHout; \
  } \
}
#define BPM_CUT_OFF(P,M,score,PEQ,enc_char,top_level,PHin,MHin,PHout,MHout) { \
  const uint8_t last = top_level-1; \
  if (gem_expect_false(score[last]<=max_distance__1)) { \
    const uint64_t last_score = score[last]+(MHin-PHin); \
    if (last_score<=max_distance && last<top && (MHin || (PEQ[BPM_PATTERN_PEQ_IDX(top_level,enc_char)] & 1))) { \
      const uint64_t Peq = PEQ[BPM_PATTERN_PEQ_IDX(top_level,enc_char)]; \
      /* Init block V */ \
      uint64_t Pv = BMP_W64_ONES; \
      uint64_t Mv = 0; \
      const uint64_t mask = level_mask[top_level]; \
      /* Compute Block */ \
      BPM_ADVANCE_BLOCK(Peq,mask,Pv,Mv,PHin,MHin,PHout,MHout); \
      /* Save Block Pv,Mv */ \
      P[top_level]=Pv; \
      M[top_level]=Mv; \
      /* Set score & increment the top level block */ \
      score[top_level] = last_score + init_score[top_level] + (PHout-MHout); \
      ++top_level; \
    } else { \
      while (score[top_level-1] > (max_distance+init_score[top_level-1])) { \
        --top_level; \
      } \
    } \
  } else { \
    while (score[top_level-1] > (max_distance+init_score[top_level-1])) { \
      --top_level; \
    } \
  } \
}

/*
 * Setup
 */
void bpm_reset_search(
    const uint64_t num_words,
    uint64_t* const P,
    uint64_t* const M,
    int64_t* const score,
    const int64_t* const init_score);
void bpm_reset_search_cutoff(
    uint8_t* const top_level,
    uint64_t* const P,
    uint64_t* const M,
    int64_t* const score,
    const int64_t* const init_score,
    const uint64_t max_distance);

/*
 * BPM (Compute edit distance)
 */
// Raw
bool bpm_compute_edit_distance_raw(
    bpm_pattern_t* const bpm_pattern,
    const uint8_t* const text,
    const uint64_t text_length,
    uint64_t* const position,
    uint64_t* const distance);
// Cut-off
bool bpm_compute_edit_distance(
    const bpm_pattern_t* const bpm_pattern,
    const uint8_t* const text,
    const uint64_t text_length,
    uint64_t* const match_distance,
    uint64_t* const match_column,
    uint64_t max_distance,
    const bool quick_abandon);
// Find all local minimums
uint64_t bpm_compute_edit_distance_all(
    pattern_t* const pattern,
    filtering_candidates_t* const filtering_candidates,
    const uint64_t begin_position,
    const uint8_t* const text,
    const uint64_t text_length,
    uint64_t max_distance);

#endif /* ALIGN_BPM_DISTANCE_H_ */
