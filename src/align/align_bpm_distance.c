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

#include "align/alignment.h"
#include "text/dna_text.h"
#include "align/pattern/pattern.h"
#include "align/align_bpm_distance.h"
#include "filtering/candidates/filtering_candidates_accessors.h"
#include "filtering/region/filtering_region.h"

/*
 * Mode/Debug
 */
#define DEBUG_BPM_TILED false

/*
 * Reset search functions
 */
void bpm_reset_search(
    const uint64_t num_words,
    uint64_t* const P,
    uint64_t* const M,
    int64_t* const score,
    const int64_t* const init_score) {
  // Reset score,P,M
  uint64_t i;
  P[0]=BMP_W64_ONES;
  M[0]=0;
  score[0] = init_score[0];
  for (i=1;i<num_words;++i) {
    P[i]=BMP_W64_ONES;
    M[i]=0;
    score[i] = score[i-1] + init_score[i];
  }
}
void bpm_reset_search_cutoff(
    uint8_t* const top_level,
    uint64_t* const P,
    uint64_t* const M,
    int64_t* const score,
    const int64_t* const init_score,
    const uint64_t max_distance) {
  // Calculate the top level (maximum bit-word for cut-off purposes)
  const uint8_t y = (max_distance>0) ? (max_distance+(BMP_W64_LENGTH-1))/BMP_W64_LENGTH : 1;
  *top_level = y;
  // Reset score,P,M
  uint64_t i;
  P[0]=BMP_W64_ONES;
  M[0]=0;
  score[0] = init_score[0];
  for (i=1;i<y;++i) {
    P[i]=BMP_W64_ONES;
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
int8_t bpm_advance_block(
    uint64_t Eq,
    const uint64_t mask,
    uint64_t Pv,
    uint64_t Mv,
    const int8_t hin,
    uint64_t* const Pv_out,
    uint64_t* const Mv_out) {
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
 * BMP
 */
bool bpm_compute_edit_distance_raw(
    bpm_pattern_t* const bpm_pattern,
    const uint8_t* const text,
    const uint64_t text_length,
    uint64_t* const position,
    uint64_t* const distance) {
  // Pattern variables
  const uint64_t* PEQ = bpm_pattern->PEQ;
  const uint64_t num_words64 = bpm_pattern->pattern_num_words64;
  uint64_t* const P = bpm_pattern->P;
  uint64_t* const M = bpm_pattern->M;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  int64_t* const score = bpm_pattern->score;
  const int64_t* const init_score = bpm_pattern->init_score;
  // Initialize search
  uint64_t min_score = ALIGN_DISTANCE_INF, min_score_column = ALIGN_COLUMN_INF;
  bpm_reset_search(num_words64,P,M,score,init_score);
  // Advance in DP-bit_encoded matrix
  uint64_t text_position;
  for (text_position=0;text_position<text_length;++text_position) {
    // Fetch next character
    const uint8_t enc_char = text[text_position];
    // Advance all blocks
    int8_t carry;
    uint64_t i;
    for (i=0,carry=0;i<num_words64;++i) {
      uint64_t* const Py = P+i;
      uint64_t* const My = M+i;
      carry = bpm_advance_block(PEQ[BPM_PATTERN_PEQ_IDX(i,enc_char)],level_mask[i],*Py,*My,carry+1,Py,My);
      score[i] += carry;
    }
    // Check match
    if (score[num_words64-1] < min_score) {
      min_score_column = text_position;
      min_score = score[num_words64-1];
    }
  }
  // Return results
  if (min_score!=ALIGN_DISTANCE_INF) {
    *distance = min_score;
    *position = min_score_column;
    return true;
  } else {
    *distance = ALIGN_DISTANCE_INF;
    *position = ALIGN_COLUMN_INF;
    return false;
  }
}
bool bpm_compute_edit_distance(
    const bpm_pattern_t* const bpm_pattern,
    const uint8_t* const text,
    const uint64_t text_length,
    uint64_t* const match_distance,
    uint64_t* const match_column,
    uint64_t max_distance,
    const bool quick_abandon) {
  PROF_START(GP_BPM_DISTANCE);
  PROF_ADD_COUNTER(GP_BPM_DISTANCE_KEY_LENGTH,bpm_pattern->pattern_length);
  PROF_ADD_COUNTER(GP_BPM_DISTANCE_TEXT_LENGTH,text_length);
  PROF_ADD_COUNTER(GP_BPM_DISTANCE_CELLS,bpm_pattern->pattern_length*text_length);
  // Pattern variables
  const uint64_t* PEQ = bpm_pattern->PEQ;
  const uint64_t num_words64 = bpm_pattern->pattern_num_words64;
  uint64_t* const P = bpm_pattern->P;
  uint64_t* const M = bpm_pattern->M;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  int64_t* const score = bpm_pattern->score;
  const int64_t* const init_score = bpm_pattern->init_score;
  const uint64_t* const pattern_left = bpm_pattern->pattern_left;
  // Initialize search
  if (max_distance >= bpm_pattern->pattern_length) {
    max_distance = bpm_pattern->pattern_length-1; // Correct max-distance
  }
  const uint64_t max_distance__1 = max_distance+1;
  const uint8_t top = num_words64-1;
  uint8_t top_level;
  uint64_t min_score = ALIGN_DISTANCE_INF, min_score_column = ALIGN_COLUMN_INF;
  bpm_reset_search_cutoff(&top_level,P,M,score,init_score,max_distance);
  // Advance in DP-bit_encoded matrix
  uint64_t text_position, text_left=text_length;
  for (text_position=0;text_position<text_length;++text_position,--text_left) {
    // Fetch next character
    const uint8_t enc_char = text[text_position];
    // Advance all blocks
    uint64_t i,PHin=0,MHin=0,PHout,MHout;
    for (i=0;i<top_level;++i) {
      uint64_t Pv = P[i];
      uint64_t Mv = M[i];
      const uint64_t mask = level_mask[i];
      const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i,enc_char)];
      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv,Mv,PHin,MHin,PHout,MHout);
      /* Save Block Pv,Mv */
      P[i]=Pv;
      M[i]=Mv;
      /* Adjust score and swap propagate Hv */
      score[i] += PHout-MHout;
      PHin=PHout;
      MHin=MHout;
    }
    // Cut-off
    const uint8_t last = top_level-1;
    if (gem_expect_false(score[last]<=max_distance__1 && last<top)) {
      const uint64_t last_score = score[last]+(MHin-PHin);
      const uint64_t Peq = PEQ[BPM_PATTERN_PEQ_IDX(top_level,enc_char)];
      if (last_score<=max_distance && (MHin || (Peq & 1))) {
        // Init block V
        uint64_t Pv = BMP_W64_ONES;
        uint64_t Mv = 0;
        const uint64_t mask = level_mask[top_level];
        /* Compute Block */
        BPM_ADVANCE_BLOCK(Peq,mask,Pv,Mv,PHin,MHin,PHout,MHout);
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
    const int64_t current_score = score[top_level-1];
    if (top_level==num_words64 && current_score<=max_distance) {
      if (current_score < min_score)  {
        min_score_column = text_position;
        min_score = current_score;
      }
    } else if (quick_abandon && min_score==ALIGN_DISTANCE_INF &&
               current_score+pattern_left[top_level] > text_left+max_distance) {
      // Quick abandon, it doesn't match (bounded by best case scenario)
      PROF_STOP(GP_BPM_DISTANCE);
      PROF_INC_COUNTER(GP_BPM_DISTANCE_QUICK_ABANDON);
      *match_distance = ALIGN_DISTANCE_INF;
      return false;
    }
  }
  // Return results
  if (min_score!=ALIGN_DISTANCE_INF) {
    *match_distance = min_score;
    *match_column = min_score_column;
    PROF_STOP(GP_BPM_DISTANCE);
    return true;
  } else {
    *match_distance = ALIGN_DISTANCE_INF;
    PROF_STOP(GP_BPM_DISTANCE);
    return false;
  }
}
/*
 * BPM all matches
 */
uint64_t bpm_compute_edit_distance_all(
    pattern_t* const pattern,
    filtering_candidates_t* const filtering_candidates,
    const uint64_t begin_position,
    const uint8_t* const text,
    const uint64_t text_length,
    uint64_t max_distance) {
  PROF_START(GP_BPM_ALL);
  // Parameters
  bpm_pattern_t* const bpm_pattern = &pattern->pattern_tiled.bpm_pattern;
  // Pattern variables
  const uint64_t* PEQ = bpm_pattern->PEQ;
  const uint64_t num_words64 = bpm_pattern->pattern_num_words64;
  const uint64_t key_length = bpm_pattern->pattern_length;
  uint64_t* const P = bpm_pattern->P;
  uint64_t* const M = bpm_pattern->M;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  int64_t* const score = bpm_pattern->score;
  const int64_t* const init_score = bpm_pattern->init_score;
  const uint64_t* const pattern_left = bpm_pattern->pattern_left;
  // Initialize search
  if (max_distance >= bpm_pattern->pattern_length) {
    max_distance = bpm_pattern->pattern_length-1; // Correct max-distance
  }
  const uint64_t max_distance__1 = max_distance+1;
  const uint8_t top = num_words64-1;
  uint8_t top_level;
  uint64_t min_score = ALIGN_DISTANCE_INF, min_score_column = ALIGN_COLUMN_INF, opt_steps_left;
  bpm_reset_search_cutoff(&top_level,P,M,score,init_score,max_distance);
  // Advance in DP-bit_encoded matrix
  const uint64_t end_position = begin_position + text_length;
  bool match_found = false;
  uint64_t text_position, text_left=text_length, num_matches_found=0;
  for (text_position=0;text_position<text_length;++text_position,--text_left) {
    // Fetch next character and advance all blocks
    uint64_t PHin=0,MHin=0,PHout,MHout;
    const uint8_t enc_char = text[text_position];
    BPM_ADVANCE_COLUMN(P,M,score,PEQ,enc_char,top_level,PHin,MHin,PHout,MHout);
    // Cut-off
    BPM_CUT_OFF(P,M,score,PEQ,enc_char,top_level,PHin,MHin,PHout,MHout);
    // Check match
    const int64_t current_score = score[top_level-1];
    if (!match_found) {
      if (top_level==num_words64 && current_score<=max_distance) {
        min_score_column = text_position;
        min_score = current_score;
        if (current_score==0) { // Don't try to optimize (exact match)
          const uint64_t text_end_offset = min_score_column+1;
          const uint64_t text_begin_offset = BOUNDED_SUBTRACTION(text_end_offset,key_length+min_score,0);
          filtering_candidates_add_region_verified(
              filtering_candidates,pattern,text_begin_offset,
              text_end_offset,begin_position,end_position,min_score);
          ++num_matches_found; // Increment the number of matches found
        } else {
          match_found = true;
          opt_steps_left = max_distance; // Setup optimization steps
        }
      }
    } else {
      if (top_level==num_words64 && current_score<min_score) { // Update minimum
        min_score_column = text_position;
        min_score = current_score;
      }
      if (opt_steps_left==0) {
        const uint64_t text_end_offset = min_score_column+1;
        const uint64_t text_begin_offset = BOUNDED_SUBTRACTION(text_end_offset,key_length+min_score,0);
        filtering_candidates_add_region_verified(
            filtering_candidates,pattern,text_begin_offset,
            text_end_offset,begin_position,end_position,min_score);
        ++num_matches_found; // Increment the number of matches found
        match_found = false;
      } else {
        --opt_steps_left;
      }
    }
    // Quick abandon
    if (min_score==ALIGN_DISTANCE_INF && current_score+pattern_left[top_level] > text_left+max_distance) {
      PROF_INC_COUNTER(GP_BPM_ALL_QUICK_ABANDON);
      break; // Quick abandon, it doesn't match (bounded by best case scenario)
    }
  }
  if (match_found) {
    const uint64_t text_end_offset = min_score_column+1;
    const uint64_t text_begin_offset = BOUNDED_SUBTRACTION(text_end_offset,key_length+min_score,0);
    filtering_candidates_add_region_verified(
        filtering_candidates,pattern,text_begin_offset,
        text_end_offset,begin_position,end_position,min_score);
    ++num_matches_found; // Increment the number of matches found
  }
  PROF_INC_COUNTER(GP_BPM_ALL_MATCHES_FOUND);
  PROF_STOP(GP_BPM_ALL);
  return num_matches_found;
}
