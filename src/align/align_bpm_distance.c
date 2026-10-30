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
/*
 * Cut-off kernels for patterns of at most 4 words (BPM_W4_MAX_WORDS)
 *   P, M and score are kept in locals (P0..P3, M0..M3, score0..score3) and
 *   written back on return; each column is unrolled for its top_level.
 *   gcc cannot tell P, M and score from PEQ and level_mask in memory, so the
 *   array loops reload and re-store them every column.
 */
#define BPM_W4_MAX_WORDS 4
/*
 * A text byte above ENC_DNA_CHAR_N indexes PEQ past the word's rows; on the
 * last word that reads P/M in memory, which the kernels keep in registers.
 * Such texts take the array version.
 */
static inline bool bpm_text_within_peq(
    const uint8_t* const text,
    const uint64_t text_length) {
  uint8_t above = 0;
  uint64_t i;
  for (i=0;i<text_length;++i) above |= (text[i] > ENC_DNA_CHAR_N);
  return above==0;
}
#define BPM_W4_SETUP() \
  const uint64_t* const level_mask = bpm_pattern->level_mask; \
  const int64_t* const init_score = bpm_pattern->init_score; \
  const uint64_t mask0 = level_mask[0]; \
  const uint64_t mask1 = (num_words64>1) ? level_mask[1] : 0; \
  const uint64_t mask2 = (num_words64>2) ? level_mask[2] : 0; \
  const uint64_t mask3 = (num_words64>3) ? level_mask[3] : 0; \
  const int64_t init_score0 = init_score[0]; \
  const int64_t init_score1 = (num_words64>1) ? init_score[1] : 0; \
  const int64_t init_score2 = (num_words64>2) ? init_score[2] : 0; \
  const int64_t init_score3 = (num_words64>3) ? init_score[3] : 0; \
  const uint64_t max_distance__1 = max_distance+1; \
  const uint8_t top = num_words64-1; \
  /* bpm_reset_search_cutoff */ \
  const uint8_t y = (max_distance>0) ? (max_distance+(BMP_W64_LENGTH-1))/BMP_W64_LENGTH : 1; \
  uint8_t top_level = y, max_level = y; \
  uint64_t P0 = BMP_W64_ONES, P1 = 0, P2 = 0, P3 = 0; \
  uint64_t M0 = 0, M1 = 0, M2 = 0, M3 = 0; \
  int64_t score0 = init_score0, score1 = 0, score2 = 0, score3 = 0; \
  if (y>1) { P1 = BMP_W64_ONES; score1 = score0 + init_score1; } \
  if (y>2) { P2 = BMP_W64_ONES; score2 = score1 + init_score2; } \
  if (y>3) { P3 = BMP_W64_ONES; score3 = score2 + init_score3; }
/* Words below max_level are those the array version would have written */
#define BPM_W4_WRITE_BACK() { \
  P[0] = P0; M[0] = M0; score[0] = score0; \
  if (num_words64>1 && max_level>1) { P[1] = P1; M[1] = M1; score[1] = score1; } \
  if (num_words64>2 && max_level>2) { P[2] = P2; M[2] = M2; score[2] = score2; } \
  if (num_words64>3 && max_level>3) { P[3] = P3; M[3] = M3; score[3] = score3; } \
}
#define BPM_W4_STEP(k) { \
  uint64_t Pv = P##k; \
  uint64_t Mv = M##k; \
  const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(k,enc_char)]; \
  BPM_ADVANCE_BLOCK(Eq,mask##k,Pv,Mv,PHin,MHin,PHout,MHout); \
  P##k = Pv; \
  M##k = Mv; \
  score##k += PHout-MHout; \
  PHin = PHout; \
  MHin = MHout; \
}
#define BPM_W4_EXTEND(k) { \
  const uint64_t Peq = PEQ[BPM_PATTERN_PEQ_IDX(k,enc_char)]; \
  uint64_t Pv = BMP_W64_ONES; \
  uint64_t Mv = 0; \
  BPM_ADVANCE_BLOCK(Peq,mask##k,Pv,Mv,PHin,MHin,PHout,MHout); \
  P##k = Pv; \
  M##k = Mv; \
  score##k = last_score + init_score##k + (PHout-MHout); \
  top_level = k+1; \
  if (max_level < k+1) max_level = k+1; \
  current_score = score##k; \
}
/* The array version's while loop; it never goes below level 1, since score0 <= init_score0 */
#define BPM_W4_DOWN_1 { top_level = 1; current_score = score0; }
#define BPM_W4_DOWN_2 { \
  if (score1 > (max_distance+init_score1)) BPM_W4_DOWN_1 \
  else { top_level = 2; current_score = score1; } \
}
#define BPM_W4_DOWN_3 { \
  if (score2 > (max_distance+init_score2)) BPM_W4_DOWN_2 \
  else { top_level = 3; current_score = score2; } \
}
#define BPM_W4_DOWN_4 { \
  if (score3 > (max_distance+init_score3)) BPM_W4_DOWN_3 \
  else { top_level = 4; current_score = score3; } \
}
/* Cut-off at top_level T (last word L=T-1); T=4 cannot extend (last<top fails) */
#define BPM_W4_CUT_OFF(T,L) \
  if (gem_expect_false(score##L<=max_distance__1)) { \
    const uint64_t last_score = score##L+(MHin-PHin); \
    if (last_score<=max_distance && L<top && (MHin || (PEQ[BPM_PATTERN_PEQ_IDX(T,enc_char)] & 1))) BPM_W4_EXTEND(T) \
    else BPM_W4_DOWN_##T \
  } else BPM_W4_DOWN_##T
#define BPM_W4_COLUMN() \
  if (num_words64==1 || top_level==1) { \
    BPM_W4_STEP(0) \
    BPM_W4_CUT_OFF(1,0) \
  } else if (num_words64==2 || top_level==2) { \
    BPM_W4_STEP(0) BPM_W4_STEP(1) \
    BPM_W4_CUT_OFF(2,1) \
  } else if (num_words64==3 || top_level==3) { \
    BPM_W4_STEP(0) BPM_W4_STEP(1) BPM_W4_STEP(2) \
    BPM_W4_CUT_OFF(3,2) \
  } else { \
    BPM_W4_STEP(0) BPM_W4_STEP(1) BPM_W4_STEP(2) BPM_W4_STEP(3) \
    BPM_W4_DOWN_4 \
  }
static inline __attribute__((always_inline)) bool bpm_compute_edit_distance_wn(
    const bpm_pattern_t* const bpm_pattern,
    const uint8_t* const text,
    const uint64_t text_length,
    uint64_t* const match_distance,
    uint64_t* const match_column,
    const uint64_t max_distance,
    const bool quick_abandon,
    const uint64_t num_words64) {
  // Pattern variables
  const uint64_t* const PEQ = bpm_pattern->PEQ;
  uint64_t* const P = bpm_pattern->P;
  uint64_t* const M = bpm_pattern->M;
  int64_t* const score = bpm_pattern->score;
  const uint64_t pattern_length = bpm_pattern->pattern_length;
  // Initialize search
  BPM_W4_SETUP();
  uint64_t min_score = ALIGN_DISTANCE_INF, min_score_column = ALIGN_COLUMN_INF;
  // Advance in DP-bit_encoded matrix
  uint64_t text_position, text_left=text_length;
  for (text_position=0;text_position<text_length;++text_position,--text_left) {
    // Fetch next character, advance all blocks and cut-off
    const uint8_t enc_char = text[text_position];
    uint64_t PHin=0,MHin=0,PHout,MHout;
    int64_t current_score;
    BPM_W4_COLUMN();
    // Check match
    if (top_level==num_words64 && current_score<=max_distance) {
      if (current_score < min_score)  {
        min_score_column = text_position;
        min_score = current_score;
      }
    } else if (quick_abandon && min_score==ALIGN_DISTANCE_INF) {
      // Rows below the active blocks; tiles have no pattern_left array
      const uint64_t active_rows = (uint64_t)top_level*BMP_W64_LENGTH;
      const uint64_t rows_left = (pattern_length > active_rows) ? pattern_length-active_rows : 0;
      if (current_score+rows_left > text_left+max_distance) {
        // Quick abandon, it doesn't match (bounded by best case scenario)
        BPM_W4_WRITE_BACK();
        PROF_STOP(GP_BPM_DISTANCE);
        PROF_INC_COUNTER(GP_BPM_DISTANCE_QUICK_ABANDON);
        *match_distance = ALIGN_DISTANCE_INF;
        return false;
      }
    }
  }
  BPM_W4_WRITE_BACK();
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
static inline __attribute__((always_inline)) uint64_t bpm_compute_edit_distance_all_wn(
    pattern_t* const pattern,
    filtering_candidates_t* const filtering_candidates,
    const uint64_t begin_position,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t max_distance,
    const uint64_t num_words64) {
  // Parameters
  bpm_pattern_t* const bpm_pattern = &pattern->pattern_tiled.bpm_pattern;
  // Pattern variables
  const uint64_t* const PEQ = bpm_pattern->PEQ;
  const uint64_t key_length = bpm_pattern->pattern_length;
  uint64_t* const P = bpm_pattern->P;
  uint64_t* const M = bpm_pattern->M;
  int64_t* const score = bpm_pattern->score;
  const uint64_t* const pattern_left = bpm_pattern->pattern_left;
  // Initialize search
  BPM_W4_SETUP();
  uint64_t min_score = ALIGN_DISTANCE_INF, min_score_column = ALIGN_COLUMN_INF, opt_steps_left = 0;
  // Advance in DP-bit_encoded matrix
  const uint64_t end_position = begin_position + text_length;
  bool match_found = false;
  uint64_t text_position, text_left=text_length, num_matches_found=0;
  for (text_position=0;text_position<text_length;++text_position,--text_left) {
    // Fetch next character, advance all blocks and cut-off
    uint64_t PHin=0,MHin=0,PHout,MHout;
    const uint8_t enc_char = text[text_position];
    int64_t current_score;
    BPM_W4_COLUMN();
    // Check match
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
  BPM_W4_WRITE_BACK();
  if (match_found) {
    const uint64_t text_end_offset = min_score_column+1;
    const uint64_t text_begin_offset = BOUNDED_SUBTRACTION(text_end_offset,key_length+min_score,0);
    filtering_candidates_add_region_verified(
        filtering_candidates,pattern,text_begin_offset,
        text_end_offset,begin_position,end_position,min_score);
    ++num_matches_found; // Increment the number of matches found
  }
  return num_matches_found;
}
/* Out of line, so the column loop gets registers of its own */
#define BPM_W4_ALL_KERNEL(NW) \
static __attribute__((noinline)) uint64_t bpm_compute_edit_distance_all_w##NW( \
    pattern_t* const pattern,filtering_candidates_t* const filtering_candidates, \
    const uint64_t begin_position,const uint8_t* const text,const uint64_t text_length, \
    const uint64_t max_distance) { \
  return bpm_compute_edit_distance_all_wn(pattern,filtering_candidates, \
      begin_position,text,text_length,max_distance,NW); \
}
BPM_W4_ALL_KERNEL(1)
BPM_W4_ALL_KERNEL(2)
BPM_W4_ALL_KERNEL(3)
BPM_W4_ALL_KERNEL(4)
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
  const uint64_t pattern_length = bpm_pattern->pattern_length;
  // Initialize search
  if (max_distance >= bpm_pattern->pattern_length) {
    max_distance = bpm_pattern->pattern_length-1; // Correct max-distance
  }
  if (num_words64 <= BPM_W4_MAX_WORDS && bpm_text_within_peq(text,text_length)) {
    switch (num_words64) {
      case 1: return bpm_compute_edit_distance_wn(bpm_pattern,text,text_length,
          match_distance,match_column,max_distance,quick_abandon,1);
      case 2: return bpm_compute_edit_distance_wn(bpm_pattern,text,text_length,
          match_distance,match_column,max_distance,quick_abandon,2);
      case 3: return bpm_compute_edit_distance_wn(bpm_pattern,text,text_length,
          match_distance,match_column,max_distance,quick_abandon,3);
      case 4: return bpm_compute_edit_distance_wn(bpm_pattern,text,text_length,
          match_distance,match_column,max_distance,quick_abandon,4);
      default: break;
    }
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
    } else if (quick_abandon && min_score==ALIGN_DISTANCE_INF) {
      // Rows below the active blocks; tiles have no pattern_left array
      const uint64_t active_rows = (uint64_t)top_level*BMP_W64_LENGTH;
      const uint64_t rows_left = (pattern_length > active_rows) ? pattern_length-active_rows : 0;
      if (current_score+rows_left > text_left+max_distance) {
        // Quick abandon, it doesn't match (bounded by best case scenario)
        PROF_STOP(GP_BPM_DISTANCE);
        PROF_INC_COUNTER(GP_BPM_DISTANCE_QUICK_ABANDON);
        *match_distance = ALIGN_DISTANCE_INF;
        return false;
      }
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
 *   A window of at least BPM_ALL_SPLIT_MIN_LENGTH columns, for a pattern of 1-4
 *   words, runs two searches in the two lanes of 2 x 64-bit vectors: A from
 *   column 0 and B, from a fresh state, from column s. The search state after
 *   a column (top_level and the active words' P/M/score; BPM_CUT_OFF
 *   overwrites a word when it activates it) decides every column after it, so
 *   when B's state equals A's after column s+X-1, A stops there and the match
 *   logic continues on B's recorded (top_level,score) columns. Otherwise A
 *   searches the rest alone. P/M/score in the pattern end as the single search
 *   leaves them.
 */
#define BPM_ALL_SPLIT_OVERLAP      256
#define BPM_ALL_SPLIT_MIN_LENGTH   (2*BPM_ALL_SPLIT_OVERLAP)
#define BPM_ALL_SPLIT_MAX_WORDS    4
#define BPM_ALL_SPLIT_MAX_RECORDS  4096
typedef uint64_t bpm_v2u_t __attribute__((vector_size(16)));
typedef int64_t bpm_v2i_t __attribute__((vector_size(16)));
#define BPM_ALL_CHECK_MATCH(top_level,current_score) \
  if (!match_found) { \
    if (top_level==num_words64 && current_score<=max_distance) { \
      min_score_column = text_position; \
      min_score = current_score; \
      if (current_score==0) { /* Don't try to optimize (exact match) */ \
        const uint64_t text_end_offset = min_score_column+1; \
        const uint64_t text_begin_offset = BOUNDED_SUBTRACTION(text_end_offset,key_length+min_score,0); \
        filtering_candidates_add_region_verified( \
            filtering_candidates,pattern,text_begin_offset, \
            text_end_offset,begin_position,end_position,min_score); \
        ++num_matches_found; /* Increment the number of matches found */ \
      } else { \
        match_found = true; \
        opt_steps_left = max_distance; /* Setup optimization steps */ \
      } \
    } \
  } else { \
    if (top_level==num_words64 && current_score<min_score) { /* Update minimum */ \
      min_score_column = text_position; \
      min_score = current_score; \
    } \
    if (opt_steps_left==0) { \
      const uint64_t text_end_offset = min_score_column+1; \
      const uint64_t text_begin_offset = BOUNDED_SUBTRACTION(text_end_offset,key_length+min_score,0); \
      filtering_candidates_add_region_verified( \
          filtering_candidates,pattern,text_begin_offset, \
          text_end_offset,begin_position,end_position,min_score); \
      ++num_matches_found; /* Increment the number of matches found */ \
      match_found = false; \
    } else { \
      --opt_steps_left; \
    } \
  }
#define BPM_ALL_QUICK_ABANDON(top_level,current_score) \
  (min_score==ALIGN_DISTANCE_INF && current_score+pattern_left[top_level] > text_left+max_distance)
// A's column text_position on the pattern's P/M/score, then its match logic
#define BPM_ALL_COLUMN_A() { \
  uint64_t PHin=0,MHin=0,PHout,MHout; \
  const uint8_t enc_char = text[text_position]; \
  BPM_ADVANCE_COLUMN(P,M,score,PEQ,enc_char,top_level,PHin,MHin,PHout,MHout); \
  BPM_CUT_OFF(P,M,score,PEQ,enc_char,top_level,PHin,MHin,PHout,MHout); \
  const int64_t current_score = score[top_level-1]; \
  BPM_ALL_CHECK_MATCH(top_level,current_score); \
  if (BPM_ALL_QUICK_ABANDON(top_level,current_score)) { \
    PROF_INC_COUNTER(GP_BPM_ALL_QUICK_ABANDON); \
    goto bpm_all_done; /* Quick abandon, it doesn't match (bounded by best case scenario) */ \
  } \
}
// Word k of both lanes (lane 0 A, lane 1 B); level masks are single bits (shift##k)
#define BPM_V2_BLOCK(k) \
  const bpm_v2u_t Eq = E##k[pair]; \
  const bpm_v2u_t Xv = Eq | M##k; \
  const bpm_v2u_t _Eq = Eq | MHin; \
  const bpm_v2u_t Xh = (((_Eq & P##k) + P##k) ^ P##k) | _Eq; \
  bpm_v2u_t Ph = M##k | ~(Xh | P##k); \
  bpm_v2u_t Mh = P##k & Xh; \
  const bpm_v2u_t PHout = (Ph >> shift##k) & one; \
  const bpm_v2u_t MHout = (Mh >> shift##k) & one; \
  Ph = (Ph << 1) | PHin; \
  Mh = (Mh << 1) | MHin; \
  const bpm_v2u_t Pn = Mh | ~(Xv | Ph); \
  const bpm_v2u_t Mn = Ph & Xv; \
  const bpm_v2i_t Sn = S##k + (bpm_v2i_t)(PHout-MHout);
#define BPM_V2_STEP(k) { \
  BPM_V2_BLOCK(k) \
  P##k = Pn; M##k = Mn; S##k = Sn; \
  PHin = PHout; MHin = MHout; \
}
// Lanes at different top levels: a word inactive in a lane keeps its values
#define BPM_V2_STEP_MIXED(k) { \
  BPM_V2_BLOCK(k) \
  P##k = (Pn & active##k) | (P##k & ~active##k); \
  M##k = (Mn & active##k) | (M##k & ~active##k); \
  S##k = (Sn & (bpm_v2i_t)active##k) | (S##k & ~(bpm_v2i_t)active##k); \
  PHsel = (PHout & last##k) | (PHsel & ~last##k); \
  MHsel = (MHout & last##k) | (MHsel & ~last##k); \
  Ssel = (Sn & (bpm_v2i_t)last##k) | (Ssel & ~(bpm_v2i_t)last##k); \
  PHin = PHout; MHin = MHout; \
}
#define BPM_V2_SET_WORD(k,L,Pv,Mv,Sv) { P##k[L] = Pv; M##k[L] = Mv; S##k[L] = Sv; }
#define BPM_V2_GET_SCORE(k,L) S##k[L]
// BPM_CUT_OFF on lane L (scalar); leaves the lane's score[top_level-1] in current_score
#define BPM_V2_CUT_OFF(L,enc_char,top_level,current_score) { \
  const uint8_t last = top_level-1; \
  const uint64_t PHl = PHin[L], MHl = MHin[L]; \
  int64_t score_last = current_score; \
  if (score_last<=max_distance__1) { \
    const uint64_t last_score = score_last+(MHl-PHl); \
    if (last_score<=max_distance && last<top && (MHl || (PEQ[BPM_PATTERN_PEQ_IDX(top_level,enc_char)] & 1))) { \
      const uint64_t Peq = PEQ[BPM_PATTERN_PEQ_IDX(top_level,enc_char)]; \
      uint64_t Pv = BMP_W64_ONES; \
      uint64_t Mv = 0; \
      uint64_t PHout, MHout; \
      const uint64_t mask = level_mask[top_level]; \
      BPM_ADVANCE_BLOCK(Peq,mask,Pv,Mv,PHl,MHl,PHout,MHout); \
      const int64_t score_top = last_score + init_score[top_level] + (PHout-MHout); \
      switch (top_level) { \
        case 1: BPM_V2_SET_WORD(1,L,Pv,Mv,score_top); break; \
        case 2: BPM_V2_SET_WORD(2,L,Pv,Mv,score_top); break; \
        default: BPM_V2_SET_WORD(3,L,Pv,Mv,score_top); break; \
      } \
      ++top_level; \
      current_score = score_top; \
    } /* else BPM_CUT_OFF's while loop cannot fire: score_last <= max_distance+1 <= max_distance+init_score */ \
  } else { \
    while (score_last > (max_distance+init_score[top_level-1])) { \
      --top_level; \
      switch (top_level) { \
        case 1: score_last = BPM_V2_GET_SCORE(0,L); break; \
        case 2: score_last = BPM_V2_GET_SCORE(1,L); break; \
        default: score_last = BPM_V2_GET_SCORE(2,L); break; \
      } \
    } \
    current_score = score_last; \
  } \
}
#define BPM_V2_LANE_MASKS(k) \
  active##k = (bpm_v2u_t){(top_level_a>k) ? BMP_W64_ONES : 0,(top_level_b>k) ? BMP_W64_ONES : 0}; \
  last##k = (bpm_v2u_t){(top_level_a==k+1) ? BMP_W64_ONES : 0,(top_level_b==k+1) ? BMP_W64_ONES : 0};
#define BPM_V2_TOP_LEVELS_CHANGED() { \
  mixed = (top_level_a!=top_level_b); \
  BPM_V2_LANE_MASKS(0) BPM_V2_LANE_MASKS(1) BPM_V2_LANE_MASKS(2) BPM_V2_LANE_MASKS(3) \
  /* A lane leaves (max_distance+1,max_distance+init_score[top_level-1]] only at BPM_CUT_OFF's branches */ \
  range_a = init_score[top_level_a-1]-1; \
  range_b = init_score[top_level_b-1]-1; \
  if (top_level_a > written_a) written_a = top_level_a; \
  if (top_level_b > written_b) written_b = top_level_b; \
}
// Both lanes' column: A's text[text_position], B's text[b_position]; leaves score_a/score_b
#define BPM_V2_COLUMN(b_position) \
  const uint8_t enc_a = text[text_position]; \
  const uint8_t enc_b = text[b_position]; \
  const uint64_t pair = enc_a*DNA__N_RANGE+enc_b; \
  bpm_v2u_t PHin = zero, MHin = zero; \
  bpm_v2i_t Ssel; \
  if (top_level_a > NW || top_level_b > NW) __builtin_unreachable(); \
  if (gem_expect_true(!mixed)) { \
    switch (top_level_a) { \
      case 1: BPM_V2_STEP(0); Ssel = S0; break; \
      case 2: BPM_V2_STEP(0); BPM_V2_STEP(1); Ssel = S1; break; \
      case 3: BPM_V2_STEP(0); BPM_V2_STEP(1); BPM_V2_STEP(2); Ssel = S2; break; \
      default: BPM_V2_STEP(0); BPM_V2_STEP(1); BPM_V2_STEP(2); BPM_V2_STEP(3); Ssel = S3; break; \
    } \
  } else { \
    bpm_v2u_t PHsel = zero, MHsel = zero; \
    Ssel = (bpm_v2i_t)zero; \
    BPM_V2_STEP_MIXED(0); \
    if (NW > 1) BPM_V2_STEP_MIXED(1); \
    if (NW > 2) BPM_V2_STEP_MIXED(2); \
    if (NW > 3) BPM_V2_STEP_MIXED(3); \
    PHin = PHsel; MHin = MHsel; \
  } \
  int64_t score_a = Ssel[0], score_b = Ssel[1]; \
  if (gem_expect_false(((uint64_t)(score_a-(int64_t)max_distance__1-1) >= range_a) | \
                       ((uint64_t)(score_b-(int64_t)max_distance__1-1) >= range_b))) { \
    const uint8_t top_level_a_old = top_level_a, top_level_b_old = top_level_b; \
    BPM_V2_CUT_OFF(0,enc_a,top_level_a,score_a); \
    BPM_V2_CUT_OFF(1,enc_b,top_level_b,score_b); \
    if (top_level_a!=top_level_a_old || top_level_b!=top_level_b_old) BPM_V2_TOP_LEVELS_CHANGED(); \
  }
// A's match logic in the paired loop; a quick abandon leaves A's lane in the pattern
#define BPM_V2_MATCH_A() \
  BPM_ALL_CHECK_MATCH(top_level_a,score_a); \
  if (BPM_ALL_QUICK_ABANDON(top_level_a,score_a)) { \
    PROF_INC_COUNTER(GP_BPM_ALL_QUICK_ABANDON); \
    BPM_V2_STORE_LANE(0,written_a,P,M,score); \
    goto bpm_all_done; /* Quick abandon, it doesn't match (bounded by best case scenario) */ \
  }
#define BPM_V2_STORE_WORD(k,L,num_words,Pd,Md,Sd) \
  if (k < num_words) { Pd[k] = P##k[L]; Md[k] = M##k[L]; Sd[k] = S##k[L]; }
#define BPM_V2_STORE_LANE(L,num_words,Pd,Md,Sd) { \
  BPM_V2_STORE_WORD(0,L,num_words,Pd,Md,Sd) \
  if (NW > 1) BPM_V2_STORE_WORD(1,L,num_words,Pd,Md,Sd) \
  if (NW > 2) BPM_V2_STORE_WORD(2,L,num_words,Pd,Md,Sd) \
  if (NW > 3) BPM_V2_STORE_WORD(3,L,num_words,Pd,Md,Sd) \
}
static inline __attribute__((always_inline)) uint64_t bpm_compute_edit_distance_all_v2(
    pattern_t* const pattern,
    filtering_candidates_t* const filtering_candidates,
    const uint64_t begin_position,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t max_distance,
    const uint64_t NW) {
  // Parameters
  bpm_pattern_t* const bpm_pattern = &pattern->pattern_tiled.bpm_pattern;
  // Pattern variables
  const uint64_t* PEQ = bpm_pattern->PEQ;
  const uint64_t num_words64 = NW;
  const uint64_t key_length = bpm_pattern->pattern_length;
  uint64_t* const P = bpm_pattern->P;
  uint64_t* const M = bpm_pattern->M;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  int64_t* const score = bpm_pattern->score;
  const int64_t* const init_score = bpm_pattern->init_score;
  const uint64_t* const pattern_left = bpm_pattern->pattern_left;
  // Initialize search
  const uint64_t max_distance__1 = max_distance+1;
  const uint8_t top = num_words64-1;
  uint8_t top_level;
  uint64_t min_score = ALIGN_DISTANCE_INF, min_score_column = ALIGN_COLUMN_INF, opt_steps_left = 0;
  bpm_reset_search_cutoff(&top_level,P,M,score,init_score,max_distance);
  const uint64_t end_position = begin_position + text_length;
  bool match_found = false;
  uint64_t text_position = 0, text_left = text_length, num_matches_found = 0;
  // Split
  const uint64_t X = BPM_ALL_SPLIT_OVERLAP;
  const uint64_t s = (text_length-X+1)/2; // B's first column
  const uint64_t num_records = text_length-s-X; // B's columns after s+X-1
  // Pair tables {PEQ of A's character, PEQ of B's character} and the lanes' state
  const bpm_v2u_t zero = {0,0}, one = {1,1};
  bpm_v2u_t E0[DNA__N_RANGE*DNA__N_RANGE], E1[DNA__N_RANGE*DNA__N_RANGE];
  bpm_v2u_t E2[DNA__N_RANGE*DNA__N_RANGE], E3[DNA__N_RANGE*DNA__N_RANGE];
  bpm_v2u_t* const E[4] = {E0,E1,E2,E3};
  uint64_t shifts[4] = {0,0,0,0};
  uint64_t k, a, b;
  for (k=0;k<4;++k) { // Tables past NW are never read
    if (k<NW) shifts[k] = __builtin_ctzll(level_mask[k]);
    for (a=0;a<DNA__N_RANGE;++a) {
      for (b=0;b<DNA__N_RANGE;++b) {
        E[k][a*DNA__N_RANGE+b] = (k<NW) ?
            (bpm_v2u_t){PEQ[BPM_PATTERN_PEQ_IDX(k,a)],PEQ[BPM_PATTERN_PEQ_IDX(k,b)]} : zero;
      }
    }
  }
  const uint64_t shift0 = shifts[0], shift1 = shifts[1], shift2 = shifts[2], shift3 = shifts[3];
  bpm_v2u_t P0 = zero, P1 = zero, P2 = zero, P3 = zero;
  bpm_v2u_t M0 = zero, M1 = zero, M2 = zero, M3 = zero;
  bpm_v2i_t S0 = (bpm_v2i_t)zero, S1 = (bpm_v2i_t)zero, S2 = (bpm_v2i_t)zero, S3 = (bpm_v2i_t)zero;
  bpm_v2u_t active0, active1, active2, active3, last0, last1, last2, last3;
  if (top_level > 0) { P0 = (bpm_v2u_t){P[0],P[0]}; M0 = (bpm_v2u_t){M[0],M[0]}; S0 = (bpm_v2i_t){score[0],score[0]}; }
  if (top_level > 1) { P1 = (bpm_v2u_t){P[1],P[1]}; M1 = (bpm_v2u_t){M[1],M[1]}; S1 = (bpm_v2i_t){score[1],score[1]}; }
  if (top_level > 2) { P2 = (bpm_v2u_t){P[2],P[2]}; M2 = (bpm_v2u_t){M[2],M[2]}; S2 = (bpm_v2i_t){score[2],score[2]}; }
  if (top_level > 3) { P3 = (bpm_v2u_t){P[3],P[3]}; M3 = (bpm_v2u_t){M[3],M[3]}; S3 = (bpm_v2i_t){score[3],score[3]}; }
  uint8_t top_level_a = top_level, top_level_b = top_level;
  uint8_t written_a = top_level, written_b = top_level;
  uint64_t range_a, range_b;
  bool mixed;
  BPM_V2_TOP_LEVELS_CHANGED();
  // A on [0,X), B on [s,s+X)
  for (;text_position<X;++text_position,--text_left) {
    BPM_V2_COLUMN(s+text_position);
    BPM_V2_MATCH_A();
  }
  // B's checkpoint after column s+X-1
  uint64_t PC[BPM_ALL_SPLIT_MAX_WORDS] = {0}, MC[BPM_ALL_SPLIT_MAX_WORDS] = {0};
  int64_t SC[BPM_ALL_SPLIT_MAX_WORDS] = {0};
  const uint8_t top_level_c = top_level_b;
  BPM_V2_STORE_LANE(1,NW,PC,MC,SC);
  // A on [X,X+num_records), B on [s+X,text_length), recorded
  uint64_t PA[BPM_ALL_SPLIT_MAX_WORDS] = {0}, MA[BPM_ALL_SPLIT_MAX_WORDS] = {0}; // B's state at A's abandon column
  int64_t SA[BPM_ALL_SPLIT_MAX_WORDS] = {0};
  uint8_t record_top[BPM_ALL_SPLIT_MAX_RECORDS];
  int64_t record_score[BPM_ALL_SPLIT_MAX_RECORDS];
  uint8_t written_b_abandon = 0;
  bool b_settled = false;
  written_b = top_level_b;
  const uint64_t paired_end = X+num_records;
  for (;text_position<paired_end;++text_position,--text_left) {
    const uint64_t b_position = s+text_position;
    BPM_V2_COLUMN(b_position);
    BPM_V2_MATCH_A();
    // Record B's column
    const uint64_t r = text_position-X;
    record_top[r] = top_level_b;
    record_score[r] = score_b;
    if (!b_settled) {
      const uint64_t b_text_left = text_length-b_position;
      if (top_level_b==num_words64 && score_b<=max_distance) {
        b_settled = true; // A match: no quick abandon from here on
      } else if (score_b+pattern_left[top_level_b] > b_text_left+max_distance) {
        b_settled = true;
        written_b_abandon = written_b;
        BPM_V2_STORE_LANE(1,written_b,PA,MA,SA);
      }
    }
  }
  // A's lane to the pattern; A alone up to s+X-1 (at most one column)
  uint64_t PB[BPM_ALL_SPLIT_MAX_WORDS] = {0}, MB[BPM_ALL_SPLIT_MAX_WORDS] = {0};
  int64_t SB[BPM_ALL_SPLIT_MAX_WORDS] = {0};
  BPM_V2_STORE_LANE(0,written_a,P,M,score);
  BPM_V2_STORE_LANE(1,written_b,PB,MB,SB);
  top_level = top_level_a;
  for (;text_position<s+X;++text_position,--text_left) {
    BPM_ALL_COLUMN_A();
  }
  // Same state after column s+X-1: B's columns are the search's
  bool same = (top_level==top_level_c);
  for (k=0;same && k<top_level_c;++k) {
    same = (P[k]==PC[k] && M[k]==MC[k] && score[k]==SC[k]);
  }
  if (same) {
    uint64_t r;
    for (r=0;r<num_records;++r,++text_position,--text_left) {
      const uint8_t top_level_r = record_top[r];
      const int64_t current_score = record_score[r];
      BPM_ALL_CHECK_MATCH(top_level_r,current_score);
      if (BPM_ALL_QUICK_ABANDON(top_level_r,current_score)) {
        PROF_INC_COUNTER(GP_BPM_ALL_QUICK_ABANDON);
        for (k=0;k<written_b_abandon;++k) { P[k]=PA[k]; M[k]=MA[k]; score[k]=SA[k]; }
        goto bpm_all_done; // Quick abandon, it doesn't match (bounded by best case scenario)
      }
    }
    for (k=0;k<written_b;++k) { P[k]=PB[k]; M[k]=MB[k]; score[k]=SB[k]; }
    goto bpm_all_done;
  }
  for (;text_position<text_length;++text_position,--text_left) {
    BPM_ALL_COLUMN_A();
  }
bpm_all_done:
  if (match_found) {
    const uint64_t text_end_offset = min_score_column+1;
    const uint64_t text_begin_offset = BOUNDED_SUBTRACTION(text_end_offset,key_length+min_score,0);
    filtering_candidates_add_region_verified(
        filtering_candidates,pattern,text_begin_offset,
        text_end_offset,begin_position,end_position,min_score);
    ++num_matches_found; // Increment the number of matches found
  }
  return num_matches_found;
}
#define BPM_ALL_V2_KERNEL(NW) \
static __attribute__((noinline)) uint64_t bpm_compute_edit_distance_all_v2_w##NW( \
    pattern_t* const pattern,filtering_candidates_t* const filtering_candidates, \
    const uint64_t begin_position,const uint8_t* const text,const uint64_t text_length, \
    const uint64_t max_distance) { \
  return bpm_compute_edit_distance_all_v2(pattern,filtering_candidates, \
      begin_position,text,text_length,max_distance,NW); \
}
BPM_ALL_V2_KERNEL(1)
BPM_ALL_V2_KERNEL(2)
BPM_ALL_V2_KERNEL(3)
BPM_ALL_V2_KERNEL(4)
// Windows the lanes take: 1-4 words of single-bit level masks, text bytes 0-4. A byte above 4
// indexes PEQ past the word's 5 rows (the next word's, or P/M), which the pair tables do not hold
static bool bpm_compute_edit_distance_all_split(
    const bpm_pattern_t* const bpm_pattern,
    const uint8_t* const text,
    const uint64_t text_length) {
  const uint64_t num_words64 = bpm_pattern->pattern_num_words64;
  if (text_length < BPM_ALL_SPLIT_MIN_LENGTH) return false;
  if (num_words64 < 1 || num_words64 > BPM_ALL_SPLIT_MAX_WORDS) return false;
  if (text_length/2 > BPM_ALL_SPLIT_MAX_RECORDS) return false;
  uint64_t i;
  for (i=0;i<num_words64;++i) {
    const uint64_t mask = bpm_pattern->level_mask[i];
    if (mask==0 || (mask & (mask-1))!=0) return false;
  }
  uint8_t above = 0;
  for (i=0;i<text_length;++i) above |= (text[i] > ENC_DNA_CHAR_N);
  return above==0;
}
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
  if (bpm_compute_edit_distance_all_split(bpm_pattern,text,text_length)) {
    uint64_t num_matches_found;
    switch (num_words64) {
      case 1: num_matches_found = bpm_compute_edit_distance_all_v2_w1(pattern,
          filtering_candidates,begin_position,text,text_length,max_distance); break;
      case 2: num_matches_found = bpm_compute_edit_distance_all_v2_w2(pattern,
          filtering_candidates,begin_position,text,text_length,max_distance); break;
      case 3: num_matches_found = bpm_compute_edit_distance_all_v2_w3(pattern,
          filtering_candidates,begin_position,text,text_length,max_distance); break;
      default: num_matches_found = bpm_compute_edit_distance_all_v2_w4(pattern,
          filtering_candidates,begin_position,text,text_length,max_distance); break;
    }
    PROF_INC_COUNTER(GP_BPM_ALL_MATCHES_FOUND);
    PROF_STOP(GP_BPM_ALL);
    return num_matches_found;
  }
  if (num_words64 >= 1 && num_words64 <= BPM_W4_MAX_WORDS && bpm_text_within_peq(text,text_length)) {
    uint64_t num_matches_found;
    switch (num_words64) {
      case 1: num_matches_found = bpm_compute_edit_distance_all_w1(pattern,
          filtering_candidates,begin_position,text,text_length,max_distance); break;
      case 2: num_matches_found = bpm_compute_edit_distance_all_w2(pattern,
          filtering_candidates,begin_position,text,text_length,max_distance); break;
      case 3: num_matches_found = bpm_compute_edit_distance_all_w3(pattern,
          filtering_candidates,begin_position,text,text_length,max_distance); break;
      default: num_matches_found = bpm_compute_edit_distance_all_w4(pattern,
          filtering_candidates,begin_position,text,text_length,max_distance); break;
    }
    PROF_INC_COUNTER(GP_BPM_ALL_MATCHES_FOUND);
    PROF_STOP(GP_BPM_ALL);
    return num_matches_found;
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
