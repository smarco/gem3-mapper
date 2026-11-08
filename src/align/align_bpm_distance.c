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
/*
 * BPM all matches in BPM_ALL_V4_LANES segments
 *   A longer window runs as BPM_ALL_V4_LANES searches of T steps each in the lanes
 *   of 64-bit vectors: lane 0 continues the search from column e (the first
 *   e < BPM_ALL_V4_LANES columns run alone) and lane j>0 starts from a fresh state
 *   at s_j = e+j*(T-X). As in the split above, lane j's columns after one where its
 *   state equals the search's are the search's. Lane j is compared after column
 *   s_j+X-1; when that fails, the search runs on alone and compares again at
 *   lane j's snapshots (every BPM_ALL_V4_SNAP steps), and then at the next lane's.
 *   The match logic then runs on the lanes' recorded (top_level,score) in column
 *   order, visiting only the columns where it can act: a score within max_distance
 *   at the last word, a quick abandon bound crossed, or a pending match.
 *   Built where 4 x 64-bit vectors are native (AVX2); elsewhere the split above runs.
 */
#if defined(__AVX2__)
#include <immintrin.h>
#define BPM_ALL_V4_LANES      4
#define BPM_ALL_V4_OVERLAP    192
#define BPM_ALL_V4_SNAP       32
#define BPM_ALL_V4_MAX_STEPS  1024
#define BPM_ALL_V4_MAX_SNAPS  (BPM_ALL_V4_MAX_STEPS/BPM_ALL_V4_SNAP)
#define BPM_ALL_V4_MIN_LENGTH (2*BPM_ALL_V4_OVERLAP)
typedef uint64_t bpm_v4u_t __attribute__((vector_size(32)));
typedef int64_t bpm_v4i_t __attribute__((vector_size(32)));
typedef uint32_t bpm_v4tl_t; // Every lane's top_level, a byte each
#define BPM_V4_FOR_LANES(F) F(0) F(1) F(2) F(3)
#define BPM_V4_LANEWISE(F,k) {F(0,k),F(1,k),F(2,k),F(3,k)}
#define BPM_V4_TOP(packed,j) ((uint8_t)((packed)>>(8*(j))))
#define BPM_V4_GT(j,k)   ((tl[j]>k) ? BMP_W64_ONES : 0)
#define BPM_V4_EQ(j,k)   ((tl[j]==k) ? BMP_W64_ONES : 0)
#define BPM_V4_DEC(j,_)  ((int64_t)(max_distance+init_score[tl[j]-1]))
#define BPM_V4_AB(j,_)   ((int64_t)(text_length-s[j]+max_distance-pattern_left[tl[j]])-(int64_t)t)
#define BPM_V4_SC(j,_)   sc[j]
// Lanes whose top_level is above k, or equal to k
#define BPM_V4_ABOVE(k)  ((bpm_v4u_t)BPM_V4_LANEWISE(BPM_V4_GT,k))
#define BPM_V4_AT(k)     ((bpm_v4u_t)BPM_V4_LANEWISE(BPM_V4_EQ,k))
// Eq of the lanes from the pair tables of lanes (0,1) and (2,3)
#if defined(__clang__) || __GNUC__ >= 12
#define BPM_V4_JOIN(E) __builtin_shufflevector(E[pr[0]],E[pr[1]],0,1,2,3)
#else
#define BPM_V4_JOIN(E) ((bpm_v4u_t){E[pr[0]][0],E[pr[0]][1],E[pr[1]][0],E[pr[1]][1]})
#endif
// Bit j: lane j's sign bit
#define BPM_V4_SIGN_BITS(v) ((uint64_t)_mm256_movemask_pd(_mm256_castsi256_pd((__m256i)(v))))
#define BPM_V4_BLOCK(k) \
  const bpm_v4u_t Xv = Eq##k | M##k; \
  const bpm_v4u_t _Eq = Eq##k | MHin; \
  const bpm_v4u_t Xh = (((_Eq & P##k) + P##k) ^ P##k) | _Eq; \
  bpm_v4u_t Ph = M##k | ~(Xh | P##k); \
  bpm_v4u_t Mh = P##k & Xh; \
  const bpm_v4u_t PHout = (k+1 < NW) ? Ph >> 63 : (Ph >> shift##k) & one; /* words below the last: bit 63 */ \
  const bpm_v4u_t MHout = (k+1 < NW) ? Mh >> 63 : (Mh >> shift##k) & one; \
  Ph = (Ph << 1) | PHin; \
  Mh = (Mh << 1) | MHin; \
  const bpm_v4u_t Pn = Mh | ~(Xv | Ph); \
  const bpm_v4u_t Mn = Ph & Xv; \
  const bpm_v4i_t Sn = S##k + (bpm_v4i_t)(PHout-MHout);
/* Every lane advances the words below the highest top_level; a word above a lane's top_level
   runs on with values nothing reads, and the lane's values from when it dropped the word are
   kept in stale_P/M/S, put back at the snapshots and at the end */
#define BPM_V4_STEP(k) { \
  BPM_V4_BLOCK(k) \
  P##k = Pn; M##k = Mn; S##k = Sn; \
  PHin = PHout; MHin = MHout; \
  PHo##k = PHout; MHo##k = MHout; \
}
// The carries and score of each lane's top word
#define BPM_V4_SELECT(k) \
  if (NW > k && hi > k) { \
    PHsel = (PHo##k & at##k) | (PHsel & ~at##k); \
    MHsel = (MHo##k & at##k) | (MHsel & ~at##k); \
    Ssel = (S##k & (bpm_v4i_t)at##k) | (Ssel & ~(bpm_v4i_t)at##k); \
  }
#define BPM_V4_SET_STALE(k,L) { stale_P[k][L] = P##k[L]; stale_M[k][L] = M##k[L]; stale_S[k][L] = S##k[L]; }
#define BPM_V4_UNSTALE(k) \
  if (NW > k) { \
    bpm_v4u_t sP, sM; bpm_v4i_t sS; \
    memcpy(&sP,stale_P[k],sizeof(sP)); memcpy(&sM,stale_M[k],sizeof(sM)); memcpy(&sS,stale_S[k],sizeof(sS)); \
    kP[k] = (P##k & above##k) | (sP & ~above##k); \
    kM[k] = (M##k & above##k) | (sM & ~above##k); \
    kS[k] = (S##k & (bpm_v4i_t)above##k) | (sS & ~(bpm_v4i_t)above##k); \
  }
#define BPM_V4_SET_WORD(k,L,Pv,Mv,Sv) { P##k[L] = Pv; M##k[L] = Mv; S##k[L] = Sv; }
// BPM_CUT_OFF on lane L (scalar) when its cut test is set; leaves the lane's score[top_level-1] in sc[L]
#define BPM_V4_CUT_OFF(L) { \
  int64_t current_score = Ssel[L]; \
  if ((cut_lanes>>L) & 1) { \
    uint8_t top_level_l = tl[L]; \
    const uint8_t enc_char = ch[L]; \
    const uint8_t last = top_level_l-1; \
    const uint64_t PHl = PHsel[L], MHl = MHsel[L]; \
    if (current_score<=max_distance__1) { \
      const uint64_t last_score = current_score+(MHl-PHl); \
      if (last_score<=max_distance && last<top && (MHl || (PEQ[BPM_PATTERN_PEQ_IDX(top_level_l,enc_char)] & 1))) { \
        const uint64_t Peq = PEQ[BPM_PATTERN_PEQ_IDX(top_level_l,enc_char)]; \
        uint64_t Pv = BMP_W64_ONES; \
        uint64_t Mv = 0; \
        uint64_t PHout, MHout; \
        const uint64_t mask = level_mask[top_level_l]; \
        BPM_ADVANCE_BLOCK(Peq,mask,Pv,Mv,PHl,MHl,PHout,MHout); \
        const int64_t score_top = last_score + init_score[top_level_l] + (PHout-MHout); \
        switch (top_level_l) { \
          case 1: BPM_V4_SET_WORD(1,L,Pv,Mv,score_top); break; \
          case 2: BPM_V4_SET_WORD(2,L,Pv,Mv,score_top); break; \
          default: BPM_V4_SET_WORD(3,L,Pv,Mv,score_top); break; \
        } \
        ++top_level_l; \
        current_score = score_top; \
      } /* else BPM_CUT_OFF's while loop cannot fire: score <= max_distance+1 <= max_distance+init_score */ \
    } else { \
      while (current_score > (max_distance+init_score[top_level_l-1])) { \
        --top_level_l; \
        switch (top_level_l) { \
          case 1: BPM_V4_SET_STALE(1,L); current_score = S0[L]; break; \
          case 2: BPM_V4_SET_STALE(2,L); current_score = S1[L]; break; \
          default: BPM_V4_SET_STALE(3,L); current_score = S2[L]; break; \
        } \
      } \
    } \
    if (top_level_l!=tl[L]) { tl[L] = top_level_l; changed = true; } \
  } \
  sc[L] = current_score; \
}
// Lane L's column can act in the match logic: a score within max_distance at the last word, or past the quick abandon bound
#define BPM_V4_EVENT(L) \
  events |= (uint8_t)(((tl[L]==num_words64 && sc[L]<=(int64_t)max_distance) || \
      (sc[L]+(int64_t)pattern_left[tl[L]] > (int64_t)(text_length-s[L]-t+max_distance))) << L);
#define BPM_V4_TOP_LEVELS_CHANGED() { \
  uint64_t j; \
  lo = tl[0]; hi = tl[0]; \
  for (j=1;j<BPM_ALL_V4_LANES;++j) { \
    if (tl[j] < lo) lo = tl[j]; \
    if (tl[j] > hi) hi = tl[j]; \
  } \
  mixed = (lo!=hi); \
  above1 = BPM_V4_ABOVE(1); above2 = BPM_V4_ABOVE(2); above3 = BPM_V4_ABOVE(3); \
  at0 = BPM_V4_AT(1); at1 = BPM_V4_AT(2); at2 = BPM_V4_AT(3); at3 = BPM_V4_AT(4); \
  below_nw = (bpm_v4i_t)~BPM_V4_AT(NW); \
  at_nw = (bpm_v4i_t)BPM_V4_AT(NW); \
  dec_bound = (bpm_v4i_t)BPM_V4_LANEWISE(BPM_V4_DEC,0); \
  abandon_bound = (bpm_v4i_t)BPM_V4_LANEWISE(BPM_V4_AB,0); \
  tl_packed = 0; \
  for (j=0;j<BPM_ALL_V4_LANES;++j) { \
    tl_packed |= (bpm_v4tl_t)tl[j] << (8*j); \
    if (tl[j] > written[j]) written[j] = tl[j]; \
  } \
}
#define BPM_V4_KEEP(dP,dM,dS) { \
  bpm_v4u_t* const kP = dP; bpm_v4u_t* const kM = dM; bpm_v4i_t* const kS = dS; \
  BPM_V4_UNSTALE(0) BPM_V4_UNSTALE(1) BPM_V4_UNSTALE(2) BPM_V4_UNSTALE(3) \
}
// Lane m's state after step t_end, from its last snapshot (lane 0 before its first: the pattern); then quick abandon
#define BPM_V4_ABANDON(t_end,written_m) { \
  PROF_INC_COUNTER(GP_BPM_ALL_QUICK_ABANDON); \
  const uint64_t n_snap = (t_end+1)/BPM_ALL_V4_SNAP; \
  uint64_t c; \
  if (n_snap==0) { \
    for (c=s[0];c<=s[0]+t_end;++c) { \
      uint64_t PHin=0,MHin=0,PHout,MHout; \
      const uint8_t enc_char = text[c]; \
      BPM_ADVANCE_COLUMN(P,M,score,PEQ,enc_char,top_level,PHin,MHin,PHout,MHout); \
      BPM_CUT_OFF(P,M,score,PEQ,enc_char,top_level,PHin,MHin,PHout,MHout); \
    } \
  } else { \
    uint64_t rP[4], rM[4]; \
    int64_t rS[4]; \
    uint8_t top_level_m = BPM_V4_TOP(snap_top[n_snap-1],m); \
    for (k=0;k<NW;++k) { rP[k] = snap_P[n_snap-1][k][m]; rM[k] = snap_M[n_snap-1][k][m]; rS[k] = snap_S[n_snap-1][k][m]; } \
    for (c=s[m]+n_snap*BPM_ALL_V4_SNAP;c<=s[m]+t_end;++c) { \
      uint64_t PHin=0,MHin=0,PHout,MHout; \
      const uint8_t enc_char = text[c]; \
      BPM_ADVANCE_COLUMN(rP,rM,rS,PEQ,enc_char,top_level_m,PHin,MHin,PHout,MHout); \
      BPM_CUT_OFF(rP,rM,rS,PEQ,enc_char,top_level_m,PHin,MHin,PHout,MHout); \
    } \
    for (k=0;k<written_m && k<NW;++k) { P[k]=rP[k]; M[k]=rM[k]; score[k]=rS[k]; } \
  } \
  goto bpm_all_done; /* Quick abandon, it doesn't match (bounded by best case scenario) */ \
}
static inline __attribute__((always_inline)) uint64_t bpm_compute_edit_distance_all_v4(
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
  const uint8_t y = top_level;
  const uint64_t end_position = begin_position + text_length;
  bool match_found = false;
  uint64_t text_position = 0, text_left = text_length, num_matches_found = 0;
  uint64_t j, k, a, b, t = 0;
  // Segments: T steps per lane, lane j from column s[j]
  const uint64_t X = BPM_ALL_V4_OVERLAP;
  const uint64_t e = (text_length+(BPM_ALL_V4_LANES-1)*X) % BPM_ALL_V4_LANES;
  const uint64_t T = (text_length-e+(BPM_ALL_V4_LANES-1)*X)/BPM_ALL_V4_LANES;
  uint64_t s[BPM_ALL_V4_LANES];
  const uint8_t* lane_text[BPM_ALL_V4_LANES];
  for (j=0;j<BPM_ALL_V4_LANES;++j) { s[j] = e+j*(T-X); lane_text[j] = text+s[j]; }
  // Lane 0's first e columns alone
  uint8_t written_a = top_level;
  for (;text_position<e;++text_position,--text_left) {
    BPM_ALL_COLUMN_A();
    if (top_level > written_a) written_a = top_level;
  }
  // Pair tables {PEQ of lane 2i's character, PEQ of lane 2i+1's character}
  const bpm_v4u_t zero = {0}, one = zero+1;
  const bpm_v4i_t md1 = (bpm_v4i_t)zero+(int64_t)max_distance__1;
  bpm_v2u_t E0[DNA__N_RANGE*DNA__N_RANGE], E1[DNA__N_RANGE*DNA__N_RANGE];
  bpm_v2u_t E2[DNA__N_RANGE*DNA__N_RANGE], E3[DNA__N_RANGE*DNA__N_RANGE];
  bpm_v2u_t* const E[4] = {E0,E1,E2,E3};
  uint64_t shifts[4] = {0,0,0,0};
  for (k=0;k<4;++k) { // Tables past NW are never read
    if (k<NW) shifts[k] = __builtin_ctzll(level_mask[k]);
    for (a=0;a<DNA__N_RANGE;++a) {
      for (b=0;b<DNA__N_RANGE;++b) {
        E[k][a*DNA__N_RANGE+b] = (k<NW) ?
            (bpm_v2u_t){PEQ[BPM_PATTERN_PEQ_IDX(k,a)],PEQ[BPM_PATTERN_PEQ_IDX(k,b)]} : (bpm_v2u_t){0,0};
      }
    }
  }
  const uint64_t shift0 = shifts[0], shift1 = shifts[1], shift2 = shifts[2], shift3 = shifts[3];
  // Lane 0 from the pattern (words written so far), the others from bpm_reset_search_cutoff
  bpm_v4u_t P0 = zero, P1 = zero, P2 = zero, P3 = zero;
  bpm_v4u_t M0 = zero, M1 = zero, M2 = zero, M3 = zero;
  bpm_v4i_t S0 = (bpm_v4i_t)zero, S1 = (bpm_v4i_t)zero, S2 = (bpm_v4i_t)zero, S3 = (bpm_v4i_t)zero;
  {
    int64_t reset_score = 0;
#define BPM_V4_INIT_WORD(k) \
    if (k < NW) { \
      reset_score += init_score[k]; \
      if (k < y) { P##k = zero+BMP_W64_ONES; S##k = (bpm_v4i_t)zero+reset_score; } \
      if (k < written_a) { P##k[0] = P[k]; M##k[0] = M[k]; S##k[0] = score[k]; } \
      else { P##k[0] = 0; M##k[0] = 0; S##k[0] = 0; } \
    }
    BPM_V4_INIT_WORD(0) BPM_V4_INIT_WORD(1) BPM_V4_INIT_WORD(2) BPM_V4_INIT_WORD(3)
  }
  // Per lane: top_level, and the highest since the column it is compared at (written back from there)
  uint8_t tl[BPM_ALL_V4_LANES], written[BPM_ALL_V4_LANES];
  tl[0] = top_level; written[0] = written_a;
  for (j=1;j<BPM_ALL_V4_LANES;++j) { tl[j] = y; written[j] = y; }
  const bpm_v4u_t above0 = ~zero;
  bpm_v4u_t above1, above2, above3, at0, at1, at2, at3;
  // Lane 0's words from the pattern above its top_level are stale already
  uint64_t stale_P[4][BPM_ALL_V4_LANES], stale_M[4][BPM_ALL_V4_LANES];
  int64_t stale_S[4][BPM_ALL_V4_LANES];
  for (k=0;k<4;++k) {
    for (j=0;j<BPM_ALL_V4_LANES;++j) { stale_P[k][j] = 0; stale_M[k][j] = 0; stale_S[k][j] = 0; }
    if (k < written_a) { stale_P[k][0] = P[k]; stale_M[k][0] = M[k]; stale_S[k][0] = score[k]; }
  }
  bpm_v4i_t below_nw, at_nw, dec_bound, abandon_bound;
  bpm_v4tl_t tl_packed;
  uint8_t lo, hi;
  bool mixed;
  BPM_V4_TOP_LEVELS_CHANGED();
  // Every lane's (top_level,score) per step, its state every BPM_ALL_V4_SNAP steps, and its match-logic columns
  bpm_v4i_t rec_score[BPM_ALL_V4_MAX_STEPS];
  bpm_v4tl_t rec_top[BPM_ALL_V4_MAX_STEPS];
  bpm_v4u_t snap_P[BPM_ALL_V4_MAX_SNAPS][4], snap_M[BPM_ALL_V4_MAX_SNAPS][4];
  bpm_v4i_t snap_S[BPM_ALL_V4_MAX_SNAPS][4];
  bpm_v4tl_t snap_top[BPM_ALL_V4_MAX_SNAPS];
  uint8_t rec_event[BPM_ALL_V4_MAX_STEPS+8]; // Bit j: lane j's column can act in the match logic
  for (t=0;t<T;++t) {
    uint8_t ch[BPM_ALL_V4_LANES];
    uint64_t pr[BPM_ALL_V4_LANES/2];
    for (j=0;j<BPM_ALL_V4_LANES;++j) ch[j] = lane_text[j][t];
    for (j=0;j<BPM_ALL_V4_LANES/2;++j) pr[j] = ch[2*j]*DNA__N_RANGE+ch[2*j+1];
    // Eq of the words advanced, and of the word above each lane's top
    const bpm_v4u_t Eq0 = BPM_V4_JOIN(E0);
    const bpm_v4u_t Eq1 = (NW > 1) ? BPM_V4_JOIN(E1) : zero;
    bpm_v4u_t Eq2 = zero, Eq3 = zero;
    if (NW > 2 && hi >= 2) Eq2 = BPM_V4_JOIN(E2);
    if (NW > 3 && hi >= 3) Eq3 = BPM_V4_JOIN(E3);
    bpm_v4u_t PHin = zero, MHin = zero, PHsel, MHsel, Eq_next;
    bpm_v4u_t PHo0 = zero, MHo0 = zero, PHo1 = zero, MHo1 = zero, PHo2 = zero, MHo2 = zero, PHo3 = zero, MHo3 = zero;
    bpm_v4i_t Ssel;
    if (hi > NW || lo < 1) __builtin_unreachable();
    switch (hi) {
      case 1: BPM_V4_STEP(0); Ssel = S0; Eq_next = Eq1; break;
      case 2: BPM_V4_STEP(0); BPM_V4_STEP(1); Ssel = S1; Eq_next = Eq2; break;
      case 3: BPM_V4_STEP(0); BPM_V4_STEP(1); BPM_V4_STEP(2); Ssel = S2; Eq_next = Eq3; break;
      default: BPM_V4_STEP(0); BPM_V4_STEP(1); BPM_V4_STEP(2); BPM_V4_STEP(3); Ssel = S3; Eq_next = zero; break;
    }
    PHsel = PHin; MHsel = MHin;
    if (gem_expect_false(mixed)) {
      Ssel = S0; PHsel = PHo0; MHsel = MHo0;
      BPM_V4_SELECT(1) BPM_V4_SELECT(2) BPM_V4_SELECT(3)
      Eq_next = (Eq1 & at0) | (Eq2 & at1) | (Eq3 & at2);
    }
    // BPM_CUT_OFF acts on a lane (sign bits; the signed tests include the unsigned ones): it extends
    // (score+MHin-PHin <= max_distance, which implies score <= max_distance+1), or drops a level
    const bpm_v4i_t Sm = Ssel-md1;
    const bpm_v4i_t extend = (Sm+(bpm_v4i_t)(MHsel-PHsel)) & (bpm_v4i_t)((MHsel|Eq_next) << 63) & below_nw;
    const bpm_v4i_t cut = extend | (dec_bound-Ssel);
    // The match logic can act on the lane's column
    const bpm_v4i_t event = (Sm & at_nw) | (abandon_bound-Ssel);
    const uint64_t cut_lanes = BPM_V4_SIGN_BITS(cut);
    rec_event[t] = (uint8_t)BPM_V4_SIGN_BITS(event);
    if (gem_expect_false(cut_lanes!=0)) {
      int64_t sc[BPM_ALL_V4_LANES];
      bool changed = false;
      BPM_V4_FOR_LANES(BPM_V4_CUT_OFF)
      if (changed) {
        BPM_V4_TOP_LEVELS_CHANGED();
        Ssel = (bpm_v4i_t)BPM_V4_LANEWISE(BPM_V4_SC,0);
        uint8_t events = 0;
        BPM_V4_FOR_LANES(BPM_V4_EVENT)
        rec_event[t] = events;
      }
    }
    rec_score[t] = Ssel;
    rec_top[t] = tl_packed;
    if (((t+1)%BPM_ALL_V4_SNAP)==0) {
      const uint64_t n = t/BPM_ALL_V4_SNAP;
      BPM_V4_KEEP(snap_P[n],snap_M[n],snap_S[n])
      snap_top[n] = tl_packed;
      if (t==X-1) {
        for (j=1;j<BPM_ALL_V4_LANES;++j) written[j] = tl[j];
      }
    }
    abandon_bound -= (bpm_v4i_t)one;
  }
  memset(rec_event+T,0,8);
  bpm_v4u_t last_P[4], last_M[4];
  bpm_v4i_t last_S[4];
  BPM_V4_KEEP(last_P,last_M,last_S)
  // The lanes in column order; the pattern holds the search's state after column text_position-1
  uint64_t m = 0, t_done = (uint64_t)-1;
  bool same = true; // Lane 0 continues the search
  for (;;) {
    if (same) {
      uint8_t written_m;
      if (t_done+1==((m==0) ? 0 : X)) {
        // Only the columns with events, and the last column of a pending match
        const uint64_t lane_bits = 0x0101010101010101ull << m;
        t = t_done+1;
        while (t<T) {
          uint64_t t_event = t, w = 0; // Lane m's next column with an event (rec_event is zero past T)
          while (t_event<T) {
            memcpy(&w,rec_event+t_event,8);
            w &= lane_bits;
            if (w!=0) break;
            t_event += 8;
          }
          t_event = (w!=0) ? t_event+__builtin_ctzll(w)/8 : T;
          if (match_found) {
            // Columns before the pending match's last one only count down opt_steps_left
            const uint64_t t_last = t+opt_steps_left;
            if (t_event > t_last) t_event = t_last;
            if (t_event >= T) { opt_steps_left -= T-t; break; }
            opt_steps_left -= t_event-t;
          } else if (t_event >= T) {
            break;
          }
          t = t_event;
          text_position = s[m]+t;
          text_left = text_length-text_position;
          const uint8_t top_level_r = BPM_V4_TOP(rec_top[t],m);
          const int64_t current_score = rec_score[t][m];
          BPM_ALL_CHECK_MATCH(top_level_r,current_score);
          if (BPM_ALL_QUICK_ABANDON(top_level_r,current_score)) {
            uint64_t u;
            written_m = (m==0) ? written_a : top_level;
            for (u=t_done+1;u<=t;++u) {
              const uint8_t top_level_u = BPM_V4_TOP(rec_top[u],m);
              if (top_level_u > written_m) written_m = top_level_u;
            }
            BPM_V4_ABANDON(t,written_m);
          }
          ++t;
        }
        written_m = written[m];
      } else {
        // Every column: lane m met the search at a later snapshot
        written_m = top_level;
        for (t=t_done+1;t<T;++t) {
          text_position = s[m]+t;
          text_left = text_length-text_position;
          const uint8_t top_level_r = BPM_V4_TOP(rec_top[t],m);
          const int64_t current_score = rec_score[t][m];
          if (top_level_r > written_m) written_m = top_level_r;
          BPM_ALL_CHECK_MATCH(top_level_r,current_score);
          if (BPM_ALL_QUICK_ABANDON(top_level_r,current_score)) BPM_V4_ABANDON(t,written_m);
        }
      }
      for (k=0;k<written_m && k<NW;++k) { P[k]=last_P[k][m]; M[k]=last_M[k][m]; score[k]=last_S[k][m]; }
      top_level = tl[m];
      text_position = s[m]+T;
      text_left = text_length-text_position;
      if (++m==BPM_ALL_V4_LANES) break;
      t_done = X-1;
    } else {
      // The search alone up to lane m's next snapshot, or its last column
      uint64_t t_next = ((t_done+1)/BPM_ALL_V4_SNAP)*BPM_ALL_V4_SNAP+BPM_ALL_V4_SNAP-1;
      const bool lane_end = (t_next >= T-1);
      if (lane_end) t_next = T-1;
      const uint64_t column_end = text_position+(t_next-t_done);
      for (;text_position<column_end;++text_position,--text_left) {
        BPM_ALL_COLUMN_A();
      }
      if (lane_end) {
        if (++m==BPM_ALL_V4_LANES) break;
        t_done = X-1;
      } else {
        t_done = t_next;
      }
    }
    // Compare the search with lane m's snapshot after step t_done
    const uint64_t n = t_done/BPM_ALL_V4_SNAP;
    same = (BPM_V4_TOP(snap_top[n],m)==top_level);
    for (k=0;same && k<top_level;++k) {
      same = (snap_P[n][k][m]==P[k] && snap_M[n][k][m]==M[k] && snap_S[n][k][m]==score[k]);
    }
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
#define BPM_ALL_V4_KERNEL(NW) \
static __attribute__((noinline)) uint64_t bpm_compute_edit_distance_all_v4_w##NW( \
    pattern_t* const pattern,filtering_candidates_t* const filtering_candidates, \
    const uint64_t begin_position,const uint8_t* const text,const uint64_t text_length, \
    const uint64_t max_distance) { \
  return bpm_compute_edit_distance_all_v4(pattern,filtering_candidates, \
      begin_position,text,text_length,max_distance,NW); \
}
BPM_ALL_V4_KERNEL(1)
BPM_ALL_V4_KERNEL(2)
BPM_ALL_V4_KERNEL(3)
BPM_ALL_V4_KERNEL(4)
// Windows of the split long enough for the segments, whose words below the last end at bit 63
static bool bpm_compute_edit_distance_all_v4_fits(
    const bpm_pattern_t* const bpm_pattern,
    const uint64_t text_length) {
  if (text_length < BPM_ALL_V4_MIN_LENGTH) return false;
  if (text_length+(BPM_ALL_V4_LANES-1)*BPM_ALL_V4_OVERLAP > BPM_ALL_V4_LANES*BPM_ALL_V4_MAX_STEPS) return false;
  uint64_t i;
  for (i=0;i+1<bpm_pattern->pattern_num_words64;++i) {
    if (bpm_pattern->level_mask[i]!=BMP_W64_MASK) return false;
  }
  return true;
}
#endif
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
#ifdef BPM_ALL_V4_LANES
    if (bpm_compute_edit_distance_all_v4_fits(bpm_pattern,text_length)) {
      switch (num_words64) {
        case 1: num_matches_found = bpm_compute_edit_distance_all_v4_w1(pattern,
            filtering_candidates,begin_position,text,text_length,max_distance); break;
        case 2: num_matches_found = bpm_compute_edit_distance_all_v4_w2(pattern,
            filtering_candidates,begin_position,text,text_length,max_distance); break;
        case 3: num_matches_found = bpm_compute_edit_distance_all_v4_w3(pattern,
            filtering_candidates,begin_position,text,text_length,max_distance); break;
        default: num_matches_found = bpm_compute_edit_distance_all_v4_w4(pattern,
            filtering_candidates,begin_position,text,text_length,max_distance); break;
      }
    } else
#endif
    {
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
/*
 * BPM cut-off with quick abandon, several texts against one pattern
 *   Each text gets what bpm_compute_edit_distance(...,quick_abandon=true) gives it,
 *   and the pattern's P/M/score end as those calls in text order leave them. For 2-4
 *   word patterns two texts run at a time, in the two lanes of 2 x 64-bit vectors, each
 *   with its own top_level: a column computes the words below the larger top_level in
 *   both lanes, so a lane's words at or above its own top_level take garbage, which
 *   activation overwrites (a word's values at deactivation are kept for the write-back).
 *   A lane at the last word tracks its minimum in the vectors; the cut-off and abandon
 *   tests run per lane in scalar code, entered only when a lane's top level changes or
 *   it crosses the abandon bound. A lane that finishes takes the next text; the last
 *   text left runs alone in the scalar kernel. A text holding a byte above ENC_DNA_CHAR_N (it indexes
 *   PEQ past its word's rows, into P/M on the last word) goes through
 *   bpm_compute_edit_distance once the texts before it are written back.
 */
typedef struct {
  const uint8_t* text;
  uint64_t text_length;
  uint64_t position;          // Columns done
  uint64_t text_idx;
  uint64_t min_score;
  uint64_t min_score_column;
  uint8_t top_level;
  uint8_t max_level;
  bool abandoned;
  uint64_t P[BPM_BATCH_MAX_WORDS]; // Words at or above top_level: their values at deactivation
  uint64_t M[BPM_BATCH_MAX_WORDS];
  int64_t score[BPM_BATCH_MAX_WORDS];
} bpm_batch_lane_t;
typedef struct {
  int64_t text_idx[BPM_BATCH_MAX_WORDS]; // The last text (in order) that wrote each word, or -1
  uint64_t P[BPM_BATCH_MAX_WORDS];
  uint64_t M[BPM_BATCH_MAX_WORDS];
  int64_t score[BPM_BATCH_MAX_WORDS];
} bpm_batch_state_t;
static void bpm_batch_lane_start(
    bpm_batch_lane_t* const lane,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t text_idx,
    const uint8_t y,
    const int64_t* const init_score) {
  lane->text = text;
  lane->text_length = text_length;
  lane->position = 0;
  lane->text_idx = text_idx;
  lane->min_score = ALIGN_DISTANCE_INF;
  lane->min_score_column = ALIGN_COLUMN_INF;
  lane->top_level = y;
  lane->max_level = y;
  lane->abandoned = false;
  // bpm_reset_search_cutoff
  int64_t score = 0;
  uint64_t k;
  for (k=0;k<BPM_BATCH_MAX_WORDS;++k) {
    if (k < y) {
      score += init_score[k];
      lane->P[k] = BMP_W64_ONES; lane->M[k] = 0; lane->score[k] = score;
    } else {
      lane->P[k] = 0; lane->M[k] = 0; lane->score[k] = 0;
    }
  }
}
static void bpm_batch_lane_record(
    const bpm_batch_lane_t* const lane,
    bpm_batch_state_t* const state,
    uint64_t* const match_distances,
    uint64_t* const match_columns) {
  const uint64_t t = lane->text_idx;
  if (!lane->abandoned && lane->min_score!=ALIGN_DISTANCE_INF) {
    match_distances[t] = lane->min_score;
    match_columns[t] = lane->min_score_column;
  } else {
    match_distances[t] = ALIGN_DISTANCE_INF;
  }
  // BPM_W4_WRITE_BACK: the words below max_level
  uint64_t k;
  for (k=0;k<lane->max_level;++k) {
    if ((int64_t)t > state->text_idx[k]) {
      state->text_idx[k] = t;
      state->P[k] = lane->P[k];
      state->M[k] = lane->M[k];
      state->score[k] = lane->score[k];
    }
  }
}
/* The rest of a lane's text, alone (bpm_compute_edit_distance_wn from the lane's state) */
static inline __attribute__((always_inline)) void bpm_batch_lane_finish(
    const bpm_pattern_t* const bpm_pattern,
    const uint64_t max_distance,
    bpm_batch_lane_t* const lane,
    const uint64_t num_words64) {
  const uint64_t* const PEQ = bpm_pattern->PEQ;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  const int64_t* const init_score = bpm_pattern->init_score;
  const uint64_t pattern_length = bpm_pattern->pattern_length;
  const uint64_t mask0 = level_mask[0];
  const uint64_t mask1 = (num_words64>1) ? level_mask[1] : 0;
  const uint64_t mask2 = (num_words64>2) ? level_mask[2] : 0;
  const uint64_t mask3 = (num_words64>3) ? level_mask[3] : 0;
  const int64_t init_score1 = (num_words64>1) ? init_score[1] : 0;
  const int64_t init_score2 = (num_words64>2) ? init_score[2] : 0;
  const int64_t init_score3 = (num_words64>3) ? init_score[3] : 0;
  const uint64_t max_distance__1 = max_distance+1;
  const uint8_t top = num_words64-1;
  uint8_t top_level = lane->top_level, max_level = lane->max_level;
  uint64_t P0 = lane->P[0], P1 = lane->P[1], P2 = lane->P[2], P3 = lane->P[3];
  uint64_t M0 = lane->M[0], M1 = lane->M[1], M2 = lane->M[2], M3 = lane->M[3];
  int64_t score0 = lane->score[0], score1 = lane->score[1], score2 = lane->score[2], score3 = lane->score[3];
  uint64_t min_score = lane->min_score, min_score_column = lane->min_score_column;
  const uint8_t* const text = lane->text;
  const uint64_t text_length = lane->text_length;
  uint64_t text_position, text_left = text_length-lane->position;
  for (text_position=lane->position;text_position<text_length;++text_position,--text_left) {
    const uint8_t enc_char = text[text_position];
    uint64_t PHin=0,MHin=0,PHout,MHout;
    int64_t current_score;
    BPM_W4_COLUMN();
    if (top_level==num_words64 && current_score<=max_distance) {
      if (current_score < min_score)  {
        min_score_column = text_position;
        min_score = current_score;
      }
    } else if (min_score==ALIGN_DISTANCE_INF) {
      const uint64_t active_rows = (uint64_t)top_level*BMP_W64_LENGTH;
      const uint64_t rows_left = (pattern_length > active_rows) ? pattern_length-active_rows : 0;
      if (current_score+rows_left > text_left+max_distance) {
        lane->abandoned = true; // Quick abandon
        break;
      }
    }
  }
  lane->top_level = top_level;
  lane->max_level = max_level;
  lane->min_score = min_score;
  lane->min_score_column = min_score_column;
  lane->P[0] = P0; lane->P[1] = P1; lane->P[2] = P2; lane->P[3] = P3;
  lane->M[0] = M0; lane->M[1] = M1; lane->M[2] = M2; lane->M[3] = M3;
  lane->score[0] = score0; lane->score[1] = score1; lane->score[2] = score2; lane->score[3] = score3;
}
// Word k of both lanes; words below the last have level mask 1<<63
#define BPM_B2_WORD(k) { \
  const bpm_v2u_t Eq = Ep[k]; \
  const bpm_v2u_t Xv = Eq | VM##k; \
  const bpm_v2u_t _Eq = Eq | MHin; \
  const bpm_v2u_t Xh = (((_Eq & VP##k) + VP##k) ^ VP##k) | _Eq; \
  bpm_v2u_t Ph = VM##k | ~(Xh | VP##k); \
  bpm_v2u_t Mh = VP##k & Xh; \
  if (k+1<NW) { PHo##k = Ph >> 63; MHo##k = Mh >> 63; } \
  else { PHo##k = (Ph >> shift_last) & one; MHo##k = (Mh >> shift_last) & one; } \
  Ph = (Ph << 1) | PHin; \
  Mh = (Mh << 1) | MHin; \
  VP##k = Mh | ~(Xv | Ph); \
  VM##k = Ph & Xv; \
  VS##k += (bpm_v2i_t)(PHo##k - MHo##k); \
  PHin = PHo##k; MHin = MHo##k; \
}
/*
 * Word 0 with a shorter dependency chain (5 operations from P to P instead of 6), for
 * columns of one word: with Xs = Xh<<1, Ps = P<<1 and T = Xv | M<<1, the new P is
 * MAJ(Xs, Ps | (~T & 1), ~T) (bit 0 of Xs and Ps is 0, where the new P is ~T's bit 0)
 */
#define BPM_B2_WORD0_SHORT() { \
  const bpm_v2u_t Eq = Ep[0]; \
  const bpm_v2u_t Xv = Eq | VM0; \
  const bpm_v2u_t T = Xv | (VM0 << 1); \
  const bpm_v2u_t Ps = (VP0 << 1) | (~T & one); \
  const bpm_v2u_t Xh = (((Eq & VP0) + VP0) ^ VP0) | Eq; \
  const bpm_v2u_t Ph = VM0 | ~(Xh | VP0); \
  const bpm_v2u_t Mh = VP0 & Xh; \
  const bpm_v2u_t Xs = Xh << 1; \
  if (1<NW) { PHo0 = Ph >> 63; MHo0 = Mh >> 63; } \
  else { PHo0 = (Ph >> shift_last) & one; MHo0 = (Mh >> shift_last) & one; } \
  VP0 = (Xs & Ps) | (~T & (Xs | Ps)); \
  VM0 = (Ph << 1) & Xv; \
  VS0 += (bpm_v2i_t)(PHo0 - MHo0); \
}
// Lane 0 from word TA-1 (or TA), lane 1 from word TB-1 (or TB)
#define BPM_B2_PICK(V0,V1,V2,V3,T) (((T)==1) ? V0 : ((T)==2) ? V1 : ((T)==3) ? V2 : V3)
#define BPM_B2_SELECT(V0,V1,V2,V3,TA,TB,type) \
  (((TA)==(TB)) ? BPM_B2_PICK(V0,V1,V2,V3,TA) : \
   ((BPM_B2_PICK(V0,V1,V2,V3,TA) & (type)lane0) | (BPM_B2_PICK(V0,V1,V2,V3,TB) & (type)lane1)))
#define BPM_B2_LANES(TA,TB,T) \
  (((TA)==(T) && (TB)==(T)) ? all_lanes : ((TA)==(T)) ? lane0 : ((TB)==(T)) ? lane1 : (bpm_v2i_t)zero)
/*
 * Columns of both lanes at top levels (TA,TB) until a lane's top level would change,
 * a lane without a match crosses the abandon bound, or the shorter text ends. A lane
 * at the last word tracks its minimum here; a lane below it leaves only to extend
 * (BPM_W4_CUT_OFF's test) or drop (score above max_distance+init_score[top_level-1]).
 */
#define BPM_B2_LOOP(TA,TB) \
  if ((TA)<=NW && (TB)<=NW) { \
    const bool full_a = (TA)==NW, full_b = (TB)==NW; \
    const bpm_v2i_t full = BPM_B2_LANES((full_a ? (TA) : 0),(full_b ? (TB) : 0),NW); \
    const bpm_v2i_t partial = ~full; \
    for (;;) { \
      const bpm_v2u_t* const Ep = E + ((uint64_t)ta[j]*DNA__N_RANGE+tb[j])*NW; \
      bpm_v2u_t PHin = zero, MHin = zero; \
      if ((TA)==1 && (TB)==1) BPM_B2_WORD0_SHORT() else BPM_B2_WORD(0) \
      if ((TA)>1 || (TB)>1) BPM_B2_WORD(1) \
      if ((TA)>2 || (TB)>2) BPM_B2_WORD(2) \
      if ((TA)>3 || (TB)>3) BPM_B2_WORD(3) \
      Ssel = BPM_B2_SELECT(VS0,VS1,VS2,VS3,TA,TB,bpm_v2i_t); \
      ++j; \
      bpm_v2i_t event_v = (bpm_v2i_t)(Ssel > lim); \
      if ((TA)>1 || (TB)>1) event_v |= (bpm_v2i_t)(Ssel > hi); \
      if (full_a || full_b) { \
        const bpm_v2i_t lower = (bpm_v2i_t)(Ssel < VMIN) & full; \
        VMIN = (Ssel & lower) | (VMIN & ~lower); \
        VMINC = (colv & (bpm_v2u_t)lower) | (VMINC & ~(bpm_v2u_t)lower); \
        lim = (no_bound & lower) | (lim & ~lower); \
        colv += one; \
      } \
      if (!full_a || !full_b) { \
        const bpm_v2u_t PHs = BPM_B2_SELECT(PHo0,PHo1,PHo2,PHo3,TA,TB,bpm_v2u_t); \
        const bpm_v2u_t MHs = BPM_B2_SELECT(MHo0,MHo1,MHo2,MHo3,TA,TB,bpm_v2u_t); \
        const bpm_v2u_t EqT = BPM_B2_SELECT(Ep[0],Ep[1],Ep[2],Ep[3], \
            (full_a ? (TB)+1 : (TA)+1),(full_b ? (TA)+1 : (TB)+1),bpm_v2u_t); \
        const bpm_v2i_t extend = \
            (bpm_v2i_t)(((bpm_v2u_t)Ssel + MHs - PHs) <= max_distance_v) & \
            (bpm_v2i_t)(((MHs | EqT) & one) != zero) & partial; \
        event_v |= extend; \
      } \
      lim -= ones_i; \
      if (gem_expect_false((event_v[0] | event_v[1])!=0)) { event = true; break; } \
      if (j==run) break; \
    } \
  }
#define BPM_B2_GET(k,L,Pd,Md,Sd) \
  switch (k) { \
    case 0: Pd = VP0[L]; Md = VM0[L]; Sd = VS0[L]; break; \
    case 1: Pd = VP1[L]; Md = VM1[L]; Sd = VS1[L]; break; \
    case 2: Pd = VP2[L]; Md = VM2[L]; Sd = VS2[L]; break; \
    default: Pd = VP3[L]; Md = VM3[L]; Sd = VS3[L]; break; \
  }
#define BPM_B2_SET(k,L,Pv,Mv,Sv) \
  switch (k) { \
    case 0: VP0[L] = Pv; VM0[L] = Mv; VS0[L] = Sv; break; \
    case 1: VP1[L] = Pv; VM1[L] = Mv; VS1[L] = Sv; break; \
    case 2: VP2[L] = Pv; VM2[L] = Mv; VS2[L] = Sv; break; \
    default: VP3[L] = Pv; VM3[L] = Mv; VS3[L] = Sv; break; \
  }
#define BPM_B2_GET_S(k,L,Sd) \
  switch (k) { \
    case 0: Sd = VS0[L]; break; \
    case 1: Sd = VS1[L]; break; \
    case 2: Sd = VS2[L]; break; \
    default: Sd = VS3[L]; break; \
  }
#define BPM_B2_GET_H(k,L,PHd,MHd) \
  switch (k) { \
    case 0: PHd = PHo0[L]; MHd = MHo0[L]; break; \
    case 1: PHd = PHo1[L]; MHd = MHo1[L]; break; \
    case 2: PHd = PHo2[L]; MHd = MHo2[L]; break; \
    default: PHd = PHo3[L]; MHd = MHo3[L]; break; \
  }
/*
 * BPM_W4_CUT_OFF and the match/abandon tests on lane L's last column (position+j-1);
 * the minimum is VMIN[L] (max_distance+1 for none) and VMINC[L]
 */
#define BPM_B2_HANDLE(L) { \
  bpm_batch_lane_t* const lane = lanes+L; \
  const uint64_t column = lane->position+j-1; \
  const uint8_t enc_char = lane->text[column]; \
  uint8_t top_level = lane->top_level; \
  int64_t current_score = Ssel[L]; \
  uint64_t PHl = 0, MHl = 0; \
  BPM_B2_GET_H(top_level-1,L,PHl,MHl); \
  bool extended = false; \
  if (current_score<=max_distance__1) { \
    const uint64_t last_score = current_score+(MHl-PHl); \
    if (last_score<=max_distance && top_level-1<top && \
        (MHl || (PEQ[BPM_PATTERN_PEQ_IDX(top_level,enc_char)] & 1))) { \
      const uint64_t Peq = PEQ[BPM_PATTERN_PEQ_IDX(top_level,enc_char)]; \
      const uint64_t mask = level_mask[top_level]; \
      uint64_t Pv = BMP_W64_ONES, Mv = 0, PHin = PHl, MHin = MHl, PHout, MHout; \
      BPM_ADVANCE_BLOCK(Peq,mask,Pv,Mv,PHin,MHin,PHout,MHout); \
      const int64_t score_top = last_score + init_score[top_level] + (PHout-MHout); \
      BPM_B2_SET(top_level,L,Pv,Mv,score_top); \
      ++top_level; \
      if (lane->max_level < top_level) lane->max_level = top_level; \
      current_score = score_top; \
      extended = true; \
    } \
  } \
  if (!extended) { \
    while (top_level>1 && current_score > (max_distance+init_score[top_level-1])) { \
      --top_level; /* Word top_level leaves the search: keep its values */ \
      BPM_B2_GET(top_level,L,lane->P[top_level],lane->M[top_level],lane->score[top_level]); \
      BPM_B2_GET_S(top_level-1,L,current_score); \
    } \
  } \
  if (top_level==num_words64 && current_score<=max_distance) { \
    if (current_score < VMIN[L])  { \
      VMINC[L] = column; \
      VMIN[L] = current_score; \
    } \
  } else if (VMIN[L] > (int64_t)max_distance) { \
    const uint64_t active_rows = (uint64_t)top_level*BMP_W64_LENGTH; \
    const uint64_t rows_left = (pattern_length > active_rows) ? pattern_length-active_rows : 0; \
    const uint64_t text_left = lane->text_length-column; \
    if (current_score+rows_left > text_left+max_distance) lane->abandoned = true; \
  } \
  lane->top_level = top_level; \
}
// Lane L's active words and minimum to its record
#define BPM_B2_SAVE(L) { \
  uint64_t k; \
  for (k=0;k<lanes[L].top_level;++k) { \
    BPM_B2_GET(k,L,lanes[L].P[k],lanes[L].M[k],lanes[L].score[k]); \
  } \
  if (VMIN[L] <= (int64_t)max_distance) { \
    lanes[L].min_score = VMIN[L]; \
    lanes[L].min_score_column = VMINC[L]; \
  } \
}
// Lane L takes the next text (texts of length 0 are recorded on the way)
#define BPM_B2_TAKE(L) { \
  active##L = false; \
  while (next_text < num_texts) { \
    bpm_batch_lane_start(lanes+L,texts[next_text],text_lengths[next_text],next_text,y,init_score); \
    ++next_text; \
    if (lanes[L].text_length > 0) { \
      uint64_t k; \
      for (k=0;k<y;++k) { \
        BPM_B2_SET(k,L,lanes[L].P[k],lanes[L].M[k],lanes[L].score[k]); \
      } \
      VMIN[L] = max_distance__1; \
      VMINC[L] = ALIGN_COLUMN_INF; \
      active##L = true; \
      break; \
    } \
    bpm_batch_lane_record(lanes+L,state,match_distances,match_columns); \
  } \
}
#define BPM_B2_NEXT(L) \
  if (lanes[L].abandoned || lanes[L].position==lanes[L].text_length) { \
    BPM_B2_SAVE(L); \
    bpm_batch_lane_record(lanes+L,state,match_distances,match_columns); \
    BPM_B2_TAKE(L); \
  }
#define BPM_B2_ROWS_LEFT(tl) \
  ((pattern_length > (uint64_t)(tl)*BMP_W64_LENGTH) ? pattern_length-(uint64_t)(tl)*BMP_W64_LENGTH : 0)
// Lane L abandons at column position+j when its score exceeds this minus j
#define BPM_B2_LIM(L) \
  ((VMIN[L] <= (int64_t)max_distance) ? BPM_B2_NO_BOUND : \
   (int64_t)(lanes[L].text_length+max_distance-BPM_B2_ROWS_LEFT(lanes[L].top_level))-(int64_t)lanes[L].position)
#define BPM_B2_NO_BOUND (INT64_MAX/4)
static inline __attribute__((always_inline)) void bpm_compute_edit_distance_batch_v2(
    const bpm_pattern_t* const bpm_pattern,
    const uint8_t* const* const texts,
    const uint64_t* const text_lengths,
    const uint64_t num_texts,
    const uint64_t max_distance,
    uint64_t* const match_distances,
    uint64_t* const match_columns,
    bpm_batch_state_t* const state,
    const uint64_t NW) {
  // Pattern variables
  const uint64_t* const PEQ = bpm_pattern->PEQ;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  const int64_t* const init_score = bpm_pattern->init_score;
  const uint64_t pattern_length = bpm_pattern->pattern_length;
  const uint64_t num_words64 = NW;
  const uint8_t top = NW-1;
  const uint64_t max_distance__1 = max_distance+1;
  const uint8_t y = (max_distance>0) ? (max_distance+(BMP_W64_LENGTH-1))/BMP_W64_LENGTH : 1;
  const uint64_t shift_last = __builtin_ctzll(level_mask[NW-1]);
  // Pair tables {PEQ of lane 0's character, PEQ of lane 1's character}, NW words per pair
  bpm_v2u_t E[DNA__N_RANGE*DNA__N_RANGE*BPM_BATCH_MAX_WORDS];
  uint64_t a, b, k;
  for (a=0;a<DNA__N_RANGE;++a) {
    for (b=0;b<DNA__N_RANGE;++b) {
      for (k=0;k<NW;++k) {
        E[(a*DNA__N_RANGE+b)*NW+k] =
            (bpm_v2u_t){PEQ[BPM_PATTERN_PEQ_IDX(k,a)],PEQ[BPM_PATTERN_PEQ_IDX(k,b)]};
      }
    }
  }
  const bpm_v2u_t zero = {0,0}, one = {1,1};
  const bpm_v2u_t max_distance_v = {max_distance,max_distance};
  const bpm_v2i_t lane0 = {-1,0}, lane1 = {0,-1}, all_lanes = {-1,-1}, ones_i = {1,1};
  const bpm_v2i_t no_bound = {BPM_B2_NO_BOUND,BPM_B2_NO_BOUND};
  bpm_v2u_t VP0 = zero, VP1 = zero, VP2 = zero, VP3 = zero;
  bpm_v2u_t VM0 = zero, VM1 = zero, VM2 = zero, VM3 = zero;
  bpm_v2i_t VS0 = (bpm_v2i_t)zero, VS1 = (bpm_v2i_t)zero, VS2 = (bpm_v2i_t)zero, VS3 = (bpm_v2i_t)zero;
  bpm_v2u_t PHo0 = zero, PHo1 = zero, PHo2 = zero, PHo3 = zero;
  bpm_v2u_t MHo0 = zero, MHo1 = zero, MHo2 = zero, MHo3 = zero;
  bpm_v2i_t VMIN = (bpm_v2i_t)zero;
  bpm_v2u_t VMINC = zero;
  // Lanes
  bpm_batch_lane_t lanes[2];
  uint64_t next_text = 0;
  bool active0, active1;
  BPM_B2_TAKE(0);
  BPM_B2_TAKE(1);
  while (active0 && active1) {
    const uint64_t left0 = lanes[0].text_length-lanes[0].position;
    const uint64_t left1 = lanes[1].text_length-lanes[1].position;
    const uint64_t run = MIN(left0,left1);
    const uint8_t* const ta = lanes[0].text+lanes[0].position;
    const uint8_t* const tb = lanes[1].text+lanes[1].position;
    const uint8_t tl0 = lanes[0].top_level, tl1 = lanes[1].top_level;
    const bpm_v2i_t hi = {(int64_t)max_distance+init_score[tl0-1],(int64_t)max_distance+init_score[tl1-1]};
    bpm_v2i_t lim = {BPM_B2_LIM(0),BPM_B2_LIM(1)};
    bpm_v2u_t colv = {lanes[0].position,lanes[1].position};
    bpm_v2i_t Ssel = (bpm_v2i_t)zero;
    uint64_t j = 0;
    bool event = false;
    switch ((tl0-1)*4+(tl1-1)) {
      case  0: BPM_B2_LOOP(1,1); break;
      case  1: BPM_B2_LOOP(1,2); break;
      case  2: BPM_B2_LOOP(1,3); break;
      case  3: BPM_B2_LOOP(1,4); break;
      case  4: BPM_B2_LOOP(2,1); break;
      case  5: BPM_B2_LOOP(2,2); break;
      case  6: BPM_B2_LOOP(2,3); break;
      case  7: BPM_B2_LOOP(2,4); break;
      case  8: BPM_B2_LOOP(3,1); break;
      case  9: BPM_B2_LOOP(3,2); break;
      case 10: BPM_B2_LOOP(3,3); break;
      case 11: BPM_B2_LOOP(3,4); break;
      case 12: BPM_B2_LOOP(4,1); break;
      case 13: BPM_B2_LOOP(4,2); break;
      case 14: BPM_B2_LOOP(4,3); break;
      default: BPM_B2_LOOP(4,4); break;
    }
    if (event) {
      BPM_B2_HANDLE(0);
      BPM_B2_HANDLE(1);
    }
    lanes[0].position += j;
    lanes[1].position += j;
    BPM_B2_NEXT(0);
    BPM_B2_NEXT(1);
  }
  // The lane left runs alone
  if (active0) {
    BPM_B2_SAVE(0);
    bpm_batch_lane_finish(bpm_pattern,max_distance,lanes+0,NW);
    bpm_batch_lane_record(lanes+0,state,match_distances,match_columns);
  }
  if (active1) {
    BPM_B2_SAVE(1);
    bpm_batch_lane_finish(bpm_pattern,max_distance,lanes+1,NW);
    bpm_batch_lane_record(lanes+1,state,match_distances,match_columns);
  }
}
#define BPM_BATCH_KERNEL(NW) \
static __attribute__((noinline)) void bpm_compute_edit_distance_batch_w##NW( \
    const bpm_pattern_t* const bpm_pattern,const uint8_t* const* const texts, \
    const uint64_t* const text_lengths,const uint64_t num_texts,const uint64_t max_distance, \
    uint64_t* const match_distances,uint64_t* const match_columns,bpm_batch_state_t* const state) { \
  bpm_compute_edit_distance_batch_v2(bpm_pattern,texts,text_lengths,num_texts, \
      max_distance,match_distances,match_columns,state,NW); \
}
BPM_BATCH_KERNEL(2)
BPM_BATCH_KERNEL(3)
BPM_BATCH_KERNEL(4)
// Patterns the lanes take: 2-4 words, level masks 1<<63 below the last word, a single bit on it
static bool bpm_compute_edit_distance_batch_lanes(
    const bpm_pattern_t* const bpm_pattern) {
  const uint64_t num_words64 = bpm_pattern->pattern_num_words64;
  if (num_words64 < 2 || num_words64 > BPM_BATCH_MAX_WORDS) return false;
  uint64_t k;
  for (k=0;k+1<num_words64;++k) {
    if (bpm_pattern->level_mask[k] != (1ull<<63)) return false;
  }
  const uint64_t mask = bpm_pattern->level_mask[num_words64-1];
  return mask!=0 && (mask & (mask-1))==0;
}
/*
 * 1-word patterns (keys of 1-64 bp): top_level stays 1 (score[0] never exceeds
 * init_score[0]), so a column is one block and a minimum update. Up to
 * BPM_1W_LANES texts run together in vector lanes. Once a lane without a match
 * crosses the quick-abandon bound no later column of its text matches (rows_left
 * is 0 and the score falls at most 1 per column), so the lanes run on and a group
 * stops when every lane has crossed the bound or ended. The batch's last text sets
 * the pattern's P/M/score: its lane's state is taken at the column where it
 * crosses the bound, or at its end, as its call leaves it.
 */
#define BPM_1W_LANES 4
typedef uint64_t bpm_v1wu_t __attribute__((vector_size(8*BPM_1W_LANES)));
typedef int64_t bpm_v1wi_t __attribute__((vector_size(8*BPM_1W_LANES)));
/*
 * One column of word 0 (no carry in), with a shorter dependency chain (5 operations from
 * P to P instead of 6): with Xs = Xh<<1, Ps = P<<1 and T = Xv | M<<1, the new P is
 * MAJ(Xs, Ps | (~T & 1), ~T) (bit 0 of Xs and Ps is 0, where the new P is ~T's bit 0)
 */
#define BPM_1W_STEP(Eq) { \
  const bpm_v1wu_t Xv = Eq | VM; \
  const bpm_v1wu_t T = Xv | (VM << 1); \
  const bpm_v1wu_t Ps = (VP << 1) | (~T & one); \
  const bpm_v1wu_t Xh = (((Eq & VP) + VP) ^ VP) | Eq; \
  const bpm_v1wu_t Ph = VM | ~(Xh | VP); \
  const bpm_v1wu_t Mh = VP & Xh; \
  const bpm_v1wu_t Xs = Xh << 1; \
  PHo = (Ph >> shift_last) & one; \
  MHo = (Mh >> shift_last) & one; \
  VP = (Xs & Ps) | (~T & (Xs | Ps)); \
  VM = (Ph << 1) & Xv; \
}
static __attribute__((noinline)) void bpm_1w_lanes(
    const bpm_pattern_t* const bpm_pattern,
    const uint64_t max_distance,
    const uint8_t* const* const texts,
    const uint64_t* const text_lengths,
    const bool watch,
    uint64_t* const match_distances,
    uint64_t* const match_columns,
    uint64_t* const state) {
  const uint64_t* const PEQ = bpm_pattern->PEQ;
  const uint64_t shift_last = __builtin_ctzll(bpm_pattern->level_mask[0]);
  const int64_t init_score = bpm_pattern->init_score[0];
  const uint8_t* const t0 = texts[0];
  const uint8_t* const t1 = texts[1];
  const uint8_t* const t2 = texts[2];
  const uint8_t* const t3 = texts[3];
  const bpm_v1wu_t zero = {0,0,0,0}, one = {1,1,1,1};
  const bpm_v1wi_t no_bound = {INT64_MAX/4,INT64_MAX/4,INT64_MAX/4,INT64_MAX/4};
  bpm_v1wu_t VP = ~zero, VM = zero, VMINC = zero, colv = zero, lengths = zero, PHo = zero, MHo = zero;
  bpm_v1wi_t VS = (bpm_v1wi_t)zero, VMIN = (bpm_v1wi_t)zero, LIM = (bpm_v1wi_t)zero;
  uint64_t min_length = UINT64_MAX, max_length = 0, l, c;
  for (l=0;l<BPM_1W_LANES;++l) {
    VS[l] = init_score;
    VMIN[l] = max_distance+1;
    LIM[l] = text_lengths[l]+max_distance; // Abandon at column c: score > length-c+max_distance
    lengths[l] = text_lengths[l];
    min_length = MIN(min_length,text_lengths[l]);
    max_length = MAX(max_length,text_lengths[l]);
  }
  bool watching = watch;
  // Columns every lane has
  for (c=0;c<min_length;++c) {
    const bpm_v1wu_t Eq = {PEQ[t0[c]],PEQ[t1[c]],PEQ[t2[c]],PEQ[t3[c]]};
    BPM_1W_STEP(Eq);
    VS += (bpm_v1wi_t)(PHo - MHo);
    const bpm_v1wi_t lower = (bpm_v1wi_t)(VS < VMIN);
    VMIN = (VS & lower) | (VMIN & ~lower);
    VMINC = (colv & (bpm_v1wu_t)lower) | (VMINC & ~(bpm_v1wu_t)lower);
    LIM = (no_bound & lower) | (LIM & ~lower);
    const bpm_v1wi_t done = (bpm_v1wi_t)(VS > LIM);
    LIM -= (bpm_v1wi_t)one;
    colv += one;
    if (watching && done[3]) {
      state[0] = VP[3]; state[1] = VM[3]; state[2] = (uint64_t)VS[3];
      watching = false;
    }
    if (gem_expect_false((done[0] & done[1] & done[2] & done[3])!=0)) goto lanes_done;
  }
  // Columns past the shortest text: a lane past its end keeps its state
  for (;c<max_length;++c) {
    bpm_v1wu_t Eq;
    Eq[0] = (c<text_lengths[0]) ? PEQ[t0[c]] : 0;
    Eq[1] = (c<text_lengths[1]) ? PEQ[t1[c]] : 0;
    Eq[2] = (c<text_lengths[2]) ? PEQ[t2[c]] : 0;
    Eq[3] = (c<text_lengths[3]) ? PEQ[t3[c]] : 0;
    const bpm_v1wu_t active = (bpm_v1wu_t)(colv < lengths);
    const bpm_v1wu_t P_old = VP, M_old = VM;
    BPM_1W_STEP(Eq);
    VP = (VP & active) | (P_old & ~active);
    VM = (VM & active) | (M_old & ~active);
    VS += (bpm_v1wi_t)((PHo - MHo) & active);
    const bpm_v1wi_t lower = (bpm_v1wi_t)(VS < VMIN) & (bpm_v1wi_t)active;
    VMIN = (VS & lower) | (VMIN & ~lower);
    VMINC = (colv & (bpm_v1wu_t)lower) | (VMINC & ~(bpm_v1wu_t)lower);
    LIM = (no_bound & lower) | (LIM & ~lower);
    const bpm_v1wi_t crossed = (bpm_v1wi_t)(VS > LIM) & (bpm_v1wi_t)active;
    const bpm_v1wi_t done = crossed | ~(bpm_v1wi_t)active;
    LIM -= (bpm_v1wi_t)one;
    colv += one;
    if (watching && crossed[3]) {
      state[0] = VP[3]; state[1] = VM[3]; state[2] = (uint64_t)VS[3];
      watching = false;
    }
    if ((done[0] & done[1] & done[2] & done[3])!=0) goto lanes_done;
  }
lanes_done:
  for (l=0;l<BPM_1W_LANES;++l) {
    if (VMIN[l] <= (int64_t)max_distance) {
      match_distances[l] = VMIN[l];
      match_columns[l] = VMINC[l];
    } else {
      match_distances[l] = ALIGN_DISTANCE_INF;
    }
  }
  if (watching) { // The watched text ended without crossing the bound
    state[0] = VP[3]; state[1] = VM[3]; state[2] = (uint64_t)VS[3];
  }
}
/* A group of up to BPM_1W_LANES texts (the last lane repeated); one text takes the call */
static void bpm_1w_group(
    const bpm_pattern_t* const bpm_pattern,
    const uint64_t max_distance,
    const uint8_t* const* const texts,
    const uint64_t* const text_lengths,
    const uint64_t* const text_idx,
    const uint64_t num_lanes,
    const uint64_t last_text_idx,
    uint64_t* const match_distances,
    uint64_t* const match_columns) {
  if (num_lanes==1) {
    bpm_compute_edit_distance(bpm_pattern,texts[0],text_lengths[0],
        match_distances+text_idx[0],match_columns+text_idx[0],max_distance,true);
    return;
  }
  const uint8_t* lane_texts[BPM_1W_LANES];
  uint64_t lane_lengths[BPM_1W_LANES], lane_distances[BPM_1W_LANES], lane_columns[BPM_1W_LANES];
  uint64_t state[3], l;
  for (l=0;l<BPM_1W_LANES;++l) {
    const uint64_t k = (l<num_lanes) ? l : num_lanes-1;
    lane_texts[l] = texts[k];
    lane_lengths[l] = text_lengths[k];
  }
  const bool has_last = (text_idx[num_lanes-1]==last_text_idx);
  bpm_1w_lanes(bpm_pattern,max_distance,lane_texts,lane_lengths,
      has_last,lane_distances,lane_columns,state);
  for (l=0;l<num_lanes;++l) {
    match_distances[text_idx[l]] = lane_distances[l];
    if (lane_distances[l]!=ALIGN_DISTANCE_INF) match_columns[text_idx[l]] = lane_columns[l];
  }
  if (has_last) {
    bpm_pattern->P[0] = state[0];
    bpm_pattern->M[0] = state[1];
    bpm_pattern->score[0] = (int64_t)state[2];
  }
}
static void bpm_compute_edit_distance_batch_1w(
    const bpm_pattern_t* const bpm_pattern,
    const uint8_t* const* const texts,
    const uint64_t* const text_lengths,
    const uint64_t num_texts,
    const uint64_t max_distance,
    uint64_t* const match_distances,
    uint64_t* const match_columns) {
  const uint8_t* group_texts[BPM_1W_LANES];
  uint64_t group_lengths[BPM_1W_LANES], group_idx[BPM_1W_LANES], num_lanes = 0, n;
  for (n=0;n<num_texts;++n) {
    if (!bpm_text_within_peq(texts[n],text_lengths[n])) {
      // Texts before it first, so that the pattern's state ends as the calls in order leave it
      if (num_lanes > 0) {
        bpm_1w_group(bpm_pattern,max_distance,group_texts,group_lengths,group_idx,
            num_lanes,num_texts-1,match_distances,match_columns);
        num_lanes = 0;
      }
      bpm_compute_edit_distance(bpm_pattern,texts[n],text_lengths[n],
          match_distances+n,match_columns+n,max_distance,true);
      continue;
    }
    group_texts[num_lanes] = texts[n];
    group_lengths[num_lanes] = text_lengths[n];
    group_idx[num_lanes] = n;
    if (++num_lanes==BPM_1W_LANES) {
      bpm_1w_group(bpm_pattern,max_distance,group_texts,group_lengths,group_idx,
          num_lanes,num_texts-1,match_distances,match_columns);
      num_lanes = 0;
    }
  }
  if (num_lanes > 0) {
    bpm_1w_group(bpm_pattern,max_distance,group_texts,group_lengths,group_idx,
        num_lanes,num_texts-1,match_distances,match_columns);
  }
}
void bpm_compute_edit_distance_batch(
    const bpm_pattern_t* const bpm_pattern,
    const uint8_t* const* const texts,
    const uint64_t* const text_lengths,
    const uint64_t num_texts,
    uint64_t max_distance,
    uint64_t* const match_distances,
    uint64_t* const match_columns) {
  uint64_t n;
  if (num_texts >= 2 && bpm_pattern->pattern_num_words64==1 &&
      bpm_pattern->pattern_length > 0 && bpm_pattern->level_mask[0]!=0) {
    if (max_distance >= bpm_pattern->pattern_length) {
      max_distance = bpm_pattern->pattern_length-1; // Correct max-distance
    }
    bpm_compute_edit_distance_batch_1w(bpm_pattern,texts,text_lengths,
        num_texts,max_distance,match_distances,match_columns);
    return;
  }
  if (num_texts < 2 || !bpm_compute_edit_distance_batch_lanes(bpm_pattern)) {
    for (n=0;n<num_texts;++n) {
      bpm_compute_edit_distance(bpm_pattern,texts[n],text_lengths[n],
          match_distances+n,match_columns+n,max_distance,true);
    }
    return;
  }
  // Initialize search
  const uint64_t num_words64 = bpm_pattern->pattern_num_words64;
  if (max_distance >= bpm_pattern->pattern_length) {
    max_distance = bpm_pattern->pattern_length-1; // Correct max-distance
  }
  // Runs of texts within PEQ go to the lanes; the others alone, in order
  uint64_t begin = 0;
  for (n=0;n<=num_texts;++n) {
    if (n<num_texts && bpm_text_within_peq(texts[n],text_lengths[n])) continue;
    const uint64_t num_run = n-begin;
    if (num_run==1) {
      bpm_compute_edit_distance(bpm_pattern,texts[begin],text_lengths[begin],
          match_distances+begin,match_columns+begin,max_distance,true);
    } else if (num_run > 1) {
      bpm_batch_state_t state = {{-1,-1,-1,-1},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
      switch (num_words64) {
        case 2: bpm_compute_edit_distance_batch_w2(bpm_pattern,texts+begin,text_lengths+begin,
            num_run,max_distance,match_distances+begin,match_columns+begin,&state); break;
        case 3: bpm_compute_edit_distance_batch_w3(bpm_pattern,texts+begin,text_lengths+begin,
            num_run,max_distance,match_distances+begin,match_columns+begin,&state); break;
        default: bpm_compute_edit_distance_batch_w4(bpm_pattern,texts+begin,text_lengths+begin,
            num_run,max_distance,match_distances+begin,match_columns+begin,&state); break;
      }
      // Leave the search state as the calls in text order leave it
      uint64_t k;
      for (k=0;k<num_words64;++k) {
        if (state.text_idx[k] < 0) continue;
        bpm_pattern->P[k] = state.P[k];
        bpm_pattern->M[k] = state.M[k];
        bpm_pattern->score[k] = state.score[k];
      }
    }
    if (n<num_texts) { // A byte above ENC_DNA_CHAR_N
      bpm_compute_edit_distance(bpm_pattern,texts[n],text_lengths[n],
          match_distances+n,match_columns+n,max_distance,true);
    }
    begin = n+1;
  }
}
