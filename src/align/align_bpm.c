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
 *   Alignment module using BPM-algorithm to compute levenshtein full-alignment
 *   (Myers' Fast Bit-Vector algorithm to compute levenshtein full-alignment)
 * TODO
 *   Block (@top) can surpass 256 bits (uint8_t) size
 */

#include "align/alignment.h"
#include "align/align_bpm.h"
#include "align/align_bpm_pattern.h"
#include "align/align_bpm_distance.h"
#include "matches/matches.h"
#include "matches/matches_cigar.h"

/*
 * BPM-DP-Matrix for patterns of at most 4 words (BPM_MATRIX_MAX_WORDS)
 *   Each word's Pv/Mv and score are carried in locals (P0..P3, M0..M3,
 *   score0..score3): what the array loop reloads from the matrix is what it
 *   stored one column earlier, or the ONES/0 of an activation. Every matrix
 *   store is kept (align_bpm_backtrace_matrix reads them) and score is
 *   written back on return. A text byte above ENC_DNA_CHAR_N indexes PEQ
 *   into the pattern's P/M, which this function never writes, so unlike the
 *   distance kernels these need no text guard.
 */
#define BPM_MATRIX_MAX_WORDS 4
#define BPM_MATRIX_STEP(k) { \
  uint64_t Pv_in = P##k; \
  uint64_t Mv_in = M##k; \
  const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(k,enc_char)]; \
  BPM_ADVANCE_BLOCK(Eq,mask##k,Pv_in,Mv_in,PHin,MHin,PHout,MHout); \
  Pv_next[k] = Pv_in; \
  Mv_next[k] = Mv_in; \
  P##k = Pv_in; \
  M##k = Mv_in; \
  score##k += PHout-MHout; \
  PHin = PHout; \
  MHin = MHout; \
}
#define BPM_MATRIX_EXTEND(k) { \
  const uint64_t Peq = PEQ[BPM_PATTERN_PEQ_IDX(k,enc_char)]; \
  uint64_t Pv_in = BMP_W64_ONES; \
  uint64_t Mv_in = 0; \
  Pv_col[k] = BMP_W64_ONES; \
  Mv_col[k] = 0; \
  BPM_ADVANCE_BLOCK(Peq,mask##k,Pv_in,Mv_in,PHin,MHin,PHout,MHout); \
  Pv_next[k] = Pv_in; \
  Mv_next[k] = Mv_in; \
  P##k = Pv_in; \
  M##k = Mv_in; \
  score##k = last_score + init_score##k + (PHout-MHout); \
  top_level = k+1; \
  if (max_level < k+1) max_level = k+1; \
  current_score = score##k; \
}
/* The array loop's while; it never goes below level 1, since score0 <= init_score0 */
#define BPM_MATRIX_DOWN_1 { top_level = 1; current_score = score0; }
#define BPM_MATRIX_DOWN_2 { \
  if (score1 > (max_distance+init_score1)) BPM_MATRIX_DOWN_1 \
  else { top_level = 2; current_score = score1; } \
}
#define BPM_MATRIX_DOWN_3 { \
  if (score2 > (max_distance+init_score2)) BPM_MATRIX_DOWN_2 \
  else { top_level = 3; current_score = score2; } \
}
#define BPM_MATRIX_DOWN_4 { \
  if (score3 > (max_distance+init_score3)) BPM_MATRIX_DOWN_3 \
  else { top_level = 4; current_score = score3; } \
}
/* Cut-off at top_level T (last word L=T-1) */
#define BPM_MATRIX_CUT_OFF(T,L) \
  if (gem_expect_false(score##L<=max_distance__1 && L<top)) { \
    const uint64_t last_score = score##L+(MHin-PHin); \
    if (last_score<=max_distance && (MHin || (PEQ[BPM_PATTERN_PEQ_IDX(T,enc_char)] & 1))) BPM_MATRIX_EXTEND(T) \
    else BPM_MATRIX_DOWN_##T \
  } else BPM_MATRIX_DOWN_##T
static inline __attribute__((always_inline)) void align_bpm_compute_matrix_wn(
    bpm_pattern_t* const bpm_pattern,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t max_distance,
    uint64_t* const Pv,
    uint64_t* const Mv,
    bpm_align_matrix_t* const bpm_align_matrix,
    const uint64_t num_words64) {
  // Pattern variables
  const uint64_t* const PEQ = bpm_pattern->PEQ;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  int64_t* const score = bpm_pattern->score;
  const int64_t* const init_score = bpm_pattern->init_score;
  const uint64_t mask0 = level_mask[0];
  const uint64_t mask1 = (num_words64>1) ? level_mask[1] : 0;
  const uint64_t mask2 = (num_words64>2) ? level_mask[2] : 0;
  const uint64_t mask3 = (num_words64>3) ? level_mask[3] : 0;
  const int64_t init_score0 = init_score[0];
  const int64_t init_score1 = (num_words64>1) ? init_score[1] : 0;
  const int64_t init_score2 = (num_words64>2) ? init_score[2] : 0;
  const int64_t init_score3 = (num_words64>3) ? init_score[3] : 0;
  const uint64_t max_distance__1 = max_distance+1;
  const uint8_t top = num_words64-1;
  // Initialize search (bpm_reset_search_cutoff)
  const uint8_t y = (max_distance>0) ? (max_distance+(BMP_W64_LENGTH-1))/BMP_W64_LENGTH : 1;
  uint8_t top_level = y, max_level = y;
  uint64_t P0 = BMP_W64_ONES, P1 = 0, P2 = 0, P3 = 0;
  uint64_t M0 = 0, M1 = 0, M2 = 0, M3 = 0;
  int64_t score0 = init_score0, score1 = 0, score2 = 0, score3 = 0;
  Pv[0] = BMP_W64_ONES; Mv[0] = 0;
  if (y>1) { P1 = BMP_W64_ONES; Pv[1] = BMP_W64_ONES; Mv[1] = 0; score1 = score0 + init_score1; }
  if (y>2) { P2 = BMP_W64_ONES; Pv[2] = BMP_W64_ONES; Mv[2] = 0; score2 = score1 + init_score2; }
  if (y>3) { P3 = BMP_W64_ONES; Pv[3] = BMP_W64_ONES; Mv[3] = 0; score3 = score2 + init_score3; }
  uint64_t min_score = ALIGN_DISTANCE_INF, min_score_column = ALIGN_COLUMN_INF;
  // Advance in DP-bit_encoded matrix
  uint64_t* Pv_col = Pv;
  uint64_t* Mv_col = Mv;
  uint64_t text_position;
  for (text_position=0;text_position<text_length;++text_position) {
    uint64_t* const Pv_next = Pv_col + num_words64;
    uint64_t* const Mv_next = Mv_col + num_words64;
    const uint8_t enc_char = text[text_position];
    uint64_t PHin=0,MHin=0,PHout,MHout;
    int64_t current_score;
    if (num_words64==1 || top_level==1) {
      BPM_MATRIX_STEP(0)
      BPM_MATRIX_CUT_OFF(1,0)
    } else if (num_words64==2 || top_level==2) {
      BPM_MATRIX_STEP(0) BPM_MATRIX_STEP(1)
      BPM_MATRIX_CUT_OFF(2,1)
    } else if (num_words64==3 || top_level==3) {
      BPM_MATRIX_STEP(0) BPM_MATRIX_STEP(1) BPM_MATRIX_STEP(2)
      BPM_MATRIX_CUT_OFF(3,2)
    } else {
      BPM_MATRIX_STEP(0) BPM_MATRIX_STEP(1) BPM_MATRIX_STEP(2) BPM_MATRIX_STEP(3)
      BPM_MATRIX_DOWN_4
    }
    // Check match
    if (top_level==num_words64 && current_score<=max_distance) {
      if (current_score < min_score)  {
        min_score_column = text_position;
        min_score = current_score;
      }
    }
    Pv_col = Pv_next;
    Mv_col = Mv_next;
  }
  // Write back the scores the array loop would have left
  score[0] = score0;
  if (num_words64>1 && max_level>1) score[1] = score1;
  if (num_words64>2 && max_level>2) score[2] = score2;
  if (num_words64>3 && max_level>3) score[3] = score3;
  // Return optimal column/distance
  bpm_align_matrix->min_score = min_score;
  bpm_align_matrix->min_score_column = min_score_column;
}
/* Out of line, so each column loop gets the registers to itself */
#define BPM_MATRIX_KERNEL(NW) \
static __attribute__((noinline)) void align_bpm_compute_matrix_w##NW( \
    bpm_pattern_t* const bpm_pattern,const uint8_t* const text,const uint64_t text_length, \
    const uint64_t max_distance,uint64_t* const Pv,uint64_t* const Mv, \
    bpm_align_matrix_t* const bpm_align_matrix) { \
  align_bpm_compute_matrix_wn(bpm_pattern,text,text_length, \
      max_distance,Pv,Mv,bpm_align_matrix,NW); \
}
BPM_MATRIX_KERNEL(1)
BPM_MATRIX_KERNEL(2)
BPM_MATRIX_KERNEL(3)
BPM_MATRIX_KERNEL(4)
/*
 * BPM. Compute BPM-DP-Matrix
 */
void align_bpm_compute_matrix(
    bpm_pattern_t* const bpm_pattern,
    uint8_t* const text,
    const uint64_t text_length,
    uint64_t max_distance,
    bpm_align_matrix_t* const bpm_align_matrix,
    mm_allocator_t* const mm_allocator) {
  PROF_START(GP_BPM_COMPUTE_MATRIX);
  // Pattern variables
  const uint64_t* PEQ = bpm_pattern->PEQ;
  const uint64_t num_words64 = bpm_pattern->pattern_num_words64;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  int64_t* const score = bpm_pattern->score;
  const int64_t* const init_score = bpm_pattern->init_score;
  // Allocate auxiliary matrix
  const uint64_t aux_matrix_size = num_words64*UINT64_SIZE*(text_length+1); /* (+1 base-column) */
  uint64_t* const Pv = (uint64_t*)mm_allocator_malloc(mm_allocator,aux_matrix_size);
  uint64_t* const Mv = (uint64_t*)mm_allocator_malloc(mm_allocator,aux_matrix_size);
  bpm_align_matrix->Mv = Mv;
  bpm_align_matrix->Pv = Pv;
  // Initialize search
  if (max_distance >= bpm_pattern->pattern_length) {
    max_distance = bpm_pattern->pattern_length-1; // Correct max-distance
  }
  switch (num_words64) {
    case 1: align_bpm_compute_matrix_w1(bpm_pattern,text,text_length,max_distance,Pv,Mv,bpm_align_matrix); PROF_STOP(GP_BPM_COMPUTE_MATRIX); return;
    case 2: align_bpm_compute_matrix_w2(bpm_pattern,text,text_length,max_distance,Pv,Mv,bpm_align_matrix); PROF_STOP(GP_BPM_COMPUTE_MATRIX); return;
    case 3: align_bpm_compute_matrix_w3(bpm_pattern,text,text_length,max_distance,Pv,Mv,bpm_align_matrix); PROF_STOP(GP_BPM_COMPUTE_MATRIX); return;
    case 4: align_bpm_compute_matrix_w4(bpm_pattern,text,text_length,max_distance,Pv,Mv,bpm_align_matrix); PROF_STOP(GP_BPM_COMPUTE_MATRIX); return;
    default: break;
  }
  const uint64_t max_distance__1 = max_distance+1;
  const uint8_t top = num_words64-1;
  uint64_t min_score = ALIGN_DISTANCE_INF, min_score_column = ALIGN_COLUMN_INF;
  uint8_t top_level;
  bpm_reset_search_cutoff(&top_level,Pv,Mv,score,init_score,max_distance);
  // Advance in DP-bit_encoded matrix
  uint64_t text_position;
  for (text_position=0;text_position<text_length;++text_position) {
    // Fetch next character
    const uint8_t enc_char = text[text_position];
    // Advance all blocks
    uint64_t i,PHin=0,MHin=0,PHout,MHout;
    for (i=0;i<top_level;++i) {
      /* Calculate Step Data */
      const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words64,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words64;
      uint64_t Pv_in = Pv[bdp_idx];
      uint64_t Mv_in = Mv[bdp_idx];
      const uint64_t mask = level_mask[i];
      const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i,enc_char)];
      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);
      /* Adjust score and swap propagate Hv */
      score[i] += PHout-MHout;
      Pv[next_bdp_idx] = Pv_in;
      Mv[next_bdp_idx] = Mv_in;
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
        const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words64,top_level);
        const uint64_t next_bdp_idx = bdp_idx+num_words64;
        uint64_t Pv_in = BMP_W64_ONES;
        uint64_t Mv_in = 0;
        Pv[bdp_idx] = BMP_W64_ONES;
        Mv[bdp_idx] = 0;
        const uint64_t mask = level_mask[top_level];
        /* Compute Block */
        BPM_ADVANCE_BLOCK(Peq,mask,Pv_in,Mv_in,PHin,MHin,PHout,MHout);
        /* Save Block Pv,Mv */
        Pv[next_bdp_idx]=Pv_in;
        Mv[next_bdp_idx]=Mv_in;
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
    }
  }
  // Return optimal column/distance
  bpm_align_matrix->min_score = min_score;
  bpm_align_matrix->min_score_column = min_score_column;
  PROF_STOP(GP_BPM_COMPUTE_MATRIX);
}
/*
 * BPM. Recover CIGAR from a matching string
 */
const cigar_t align_bpm_backtrace_operation[16] = {
    /*
     * left_gap_alignment x deletion x insertion x match
     */
//  if (!left_gap_alignment)
//  if (match) {
//    operation = cigar_match;
//  } else if (deletion) {
//    operation = cigar_del;
//  } else if (insertion) {
//    operation = cigar_ins;
//  } else {
//    operation = cigar_mismatch;
//  }
    /* 0000 */ cigar_mismatch,
    /* 0001 */ cigar_match,
    /* 0010 */ cigar_ins,
    /* 0011 */ cigar_match,
    /* 0100 */ cigar_del,
    /* 0101 */ cigar_match,
    /* 0110 */ cigar_del,
    /* 0111 */ cigar_match,
//  if (left_gap_alignment)
//  if (deletion) {
//    operation = cigar_del;
//  } else if (insertion) {
//    operation = cigar_ins;
//  } else if (match) {
//    operation = cigar_match;
//  } else {
//    operation = cigar_mismatch;
//  }
    /* 1000 */ cigar_mismatch,
    /* 1001 */ cigar_match,
    /* 1010 */ cigar_ins,
    /* 1011 */ cigar_ins,
    /* 1100 */ cigar_del,
    /* 1101 */ cigar_del,
    /* 1110 */ cigar_del,
    /* 1111 */ cigar_del,
};
void align_bpm_backtrace_matrix(
    const bpm_pattern_t* const bpm_pattern,
    const uint8_t* const key,
    uint8_t* const text,
    const bool left_gap_alignment,
    bpm_align_matrix_t* const bpm_align_matrix,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector) {
  // Parameters
  const uint64_t pattern_length = bpm_pattern->pattern_length;
  const uint64_t* const Pv = bpm_align_matrix->Pv;
  const uint64_t* const Mv = bpm_align_matrix->Mv;
  // Allocate CIGAR string memory (worst case)
  match_alignment->cigar_offset = vector_get_used(cigar_vector); // Set CIGAR offset
  vector_reserve_additional(cigar_vector,MIN(pattern_length,2*bpm_align_matrix->min_score+1)); // Reserve
  cigar_element_t* cigar_buffer = vector_get_free_elm(cigar_vector,cigar_element_t); // Sentinel
  cigar_element_t* const cigar_buffer_base = cigar_buffer;
  cigar_buffer->type = cigar_null; // Trick
  // Retrieve the alignment. Store the match
  const uint64_t num_words64 = bpm_pattern->pattern_num_words64;
  int64_t match_effective_length = pattern_length;
  int64_t h = bpm_align_matrix->min_score_column;
  int64_t v = pattern_length - 1;
  while (v >= 0 && h >= 0) {
    const uint8_t block = v / UINT64_LENGTH;
    const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(h+1,num_words64,block);
    const uint64_t mask = 1L << (v % UINT64_LENGTH);
    // Select CIGAR operation
    //    const uint64_t left_gap_alignment_b = (left_gap_alignment==true);
    //    const uint64_t deletion = (Pv[bdp_idx] & mask) != 0;
    //    const uint64_t insertion = (Mv[(bdp_idx-num_words64)] & mask) != 0;
    //    const uint64_t match = (text[h]==key[v] && key[v]!=ENC_DNA_CHAR_N) != 0; // N's Inequality
    //    const uint64_t operation_idx = left_gap_alignment_b<<3 | deletion<<2 | insertion<<1 | match;
    //    const cigar_t operation = align_bpm_backtrace_operation[operation_idx];
    // CIGAR operation Test
    #define IS_DELETION   (Pv[bdp_idx] & mask)
    #define IS_INSERTION  (Mv[(bdp_idx-num_words64)] & mask)
    #define IS_MATCH      ((text[h]==key[v] && key[v]!=ENC_DNA_CHAR_N))
    cigar_t operation;
    if (!left_gap_alignment) {
      if (IS_DELETION) {
        operation = cigar_del;
      } else if (IS_INSERTION) {
        operation = cigar_ins;
      } else if (IS_MATCH) {
        operation = cigar_match;
      } else {
        operation = cigar_mismatch;
      }
    } else {
      if (IS_MATCH) {
        operation = cigar_match;
      } else if (IS_DELETION) {
        operation = cigar_del;
      } else if (IS_INSERTION) {
        operation = cigar_ins;
      } else {
        operation = cigar_mismatch;
      }
    }
    // Add selected CIGAR operation
    switch (operation) {
      case cigar_del:
        matches_cigar_buffer_add_cigar_element(&cigar_buffer,cigar_del,1); // Deletion <-1>@v
        --v; --match_effective_length;
        break;
      case cigar_ins:
        matches_cigar_buffer_add_cigar_element(&cigar_buffer,cigar_ins,1); // Insertion <+1>@v
        --h; ++match_effective_length;
        break;
      case cigar_mismatch:
        matches_cigar_buffer_add_mismatch(&cigar_buffer,text[h]);
        --h; --v;
        break;
      case cigar_match:
        matches_cigar_buffer_add_cigar_element(&cigar_buffer,cigar_match,1); // Match
        --h; --v;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  if (v >= 0) {
    matches_cigar_buffer_add_cigar_element(&cigar_buffer,cigar_del,v+1); // <-(@v+1)>@v
    match_effective_length -= v+1;
  }
  if (h >= 0) {
    match_alignment->match_position += h+1; // We need to correct the matching_position
  }
  // Set effective length
  match_alignment->effective_length = match_effective_length;
  // Set CIGAR buffer used
  if (cigar_buffer->type!=cigar_null) ++(cigar_buffer);
  const uint64_t num_cigar_elements = cigar_buffer - cigar_buffer_base;
  match_alignment->cigar_length = num_cigar_elements; // Set CIGAR length
  // Reverse CIGAR Elements
  if (num_cigar_elements > 0) {
    const uint64_t middle_point = num_cigar_elements/2;
    uint64_t i;
    for (i=0;i<middle_point;++i) {
      SWAP(cigar_buffer_base[i],cigar_buffer_base[num_cigar_elements-i-1]);
    }
  }
  // Set used
  vector_add_used(cigar_vector,num_cigar_elements);
}
/*
 * BPM Align match
 */
void align_bpm_match(
    bpm_pattern_t* const bpm_pattern,
    uint8_t* const key,
    uint8_t* const text,
    const uint64_t text_length,
    const uint64_t max_distance,
    const bool left_gap_alignment,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_allocator_t* const mm_allocator) {
  // Fill Matrix (Pv,Mv)
  mm_allocator_push_state(mm_allocator); // Save allocator state
  bpm_align_matrix_t bpm_align_matrix;
  align_bpm_compute_matrix(
      bpm_pattern,text,text_length,
      max_distance,&bpm_align_matrix,mm_allocator);
  // Set distance
  match_alignment->score = bpm_align_matrix.min_score;
  if (bpm_align_matrix.min_score == ALIGN_DISTANCE_INF) {
    match_alignment->cigar_length = 0;
    mm_allocator_pop_state(mm_allocator); // Free
    return;
  }
  // Backtrace and generate CIGAR
  align_bpm_backtrace_matrix(
      bpm_pattern,key,text,left_gap_alignment,
      &bpm_align_matrix,match_alignment,cigar_vector);
  // Free
  mm_allocator_pop_state(mm_allocator);
}
