/*
 * PROJECT: GEMMapper
 * FILE: align_bpm.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 * TODO
 *   - Block can surpass 256->uint8_t size
 */

#include "align/align_bpm.h"
#include "align/align_bpm_pattern.h"
#include "align/align_bpm_distance.h"
#include "align/align.h"
#include "matches/matches.h"
#include "matches/matches_cigar.h"

/*
 * BPM. Compute BPM-DP-Matrix
 *   @align_input->bpm_pattern
 *   @align_input->text
 *   @align_input->text_length
 */
void align_bpm_compute_matrix(
    match_align_input_t* const align_input,
    uint64_t max_distance,
    bpm_align_matrix_t* const bpm_align_matrix,
    mm_stack_t* const mm_stack) {
  // Parameters
  const bpm_pattern_t* const bpm_pattern = align_input->bpm_pattern;
  uint8_t* const text = align_input->text;
  const uint64_t text_length = align_input->text_length;
  // Pattern variables
  const uint64_t* PEQ = bpm_pattern->PEQ;
  const uint64_t num_words64 = bpm_pattern->pattern_num_words64;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  int64_t* const score = bpm_pattern->score;
  const int64_t* const init_score = bpm_pattern->init_score;
  // Allocate auxiliary matrix
  const uint64_t aux_matrix_size = num_words64*UINT64_SIZE*(text_length+1); /* (+1 base-column) */
  uint64_t* const Pv = (uint64_t*)mm_stack_malloc(mm_stack,aux_matrix_size);
  uint64_t* const Mv = (uint64_t*)mm_stack_malloc(mm_stack,aux_matrix_size);
  bpm_align_matrix->Mv = Mv;
  bpm_align_matrix->Pv = Pv;
  // Initialize search
  if (max_distance >= bpm_pattern->pattern_length) {
    max_distance = bpm_pattern->pattern_length-1; // Correct max-distance
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
}
/*
 * BPM. Recover CIGAR from a matching string
 *   @align_input->key
 *   @align_input->bpm_pattern
 *   @align_input->text
 *   @bpm_align_matrix->Pv
 *   @bpm_align_matrix->Mv
 *   @bpm_align_matrix->min_score
 *   @match_alignment->match_position (Adjusted)
 */
void align_bpm_backtrace_matrix(
    match_align_input_t* const align_input,
    const bool left_gap_alignment,
    bpm_align_matrix_t* const bpm_align_matrix,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector) {
  // Parameters
  const uint8_t* const key = align_input->key;
  const bpm_pattern_t* const bpm_pattern = align_input->bpm_pattern;
  const uint64_t pattern_length = bpm_pattern->pattern_length;
  uint8_t* const text = align_input->text;
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
    const bool deletion = Pv[bdp_idx] & mask;
    const bool insertion = Mv[(bdp_idx-num_words64)] & mask;
    const bool match = text[h]==key[v];
    cigar_t operation;
    if (left_gap_alignment) {
      if (deletion) {
        operation = cigar_del;
      } else if (insertion) {
        operation = cigar_ins;
      } else if (match) {
        operation = cigar_match;
      } else {
        operation = cigar_mismatch;
      }
    } else {
      if (match) {
        operation = cigar_match;
      } else if (deletion) {
        operation = cigar_del;
      } else if (insertion) {
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
 *   @align_input->key
 *   @align_input->bpm_pattern
 *   @align_input->text
 *   @align_input->text_length
 *   @match_alignment->match_position (Adjusted)
 */
void align_bpm_match(
    match_align_input_t* const align_input,
    const uint64_t max_distance,
    const bool left_gap_alignment,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_stack_t* const mm_stack) {
  // Fill Matrix (Pv,Mv)
  mm_stack_push_state(mm_stack); // Save stack state
  bpm_align_matrix_t bpm_align_matrix;
  align_bpm_compute_matrix(align_input,max_distance,&bpm_align_matrix,mm_stack);
  // Set distance
  match_alignment->score = bpm_align_matrix.min_score;
  if (bpm_align_matrix.min_score == ALIGN_DISTANCE_INF) {
    match_alignment->cigar_length = 0;
    mm_stack_pop_state(mm_stack); // Free
    return;
  }
  // Backtrace and generate CIGAR
  align_bpm_backtrace_matrix(align_input,left_gap_alignment,&bpm_align_matrix,match_alignment,cigar_vector);
  // Free
  mm_stack_pop_state(mm_stack);
}
