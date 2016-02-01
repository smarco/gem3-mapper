/*
 * PROJECT: GEMMapper
 * FILE: align_bpm_distance.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 * TODO
 *   - Block can surpass 256->uint8_t size
 */

#include "align_bpm_distance.h"
#include "align.h"
#include "pattern.h"
#include "filtering_region.h"

/*
 * Mode/Debug
 */
#define DEBUG_BPM_TILED false

/*
 * Reset search functions
 */
void bpm_reset_search(
    const uint64_t num_words,uint64_t* const P,uint64_t* const M,
    int64_t* const score,const int64_t* const init_score) {
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
    uint8_t* const top_level,uint64_t* const P,uint64_t* const M,int64_t* const score,
    const int64_t* const init_score,const uint64_t max_distance) {
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
 * BMP
 */
bool bpm_compute_edit_distance_raw(
    bpm_pattern_t* const bpm_pattern,const uint8_t* const text,const uint64_t text_length,
    uint64_t* const position,uint64_t* const distance) {
  // Pattern variables
  const uint64_t* PEQ = bpm_pattern->PEQ;
  const uint64_t num_words = bpm_pattern->pattern_num_words;
  uint64_t* const P = bpm_pattern->P;
  uint64_t* const M = bpm_pattern->M;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  int64_t* const score = bpm_pattern->score;
  const int64_t* const init_score = bpm_pattern->init_score;

  // Initialize search
  uint64_t min_score = ALIGN_DISTANCE_INF, min_score_column = ALIGN_COLUMN_INF;
  bpm_reset_search(num_words,P,M,score,init_score);

  // Advance in DP-bit_encoded matrix
  uint64_t text_position;
  for (text_position=0;text_position<text_length;++text_position) {
    // Fetch next character
    const uint8_t enc_char = text[text_position];

    // Advance all blocks
    int8_t carry;
    uint64_t i;
    for (i=0,carry=0;i<num_words;++i) {
      uint64_t* const Py = P+i;
      uint64_t* const My = M+i;
      carry = bpm_advance_block(PEQ[BPM_PATTERN_PEQ_IDX(i,enc_char)],level_mask[i],*Py,*My,carry+1,Py,My);
      score[i] += carry;
    }

    // Check match
    if (score[num_words-1] < min_score) {
      min_score_column = text_position;
      min_score = score[num_words-1];
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
    const uint8_t* const text,const uint64_t text_length,
    uint64_t* const match_distance,uint64_t* const match_column,
    const uint64_t max_distance,const bool quick_abandon) {
  PROF_START(GP_BPM_DISTANCE);
  PROF_ADD_COUNTER(GP_BPM_DISTANCE_CELLS,bpm_pattern->pattern_length*text_length);
  // Pattern variables
  const uint64_t* PEQ = bpm_pattern->PEQ;
  const uint64_t num_words = bpm_pattern->pattern_num_words;
  uint64_t* const P = bpm_pattern->P;
  uint64_t* const M = bpm_pattern->M;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  int64_t* const score = bpm_pattern->score;
  const int64_t* const init_score = bpm_pattern->init_score;
  const uint64_t* const pattern_left = bpm_pattern->pattern_left;

  // Initialize search
  const uint64_t max_distance__1 = max_distance+1;
  const uint8_t top = num_words-1;
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
    if (gem_expect_false(score[last]<=max_distance__1)) {
      const uint64_t last_score = score[last]+(MHin-PHin);
      const uint64_t Peq = PEQ[BPM_PATTERN_PEQ_IDX(top_level,enc_char)];
      if (last_score<=max_distance && last<top && (MHin || (Peq & 1))) {
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
    if (top_level==num_words && current_score<=max_distance) {
      if (current_score < min_score)  {
        min_score_column = text_position;
        min_score = current_score;
      }
    } else if (quick_abandon && min_score==ALIGN_DISTANCE_INF &&
        current_score+pattern_left[top_level] > text_left+max_distance) {
      // Quick abandon, it doesn't match (bounded by best case scenario)
      // TODO Test if (abandon_cond() && min_score!=ALIGN_DISTANCE_INF) return best_distace_score;
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
 * BPM Tiled (bound)
 */
void bpm_compute_edit_distance_tiled(
    bpm_pattern_t* const bpm_pattern,const uint8_t* const text,
    const uint64_t text_length,uint64_t* const match_distance,
    uint64_t* const match_column,const uint64_t max_error,
    mm_stack_t* const mm_stack) {
  PROF_START(GP_BPM_DISTANCE_TILED);
  // Fetch pattern dimensions
  const uint64_t num_pattern_tiles = bpm_pattern->num_pattern_tiles;
  const uint64_t words_per_tile =  bpm_pattern->words_per_tile;
  PROF_ADD_COUNTER(GP_BMP_DISTANCE_TILED_NUM_TILES,num_pattern_tiles);
  // Calculate tile dimensions
  pattern_tiled_t pattern_tiled;
  const bool pattern_can_align = pattern_tiled_init(&pattern_tiled,
      bpm_pattern->pattern_length,words_per_tile*BPM_ALIGN_WORD_LENGTH,text_length,max_error);
  if (!pattern_can_align) {
    *match_distance = ALIGN_DISTANCE_INF;
    *match_column = ALIGN_COLUMN_INF;
    PROF_STOP(GP_BPM_DISTANCE_TILED);
    return;
  }
  // Initialize current tile variables
  bpm_pattern_t* const bpm_pattern_tiles = bpm_pattern_compile_tiles(bpm_pattern,1,max_error,mm_stack);
  uint64_t tile_offset, global_distance = 0, distance_link_tiles = 0;
  for (tile_offset=0;tile_offset<num_pattern_tiles;++tile_offset) {
    PROF_ADD_COUNTER(GP_BMP_DISTANCE_TILED_NUM_TILES_VERIFIED,1);
    // BPM Cut-off
    bpm_compute_edit_distance(bpm_pattern_tiles+tile_offset,
        text+pattern_tiled.tile_offset,pattern_tiled.tile_wide,
        &pattern_tiled.tile_distance,&pattern_tiled.tile_match_column,
        max_error,false);
    // Update global distance
    global_distance += pattern_tiled.tile_distance;
    if (global_distance > max_error) {
      *match_distance = ALIGN_DISTANCE_INF;
      *match_column = ALIGN_COLUMN_INF;
      PROF_STOP(GP_BPM_DISTANCE_TILED);
      return;
    }
    // Bound the joint distance by estimating the cost of connecting the path through the tiles
    distance_link_tiles += pattern_tiled_bound_matching_path(&pattern_tiled);
    // Calculate next tile
    pattern_tiled_calculate_next(&pattern_tiled);
  }
  // Assign distance
  *match_distance = global_distance + distance_link_tiles; // Bound
  if (*match_distance==0) {
    *match_column = pattern_tiled.prev_tile_match_position;
  } else {
    if (*match_distance > max_error) {
      *match_distance = max_error;
    }
    *match_column = BOUNDED_ADDITION(pattern_tiled.prev_tile_match_position,*match_distance,text_length-1);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_BPM_TILED) {
    uint64_t check_pos, check_distance;
    bpm_compute_edit_distance(bpm_pattern,text,text_length,&check_distance,&check_pos,max_error,true);
    gem_cond_fatal_error_msg(check_distance < *match_distance,
        "BPM.Compute.Edit.Distance.Tiled :: Verification failed ("
        "Single.Tile.Distance=%"PRIu64" < Tiles.Bound.Distance=%"PRIu64")",
        check_distance,*match_distance);
  }
  PROF_STOP(GP_BPM_DISTANCE_TILED);
}
/*
 * BPM all matches
 */
uint64_t bpm_compute_edit_distance_all(
    const bpm_pattern_t* const bpm_pattern,vector_t* const filtering_regions,
    const uint64_t text_trace_offset,const uint64_t begin_position,
    const uint8_t* const text,const uint64_t text_length,const uint64_t max_distance) {
  PROF_START(GP_BPM_ALL);
  // Pattern variables
  const uint64_t* PEQ = bpm_pattern->PEQ;
  const uint64_t num_words = bpm_pattern->pattern_num_words;
  uint64_t* const P = bpm_pattern->P;
  uint64_t* const M = bpm_pattern->M;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  int64_t* const score = bpm_pattern->score;
  const int64_t* const init_score = bpm_pattern->init_score;
  const uint64_t* const pattern_left = bpm_pattern->pattern_left;
  // Initialize search
  const uint64_t max_distance__1 = max_distance+1;
  const uint8_t top = num_words-1;
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
      if (top_level==num_words && current_score<=max_distance) {
        min_score_column = text_position;
        min_score = current_score;
        if (current_score==0) { // Don't try to optimize (exact match)
          filtering_region_add(filtering_regions,text_trace_offset,
              begin_position,end_position,min_score,min_score_column);
          ++num_matches_found; // Increment the number of matches found
        } else {
          match_found = true;
          opt_steps_left = max_distance; // Setup optimization steps
        }
      }
    } else {
      if (top_level==num_words && current_score<min_score) { // Update minimum
        min_score_column = text_position;
        min_score = current_score;
      }
      if (opt_steps_left==0) {
        filtering_region_add(filtering_regions,text_trace_offset,
            begin_position,end_position,min_score,min_score_column);
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
    filtering_region_add(filtering_regions,text_trace_offset,
        begin_position,end_position,min_score,min_score_column);
    ++num_matches_found; // Increment the number of matches found
  }
  PROF_INC_COUNTER(GP_BPM_ALL_MATCHES_FOUND);
  PROF_STOP(GP_BPM_ALL);
  return num_matches_found;
}
