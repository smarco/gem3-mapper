/*
 * PROJECT: GEMMapper
 * FILE: bpm_align.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 * TODO
 *   - Block can surpass 256->uint8_t size
 */

#include "pattern.h"
#include "bpm_align.h"
#include "swg_align.h"
#include "filtering_region.h"
#include "match_align_dto.h"
#include "matches.h"

/*
 * Checks
 */
#define BPM_PATTERN_CHECK(bpm_pattern) GEM_CHECK_NULL(bpm_pattern->PEQ)

/*
 * Constants
 */
// Mode
//#define BPM_TILED
// Debug
#define DEBUG_BPM_TILED false
// Constants
#define BMP_W64_LENGTH UINT64_LENGTH
#define BMP_W64_ONES   UINT64_MAX
#define BMP_W64_MASK   (1ull<<63)
// Pattern Accessors
#define BPM_PATTERN_PEQ_IDX(word_pos,encoded_character)   ((word_pos*DNA__N_RANGE)+(encoded_character))
#define BPM_PATTERN_BDP_IDX(position,num_words,word_pos)  ((position)*(num_words)+(word_pos))

GEM_INLINE void bpm_pattern_compile_chunks(
    bpm_pattern_t* const bpm_pattern,uint8_t* const pattern,
    const uint64_t pattern_length,const uint64_t max_error,mm_stack_t* const mm_stack) {
  // Init BPM chunks
  const uint64_t words_per_chunk = DIV_CEIL(max_error,BPM_ALIGN_WORD_LENGTH);
  const uint64_t num_chunks = DIV_CEIL(bpm_pattern->pattern_num_words,words_per_chunk);
  bpm_pattern->words_per_chunk = words_per_chunk;
  bpm_pattern->num_pattern_chunks = num_chunks;
  bpm_pattern->bpm_pattern_chunks = mm_stack_calloc(mm_stack,num_chunks,bpm_pattern_t,false);
  uint64_t offset_words = 0, i;
  for (i=0;i<num_chunks;++i) {
    // Initialize pattern-chunk variables
    bpm_pattern_t* const bpm_pattern_chunk = bpm_pattern->bpm_pattern_chunks + i;
    bpm_pattern_chunk->P = bpm_pattern->P;
    bpm_pattern_chunk->M = bpm_pattern->M;
    bpm_pattern_chunk->PEQ = bpm_pattern->PEQ + offset_words*DNA__N_RANGE;
    bpm_pattern_chunk->pattern_num_words = words_per_chunk;
    bpm_pattern_chunk->level_mask = bpm_pattern->level_mask + offset_words;
    bpm_pattern_chunk->score = bpm_pattern->score + offset_words;
    bpm_pattern_chunk->init_score = bpm_pattern->init_score + offset_words;
    bpm_pattern_chunk->pattern_left = NULL;
    offset_words += words_per_chunk;
  }
}
GEM_INLINE void bpm_pattern_compile(
    bpm_pattern_t* const bpm_pattern,uint8_t* const pattern,
    const uint64_t pattern_length,const uint64_t max_error,mm_stack_t* const mm_stack) {
  GEM_CHECK_NULL(bpm_pattern);
  GEM_CHECK_NULL(pattern);
  GEM_CHECK_ZERO(pattern_length);
  // Calculate dimensions
  const uint64_t word_length = BPM_ALIGN_WORD_LENGTH;
  const uint64_t word_size = BPM_ALIGN_WORD_SIZE;
  const uint64_t pattern_num_words = DIV_CEIL(pattern_length,BPM_ALIGN_WORD_LENGTH);
  const uint64_t PEQ_length = pattern_num_words*word_length;
  const uint64_t pattern_mod = pattern_length%word_length;
  // Init fields
  bpm_pattern->pattern_word_size = word_size;
  bpm_pattern->pattern_length = pattern_length;
  bpm_pattern->pattern_num_words = pattern_num_words;
  bpm_pattern->pattern_mod = pattern_mod;
  bpm_pattern->PEQ_length = PEQ_length;
  // Allocate memory
  const uint64_t aux_vector_size = pattern_num_words*word_size;
  const uint64_t PEQ_size = DNA__N_RANGE*aux_vector_size;
  const uint64_t score_size = pattern_num_words*UINT64_SIZE;
  bpm_pattern->PEQ = mm_stack_malloc(mm_stack,PEQ_size);
  bpm_pattern->P = mm_stack_malloc(mm_stack,aux_vector_size);
  bpm_pattern->M = mm_stack_malloc(mm_stack,aux_vector_size);
  bpm_pattern->level_mask = mm_stack_malloc(mm_stack,aux_vector_size);
  bpm_pattern->score = mm_stack_malloc(mm_stack,score_size);
  bpm_pattern->init_score = mm_stack_malloc(mm_stack,score_size);
  bpm_pattern->pattern_left = mm_stack_malloc(mm_stack,(pattern_num_words+1)*UINT64_SIZE);
  // Init PEQ
  memset(bpm_pattern->PEQ,0,PEQ_size);
  uint64_t i;
  for (i=0;i<pattern_length;++i) {
    const uint8_t enc_char = pattern[i];
    const uint64_t block = i/BPM_ALIGN_WORD_LENGTH;
    const uint64_t mask = 1ull<<(i%BPM_ALIGN_WORD_LENGTH);
    bpm_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block,enc_char)] |= mask;
  }
  for (;i<PEQ_length;++i) {
    const uint64_t block = i/BPM_ALIGN_WORD_LENGTH;
    const uint64_t mask = 1ull<<(i%BPM_ALIGN_WORD_LENGTH);
    uint64_t j;
    for (j=0;j<DNA__N_RANGE;++j) {
      bpm_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block,j)] |= mask;
    }
  }
  // Init auxiliary data
  uint64_t pattern_left = pattern_length;
  const uint64_t top = pattern_num_words-1;
  memset(bpm_pattern->level_mask,0,aux_vector_size);
  for (i=0;i<top;++i) {
    bpm_pattern->level_mask[i] = BMP_W64_MASK;
    bpm_pattern->init_score[i] = word_length;
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > word_length) ? pattern_left-word_length : 0;
  }
  for (;i<=pattern_num_words;++i) {
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > word_length) ? pattern_left-word_length : 0;
  }
  if (pattern_mod>0) {
    const uint64_t mask_shift = pattern_mod-1;
    bpm_pattern->level_mask[top] = 1ull<<(mask_shift);
    bpm_pattern->init_score[top] = pattern_mod;
  } else {
    bpm_pattern->level_mask[top] = BMP_W64_MASK;
    bpm_pattern->init_score[top] = word_length;
  }
#ifdef BPM_TILED
  // Init BPM chunks
  bpm_pattern_compile_chunks(bpm_pattern,pattern,pattern_length,max_error,mm_stack);
#endif
  // Init BPM-GPU Dimensions
  bpm_pattern->gpu_num_entries = DIV_CEIL(pattern_length,BPM_GPU_PATTERN_ENTRY_LENGTH);
  bpm_pattern->gpu_entries_per_chunk = DIV_CEIL(max_error,BPM_GPU_PATTERN_ENTRY_LENGTH);
  bpm_pattern->gpu_num_chunks = DIV_CEIL(bpm_pattern->gpu_num_entries,bpm_pattern->gpu_entries_per_chunk);
}
/*
 * Reset search functions
 */
GEM_INLINE void bpm_reset_search(
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
GEM_INLINE void bpm_reset_search_cutoff(
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
GEM_INLINE int8_t bpm_advance_block(
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
GEM_INLINE bool bpm_get_distance_raw(
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
GEM_INLINE bool bpm_get_distance_cutoff(
    const bpm_pattern_t* const bpm_pattern,
    const uint8_t* const text,const uint64_t text_length,
    uint64_t* const match_end_column,uint64_t* const distance,
    const uint64_t max_distance,const bool quick_abandon) {
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
      PROF_INC_COUNTER(GP_BPM_QUICK_ABANDON);
      *distance = ALIGN_DISTANCE_INF;
      return false;
    }
  }
  // Return results
  if (min_score!=ALIGN_DISTANCE_INF) {
    *distance = min_score;
    *match_end_column = min_score_column; // Exact
    return true;
  } else {
    *distance = ALIGN_DISTANCE_INF;
    return false;
  }
}
/*
 * BPM Tiled (bound)
 */
GEM_INLINE void bpm_get_distance_cutoff_tiled(
    bpm_pattern_t* const bpm_pattern,const uint8_t* const text,const uint64_t text_length,
    uint64_t* const levenshtein_distance,uint64_t* const levenshtein_match_end_column,
    const uint64_t max_error) {
#ifndef BPM_TILED
  PROF_START(GP_BPM_TILED);
  // BPM Cut-off + Quick-abandon
  bpm_get_distance_cutoff(bpm_pattern,text,text_length,
      levenshtein_match_end_column,levenshtein_distance,max_error,true);
  PROF_ADD_COUNTER(GP_BMP_TILED_NUM_TILES,1);
  PROF_ADD_COUNTER(GP_BMP_TILED_NUM_TILES_VERIFIED,1);
  PROF_STOP(GP_BPM_TILED);
#else
  PROF_START(GP_BPM_TILED);
  // Fetch pattern dimensions
  const uint64_t num_pattern_chunks = bpm_pattern->num_pattern_chunks;
  const uint64_t words_per_chunk =  bpm_pattern->words_per_chunk;
  PROF_ADD_COUNTER(GP_BMP_TILED_NUM_TILES,num_pattern_chunks);
  // Calculate tile dimensions
  pattern_tiled_t pattern_tiled;
  const bool pattern_can_align = pattern_tiled_init(&pattern_tiled,
      bpm_pattern->pattern_length,words_per_chunk*BPM_ALIGN_WORD_LENGTH,text_length,max_error);
  if (!pattern_can_align) {
    *levenshtein_distance = ALIGN_DISTANCE_INF;
    *levenshtein_match_end_column = ALIGN_COLUMN_INF; // FIXME Needed?
    PROF_STOP(GP_BPM_TILED);
    return;
  }
  // Initialize current tile variables
  uint64_t pattern_chunk, global_distance = 0, distance_link_tiles = 0;
  for (pattern_chunk=0;pattern_chunk<num_pattern_chunks;++pattern_chunk) {
    PROF_ADD_COUNTER(GP_BMP_TILED_NUM_TILES_VERIFIED,1);
    // BPM Cut-off
    bpm_get_distance_cutoff(bpm_pattern->bpm_pattern_chunks+pattern_chunk,
        text+pattern_tiled.tile_offset,pattern_tiled.tile_wide,
        &pattern_tiled.tile_match_column,&pattern_tiled.tile_distance,max_error,false);
    // Update global distance
    global_distance += pattern_tiled.tile_distance;
    if (global_distance > max_error) {
      *levenshtein_distance = ALIGN_DISTANCE_INF;
      *levenshtein_match_end_column = ALIGN_COLUMN_INF; // FIXME Needed?
      PROF_STOP(GP_BPM_TILED);
      return;
    }
    // Bound the joint distance by estimating the cost of connecting the path through the tiles
    distance_link_tiles += pattern_tiled_bound_matching_path(&pattern_tiled);
    // Calculate next tile
    pattern_tiled_calculate_next(&pattern_tiled);
  }
  *levenshtein_distance = global_distance + distance_link_tiles; // Bound
  if (*levenshtein_distance==0) {
    *levenshtein_match_end_column = pattern_tiled.prev_tile_match_position;
  } else { // No reduction, we use the whole text
    *levenshtein_match_end_column = text_length-1;
    // if (*levenshtein_distance > max_filtering_error) { // FIXME: Remove me
    //  *levenshtein_distance = max_filtering_error; // Adjust bound overestimations
    // }
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_BPM_TILED) {
    uint64_t check_pos, check_distance;
    bpm_get_distance_cutoff(bpm_pattern,text,text_length,&check_pos,&check_distance,max_error,true);
    gem_cond_fatal_error_msg(check_distance < *levenshtein_distance,
        "Pattern-chunk verification failed (SUM(d(P_i))=%lu <= d(P)=%lu)",check_distance,*levenshtein_distance);
  }
  PROF_STOP(GP_BPM_TILED);
#endif
}
/*
 * BPM all matches
 */
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
    const uint64_t Peq = PEQ[BPM_PATTERN_PEQ_IDX(top_level,enc_char)]; \
    if (last_score<=max_distance && last<top && (MHin || (Peq & 1))) { \
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
#define BPM_ADD_FILTERING_REGION(text_trace_offset,begin_position,end_position,min_score_column,min_score,num_matches_found) { \
  filtering_region_t* filtering_region; \
  vector_alloc_new(filtering_regions,filtering_region_t,filtering_region); \
  /* State */ \
  filtering_region->status = filtering_region_accepted; \
  /* Text-trace */ \
  filtering_region->text_trace_offset = text_trace_offset; \
  /* Location */ \
  filtering_region->begin_position = begin_position; \
  filtering_region->end_position = end_position; \
  filtering_region->base_position_offset = 0; \
  /* Regions Matching */ \
  filtering_region->match_scaffold.scaffold_regions = NULL; \
  filtering_region->match_scaffold.num_scaffold_regions = 0; \
  filtering_region->match_scaffold.scaffolding_coverage = 0; \
  /* Alignment distance */ \
  filtering_region->align_distance = min_score; \
  filtering_region->align_match_end_column = min_score_column; \
  /* Increment the number of matches found */ \
  ++num_matches_found; \
}
GEM_INLINE uint64_t bpm_search_all(
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
        if (current_score < min_score)  { // Founded match, Update minimum
          min_score_column = text_position;
          min_score = current_score;
          if (current_score==0) { // Don't try to optimize, carry on
            match_found = false;
            BPM_ADD_FILTERING_REGION(text_trace_offset,begin_position,end_position,min_score_column,min_score,num_matches_found);
          } else {
            match_found = true;
            opt_steps_left = max_distance; // Setup optimization steps
          }
        }
      }
    } else {
      if (top_level==num_words && current_score<min_score) { // Update minimum
        min_score_column = text_position;
        min_score = current_score;
      }
      if (opt_steps_left==0) {
        match_found = false;
        BPM_ADD_FILTERING_REGION(text_trace_offset,begin_position,end_position,min_score_column,min_score,num_matches_found);
      } else {
        --opt_steps_left;
      }
    }
    // Quick abandon
    if (min_score==ALIGN_DISTANCE_INF && current_score+pattern_left[top_level] > text_left+max_distance) {
      // Quick abandon, it doesn't match (bounded by best case scenario)
      PROF_INC_COUNTER(GP_BPM_ALL_QUICK_ABANDON);
      break;
    }
  }
  if (match_found) {
    BPM_ADD_FILTERING_REGION(text_trace_offset,begin_position,end_position,min_score_column,min_score,num_matches_found);
  }
  PROF_INC_COUNTER(GP_BPM_ALL_MATCHES_FOUND);
  PROF_STOP(GP_BPM_ALL);
  return num_matches_found;
}
/*
 * BPM. Compute BPM-DP-Matrix
 *   @align_input->key
 *   @align_input->bpm_pattern
 *   @align_input->text
 *   @align_input->text_length
 */
GEM_INLINE void bpm_align_compute_matrix(
    match_align_input_t* const align_input,const uint64_t max_distance,
    bpm_align_matrix_t* const bpm_align_matrix,mm_stack_t* const mm_stack) {
  // Parameters
  const bpm_pattern_t* const bpm_pattern = align_input->bpm_pattern;
  uint8_t* const text = align_input->text;
  const uint64_t text_length = align_input->text_length;
  // Pattern variables
  const uint64_t* PEQ = bpm_pattern->PEQ;
  const uint64_t num_words = bpm_pattern->pattern_num_words;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  int64_t* const score = bpm_pattern->score;
  const int64_t* const init_score = bpm_pattern->init_score;
  // Allocate auxiliary matrix
  const uint64_t aux_matrix_size = num_words*bpm_pattern->pattern_word_size*(text_length+1); /* (+1 base-column) */
  uint64_t* const Pv = (uint64_t*)mm_stack_malloc(mm_stack,aux_matrix_size);
  uint64_t* const Mv = (uint64_t*)mm_stack_malloc(mm_stack,aux_matrix_size);
  bpm_align_matrix->Mv = Mv;
  bpm_align_matrix->Pv = Pv;
  // Initialize search
  const uint64_t max_distance__1 = max_distance+1;
  const uint8_t top = num_words-1;
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
      const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words,i);
      const uint64_t next_bdp_idx = bdp_idx+num_words;
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
    if (gem_expect_false(score[last]<=max_distance__1)) {
      const uint64_t last_score = score[last]+(MHin-PHin);
      const uint64_t Peq = PEQ[BPM_PATTERN_PEQ_IDX(top_level,enc_char)];
      if (last_score<=max_distance && last<top && (MHin || (Peq & 1))) {
        // Init block V
        const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(text_position,num_words,top_level);
        const uint64_t next_bdp_idx = bdp_idx+num_words;
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
    if (top_level==num_words && current_score<=max_distance) {
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
 *   @align_input->key_length
 *   @align_input->bpm_pattern
 *   @align_input->text
 *   @bpm_align_matrix->Pv
 *   @bpm_align_matrix->Mv
 *   @bpm_align_matrix->min_score
 *   @match_alignment->match_position (Adjusted)
 */
GEM_INLINE void bpm_align_backtrack_matrix(
    match_align_input_t* const align_input,const bool left_gap_alignment,
    bpm_align_matrix_t* const bpm_align_matrix,match_alignment_t* const match_alignment,
    vector_t* const cigar_vector) {
  // Parameters
  const uint8_t* const key = align_input->key;
  const uint64_t key_length = align_input->key_length;
  const bpm_pattern_t* const bpm_pattern = align_input->bpm_pattern;
  uint8_t* const text = align_input->text;
  const uint64_t* const Pv = bpm_align_matrix->Pv;
  const uint64_t* const Mv = bpm_align_matrix->Mv;
  // Allocate CIGAR string memory (worst case)
  match_alignment->cigar_offset = vector_get_used(cigar_vector); // Set CIGAR offset
  vector_reserve_additional(cigar_vector,MIN(key_length,2*bpm_align_matrix->min_score+1)); // Reserve
  cigar_element_t* cigar_buffer = vector_get_free_elm(cigar_vector,cigar_element_t); // Sentinel
  cigar_element_t* const cigar_buffer_base = cigar_buffer;
  cigar_buffer->type = cigar_null; // Trick
  // Retrieve the alignment. Store the match
  const uint64_t num_words = bpm_pattern->pattern_num_words;
  int64_t match_effective_length = key_length;
  int64_t h = bpm_align_matrix->min_score_column;
  int64_t v = key_length - 1;
  while (v >= 0 && h >= 0) {
    const uint8_t block = v / UINT64_LENGTH;
    const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(h+1,num_words,block);
    const uint64_t mask = 1L << (v % UINT64_LENGTH);
    // Select CIGAR operation
    const bool deletion = Pv[bdp_idx] & mask;
    const bool insertion = Mv[(bdp_idx-num_words)] & mask;
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
        matches_cigar_buffer_add_cigar_element(&cigar_buffer,cigar_del,1,NULL); // Deletion <-1>@v
        --v; --match_effective_length;
        break;
      case cigar_ins:
        matches_cigar_buffer_add_cigar_element(&cigar_buffer,cigar_ins,1,text+v); // Insertion <+1>@v
        --h; ++match_effective_length;
        break;
      case cigar_mismatch:
        if (cigar_buffer->type!=cigar_null) ++(cigar_buffer);
        cigar_buffer->type = cigar_mismatch;
        cigar_buffer->mismatch = text[h]; // Mismatch
        --h; --v;
        break;
      case cigar_match:
        matches_cigar_buffer_add_cigar_element(&cigar_buffer,cigar_match,1,NULL); // Match
        --h; --v;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  if (v >= 0) {
    matches_cigar_buffer_add_cigar_element(&cigar_buffer,cigar_del,v+1,NULL); // <-(@v+1)>@v
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
GEM_INLINE void bpm_align_match(
    match_align_input_t* const align_input,const uint64_t max_distance,
    const bool left_gap_alignment,match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,mm_stack_t* const mm_stack) {
  // Fill Matrix (Pv,Mv)
  mm_stack_push_state(mm_stack); // Save stack state
  bpm_align_matrix_t bpm_align_matrix;
  bpm_align_compute_matrix(align_input,max_distance,&bpm_align_matrix,mm_stack);
  // Set distance
  match_alignment->score = bpm_align_matrix.min_score;
  if (bpm_align_matrix.min_score == ALIGN_DISTANCE_INF) {
    mm_stack_pop_state(mm_stack,false); // Free
    return;
  }
  // Backtrace and generate CIGAR
  bpm_align_backtrack_matrix(align_input,left_gap_alignment,&bpm_align_matrix,match_alignment,cigar_vector);
  // Free
  mm_stack_pop_state(mm_stack,false);
}
