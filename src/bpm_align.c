/*
 * PROJECT: GEMMapper
 * FILE: bpm_align.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 * TODO
 *   - Block can surpass 256->uint8_t size
 */

#include "bpm_align.h"
#include "matches.h"

#define BPM_PATTERN_CHECK(bpm_pattern) GEM_CHECK_NULL(bpm_pattern->PEQ)

/*
 * (Re)alignment Basic: Dynamic Programming - LEVENSHTEIN/EDIT
 */
#define GT_DP(i,j) dp_array[(i)*pattern_len+(j)]
GEM_INLINE void align_levenshtein_dp_matrix_print(
    uint64_t* const dp_array,const uint64_t pattern_len,const uint64_t sequence_len,
    const uint64_t pattern_limit,const uint64_t sequence_limit) {
  uint64_t i, j;
  for (j=0;j<pattern_limit;++j) {
    for (i=0;i<sequence_limit;++i) {
      fprintf(stdout,"%02"PRIu64" ",GT_DP(i,j));
    }
    fprintf(stdout,"\n");
  }
  fprintf(stdout,"\n");
}
GEM_INLINE int64_t align_levenshtein_get_distance(
    const char* const pattern,const uint64_t pattern_length,
    const char* const sequence,const uint64_t sequence_length,
    const bool ends_free,uint64_t* const position) {
  GEM_CHECK_NULL(pattern); GEM_CHECK_ZERO(pattern_length);
  GEM_CHECK_NULL(sequence); GEM_CHECK_ZERO(sequence_length);
  // Allocate DP-matrix
  const uint64_t pattern_len = pattern_length+1;
  const uint64_t sequence_len = sequence_length+1;
  uint64_t* dp_array[2];
  dp_array[0] = mm_calloc(2*pattern_len,uint64_t,false);
  dp_array[1] = dp_array[0] + pattern_len;
  // Init DP-Matrix
  uint64_t min_val = UINT64_MAX, i_pos = UINT64_MAX;
  uint64_t i, j, idx_a=0, idx_b=0;
  for (j=0;j<pattern_len;++j) dp_array[0][j]=j;
  // Calculate DP-Matrix
  for (i=1;i<sequence_len;++i) {
    // Fix indexes
    idx_a = idx_b;
    idx_b = i % 2;
    // Fix first cell
    dp_array[idx_b][0] = (ends_free) ? 0 : dp_array[idx_a][0]+1;
    // Develop row
    for (j=1;j<pattern_len;++j) {
      const uint64_t ins = dp_array[idx_a][j]   + 1;
      const uint64_t del = dp_array[idx_b][j-1] + 1;
      const uint64_t sub = dp_array[idx_a][j-1] + ((sequence[i-1]==pattern[j-1]) ? 0 : 1);
      dp_array[idx_b][j] = MIN(sub,MIN(ins,del));
    }
    // Check last cell value
    if (ends_free && dp_array[idx_b][pattern_length] < min_val) {
      min_val = dp_array[idx_b][pattern_length];
      i_pos = i;
    }
  }
  // Return results & Free
  int64_t distance = INT64_MAX;
  if (ends_free) {
    *position = i_pos-1;
    distance = min_val;
  } else {
    *position = pattern_length;
    distance = dp_array[idx_b][pattern_length];
  }
  mm_free(dp_array[0]);
  return distance;
}
//GEM_INLINE int64_t gt_map_block_realign_levenshtein(
//    char* const pattern,const uint64_t pattern_length,
//    char* const sequence,const uint64_t sequence_length,const bool ends_free) {
//  GT_NULL_CHECK(pattern); GT_ZERO_CHECK(pattern_length);
//  GT_NULL_CHECK(sequence); GT_ZERO_CHECK(sequence_length);
//  // Allocate DP matrix
//  const uint64_t pattern_len = pattern_length+1;
//  const uint64_t sequence_len = sequence_length+1;
//  uint64_t* dp_array;
//  dp_array = mm_calloc(pattern_len*sequence_len,uint64_t,false);
//  // Init DP-Matrix
//  uint64_t min_val = UINT64_MAX, i_pos = UINT64_MAX;
//  uint64_t i, j;
//  for (i=0;i<sequence_len;++i) GT_DP(i,0)=(ends_free)?0:i;
//  for (j=0;j<pattern_len;++j) GT_DP(0,j)=j;
//  // Calculate DP-Matrix
//  for (i=1;i<sequence_len;++i) {
//    for (j=1;j<pattern_len;++j) {
//      const uint64_t ins = GT_DP(i-1,j) + 1;
//      const uint64_t del = GT_DP(i,j-1) + 1;
//      const uint64_t sub = GT_DP(i-1,j-1) + ((sequence[i-1]==pattern[j-1]) ? 0 : 1);
//      GT_DP(i,j) = GT_MIN(sub,GT_MIN(ins,del));
//    }
//    // Check last cell value
//    if (ends_free && GT_DP(i,pattern_length) < min_val) {
//      min_val = GT_DP(i,pattern_length);
//      i_pos = i;
//    }
//  }
//  // DEBUG
//  align_levenshtein_dp_matrix_print(dp_array,pattern_len,sequence_len,30,30);
//  // Return results & Free
//  const int64_t distance = (ends_free) ? min_val : GT_DP(i_pos,pattern_len-1);
//  mm_free(dp_array);
//  return distance;
//}

/*
 * Bit-compressed (Re)alignment
 *   BPM[BitParalellMyers] - Myers' Fast Bit-Vector algorithm (Levenshtein)
 */
// Constants
#define BMP_W64_LENGTH UINT64_LENGTH
#define BMP_W64_ONES   UINT64_MAX

#define BMP_INITIAL_BUFFER_SIZE ((16*DNA_RANGE+16*5)*64)
#define BMP_W8_LENGTH 8
#define BMP_W8_MASK (1<<7)

GEM_INLINE void bpm_pattern_compile(
    bpm_pattern_t* const bpm_pattern,const uint64_t word_size,
    uint8_t* const pattern,const uint64_t pattern_length,mm_stack_t* const mm_stack) {
  GEM_CHECK_NULL(bpm_pattern);
  GEM_CHECK_ZERO(word_size);
  GEM_CHECK_NULL(pattern);
  GEM_CHECK_ZERO(pattern_length);
  // Calculate dimensions
  const uint64_t word_length = word_size*BMP_W8_LENGTH;
  const uint64_t pattern_num_words = (pattern_length+(word_length-1))/word_length;
  const uint64_t PEQ_length = pattern_num_words*word_length;
  const uint64_t pattern_mod = pattern_length%word_length;
  // Init fields
  bpm_pattern->pattern_word_size = word_size;
  bpm_pattern->pattern_length = pattern_length;
  bpm_pattern->pattern_num_words = pattern_num_words;
  bpm_pattern->pattern_mod = pattern_mod;
  bpm_pattern->PEQ_length = PEQ_length;
  // Allocate memory
  const uint64_t PEQ_size = DNA__N_RANGE*pattern_num_words*word_size;
  const uint64_t aux_vector_size = pattern_num_words*word_size;
  const uint64_t score_size = pattern_num_words*UINT64_SIZE;
  bpm_pattern->PEQ = mm_stack_malloc(mm_stack,PEQ_size);
  bpm_pattern->P = mm_stack_malloc(mm_stack,aux_vector_size);
  bpm_pattern->M = mm_stack_malloc(mm_stack,aux_vector_size);
  bpm_pattern->level_mask = mm_stack_malloc(mm_stack,aux_vector_size);
  bpm_pattern->score = mm_stack_malloc(mm_stack,score_size);
  bpm_pattern->init_score = mm_stack_malloc(mm_stack,score_size);
  // Init PEQ
  memset(bpm_pattern->PEQ,0,PEQ_size);
  uint64_t i;
  for (i=0;i<pattern_length;++i) {
    const uint8_t enc_char = pattern[i];
    const uint64_t block = i/BMP_W8_LENGTH;
    const uint8_t mask = 1<<(i%BMP_W8_LENGTH);
    bpm_pattern->PEQ[BPM_PATTERN_PEQ_IDX(enc_char,block,pattern_num_words*word_size)] |= mask;
  }
  for (;i<PEQ_length;++i) {
    const uint64_t block = i/BMP_W8_LENGTH;
    const uint8_t mask = 1<<(i%BMP_W8_LENGTH);
    uint64_t j;
    for (j=0;j<DNA_RANGE;++j) {
      bpm_pattern->PEQ[BPM_PATTERN_PEQ_IDX(j,block,pattern_num_words*word_size)] |= mask;
    }
  }
  // Init auxiliary data
  const uint64_t top = pattern_num_words-1;
  memset(bpm_pattern->level_mask,0,aux_vector_size);
  for (i=0;i<top;++i) {
    bpm_pattern->level_mask[i*word_size+(word_size-1)] = BMP_W8_MASK;
    bpm_pattern->init_score[i] = word_length;
  }
  if (pattern_mod>0) {
    const uint64_t mask_shift = pattern_mod-1;
    const uint64_t block = mask_shift/BMP_W8_LENGTH;
    const uint8_t mask = 1<<(mask_shift%BMP_W8_LENGTH);
    bpm_pattern->level_mask[top*word_size+block] = mask;
    bpm_pattern->init_score[top] = pattern_mod;
  } else {
    bpm_pattern->level_mask[top*word_size+(word_size-1)] = BMP_W8_MASK;
    bpm_pattern->init_score[top] = word_length;
  }
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
GEM_INLINE void bpm_reset_search__cutoff(
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

  hout += T_hout_64[(Ph & mask)!=0][(Mh & mask)!=0]; // FIXME

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
GEM_INLINE bool bpm_get_distance(
    bpm_pattern_t* const bpm_pattern,const uint8_t* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance) {
  // Pattern variables
  const uint64_t* PEQ = (uint64_t*)bpm_pattern->PEQ;
  const uint64_t num_words = bpm_pattern->pattern_num_words;
  uint64_t* const P = (uint64_t*)bpm_pattern->P;
  uint64_t* const M = (uint64_t*)bpm_pattern->M;
  const uint64_t* const level_mask = (uint64_t*)bpm_pattern->level_mask;
  int64_t* const score = (int64_t*)bpm_pattern->score;
  const int64_t* const init_score = (int64_t*)bpm_pattern->init_score;

  // Initialize search
  uint64_t min_score = UINT64_MAX, min_score_position = UINT64_MAX;
  bpm_reset_search(num_words,P,M,score,init_score);

  // Advance in DP-bit_encoded matrix
  uint64_t sequence_position;
  for (sequence_position=0;sequence_position<sequence_length;++sequence_position) {
    // Fetch next character
    const uint8_t enc_char = sequence[sequence_position];

    // Advance all blocks
    int8_t carry;
    uint64_t i;
    for (i=0,carry=0;i<num_words;++i) {
      uint64_t* const Py = P+i;
      uint64_t* const My = M+i;
      carry = bpm_advance_block(
          PEQ[BPM_PATTERN_PEQ_IDX(enc_char,i,num_words)],level_mask[i],*Py,*My,carry+1,Py,My);
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
GEM_INLINE bool bpm_get_distance__cutoff(
    const bpm_pattern_t* const bpm_pattern,
    const uint8_t* const sequence,const uint64_t sequence_length,
    uint64_t* const position,uint64_t* const distance,const uint64_t max_distance) {
  // Pattern variables
  const uint64_t* PEQ = (uint64_t*)bpm_pattern->PEQ;
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
  bpm_reset_search__cutoff(&top_level,P,M,score,init_score,max_distance);

  // Advance in DP-bit_encoded matrix
  uint64_t sequence_position;
  for (sequence_position=0;sequence_position<sequence_length;++sequence_position) {
    // Fetch next character
    const uint8_t enc_char = sequence[sequence_position];

    // Advance all blocks
    int8_t carry;
    uint64_t i;
    for (i=0,carry=0;i<top_level;++i) {
      uint64_t* const Py = P+i;
      uint64_t* const My = M+i;
      carry = bpm_advance_block(
          PEQ[BPM_PATTERN_PEQ_IDX(enc_char,i,num_words)],level_mask[i],*Py,*My,carry+1,Py,My);
      score[i] += carry;
    }

    // Cut-off
    const uint8_t last = top_level-1;
    if ((score[last]-carry)<=max_distance && last<top &&
        ( (PEQ[BPM_PATTERN_PEQ_IDX(enc_char,top_level,num_words)] & 1) || (carry<0) )  ) {
      // Init block V
      P[top_level]=BMP_W64_ONES;
      M[top_level]=0;

      uint64_t* const Py = P+top_level;
      uint64_t* const My = M+top_level;
      score[top_level] = score[top_level-1] + init_score[top_level] - carry +
          bpm_advance_block(PEQ[BPM_PATTERN_PEQ_IDX(enc_char,top_level,num_words)],
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
//GEM_INLINE bool gt_map_block_bpm_get_distance_simd128_unfold64(
//    gt_bpm_pattern* const bpm_pattern,char* const sequence,const uint64_t sequence_length,
//    uint64_t* const position,uint64_t* const distance) {
//  // Pattern variables
//  const uint64_t* PEQ = (uint64_t*)bpm_pattern->PEQ;
//  const uint64_t num_words = bpm_pattern->pattern_num_words;
//  uint64_t* const P = (uint64_t*)bpm_pattern->P;
//  uint64_t* const M = (uint64_t*)bpm_pattern->M;
//  const uint64_t* const level_mask = (uint64_t*)bpm_pattern->level_mask;
//  int64_t* const score = (int64_t*)bpm_pattern->score;
//  const int64_t* const init_score = (int64_t*)bpm_pattern->init_score;
//
//  // Initialize search
//  uint64_t min_score = UINT64_MAX, min_score_position = UINT64_MAX;
//  uint64_t i;
//  P[0]=UINT64_MAX;
//  P[1]=UINT64_MAX;
//  M[0]=0;
//  M[1]=0;
//  score[0] = init_score[0];
//  for (i=1;i<num_words;++i) {
//    const uint64_t idx = i<<1;
//    P[idx]=UINT64_MAX;
//    P[idx+1]=UINT64_MAX;
//    M[idx]=0;
//    M[idx+1]=0;
//    score[i] = score[i-1] + init_score[i];
//  }
//
//  // Advance in DP-bit_encoded matrix
//  uint64_t sequence_position;
//  for (sequence_position=0;sequence_position<sequence_length;++sequence_position) {
//    // Fetch next character
//    const uint8_t enc_char = gt_cdna_encode(sequence[sequence_position]);
//
//    // Advance all blocks
//    uint64_t PHin=0,MHin=0,PHout,MHout;
//    for (i=0;i<num_words;++i) {
//      uint64_t* const P_mem_addr = P+(i<<1);
//      uint64_t* const M_mem_addr = M+(i<<1);
//      const uint64_t* const mask_mem_addr = (level_mask+(i<<1));
//      const uint64_t* const Eq_mem_addr = PEQ+((enc_char*num_words+i)<<1);
//
//      /*
//       * Compute Block
//       */
//      // Xv = Eq | Mv;
//      uint64_t Xv[2];
//      Xv[0] = Eq_mem_addr[0] | M_mem_addr[0];
//      Xv[1] = Eq_mem_addr[1] | M_mem_addr[1];
//      // _Eq = Eq | MHin;
//      uint64_t _Eq[2];
//      _Eq[0] = Eq_mem_addr[0] | MHin;
//      _Eq[1] = Eq_mem_addr[1];
//      // Xh = (((_Eq & Pv) + Pv) ^ Pv) | _Eq;
//      uint64_t Xh[2];
//      Xh[0] = (_Eq[0] & P_mem_addr[0]);
//      Xh[1] = (_Eq[1] & P_mem_addr[1]);
//      Xh[0] = Xh[0] + P_mem_addr[0];
//      Xh[1] = Xh[1] + P_mem_addr[1] + (Xh[0]<P_mem_addr[0]);
//      Xh[0] = (Xh[0] ^ P_mem_addr[0]) | _Eq[0];
//      Xh[1] = (Xh[1] ^ P_mem_addr[1]) | _Eq[1];
//
//      /* Calculate Hout */
//      // Ph = Mv | ~(Xh | Pv);
//      uint64_t Ph[2];
//      Ph[0] = M_mem_addr[0] | ~(Xh[0] | P_mem_addr[0]);
//      Ph[1] = M_mem_addr[1] | ~(Xh[1] | P_mem_addr[1]);
//      // Mh = Pv & Xh;
//      uint64_t Mh[2];
//      Mh[0] = P_mem_addr[0] & Xh[0];
//      Mh[1] = P_mem_addr[1] & Xh[1];
//
//      /* Account Hout that propagates for the next block */
//      // PHout = (Ph & mask)!=0;
//      // MHout = (Mh & mask)!=0;
//      PHout = ((Ph[0] & mask_mem_addr[0]) | (Ph[1] & mask_mem_addr[1])) !=0;
//      MHout = ((Mh[0] & mask_mem_addr[0]) | (Mh[1] & mask_mem_addr[1])) !=0;
//
//      /* Hout become the Hin of the next cell */
//      // Ph <<= 1;
//      // Mh <<= 1;
//      /* Account Hin coming from the previous block */
//      // Ph |= PHin;
//      // Mh |= MHin;
//      Ph[1] = (Ph[1]<<1) | ((Ph[0] & GT_BMP_W64_MASK)!=0);
//      Ph[0] = (Ph[0]<<1) | PHin;
//      Mh[1] = (Mh[1]<<1) | ((Mh[0] & GT_BMP_W64_MASK)!=0);
//      Mh[0] = (Mh[0]<<1) | MHin;
//
//      /* Finally, generate the Vout */
//      // Pv = Mh | ~(Xv | Ph);
//      P_mem_addr[0] = Mh[0] | ~(Xv[0] | Ph[0]);
//      P_mem_addr[1] = Mh[1] | ~(Xv[1] | Ph[1]);
//      // Mv = Ph & Xv
//      M_mem_addr[0] = Ph[0] & Xv[0];
//      M_mem_addr[1] = Ph[1] & Xv[1];
//
//      /* Adjust score and swap propagate Hv */
//      score[i] += PHout-MHout;
//      PHin=PHout;
//      MHin=MHout;
//    }
//
//    // Check match
//    if (score[num_words-1] < min_score) {
//      min_score_position = sequence_position;
//      min_score = score[num_words-1];
//    }
//  }
//  // Return results
//  if (min_score!=UINT64_MAX) {
//    *distance = min_score;
//    *position = min_score_position;
//    return true;
//  } else {
//    *distance = UINT64_MAX;
//    *position = UINT64_MAX;
//    return false;
//  }
//}



/*
 * Recover CIGAR from a matching string
 */
#define BPM_ALIGN_ADD_CIGAR_ELEMENT(cigar_buffer,cigar_element_type,element_length) \
  if (cigar_buffer->type == cigar_element_type) { \
    cigar_buffer->length += element_length; \
  } else { \
    if (cigar_buffer->type!=cigar_null) ++(cigar_buffer); \
    cigar_buffer->type = cigar_element_type; \
    cigar_buffer->length = element_length; \
  }
GEM_INLINE void bpm_align_match(
    const uint8_t* const key,const bpm_pattern_t* const bpm_pattern,
    const uint8_t* const sequence,uint64_t* const match_position,
    const uint64_t matching_distance,const uint64_t matching_column,
    vector_t* const cigar_vector,uint64_t* const cigar_vector_offset,uint64_t* const cigar_length,
    mm_stack_t* const mm_stack) {
  // Pattern variables
  const uint64_t* PEQ = (uint64_t*)bpm_pattern->PEQ;
  const uint64_t num_words = bpm_pattern->pattern_num_words;
  const uint64_t* const level_mask = (uint64_t*)bpm_pattern->level_mask;
  int64_t* const score = (int64_t*)bpm_pattern->score;
  const int64_t* const init_score = (int64_t*)bpm_pattern->init_score;
  // Allocate auxiliary matrix
  const uint64_t aux_matrix_size =
      bpm_pattern->pattern_num_words*bpm_pattern->pattern_word_size*(matching_column+2); /* (+1 length) (+1 base-column)*/
  mm_stack_push_state(mm_stack); // Save stack state
  uint64_t* const Pv = (uint64_t*)mm_stack_malloc(mm_stack,aux_matrix_size);
  uint64_t* const Mv = (uint64_t*)mm_stack_malloc(mm_stack,aux_matrix_size);
  // Initialize search
  const uint8_t top = num_words-1;
  uint8_t top_level;
  bpm_reset_search__cutoff(&top_level,Pv,Mv,score,init_score,matching_distance);

  // Advance in DP-bit_encoded matrix
  uint64_t sequence_position;
  for (sequence_position=0;sequence_position<=matching_column;++sequence_position) {
    // Fetch next character
    const uint8_t enc_char = sequence[sequence_position];
    // Advance all blocks
    int8_t carry;
    uint64_t i;
    for (i=0,carry=0;i<top_level;++i) {
      const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(sequence_position,i,num_words);
      const uint64_t next_bdp_idx = bdp_idx+num_words;
      uint64_t* const Pv_in = Pv + bdp_idx;
      uint64_t* const Mv_in = Mv + bdp_idx;
      uint64_t* const Pv_out = Pv + next_bdp_idx;
      uint64_t* const Mv_out = Mv + next_bdp_idx;
      carry = bpm_advance_block(
          PEQ[BPM_PATTERN_PEQ_IDX(enc_char,i,num_words)],
          level_mask[i],*Pv_in,*Mv_in,carry+1,Pv_out,Mv_out);
      score[i] += carry;
    }
    // Cut-off
    const uint8_t last = top_level-1;
    if ((score[last]-carry)<=matching_distance && last<top &&
        ( (PEQ[BPM_PATTERN_PEQ_IDX(enc_char,top_level,num_words)] & 1) || (carry<0) )  ) {
      // Init block V
      const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(sequence_position,top_level,num_words);
      const uint64_t next_bdp_idx = bdp_idx+num_words;
      Pv[bdp_idx]=BMP_W64_ONES;
      Mv[bdp_idx]=0;
      // Advance block
      uint64_t* const Pv_in = Pv + bdp_idx;
      uint64_t* const Mv_in = Mv + bdp_idx;
      uint64_t* const Pv_out = Pv + next_bdp_idx;
      uint64_t* const Mv_out = Mv + next_bdp_idx;
      score[top_level] = score[top_level-1] + init_score[top_level] - carry +
          bpm_advance_block(PEQ[BPM_PATTERN_PEQ_IDX(enc_char,top_level,num_words)],
              level_mask[top_level],*Pv_in,*Mv_in,carry+1,Pv_out,Mv_out);
      ++top_level;
    } else {
      while (score[top_level-1]>matching_distance+init_score[top_level-1]) {
        --top_level;
      }
    }
  }
  gem_fatal_check_msg(top_level!=num_words || score[top_level-1]!=matching_distance,
      "BPM-Align. No match found at given position/column");
  /*
   * Retrieve the alignment. Store the match (Backtrace and generate CIGAR)
   */
  // Allocate CIGAR string memory (worst case)
  const uint64_t pattern_length = bpm_pattern->pattern_length;
  *cigar_vector_offset = vector_get_used(cigar_vector); // Set CIGAR offset
  vector_reserve_additional(cigar_vector,pattern_length); // Reserve
  cigar_element_t* cigar_buffer = vector_get_free_elm(cigar_vector,cigar_element_t); // Sentinel
  cigar_element_t* const cigar_buffer_base = cigar_buffer;
  cigar_buffer->type = cigar_null; // Trick
  // Start Backtrace
  int64_t h = matching_column;
  int64_t v = pattern_length - 1;
  while (v >= 0 && h >= 0) {
    const uint8_t block = v / UINT64_LENGTH;
    const uint64_t bdp_idx = BPM_PATTERN_BDP_IDX(h+1,block,num_words);
    const uint64_t mask = 1L << (v % UINT64_LENGTH);
    if (Pv[bdp_idx] & mask) {
      BPM_ALIGN_ADD_CIGAR_ELEMENT(cigar_buffer,cigar_del,1); // Deletion <-1>@v
      --v;
    } else if (Mv[(bdp_idx-num_words)] & mask) {
      BPM_ALIGN_ADD_CIGAR_ELEMENT(cigar_buffer,cigar_ins,1); // Insertion <+1>@v
      --h;
    } else if (sequence[h] != key[v]) {
      // Mismatch
      if (cigar_buffer->type!=cigar_null) ++(cigar_buffer);
      cigar_buffer->type = cigar_mismatch;
      cigar_buffer->mismatch = sequence[h];
      --h; --v;
    } else {
      BPM_ALIGN_ADD_CIGAR_ELEMENT(cigar_buffer,cigar_match,1); // Match
      --h; --v;
    }
  }
  if (v >= 0) {
    BPM_ALIGN_ADD_CIGAR_ELEMENT(cigar_buffer,cigar_del,v+1); // <-(@v+1)>@v
  }
  if (h >= 0) {
    *match_position += h+1; // We need to correct the matching_position
  }
  // Set CIGAR buffer used
  if (cigar_buffer->type!=cigar_null) ++(cigar_buffer);
  const uint64_t num_cigar_elements = cigar_buffer - cigar_buffer_base;
  vector_add_used(cigar_vector,num_cigar_elements);
  *cigar_length = num_cigar_elements; // Set CIGAR length
  // Reverse CIGAR Elements
  if (num_cigar_elements > 0) {
    const uint64_t middle_point = num_cigar_elements/2;
    uint64_t i;
    for (i=0;i<middle_point;++i) {
      SWAP(cigar_buffer_base[i],cigar_buffer_base[num_cigar_elements-i-1]);
    }
  }
  mm_stack_pop_state(mm_stack,false); // Free
//  // Check // TODO
//  gem_check(!align_check(
//      key,bpm_pattern->pattern_length,
//      sequence+(h+1),matching_column-(h+1),
//      cigar_buffer,*cigar_buffer_offset,*cigar_length);
}










//// Myers DEBUG
//bool fmi_matches_check_alignment(
//    const ch_t* const key,const uint64_t key_length,
//    const ch_t* const text,uint64_t alg_length,
//    mismatch* misms,const int64_t num_misms) {
//  register bool correct = true;
//  register int64_t p_ref = 0, p_read = 0, misms_pos=0;
//  while (correct && p_read<key_length && p_ref<alg_length) {
//    if (misms_pos<num_misms && misms[misms_pos].position == p_read) {
//      if (misms[misms_pos].mismatch>=256) {
//        register const uint64_t tmp=misms[misms_pos].mismatch/256;
//        register const bool is_del = tmp%2;
//        register const uint64_t size = tmp/2;
//        if (!is_del) p_ref+=size;
//        else p_read+=size;
//      } else {
//        correct = (misms[misms_pos].mismatch == text[p_ref] && key[p_read] != text[p_ref]);
//        p_read++; p_ref++;
//      }
//      ++misms_pos;
//    } else {
//      correct = key[p_read] == text[p_ref];
//      p_read++; p_ref++;
//    }
//  }
//  while (correct && p_read<key_length) {
//    if (misms_pos>=0 && misms[misms_pos].position == p_read) {
//      if (misms[misms_pos].mismatch>=256) {
//        register const uint64_t tmp=misms[misms_pos].mismatch/256;
//        register const bool is_del = tmp%2;
//        register const uint64_t size = tmp/2;
//        if (!is_del) correct = false;
//        else p_read+=size;
//      } else {
//        correct = false;
//      }
//      ++misms_pos;
//    } else {
//      correct = key[p_read] == text[p_ref];
//      p_read++; p_ref++;
//    }
//  }
//  return correct && p_read==key_length && misms_pos==num_misms && p_ref==alg_length;
//}











