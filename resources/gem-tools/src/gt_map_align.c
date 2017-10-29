/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_align.c
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_map_align.h"
#include "gt_output_map.h"

// Factor to multiply the read length as to allow expansion at realignment
#define GT_MAP_REALIGN_EXPANSION_FACTOR (0.20)

// Types of misms (as to backtrack the mismatches/indels)
#define GT_MAP_ALG_MISMS_NONE 0
#define GT_MAP_ALG_MISMS_MISMS 1
#define GT_MAP_ALG_MISMS_INS 2
#define GT_MAP_ALG_MISMS_DEL 3

/*
 * Map check/recover operators
 */
GT_INLINE gt_status gt_map_block_check_alignment(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern); GT_NULL_CHECK(sequence);
  // Check Alignment
  uint64_t pattern_centinel=0, sequence_centinel=0;
  // Init misms
  const uint64_t num_misms = gt_map_get_num_misms(map);
  uint64_t misms_offset=0;
  gt_misms* misms;
  GT_MAP_CHECK__RELOAD_MISMS_PTR(map,misms_offset,misms,num_misms);
  // Traverse the sequence
  while (pattern_centinel<pattern_length || sequence_centinel<sequence_length) {
    if (misms!=NULL && misms->position==pattern_centinel) { // Misms
      switch (misms->misms_type) {
        case MISMS:
          if (pattern_centinel>=pattern_length || sequence_centinel>=sequence_length) return GT_MAP_CHECK_ALG_MISMS_OUT_OF_SEQ;
          if (pattern[pattern_centinel]==sequence[sequence_centinel]) return GT_MAP_CHECK_ALG_NO_MISMS;
          if (misms->base!=sequence[sequence_centinel]) return GT_MAP_CHECK_ALG_BAD_MISMS;
          ++pattern_centinel; ++sequence_centinel;
          break;
        case INS:
          if (sequence_centinel+misms->size>sequence_length) return GT_MAP_CHECK_ALG_INS_OUT_OF_SEQ;
          sequence_centinel+=misms->size;
          break;
        case DEL:
          if (pattern_centinel+misms->size>pattern_length) return GT_MAP_CHECK_ALG_DEL_OUT_OF_SEQ;
          pattern_centinel+=misms->size;
          break;
      }
      ++misms_offset;
      GT_MAP_CHECK__RELOAD_MISMS_PTR(map,misms_offset,misms,num_misms);
    } else { // Match
      if (pattern_centinel>=pattern_length || sequence_centinel>=sequence_length) return GT_MAP_CHECK_ALG_MATCH_OUT_OF_SEQ;
      if (pattern[pattern_centinel]!=sequence[sequence_centinel]) return GT_MAP_CHECK_ALG_MISMATCH;
      ++pattern_centinel;
      ++sequence_centinel;
    }
  }
  // Ok, no complains
  return 0;
}
GT_INLINE gt_status gt_map_block_check_alignment_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Retrieve the sequence
  const uint64_t sequence_length = gt_map_get_length(map);
  gt_string* const sequence = gt_string_new(sequence_length+1);
  gt_status error_code;
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      sequence_length,0,sequence))) return error_code;
  // Check Alignment
  error_code = gt_map_block_check_alignment(map,
      gt_string_get_string(pattern),gt_string_get_length(pattern),
      gt_string_get_string(sequence),sequence_length);
  gt_string_delete(sequence);
  return error_code;
}
GT_INLINE gt_status gt_map_check_alignment_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Handle SMs
  gt_status error_code;
  if (gt_map_get_num_blocks(map)==1) { // Single-Block
    return gt_map_block_check_alignment_sa(map,pattern,sequence_archive);
  } else { // Realigning Slit-Maps (let's try not to spoil the splice-site consensus)
    gt_string* read_chunk = gt_string_new(0);
    uint64_t offset = 0;
    GT_MAP_ITERATE(map,map_block) {
      gt_string_set_nstring(read_chunk,gt_string_get_string(pattern)+offset,gt_map_get_base_length(map_block));
      if ((error_code=gt_map_block_check_alignment_sa(map_block,read_chunk,sequence_archive))) {
        gt_string_delete(read_chunk);
        return error_code;
      }
      offset += gt_map_get_base_length(map_block);
    }
    gt_string_delete(read_chunk);
    return 0;
  }
}
GT_INLINE gt_status gt_map_block_recover_mismatches(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern); GT_NULL_CHECK(sequence);
  GT_ZERO_CHECK(pattern_length); GT_ZERO_CHECK(sequence_length);
  // Set misms pattern
  gt_misms new_misms;
  new_misms.misms_type = MISMS;
  const uint64_t num_base_misms = gt_map_get_num_misms(map);
  // Traverse pattern & annotate mismatches
  uint64_t sequence_centinel=0, pattern_centinel=0;
  uint64_t misms_offset=0;
  gt_misms misms = {.position = UINT64_MAX};
  GT_MAP_CHECK__RELOAD_MISMS(map,misms_offset,misms,num_base_misms);
  while (sequence_centinel<sequence_length || pattern_centinel<pattern_length) {
    if (misms.position==pattern_centinel) { // Misms annotated
      switch (misms.misms_type) {
        case MISMS:
          if (pattern_centinel>=pattern_length || sequence_centinel>=sequence_length) {
            return GT_MAP_CHECK_ALG_MATCH_OUT_OF_SEQ;
          }
          if (pattern[pattern_centinel]==sequence[sequence_centinel]) return GT_MAP_CHECK_ALG_NO_MISMS;
          ++pattern_centinel; ++sequence_centinel;
          break;
        case INS:
          if (sequence_centinel+misms.size>sequence_length) return GT_MAP_CHECK_ALG_INS_OUT_OF_SEQ;
          sequence_centinel+=misms.size;
          break;
        case DEL:
          if (pattern_centinel+misms.size>pattern_length) return GT_MAP_CHECK_ALG_DEL_OUT_OF_SEQ;
          pattern_centinel+=misms.size;
          break;
      }
      // Add the misms + reload old misms
      gt_map_add_misms(map,&misms);
      ++misms_offset;
      GT_MAP_CHECK__RELOAD_MISMS(map,misms_offset,misms,num_base_misms);
    } else { // Nothing annotated
      if (pattern_centinel>=pattern_length || sequence_centinel>=sequence_length) {
        return GT_MAP_CHECK_ALG_MATCH_OUT_OF_SEQ;
      }
      if (pattern[pattern_centinel]!=sequence[sequence_centinel]) {
        new_misms.position = pattern_centinel;
        new_misms.base = sequence[sequence_centinel];
        gt_map_add_misms(map,&new_misms);
      }
      ++pattern_centinel; ++sequence_centinel;
    }
  }
  // Set new misms vector
  const uint64_t total_misms = gt_map_get_num_misms(map);
  uint64_t i, j;
  for (i=0,j=num_base_misms;j<total_misms;++j,++i) {
    gt_map_set_misms(map,gt_map_get_misms(map,j),i);
  }
  gt_map_set_num_misms(map,total_misms-num_base_misms);
  return 0;
}
GT_INLINE gt_status gt_map_block_recover_mismatches_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Retrieve the sequence
  gt_status error_code;
  const uint64_t pattern_length = gt_string_get_length(pattern);
  const uint64_t sequence_length = gt_map_get_length(map);
  gt_string* const sequence = gt_string_new(pattern_length+1);
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      sequence_length,0,sequence))) return error_code;
  // Recover mismatches
  const uint64_t num_misms = gt_map_get_num_misms(map);
  if ((error_code=gt_map_block_recover_mismatches(map,gt_string_get_string(pattern),
      pattern_length,gt_string_get_string(sequence),sequence_length))) {
    gt_map_set_num_misms(map,num_misms); // Restore state
    gt_error(MAP_RECOVER_MISMS_WRONG_BASE_ALG);
    gt_string_delete(sequence);
    return error_code;
  }
  gt_string_delete(sequence);
  return 0;
}
GT_INLINE gt_status gt_map_recover_mismatches_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Handle SMs
  gt_status error_code;
  if (gt_map_get_num_blocks(map)==1) { // Single-Block
    if ((error_code=gt_map_block_recover_mismatches_sa(map,pattern,sequence_archive))) {
      return error_code;
    }
    return 0;
  } else { // Split-Map (let's try not to spoil the splice-site consensus)
    uint64_t offset = 0;
    gt_string* read_chunk = gt_string_new(0);
    GT_MAP_ITERATE(map,map_block) {
      gt_string_set_nstring(read_chunk,gt_string_get_string(pattern)+offset,gt_map_get_base_length(map_block));
      if ((error_code=gt_map_block_recover_mismatches_sa(map_block,read_chunk,sequence_archive))) {
        gt_string_delete(read_chunk);
        return error_code;
      }
      offset += gt_map_get_base_length(map_block);
    }
    gt_string_delete(read_chunk);
    return 0;
  }
}
/*
 * Map (Re)alignment operators: HAMMING
 */
GT_INLINE gt_status gt_map_block_realign_hamming(
    gt_map* const map,char* const pattern,char* const sequence,const uint64_t length) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern);
  GT_NULL_CHECK(sequence);
  GT_ZERO_CHECK(length);
  // Set misms pattern & clear map mismatches
  gt_misms misms;
  misms.misms_type = MISMS;
  gt_map_clear_misms(map);
  // Traverse pattern & annotate mismatches
  uint64_t i;
  for (i=0;i<length;++i) {
    if (pattern[i]!=sequence[i]) {
      misms.position = i;
      misms.base = sequence[i];
      gt_map_add_misms(map,&misms);
    }
  }
  return 0;
}
GT_INLINE gt_status gt_map_block_realign_hamming_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Retrieve the sequence
  gt_status error_code;
  const uint64_t pattern_length = gt_string_get_length(pattern);
  gt_string* const sequence = gt_string_new(pattern_length+1);
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      pattern_length,0,sequence))) {
    gt_string_delete(sequence);
    return error_code;
  }
  // Realign Hamming
  error_code=gt_map_block_realign_hamming(map,gt_string_get_string(pattern),gt_string_get_string(sequence),pattern_length);
  gt_string_delete(sequence);
  return error_code;
}
GT_INLINE gt_status gt_map_realign_hamming_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Handle SMs
  gt_status error_code;
  if (gt_map_get_num_blocks(map)==1) { // Single-Block
    return gt_map_block_realign_hamming_sa(map,pattern,sequence_archive);
  } else { // Realigning Slit-Maps (let's try not to spoil the splice-site consensus)
    gt_string* read_chunk = gt_string_new(0);
    uint64_t offset = 0;
    GT_MAP_ITERATE(map,map_block) {
      gt_string_set_nstring(read_chunk,gt_string_get_string(pattern)+offset,gt_map_get_base_length(map_block));
      if ((error_code=gt_map_block_realign_hamming_sa(map_block,read_chunk,sequence_archive))) {
        gt_string_delete(read_chunk);
        return error_code;
      }
      offset += gt_map_get_base_length(map_block);
    }
    gt_string_delete(read_chunk);
    return 0;
  }
}
/*
 * Map (Re)alignment operators: LEVENSHTEIN
 */
#define GT_DP(i,j) dp_array[(i)*pattern_len+(j)]
#define GT_DP_SET_MISMS(misms,position_pattern,position_sequence,prev_misms,num_misms) { \
  misms.misms_type = MISMS; \
  misms.position = position_pattern; \
  misms.base = sequence[position_sequence]; \
  gt_map_add_misms(map,&misms); ++num_misms; \
  prev_misms = GT_MAP_ALG_MISMS_MISMS; \
}
#define GT_DP_SET_INS(map,misms,pos,length,prev_misms,num_misms) { \
  if (prev_misms != GT_MAP_ALG_MISMS_INS) { \
    misms.position = pos-length+2; \
    misms.misms_type = INS; \
    misms.size = length; \
    gt_map_add_misms(map,&misms); ++num_misms; \
  } else { \
    gt_misms* _misms = gt_map_get_misms(map,num_misms-1); \
    _misms->size+=length; \
  } \
  prev_misms = GT_MAP_ALG_MISMS_INS; \
}
#define GT_DP_SET_DEL(map,misms,pos,length,prev_misms,num_misms) { \
    if (prev_misms != GT_MAP_ALG_MISMS_DEL) { \
    misms.position = pos-length+1; \
    misms.misms_type = DEL; \
    misms.size = length; \
    gt_map_add_misms(map,&misms); ++num_misms; \
  } else { \
    gt_misms* _misms = gt_map_get_misms(map,num_misms-1); \
    _misms->position-=length; \
    _misms->size+=length; \
  } \
  prev_misms = GT_MAP_ALG_MISMS_DEL; \
}
GT_INLINE void gt_map_realign_dp_matrix_print(
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
GT_INLINE void gt_map_realign_levenshtein_get_distance(
    char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length,
    const bool ends_free,uint64_t* const position,uint64_t* const distance,
    gt_vector* const buffer) {
  GT_NULL_CHECK(pattern); GT_ZERO_CHECK(pattern_length);
  GT_NULL_CHECK(sequence); GT_ZERO_CHECK(sequence_length);
  // Allocate DP-matrix
  const uint64_t pattern_len = pattern_length+1;
  const uint64_t sequence_len = sequence_length+1;
  uint64_t* dp_array[2];
  if (buffer!=NULL) {
    gt_vector_prepare(buffer,uint64_t,2*pattern_len);
    dp_array[0] = gt_vector_get_mem(buffer,uint64_t);
  } else {
    dp_array[0] = gt_calloc(2*pattern_len,uint64_t,false);
  }
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
      dp_array[idx_b][j] = GT_MIN(sub,GT_MIN(ins,del));
    }
    // Check last cell value
    if (ends_free && dp_array[idx_b][pattern_length] < min_val) {
      min_val = dp_array[idx_b][pattern_length];
      i_pos = i;
    }
  }
  // Return results
  if (ends_free) {
    *position=i_pos-1;
    *distance=min_val;
  } else {
    *position=pattern_length;
    *distance=dp_array[idx_b][pattern_length];
  }
  // Free
  if (buffer==NULL) gt_free(dp_array[0]);
}
GT_INLINE gt_status gt_map_block_realign_levenshtein(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length,
    const bool ends_free,gt_vector* const buffer) {
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern); GT_ZERO_CHECK(pattern_length);
  GT_NULL_CHECK(sequence); GT_ZERO_CHECK(sequence_length);
  // Clear map misms
  gt_map_clear_misms(map);
  // Allocate DP matrix
  const uint64_t pattern_len = pattern_length+1;
  const uint64_t sequence_len = sequence_length+1;
  uint64_t* dp_array;
  if (buffer!=NULL) {
    gt_vector_prepare(buffer,uint64_t,pattern_len*sequence_len);
    dp_array = gt_vector_get_mem(buffer,uint64_t);
  } else {
    dp_array = gt_calloc(pattern_len*sequence_len,uint64_t,false);
  }
  // Init DP-Matrix
  uint64_t min_val = UINT64_MAX, i_pos = UINT64_MAX;
  uint64_t i, j;
  for (i=0;i<sequence_len;++i) GT_DP(i,0)=(ends_free)?0:i;
  for (j=0;j<pattern_len;++j) GT_DP(0,j)=j;
  // Calculate DP-Matrix
  for (i=1;i<sequence_len;++i) {
    for (j=1;j<pattern_len;++j) {
      const uint64_t ins = GT_DP(i-1,j) + 1;
      const uint64_t del = GT_DP(i,j-1) + 1;
      const uint64_t sub = GT_DP(i-1,j-1) + ((sequence[i-1]==pattern[j-1]) ? 0 : 1);
      GT_DP(i,j) = GT_MIN(sub,GT_MIN(ins,del));
    }
    // Check last cell value
    if (ends_free && GT_DP(i,pattern_length) < min_val) {
      min_val = GT_DP(i,pattern_length);
      i_pos = i;
    }
  }
  // DEBUG
  gt_map_realign_dp_matrix_print(dp_array,pattern_len,sequence_len,30,30);
  // Backtrack all edit operations
  uint64_t num_misms = 0, prev_misms = GT_MAP_ALG_MISMS_NONE;
  gt_misms misms;
  for (i=i_pos,j=pattern_len-1;i>0 && j>0;) {
    const uint32_t current_cell = GT_DP(i,j);
    if (sequence[i-1]==pattern[j-1]) { // Match
      prev_misms = GT_MAP_ALG_MISMS_NONE;
      --i; --j;
    } else {
      if (GT_DP(i-1,j)+1 == current_cell) { // Ins
        GT_DP_SET_INS(map,misms,j-1,1,prev_misms,num_misms);
        --i;
      } else if (GT_DP(i,j-1)+1 == current_cell) { // Del
        GT_DP_SET_DEL(map,misms,j-1,1,prev_misms,num_misms);
        --j;
      } else if (GT_DP(i-1,j-1)+1 == current_cell) { // Misms
        GT_DP_SET_MISMS(misms,j-1,i-1,prev_misms,num_misms);
        --i; --j;
      }
    }
  }
  if (i>0) {
    if (!ends_free) { // Insert the rest of the pattern
      GT_DP_SET_INS(map,misms,i-2,i,prev_misms,num_misms);
    } else {
      map->position+=i;
    }
  }
  if (j>0) { // Delete the rest of the sequence
    GT_DP_SET_DEL(map,misms,j-1,j,prev_misms,num_misms);
  }
  // Flip all mismatches
  uint64_t z;
  const uint64_t mid_point = num_misms/2;
  for (z=0;z<mid_point;++z) {
    misms = *gt_map_get_misms(map,z);
    gt_map_set_misms(map,gt_map_get_misms(map,num_misms-1-z),z);
    gt_map_set_misms(map,&misms,num_misms-1-z);
  }
  // Set map base length
  gt_map_set_base_length(map,pattern_length);
  // Safe check
  gt_debug_block(gt_map_block_check_alignment(map,pattern,pattern_length,
      sequence+((ends_free)?i:0),gt_map_get_length(map))!=0) {
    gt_output_map_fprint_map_block_pretty(stderr,map,pattern,pattern_length,
       sequence+((ends_free)?i:0),gt_map_get_length(map));
    gt_cond_fatal_error(gt_map_block_check_alignment(map,pattern,pattern_length,
      sequence+((ends_free)?i:0),gt_map_get_length(map))!=0,MAP_ALG_WRONG_ALG);
  }
  // Free
  if (buffer==NULL) gt_free(dp_array);
  return 0;
}
GT_INLINE gt_status gt_map_block_realign_levenshtein_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive,
    const float extra_length_factor,const bool ends_free) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  gt_status error_code;
  // Retrieve the sequence
  const uint64_t decode_length = (ends_free) ? gt_string_get_length(pattern) : gt_map_get_length(map);
  const uint64_t extra_length = (ends_free) ? ((float)decode_length*extra_length_factor) : 0;
  gt_string* const sequence = gt_string_new(decode_length+extra_length+1);
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      decode_length,extra_length,sequence))) {
    gt_string_delete(sequence); // Free
    return error_code;
  }
  // Realign Levenshtein
  error_code = gt_map_block_realign_levenshtein(map,
      gt_string_get_string(pattern),gt_string_get_length(pattern),
      gt_string_get_string(sequence),gt_string_get_length(sequence),ends_free,NULL);
  gt_string_delete(sequence); // Free
  return error_code;
}
GT_INLINE gt_status gt_map_realign_levenshtein_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Handle SMs
  gt_status error_code;
  if (gt_map_get_num_blocks(map)==1) {
    return gt_map_block_realign_levenshtein_sa(map,pattern,sequence_archive,GT_MAP_REALIGN_EXPANSION_FACTOR,true);
  } else { // Realigning MultipleBlocks (don't spoil the boundaries)
    gt_string* read_chunk = gt_string_new(0);
    uint64_t offset = 0;
    GT_MAP_ITERATE(map,map_block) {
      gt_string_set_nstring(read_chunk,gt_string_get_string(pattern)+offset,gt_map_get_base_length(map_block));
      if ((error_code=gt_map_block_realign_levenshtein_sa(map_block,read_chunk,sequence_archive,0.0,false))) {
        return error_code;
      }
      offset += gt_map_get_base_length(map_block);
    }
    gt_string_delete(read_chunk);
    return 0;
  }
}
GT_INLINE gt_status gt_map_block_realign_weighted(
    gt_map* const map,char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length,int32_t (*gt_weigh_fx)(char*,char*)) {
  GT_NOT_IMPLEMENTED(); // TODO
  return 0;
}
GT_INLINE gt_status gt_map_block_realign_weighted_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive,
    const float extra_length_factor,int32_t (*gt_weigh_fx)(char*,char*)) {
  GT_NOT_IMPLEMENTED(); // TODO
  return 0;
}
GT_INLINE gt_status gt_map_realign_weighted_sa(
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive,
    int32_t (*gt_weigh_fx)(char*,char*)) {
  GT_NOT_IMPLEMENTED(); // TODO
  return 0;
}
/*
 * Bit-compressed (Re)alignment
 *   BMP[BitParalellMayers] - Myers' Fast Bit-Vector algorithm (Levenshtein)
 */

#define GT_BMP_INITIAL_BUFFER_SIZE ((16*GT_DNA_RANGE+16*5)*64)
#define GT_BMP_W8_LENGTH 8
#define GT_BMP_W8_MASK (1<<7)
#define GT_PEQ_IDX(encoded_character,block_id,num_blocks) (encoded_character*num_blocks+block_id)

GT_INLINE gt_bpm_pattern* gt_map_bpm_pattern_new() {
  gt_bpm_pattern* const bpm_pattern = gt_alloc(gt_bpm_pattern);
  // Peq
  bpm_pattern->peq = NULL;
  // Auxiliary data
  bpm_pattern->P = NULL;
  bpm_pattern->M = NULL;
  bpm_pattern->level_mask = NULL;
  bpm_pattern->score = NULL;
  bpm_pattern->init_score = NULL;
  // Auxiliary buffer // TODO: Replace with slabs
  bpm_pattern->internal_buffer=gt_vector_new(GT_BMP_INITIAL_BUFFER_SIZE,1);
  // Return
  return bpm_pattern;
}
GT_INLINE void gt_map_bpm_pattern_delete(gt_bpm_pattern* const bpm_pattern) {
  gt_vector_delete(bpm_pattern->internal_buffer);
  gt_free(bpm_pattern);
}
GT_INLINE void gt_map_bpm_pattern_compile(
    gt_bpm_pattern* const bpm_pattern,const uint64_t word_size,
    char* const pattern,const uint64_t pattern_length) {
  GT_NULL_CHECK(bpm_pattern);
  GT_ZERO_CHECK(word_size);
  GT_NULL_CHECK(pattern);
  GT_ZERO_CHECK(pattern_length);
  // Calculate dimensions
  const uint64_t word_length = word_size*GT_BMP_W8_LENGTH;
  const uint64_t pattern_num_words = (pattern_length+(word_length-1))/word_length;
  const uint64_t peq_length = pattern_num_words*word_length;
  const uint64_t pattern_mod = pattern_length%word_length;
  // Init fields
  bpm_pattern->pattern_length = pattern_length;
  bpm_pattern->pattern_num_words = pattern_num_words;
  bpm_pattern->pattern_mod = pattern_mod;
  bpm_pattern->peq_length = peq_length;
  // Allocate memory
  const uint64_t peq_size = GT_DNA_RANGE*pattern_num_words*word_size;
  const uint64_t aux_vector_size = pattern_num_words*word_size;
  const uint64_t score_size = pattern_num_words*8;
  gt_vector_reserve(bpm_pattern->internal_buffer,peq_size+3*aux_vector_size+2*score_size,false);
  bpm_pattern->peq = gt_vector_get_mem(bpm_pattern->internal_buffer,void);
  bpm_pattern->P = bpm_pattern->peq + peq_size;
  bpm_pattern->M = bpm_pattern->P + aux_vector_size;
  bpm_pattern->level_mask = bpm_pattern->M + aux_vector_size;
  bpm_pattern->score = (int64_t*) (bpm_pattern->level_mask + aux_vector_size);
  bpm_pattern->init_score = bpm_pattern->score + pattern_num_words;
  // Init peq
  memset(bpm_pattern->peq,0,peq_size);
  uint64_t i;
  for (i=0;i<pattern_length;++i) {
    const uint8_t enc_char = gt_cdna_encode(pattern[i]);
    const uint64_t block = i/GT_BMP_W8_LENGTH;
    const uint8_t mask = 1<<(i%GT_BMP_W8_LENGTH);
    bpm_pattern->peq[GT_PEQ_IDX(enc_char,block,pattern_num_words*word_size)] |= mask;
  }
  for (;i<peq_length;++i) {
    const uint64_t block = i/GT_BMP_W8_LENGTH;
    const uint8_t mask = 1<<(i%GT_BMP_W8_LENGTH);
    uint64_t j;
    for (j=0;j<GT_DNA_RANGE;++j) {
      bpm_pattern->peq[GT_PEQ_IDX(j,block,pattern_num_words*word_size)] |= mask;
    }
  }
  // Init auxiliary data
  const uint64_t top = pattern_num_words-1;
  memset(bpm_pattern->level_mask,0,aux_vector_size);
  for (i=0;i<top;++i) {
    bpm_pattern->level_mask[i*word_size+(word_size-1)] = GT_BMP_W8_MASK;
    bpm_pattern->init_score[i] = word_length;
  }
  if (pattern_mod>0) {
    const uint64_t mask_shift = pattern_mod-1;
    const uint64_t block = mask_shift/GT_BMP_W8_LENGTH;
    const uint8_t mask = 1<<(mask_shift%GT_BMP_W8_LENGTH);
    bpm_pattern->level_mask[top*word_size+block] = mask;
    bpm_pattern->init_score[top] = pattern_mod;
  } else {
    bpm_pattern->level_mask[top*word_size+(word_size-1)] = GT_BMP_W8_MASK;
    bpm_pattern->init_score[top] = word_length;
  }
}

/*
 * Generic X-vector
 */
#define gt_bmp_vector_t uint8_t
#include "gt_map_align_bpm_generic.c"
#define gt_bmp_vector_t uint16_t
#include "gt_map_align_bpm_generic.c"
#define gt_bmp_vector_t uint32_t
#include "gt_map_align_bpm_generic.c"
#define gt_bmp_vector_t uint64_t
#include "gt_map_align_bpm_generic.c"
#define gt_bmp_vector_t __int128
#include "gt_map_align_bpm_generic.c"
