/*
 * PROJECT: GEMMapper
 * FILE: align.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */


#include "align/align.h"
#include "data_structures/dna_text.h"
#include "matches/match_alignment.h"

/*
 * Check matches (CIGAR string against text & pattern)
 */
bool align_check(
    FILE* const stream,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length,
    vector_t* const cigar_vector,
    uint64_t const cigar_offset,
    uint64_t const cigar_length,
    const bool verbose) {
  // Traverse CIGAR
  cigar_element_t* const cigar_base = vector_get_elm(cigar_vector,cigar_offset,cigar_element_t);
  uint64_t read_pos=0, text_pos=0;
  uint64_t i;
  for (i=0;i<cigar_length;++i) {
    cigar_element_t* const cigar_element = cigar_base + i;
    switch (cigar_element->type) {
      case cigar_match: {
        // Check all matching characters
        uint64_t j;
        for (j=0;j<cigar_element->length;++j) {
          if (key[read_pos] != text[text_pos]) {
            if (verbose) {
              fprintf(stream,"Align Check. Alignment not matching "
                  "(key[%"PRIu64"]=%c != text[%"PRIu64"]=%c)\n",
                  read_pos,dna_decode(key[read_pos]),text_pos,
                  dna_decode(text[text_pos]));
            }
            return false;
          }
          ++read_pos;
          ++text_pos;
        }
        break;
      }
      case cigar_mismatch:
        // Check mismatch
        if (key[read_pos] == text[text_pos]) {
          if (verbose) {
            fprintf(stream,"Align Check. Alignment not mismatching "
                "(key[%"PRIu64"]=%c == text[%"PRIu64"]=%c, CIGAR=%c)\n",
              read_pos,dna_decode(key[read_pos]),text_pos,dna_decode(text[text_pos]),
              dna_decode(cigar_element->mismatch));
          }
          return false;
        } else if (cigar_element->mismatch != text[text_pos]) {
          if (verbose) {
            fprintf(stream,"Align Check. Alignment not mismatching as CIGAR states "
                "(key[%"PRIu64"]=%c == text[%"PRIu64"]=%c, CIGAR=%c)\n",
              read_pos,dna_decode(key[read_pos]),text_pos,dna_decode(text[text_pos]),
              dna_decode(cigar_element->mismatch));
          }
          return false;
        }
        ++read_pos;
        ++text_pos;
        break;
      case cigar_ins:
        text_pos += cigar_element->length;
        break;
      case cigar_del:
        read_pos += cigar_element->length;
        break;
      case cigar_null:
        gem_cond_error_msg(verbose,"Align Check. CIGAR Null");
        return false;
        break;
      default:
        break;
    }
  }
  // Check alignment length
  if (read_pos != key_length) {
    if (verbose) {
      fprintf(stream,"Align Check. Alignment incorrect length "
          "(key-aligned=%"PRIu64",key-length=%"PRIu64")\n",read_pos,key_length);
    }
    return false;
  }
  if (text_pos != text_length) {
    if (verbose) {
      fprintf(stream,"Align Check. Alignment incorrect length "
          "(text-aligned=%"PRIu64",text-length=%"PRIu64")\n",text_pos,text_length);
    }
    return false;
  }
  return true;
}
/*
 * Compute edit distance (Basic DP-Matrix Alignment)
 */
int64_t align_dp_compute_edit_distance(
    const char* const key,
    const uint64_t key_length,
    const char* const text,
    const uint64_t text_length,
    const bool ends_free,
    uint64_t* const position) {
  GEM_CHECK_NULL(key); GEM_CHECK_ZERO(key_length);
  GEM_CHECK_NULL(text); GEM_CHECK_ZERO(text_length);
  // Allocate DP-matrix
  const uint64_t key_len = key_length+1;
  const uint64_t text_len = text_length+1;
  uint64_t* dp_array[2];
  dp_array[0] = mm_calloc(2*key_len,uint64_t,false);
  dp_array[1] = dp_array[0] + key_len;
  // Init DP-Matrix
  uint64_t min_val = UINT64_MAX, i_pos = UINT64_MAX;
  uint64_t i, j, idx_a=0, idx_b=0;
  for (j=0;j<key_len;++j) dp_array[0][j]=j;
  // Calculate DP-Matrix
  for (i=1;i<text_len;++i) {
    // Fix indexes
    idx_a = idx_b;
    idx_b = i % 2;
    // Fix first cell
    dp_array[idx_b][0] = (ends_free) ? 0 : dp_array[idx_a][0]+1;
    // Develop row
    for (j=1;j<key_len;++j) {
      const uint64_t ins = dp_array[idx_a][j]   + 1;
      const uint64_t del = dp_array[idx_b][j-1] + 1;
      const uint64_t sub = dp_array[idx_a][j-1] + ((text[i-1]==key[j-1]) ? 0 : 1);
      dp_array[idx_b][j] = MIN(sub,MIN(ins,del));
    }
    // Check last cell value
    if (ends_free && dp_array[idx_b][key_length] < min_val) {
      min_val = dp_array[idx_b][key_length];
      i_pos = i;
    }
  }
  // Return results & Free
  int64_t distance = INT64_MAX;
  if (ends_free) {
    *position = i_pos-1;
    distance = min_val;
  } else {
    *position = key_length;
    distance = dp_array[idx_b][key_length];
  }
  mm_free(dp_array[0]);
  return distance;
}
