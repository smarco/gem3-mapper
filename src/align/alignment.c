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
 *   Basic alignment module to check edit-distance based alignments/CIGARs
 *   and to generate the best correct CIGAR by means of a simple algorithm
 *   as to double-check
 */

#include "text/dna_text.h"
#include "align/alignment.h"
#include "align/align_bpm_distance.h"
#include "matches/align/match_alignment.h"
#include "matches/matches_cigar.h"

/*
 * Check matches (CIGAR string against text & pattern)
 */
bool alignment_check(
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
int64_t alignment_dp_compute_edit_distance(
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
/*
 * Verify levenshtein using BPM
 */
void alignment_verify_levenshtein_bpm(
    alignment_t* const alignment,
    const uint64_t filtering_max_error,
    bpm_pattern_t* const bpm_pattern,
    bpm_pattern_t* const bpm_pattern_tiles,
    text_trace_t* const text_trace) {
  // Parameters
  const uint8_t* const text = text_trace->text;
  const uint64_t num_pattern_tiles = bpm_pattern_tiles->num_pattern_tiles;
  alignment_tile_t* const alignment_tiles = alignment->alignment_tiles;
  const uint64_t max_error = MIN(filtering_max_error,bpm_pattern->pattern_length);
  // Align tiles
  uint64_t max_remaining_error = max_error;
  uint64_t tile_pos, global_distance=0;
  PROF_ADD_COUNTER(GP_BMP_DISTANCE_NUM_TILES,num_pattern_tiles);
  for (tile_pos=0;tile_pos<num_pattern_tiles;++tile_pos) {
    alignment_tile_t* const alignment_tile = alignment_tiles + tile_pos;
    if (alignment_tile->distance!=ALIGN_DISTANCE_INF) {
      global_distance += alignment_tile->distance;
    } else {
      bpm_pattern_t* const bpm_pattern_tile = bpm_pattern_tiles+tile_pos;
      const uint64_t max_tile_error = MIN(max_remaining_error,bpm_pattern->pattern_length);
      const uint64_t tile_offset = alignment_tile->text_begin_offset;
      const uint64_t tile_wide = alignment_tile->text_end_offset-alignment_tile->text_begin_offset;
      uint64_t tile_distance, tile_match_column;
      bpm_compute_edit_distance(bpm_pattern_tile,text+tile_offset,
          tile_wide,&tile_distance,&tile_match_column,max_tile_error,false);
      PROF_INC_COUNTER(GP_BMP_DISTANCE_NUM_TILES_VERIFIED);
      // Store tile alignment
      if (tile_distance!=ALIGN_DISTANCE_INF) {
        const uint64_t tile_end_offset = tile_match_column+1;
        const uint64_t tile_tall = bpm_pattern_tile->pattern_length;
        const uint64_t tile_begin_offset = BOUNDED_SUBTRACTION(tile_end_offset,tile_tall+tile_distance,0);
        alignment_tile->distance = tile_distance;
        alignment_tile->text_end_offset = tile_offset + tile_end_offset;
        alignment_tile->text_begin_offset = tile_offset + tile_begin_offset;
        // Update distance
        max_remaining_error = BOUNDED_SUBTRACTION(max_remaining_error,tile_distance,0);
        global_distance += tile_distance;
      } else {
        global_distance = ALIGN_DISTANCE_INF;
        break; // Stop verify
      }
    }
    // Check global-distance
    if (global_distance > max_error) {
      global_distance = ALIGN_DISTANCE_INF;
      break; // Stop verify
    }
  }
  // Setup alignment result
  alignment->distance_min_bound = global_distance;
}
