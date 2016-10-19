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
 * Debug
 */
#define FILTERING_REGION_VERIFY_CHECK_KMER_FILTER false

/*
 * Setup
 */
void alignment_init(
    alignment_t* const alignment,
    const uint64_t key_length,
    const uint64_t text_begin_offset,
    const uint64_t text_end_offset,
    const uint64_t max_error,
    const uint64_t num_tiles,
    const uint64_t tile_length,
    mm_stack_t* const mm_stack) {
  // Allocate alignment
  alignment->num_tiles = num_tiles;
  alignment->distance_min_bound = ALIGN_DISTANCE_UNKNOWN;
  if (alignment->alignment_tiles==NULL) {
    alignment->alignment_tiles = mm_stack_calloc(mm_stack,num_tiles,alignment_tile_t,false);
  }
  // Init all tiles
  const uint64_t text_length = text_end_offset-text_begin_offset;
  alignment_tile_t* const alignment_tiles = alignment->alignment_tiles;
  if (num_tiles==1) {
    alignment_tiles->distance = ALIGN_DISTANCE_UNKNOWN;
    alignment_tiles->text_end_offset = text_end_offset;
    alignment_tiles->text_begin_offset = text_begin_offset;
  } else {
    // Error parameters
    const uint64_t adjusted_max_error = MIN(max_error,tile_length);
    // Calculate tile dimensions
    pattern_tiled_t pattern_tiled;
    pattern_tiled_init(&pattern_tiled,key_length,tile_length,text_length,adjusted_max_error);
    uint64_t tile_pos;
    for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
      // Init Tile
      alignment_tiles[tile_pos].distance = ALIGN_DISTANCE_UNKNOWN;
      alignment_tiles[tile_pos].text_end_offset = text_begin_offset+pattern_tiled.tile_offset+pattern_tiled.tile_wide;
      alignment_tiles[tile_pos].text_begin_offset = text_begin_offset+pattern_tiled.tile_offset;
      // Calculate next tile
      pattern_tiled_calculate_next(&pattern_tiled);
    }
  }
}
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
 * Verify levenshtein using Filters (BPM + kmer-counting)
 */
uint64_t alignment_verify_levenshtein_kmer_filter(
    alignment_tile_t* const alignment_tile,
    alignment_filters_tile_t* const filters_tiles,
    uint8_t* const key,
    uint8_t* const text,
    mm_stack_t* const mm_stack) {
  // Parameters
  const uint64_t tile_offset = alignment_tile->text_begin_offset;
  const uint64_t tile_wide = alignment_tile->text_end_offset-alignment_tile->text_begin_offset;
  // Check kmer-filter compiled
  if (filters_tiles->kmer_filter_tile==NULL) {
    filters_tiles->kmer_filter_tile = mm_stack_alloc(mm_stack,kmer_counting_t);
    kmer_counting_compile(
        filters_tiles->kmer_filter_tile,
        key + filters_tiles->tile_offset,
        filters_tiles->tile_length,
        filters_tiles->max_error,
        mm_stack);
  }
  // Kmer Filter
  const uint64_t test_positive = kmer_counting_filter(
      filters_tiles->kmer_filter_tile,text+tile_offset,tile_wide);
  if (test_positive==ALIGN_DISTANCE_INF) {
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_DISCARDED);
    // DEBUG
    gem_cond_debug_block(FILTERING_REGION_VERIFY_CHECK_KMER_FILTER)  {
      const bpm_pattern_t* const bpm_pattern = filters_tiles->bpm_pattern_tile;
      uint64_t distance, match_column;
      bpm_compute_edit_distance(bpm_pattern,text+tile_offset,tile_wide,
          &distance,&match_column,filters_tiles->max_error,false);
      gem_cond_error_msg(distance != ALIGN_DISTANCE_INF,
          "Filtering.Region.Verify: K-mer filtering wrong discarding (edit-distance=%lu)",distance);
    }
    return ALIGN_DISTANCE_INF;
  } else {
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_ACCEPTED);
    return ALIGN_DISTANCE_UNKNOWN;
  }
}
void alignment_verify_levenshtein_bpm(
    alignment_tile_t* const alignment_tile,
    alignment_filters_tile_t* const filters_tiles,
    const uint8_t* const text,
    uint64_t* const tile_distance,
    uint64_t* const tile_match_column) {
  // Parameters
  bpm_pattern_t* const bpm_pattern_tile = filters_tiles->bpm_pattern_tile;
  // Align BPM
  const uint64_t max_tile_error = filters_tiles->max_error;
  const uint64_t tile_offset = alignment_tile->text_begin_offset;
  const uint64_t tile_wide = alignment_tile->text_end_offset-alignment_tile->text_begin_offset;
  bpm_compute_edit_distance(bpm_pattern_tile,text+tile_offset,
      tile_wide,tile_distance,tile_match_column,max_tile_error,false); // TODO Activate quick-abandon at compile tiles
  PROF_INC_COUNTER(GP_BMP_DISTANCE_NUM_TILES_VERIFIED);
}
void alignment_verify_levenshtein(
    alignment_t* const alignment,
    alignment_filters_t* const filters,
    uint8_t* const key,
    uint8_t* const text,
    const uint64_t max_error) {
  // Parameters
  alignment_tile_t* const alignment_tiles = alignment->alignment_tiles;
  // Parameters pattern
  alignment_filters_tile_t* const alignment_filters_tiles = filters->tiles;
  const uint64_t num_pattern_tiles = filters->num_tiles;
  // Align tiles
  uint64_t tile_pos, global_distance=0;
  PROF_ADD_COUNTER(GP_CANDIDATE_TILES,num_pattern_tiles);
  for (tile_pos=0;tile_pos<num_pattern_tiles;++tile_pos) {
    alignment_tile_t* const alignment_tile = alignment_tiles + tile_pos;
    alignment_filters_tile_t* const alignment_filters_tile = alignment_filters_tiles + tile_pos;
    if (alignment_tile->distance!=ALIGN_DISTANCE_UNKNOWN) {
      if (alignment_tile->distance==ALIGN_DISTANCE_INF) {
        global_distance = ALIGN_DISTANCE_INF;
      } else {
        global_distance += alignment_tile->distance;
      }
    } else {
      uint64_t tile_distance, tile_match_column;
      // Use kmer filter
      tile_distance = alignment_verify_levenshtein_kmer_filter(
          alignment_tile,alignment_filters_tile,key,text,filters->mm_stack);
      // Check BPM
      if (tile_distance==ALIGN_DISTANCE_UNKNOWN) {
        alignment_verify_levenshtein_bpm(
            alignment_tile,alignment_filters_tile,
            text,&tile_distance,&tile_match_column);
      }
      alignment_tile->distance = tile_distance;
      // Store tile alignment
      if (tile_distance==ALIGN_DISTANCE_INF) {
        global_distance = ALIGN_DISTANCE_INF;
        break; // Stop verify
      } else {
        const uint64_t tile_end_offset = tile_match_column+1;
        const uint64_t tile_tall = alignment_filters_tile->tile_length;
        const uint64_t tile_begin_offset = BOUNDED_SUBTRACTION(tile_end_offset,tile_tall+tile_distance,0);
        const uint64_t tile_offset = alignment_tile->text_begin_offset;
        alignment_tile->distance = tile_distance;
        alignment_tile->text_end_offset = tile_offset + tile_end_offset;
        alignment_tile->text_begin_offset = tile_offset + tile_begin_offset;
        // Update distance
        global_distance += tile_distance;
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
