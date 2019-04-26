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
 * Setup
 */
void alignment_init(
    alignment_t* const alignment,
    const uint64_t key_length,
    const uint64_t text_begin_offset,
    const uint64_t text_end_offset,
    const uint64_t max_error,
    const uint64_t num_tiles,
    const uint64_t tile_length) {
  // Allocate alignment
  alignment->num_tiles = num_tiles;
  alignment->distance_min_bound = ALIGN_DISTANCE_UNKNOWN;
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
    pattern_tiling_t pattern_tiling;
    pattern_tiling_init(&pattern_tiling,key_length,tile_length,text_length,adjusted_max_error);
    uint64_t tile_pos;
    for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
      // Init Tile
      alignment_tiles[tile_pos].distance = ALIGN_DISTANCE_UNKNOWN;
      alignment_tiles[tile_pos].text_end_offset = text_begin_offset+pattern_tiling.tile_offset+pattern_tiling.tile_wide;
      alignment_tiles[tile_pos].text_begin_offset = text_begin_offset+pattern_tiling.tile_offset;
      // Calculate next tile
      pattern_tiling_next(&pattern_tiling);
    }
  }
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
 * Verify edit distance using Kmer
 */
void alignment_verify_edit_kmer(
    alignment_t* const alignment,
    pattern_tiled_t* const pattern_tiled,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length,
    const uint64_t max_error,
    const uint64_t kmer_tiles,
    const uint64_t kmer_length) {
  PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_CANDIDATES);
  // Compute minimum edit distance bound
  const uint64_t distance_min_bound =
      kmer_counting_min_bound_nway(&pattern_tiled->kmer_filter_nway,
          text,text_length,max_error,pattern_tiled->mm_allocator);
  //  // DEBUG
  //  const bpm_pattern_t* const bpm_pattern = pattern_tiled->bpm_pattern;
  //  uint64_t distance, match_column;
  //  bpm_compute_edit_distance(bpm_pattern,text,text_length,
  //      &distance,&match_column,1000,false);
  //  if (distance < distance_min_bound) { fprintf(stderr,"error\n"); }
  if (distance_min_bound > max_error) {
    alignment->distance_min_bound = ALIGN_DISTANCE_INF;
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_DISCARDED);
  } else {
    alignment->distance_min_bound = ALIGN_DISTANCE_UNKNOWN;
    PROF_INC_COUNTER(GP_FC_KMER_COUNTER_FILTER_ACCEPTED);
  }
}
/*
 * Verify edit distance using BPM
 */
void alignment_verify_edit_bpm_tile(
    alignment_tile_t* const alignment_tile,
    pattern_tile_t* const pattern_tile,
    const uint8_t* const text,
    uint64_t* const tile_distance,
    uint64_t* const tile_match_column) {
  // Parameters
  bpm_pattern_t* const bpm_pattern_tile = &pattern_tile->bpm_pattern_tile;
  // Align BPM
  const uint64_t max_tile_error = pattern_tile->max_error;
  const uint64_t tile_offset = alignment_tile->text_begin_offset;
  const uint64_t tile_wide = alignment_tile->text_end_offset-alignment_tile->text_begin_offset;
  bpm_compute_edit_distance(bpm_pattern_tile,text+tile_offset,
      tile_wide,tile_distance,tile_match_column,max_tile_error,false); // TODO Activate quick-abandon at compile tiles
  PROF_INC_COUNTER(GP_BMP_DISTANCE_NUM_TILES_VERIFIED);
}
void alignment_verify_edit_bpm(
    alignment_t* const alignment,
    pattern_tiled_t* const pattern_tiled,
    uint8_t* const key,
    uint8_t* const text,
    const uint64_t max_error) {
  // Parameters
  const uint64_t num_pattern_tiles = pattern_tiled->num_tiles;
  // Align tiles
  uint64_t tile_pos, distance_min_bound=0;
  PROF_ADD_COUNTER(GP_CANDIDATE_TILES,num_pattern_tiles);
  for (tile_pos=0;tile_pos<num_pattern_tiles;++tile_pos) {
    alignment_tile_t* const alignment_tile = alignment->alignment_tiles + tile_pos;
    pattern_tile_t* const pattern_tile = pattern_tiled->tiles + tile_pos;
    if (alignment_tile->distance!=ALIGN_DISTANCE_UNKNOWN) {
      if (alignment_tile->distance==ALIGN_DISTANCE_INF) {
        distance_min_bound = ALIGN_DISTANCE_INF;
      } else {
        distance_min_bound += alignment_tile->distance;
      }
    } else {
      uint64_t tile_distance, tile_match_column;
      // Check BPM
      alignment_verify_edit_bpm_tile(
          alignment_tile,pattern_tile,
          text,&tile_distance,&tile_match_column);
      alignment_tile->distance = tile_distance;
      // Store tile alignment
      if (tile_distance==ALIGN_DISTANCE_INF) {
        distance_min_bound = ALIGN_DISTANCE_INF;
        break; // Stop verify
      } else {
        const uint64_t tile_end_offset = tile_match_column+1;
        const uint64_t tile_tall = pattern_tile->tile_length;
        const uint64_t tile_begin_offset =
        	BOUNDED_SUBTRACTION(tile_end_offset,tile_tall+tile_distance,0);
        const uint64_t tile_offset = alignment_tile->text_begin_offset;
        alignment_tile->distance = tile_distance;
        alignment_tile->text_end_offset = tile_offset + tile_end_offset;
        alignment_tile->text_begin_offset = tile_offset + tile_begin_offset;
        // Update distance
        distance_min_bound += tile_distance;
      }
    }
    // Check global-distance
    if (distance_min_bound > max_error) {
      distance_min_bound = ALIGN_DISTANCE_INF;
      break; // Stop verify
    }
  }
  // Setup alignment result
  alignment->distance_min_bound = distance_min_bound;
  // PROFILE
  PROF_BLOCK() {
    PROF_INC_COUNTER(GP_FC_BPM_FILTER_CANDIDATES);
    if (distance_min_bound==ALIGN_DISTANCE_INF) {
      PROF_INC_COUNTER(GP_FC_BPM_FILTER_DISCARDED);
    } else {
      PROF_INC_COUNTER(GP_FC_BPM_FILTER_ACCEPTED);
    }
  }
}
/*
 * Check alignment (CIGAR string against text & pattern)
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
    const bool local_alignment,
    const bool verbose) {
  // Traverse CIGAR
  uint64_t read_pos=0, text_pos=0;
  uint64_t i;
  for (i=0;i<cigar_length;++i) {
    cigar_element_t* const cigar_element = vector_get_elm(cigar_vector,cigar_offset+i,cigar_element_t);
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
        if (key[read_pos] == text[text_pos] && key[read_pos] != ENC_DNA_CHAR_N) {
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
  // Local Alignment
  if (local_alignment) return true;
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
 * Display
 */
void alignment_print_pretty(
    FILE* const stream,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint8_t* const text,
    const uint64_t text_length,
    vector_t* const cigar_vector,
    uint64_t const cigar_offset,
    uint64_t const cigar_length,
    mm_allocator_t* const mm_allocator) {
  mm_allocator_push_state(mm_allocator);
  tab_fprintf(stream,"[GEM]>Match.alignment\n");
  tab_global_inc();
  tab_fprintf(stream,"=> CIGAR  ");
  const uint64_t max_buffer_length = text_length+key_length+1;
  char* const key_alg = mm_allocator_calloc(mm_allocator,max_buffer_length,char,true);
  char* const ops_alg = mm_allocator_calloc(mm_allocator,max_buffer_length,char,true);
  char* const text_alg = mm_allocator_calloc(mm_allocator,max_buffer_length,char,true);
  cigar_element_t* cigar_element = vector_get_elm(cigar_vector,cigar_offset,cigar_element_t);
  uint64_t i, j, alg_pos = 0, read_pos = 0, text_pos = 0;
  for (i=0;i<cigar_length;++i,++cigar_element) {
    switch (cigar_element->type) {
      case cigar_match:
        fprintf(stream,"%d",(uint32_t)cigar_element->length);
        for (j=0;j<cigar_element->length;++j) {
          if (key[read_pos] != text[text_pos]) {
            key_alg[alg_pos] = dna_decode(key[read_pos]);
            ops_alg[alg_pos] = 'X';
            text_alg[alg_pos++] = dna_decode(text[text_pos]);
          } else {
            key_alg[alg_pos] = dna_decode(key[read_pos]);
            ops_alg[alg_pos] = '|';
            text_alg[alg_pos++] = dna_decode(text[text_pos]);
          }
          read_pos++; text_pos++;
        }
        break;
      case cigar_mismatch:
        fprintf(stream,"%c",dna_decode(cigar_element->mismatch));
        if (key[read_pos] != text[text_pos]) {
          key_alg[alg_pos] = dna_decode(key[read_pos++]);
          ops_alg[alg_pos] = 'M';
          text_alg[alg_pos++] = dna_decode(text[text_pos++]);
        } else {
          key_alg[alg_pos] = dna_decode(key[read_pos++]);
          ops_alg[alg_pos] = '!';
          text_alg[alg_pos++] = dna_decode(text[text_pos++]);
        }
        break;
      case cigar_ins:
        fprintf(stream,">%u+",cigar_element->length);
        for (j=0;j<cigar_element->length;++j) {
          key_alg[alg_pos] = '-';
          ops_alg[alg_pos] = ' ';
          text_alg[alg_pos++] = dna_decode(text[text_pos++]);
        }
        break;
      case cigar_del:
        for (j=0;j<cigar_element->length;++j) {
          key_alg[alg_pos] = dna_decode(key[read_pos++]);
          ops_alg[alg_pos] = ' ';
          text_alg[alg_pos++] = '-';
        }
        if (cigar_element->attributes==cigar_attr_trim) {
          fprintf(stream,"(%u)",cigar_element->length);
        } else {
          fprintf(stream,">%u-",cigar_element->length);
        }
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  key_alg[alg_pos] = '\0';
  ops_alg[alg_pos] = '\0';
  text_alg[alg_pos] = '\0';
  fprintf(stream,"\n");
  tab_fprintf(stream,"=> Pretty.Alignment\n");
  tab_fprintf(stream,"   KEY--%s--\n",key_alg);
  tab_fprintf(stream,"        %s  \n",ops_alg);
  tab_fprintf(stream,"   TXT--%s--\n",text_alg);
  tab_global_dec();
  mm_allocator_pop_state(mm_allocator);
}
