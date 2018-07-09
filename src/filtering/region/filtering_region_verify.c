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
 *   Filtering region module provides functions to verify a
 *   filtering-region against its corresponding text-region
 *   (compute the distance of the alignment between both)
 */

#include "align/alignment.h"
#include "filtering/region/filtering_region_verify.h"
#include "align/align_bpm_distance.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_REGION  GEM_DEEP_DEBUG
#define DEBUG_FILTERING_REGION_VERIFY_KMER_FILTER false

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Verify Hamming
 */
void filtering_region_verify_hamming_text(
    alignment_t* const alignment,
    const uint8_t* const text,
    const uint8_t* const key,
    const uint64_t key_length,
    const uint64_t max_mismatches) {
  // Check candidate
  uint64_t i, mismatches;
  for (i=0,mismatches=0;i<key_length;++i) {
    const uint8_t candidate_enc = text[i];
    // Check Mismatch
    if (candidate_enc==ENC_DNA_CHAR_N || candidate_enc != key[i]) {
      // Check Real Mismatch
      if (++mismatches > max_mismatches) {
        alignment->distance_min_bound = ALIGN_DISTANCE_INF;
        return;
      }
    }
  }
  alignment->distance_min_bound = mismatches;
}
void filtering_region_verify_hamming(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern) {
  // Parameters
  const uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  const uint64_t max_error = pattern->max_effective_filtering_error;
  text_trace_t* const text_trace = &filtering_region->text_trace;
  const uint64_t text_base_offset =
      filtering_region->text_source_region_offset -
      filtering_region->key_source_region_offset;
  const uint8_t* const text = text_trace->text + text_base_offset;
  const uint64_t text_length =
      filtering_region->text_end_position -
      (filtering_region->text_begin_position+text_base_offset);
  alignment_t* const alignment = &filtering_region->alignment;
  // Check length
  if (text_length >= key_length) {
    // Verify Hamming
    filtering_region_verify_hamming_text(alignment,text,key,key_length,max_error);
    if (alignment->distance_min_bound != ALIGN_DISTANCE_INF) {
      alignment->num_tiles = 1;
      alignment_tile_t* const alignment_tiles =
          filtering_candidates_allocate_alignment_tiles(filtering_candidates,1);
      alignment->alignment_tiles = alignment_tiles;
      alignment_tiles->distance = alignment->distance_min_bound;
      alignment_tiles->text_begin_offset = text_base_offset;
      alignment_tiles->text_end_offset = text_base_offset + key_length;
    }
  } else {
    alignment->distance_min_bound = ALIGN_DISTANCE_INF;
  }
}
/*
 * Verify Leveshtein
 */
void filtering_region_verify_levenshtein(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    const bool kmer_count_filter,
    pattern_t* const pattern) {
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  const uint64_t max_error = pattern->max_effective_filtering_error;
  text_trace_t* const text_trace = &filtering_region->text_trace;
  const uint64_t text_length = text_trace->text_padded_length;
  uint8_t* const text = text_trace->text_padded;
  alignment_t* const alignment = &filtering_region->alignment;
  // Prepare Alignment
  filtering_candidates_init_alignment(filtering_candidates,alignment,pattern,text_length,max_error);
  // Filter using kmer-counting
  if (kmer_count_filter) {
    alignment_verify_edit_kmer(
        &filtering_region->alignment,&pattern->pattern_tiled,
        pattern->key,pattern->key_length,text,text_length,max_error,
        search_parameters->candidate_verification.kmer_tiles,
        search_parameters->candidate_verification.kmer_length);
    if (alignment->distance_min_bound==ALIGN_DISTANCE_INF) {
      // DEBUG
      gem_cond_debug_block(DEBUG_FILTERING_REGION_VERIFY_KMER_FILTER)  {
        const bpm_pattern_t* const bpm_pattern = &pattern->pattern_tiled.bpm_pattern;
        uint64_t distance, match_column;
        bpm_compute_edit_distance(bpm_pattern,text,text_length,
            &distance,&match_column,max_error,false);
        gem_cond_error_msg(distance != ALIGN_DISTANCE_INF,
            "Filtering.Region.Verify: K-mer filtering wrong discarding (edit-distance=%"PRIu64")",distance);
      }
      // Return Discarded
      return;
    }
  }
  // Filter using BPM
  alignment_verify_edit_bpm(&filtering_region->alignment,
      &pattern->pattern_tiled,pattern->key,text,max_error);
}
/*
 * Verify (Switch)
 */
bool filtering_region_verify(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    const bool kmer_count_filter,
    pattern_t* const pattern) {
  PROF_START(GP_FC_VERIFY_CANDIDATES_REGION);
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  alignment_t* const alignment = &filtering_region->alignment;
  // Check align-distance (already known or verified)
  if (alignment->distance_min_bound == ALIGN_DISTANCE_UNKNOWN) {
    // Retrieve text-candidate
    filtering_region_retrieve_text(filtering_region,
        pattern,archive_text,filtering_candidates->mm_allocator);
    // Select alignment model
    switch (search_parameters->match_alignment_model) {
      case match_alignment_model_hamming: {
        // Verify Hamming
        filtering_region_verify_hamming(
            filtering_candidates,filtering_region,pattern);
        break;
      }
      case match_alignment_model_levenshtein:
      case match_alignment_model_gap_affine:
      case match_alignment_model_none:
        // Verify Levenshtein
        filtering_region_verify_levenshtein(
            filtering_candidates,filtering_region,
            kmer_count_filter,pattern);
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  // Check distance
  if (alignment->distance_min_bound != ALIGN_DISTANCE_INF) {
    filtering_region->status = filtering_region_accepted;
    PROF_STOP(GP_FC_VERIFY_CANDIDATES_REGION);
    PROF_INC_COUNTER(GP_ACCEPTED_REGIONS);
    return true;
  } else {
    filtering_region->status = filtering_region_verified_discarded;
    PROF_STOP(GP_FC_VERIFY_CANDIDATES_REGION);
    PROF_INC_COUNTER(GP_DISCARDED_REGIONS);
    return false;
  }
}
uint64_t filtering_region_verify_extension(
    filtering_candidates_t* const filtering_candidates,
    text_trace_t* const text_trace,
    const uint64_t index_position,
    pattern_t* const pattern,
    const uint64_t max_filtering_error) {
  PROFILE_START(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // Text (candidate)
  const uint8_t* const text = text_trace->text;
  const uint64_t text_length = text_trace->text_length;
  // Myers's BPM algorithm
  const uint64_t num_matches_found =
      bpm_compute_edit_distance_all(
          pattern,filtering_candidates,index_position,
          text,text_length,max_filtering_error);
  PROF_ADD_COUNTER(GP_ACCEPTED_REGIONS,num_matches_found);
  PROF_ADD_COUNTER(GP_FC_EXTEND_VERIFY_CANDIDATES_LENGTH,text_length);
  // Return number of filtering regions added (accepted)
  PROF_ADD_COUNTER(GP_FC_EXTEND_VERIFY_CANDIDATES_FOUND,num_matches_found);
  PROFILE_STOP(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS,PROFILE_LEVEL);
  return num_matches_found;
}
