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
    const bool* const allowed_enc,
    const uint64_t max_mismatches) {
  // Check candidate
  uint64_t i, mismatches;
  for (i=0,mismatches=0;i<key_length;++i) {
    const uint8_t candidate_enc = text[i];
    // Check Mismatch
    if (!allowed_enc[candidate_enc] || candidate_enc != key[i]) {
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
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const bool* const allowed_enc,
    mm_stack_t* const mm_stack) {
  // Parameters
  const uint8_t* const key = pattern->key;
  const uint64_t key_length = pattern->key_length;
  const uint64_t max_error = pattern->max_effective_filtering_error;
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
    filtering_region_verify_hamming_text(alignment,text,key,key_length,allowed_enc,max_error);
    if (alignment->distance_min_bound != ALIGN_DISTANCE_INF) {
      alignment->num_tiles = 1;
      alignment_tile_t* const alignment_tiles = mm_stack_calloc(mm_stack,1,alignment_tile_t,false);
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
    pattern_t* const pattern,
    text_trace_t* const text_trace,
    const bool kmer_filter,
    mm_stack_t* const mm_stack) {
  // Parameters
  const uint64_t max_error = pattern->max_effective_filtering_error;
  const uint64_t text_length = text_trace->text_padded_length;
  uint8_t* const text = text_trace->text_padded;
  // Prepare Alignment
  filtering_candidates_init_alignment(
      filtering_candidates,&filtering_region->alignment,
      pattern,text_length,max_error,false);
  // Verify Levenshtein
  alignment_verify_levenshtein(
      &filtering_region->alignment,
      &pattern->alignment_filters,
      pattern->key,text,max_error);
}
/*
 * Verify (Switch)
 */
bool filtering_region_verify(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern,
    const bool kmer_filter) {
  PROF_START(GP_FC_VERIFY_CANDIDATES_REGION);
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  text_collection_t* const text_collection = &filtering_candidates->text_collection;
  alignment_t* const alignment = &filtering_region->alignment;
  // Check align-distance (already known or verified)
  if (alignment->distance_min_bound == ALIGN_DISTANCE_UNKNOWN) {
    // Retrieve text-candidate
    filtering_region_retrieve_text(filtering_region,pattern,archive_text,text_collection);
    text_trace_t* const text_trace = text_collection_get_trace(
        text_collection,filtering_region->text_trace_offset);
    // Select alignment model
    switch (search_parameters->match_alignment_model) {
      case match_alignment_model_hamming: {
        // Verify Hamming
        filtering_region_verify_hamming(filtering_region,pattern,text_trace,
            search_parameters->allowed_enc,filtering_candidates->mm->mm_general);
        break;
      }
      case match_alignment_model_levenshtein:
      case match_alignment_model_gap_affine:
        // Verify Levenshtein
        filtering_region_verify_levenshtein(filtering_candidates,filtering_region,
            pattern,text_trace,kmer_filter,filtering_candidates->mm->mm_general);
        break;
      case match_alignment_model_none:
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
    const uint64_t text_trace_offset,
    const uint64_t index_position,
    pattern_t* const pattern) {
  PROFILE_START(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS,PROFILE_LEVEL);
  // Text (candidate)
  text_trace_t* const text_trace = text_collection_get_trace(
      &filtering_candidates->text_collection,text_trace_offset);
  const uint8_t* const text = text_trace->text;
  const uint64_t text_length = text_trace->text_length;
  // Pattern
  const uint64_t max_filtering_error = pattern->max_effective_filtering_error;
  // Myers's BPM algorithm
  const uint64_t num_matches_found =
      bpm_compute_edit_distance_all(
          pattern,filtering_candidates,
          text_trace_offset,index_position,
          text,text_length,max_filtering_error);
  PROF_ADD_COUNTER(GP_ACCEPTED_REGIONS,num_matches_found);
  PROF_ADD_COUNTER(GP_FC_EXTEND_VERIFY_CANDIDATES_LENGTH,text_length);
  // Return number of filtering regions added (accepted)
  PROF_ADD_COUNTER(GP_FC_EXTEND_VERIFY_CANDIDATES_FOUND,num_matches_found);
  PROFILE_STOP(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS,PROFILE_LEVEL);
  return num_matches_found;
}
