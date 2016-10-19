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
 *   Filtering candidates module provides functions to process all
 *   bwt-encoded positions (originated from a candidate generation process)
 *   and compose them into filtering-regions (using decoded text coordinates)
 *   This "buffered" module operates in batches of filtering-candidates and
 *   makes use of GPU-buffers to offloads the decoding of positions to a GPU
 */

#include "align/alignment.h"
#include "filtering/candidates/filtering_candidates_process_buffered.h"
#include "filtering/candidates/filtering_candidates_process.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Decode Candidates Helpers (Buffered/Batch)
 */
void filtering_candidates_decode_retrieve_text_sample(
    filtering_candidates_t* const filtering_candidates,
    filtering_position_t* const filtering_position,
    const uint64_t filtering_position_offset,
    const uint64_t region_interval_lo,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_offset_begin) {
  // Parameters
  fm_index_t* const fm_index = filtering_candidates->archive->fm_index;
  // Retrieve decoded position & fetch sample
  uint64_t text_position;
  gpu_buffer_fmi_decode_get_position_text(gpu_buffer_fmi_decode,
      buffer_offset_begin+filtering_position_offset,&text_position);
  if (text_position == -1) {
    // Re-Decode (GPU decode failed)
    PROF_INC_COUNTER(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL_TOTAL);
    PROFILE_START(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL,PROFILE_LEVEL);
    filtering_position->region_text_position =
        fm_index_decode(fm_index,region_interval_lo+filtering_position_offset);
    PROFILE_STOP(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL,PROFILE_LEVEL);
  } else {
    filtering_position->region_text_position = text_position;
  }
  // DEBUG
#ifdef CUDA_CHECK_BUFFERED_DECODE_POSITIONS
  const uint64_t region_text_position = fm_index_decode(fm_index,region_interval_lo+filtering_position_offset);
  gem_cond_fatal_error_msg(filtering_position->region_text_position!=region_text_position,
      "Filtering.Candidates.Process.Buffered. Check decoded position failed (%lu!=%lu)",
      filtering_position->region_text_position,region_text_position);
#endif
}
void filtering_candidates_decode_batch_retrieve_sampled_position(
    filtering_candidates_t* const filtering_candidates,
    fc_batch_decode_candidate* const batch_candidate,
    const uint64_t region_interval_lo,
    const uint64_t region_interval_position,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_offset_begin) {
  // Parameters
  fm_index_t* const fm_index = filtering_candidates->archive->fm_index;
  const bwt_t* const bwt = fm_index->bwt;
  // Retrieve sampled-position from buffer
  gpu_buffer_fmi_decode_get_position_sa(
      gpu_buffer_fmi_decode,buffer_offset_begin+region_interval_position,
      &batch_candidate->index_position,&batch_candidate->distance);
  // If buffered-decode failed, retrieve sampled-position from index
  if (batch_candidate->index_position == -1) {
    PROF_INC_COUNTER(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL_TOTAL);
    PROFILE_START(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL,PROFILE_LEVEL);
    fm_index_retrieve_bwt_sampled(fm_index,region_interval_lo+region_interval_position,
        &batch_candidate->index_position,&batch_candidate->distance);
    PROFILE_STOP(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL,PROFILE_LEVEL);
  }
  // Prefetch BM-sampled
  bwt_prefetch(bwt,batch_candidate->index_position,&batch_candidate->bwt_block_locator);
}
void filtering_candidates_decode_batch_retrieve_bm_sampled(
    filtering_candidates_t* const filtering_candidates,
    fc_batch_decode_candidate* const batch_candidate) {
  // Parameters
  fm_index_t* const fm_index = filtering_candidates->archive->fm_index;
  const bwt_t* const bwt = fm_index->bwt;
  const sampled_sa_t* const sampled_sa = fm_index->sampled_sa;
  // Retrieve BM-sampled (position of the SA-sample)
  bool is_sampled;
  batch_candidate->index_position =
      bwt_prefetched_LF(bwt,batch_candidate->index_position,
          &is_sampled,&batch_candidate->bwt_block_locator);
  // Prefetch sampled-position
  sampled_sa_prefetch_sample(sampled_sa,batch_candidate->index_position);
}
void filtering_candidates_decode_batch_retrieve_sa_sample(
    filtering_candidates_t* const filtering_candidates,
    fc_batch_decode_candidate* const batch_candidate,
    filtering_position_t* const filtering_position,
    const uint64_t region_interval_lo,
    const uint64_t region_interval_position) {
  // Parameters
  fm_index_t* const fm_index = filtering_candidates->archive->fm_index;
  const sampled_sa_t* const sampled_sa = fm_index->sampled_sa;
  const uint64_t bwt_length = fm_index_get_length(fm_index);
  // Retrieve sampled-position
  const uint64_t sampled_position = sampled_sa_get_sample(sampled_sa,batch_candidate->index_position);
  filtering_position->region_text_position = (sampled_position + batch_candidate->distance) % bwt_length;
  // DEBUG
#ifdef CUDA_CHECK_BUFFERED_DECODE_POSITIONS
  const uint64_t region_text_position = fm_index_decode(fm_index,region_interval_lo+region_interval_position);
  gem_cond_fatal_error_msg(filtering_position->region_text_position!=region_text_position,
      "Filtering.Candidates.Process.Buffered. Check decoded position failed (%lu!=%lu)",
      filtering_position->region_text_position,region_text_position);
#endif
}
/*
 * Decode Candidates Buffered (from GPU-Buffer)
 */
void filtering_candidates_decode_filtering_positions_buffered(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    const uint64_t buffered_positions_offset,
    pattern_t* const pattern,
    region_search_t* const region_search) {
  PROFILE_START(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  // Parameters
  filtering_position_buffered_t* const positions_buffered =
      filtering_candidates_buffered->positions_buffered + buffered_positions_offset;
  const bool exact_match = (region_search->begin==0 && region_search->end==pattern->key_length);
  const uint64_t region_lo = region_search->lo;
  const uint64_t region_hi = region_search->hi;
  // Add all candidate positions
  const uint64_t num_candidates = region_hi-region_lo;
  uint64_t i;
  for (i=0;i<num_candidates;++i) {
    // Allocate
    filtering_position_t* const fposition = filtering_candidates_allocate_position(filtering_candidates);
    // Configure
    fposition->source_region_begin = positions_buffered[i].source_region_begin;
    fposition->source_region_end = positions_buffered[i].source_region_end;
    fposition->source_region_error = 0;
    // Decode position
    fposition->region_text_position =
        fm_index_decode(filtering_candidates->archive->fm_index,region_lo+i);
    // Adjust Position
    filtering_candidates_compute_text_coordinates(filtering_candidates,fposition,pattern);
    fposition->align_distance = exact_match ? 0 : ALIGN_DISTANCE_UNKNOWN;
  }
  PROF_ADD_COUNTER(GP_FC_DECODE_POSITIONS,num_candidates);
  PROF_ADD_COUNTER(GP_ASSW_DECODE_CANDIDATES_RETRIVED,num_candidates);
  PROFILE_STOP(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
}
void filtering_candidates_decode_sa_filtering_positions_buffered(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    const uint64_t buffered_positions_offset,
    pattern_t* const pattern,
    region_search_t* const region_search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_offset_begin) {
  PROFILE_START(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  // Parameters
  filtering_position_buffered_t* const positions_buffered =
      filtering_candidates_buffered->positions_buffered + buffered_positions_offset;
  const bool exact_match = (region_search->begin==0 && region_search->end==pattern->key_length);
  const uint64_t region_lo = region_search->lo;
  const uint64_t region_hi = region_search->hi;
  const uint64_t num_candidates = region_hi-region_lo;
  // Retrieve decoded positions
  fc_batch_decode_candidate batch[DECODE_NUM_POSITIONS_PREFETCHED];
  uint64_t num_left_positions = num_candidates;
  uint64_t i, current_position = 0;
  while (num_left_positions > 0) {
    const uint64_t batch_size = MIN(num_left_positions,DECODE_NUM_POSITIONS_PREFETCHED);
    // Retrieve sampled-position & prefetch BM-sampled
    for (i=0;i<batch_size;++i) {
      filtering_candidates_decode_batch_retrieve_sampled_position(
          filtering_candidates,batch+i,region_lo,current_position+i,
          gpu_buffer_fmi_decode,buffer_offset_begin);
    }
    // Retrieve BM-sampled & prefetch SA-sample
    for (i=0;i<batch_size;++i) {
      filtering_candidates_decode_batch_retrieve_bm_sampled(filtering_candidates,batch+i);
    }
    // Retrieve SA-sample, locate position & adjust
    for (i=0;i<batch_size;++i) {
      // Allocate position
      filtering_position_t* const fposition = filtering_candidates_allocate_position(filtering_candidates);
      // Retrieve SA-sample
      fposition->source_region_begin = positions_buffered[current_position+i].source_region_begin;
      fposition->source_region_end = positions_buffered[current_position+i].source_region_end;
      fposition->source_region_error = 0;
      filtering_candidates_decode_batch_retrieve_sa_sample(
          filtering_candidates,batch+i,fposition,region_lo,current_position+i);
      // Adjust Position
      filtering_candidates_compute_text_coordinates(filtering_candidates,fposition,pattern);
      fposition->align_distance = exact_match ? 0 : ALIGN_DISTANCE_UNKNOWN;
    }
    // Next batch
    current_position = current_position + batch_size;
    num_left_positions -= batch_size;
  }
  PROF_ADD_COUNTER(GP_FC_DECODE_POSITIONS,num_candidates);
  PROF_ADD_COUNTER(GP_ASSW_DECODE_CANDIDATES_RETRIVED,num_candidates);
  PROFILE_STOP(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
}
void filtering_candidates_decode_text_filtering_positions_buffered(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    const uint64_t buffered_positions_offset,
    pattern_t* const pattern,
    region_search_t* const region_search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_offset_begin) {
  PROFILE_START(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  // Parameters
  filtering_position_buffered_t* const positions_buffered =
      filtering_candidates_buffered->positions_buffered + buffered_positions_offset;
  const bool exact_match = (region_search->begin==0 && region_search->end==pattern->key_length);
  const uint64_t region_lo = region_search->lo;
  const uint64_t region_hi = region_search->hi;
  const uint64_t num_candidates = region_hi-region_lo;
  // Add all candidate positions
  uint64_t i;
  for (i=0;i<num_candidates;++i) {
    // Allocate position
    filtering_position_t* const fposition = filtering_candidates_allocate_position(filtering_candidates);
    // Configure
    fposition->source_region_begin = positions_buffered[i].source_region_begin;
    fposition->source_region_end = positions_buffered[i].source_region_end;
    fposition->source_region_error = 0;
    // Retrieve SA-sample
    filtering_candidates_decode_retrieve_text_sample(filtering_candidates,
        fposition,i,region_lo,gpu_buffer_fmi_decode,buffer_offset_begin);
    // Adjust Position
    filtering_candidates_compute_text_coordinates(filtering_candidates,fposition,pattern);
    fposition->align_distance = exact_match ? 0 : ALIGN_DISTANCE_UNKNOWN;
  }
  PROF_ADD_COUNTER(GP_FC_DECODE_POSITIONS,num_candidates);
  PROF_ADD_COUNTER(GP_ASSW_DECODE_CANDIDATES_RETRIVED,num_candidates);
  PROFILE_STOP(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
}
/*
 * Process Candidates Buffered (from GPU-Buffer)
 */
void filtering_candidates_process_candidates_buffered(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const bool compose_alignment_regions) {
  PROFILE_START(GP_FC_PROCESS_CANDIDATES,PROFILE_LEVEL);
  // Retrieve total candidate positions
  PROF_ADD_COUNTER(GP_CANDIDATE_POSITIONS,filtering_candidates_get_num_positions(filtering_candidates));
  // Compose matching-regions into candidate regions (also filter out duplicated positions or already checked)
  PROFILE_START(GP_FC_COMPOSE_REGIONS,PROFILE_LEVEL);
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  filtering_candidates_compose_filtering_regions(
      filtering_candidates,pattern,compose_alignment_regions && !search_parameters->force_full_swg);
  PROFILE_STOP(GP_FC_COMPOSE_REGIONS,PROFILE_LEVEL);
  PROF_ADD_COUNTER(GP_CANDIDATE_REGIONS,filtering_candidates_get_num_regions(filtering_candidates));
  // Return total candidate regions
  PROFILE_STOP(GP_FC_PROCESS_CANDIDATES,PROFILE_LEVEL);
}
