/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2013-2017 by Alejandro Chacon <alejandro.chacond@gmail.com>
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
 *            Alejandro Chacon <alejandro.chacond@gmail.com>
 * DESCRIPTION:
 *   Filtering candidates module provides functions to process all
 *   bwt-encoded positions (originated from a candidate generation process)
 *   and compose them into filtering-regions (using decoded text coordinates)
 *   This "buffered" module operates in batches of filtering-candidates and
 *   makes use of GPU-buffers to offloads the decoding of positions to a GPU
 */

#include "filtering/candidates/filtering_candidates_buffered_process.h"
#include "align/alignment.h"
#include "filtering/candidates/filtering_candidates_process.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Decode SA-Positions Buffered (from GPU-Buffer)
 */
void filtering_candidates_buffered_decode_sa_positions(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern) {
  PROFILE_START(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  // Parameters
  const uint64_t num_positions_buffered = filtering_candidates_buffered->num_positions;
  filtering_position_t** const positions_buffered = filtering_candidates_buffered->positions;
  fm_index_t* const fm_index = filtering_candidates->archive->fm_index;
  // Add all candidate positions
  uint64_t i;
  for (i=0;i<num_positions_buffered;++i) {
    // Add position
    filtering_position_t* const fposition = positions_buffered[i];
    vector_insert(filtering_candidates->filtering_positions,fposition,filtering_position_t*);
    // Decode position
    fposition->region_text_position = fm_index_decode(fm_index,fposition->region_index_position);
    // Adjust Position
    filtering_candidates_compute_text_coordinates(filtering_candidates,fposition,pattern);
  }
  PROF_ADD_COUNTER(GP_FC_DECODE_POSITIONS,num_positions_buffered);
  PROF_ADD_COUNTER(GP_ASSW_DECODE_CANDIDATES_RETRIVED,num_positions_buffered);
  PROFILE_STOP(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
}
/*
 * Decode sampled SA-Positions Buffered (from GPU-Buffer)
 */
void filtering_candidates_decode_batch_retrieve_sampled_position(
    filtering_candidates_t* const filtering_candidates,
    filtering_position_t* const filtering_position,
    fc_batch_decode_candidate* const batch_candidate,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t gpu_buffer_fmi_decode_offset) {
  // Parameters
  fm_index_t* const fm_index = filtering_candidates->archive->fm_index;
  const bwt_t* const bwt = fm_index->bwt;
  // Retrieve sampled-position from buffer
  gpu_buffer_fmi_decode_get_position_sa(
      gpu_buffer_fmi_decode,gpu_buffer_fmi_decode_offset,
      &batch_candidate->index_position,&batch_candidate->distance);
  // If buffered-decode failed, retrieve sampled-position from index
  if (batch_candidate->index_position == -1) {
    PROF_INC_COUNTER(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL_TOTAL);
    PROFILE_START(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL,PROFILE_LEVEL);
    fm_index_retrieve_bwt_sampled(fm_index,filtering_position->region_index_position,
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
    filtering_position_t* const filtering_position,
    fc_batch_decode_candidate* const batch_candidate) {
  // Parameters
  fm_index_t* const fm_index = filtering_candidates->archive->fm_index;
  const sampled_sa_t* const sampled_sa = fm_index->sampled_sa;
  const uint64_t bwt_length = fm_index_get_length(fm_index);
  // Retrieve sampled-position
  const uint64_t sampled_position = sampled_sa_get_sample(sampled_sa,batch_candidate->index_position);
  filtering_position->region_text_position = (sampled_position + batch_candidate->distance) % bwt_length;
  // DEBUG
#ifdef GPU_CHECK_DECODE_POSITIONS
  const uint64_t region_text_position = fm_index_decode(fm_index,filtering_position->region_index_position);
  gem_cond_fatal_error_msg(filtering_position->region_text_position!=region_text_position,
      "Filtering.Candidates.Process.Buffered. Check decoded position failed (%lu!=%lu)",
      filtering_position->region_text_position,region_text_position);
#endif
}
void filtering_candidates_buffered_decode_sampled_sa_positions(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_offset_begin) {
  PROFILE_START(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  // Parameters
  const uint64_t num_positions_buffered = filtering_candidates_buffered->num_positions;
  filtering_position_t** const positions_buffered = filtering_candidates_buffered->positions;
  // Retrieve decoded positions
  fc_batch_decode_candidate batch[DECODE_NUM_POSITIONS_PREFETCHED];
  uint64_t num_left_positions = num_positions_buffered;
  uint64_t i, current_position = 0;
  while (num_left_positions > 0) {
    const uint64_t batch_size = MIN(num_left_positions,DECODE_NUM_POSITIONS_PREFETCHED);
    // Retrieve sampled-position & prefetch BM-sampled
    for (i=0;i<batch_size;++i) {
      filtering_candidates_decode_batch_retrieve_sampled_position(
          filtering_candidates,positions_buffered[current_position+i],
          batch+i,gpu_buffer_fmi_decode,buffer_offset_begin+current_position+i);
    }
    // Retrieve BM-sampled & prefetch SA-sample
    for (i=0;i<batch_size;++i) {
      filtering_candidates_decode_batch_retrieve_bm_sampled(filtering_candidates,batch+i);
    }
    // Retrieve SA-sample, locate position & adjust
    for (i=0;i<batch_size;++i) {
      // Add position
      filtering_position_t* const fposition = positions_buffered[current_position+i];
      vector_insert(filtering_candidates->filtering_positions,fposition,filtering_position_t*);
      // Retrieve SA-sample
      filtering_candidates_decode_batch_retrieve_sa_sample(filtering_candidates,fposition,batch+i);
      // Adjust Position
      filtering_candidates_compute_text_coordinates(filtering_candidates,fposition,pattern);
    }
    // Next batch
    current_position = current_position + batch_size;
    num_left_positions -= batch_size;
  }
  PROF_ADD_COUNTER(GP_FC_DECODE_POSITIONS,num_positions_buffered);
  PROF_ADD_COUNTER(GP_ASSW_DECODE_CANDIDATES_RETRIVED,num_positions_buffered);
  PROFILE_STOP(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
}
/*
 * Decode Text-Positions Buffered (from GPU-Buffer)
 */
void filtering_candidates_decode_retrieve_text_sample(
    filtering_candidates_t* const filtering_candidates,
    filtering_position_t* const filtering_position,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t gpu_buffer_fmi_decode_offset) {
  // Parameters
  fm_index_t* const fm_index = filtering_candidates->archive->fm_index;
  // Retrieve decoded position & fetch sample
  uint64_t text_position;
  gpu_buffer_fmi_decode_get_position_text(
      gpu_buffer_fmi_decode,gpu_buffer_fmi_decode_offset,&text_position);
  if (text_position == -1) {
    // Re-Decode (GPU decode failed)
    PROF_INC_COUNTER(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL_TOTAL);
    PROFILE_START(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL,PROFILE_LEVEL);
    filtering_position->region_text_position =
        fm_index_decode(fm_index,filtering_position->region_index_position);
    PROFILE_STOP(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL,PROFILE_LEVEL);
  } else {
    filtering_position->region_text_position = text_position;
  }
  // DEBUG
#ifdef GPU_CHECK_DECODE_POSITIONS
  const uint64_t region_text_position = fm_index_decode(fm_index,filtering_position->region_index_position);
  gem_cond_fatal_error_msg(filtering_position->region_text_position!=region_text_position,
      "Filtering.Candidates.Process.Buffered. Check decoded position failed (%lu!=%lu)",
      filtering_position->region_text_position,region_text_position);
#endif
}
void filtering_candidates_buffered_decode_text_positions(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_offset_begin) {
  PROFILE_START(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  // Parameters
  const uint64_t num_positions_buffered = filtering_candidates_buffered->num_positions;
  filtering_position_t** const positions_buffered = filtering_candidates_buffered->positions;
  // Add all candidate positions
  uint64_t i;
  for (i=0;i<num_positions_buffered;++i) {
    // Add position
    filtering_position_t* const fposition = positions_buffered[i];
    vector_insert(filtering_candidates->filtering_positions,fposition,filtering_position_t*);
    // Retrieve SA-sample
    filtering_candidates_decode_retrieve_text_sample(
        filtering_candidates,fposition,
        gpu_buffer_fmi_decode,buffer_offset_begin+i);
    // Adjust Position
    filtering_candidates_compute_text_coordinates(filtering_candidates,fposition,pattern);
  }
  PROF_ADD_COUNTER(GP_FC_DECODE_POSITIONS,num_positions_buffered);
  PROF_ADD_COUNTER(GP_ASSW_DECODE_CANDIDATES_RETRIVED,num_positions_buffered);
  PROFILE_STOP(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
}
/*
 * Process Candidates Buffered (from GPU-Buffer)
 */
void filtering_candidates_buffered_process_candidates(
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
      filtering_candidates,pattern,
      compose_alignment_regions && !search_parameters->alignment_force_full_swg);
  PROFILE_STOP(GP_FC_COMPOSE_REGIONS,PROFILE_LEVEL);
  PROF_ADD_COUNTER(GP_CANDIDATE_REGIONS,filtering_candidates_get_num_regions(filtering_candidates));
  // Return total candidate regions
  PROFILE_STOP(GP_FC_PROCESS_CANDIDATES,PROFILE_LEVEL);
}
