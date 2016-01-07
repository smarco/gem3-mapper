/*
 * PROJECT: GEMMapper
 * FILE: filtering_candidates_process_buffered.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "filtering_candidates_process_buffered.h"
#include "filtering_candidates_process.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Decode Candidates Helpers (Buffered/Batch)
 */
void filtering_candidates_decode_retrieve_sa_sample(
    filtering_position_t* const filtering_position,const uint64_t position,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,const uint64_t buffer_offset_begin,
    const fm_index_t* const fm_index,const uint64_t region_lo) {
  // Retrieve decoded position & fetch sample
  uint64_t bwt_sampled_position, lf_steps;
  gpu_buffer_fmi_decode_get_result(
      gpu_buffer_fmi_decode,buffer_offset_begin+position,&bwt_sampled_position,&lf_steps);
  if (bwt_sampled_position != -1) {
    // Recover Sample
    fm_index_retrieve_sa_sample(fm_index,bwt_sampled_position,
        lf_steps,&filtering_position->region_text_position);
  } else {
    // Re-Decode (GPU decode failed)
    PROF_INC_COUNTER(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL_TOTAL);
    PROFILE_START(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL,PROFILE_LEVEL);
    filtering_position->region_text_position = fm_index_decode(fm_index,region_lo+position);
    PROFILE_STOP(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL,PROFILE_LEVEL);
  }
  // DEBUG
#ifdef CUDA_CHECK_BUFFERED_DECODE_POSITIONS
  const uint64_t region_text_position = fm_index_decode(fm_index,region_lo+position);
  gem_cond_fatal_error_msg(filtering_position->region_text_position!=region_text_position,
      "Filtering.Candidates.Process.Buffered. Check decoded position failed (%lu!=%lu)",
      filtering_position->region_text_position,region_text_position);
#endif
}
void filtering_candidates_decode_batch_retrieve_sampled_position(
    fc_batch_decode_candidate* const batch_candidate,const uint64_t current_position,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,const uint64_t buffer_offset_begin,
    const fm_index_t* const fm_index,const uint64_t region_lo) {
  // Parameters
  const bwt_t* const bwt = fm_index->bwt;
  // Retrieve sampled-position from buffer
  gpu_buffer_fmi_decode_get_result(
      gpu_buffer_fmi_decode,buffer_offset_begin+current_position,
      &batch_candidate->index_position,&batch_candidate->distance);
  // If buffered-decode failed, retrieve sampled-position from index
  if (batch_candidate->index_position == -1) {
    PROF_INC_COUNTER(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL_TOTAL);
    PROFILE_START(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL,PROFILE_LEVEL);
    fm_index_retrieve_bwt_sampled(fm_index,region_lo+current_position,
        &batch_candidate->index_position,&batch_candidate->distance);
    PROFILE_STOP(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL,PROFILE_LEVEL);
  }
  // Prefetch BM-sampled
  bwt_prefetch(bwt,batch_candidate->index_position,&batch_candidate->bwt_block_locator);
}
void filtering_candidates_decode_batch_retrieve_bm_sampled(
    fc_batch_decode_candidate* const batch_candidate,const fm_index_t* const fm_index) {
  // Parameters
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
    fc_batch_decode_candidate* const batch_candidate,const uint64_t current_position,
    filtering_position_t* const filtering_position,const fm_index_t* const fm_index,
    const uint64_t region_lo) {
  // Parameters
  const sampled_sa_t* const sampled_sa = fm_index->sampled_sa;
  const uint64_t bwt_length = fm_index_get_length(fm_index);
  // Retrieve sampled-position
  const uint64_t sampled_position = sampled_sa_get_sample(sampled_sa,batch_candidate->index_position);
  filtering_position->region_text_position = (sampled_position + batch_candidate->distance) % bwt_length;
  // DEBUG
#ifdef CUDA_CHECK_BUFFERED_DECODE_POSITIONS
  const uint64_t region_text_position = fm_index_decode(fm_index,region_lo+current_position);
  gem_cond_fatal_error_msg(filtering_position->region_text_position!=region_text_position,
      "Filtering.Candidates.Process.Buffered. Check decoded position failed (%lu!=%lu)",
      filtering_position->region_text_position,region_text_position);
#endif
}
/*
 * Decode Candidates Buffered (from GPU-Buffer)
 */
void filtering_candidates_decode_filtering_positions_buffered(
    filtering_candidates_t* const filtering_candidates,const locator_t* const locator,
    archive_text_t* const archive_text,const fm_index_t* const fm_index,
    const uint64_t region_begin,const uint64_t region_end,
    const uint64_t region_lo,const uint64_t region_hi,
    const uint64_t key_length,const uint64_t boundary_error,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,const uint64_t buffer_offset_begin) {
  PROFILE_START(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  // Reserve
  vector_t* const filtering_positions = filtering_candidates->filtering_positions;
  const uint64_t num_candidates = region_hi-region_lo;
  vector_reserve_additional(filtering_positions,num_candidates);
  filtering_position_t* const filtering_position = vector_get_free_elm(filtering_positions,filtering_position_t);
  // Add all candidate positions
  uint64_t i;
  for (i=0;i<num_candidates;++i) {
    filtering_position_t* const fposition = filtering_position+i;
    // Retrieve SA-sample
    filtering_candidates_decode_retrieve_sa_sample(fposition,i,
        gpu_buffer_fmi_decode,buffer_offset_begin,fm_index,region_lo);
    // Locate Position
    fposition->locator_interval =
        locator_lookup_interval(locator,fposition->region_text_position);
    // Adjust Position
    filtering_candidates_adjust_filtering_position(fposition,
        archive_text,region_begin,key_length-region_end,boundary_error);
  }
  // Add used
  PROF_ADD_COUNTER(GP_FC_DECODE_POSITIONS,num_candidates);
  vector_add_used(filtering_positions,num_candidates);
  PROFILE_STOP(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
}
void filtering_candidates_decode_filtering_positions_buffered_prefetched(
    filtering_candidates_t* const filtering_candidates,const locator_t* const locator,
    archive_text_t* const archive_text,const fm_index_t* const fm_index,
    const uint64_t region_begin,const uint64_t region_end,
    const uint64_t region_lo,const uint64_t region_hi,
    const uint64_t key_length,const uint64_t boundary_error,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,const uint64_t buffer_offset_begin) {
  PROFILE_START(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
  // Reserve candidate positions
  vector_t* const filtering_positions = filtering_candidates->filtering_positions;
  const uint64_t num_candidates = region_hi-region_lo;
  vector_reserve_additional(filtering_positions,num_candidates);
  filtering_position_t* const filtering_position = vector_get_free_elm(filtering_positions,filtering_position_t);
  // Retrieve decoded positions
  fc_batch_decode_candidate batch[DECODE_NUM_POSITIONS_PREFETCHED];
  uint64_t num_left_positions = num_candidates;
  uint64_t i, current_position = 0;
  while (num_left_positions > 0) {
    const uint64_t batch_size = MIN(num_left_positions,DECODE_NUM_POSITIONS_PREFETCHED);
    // Retrieve sampled-position & prefetch BM-sampled
    for (i=0;i<batch_size;++i) {
      filtering_candidates_decode_batch_retrieve_sampled_position(
          batch+i,current_position+i,gpu_buffer_fmi_decode,
          buffer_offset_begin,fm_index,region_lo);
    }
    // Retrieve BM-sampled & prefetch SA-sample
    for (i=0;i<batch_size;++i) {
      filtering_candidates_decode_batch_retrieve_bm_sampled(batch+i,fm_index);
    }
    // Retrieve SA-sample, locate position & adjust
    for (i=0;i<batch_size;++i) {
      // Retrieve SA-sample
      filtering_position_t* const fposition = filtering_position+current_position+i;
      filtering_candidates_decode_batch_retrieve_sa_sample(
          batch+i,current_position+i,fposition,fm_index,region_lo);
      // Locate Position
      fposition->locator_interval = locator_lookup_interval(locator,fposition->region_text_position);
      // Adjust Position
      filtering_candidates_adjust_filtering_position(fposition,
          archive_text,region_begin,key_length-region_begin,boundary_error);
    }
    // Next batch
    current_position = current_position + batch_size;
    num_left_positions -= batch_size;
  }
  // Add used
  PROF_ADD_COUNTER(GP_FC_DECODE_POSITIONS,num_candidates);
  vector_add_used(filtering_positions,num_candidates);
  PROFILE_STOP(GP_FC_DECODE_CANDIDATES_BUFFERED,PROFILE_LEVEL);
}
/*
 * Process Candidates Buffered (from GPU-Buffer)
 */
void filtering_candidates_process_candidates_buffered(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,const pattern_t* const pattern,
    const as_parameters_t* const as_parameters,
    const bool compose_region_chaining,mm_stack_t* const mm_stack) {
  PROFILE_START(GP_FC_PROCESS_CANDIDATES,PROFILE_LEVEL);
  // Retrieve total candidate positions
  PROF_ADD_COUNTER(GP_CANDIDATE_POSITIONS,vector_get_used(filtering_candidates->filtering_positions));
  // Compose matching regions into candidate regions
  //   (also filter out duplicated positions or already checked)
  PROFILE_START(GP_FC_COMPOSE_REGIONS,PROFILE_LEVEL);
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  const uint64_t key_length = pattern->key_length;
  filtering_candidates_compose_filtering_regions(
      filtering_candidates,key_length,pattern->max_effective_bandwidth,
      compose_region_chaining && search_parameters->alignment_scaffolding,mm_stack);
  PROFILE_STOP(GP_FC_COMPOSE_REGIONS,PROFILE_LEVEL);
  PROF_ADD_COUNTER(GP_CANDIDATE_REGIONS,vector_get_used(filtering_candidates->filtering_regions));
  // Return total candidate regions
  PROFILE_STOP(GP_FC_PROCESS_CANDIDATES,PROFILE_LEVEL);
}
