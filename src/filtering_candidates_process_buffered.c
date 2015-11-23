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
 * Add all decoded candidate-positions from gpu-buffer
 */
GEM_INLINE void filtering_candidates_decode_filtering_positions_buffered(
    filtering_candidates_t* const filtering_candidates,const locator_t* const locator,
    archive_text_t* const archive_text,const fm_index_t* const fm_index,
    const uint64_t region_begin,const uint64_t region_end,
    const uint64_t region_lo,const uint64_t region_hi,
    const uint64_t key_length,const uint64_t boundary_error,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,const uint64_t buffer_offset_begin) {
  // Reserve
  vector_t* const filtering_positions = filtering_candidates->filtering_positions;
  const uint64_t num_candidates = region_hi-region_lo;
  vector_reserve_additional(filtering_positions,num_candidates);
  filtering_position_t* filtering_position = vector_get_free_elm(filtering_positions,filtering_position_t);
  // Add all candidate positions
  uint64_t i;
  for (i=0;i<num_candidates;++i,++filtering_position) {
    // Retrieve decoded position & fetch sample
    uint64_t bwt_sampled_position, lf_steps;
    gpu_buffer_fmi_decode_get_result(gpu_buffer_fmi_decode,buffer_offset_begin+i,&bwt_sampled_position,&lf_steps);
    if (bwt_sampled_position != -1) {
      // Recover Sample
      fm_index_retrieve_sa_sample(fm_index,bwt_sampled_position,lf_steps,&filtering_position->region_text_position);
    } else {
      // Re-Decode (GPU decode failed)
      filtering_position->region_text_position = fm_index_decode(fm_index,region_lo+i);
    }
    // Locate Position
    filtering_position->locator_interval = locator_lookup_interval(locator,filtering_position->region_text_position);
    // Adjust Position
    filtering_candidates_adjust_filtering_position(filtering_position,
        archive_text,region_begin,key_length-region_end,boundary_error);
  }
  // Add used
  vector_add_used(filtering_positions,num_candidates);
}
/*
 * Process Candidates from GPU-Buffer
 */
GEM_INLINE void filtering_candidates_process_candidates_buffered(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,const pattern_t* const pattern,
    const as_parameters_t* const as_parameters,
    const bool compose_region_chaining,mm_stack_t* const mm_stack) {
  PROFILE_START(GP_FC_PROCESS_CANDIDATES,PROFILE_LEVEL);
  // Retrieve total candidate positions
  uint64_t pending_candidates = vector_get_used(filtering_candidates->filtering_positions);
  PROF_ADD_COUNTER(GP_CANDIDATE_POSITIONS,pending_candidates);
  // Compose matching regions into candidate regions
  //   (also filter out duplicated positions or already checked)
  PROFILE_START(GP_FC_COMPOSE_REGIONS,PROFILE_LEVEL);
  search_parameters_t* const search_parameters = as_parameters->search_parameters;
  const uint64_t key_length = pattern->key_length;
  pending_candidates = filtering_candidates_compose_filtering_regions(
      filtering_candidates,key_length,pattern->max_effective_bandwidth,
      compose_region_chaining && search_parameters->alignment_scaffolding,mm_stack);
  PROFILE_STOP(GP_FC_COMPOSE_REGIONS,PROFILE_LEVEL);
  PROF_ADD_COUNTER(GP_CANDIDATE_REGIONS,pending_candidates);
  // Return total candidate regions
  PROFILE_STOP(GP_FC_PROCESS_CANDIDATES,PROFILE_LEVEL);
}
