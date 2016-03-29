/*
 * PROJECT: GEMMapper
 * FILE: mapper_cuda_se.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper/mapper_cuda_se.h"
#include "mapper/mapper.h"
#include "archive/archive_search_se.h"
#include "io/input_file.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * Check Occupancy
 */
bool mapper_se_cuda_stage_read_input_sequences_exhausted(mapper_cuda_search_t* const mapper_search) {
  // Check pending search
  if (mapper_search->pending_search_region_profile_end1!=NULL) return false;
  // Check end_of_block
  if (!buffered_input_file_eob(mapper_search->buffered_fasta_input_end1)) return false;
  // Reload buffer
  if (buffered_input_file_reload(mapper_search->buffered_fasta_input_end1,0)==INPUT_STATUS_EOF) return true;
  // Clear pipeline (release intermediate memory & start pipeline fresh)
  search_pipeline_clear(mapper_search->search_pipeline);
  return false;
}
bool mapper_se_cuda_stage_region_profile_output_exhausted(mapper_cuda_search_t* const mapper_search) {
  // Check Stage Region-Profile
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  if (!search_stage_region_profile_retrieve_finished(search_pipeline->stage_region_profile)) return false;
  // Check pending search
  if (mapper_search->pending_search_decode_candidates_end1!=NULL) return false;
  // Exhausted
  return true;
}
bool mapper_se_cuda_stage_decode_candidates_output_exhausted(mapper_cuda_search_t* const mapper_search) {
  // Check Stage Decode-Candidates
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  if (!search_stage_decode_candidates_retrieve_finished(search_pipeline->stage_decode_candidates)) return false;
  // Check pending search
  if (mapper_search->pending_search_verify_candidates_end1!=NULL) return false;
  // Exhausted
  return true;
}
/*
 * Region Profile (until all buffers are filled or input-block exhausted)
 *   Read from input-block,
 *   Generate region-profile partition
 *   Send to CUDA-region-profile
 */
void mapper_se_cuda_region_profile(mapper_cuda_search_t* const mapper_search) {
  PROFILE_START(GP_MAPPER_CUDA_SE_REGION_PROFILE,PROFILE_LEVEL);
  // Parameters
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  search_stage_region_profile_t* const stage_region_profile = search_pipeline->stage_region_profile;
  archive_search_t* archive_search = NULL;
  // Reschedule search (that couldn't fit into the buffer)
  if (mapper_search->pending_search_region_profile_end1!=NULL) {
    search_stage_region_profile_send_se_search(stage_region_profile,mapper_search->pending_search_region_profile_end1);
    mapper_search->pending_search_region_profile_end1 = NULL;
  }
  // Generation. Keep processing the current input-block
  while (!buffered_input_file_eob(mapper_search->buffered_fasta_input_end1)) {
    // Request a clean archive-search
    search_pipeline_allocate_se(search_pipeline,&archive_search);
    // Parse Sequence
    const error_code_t error_code = input_fasta_parse_sequence(
        mapper_search->buffered_fasta_input_end1,archive_search_get_sequence(archive_search),
        parameters->io.fastq_strictly_normalized,false);
    gem_cond_fatal_error(error_code==INPUT_STATUS_FAIL,MAPPER_CUDA_ERROR_PARSING);
    PROF_INC_COUNTER(GP_MAPPER_NUM_READS);
    // Generate Candidates (Search into the archive)
    archive_search_se_stepwise_init_search(archive_search);
    archive_search_se_stepwise_region_profile_generate_static(archive_search);
    // Send search to GPU region-profile
    const bool search_sent = search_stage_region_profile_send_se_search(stage_region_profile,archive_search);
    if (!search_sent) {
      mapper_search->pending_search_region_profile_end1 = archive_search; // Pending Search
      break;
    }
  }
  PROFILE_STOP(GP_MAPPER_CUDA_SE_REGION_PROFILE,PROFILE_LEVEL);
}
/*
 * Decode Candidates (until all buffers are filled or stage region-profile exhausted)
 *   Read from stage Region-Profile
 *   Generate Decode-Candidates
 *   Send to CUDA Decode-Candidates
 */
void mapper_se_cuda_decode_candidates(mapper_cuda_search_t* const mapper_search) {
  PROFILE_START(GP_MAPPER_CUDA_SE_DECODE_CANDIDATES,PROFILE_LEVEL);
  // Parameters
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  search_stage_region_profile_t* const stage_region_profile = search_pipeline->stage_region_profile;
  search_stage_decode_candidates_t* const stage_decode_candidates = search_pipeline->stage_decode_candidates;
  archive_search_t* archive_search = NULL;
  // Reschedule search (that couldn't fit into the buffer)
  if (mapper_search->pending_search_decode_candidates_end1!=NULL) {
    search_stage_decode_candidates_send_se_search(
        stage_decode_candidates,mapper_search->pending_search_decode_candidates_end1);
    mapper_search->pending_search_decode_candidates_end1 = NULL;
  }
  // Read from stage Region-Profile
  bool pending_searches;
  while ((pending_searches=search_stage_region_profile_retrieve_se_search(stage_region_profile,&archive_search))) {
    // Send to CUDA Decode-Candidates
    const bool search_sent = search_stage_decode_candidates_send_se_search(stage_decode_candidates,archive_search);
    if (!search_sent) {
      mapper_search->pending_search_decode_candidates_end1 = archive_search; // Pending Search
      break;
    }
  }
  // Clean
  if (!pending_searches) {
    search_stage_region_profile_clear(stage_region_profile,NULL);
  }
  PROFILE_STOP(GP_MAPPER_CUDA_SE_DECODE_CANDIDATES,PROFILE_LEVEL);
}
/*
 * Generate Candidates (until all buffers are filled or input-block exhausted)
 *   Read from input-block,
 *   Generate region-profile partition (adaptive) + decode-candidates
 *   Send to CUDA Decode-Candidates
 */
void mapper_se_cuda_generate_candidates(mapper_cuda_search_t* const mapper_search) {
  PROFILE_START(GP_MAPPER_CUDA_SE_DECODE_CANDIDATES,PROFILE_LEVEL);
  // Parameters
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  search_stage_decode_candidates_t* const stage_decode_candidates = search_pipeline->stage_decode_candidates;
  archive_search_t* archive_search = NULL;
  // Reschedule search (that couldn't fit into the buffer)
  if (mapper_search->pending_search_decode_candidates_end1!=NULL) {
    search_stage_decode_candidates_send_se_search(
        stage_decode_candidates,mapper_search->pending_search_decode_candidates_end1);
    mapper_search->pending_search_decode_candidates_end1 = NULL;
  }
  // Generation. Keep processing the current input-block
  while (!buffered_input_file_eob(mapper_search->buffered_fasta_input_end1)) {
    // Request a clean archive-search
    search_pipeline_allocate_se(search_pipeline,&archive_search);
    search_stage_decode_candidates_prepare_se_search(stage_decode_candidates,archive_search);
    // Parse Sequence
    const error_code_t error_code = input_fasta_parse_sequence(
        mapper_search->buffered_fasta_input_end1,archive_search_get_sequence(archive_search),
        parameters->io.fastq_strictly_normalized,false);
    gem_cond_fatal_error(error_code==INPUT_STATUS_FAIL,MAPPER_CUDA_ERROR_PARSING);
    PROF_INC_COUNTER(GP_MAPPER_NUM_READS);
    // Generate Candidates (Search into the archive)
    archive_search_se_stepwise_init_search(archive_search);
    archive_search_se_stepwise_region_profile_generate_adaptive(archive_search);
    // Send search to GPU region-profile
    const bool search_sent = search_stage_decode_candidates_send_se_search(stage_decode_candidates,archive_search);
    if (!search_sent) {
      mapper_search->pending_search_decode_candidates_end1 = archive_search; // Pending Search
      break;
    }
  }
  PROFILE_STOP(GP_MAPPER_CUDA_SE_DECODE_CANDIDATES,PROFILE_LEVEL);
}
/*
 * Verify Candidates (until all buffers are filled or stage Decode-Candidates exhausted)
 *   Read from stage Decode-Candidates
 *   Generate Verify-Candidates
 *   Send to CUDA Verify-Candidates
 */
void mapper_se_cuda_verify_candidates(mapper_cuda_search_t* const mapper_search) {
  PROFILE_START(GP_MAPPER_CUDA_SE_VERIFY_CANDIDATES,PROFILE_LEVEL);
  // Parameters
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  search_stage_decode_candidates_t* const stage_decode_candidates = search_pipeline->stage_decode_candidates;
  search_stage_verify_candidates_t* const stage_verify_candidates = search_pipeline->stage_verify_candidates;
  archive_search_t* archive_search = NULL;
  // Reschedule search (that couldn't fit into the buffer)
  if (mapper_search->pending_search_verify_candidates_end1!=NULL) {
    search_stage_verify_candidates_send_se_search(
        stage_verify_candidates,mapper_search->pending_search_verify_candidates_end1);
    mapper_search->pending_search_verify_candidates_end1 = NULL;
  }
  // Read from stage Decode-Candidates,
  bool pending_searches;
  while ((pending_searches=search_stage_decode_candidates_retrieve_se_search(stage_decode_candidates,&archive_search))) {
    // Generate Verify-Candidates
    // /* NOP */ archive_search_se_stepwise_verify_candidates_generate(archive_search);
    // Send to CUDA Verify-Candidates
    const bool search_sent = search_stage_verify_candidates_send_se_search(stage_verify_candidates,archive_search);
    if (!search_sent) {
      mapper_search->pending_search_verify_candidates_end1 = archive_search; // Pending Search
      break;
    }
  }
  // Clean
  if (!pending_searches) {
    search_stage_decode_candidates_clear(stage_decode_candidates,NULL);
  }
  PROFILE_STOP(GP_MAPPER_CUDA_SE_VERIFY_CANDIDATES,PROFILE_LEVEL);
}
/*
 * Finish Search (until stage Verify-Candidates exhausted)
 *   Read from stage Verify-Candidates
 *   Finish Search
 *   Output Matches
 */
void mapper_se_cuda_finish_search(mapper_cuda_search_t* const mapper_search) {
  PROFILE_START(GP_MAPPER_CUDA_SE_FINISH_SEARCH,PROFILE_LEVEL);
  // Parameters
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  search_stage_verify_candidates_t* const stage_verify_candidates = search_pipeline->stage_verify_candidates;
  archive_search_t* archive_search = NULL;
  // Process all search-groups generated
  while (search_stage_verify_candidates_retrieve_se_search(stage_verify_candidates,&archive_search)) {
    // Finish Search
    archive_search_se_stepwise_finish_search(archive_search,stage_verify_candidates->matches); // Finish search
    // Output Matches
    mapper_SE_output_matches(parameters,mapper_search->buffered_output_file,
        archive_search,stage_verify_candidates->matches,mapper_search->mapping_stats);
    // Update processed
    if (++mapper_search->reads_processed == MAPPER_TICKER_STEP) {
      ticker_update_mutex(mapper_search->ticker,mapper_search->reads_processed);
      mapper_search->reads_processed=0;
    }
  }
  // Clean
  search_stage_verify_candidates_clear(stage_verify_candidates,search_pipeline->archive_search_cache);
  PROFILE_STOP(GP_MAPPER_CUDA_SE_FINISH_SEARCH,PROFILE_LEVEL);
}
/*
 * Mapper SE-CUDA
 */
void* mapper_cuda_se_thread(mapper_cuda_search_t* const mapper_search) {
  // GEM-thread error handler
  gem_thread_register_id(mapper_search->thread_id+1);
  PROFILE_START(GP_MAPPER_CUDA_SE,PROFILE_LEVEL);
  // Parameters
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  const mapper_parameters_cuda_t* const cuda_parameters = &parameters->cuda;
  // Create new buffered reader/writer
  mapper_SE_prepare_io_buffers(parameters,cuda_parameters->input_buffer_size,
      &mapper_search->buffered_fasta_input_end1,&mapper_search->buffered_output_file);
  // Create search-pipeline & initialize matches
  mapper_search->search_pipeline = search_pipeline_new(parameters,
      mapper_search->gpu_buffer_collection,mapper_search->gpu_buffers_offset,false);
  mapper_search->pending_search_region_profile_end1 = NULL;
  mapper_search->pending_search_decode_candidates_end1 = NULL;
  mapper_search->pending_search_verify_candidates_end1 = NULL;
  // FASTA/FASTQ reading loop
  mm_stack_t* const mm_stack = mapper_search->search_pipeline->mm_stack;
  mapper_search->reads_processed = 0;
  while (!mapper_se_cuda_stage_read_input_sequences_exhausted(mapper_search)) {
#ifdef MAPPER_CUDA_ADAPTIVE_REGION_PROFILE
    // Decode Candidates
    mapper_se_cuda_generate_candidates(mapper_search);
    mm_stack_push_state(mm_stack);
    do {
      // Verify Candidates
      mapper_se_cuda_verify_candidates(mapper_search);
      // Finish Search
      mm_stack_push_state(mm_stack);
      mapper_se_cuda_finish_search(mapper_search);
      mm_stack_pop_state(mm_stack);
    } while (!mapper_se_cuda_stage_decode_candidates_output_exhausted(mapper_search));
    mm_stack_pop_state(mm_stack);
#else
    // Region Profile
    mapper_se_cuda_region_profile(mapper_search);
    mm_stack_push_state(mm_stack);
    do {
      // Decode Candidates
      mapper_se_cuda_decode_candidates(mapper_search);
      mm_stack_push_state(mm_stack);
      do {
        // Verify Candidates
        mapper_se_cuda_verify_candidates(mapper_search);
        // Finish Search
        mm_stack_push_state(mm_stack);
        mapper_se_cuda_finish_search(mapper_search);
        mm_stack_pop_state(mm_stack);
      } while (!mapper_se_cuda_stage_decode_candidates_output_exhausted(mapper_search));
      mm_stack_pop_state(mm_stack);
    } while (!mapper_se_cuda_stage_region_profile_output_exhausted(mapper_search));
    mm_stack_pop_state(mm_stack);
#endif
  }
  // Clean up
  ticker_update_mutex(mapper_search->ticker,mapper_search->reads_processed); // Update processed
  search_pipeline_delete(mapper_search->search_pipeline);
  buffered_input_file_close(mapper_search->buffered_fasta_input_end1);
  buffered_output_file_close(mapper_search->buffered_output_file);
  PROFILE_STOP(GP_MAPPER_CUDA_SE,PROFILE_LEVEL);
  pthread_exit(0);
}
