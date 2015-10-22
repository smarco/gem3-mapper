/*
 * PROJECT: GEMMapper
 * FILE: mapper_cuda_se.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper_cuda_se.h"
#include "mapper.h"
#include "archive_search_se.h"
#include "input_file.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * Check the end_of_block (We cannot reload input-buffer until all
 *   the searches of the previous block are solved)
 */
GEM_INLINE bool mapper_input_sequences_exhausted(mapper_cuda_search_t* const mapper_search) {
  if (!buffered_input_file_eob(mapper_search->buffered_fasta_input)) return false;
  if (search_group_verify_candidates_buffer_is_empty(mapper_search->search_group->search_group_vc)) {
    if (buffered_input_file_reload__dump_attached(mapper_search->buffered_fasta_input)==INPUT_STATUS_EOF) {
      return true;
    }
  }
  return false;
}
/*
 * Read from input-block, generate candidates & send to CUDA-verification
 *  (until all buffers are filled or input-block exhausted)
 */
GEM_INLINE void mapper_se_cuda_generate_candidates(mapper_cuda_search_t* const mapper_search) {
  PROFILE_START(GP_MAPPER_CUDA_THREAD_GENERATING,PROFILE_LEVEL);
  // Parameters
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  search_group_t* const search_group = mapper_search->search_group;
  search_group_verify_candidates_t* const search_group_vc = search_group->search_group_vc;
  // Generation. Keep processing the current input-block
  while (!buffered_input_file_eob(mapper_search->buffered_fasta_input)) {
    // Request a clean archive-search
    mapper_search->archive_search_generate_candidates = search_group_allocate_se(search_group);
    // Parse Sequence
    const error_code_t error_code = input_fasta_parse_sequence(
        mapper_search->buffered_fasta_input,archive_search_get_sequence(mapper_search->archive_search_generate_candidates),
        parameters->io.fastq_strictly_normalized,parameters->io.fastq_try_recovery,false);
    gem_cond_fatal_error(error_code==INPUT_STATUS_FAIL,MAPPER_CUDA_ERROR_PARSING);
    // Generate Candidates (Search into the archive)
    archive_search_se_stepwise_init_search(mapper_search->archive_search_generate_candidates);
    archive_search_se_stepwise_generate_candidates(mapper_search->archive_search_generate_candidates);
    // Add archive-search to group (Put candidates in buffer)
    if (!search_group_verify_candidates_se_add_search(search_group_vc,mapper_search->archive_search_generate_candidates)) {
      break; // Go to select-candidates
    }
    mapper_search->archive_search_generate_candidates = NULL; // Last archive-search fitted into the BPM-buffer
  }
  PROFILE_STOP(GP_MAPPER_CUDA_THREAD_GENERATING,PROFILE_LEVEL);
}
/*
 *
 */
GEM_INLINE void mapper_se_cuda_finish_search(mapper_cuda_search_t* const mapper_search) {
  PROFILE_START(GP_MAPPER_CUDA_THREAD_SELECTING,PROFILE_LEVEL);
  // Parameters
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  search_group_t* const search_group = mapper_search->search_group;
  search_group_verify_candidates_t* const search_group_vc = search_group->search_group_vc;
  text_collection_t* const text_collection = search_group_get_text_collection(search_group);
  matches_t* const matches = mapper_search->matches;
  // Process all search-groups generated
  search_group_verify_candidates_retrieve_begin(search_group->search_group_vc); // Start retrieving
  while (search_group_verify_candidates_se_get_search(
      search_group_vc,&mapper_search->archive_search_verify_candidates,text_collection,matches)) {
    // Finish Search
    archive_search_se_stepwise_finish_search(mapper_search->archive_search_verify_candidates,matches);
    archive_select_se_matches(mapper_search->archive_search_verify_candidates,false,matches);
    // Output matches
    mapper_SE_output_matches(parameters,mapper_search->buffered_output_file,
        mapper_search->archive_search_verify_candidates,matches,mapper_search->mapping_stats);
    // Update processed
    if (++mapper_search->reads_processed == MAPPER_TICKER_STEP) {
      ticker_update_mutex(mapper_search->ticker,mapper_search->reads_processed);
      mapper_search->reads_processed=0;
    }
  }
  // Reset search-group
  search_group_clear(search_group);
  PROFILE_STOP(GP_MAPPER_CUDA_THREAD_SELECTING,PROFILE_LEVEL);
}
/*
 * Reschedule pending searches
 */
GEM_INLINE void mapper_se_cuda_reschedule_pending_searches(mapper_cuda_search_t* const mapper_search) {
  // Check if the last archive-search couldn't fit into the BPM-buffer
  if (mapper_search->archive_search_generate_candidates!=NULL) {
    PROFILE_START(GP_MAPPER_CUDA_THREAD_RESTART_UNFIT,PROFILE_LEVEL);
    archive_search_se_stepwise_init_search(mapper_search->archive_search_generate_candidates);
    archive_search_se_stepwise_generate_candidates(mapper_search->archive_search_generate_candidates);
    search_group_verify_candidates_se_add_search(
        mapper_search->search_group->search_group_vc,mapper_search->archive_search_generate_candidates);
    mapper_search->archive_search_generate_candidates = NULL;
    PROFILE_STOP(GP_MAPPER_CUDA_THREAD_RESTART_UNFIT,PROFILE_LEVEL);
  }
}
/*
 * Mapper SE-CUDA
 */
void* mapper_cuda_se_thread(mapper_cuda_search_t* const mapper_search) {
  // GEM-thread error handler
  gem_thread_register_id(mapper_search->thread_id+1);
  PROFILE_START(GP_MAPPER_CUDA_THREAD,PROFILE_LEVEL);
  // Parameters
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  const mapper_parameters_cuda_t* const cuda_parameters = &parameters->cuda;
  search_group_t* const search_group = mapper_search->search_group;
  // Create new buffered reader/writer
  mapper_search->buffered_fasta_input = buffered_input_file_new(parameters->input_file,cuda_parameters->input_buffer_lines);
  mapper_search->buffered_output_file = buffered_output_file_new(parameters->output_file);
  buffered_input_file_attach_buffered_output(mapper_search->buffered_fasta_input,mapper_search->buffered_output_file);
  // Init search-group & matches
  search_group_init(search_group);
  mapper_search->matches = matches_new();
  matches_configure(mapper_search->matches,search_group_get_text_collection(search_group));
  // FASTA/FASTQ reading loop
  mapper_search->reads_processed = 0;
  while (true) {
    // Check input-block
    if (mapper_input_sequences_exhausted(mapper_search)) break;
    // Read sequence, generate candidates & send to CUDA-verification
    mapper_se_cuda_generate_candidates(mapper_search);
    // Retrieve CUDA-verification results & finish the search
    mapper_se_cuda_finish_search(mapper_search);
    // Reschedule searches that couldn't fit into any stage-buffer
    mapper_se_cuda_reschedule_pending_searches(mapper_search);
  }
  // Clean up
  ticker_update_mutex(mapper_search->ticker,mapper_search->reads_processed); // Update processed
  buffered_input_file_close(mapper_search->buffered_fasta_input);
  buffered_output_file_close(mapper_search->buffered_output_file);
  matches_delete(mapper_search->matches);
  PROFILE_STOP(GP_MAPPER_CUDA_THREAD,PROFILE_LEVEL);
  pthread_exit(0);
}
