/*
 * PROJECT: GEMMapper
 * FILE: mapper_cuda_pe.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper_cuda_pe.h"
#include "archive_search_pe.h"

/*
 * Mapper PE-CUDA
 */
void* mapper_cuda_pe_thread(mapper_cuda_search_t* const mapper_search) {
//  // GEM-thread error handler
//  gem_thread_register_id(mapper_search->thread_id+1);
//  PROFILE_START(GP_MAPPER_CUDA_THREAD);
//  // Parameters
//  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
//  const mapper_parameters_cuda_t* const cuda_parameters = &parameters->cuda;
//  // Create new buffered reader/writer
//  buffered_output_file_t* const buffered_output_file = buffered_output_file_new(parameters->output_file);
//  mapper_PE_prepare_io_buffers(parameters,cuda_parameters->input_buffer_lines,
//      &mapper_search->buffered_fasta_input_end1,&mapper_search->buffered_fasta_input_end2,buffered_output_file);
//  buffered_input_file_t* const buffered_fasta_input_end1 = mapper_search->buffered_fasta_input_end1;
//  buffered_input_file_t* const buffered_fasta_input_end2 = mapper_search->buffered_fasta_input_end2;
//  // Create search-group, archive_search pointers & paired-end matches
//  archive_search_group_t* const search_group = mapper_search->search_group;
//  archive_search_group_init_bpm_buffers(search_group); // Init BPM-buffers
//  text_collection_t* const text_collection = archive_search_group_get_text_collection(search_group);
//  paired_matches_t* const paired_matches = paired_matches_new();
//  paired_matches_configure(paired_matches,text_collection);
//  archive_search_t *archive_search_generate_end1 = NULL, *archive_search_generate_end2 = NULL;
//  archive_search_t *archive_search_select_end1, *archive_search_select_end2;
//  bpm_gpu_buffer_t *bpm_gpu_buffer_end1, *bpm_gpu_buffer_end2;
//  // FASTA/FASTQ reading loop
//  error_code_t error_code = 0;
//  uint64_t reads_processed = 0;
//  while (true) {
//    PROFILE_START(GP_MAPPER_CUDA_THREAD_GENERATING);
//    // Check the end_of_block (We cannot reload input-buffer until all the searches of the previous block are solved)
//    if (buffered_input_file_eob(mapper_search->buffered_fasta_input_end1) && archive_search_group_is_empty(search_group)) {
//      error_code = mapper_PE_reload_buffers(parameters,buffered_fasta_input_end1,buffered_fasta_input_end2);
//      if (error_code==INPUT_STATUS_EOF) break;
//    }
//    /*
//     * First-stage: Generation. Keep processing the current input-block
//     */
//    while (!buffered_input_file_eob(mapper_search->buffered_fasta_input_end1)) {
//      if (buffered_input_file_eob(mapper_search->buffered_fasta_input_end2)) {
//        MAPPER_ERROR_PE_PARSE_UNSYNCH_INPUT_FILES(parameters);
//      }
//      // Request a clean archive-search
//      archive_search_group_allocate_pe(search_group,&archive_search_generate_end1,&archive_search_generate_end2);
//      // Parse Sequence
//      error_code = mapper_PE_parse_paired_sequences(parameters,
//          buffered_fasta_input_end1,buffered_fasta_input_end2,
//          archive_search_generate_end1,archive_search_generate_end2);
//      gem_cond_fatal_error(error_code==INPUT_STATUS_FAIL,MAPPER_CUDA_ERROR_PARSING);
//      // Begin Search
//      archive_search_pe_generate_candidates(archive_search_generate_end1,archive_search_generate_end2,paired_matches);
//      // Add archive-search to group (Put candidates in buffer)
//      if (!archive_search_group_add_paired_search(search_group,
//          archive_search_generate_end1,archive_search_generate_end2)) break; // Go to select-candidates
//      archive_search_generate_end1 = NULL; // Last archive-search is in BPM-buffer
//    }
//    PROFILE_STOP(GP_MAPPER_CUDA_THREAD_GENERATING);
//    /*
//     * Second-stage: Retrieval (Select & Output)
//     */
//    PROFILE_START(GP_MAPPER_CUDA_THREAD_SELECTING);
//    // Start retrieving
//    archive_search_group_retrieve_begin(search_group);
//    // Process all search-groups generated
//    while (archive_search_group_get_paired_search(search_group,
//        &archive_search_select_end1,&bpm_gpu_buffer_end1,
//        &archive_search_select_end2,&bpm_gpu_buffer_end2)) {
//      // Retrieve candidates
//      paired_matches_clear(paired_matches); // Clear Paired Matches
//      text_collection_clear(text_collection); // Clear text-collection
//      archive_search_retrieve_candidates(archive_search_select_end1,bpm_gpu_buffer_end1,paired_matches->matches_end1);
//      archive_search_retrieve_candidates(archive_search_select_end2,bpm_gpu_buffer_end2,paired_matches->matches_end2);
//      // Finish Search
//      archive_search_pe_finish_search(archive_search_select_end1,archive_search_select_end2,paired_matches);
//      archive_select_paired_matches(archive_search_select_end1,archive_search_select_end2,paired_matches);
//      // Output matches
//      mapper_PE_output_matches(parameters,buffered_output_file,
//          archive_search_select_end1,archive_search_select_end2,paired_matches,mapper_search->mapping_stats);
//      // Update processed
//      if (++reads_processed == MAPPER_TICKER_STEP) {
//        ticker_update_mutex(mapper_search->ticker,reads_processed);
//        reads_processed=0;
//      }
//    }
//    archive_search_group_clear(search_group); // Reset search-group
//    // Check if the last archive-search couldn't fit into the BPM-buffer
//    if (archive_search_generate_end1!=NULL) {
//      PROFILE_START(GP_MAPPER_CUDA_THREAD_RESTART_UNFIT);
//      archive_search_pe_generate_candidates(archive_search_generate_end1,archive_search_generate_end2,paired_matches);
//      archive_search_group_add_paired_search(search_group,archive_search_generate_end1,archive_search_generate_end2);
//      archive_search_generate_end1 = NULL;
//      PROFILE_STOP(GP_MAPPER_CUDA_THREAD_RESTART_UNFIT);
//    }
//    PROFILE_STOP(GP_MAPPER_CUDA_THREAD_SELECTING);
//  }
//  // Clean up & Quit
//  ticker_update_mutex(mapper_search->ticker,reads_processed); // Update processed
//  buffered_input_file_close(mapper_search->buffered_fasta_input_end1);
//  if (parameters->io.separated_input_files) buffered_input_file_close(mapper_search->buffered_fasta_input_end2);
//  buffered_output_file_close(buffered_output_file);
//  paired_matches_delete(paired_matches);
//  PROFILE_STOP(GP_MAPPER_CUDA_THREAD);
  pthread_exit(0);
}
