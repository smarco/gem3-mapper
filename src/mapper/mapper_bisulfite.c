/*
 * PROJECT: GEMMapper
 * FILE: mapper_bisulfite.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper/mapper_bisulfite.h"
#include "archive/archive_search_se.h"
#include "archive/archive_search_pe.h"

#define SEQUENCE_INITIAL_LENGTH 200

string_t* mapper_bisulfite_process(string_t* orig,const archive_search_t* end,const char* BS_table) {
  string_t *seq = &archive_search_get_sequence(end)->read;
  string_copy(orig,seq);
  uint64_t len = string_get_length(seq);
  char* seq_buffer = string_get_buffer(seq);
  int64_t pos;
  for (pos=0;pos<len;pos++) {
    seq_buffer[pos]=BS_table[(int)seq_buffer[pos]];
  }
  return seq;
}
/*
 * SE Bisulfite Mapper Thread
 */
void* mapper_SE_bisulfite_thread(mapper_search_t* const mapper_search) {
  // GEM-thread error handler
  gem_thread_register_id(mapper_search->thread_id+1);

  // Create new buffered reader/writer
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  mapper_search->buffered_fasta_input =
      buffered_input_file_new(parameters->input_file,parameters->io.input_buffer_size);
  buffered_output_file_t* const buffered_output_file = buffered_output_file_new(parameters->output_file);
  buffered_input_file_attach_buffered_output(mapper_search->buffered_fasta_input,buffered_output_file);

  // Create an Archive-Search
  mm_search_t* const mm_search = mm_search_new(mm_pool_get_slab(mm_pool_32MB));
  search_parameters_t* const search_parameters = &parameters->base_search_parameters;
  archive_search_se_new(parameters->archive,search_parameters,
      false,NULL,&mapper_search->archive_search);
  archive_search_se_inject_mm(mapper_search->archive_search,mm_search);
  matches_t* const matches = matches_new();
  matches_configure(matches,mapper_search->archive_search->text_collection);

  // Temporary storage for original reads before bisulfite conversion
  uint64_t reads_processed = 0;
  string_t orig_end;
  string_init(&orig_end,SEQUENCE_INITIAL_LENGTH);
  bisulfite_read_t bs_read_mode = parameters->base_search_parameters.bisulfite_read;
  sequence_end_t read_end = (bs_read_mode==bisulfite_read_2) ? paired_end2:paired_end1;

  while (mapper_SE_read_single_sequence(mapper_search)) {
    PROF_INC_COUNTER(GP_MAPPER_NUM_READS);
    // Fully convert reads before searching into archive, making a copy of the original
    if (bs_read_mode==bisulfite_read_inferred) {
      read_end = sequence_get_end_info(&mapper_search->archive_search->sequence);
    }
    string_t* seq_end = NULL;
    switch(read_end) {
    case paired_end1:
      seq_end = mapper_bisulfite_process(&orig_end,mapper_search->archive_search,dna_bisulfite_C2T_table);
      break;
    case paired_end2:
      seq_end = mapper_bisulfite_process(&orig_end,mapper_search->archive_search,dna_bisulfite_G2A_table);
      break;
    default:
      break;
   }
   if (bs_read_mode==bisulfite_read_interleaved) {
     read_end = (read_end==paired_end1) ? paired_end2 : paired_end1;
   }
   // Search into the archive
   archive_search_se(mapper_search->archive_search,matches); // Search matches

   // Copy back original read
   if (seq_end != NULL) string_copy(seq_end,&orig_end);

   // Output matches
   mapper_SE_output_matches(parameters,buffered_output_file,
       mapper_search->archive_search,matches,mapper_search->mapping_stats);

   // Update processed
   if (++reads_processed == MAPPER_TICKER_STEP) {
    ticker_update_mutex(mapper_search->ticker,reads_processed);
    reads_processed=0;
   }

   // Clear
   mm_search_clear(mm_search);
   matches_clear(matches);
  }

  // Update processed
  ticker_update_mutex(mapper_search->ticker,reads_processed);

  // Clean up
  string_destroy(&orig_end); // Free up bisulfite temporary storage
  buffered_input_file_close(mapper_search->buffered_fasta_input);
  buffered_output_file_close(buffered_output_file);
  archive_search_delete(mapper_search->archive_search);
  matches_delete(matches);
  mm_search_delete(mm_search);
  pthread_exit(0);
}

/*
 * PE Bisulfite Mapper Thread
 */
void* mapper_PE_bisulfite_thread(mapper_search_t* const mapper_search) {
  // GEM-thread error handler
  gem_thread_register_id(mapper_search->thread_id+1);

  // Create new buffered reader/writer
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  mapper_PE_prepare_io_buffers(
      parameters,parameters->io.input_buffer_size,&mapper_search->buffered_fasta_input_end1,
      &mapper_search->buffered_fasta_input_end2,&mapper_search->buffered_output_file);

  // Create an Archive-Search
  mm_search_t* const mm_search = mm_search_new(mm_pool_get_slab(mm_pool_32MB));
  search_parameters_t* const search_parameters = &parameters->base_search_parameters;
  archive_search_pe_new(parameters->archive,search_parameters,false,NULL,
      &mapper_search->archive_search_end1,&mapper_search->archive_search_end2);
  archive_search_pe_inject_mm(mapper_search->archive_search_end1,mapper_search->archive_search_end2,mm_search);
  mapper_search->paired_matches = paired_matches_new(mm_search->text_collection);
  paired_matches_configure(mapper_search->paired_matches,mapper_search->archive_search_end1->text_collection);

  // Temporary storage for original reads before bisulfite conversion
  uint64_t reads_processed = 0;
  string_t orig_end1,orig_end2;
  string_init(&orig_end1,SEQUENCE_INITIAL_LENGTH);
  string_init(&orig_end2,SEQUENCE_INITIAL_LENGTH);

  archive_search_t* const archive_search_end1 = mapper_search->archive_search_end1;
  archive_search_t* const archive_search_end2 = mapper_search->archive_search_end2;
  paired_matches_t* const paired_matches = mapper_search->paired_matches;
  while (mapper_PE_read_paired_sequences(mapper_search)) {
    PROF_INC_COUNTER(GP_MAPPER_NUM_READS);

    // Fully convert reads before searching into archive, making a copy of the original
    string_t* seq_end1 = mapper_bisulfite_process(&orig_end1,archive_search_end1,dna_bisulfite_C2T_table);
    string_t* seq_end2 = mapper_bisulfite_process(&orig_end2,archive_search_end2,dna_bisulfite_G2A_table);

    // Search into the archive
    archive_search_pe(archive_search_end1,archive_search_end2,paired_matches);

    // Copy back original read
    string_copy(seq_end1,&orig_end1);
    string_copy(seq_end2,&orig_end2);

    // Output matches
    mapper_PE_output_matches(parameters,mapper_search->buffered_output_file,
        archive_search_end1,archive_search_end2,paired_matches,mapper_search->mapping_stats);

    // Update processed
    if (++reads_processed == MAPPER_TICKER_STEP) {
      ticker_update_mutex(mapper_search->ticker,reads_processed);
      reads_processed=0;
    }

    // Clear
    mm_search_clear(mm_search);
    paired_matches_clear(paired_matches,true);
  }

  // Update processed
  ticker_update_mutex(mapper_search->ticker,reads_processed);

  // Clean up
  string_destroy(&orig_end1);
  string_destroy(&orig_end2);
  buffered_input_file_close(mapper_search->buffered_fasta_input_end1);
  if (parameters->io.separated_input_files) buffered_input_file_close(mapper_search->buffered_fasta_input_end2);
  buffered_output_file_close(mapper_search->buffered_output_file);
  archive_search_delete(mapper_search->archive_search_end1);
  archive_search_delete(mapper_search->archive_search_end2);
  paired_matches_delete(mapper_search->paired_matches);
  mm_search_delete(mm_search);
  pthread_exit(0);
}
