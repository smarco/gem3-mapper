/*
 * PROJECT: GEMMapper
 * FILE: archive_builder.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive/archive_builder.h"
#include "gpu/gpu_structures.h"

/*
 * Constructor
 */
uint64_t sa_sort_length_cmp_values[] = {0,1,5,10,100,1000,10000};
#define sa_sort_length_cmp_num_ranges 6
archive_builder_t* archive_builder_new(
    fm_t* const output_file,
    char* const output_file_name_prefix,
    const archive_type type,
    const uint64_t ns_threshold,
    const sampling_rate_t sa_sampling_rate,
    const sampling_rate_t text_sampling_rate,
    const uint64_t num_threads,
    const uint64_t max_memory) {
  // Allocate
  archive_builder_t* const archive_builder = mm_alloc(archive_builder_t);
  /*
   * Meta-information
   */
  archive_builder->type = type;
  archive_builder->ns_threshold = ns_threshold;
  archive_builder->sa_sampling_rate = sa_sampling_rate;
  archive_builder->text_sampling_rate = text_sampling_rate;
  /*
   * Misc
   */
  // Build Parameters
  archive_builder->num_threads = num_threads;
  archive_builder->max_memory = max_memory;
  /*
   * Archive Components
   */
  // MFASTA Input parsing
  input_multifasta_state_clear(&(archive_builder->parsing_state));
  // Locator
  archive_builder->mm_slab_8MB = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_32M,MM_UNLIMITED_MEM,"archive_builder.8MB");
  archive_builder->mm_slab_32MB = mm_slab_new_(BUFFER_SIZE_32M,BUFFER_SIZE_256M,MM_UNLIMITED_MEM,"archive_builder.32MB");
  archive_builder->locator = locator_builder_new(archive_builder->mm_slab_8MB);
  // Text
  archive_builder->character_occurrences = mm_calloc(DNA_EXT_RANGE*DNA_EXT_RANGE,uint64_t,true);
  archive_builder->enc_rl_text = NULL; // RL-Text is optional
  archive_builder->sampled_rl = NULL; // RL-Text is optional
  // Output
  archive_builder->output_file_manager = output_file;
  archive_builder->output_file_name_prefix = output_file_name_prefix;
  // Return
  return archive_builder;
}
void archive_builder_delete(archive_builder_t* const archive_builder) {
  // Close FM
  fm_close(archive_builder->output_file_manager);
  // Free Text(s)
  if (archive_builder->enc_text!=NULL) {
    dna_text_delete(archive_builder->enc_text);
  }
  if (archive_builder->enc_rl_text!=NULL) {
    dna_text_delete(archive_builder->enc_rl_text);
  }
  dna_text_delete(archive_builder->enc_bwt);
  // Free Occ
  mm_free(archive_builder->character_occurrences);
  // Free Sampled-SA
  sampled_sa_builder_delete(archive_builder->sampled_sa);
  // Free mm-slab
  mm_slab_delete(archive_builder->mm_slab_8MB);
  mm_slab_delete(archive_builder->mm_slab_32MB);
  // Free handler
  mm_free(archive_builder);
}
/*
 * Writers
 */
void archive_builder_write_header(archive_builder_t* const archive_builder) {
  // Write Header
  fm_write_uint64(archive_builder->output_file_manager,ARCHIVE_MODEL_NO);
  fm_write_uint64(archive_builder->output_file_manager,archive_builder->type);
  fm_write_uint64(archive_builder->output_file_manager,1ull); // TODO Remove next index upgrade
  fm_write_uint64(archive_builder->output_file_manager,archive_builder->ns_threshold);
  fm_write_uint64(archive_builder->output_file_manager,0ull); // TODO Remove next index upgrade
}
void archive_builder_write_locator(archive_builder_t* const archive_builder) {
  // Write Locator
  locator_builder_write(archive_builder->output_file_manager,archive_builder->locator);
}
void archive_builder_write_index(
    archive_builder_t* const archive_builder,
    const bool gpu_index,
    const bool check_index,
    const bool verbose) {
  // Select proper text
  dna_text_t* const enc_text = (archive_builder->enc_rl_text==NULL) ?
      archive_builder->enc_text : archive_builder->enc_rl_text;
  // Write Text
  archive_text_write(archive_builder->output_file_manager,
      enc_text,false,archive_builder->forward_text_length,
      archive_builder->sampled_rl,verbose);
  if (archive_builder->sampled_rl!=NULL) sampled_rl_delete(archive_builder->sampled_rl); // Free
  // Create & write the FM-index
  bwt_builder_t* const bwt_builder = fm_index_write(
      archive_builder->output_file_manager,archive_builder->enc_bwt,
      archive_builder->character_occurrences,archive_builder->sampled_sa,check_index,verbose);
  // Create & write the GPU FM-Index
  if (gpu_index) {
    sampled_sa_builder_t* const sampled_sa = archive_builder->sampled_sa;
    const uint32_t sa_sampling_rate = sampled_sa_builder_get_sa_sampling_rate(sampled_sa);
    gpu_structures_write(
        archive_builder->output_file_name_prefix,enc_text,
        archive_builder->forward_text_length,bwt_builder,
        sampled_sa->sa_raw_samples,sa_sampling_rate);
  }
  // Free
  bwt_builder_delete(bwt_builder);
}

