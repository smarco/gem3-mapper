/*
 * PROJECT: GEMMapper
 * FILE: archive_builder.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_builder.h"

/*
 * Constructor
 */
uint64_t sa_sort_length_cmp_values[] = {0,1,5,10,100,1000,10000};
#define sa_sort_length_cmp_num_ranges 6
GEM_INLINE archive_builder_t* archive_builder_new(
    fm_t* const output_file,char* const output_file_name_prefix,
    const index_t index_type,const filter_t filter_type,
    const indexed_complement_t indexed_complement,const uint64_t complement_size_threshold,
    const uint64_t ns_threshold,const sampling_rate_t sampling_rate,
    const uint64_t num_threads,const uint64_t max_memory) {
  // Allocate
  archive_builder_t* const archive_builder = mm_alloc(archive_builder_t);
  /*
   * Meta-information
   */
  archive_builder->index_type = index_type;
  archive_builder->filter_type = filter_type;
  archive_builder->indexed_complement = indexed_complement;
  archive_builder->complement_size_threshold = complement_size_threshold;
  archive_builder->ns_threshold = ns_threshold;
  archive_builder->sampling_rate = sampling_rate;
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
  // Graph Components
  archive_builder->graph = NULL; // Graph is optional
  // Locator
  archive_builder->locator = locator_builder_new(mm_pool_get_slab(mm_pool_2MB));
  // Text
  archive_builder->character_occurrences = mm_calloc(DNA_EXT_RANGE*DNA_EXT_RANGE,uint64_t,true);
  // Output
  archive_builder->output_file_manager = output_file;
  archive_builder->output_file_name_prefix = output_file_name_prefix;
  // Return
  return archive_builder;
}
GEM_INLINE void archive_builder_delete(archive_builder_t* const archive_builder) {
  ARCHIVE_BUILDER_CHECK(archive_builder);
  // Archive Components
  mm_free(archive_builder->character_occurrences);
  // Free handler
  mm_free(archive_builder);
}
/*
 * Archive Write
 *   // 1. archive_builder_write_header()
 *   // 2. archive_builder_write_locator()
 *   // 3. archive_builder_write_graph__jump_table()
 *   // 4. archive_builder_build_index()
 */
GEM_INLINE void archive_builder_write_header(archive_builder_t* const archive_builder) {
  ARCHIVE_BUILDER_CHECK(archive_builder);
  // Write Header
  fm_write_uint64(archive_builder->output_file_manager,archive_builder->index_type);
  fm_write_uint64(archive_builder->output_file_manager,archive_builder->filter_type);
  fm_write_uint64(archive_builder->output_file_manager,archive_builder->indexed_complement);
  fm_write_uint64(archive_builder->output_file_manager,archive_builder->ns_threshold);
}
GEM_INLINE void archive_builder_write_locator(archive_builder_t* const archive_builder) {
  // Write Locator
  locator_builder_write(archive_builder->output_file_manager,archive_builder->locator);
}
GEM_INLINE void archive_builder_write_graph__jump_table(archive_builder_t* const archive_builder,const bool display_links) {
  // Write Graph & Jump-Locator (if any)
  if (archive_builder->index_type == fm_dna_graph) {
    // Write Graph
    graph_text_builder_write(archive_builder->output_file_manager,archive_builder->graph,archive_builder->locator);
    graph_text_builder_link_table_print(gem_info_get_stream(),archive_builder->graph,archive_builder->locator,true); // DEBUG: Print sorted links
    graph_text_builder_delete(archive_builder->graph); // Free
  }
}
/*
 * Archive Build STEP4 :: Create Index (FM-Index)
 *   1. Generate archive
 *     1.1 Write IndexText
 *     1.2 Write Sampled-SA
 *     1.3 FM-Index
 *       1.3.1 FM-Index Structure
 *       1.3.2 BWT Structure
 *       1.3.3 Memoization Table (Rank calls)
 */
GEM_INLINE void archive_builder_build_index(
    archive_builder_t* const archive_builder,
    const bool check_index,const bool verbose) {
  /*
   * Write Text & Free
   */
  // RL-Text
  if (archive_builder->index_type == fm_dna_run_length) {
    // Swap back text and RL-text
    SWAP(archive_builder->enc_text,archive_builder->enc_rl_text);
    // Store Sampled-RL
    // TODO
    // Delete RL-Text
    dna_text_delete(archive_builder->enc_rl_text); // Free
  }
  // Text
  dna_text_write(archive_builder->output_file_manager,archive_builder->enc_text);
  if (verbose) dna_text_print(gem_info_get_stream(),archive_builder->enc_text);
  dna_text_delete(archive_builder->enc_text); // Free
  /*
   * Create & write the FM-index
   */
  fm_index_builder(archive_builder->output_file_manager,
      archive_builder->enc_bwt,archive_builder->character_occurrences,archive_builder->sampled_sa,
      check_index,verbose,archive_builder->num_threads);
  fm_close(archive_builder->output_file_manager); // Close FM
}

