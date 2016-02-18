/*
 * PROJECT: GEMMapper
 * FILE: archive_builder_index.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive/archive_builder_index.h"
#include "archive/archive_builder_text.h"
#include "data_structures/dna_text.h"

/*
 * Build BWT (SA)
 */
void archive_builder_index_build_bwt(
    archive_builder_t* const archive_builder,
    const bool gpu_index,
    const bool dump_bwt,
    const bool dump_explicit_sa,
    const bool verbose) {
  // Allocate BWT-text
  const uint64_t text_length = dna_text_get_length(archive_builder->enc_text);
  archive_builder->enc_bwt = dna_text_new(text_length);
  dna_text_set_length(archive_builder->enc_bwt,text_length);
  // SA-Builder (SA-sorting)
  archive_builder->sa_builder = sa_builder_new(
      archive_builder->output_file_name_prefix,archive_builder->enc_text,
      archive_builder->num_threads,archive_builder->max_memory);
  // Count k-mers
  sa_builder_count_suffixes(archive_builder->sa_builder,archive_builder->character_occurrences,verbose);
  // Write SA-positions
  sa_builder_store_suffixes(archive_builder->sa_builder,verbose);
  // Sort suffixes & sample SA
  archive_builder->sampled_sa = sampled_sa_builder_new(
      text_length,archive_builder->num_threads,archive_builder->sa_sampling_rate,
      archive_builder->text_sampling_rate,gpu_index,mm_pool_get_slab(mm_pool_32MB)); // Allocate Sampled-SA
  sa_builder_sort_suffixes(archive_builder->sa_builder,
      archive_builder->enc_bwt,archive_builder->sampled_sa,verbose);
  // DEBUG
  if (dump_bwt) archive_builder_index_print_bwt(archive_builder,".bwt",true);
  // if (dump_explicit_sa) archive_builder_index_print_explicit_sa(archive_builder,".sa"); // No full sort (only hash-sort)
  // Free
  sa_builder_delete(archive_builder->sa_builder); // Delete SA-Builder
}
void archive_builder_index_build_bwt_reverse(
    archive_builder_t* const archive_builder,
    const bool dump_reverse_indexed_text,
    const bool dump_bwt,
    const bool dump_explicit_sa,
    const bool verbose) {
  // Reverse BWT-text
  const uint64_t text_length = dna_text_get_length(archive_builder->enc_text);
  uint8_t* const enc_text_buffer = dna_text_get_text(archive_builder->enc_text);
  const uint64_t text_length_half = text_length/2;
  uint64_t i;
  for (i=0;i<text_length_half;++i) {
    SWAP(enc_text_buffer[i],enc_text_buffer[text_length-i-1]);
  }
  // DEBUG index_text
  if (dump_reverse_indexed_text) archive_builder_text_dump(archive_builder,".rev.text");
  memset(archive_builder->character_occurrences,0,DNA_EXT_RANGE*DNA_EXT_RANGE*UINT64_SIZE);
  // SA-Builder (SA-sorting)
  archive_builder->sa_builder = sa_builder_new(
      archive_builder->output_file_name_prefix,archive_builder->enc_text,
      archive_builder->num_threads,archive_builder->max_memory);
  // Count k-mers
  sa_builder_count_suffixes(archive_builder->sa_builder,archive_builder->character_occurrences,verbose);
  // Write SA-positions
  sa_builder_store_suffixes(archive_builder->sa_builder,verbose);
  // Sort suffixes & sample SA
  sa_builder_sort_suffixes(archive_builder->sa_builder,archive_builder->enc_bwt,NULL,verbose);
  // DEBUG
  if (dump_bwt) archive_builder_index_print_bwt(archive_builder,".rev.bwt",true);
//  if (dump_explicit_sa) archive_builder_index_print_explicit_sa(archive_builder,".rev.sa"); // No full sort (only hash-sort)
  // Free
  sa_builder_delete(archive_builder->sa_builder); // Delete SA-Builder
}
/*
 * Display
 */
void archive_builder_index_print_explicit_sa(
    archive_builder_t* const archive_builder,
    const char* const extension) {
  // Open File
  char* const debug_explicit_sa_file_name = gem_strcat(archive_builder->output_file_name_prefix,extension);
  FILE* const debug_explicit_sa_file = fopen(debug_explicit_sa_file_name,"w");
  // Traverse SA & Print Explicit SA => (SApos,SA[SApos...SApos+SAFixLength])
  sa_builder_t* const sa_builder = archive_builder->sa_builder;
  const uint64_t sa_length = dna_text_get_length(sa_builder->enc_text);
  fm_t* const sa_positions_file = fm_open_file(sa_builder->sa_positions_file_name,FM_READ);
  uint64_t i;
  for (i=0;i<sa_length;++i) {
    const uint64_t sa_position = fm_read_uint64(sa_positions_file);
    sa_builder_debug_print_sa(debug_explicit_sa_file,sa_builder,sa_position,100);
  }
  // Free
  fm_close(sa_positions_file);
  fclose(debug_explicit_sa_file);
  mm_free(debug_explicit_sa_file_name);
}
void archive_builder_index_print_bwt(
    archive_builder_t* const archive_builder,
    const char* const extension,
    const bool verbose) {
  ticker_t ticker;
  ticker_percentage_reset(&ticker,verbose,"Building-BWT::Dumping BWT",1,1,true);
  char* const file_name = gem_strcat(archive_builder->output_file_name_prefix,extension);
  FILE* const bwt_file = fopen(file_name,"w"); // Open file
  dna_text_pretty_print_content(bwt_file,archive_builder->enc_bwt,80); // Dump BWT
  fclose(bwt_file);
  mm_free(file_name);
  ticker_finish(&ticker);
}
