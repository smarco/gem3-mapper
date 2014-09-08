/*
 * PROJECT: GEMMapper
 * FILE: archive_builder.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *
 *   TODO:
 *     Loop-Peeling CUDA-version. Align index-sequences to 8 chars (mod8=0)
 */

#ifndef ARCHIVE_BUILDER_H_
#define ARCHIVE_BUILDER_H_

#include "essentials.h"

#include "archive.h"
#include "input_file.h"
#include "input_multifasta_parser.h"
#include "sa_builder.h"
#include "graph_text_builder.h"

/*
 * Debug
 */
#define ARCHIVE_BUILDER_DEBUG_SUFFIX_LENGTH   100
#define ARCHIVE_BUILDER_DEBUG_RECORD_STATS false

/*
 * Checkers
 */
#define ARCHIVE_BUILDER_CHECK(archive_builder) GEM_CHECK_NULL(archive_builder)

/*
 * Archive Builder (Indexer)
 */
typedef enum { index_complement_no=0, index_complement_yes=1, index_complement_auto=2 } indexed_complement_t;
typedef struct {
  /*
   * Meta-information
   */
  index_t index_type;                      // Index/Archive type (architecture)
  filter_t filter_type;                    // Filter applied to the original text (MFasta)
  indexed_complement_t indexed_complement; // Forces the storage of the RC
  uint64_t complement_size_threshold;      // Maximum text size allowed to store the RC
  uint64_t ns_threshold;                   // Minimum length of a stretch of Ns to be removed
  /*
   * Archive Components
   */
  /* MFASTA Input parsing */
  input_multifasta_state_t parsing_state;   // Text-Building state
  /* SA-Builder */
  sa_builder_t* sa_builder;                 // SA-Builder
  /* Graph Components */
  graph_text_builder_t* graph;              // Graph structure
  /* Locator */
  locator_builder_t* locator;               // Sequence locator (from MultiFASTA)
  /* Text */
  dna_text_t* enc_text;                     // Encoded Input Text (from MultiFASTA)
  dna_text_t* enc_rl_text;                  // Encoded Run-Length Compacted Input Text (from MultiFASTA)
  uint64_t* character_occurrences;          // Total occurrences of each character
  uint64_t* sampled_rl_text;                // Sampled RL-text (to retrieve approximate source text-positions)
  /* FM-Index */
  dna_text_t* enc_bwt;                      // BWT text
  sampling_rate_t sampling_rate;            // Sampling Rate
  sampled_sa_builder_t* sampled_sa;         // Sampled SA
  /* Output */
  char* output_file_name_prefix;            // Output Text FileName Prefix
  fm_t* output_file_manager;                // Output Manager
  /*
   * Misc
   */
  /* Build Parameters */
  uint64_t num_threads;                     // Total number threads to split the work across
  uint64_t max_memory;
  /* Misc */
  ticker_t ticker;                          // Index Builder Ticker
} archive_builder_t;

/*
 * Archive Builder
 */
GEM_INLINE archive_builder_t* archive_builder_new(
    fm_t* const output_file,char* const output_file_name_prefix,
    const index_t index_type,const filter_t filter_type,
    const indexed_complement_t indexed_complement,const uint64_t complement_size_threshold,
    const uint64_t ns_threshold,const sampling_rate_t sampling_rate,
    const uint64_t num_threads,const uint64_t max_memory);
GEM_INLINE void archive_builder_delete(archive_builder_t* const archive_builder);

/*
 * Archive Write
 *   // 1. archive_builder_write_header()
 *   // 2. archive_builder_write_locator()
 *   // 3. archive_builder_build_index()
 */
GEM_INLINE void archive_builder_write_header(archive_builder_t* const archive_builder);
GEM_INLINE void archive_builder_write_locator(archive_builder_t* const archive_builder);
GEM_INLINE void archive_builder_write_graph__jump_table(archive_builder_t* const archive_builder,const bool display_links);

/*
 * STEP1 Archive Build :: Process MultiFASTA,GRAPH file(s)
 */
GEM_INLINE void archive_builder_process_graph(
    archive_builder_t* const archive_builder,input_file_t* const input_graph,
    const bool dump_graph_links,const bool verbose);
GEM_INLINE void archive_builder_process_multifasta(
    archive_builder_t* const archive_builder,input_file_t* const input_multifasta,
    const bool dump_locator_intervals,const bool dump_indexed_text,const bool verbose);
GEM_INLINE void archive_builder_process_run_length_text(
    archive_builder_t* const archive_builder,const bool dump_run_length_text,const bool verbose);
GEM_INLINE void archive_builder_process_multifasta__graph(
    archive_builder_t* const archive_builder,input_file_t* const input_multifasta,
    const bool dump_locator_intervals,const bool dump_indexed_text,const bool dump_graph_links,const bool verbose);

/*
 * STEP2 Archive Build :: Build BWT (SA)
 */
GEM_INLINE void archive_builder_build_bwt(
    archive_builder_t* const archive_builder,
    const bool dump_explicit_sa,const bool dump_bwt,const bool verbose);

/*
 * STEP3 Archive Build :: Create Index (FM-Index)
 */
GEM_INLINE void archive_builder_build_index(
    archive_builder_t* const archive_builder,const bool check_index,const bool verbose);

#endif /* ARCHIVE_BUILDER_H_ */
