/*
 * PROJECT: GEMMapper
 * FILE: archive_builder.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef ARCHIVE_BUILDER_H_
#define ARCHIVE_BUILDER_H_

#include "utils/essentials.h"
#include "io/input_file.h"
#include "io/input_multifasta_parser.h"
#include "archive/archive.h"
#include "archive/locator_builder.h"
#include "fm_index/sa_builder.h"

/*
 * Debug
 */
#define ARCHIVE_BUILDER_DEBUG_SUFFIX_LENGTH   100
#define ARCHIVE_BUILDER_DEBUG_RECORD_STATS false

/*
 * Archive Builder (Indexer)
 */
typedef enum { index_complement_no=0, index_complement_yes=1, index_complement_auto=2 } indexed_complement_t;
typedef struct {
  /*
   * Meta-information
   */
  archive_type type;                       // Archive type
  indexed_complement_t indexed_complement; // Forces the storage of the RC
  uint64_t complement_size_threshold;      // Maximum text size allowed to store the RC
  uint64_t ns_threshold;                   // Minimum length of a stretch of Ns to be removed
  bool indexed_reverse_text;               // Indexed reverse text (backwards text)
  /*
   * Archive Components
   */
  /* MFASTA Input parsing */
  input_multifasta_state_t parsing_state;   // Text-Building state
  /* SA-Builder */
  sa_builder_t* sa_builder;                 // SA-Builder
  /* Locator */
  locator_builder_t* locator;               // Sequence locator (from MultiFASTA)
  /* Text */
  uint64_t forward_text_length;             // Length of the forward text
  dna_text_t* enc_text;                     // Encoded Input Text (from MultiFASTA)
  dna_text_t* enc_rl_text;                  // Encoded Run-Length Compacted Input Text (from MultiFASTA)
  uint64_t* character_occurrences;          // Total occurrences of each character
  sampled_rl_t* sampled_rl;                 // Sampled RL-text (to retrieve approximate source text-positions)
  /* FM-Index */
  dna_text_t* enc_bwt;                      // BWT text
  sampling_rate_t sa_sampling_rate;         // SA Sampling Rate
  sampling_rate_t text_sampling_rate;       // Text Sampling Rate
  sampled_sa_builder_t* sampled_sa;         // Sampled SA Builder
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
archive_builder_t* archive_builder_new(
    fm_t* const output_file,
    char* const output_file_name_prefix,
    const archive_type type,
    const indexed_complement_t indexed_complement,
    const uint64_t complement_size_threshold,
    const uint64_t ns_threshold,
    const sampling_rate_t sa_sampling_rate,
    const sampling_rate_t text_sampling_rate,
    const bool indexed_reverse_text,
    const uint64_t num_threads,
    const uint64_t max_memory);
void archive_builder_delete(archive_builder_t* const archive_builder);

/*
 * Writers
 */
void archive_builder_write_header(archive_builder_t* const archive_builder);
void archive_builder_write_locator(archive_builder_t* const archive_builder);
void archive_builder_write_index(
    archive_builder_t* const archive_builder,
    const bool gpu_index,
    const bool check_index,
    const bool verbose);
void archive_builder_write_index_reverse(
    archive_builder_t* const archive_builder,
    const bool check_index,
    const bool verbose);

#endif /* ARCHIVE_BUILDER_H_ */
