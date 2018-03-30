/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Archive builder main module
 */

#ifndef ARCHIVE_BUILDER_H_
#define ARCHIVE_BUILDER_H_

#include "utils/essentials.h"
#include "io/input_multifasta.h"
#include "archive/archive.h"
#include "archive/locator_builder.h"
#include "fm_index/sa_builder/sa_builder.h"
#include "fm_index/rank_mtable_builder.h"

/*
 * Debug
 */
#define ARCHIVE_BUILDER_DEBUG_SUFFIX_LENGTH   100
#define ARCHIVE_BUILDER_DEBUG_RECORD_STATS    false

/*
 * Archive Builder (Indexer)
 */
typedef struct {
  /* Meta-information */
  archive_type type;                             // Archive type
  uint64_t ns_threshold;                         // Minimum length of a stretch of Ns to be removed
  bool gpu_index;                                // Index generated used GPU compiled GEM
  /* Input Multi-FASTA */
  input_multifasta_file_t input_multifasta_file; // Input multi-FASTA file
  /* Archive */
  sa_builder_t* sa_builder;                      // SA-Builder
  locator_builder_t* locator;                    // Sequence locator (from MultiFASTA)
  /* Text */
  uint64_t forward_text_length;                  // Length of the forward text
  dna_text_t* enc_text;                          // Encoded Input Text (from MultiFASTA)
  dna_text_t* enc_rl_text;                       // Encoded Run-Length Compacted Input Text (from MultiFASTA)
  uint64_t* character_occurrences;               // Total occurrences of each character
  sampled_rl_t* sampled_rl;                      // Sampled RL-text (to retrieve approximate source text-positions)
  /* FM-Index */
  dna_text_t* enc_bwt;                           // BWT text
  sampling_rate_t sa_sampling_rate;              // SA Sampling Rate
  sampling_rate_t text_sampling_rate;            // Text Sampling Rate
  sampled_sa_builder_t* sampled_sa;              // Sampled SA Builder
  /* Output */
  char* output_file_name_prefix;                 // Output Text FileName Prefix
  fm_t* output_file_manager;                     // Output Manager
  /* Misc */
  uint64_t num_threads;                          // Total number threads to split the work across
  uint64_t max_memory;                           // Max. memory to use
  ticker_t ticker;                               // Index Builder Ticker
  FILE* info_file;                               // Index info file
  /* MM */
  mm_slab_t* mm_slab_8MB;                        // MM-Slab
  mm_slab_t* mm_slab_32MB;                       // MM-Slab
} archive_builder_t;

/*
 * Archive Builder
 */
archive_builder_t* archive_builder_new(
    char* const input_multifasta_file_name,
    fm_t* const output_file,
    char* const output_file_name_prefix,
    const archive_type type,
    const bool gpu_index,
    const uint64_t ns_threshold,
    const sampling_rate_t sa_sampling_rate,
    const sampling_rate_t text_sampling_rate,
    const uint64_t num_threads,
    const uint64_t max_memory,
    FILE* const info_file);
void archive_builder_delete(
    archive_builder_t* const archive_builder);

/*
 * Writers
 */
void archive_builder_write_header(
    archive_builder_t* const archive_builder);
void archive_builder_write_locator(
    archive_builder_t* const archive_builder);
void archive_builder_write_index(
    archive_builder_t* const archive_builder,
    const bool gpu_index,
    const bool check_index,
    const bool verbose);

#endif /* ARCHIVE_BUILDER_H_ */
