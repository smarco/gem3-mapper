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
 *   Mapper module encapsulates and provides accessors to all
 *   the parameters used by the mapper
 */

#ifndef MAPPER_PARAMETERS_H_
#define MAPPER_PARAMETERS_H_

#include "utils/essentials.h"
#include "io/output_sam.h"
#include "io/output_map.h"
#include "io/input_text.h"
#include "io/buffered_input_file.h"
#include "io/buffered_output_file.h"
#include "archive/search/archive_search_handlers.h"
#include "gpu/gpu_buffer_collection.h"
#include "stats/report_stats_mstats.h"

/* Mapper Mode */
typedef enum { mapper_se, mapper_pe } mapper_type;

/* I/O Parameters */
typedef struct {
  /* Input */
  char* index_file_name;
  bool separated_input_files;
  char* input_file_name;
  char* input_file_name_end1;
  char* input_file_name_end2;
  fm_type input_compression;
  uint64_t input_num_blocks;
  uint64_t input_block_size;
  uint64_t input_buffer_size;
  bool fastq_strictly_normalized;
  /* Output */
  char* output_file_name;
  char* report_file_name;
  file_format_t output_format;
  output_sam_parameters_t sam_parameters;
  output_map_parameters_t map_parameters;
  fm_type output_compression;
  uint64_t output_buffer_size;
  uint64_t output_num_buffers;
  uint64_t mapper_ticker_step;
} mapper_parameters_io_t;

/* System */
typedef struct {
  uint64_t num_threads;
  uint64_t max_memory;
  char* tmp_folder;
} mapper_parameters_system_t;

/* CUDA settings */
typedef struct {
  /* CUDA */
  bool gpu_enabled;
  uint64_t gpu_devices;                  // Bitmask with GPU devices active
  /* I/O */
  uint64_t input_block_size;
  uint64_t input_buffer_size;
  uint64_t output_buffer_size;
  uint64_t output_num_buffers;
  /* GPU Buffering */
  uint64_t gpu_buffer_size;              // Size of each GPU-buffer
  uint64_t num_fmi_bsearch_buffers;      // Number of FMI-BSearch buffers per thread
  uint64_t num_fmi_decode_buffers;       // Number of FMI-Decode buffers per thread
  uint64_t num_kmer_filter_buffers;      // Number of Kmer-filter buffers per thread
  uint64_t num_bpm_distance_buffers;     // Number of BPM-Distance buffers per thread
  uint64_t num_bpm_align_buffers;        // Number of BPM-Align buffers per thread
  /* Stages Configuration */
  bool cpu_emulation;
} mapper_parameters_cuda_t;

/*
 * Misc
 */
typedef struct {
  /* Debug */
  /* QC */
  bool quality_control;
  bool profile;
  profile_reduce_type profile_reduce_type;
  /* Verbose */
  bool verbose_user;
  bool verbose_dev;
} mapper_parameters_misc_t;

/*
 * Mapper Parameters
 */
typedef struct {
  /* CMD line */
  int argc;
  char** argv;                                    // Arguments String
  char* gem_version;                              // GEM version
  /* GEM Structures */
  archive_t* archive;                             // GEM Archive
  gpu_buffer_collection_t* gpu_buffer_collection; // GEM-GPU Index
  input_file_sliced_t* input_file;                // Single file input
  input_file_sliced_t* input_file_end1;           // Split file input (end/1)
  input_file_sliced_t* input_file_end2;           // Split file input (end/2)
  pthread_mutex_t input_file_mutex;               // Input-mutex
  FILE* output_stream;                            // Output Stream
  output_file_t* output_file;                     // Output Handler
  mapping_stats_t* global_mapping_stats;          // Stats Report
  pthread_mutex_t error_report_mutex;             // Mutex for error reporting
  /* Mapper Type */
  mapper_type mapper_type;
  /* I/O Parameters */
  mapper_parameters_io_t io;
  /* Search Parameters */
  search_parameters_t search_parameters;
  /* System */
  mapper_parameters_system_t system;
  /* CUDA settings */
  mapper_parameters_cuda_t cuda;
  /* Miscellaneous */
  mapper_parameters_misc_t misc;
  /* Profile */
  gem_timer_t mapper_time;
  gem_timer_t loading_time;
  /* Parsing Arguments Flags */
  bool min_reported_strata_set;
  bool max_reported_matches_set;
  char *bs_suffix1;
  char *bs_suffix2;
} mapper_parameters_t;

/*
 * Mapper Parameters
 */
void mapper_parameters_set_defaults(mapper_parameters_t* const mapper_parameters);
void mapper_parameters_print(FILE* const stream,mapper_parameters_t* const parameters);

#endif /* MAPPER_PARAMETERS_H_ */
