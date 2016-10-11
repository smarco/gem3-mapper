/*
 * PROJECT: GEMMapper
 * FILE: mapper.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MAPPER_PARAMETERS_H_
#define MAPPER_PARAMETERS_H_

#include "utils/essentials.h"
#include "io/output_sam.h"
#include "io/output_map.h"
#include "io/input_file.h"
#include "io/input_parser.h"
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
  bool cuda_enabled;
  /* I/O */
  uint64_t input_block_size;
  uint64_t input_buffer_size;
  uint64_t output_buffer_size;
  uint64_t output_num_buffers;
  /* GPU Buffering */
  uint64_t gpu_buffer_size;              // Size of each GPU-buffer
  uint64_t num_fmi_bsearch_buffers;      // Number of FMI-BSearch buffers per thread
  uint64_t num_fmi_decode_buffers;       // Number of FMI-Decode buffers per thread
  uint64_t num_bpm_buffers;              // Number of BPM buffers per thread
  /* Stages Configuration */
  bool cpu_emulation;
} mapper_parameters_cuda_t;

/* Hints */
//typedef struct {
//  uint64_t avg_read_length;        // Hint on the average read length
//  uint64_t std_read_length;        // Hint on the standard deviation of the read length
//  uint64_t candidates_per_query;   // Hint on the number of candidates per query
//  uint64_t dummy;
//} mapper_parameters_hints_t;

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
  char** argv;
  char* gem_version;
  /* GEM Structures */
  archive_t* archive;                             // GEM Archive
  gpu_buffer_collection_t* gpu_buffer_collection; // GEM-GPU Index
  input_file_sliced_t* input_file;
  input_file_sliced_t* input_file_end1;
  input_file_sliced_t* input_file_end2;
  pthread_mutex_t input_file_mutex;
  FILE* output_stream;
  output_file_t* output_file;
  mapping_stats_t* global_mapping_stats;          // Stats Report
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
  /* Hints*/
  // mapper_parameters_hints_t hints;
  /* Miscellaneous */
  mapper_parameters_misc_t misc;
  /* Profile */
  gem_timer_t mapper_time;
  gem_timer_t loading_time;
} mapper_parameters_t;

/*
 * Mapper Parameters
 */
void mapper_parameters_set_defaults(mapper_parameters_t* const mapper_parameters);
void mapper_parameters_print(FILE* const stream,mapper_parameters_t* const parameters);

#endif /* MAPPER_PARAMETERS_H_ */
