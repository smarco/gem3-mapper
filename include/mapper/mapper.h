/*
 * PROJECT: GEMMapper
 * FILE: mapper.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MAPPER_H_
#define MAPPER_H_

// GEM essentials
#include "utils/essentials.h"
#include "system/profiler.h"
#include "io/input_file.h"
#include "io/input_parser.h"
#include "io/input_fasta_parser.h"
#include "io/output_map.h"
#include "io/output_sam.h"
#include "data_structures/quality_model.h"
#include "archive/archive.h"
#include "archive/archive_search.h"
#include "archive/archive_search_parameters.h"
#include "archive/archive_select.h"
#include "gpu/gpu_buffer_collection.h"
#include "stats/report_stats_mstats.h"

/*
 * Constants
 */
#define MAPPER_TICKER_STEP  100000

/*
 * Mapper Parameters
 */
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
  bool cpu_emulated;
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
} mapper_parameters_cuda_t;
/* Hints */
typedef struct {
//  uint64_t avg_read_length;        // Hint on the average read length
//  uint64_t std_read_length;        // Hint on the standard deviation of the read length
//  uint64_t candidates_per_query;   // Hint on the number of candidates per query
  uint64_t dummy;
} mapper_parameters_hints_t;
/* Misc */
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
/* Mapper Parameters */
typedef struct {
  /* CMD line */
  int argc;
  char** argv;
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
  search_parameters_t base_search_parameters;
  /* System */
  mapper_parameters_system_t system;
  /* CUDA settings */
  mapper_parameters_cuda_t cuda;
  /* Hints*/
  mapper_parameters_hints_t hints;
  /* Miscellaneous */
  mapper_parameters_misc_t misc;
  /* Profile */
  gem_timer_t mapper_time;
  gem_timer_t loading_time;
} mapper_parameters_t;
/*
 * Mapper Search
 */
typedef struct {
  /* Thread Info */
  uint64_t thread_id;
  pthread_t* thread_data;
  /* I/O */
  bool paired_end;
  buffered_input_file_t* buffered_fasta_input;
  buffered_input_file_t* buffered_fasta_input_end1;
  buffered_input_file_t* buffered_fasta_input_end2;
  buffered_output_file_t* buffered_output_file;
  /* Mapper parameters */
  mapper_parameters_t* mapper_parameters;
	/* Per thread stats report structures */
	mapping_stats_t* mapping_stats;
  /* Archive-Search */
  archive_search_t* archive_search;
  archive_search_t* archive_search_end1;
  archive_search_t* archive_search_end2;
  paired_matches_t* paired_matches;
  /* Ticker */
  ticker_t* ticker;
} mapper_search_t;

/*
 * Report
 */
void mapper_display_input_state(
    FILE* stream,
    buffered_input_file_t* const buffered_fasta_input,
    const sequence_t* const sequence);

/*
 * Mapper Parameters
 */
void mapper_parameters_set_defaults(mapper_parameters_t* const mapper_parameters);
void mapper_parameters_print(FILE* const stream,mapper_parameters_t* const parameters);

/*
 * Index loader
 */
void mapper_load_index(mapper_parameters_t* const parameters);

/*
 * Input (Low-level)
 */
void mapper_SE_prepare_io_buffers(
    const mapper_parameters_t* const parameters,
    const uint64_t input_buffer_lines,
    buffered_input_file_t** const buffered_fasta_input,
    buffered_output_file_t** const buffered_output_file);
void mapper_PE_prepare_io_buffers(
    const mapper_parameters_t* const parameters,
    const uint64_t input_buffer_lines,
    buffered_input_file_t** const buffered_fasta_input_end1,
    buffered_input_file_t** const buffered_fasta_input_end2,
    buffered_output_file_t** const buffered_output_file);
uint64_t mapper_PE_reload_buffers(
    mapper_parameters_t* const parameters,
    buffered_input_file_t* const buffered_fasta_input_end1,
    buffered_input_file_t* const buffered_fasta_input_end2);
error_code_t mapper_PE_parse_paired_sequences(
    const mapper_parameters_t* const parameters,
    buffered_input_file_t* const buffered_fasta_input_end1,
    buffered_input_file_t* const buffered_fasta_input_end2,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);

/*
 * Input (High-level)
 */
error_code_t mapper_SE_read_single_sequence(mapper_search_t* const mapper_search);
error_code_t mapper_PE_read_paired_sequences(mapper_search_t* const mapper_search);

/*
 * Output
 */
void mapper_SE_output_matches(
    mapper_parameters_t* const parameters,
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search,
    matches_t* const matches,mapping_stats_t* mstats);
void mapper_PE_output_matches(
    mapper_parameters_t* const parameters,
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,
    mapping_stats_t* mstats);

/*
 * SE Mapper
 */
void mapper_SE_run(mapper_parameters_t* const mapper_parameters);

/*
 * PE Mapper
 */
void mapper_PE_run(mapper_parameters_t* const mapper_parameters);

#endif /* MAPPER_H_ */
