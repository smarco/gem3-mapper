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
#include "essentials.h"

// I/O
#include "input_file.h"
#include "input_parser.h"
#include "input_fasta_parser.h"
#include "output_map.h"
#include "output_sam.h"

// GEM Index
#include "archive.h"
#include "archive_search.h"
#include "archive_select.h"

// ASM
#include "search_parameters.h"

// Qualities
#include "quality_model.h"

/*
 * Constants
 */
#define MAPPER_TICKER_STEP 100000

/*
 * Mapper Parameters
 */
/* Mapper Mode */
typedef enum { mapper_se, mapper_pe, mapper_graph } mapper_type;
/* I/O Parameters */
typedef struct {
  /* Input */
  char* index_file_name;
  bool check_index;
  bool separated_input_files;
  char* input_file_name;
  char* input_file_name_end1;
  char* input_file_name_end2;
  fm_type input_compression;
  uint64_t input_block_size;
  uint64_t input_num_buffers;
  uint64_t input_buffer_lines;
  /* I/O Attributes (qualities, ...) */
  bool fastq_strictly_normalized;
  bool fastq_try_recovery;
  /* Output */
  char* output_file_name;
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
  /* I/O */
  uint64_t input_block_size;
  uint64_t input_num_buffers;
  uint64_t input_buffer_lines;
  uint64_t output_buffer_size;
  uint64_t output_num_buffers;
  /* BPM Buffers */
  uint64_t num_search_groups_per_thread; // Total number of search-groups deployed
  uint64_t bpm_buffer_size;              // Size of each BPM-buffer
} mapper_parameters_cuda_t;
/* Hints */
typedef struct {
  uint64_t avg_read_length;        // Hint on the average read length
  uint64_t std_read_length;        // Hint on the standard deviation of the read length
  uint64_t candidates_per_query;   // Hint on the number of candidates per query
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
  archive_t* archive;
  input_file_t* input_file;
  input_file_t* input_file_end1;
  input_file_t* input_file_end2;
  pthread_mutex_t input_file_mutex;
  FILE* output_stream;
  output_file_t* output_file;
  /* Mapper Type */
  mapper_type mapper_type;
  /* I/O Parameters */
  mapper_parameters_io_t io;
  /* Search Parameters */
  search_parameters_t search_parameters;
  /* Select Parameters */
  select_parameters_t select_parameters;
  /* System */
  mapper_parameters_system_t system;
  /* CUDA settings */
  mapper_parameters_cuda_t cuda;
  /* Hints*/
  mapper_parameters_hints_t hints;
  /* Miscellaneous */
  mapper_parameters_misc_t misc;
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
  /* Mapper parameters */
  mapper_parameters_t* mapper_parameters;
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
    FILE* stream,buffered_input_file_t* const buffered_fasta_input,
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
 * Input
 */
void mapper_PE_prepare_io_buffers(
    const mapper_parameters_t* const parameters,const uint64_t input_buffer_lines,
    buffered_input_file_t** const buffered_fasta_input_end1,
    buffered_input_file_t** const buffered_fasta_input_end2,
    buffered_output_file_t* const buffered_output_file);
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
error_code_t mapper_SE_read_single_sequence(mapper_search_t* const mapper_search);
error_code_t mapper_PE_read_paired_sequences(mapper_search_t* const mapper_search);

/*
 * Output
 */
void mapper_SE_output_matches(
    mapper_parameters_t* const parameters,buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search,matches_t* const matches);
void mapper_PE_output_matches(
    mapper_parameters_t* const parameters,buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

/*
 * SE Mapper
 */
void mapper_SE_run(mapper_parameters_t* const mapper_parameters);

/*
 * PE Mapper
 */
void mapper_PE_run(mapper_parameters_t* const mapper_parameters);

/*
 * Error Messages
 */
#define GEM_ERROR_MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_PAIR "Parsing Input Files. Files '%s,%s' doesn't contain the same number of reads (cannot pair)"
#define GEM_ERROR_MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_SINGLE "Parsing Input File. File '%s' doesn't contain a pair number of reads (cannot pair)"
#define GEM_ERROR_MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_EOF "Parsing Input Files. File '%s' could not read second end (unexpected end-of-file)"
#define GEM_ERROR_MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_NOT_EOF "Parsing Input Files. File '%s' has too many reads (expected end-of-file)"

/*
 * Error Macros
 */
#define MAPPER_ERROR_PE_PARSE_UNSYNCH_INPUT_FILES(mapper_parameters) \
  if (mapper_parameters->io.separated_input_files) { \
    gem_fatal_error(MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_PAIR, \
        input_file_get_file_name(mapper_parameters->input_file_end1), \
        input_file_get_file_name(mapper_parameters->input_file_end2)); \
  } else { \
    gem_fatal_error(MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_SINGLE, \
        input_file_get_file_name(mapper_parameters->input_file)); \
  }


#endif /* MAPPER_H_ */
