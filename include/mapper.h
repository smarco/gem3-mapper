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
#include "approximate_search_parameters.h"

// Qualities
#include "quality_model.h"

/*
 * Constants
 */
#define MAPPER_TICKER_STEP 1000

/*
 * Mapper Parameters
 */
typedef enum { mapper_se, mapper_pe, mapper_se_cuda, mapper_graph } mapper_type;
typedef struct {
  /* CMD line */
  int argc;
  char** argv;
  /* Mapper Type */
  mapper_type mapper_type;
  /* I/O Parameters */
  char *index_file_name;
  bool check_index;
  char *input_file_name;
  fm_type input_compression;
  char *output_file_name;
  fm_type output_compression;
  /* I/O Structures */
  archive_t* archive;
  input_file_t* input_file;
  FILE* output_stream;
  output_file_t* output_file;
  uint64_t max_output_buffers;
  file_format_t output_format;
  output_sam_parameters_t sam_parameters;
  /* I/O Attributes (qualities, ...) */
  bool fastq_strictly_normalized;
  bool fastq_try_recovery;
  /* Search Parameters */
  search_parameters_t search_parameters;
  /* Select Parameters */
  select_parameters_t select_parameters;
  /* System */
  uint64_t num_threads;
  uint64_t max_memory;
  char* tmp_folder;
  /* Miscellaneous */
  bool stats;
  bool verbose_user;
  bool verbose_dev;
  /* Extras */
} mapper_parameters_t;

/*
 * Mapper Parameters
 */
GEM_INLINE void mapper_parameters_set_defaults(mapper_parameters_t* const mapper_parameters);
GEM_INLINE void mapper_parameters_print(
    FILE* const stream,mapper_parameters_t* const parameters,const bool dump_index_info);

/*
 * I/O
 */
GEM_INLINE void mapper_SE_output_matches(
    const mapper_parameters_t* const parameters,
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read,matches_t* const matches);

/*
 * SE Mapper
 */
GEM_INLINE void mapper_SE_run(mapper_parameters_t* const mapper_parameters);

/*
 * PE Mapper
 */
GEM_INLINE void mapper_PE_run(const mapper_parameters_t* const mapper_parameters);

/*
 * Error Msg/Functions
 */
void mapper_error_report(FILE* stream); // Mapper error-report function

#endif /* MAPPER_H_ */
