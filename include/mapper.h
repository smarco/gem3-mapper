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

// ASM
#include "approximate_search_parameters.h"

// Qualities
#include "quality_model.h"

/*
 * Mapper Parameters
 */
typedef enum { mapper_se, mapper_pe, mapper_se_cuda, mapper_graph } mapper_type;
typedef struct {
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
  file_format_t output_format;
  /* I/O Attributes (qualities, ...) */
  bool fastq_strictly_normalized;
  bool fastq_try_recovery;
  quality_format_t quality_format;
  quality_model_t quality_model;
  uint64_t quality_threshold;
  /* Single-end Alignment */
  mapping_mode_t mapping_mode;
  float mapping_degree;
  float max_search_error;
  float max_filtering_error;
  float complete_strata_after_best;
  float min_matching_length;
  uint64_t max_search_matches;
  char* mismatch_alphabet;
  /* Paired-end Alignment */
  /* Reporting */
  uint64_t min_decoded_strata;
  uint64_t max_decoded_matches;
  uint64_t min_reported_matches;
  uint64_t max_reported_matches;
  /* System */
  uint64_t num_threads;
  uint64_t max_memory;
  char* tmp_folder;
  /* Miscellaneous */
  bool user_verbose;
  bool dev_verbose;
  /* Extras */
} mapper_parameters_t;

/*
 * Mapper Parameters
 */

GEM_INLINE void mapper_parameters_set_defaults(mapper_parameters_t* const mapper_parameters);

/*
 * SE Mapper
 */
GEM_INLINE void mapper_SE_configure_archive_search(
    archive_search_t* const archive_search,
    const mapper_parameters_t* const parameters);
GEM_INLINE void mapper_SE_output_matches(
    const mapper_parameters_t* const parameters,
    buffered_output_file_t* const buffered_output_file,
    const sequence_t* const seq_read,matches_t* const matches);
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
