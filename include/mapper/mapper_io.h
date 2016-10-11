/*
 * PROJECT: GEMMapper
 * FILE: mapper.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MAPPER_IO_H_
#define MAPPER_IO_H_

#include "utils/essentials.h"
#include "mapper/mapper_parameters.h"

/*
 * Index loader
 */
void mapper_load_index(mapper_parameters_t* const parameters);

/*
 * Input (Low-level)
 */
void mapper_se_prepare_io_buffers(
    const mapper_parameters_t* const parameters,
    const uint64_t input_buffer_lines,
    buffered_input_file_t** const buffered_fasta_input,
    buffered_output_file_t** const buffered_output_file);
void mapper_pe_prepare_io_buffers(
    const mapper_parameters_t* const parameters,
    const uint64_t input_buffer_lines,
    buffered_input_file_t** const buffered_fasta_input_end1,
    buffered_input_file_t** const buffered_fasta_input_end2,
    buffered_output_file_t** const buffered_output_file);
uint64_t mapper_pe_reload_buffers(
    mapper_parameters_t* const parameters,
    buffered_input_file_t* const buffered_fasta_input_end1,
    buffered_input_file_t* const buffered_fasta_input_end2,
    mapper_stats_t* const mapper_stats);
error_code_t mapper_pe_parse_paired_sequences(
    const mapper_parameters_t* const parameters,
    buffered_input_file_t* const buffered_fasta_input_end1,
    buffered_input_file_t* const buffered_fasta_input_end2,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);

/*
 * Input
 */
error_code_t mapper_se_read_single_sequence(
    archive_search_t* const archive_search,
    buffered_input_file_t* const buffered_fasta_input);
error_code_t mapper_pe_read_paired_sequences(
    mapper_parameters_t* const mapper_parameters,
    buffered_input_file_t* const buffered_fasta_input_end1,
    buffered_input_file_t* const buffered_fasta_input_end2,
    archive_search_handlers_t* const archive_search_handlers,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);

/*
 * Output
 */
void mapper_se_output_matches(
    mapper_parameters_t* const parameters,
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search,
    matches_t* const matches,mapping_stats_t* mstats);
void mapper_pe_output_matches(
    mapper_parameters_t* const parameters,
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,
    mapping_stats_t* mstats);

#endif /* MAPPER_IO_H_ */
