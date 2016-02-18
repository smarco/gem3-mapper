/*
 * PROJECT: GEMMapper
 * FILE: input_fasta_parser.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Input parser for FASTA/FASTQ/MULTIFASTA files
 */

#ifndef INPUT_FASTA_PARSER_H_
#define INPUT_FASTA_PARSER_H_

#include "system/commons.h"
#include "io/input_file.h"
#include "io/buffered_input_file.h"
#include "data_structures/sequence.h"

/*
 * Input FASTA/FASTQ Constants
 */
#define FASTA_TAG_BEGIN '>'
#define FASTQ_TAG_BEGIN '@'
#define FASTQ_SEP '+'

/*
 * FASTQ File basics
 */
bool input_file_test_fasta(
    input_file_t* const input_file,
    fasta_file_format_t* const fasta_file_format,
    const bool show_errors);
void input_fasta_parser_prompt_error(
    buffered_input_file_t* const buffered_fasta_input,
    uint64_t line_num,
    uint64_t column_pos,
    const error_code_t error_code);
void input_fasta_parser_next_record(
    buffered_input_file_t* const buffered_fasta_input,
    char* const line_start);

/*
 * Accessors
 */
bool input_fasta_is_fasta(input_file_t* const input_file);
bool input_fasta_is_fastq(input_file_t* const input_file);

/*
 * Read Parser
 */
error_code_t input_fasta_parse_sequence(
    buffered_input_file_t* const buffered_fasta_input,
    sequence_t* const seq_read,
    const bool strictly_normalized,
    const bool try_recovery,
    const bool check_input_buffer);

#endif /* INPUT_FASTA_PARSER_H_ */
