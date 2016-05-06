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
 * Read Parser
 */
int input_fasta_parse_sequence(
    buffered_input_file_t* const buffered_fasta_input,
    sequence_t* const sequence,
    const bool strictly_normalized,
    const bool check_input_buffer);

/*
 * Display
 */
void input_fasta_parser_prompt_error(
    buffered_input_file_t* const buffered_input,
    const error_code_t error_code);

#endif /* INPUT_FASTA_PARSER_H_ */
