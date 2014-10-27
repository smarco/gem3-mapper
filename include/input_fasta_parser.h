/*
 * PROJECT: GEMMapper
 * FILE: input_fasta_parser.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Input parser for FASTA/FASTQ/MULTIFASTA files
 */

#ifndef INPUT_FASTA_PARSER_H_
#define INPUT_FASTA_PARSER_H_

#include "commons.h"

#include "input_file.h"
#include "buffered_input_file.h"
#include "sequence.h"

/*
 * Input FASTA/FASTQ Constants
 */
#define FASTA_TAG_BEGIN '>'
#define FASTQ_TAG_BEGIN '@'
#define FASTQ_SEP '+'

/*
 * FASTQ File basics
 */
GEM_INLINE bool input_file_test_fasta(
    input_file_t* const input_file,fasta_file_format_t* const fasta_file_format,const bool show_errors);
GEM_INLINE void input_fasta_parser_prompt_error(
    buffered_input_file_t* const buffered_fasta_input,
    uint64_t line_num,uint64_t column_pos,const error_code_t error_code);
GEM_INLINE void input_fasta_parser_next_record(
    buffered_input_file_t* const buffered_fasta_input,char* const line_start);

/*
 * Accessors
 */
GEM_INLINE bool input_fasta_is_fasta(input_file_t* const input_file);
GEM_INLINE bool input_fasta_is_fastq(input_file_t* const input_file);

/*
 * Read Parser
 */
GEM_INLINE error_code_t input_fasta_parse_sequence(
    buffered_input_file_t* const buffered_fasta_input,sequence_t* const seq_read,
    const bool strictly_normalized,const bool try_recovery,const bool check_input_buffer);

/*
 * Error codes
 */
//#define FASTA_ERROR_WRONG_FILE_FORMAT 10   /**/
#define FASTA_ERROR_TAG_BEGINNING 20
#define FASTA_ERROR_READ_BAD_CHARACTER 30
#define FASTA_ERROR_SEPARATOR_BAD_CHARACTER 40
#define FASTA_ERROR_QUALITIES_BAD_CHARACTER 50
#define FASTA_ERROR_LENGTHS 60

/*
 * Error Messages
 */
#define GEM_ERROR_PARSE_FASTQ_READ_BAD_CHARACTER "Parsing FASTA/FASTQ error(%s:%"PRIu64":%"PRIu64"). Input read sequence contains bad character %s"
#define GEM_ERROR_PARSE_FASTQ_QUALITIES_BAD_CHARACTER "Parsing FASTA/FASTQ error(%s:%"PRIu64":%"PRIu64"). Input read qualities contains invalid quality value"

//#define GEM_ERROR_PARSE_FASTQ_BAD_FILE_FORMAT "Input file format is not FASTQ/FASTA"

#define GEM_ERROR_PARSE_FASTA "Parsing FASTA/FASTQ error(%s:%"PRIu64":%"PRIu64")"
//#define GEM_ERROR_WRONG_FILE_FORMAT "Parsing FASTA/FASTQ error(%s:%"PRIu64"). Wrong FASTA/FASTQ syntax"
#define GEM_ERROR_FASTA_TAG_BEGINNING "Parsing FASTA/FASTQ error(%s:%"PRIu64"). Beginning Symbol ('>' or '@') not found. Bad syntax"
#define GEM_ERROR_FASTA_SEPARATOR_BAD_CHARACTER "Parsing FASTA/FASTQ error(%s:%"PRIu64"). Wrong separator symbol ('+'). Bad syntax"
#define GEM_ERROR_FASTA_LENGTHS "Parsing FASTA/FASTQ error(%s:%"PRIu64"). Different read and qualities lengths"

#endif /* INPUT_FASTA_PARSER_H_ */
