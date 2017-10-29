/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_fasta_parser.h
 * DATE: 17/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Input parser for FASTA/FASTQ/MULTIFASTA files
 */

#ifndef GT_INPUT_FASTA_PARSER_H_
#define GT_INPUT_FASTA_PARSER_H_

#include "gt_commons.h"
#include "gt_dna_read.h"
#include "gt_sequence_archive.h"

#include "gt_input_file.h"
#include "gt_buffered_input_file.h"
#include "gt_input_parser.h"

// Codes gt_status
#define GT_IFP_OK   GT_STATUS_OK
#define GT_IFP_FAIL GT_STATUS_FAIL
#define GT_IFP_EOF  0

/*
 * Parsing error/state codes
 */
#define GT_IFP_PE_WRONG_FILE_FORMAT 10
#define GT_IFP_PE_PREMATURE_EOL 11
#define GT_IFP_PE_PREMATURE_EOB 12
#define GT_IFP_PE_BAD_CHARACTER 13

#define GT_IFP_PE_TAG_BAD_BEGINNING 20

#define GT_IFP_PE_READ_BAD_CHARACTER 30

#define GT_IFP_PE_SEPARATOR_BAD_CHARACTER 40

#define GT_IFP_PE_QUALS_BAD_CHARACTER 50
#define GT_IFP_PE_QUALS_BAD_LENGTH 51

/*
 * FASTQ File basics
 */
GT_INLINE bool gt_input_file_test_fasta(
    gt_input_file* const input_file,gt_fasta_file_format* const fasta_file_format,const bool show_errors);
GT_INLINE void gt_input_fasta_parser_prompt_error(
    gt_buffered_input_file* const buffered_fasta_input,
    uint64_t line_num,uint64_t column_pos,const gt_status error_code);
GT_INLINE void gt_input_fasta_parser_next_record(gt_buffered_input_file* const buffered_fasta_input,char* const line_start);

#define gt_input_fasta_is_fasta(input_file) (input_file->file_format==FASTA && input_file->fasta_type.fasta_format==F_FASTA)
#define gt_input_fasta_is_fastq(input_file) (input_file->file_format==FASTA && input_file->fasta_type.fasta_format==F_FASTQ)
#define gt_input_fasta_is_multifasta(input_file) (input_file->file_format==FASTA && input_file->fasta_type.fasta_format==F_MULTI_FASTA)

#define gt_input_fasta_get_format(input_file) (input_file->fasta_type.fasta_format)

/*
 * High Level Parsers
 */
GT_INLINE gt_status gt_input_fasta_parser_get_read(
    gt_buffered_input_file* const buffered_fasta_input,gt_dna_read* const dna_read);
GT_INLINE gt_status gt_input_fasta_parser_get_alignment(
    gt_buffered_input_file* const buffered_fasta_input,gt_alignment* const alignment);
GT_INLINE gt_status gt_input_fasta_parser_get_template(
    gt_buffered_input_file* const buffered_fasta_input,gt_template* const template,const bool paired_read);

GT_INLINE gt_status gt_input_multifasta_parser_get_archive(
    gt_input_file* const input_multifasta_file,gt_sequence_archive* const sequence_archive);

/*
 * Synch read of blocks
 */
GT_INLINE gt_status gt_input_fasta_parser_synch_blocks(
    gt_buffered_input_file* const buffered_input1,gt_buffered_input_file* const buffered_input2,pthread_mutex_t* const input_mutex);
GT_INLINE gt_status gt_input_fasta_parser_synch_blocks_v(
    pthread_mutex_t* const input_mutex,uint64_t num_inputs,gt_buffered_input_file* const buffered_input,va_list v_args);
GT_INLINE gt_status gt_input_fasta_parser_synch_blocks_va(
    pthread_mutex_t* const input_mutex,const uint64_t num_inputs,gt_buffered_input_file* const buffered_input,...);
GT_INLINE gt_status gt_input_fasta_parser_synch_blocks_a(
    pthread_mutex_t* const input_mutex,gt_buffered_input_file** const buffered_input,const uint64_t num_inputs);

/*
 * Error Messages
 */
// IFP (Input FASTA Parser). General
#define GT_ERROR_PARSE_FASTA "Parsing FASTA/FASTQ error(%s:%"PRIu64":%"PRIu64")"
#define GT_ERROR_PARSE_FASTA_BAD_FILE_FORMAT "Parsing FASTA/FASTQ error(%s:%"PRIu64"). Not a FASTA/FASTQ file"

#endif /* GT_INPUT_FASTA_PARSER_H_ */
