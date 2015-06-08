/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_sam_parser.h
 * DATE: 17/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */


#ifndef GT_INPUT_SAM_PARSER_H_
#define GT_INPUT_SAM_PARSER_H_

#include "gt_commons.h"
#include "gt_dna_string.h"
#include "gt_alignment_utils.h"
#include "gt_template_utils.h"

#include "gt_sequence_archive.h"

#include "gt_input_file.h"
#include "gt_buffered_input_file.h"
#include "gt_input_parser.h"
#include "gt_input_fasta_parser.h"

#include "gt_sam_attributes.h"

// Codes gt_status
#define GT_ISP_OK   GT_STATUS_OK
#define GT_ISP_FAIL GT_STATUS_FAIL
#define GT_ISP_EOF  0

/*
 * Parsing error/state codes
 */
#define GT_ISP_PE_WRONG_FILE_FORMAT 10
#define GT_ISP_PE_PREMATURE_EOL 11
#define GT_ISP_PE_EXPECTED_NUMBER 12
#define GT_ISP_PE_BAD_CHARACTER 14
#define GT_ISP_PE_WRONG_READ_CONTENT 15
/* CIGAR */
#define GT_ISP_PE_CIGAR_PREMATURE_END 20
#define GT_ISP_PE_SAM_UNMAPPED_XA 21
/* PairedEnd Parsing */
#define GT_ISP_PE_WRONG_NUM_XA 30
#define GT_ISP_PE_UNSOLVED_PENDING_MAPS 32

/*
 * SAM file format constants
 */
#define GT_SAM_HEADER_BEGIN '@'

typedef struct {
  bool parse_optional_fields; // TODO
  bool sam_soap_style;
} gt_sam_parser_attributes;
#define GT_SAM_PARSER_ATTR_DEFAULT { .sam_soap_style=false }

GT_INLINE gt_sam_parser_attributes* gt_input_sam_parser_attributes_new();
GT_INLINE void gt_input_sam_parser_attributes_delete(gt_sam_parser_attributes* const attributes);
GT_INLINE void gt_input_sam_parser_attributes_reset_defaults(gt_sam_parser_attributes* const attributes);
GT_INLINE void gt_input_sam_parser_attributes_set_soap_compilant(gt_sam_parser_attributes* const attributes);

/*
 * SAM File basics
 */
GT_INLINE bool gt_input_file_test_sam(
    gt_input_file* const input_file,gt_sam_headers* const sam_headers,const bool show_errors);
GT_INLINE void gt_input_sam_parser_prompt_error(
    gt_buffered_input_file* const buffered_map_input,
    uint64_t line_num,uint64_t column_pos,const gt_status error_code);
GT_INLINE void gt_input_sam_parser_next_record(gt_buffered_input_file* const buffered_map_input);

/*
 * High Level Parsers
 */
GT_INLINE gt_status gt_input_sam_parser_get_template(
    gt_buffered_input_file* const buffered_map_input,gt_template* const template,gt_sam_parser_attributes* const attributes);
GT_INLINE gt_status gt_input_sam_parser_get_alignment(
    gt_buffered_input_file* const buffered_map_input,gt_alignment* const alignment,gt_sam_parser_attributes* const attributes);

/*
 * Error Messages. ISP (Input SAM Parser). General
 */
#define GT_ERROR_PARSE_SAM "Parsing SAM error(%s:%"PRIu64":%"PRIu64")"
#define GT_ERROR_PARSE_SAM_BAD_FILE_FORMAT "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Not a SAM file"
#define GT_ERROR_PARSE_SAM_BAD_CHARACTER "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Bad character found"
#define GT_ERROR_PARSE_SAM_UNMAPPED_XA "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Unmapped read contains XA field (inconsistency)"
#define GT_ERROR_PARSE_SAM_PREMATURE_EOL "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Premature End-of-line found"
#define GT_ERROR_PARSE_SAM_EXPECTED_NUMBER "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Expected number"
#define GT_ERROR_PARSE_SAM_WRONG_READ_CONTENT "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Read in template doesn't match previously parse reads with same tag"
#define GT_ERROR_PARSE_SAM_CIGAR_PREMATURE_END "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Premature end of CIGAR string"
#define GT_ERROR_PARSE_SAM_WRONG_NUM_XA "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Wrong number of eXtra mAps (as to pair them)"
#define GT_ERROR_PARSE_SAM_UNSOLVED_PENDING_MAPS "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Failed to pair maps"

#endif /* GT_INPUT_SAM_PARSER_H_ */
