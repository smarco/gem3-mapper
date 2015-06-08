/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_generic_parser.h
 * DATE: 28/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Generic parser for {MAP,SAM,FASTQ}
 */

#ifndef GT_INPUT_GENERIC_PARSER_H_
#define GT_INPUT_GENERIC_PARSER_H_

#include "gt_commons.h"
#include "gt_alignment_utils.h"
#include "gt_template_utils.h"

#include "gt_input_file.h"
#include "gt_buffered_input_file.h"

#include "gt_input_parser.h"
#include "gt_input_fasta_parser.h"
#include "gt_input_map_parser.h"
#include "gt_input_sam_parser.h"

#define GT_IGP_FAIL -1
#define GT_IGP_EOF 0
#define GT_IGP_OK 1

/*
 * Parsing Attributes
 */
typedef struct {
  gt_sam_parser_attributes *sam_parser_attributes; /* SAM specific */
  gt_map_parser_attributes *map_parser_attributes; /* MAP specific */
} gt_generic_parser_attributes;

GT_INLINE gt_generic_parser_attributes* gt_input_generic_parser_attributes_new(const bool paired_read);
GT_INLINE void gt_input_generic_parser_attributes_delete(gt_generic_parser_attributes* const attributes);

GT_INLINE void gt_input_generic_parser_attributes_reset_defaults(gt_generic_parser_attributes* const attributes);
GT_INLINE bool gt_input_generic_parser_attributes_is_paired(gt_generic_parser_attributes* const attributes);
GT_INLINE void gt_input_generic_parser_attributes_set_paired(gt_generic_parser_attributes* const attributes,const bool is_paired);

/*
 * Generic Parser
 */
GT_INLINE gt_status gt_input_generic_parser_get_alignment(
    gt_buffered_input_file* const buffered_input,gt_alignment* const alignment,gt_generic_parser_attributes* const attributes);
GT_INLINE gt_status gt_input_generic_parser_get_template(
    gt_buffered_input_file* const buffered_input,gt_template* const template,gt_generic_parser_attributes* const attributes);

/*
 * Synch read of blocks
 */
GT_INLINE gt_status gt_input_generic_parser_synch_blocks_v(
    pthread_mutex_t* const input_mutex,gt_generic_parser_attributes* const attributes,uint64_t num_inputs,
    gt_buffered_input_file* const buffered_input,va_list v_args);
GT_INLINE gt_status gt_input_generic_parser_synch_blocks_va(
    pthread_mutex_t* const input_mutex,gt_generic_parser_attributes* const attributes,
    const uint64_t num_inputs,gt_buffered_input_file* const buffered_input,...);
GT_INLINE gt_status gt_input_generic_parser_synch_blocks_a(
    pthread_mutex_t* const input_mutex,gt_buffered_input_file** const buffered_input,
    const uint64_t num_inputs,gt_generic_parser_attributes* const attributes);

#endif /* GT_INPUT_GENERIC_PARSER_H_ */
