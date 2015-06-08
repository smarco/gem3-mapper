/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_generic_parser.c
 * DATE: 28/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Generic parser for {MAP,SAM,FASTQ}
 */

#include "gt_input_generic_parser.h"

/*
 * Setup
 */
GT_INLINE gt_generic_parser_attributes* gt_input_generic_parser_attributes_new(const bool paired_reads) {
  gt_generic_parser_attributes* attributes = gt_alloc(gt_generic_parser_attributes);
  /* MAP */
  attributes->map_parser_attributes = gt_input_map_parser_attributes_new(paired_reads);
  /* SAM */
  attributes->sam_parser_attributes = gt_input_sam_parser_attributes_new();
  return attributes;
}
GT_INLINE void gt_input_generic_parser_attributes_delete(gt_generic_parser_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  gt_input_map_parser_attributes_delete(attributes->map_parser_attributes);
  gt_input_sam_parser_attributes_delete(attributes->sam_parser_attributes);
  gt_free(attributes);
}

/*
 * Accessors
 */
GT_INLINE void gt_input_generic_parser_attributes_reset_defaults(gt_generic_parser_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  gt_input_map_parser_attributes_reset_defaults(attributes->map_parser_attributes);
  gt_input_sam_parser_attributes_reset_defaults(attributes->sam_parser_attributes);
}
GT_INLINE bool gt_input_generic_parser_attributes_is_paired(gt_generic_parser_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  return gt_input_map_parser_attributes_is_paired(attributes->map_parser_attributes);
}
GT_INLINE void gt_input_generic_parser_attributes_set_paired(gt_generic_parser_attributes* const attributes,const bool is_paired) {
  GT_NULL_CHECK(attributes);
  gt_input_map_parser_attributes_set_paired(attributes->map_parser_attributes,is_paired);
}

/*
 * Parsers Helpers
 */
GT_INLINE gt_status gt_input_generic_parser_get_alignment(
    gt_buffered_input_file* const buffered_input,gt_alignment* const alignment,gt_generic_parser_attributes* const attributes) {
  gt_status error_code = GT_IGP_FAIL;
  switch (buffered_input->input_file->file_format) {
    case SAM:
      return gt_input_sam_parser_get_alignment(buffered_input,alignment,attributes->sam_parser_attributes);
      break;
    case FASTA:
      return gt_input_fasta_parser_get_alignment(buffered_input,alignment);
      break;
    case MAP:
    default: // gt_fatal_error_msg("File type not supported");
      return gt_input_map_parser_get_alignment(buffered_input,alignment,attributes->map_parser_attributes);
      break;
  }
  return error_code;
}
GT_INLINE gt_status gt_input_generic_parser_get_template(
    gt_buffered_input_file* const buffered_input,gt_template* const template,gt_generic_parser_attributes* const attributes) {
  gt_status error_code = GT_IGP_FAIL;
  switch (buffered_input->input_file->file_format) {
    case SAM:
      if (gt_input_generic_parser_attributes_is_paired(attributes)) {
        error_code = gt_input_sam_parser_get_template(buffered_input,template,attributes->sam_parser_attributes);
        gt_template_get_block_dyn(template,0);
        gt_template_get_block_dyn(template,1); // Make sure is a template
        return error_code;
      } else {
        return gt_input_sam_parser_get_alignment(
            buffered_input,gt_template_get_block_dyn(template,0),attributes->sam_parser_attributes);
      }
      break;
    case FASTA:
      return gt_input_fasta_parser_get_template(buffered_input,template,gt_input_generic_parser_attributes_is_paired(attributes));
      break;
    case MAP:
    default: // gt_fatal_error_msg("File type not supported");
      return gt_input_map_parser_get_template(buffered_input,template,attributes->map_parser_attributes);
      break;
  }
  return error_code;
}


/*
 * Synch read of blocks
 */
GT_INLINE gt_status gt_input_generic_parser_synch_blocks_v(
    pthread_mutex_t* const input_mutex,gt_generic_parser_attributes* const attributes,uint64_t num_inputs,
    gt_buffered_input_file* const buffered_input,va_list v_args) {
  gt_status error_code = GT_IGP_FAIL;
  switch (buffered_input->input_file->file_format) {
    case MAP:
      return gt_input_map_parser_synch_blocks_v(input_mutex,attributes->map_parser_attributes,num_inputs,buffered_input,v_args);
      break;
    case SAM:
      GT_NOT_IMPLEMENTED();
      break;
    case FASTA:
      return gt_input_fasta_parser_synch_blocks_v(input_mutex,num_inputs,buffered_input,v_args);
      break;
    default:
      gt_fatal_error_msg("File type not supported");
      break;
  }
  return error_code;
}
GT_INLINE gt_status gt_input_generic_parser_synch_blocks_va(
    pthread_mutex_t* const input_mutex,gt_generic_parser_attributes* const attributes,
    const uint64_t num_inputs,gt_buffered_input_file* const buffered_input,...) {
  GT_NULL_CHECK(input_mutex);
  GT_NULL_CHECK(attributes);
  GT_ZERO_CHECK(num_inputs);
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_input);
  va_list v_args;
  va_start(v_args,buffered_input);
  gt_status error_code = gt_input_generic_parser_synch_blocks_v(input_mutex,attributes,num_inputs,buffered_input,v_args);
  va_end(v_args);
  return error_code;
}
GT_INLINE gt_status gt_input_generic_parser_synch_blocks_a(
    pthread_mutex_t* const input_mutex,gt_buffered_input_file** const buffered_input,
    const uint64_t num_inputs,gt_generic_parser_attributes* const attributes) {
  gt_status error_code = GT_IGP_FAIL;
  switch (buffered_input[0]->input_file->file_format) {
    case MAP:
      return gt_input_map_parser_synch_blocks_a(input_mutex,buffered_input,num_inputs,attributes->map_parser_attributes);
      break;
    case SAM:
      GT_NOT_IMPLEMENTED();
      break;
    case FASTA:
      return gt_input_fasta_parser_synch_blocks_a(input_mutex,buffered_input,num_inputs);
      break;
    default:
      gt_fatal_error_msg("File type not supported");
      break;
  }
  return error_code;
}
