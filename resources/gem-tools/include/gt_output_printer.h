/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_printer.h
 * DATE: 01/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef GT_OUTPUT_PRINTER_H_
#define GT_OUTPUT_PRINTER_H_

#include "gt_essentials.h"

#include "gt_output_buffer.h"
#include "gt_buffered_output_file.h"
#include "gt_generic_printer.h"
#include "gt_attributes.h"

// Print Segmented Read Attribute
GT_GENERIC_PRINTER_PROTOTYPE(gt_output,print_segmented_read_info,const uint64_t segment_id,const uint64_t total_segments);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output,print_segmented_read,gt_segmented_read_info* const segmented_read_info);

// Print Trim Attributes
GT_GENERIC_PRINTER_PROTOTYPE(gt_output,print_trim,
    const gt_read_trim_t trim_type,gt_string* const trimmed_read,gt_string* const trimmed_qualities,const uint64_t length);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output,print_right_trim,gt_read_trim* const right_trim);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output,print_left_trim,gt_read_trim* const right_trim);

// Prints all TAG attributes (trims,sg,casava,extras,...)
GT_GENERIC_PRINTER_PROTOTYPE(gt_output,print_tag_attributes,
    gt_attributes* const attributes,const bool print_casava_flags,const bool print_extra_tag_attributes);

#endif /* GT_OUTPUT_PRINTER_H_ */
