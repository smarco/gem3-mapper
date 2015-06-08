/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_generic_printer.h
 * DATE: 28/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Generic printer for {FASTA,FASTQ,MAP,SAM}
 */

#ifndef GT_OUTPUT_GENERIC_PRINTER_H_
#define GT_OUTPUT_GENERIC_PRINTER_H_

#include "gt_essentials.h"
#include "gt_template.h"

#include "gt_input_file.h"

#include "gt_generic_printer.h"
#include "gt_buffered_output_file.h"
#include "gt_output_fasta.h"
#include "gt_output_map.h"
#include "gt_output_sam.h"


/*
 * Generic Printer Attributes
 */
typedef struct {
  /* Format */
  gt_file_format output_format;
  /* Format Specific Attributes */
  gt_output_map_attributes *output_map_attributes;
  gt_output_sam_attributes *output_sam_attributes;
  gt_output_fasta_attributes *output_fasta_attributes;
} gt_generic_printer_attributes;

GT_INLINE gt_generic_printer_attributes* gt_generic_printer_attributes_new(const gt_file_format file_format);
GT_INLINE void gt_generic_printer_attributes_clear(gt_generic_printer_attributes* const attributes);
GT_INLINE void gt_generic_printer_attributes_delete(gt_generic_printer_attributes* const attributes);
GT_INLINE void gt_generic_printer_attributes_set_format(gt_generic_printer_attributes* const attributes,const gt_file_format file_format);

/*
 * Generic Printer
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_generic,print_alignment,gt_alignment* const alignment,gt_generic_printer_attributes* const attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_generic,print_template,gt_template* const template,gt_generic_printer_attributes* const attributes);

#endif /* GT_OUTPUT_GENERIC_PRINTER_H_ */
