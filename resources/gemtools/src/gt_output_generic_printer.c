/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_generic_printer.c
 * DATE: 28/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Generic printer for {FASTA,FASTQ,MAP,SAM}
 */

#include "gt_output_generic_printer.h"

/*
 * Generic Printer Attributes
 */
GT_INLINE gt_generic_printer_attributes* gt_generic_printer_attributes_new(const gt_file_format file_format) {
  gt_generic_printer_attributes* attributes = gt_alloc(gt_generic_printer_attributes);
  attributes->output_sam_attributes = NULL;
  attributes->output_fasta_attributes = NULL;
  attributes->output_map_attributes = NULL;
  gt_generic_printer_attributes_set_format(attributes,file_format);
  return attributes;
}
GT_INLINE void gt_generic_printer_attributes_delete(gt_generic_printer_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  if (attributes->output_sam_attributes!=NULL) gt_output_sam_attributes_delete(attributes->output_sam_attributes);
  if (attributes->output_fasta_attributes!=NULL) gt_output_fasta_attributes_delete(attributes->output_fasta_attributes);
  if (attributes->output_map_attributes!=NULL) gt_output_map_attributes_delete(attributes->output_map_attributes);
  gt_free(attributes);
}
GT_INLINE void gt_generic_printer_attributes_set_format(
    gt_generic_printer_attributes* const attributes,const gt_file_format file_format) {
  GT_NULL_CHECK(attributes);
  switch (file_format) {
    case SAM:
      attributes->output_format = SAM;
      attributes->output_sam_attributes = gt_output_sam_attributes_new();
      break;
    case FASTA:
      attributes->output_format = FASTA;
      attributes->output_fasta_attributes = gt_output_fasta_attributes_new();
      break;
    case MAP:
    default:
      attributes->output_format = MAP;
      attributes->output_map_attributes = gt_output_map_attributes_new();
      break;
  }
}

/*
 * Generic Printer
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_generic,print_alignment,gt_alignment* const alignment,gt_generic_printer_attributes* const attributes);
GT_INLINE gt_status gt_output_generic_gprint_alignment(gt_generic_printer* const gprinter,
    gt_alignment* const alignment,gt_generic_printer_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(attributes);
  switch (attributes->output_format) {
    case SAM:
      gt_output_sam_gprint_alignment(gprinter,alignment,attributes->output_sam_attributes);
      break;
    case FASTA:
      gt_output_fasta_gprint_alignment(gprinter,alignment,attributes->output_fasta_attributes);
      break;
    case MAP:
    default:
      gt_output_map_gprint_alignment(gprinter,alignment,attributes->output_map_attributes);
      break;
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_generic,print_template,gt_template* const template,gt_generic_printer_attributes* const attributes);
GT_INLINE gt_status gt_output_generic_gprint_template(gt_generic_printer* const gprinter,
    gt_template* const template,gt_generic_printer_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(attributes);
  switch (attributes->output_format) {
    case SAM:
      gt_output_sam_gprint_template(gprinter,template,attributes->output_sam_attributes);
      break;
    case FASTA:
      gt_output_fasta_gprint_template(gprinter,template,attributes->output_fasta_attributes);
      break;
    case MAP:
    default:
      gt_output_map_gprint_gem_template(gprinter,template,attributes->output_map_attributes);
      break;
  }
  return 0;
}
