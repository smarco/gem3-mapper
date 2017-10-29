/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_fasta.h
 * DATE: 01/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Output printers to FASTA/FASTQ/MULTIFASTA
 */

#ifndef GT_OUTPUT_FASTA_H_
#define GT_OUTPUT_FASTA_H_

#include "gt_essentials.h"
#include "gt_dna_read.h"
#include "gt_sequence_archive.h"
#include "gt_input_file.h"

#include "gt_alignment.h"
#include "gt_template.h"

#include "gt_input_parser.h"
#include "gt_output_buffer.h"
#include "gt_buffered_output_file.h"
#include "gt_generic_printer.h"
#include "gt_output_printer.h"
#include "gt_attributes.h"

/*
 * Output attributes
 */
typedef struct {
  bool print_extra;            // Print extra information stored in attributes
  bool print_casava;           // If available print CASAVA IDs, otherwise, appends /1 /2 for paired reads
  gt_file_fasta_format format; // File format {F_FASTA, F_FASTQ, F_MULTI_FASTA}
} gt_output_fasta_attributes;

GT_INLINE gt_output_fasta_attributes* gt_output_fasta_attributes_new();
GT_INLINE void gt_output_fasta_attributes_delete(gt_output_fasta_attributes* attributes);
GT_INLINE void gt_output_fasta_attributes_reset_defaults(gt_output_fasta_attributes* const attributes);

GT_INLINE bool gt_output_fasta_attributes_is_print_extra(gt_output_fasta_attributes* const attributes);
GT_INLINE void gt_output_fasta_attributes_set_print_extra(gt_output_fasta_attributes* const attributes,const bool print_extra);
GT_INLINE bool gt_output_fasta_attributes_is_print_casava(gt_output_fasta_attributes* const attributes);
GT_INLINE void gt_output_fasta_attributes_set_print_casava(gt_output_fasta_attributes* const attributes,const bool print_casava);
GT_INLINE gt_file_fasta_format gt_output_fasta_attributes_get_format(gt_output_fasta_attributes* const attributes);
GT_INLINE void gt_output_fasta_attributes_set_format(gt_output_fasta_attributes* const attributes,const gt_file_fasta_format format);

/*
 * FASTA building block printers
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_fasta,print_tag,const bool is_fasta,gt_string* const tag,gt_attributes* const attributes,gt_output_fasta_attributes* const output_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_fasta,print_fasta,gt_string* const tag,gt_string* const read,gt_attributes* const attributes,gt_output_fasta_attributes* const output_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_fasta,print_fastq,gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_attributes* const attributes,gt_output_fasta_attributes* const output_attributes);

/*
 * FASTA High-level Printers
 */
// GT_GENERIC_PRINTER_PROTOTYPE(gt_output_fasta,print_dna_read,gt_file_fasta_format fasta_format,gt_dna_read* const dna_read,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_gprint_dna_read(gt_generic_printer* const gprinter,gt_dna_read* const dna_read,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_fprint_dna_read(FILE* file,gt_dna_read* const dna_read,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_sprint_dna_read(gt_string* const string,gt_dna_read* const dna_read,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_bprint_dna_read(gt_output_buffer* const output_buffer,gt_dna_read* const dna_read,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_ofprint_dna_read(gt_output_file* const output_file,gt_dna_read* const dna_read,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_bofprint_dna_read(gt_buffered_output_file* const buffered_output_file,gt_dna_read* const dna_read,gt_output_fasta_attributes* const output_attributes);

// GT_GENERIC_PRINTER_PROTOTYPE(gt_output_fasta,print_alignment,gt_alignment* const alignment,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_gprint_alignment(gt_generic_printer* const gprinter,gt_alignment* const alignment,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_fprint_alignment(FILE* file,gt_alignment* const alignment,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_sprint_alignment(gt_string* const string,gt_alignment* const alignment,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_bprint_alignment(gt_output_buffer* const output_buffer,gt_alignment* const alignment,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_ofprint_alignment(gt_output_file* const output_file,gt_alignment* const alignment,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_bofprint_alignment(gt_buffered_output_file* const buffered_output_file,gt_alignment* const alignment,gt_output_fasta_attributes* const output_attributes);

// GT_GENERIC_PRINTER_PROTOTYPE(gt_output_fasta,print_template,gt_template* const template,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_gprint_template(gt_generic_printer* const gprinter,gt_template* const template,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_fprint_template(FILE* file,gt_template* const template,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_sprint_template(gt_string* const string,gt_template* const template,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_bprint_template(gt_output_buffer* const output_buffer,gt_template* const template,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_ofprint_template(gt_output_file* const output_file,gt_template* const template,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_bofprint_template(gt_buffered_output_file* const buffered_output_file,gt_template* const template,gt_output_fasta_attributes* const output_attributes);

// GT_GENERIC_PRINTER_PROTOTYPE(gt_output_fasta,print_sequence_archive,gt_sequence_archive* const sequence_archive,const uint64_t column_width,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_gprint_sequence_archive(gt_generic_printer* const gprinter,gt_sequence_archive* const sequence_archive,const uint64_t column_width,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_fprint_sequence_archive(FILE* file,gt_sequence_archive* const sequence_archive,const uint64_t column_width,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_sprint_sequence_archive(gt_string* const string,gt_sequence_archive* const sequence_archive,const uint64_t column_width,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_bprint_sequence_archive(gt_output_buffer* const output_buffer,gt_sequence_archive* const sequence_archive,const uint64_t column_width,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_ofprint_sequence_archive(gt_output_file* const output_file,gt_sequence_archive* const sequence_archive,const uint64_t column_width,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_bofprint_sequence_archive(gt_buffered_output_file* const buffered_output_file,gt_sequence_archive* const sequence_archive,const uint64_t column_width,gt_output_fasta_attributes* const output_attributes);

#endif /* GT_OUTPUT_FASTA_H_ */
