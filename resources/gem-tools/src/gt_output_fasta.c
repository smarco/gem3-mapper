/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_sam.c
 * DATE: 01/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Output printers to FASTA/FASTQ/MULTIFASTA
 */

#include "gt_output_fasta.h"

GT_INLINE gt_output_fasta_attributes* gt_output_fasta_attributes_new() {
	gt_output_fasta_attributes* attributes = gt_alloc(gt_output_fasta_attributes);
	gt_output_fasta_attributes_reset_defaults(attributes);
	return attributes;
}
GT_INLINE void gt_output_fasta_attributes_delete(gt_output_fasta_attributes* attributes) {
  GT_NULL_CHECK(attributes);
  gt_free(attributes);
}
GT_INLINE void gt_output_fasta_attributes_reset_defaults(gt_output_fasta_attributes* const attributes) {
	attributes->print_extra = true;
	attributes->print_casava = true;
	attributes->format = F_FASTQ;
}
GT_INLINE bool gt_output_fasta_attributes_is_print_extra(gt_output_fasta_attributes* const attributes) {
	return attributes->print_extra;
}
GT_INLINE void gt_output_fasta_attributes_set_print_extra(gt_output_fasta_attributes* const attributes,const bool print_extra) {
	attributes->print_extra = print_extra;
}
GT_INLINE bool gt_output_fasta_attributes_is_print_casava(gt_output_fasta_attributes* const attributes) {
	return attributes->print_casava;
}
GT_INLINE void gt_output_fasta_attributes_set_print_casava(gt_output_fasta_attributes* const attributes,const bool print_casava) {
	attributes->print_casava = print_casava;
}
GT_INLINE gt_file_fasta_format gt_output_fasta_attributes_get_format(gt_output_fasta_attributes* const attributes) {
	return attributes->format;
}
GT_INLINE void gt_output_fasta_attributes_set_format(gt_output_fasta_attributes* const attributes,const gt_file_fasta_format format) {
	attributes->format = format;
}

/*
 * FASTA building block printers
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS is_fasta,tag,attributes,output_attributes
GT_INLINE gt_status gt_output_fasta_gprint_tag(gt_generic_printer* const gprinter,
    bool const is_fasta,gt_string* const tag,gt_attributes* const attributes,gt_output_fasta_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  // Print begin character & TAG
  if (is_fasta) {
    gt_gprintf(gprinter,">"PRIgts,PRIgts_content(tag));
  } else {
    gt_gprintf(gprinter,"@"PRIgts,PRIgts_content(tag));
  }
  // Print TAG Attributes
  if (attributes!=NULL) {
    gt_output_gprint_tag_attributes(gprinter,attributes,
        gt_output_fasta_attributes_is_print_casava(output_attributes),
        gt_output_fasta_attributes_is_print_extra(output_attributes));
  }
  // Print EOL
  gt_gprintf(gprinter,"\n");
  return 0;
}

#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,attributes,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_fasta,print_fasta,
    gt_string* const tag,gt_string* const read,gt_attributes* const attributes,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_gprint_fasta(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_attributes* const attributes,gt_output_fasta_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  GT_STRING_CHECK(read);
  GT_ATTRIBUTES_CHECK(attributes);
  gt_output_fasta_gprint_tag(gprinter,true,tag,attributes,output_attributes);
  gt_gprintf(gprinter,PRIgts"\n",PRIgts_content(read));
  return 0;
}

#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities,attributes,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_fasta,print_fastq,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_attributes* const attributes,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_gprint_fastq(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_attributes* const attributes,gt_output_fasta_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  GT_STRING_CHECK(read);
  //gt_gprintf(gprinter,"@"PRIgts"\n",PRIgts_content(tag));
  gt_output_fasta_gprint_tag(gprinter, false, tag, attributes, output_attributes);
  gt_gprintf(gprinter,PRIgts"\n",PRIgts_content(read));
  if (!gt_string_is_null(qualities)) {
    gt_gprintf(gprinter,"+\n"PRIgts"\n",PRIgts_content(qualities));
  } else { // Print dummy qualities
    const uint64_t read_length = gt_string_get_length(read);
    uint64_t i;
    for (i=0;i<read_length;++i) gt_gprintf(gprinter,"X");
    gt_gprintf(gprinter,"\n");
  }
  // TODO attributes
  return 0;
}


/*
 * FASTA High-level Printers
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS dna_read,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_fasta,print_dna_read,
    gt_dna_read* const dna_read,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_gprint_dna_read(gt_generic_printer* const gprinter,
    gt_dna_read* const dna_read,gt_output_fasta_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_DNA_READ_CHECK(dna_read);
  switch (gt_output_fasta_attributes_get_format(output_attributes)) {
    case F_FASTA:
      return gt_output_fasta_gprint_fasta(gprinter,dna_read->tag,dna_read->read,dna_read->attributes,output_attributes);
      break;
    case F_FASTQ:
      return gt_output_fasta_gprint_fastq(gprinter,dna_read->tag,dna_read->read,dna_read->qualities,dna_read->attributes,output_attributes);
      break;
    default:
      gt_fatal_error(SELECTION_NOT_VALID);
      break;
  }
}

#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_fasta,print_alignment,
    gt_alignment* const alignment,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_gprint_alignment(gt_generic_printer* const gprinter,
    gt_alignment* const alignment,gt_output_fasta_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  switch (gt_output_fasta_attributes_get_format(output_attributes)) {
    case F_FASTA:
      return gt_output_fasta_gprint_fasta(gprinter,alignment->tag,alignment->read,alignment->attributes, output_attributes);
      break;
    case F_FASTQ:
      return gt_output_fasta_gprint_fastq(gprinter,alignment->tag,alignment->read,alignment->qualities,alignment->attributes, output_attributes);
      break;
    default:
      gt_fatal_error(SELECTION_NOT_VALID);
      break;
  }
}

#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_fasta,print_template,
   gt_template* const template,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_gprint_template(gt_generic_printer* const gprinter,
   gt_template* const template,gt_output_fasta_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  gt_status error_code;
  GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
    error_code=gt_output_fasta_gprint_alignment(gprinter,alignment, output_attributes);
    if (error_code) return error_code;
  }
  return 0;
}

#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS sequence_archive,column_width,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_fasta,print_sequence_archive,
    gt_sequence_archive* const sequence_archive,const uint64_t column_width,gt_output_fasta_attributes* const output_attributes);
GT_INLINE gt_status gt_output_fasta_gprint_sequence_archive(gt_generic_printer* const gprinter,
    gt_sequence_archive* const sequence_archive,const uint64_t column_width,gt_output_fasta_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  // Dump the content of the reference file (using Iterators)
  gt_segmented_sequence* seq;
  gt_sequence_archive_iterator seq_arch_it;
  gt_sequence_archive_new_iterator(sequence_archive,&seq_arch_it);
  while ((seq=gt_sequence_archive_iterator_next(&seq_arch_it))) {
    // Print TAG
    gt_gprintf(gprinter,">"PRIgts"\n",PRIgts_content(seq->seq_name));
	  // Print Sequences
    gt_segmented_sequence_iterator sequence_iterator;
    gt_segmented_sequence_new_iterator(seq,0,GT_ST_FORWARD,&sequence_iterator);
    uint64_t chars_written = 0;
    if (!gt_segmented_sequence_iterator_eos(&sequence_iterator)) {
      while (!gt_segmented_sequence_iterator_eos(&sequence_iterator)) {
        gt_gprintf(gprinter,"%c",gt_segmented_sequence_iterator_next(&sequence_iterator));
        if ((++chars_written)%column_width==0) gt_gprintf(gprinter,"\n");
      }
      gt_gprintf(gprinter,"\n");
    }
  }
  return 0;
}

