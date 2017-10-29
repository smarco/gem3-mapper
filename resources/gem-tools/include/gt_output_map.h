/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_map.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_OUTPUT_MAP_H_
#define GT_OUTPUT_MAP_H_

#include "gt_essentials.h"
#include "gt_template.h"
#include "gt_sequence_archive.h"

#include "gt_generic_printer.h"
#include "gt_buffered_output_file.h"
#include "gt_output_printer.h"

/*
 * Error/state codes (Map Output Error)
 */
#define GT_MOE_INCONSISTENT_COUNTERS 10
#define GT_MOE_ERROR_PRINTING_MISM_STRING 20
#define GT_MOE_ERROR_PRINTING_MAP_BLOCKS 30

/*
 * Output attributes
 */
typedef struct {
  /* Tag */
  bool print_extra; // print extra information stored in attributes
  bool print_casava; // if available print casava ids, otherwise, appends /1 /2 for paired reads
  /* Counters */
  bool compact;
  /* Maps */
  bool print_scores; // Print alignment scores
  bool hex_print_scores; // print alignment scores in hex
  uint64_t max_printable_maps; // Maximum number of maps printed
} gt_output_map_attributes;
#define GT_OUTPUT_MAP_ATTR_DEFAULT() { \
   /* Tag */ \
  .print_extra=true, \
  .print_casava=true, \
   /* Counters */ \
  .compact=false, \
   /* Maps */ \
  .print_scores=true, \
  .max_printable_maps=GT_ALL \
}
#define GT_OUTPUT_MAP_CHECK_ATTRIBUTES(attributes) \
  gt_output_map_attributes __##attributes; \
  if (attributes==NULL) { /* Check null output_map_attributes */ \
    gt_output_map_attributes_reset_defaults(&__##attributes); \
    attributes = &__##attributes; \
  }

GT_INLINE gt_output_map_attributes* gt_output_map_attributes_new();
GT_INLINE void gt_output_map_attributes_delete(gt_output_map_attributes* const attributes);
GT_INLINE void gt_output_map_attributes_reset_defaults(gt_output_map_attributes* const attributes);
/* Tag */
GT_INLINE bool gt_output_map_attributes_is_print_casava(gt_output_map_attributes* const attributes);
GT_INLINE void gt_output_map_attributes_set_print_casava(gt_output_map_attributes* const attributes,const bool print_casava);
GT_INLINE bool gt_output_map_attributes_is_print_extra(gt_output_map_attributes* const attributes);
GT_INLINE void gt_output_map_attributes_set_print_extra(gt_output_map_attributes* const attributes,const bool print_extra);
/* Maps */
GT_INLINE bool gt_output_map_attributes_is_print_scores(gt_output_map_attributes* const attributes);
GT_INLINE void gt_output_map_attributes_set_print_scores(gt_output_map_attributes* const attributes,const bool print_scores);
GT_INLINE uint64_t gt_output_map_attributes_get_max_printable_maps(gt_output_map_attributes* const attributes);
GT_INLINE void gt_output_map_attributes_set_max_printable_maps(gt_output_map_attributes* const attributes,const uint64_t max_printable_maps);


/*
 * TAG building block printers
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_tag,gt_string* const tag,
    gt_attributes* const attributes,gt_output_map_attributes* const output_map_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_template_tag,gt_template* const template);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_alignment_tag,gt_alignment* const alignment);

/*
 * MAP building block printers
 *   - If @gt_output_map_attributes==NULL then defaults are applied
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_mismatch_string,gt_map* const map,gt_output_map_attributes* output_map_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_counters,gt_vector* const counters,gt_attributes* const attributes,gt_output_map_attributes* output_map_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_map_block,gt_map* const map,gt_output_map_attributes* output_map_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_map,gt_map* const map,gt_output_map_attributes* output_map_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_map_list,gt_vector* const maps,gt_output_map_attributes* output_map_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_map_placeholder,gt_vector* const mmap_placeholder,gt_output_map_attributes* output_map_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_template_maps,gt_template* const template,gt_output_map_attributes* output_map_attributes); // Sorted print by mismatch number
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_alignment_maps,gt_alignment* const alignment,gt_output_map_attributes* output_map_attributes); // Sorted print by mismatch number

/*
 * High-level MAP Printers {Alignment/Template}
 */
GT_INLINE gt_status gt_output_map_gprint_template(gt_generic_printer* const gprinter,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_fprint_template(FILE* file,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_sprint_template(gt_string* const string,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_bprint_template(gt_output_buffer* const output_buffer,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_ofprint_template(gt_output_file* const output_file,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_bofprint_template(gt_buffered_output_file* const buffered_output_file,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_alignment(gt_generic_printer* const gprinter,gt_alignment* const alignment,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_fprint_alignment(FILE* file,gt_alignment* const alignment,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_sprint_alignment(gt_string* const string,gt_alignment* const alignment,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_bprint_alignment(gt_output_buffer* const output_buffer,gt_alignment* const alignment,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_ofprint_alignment(gt_output_file* const output_file,gt_alignment* const alignment,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_bofprint_alignment(gt_buffered_output_file* const buffered_output_file,gt_alignment* const alignment,gt_output_map_attributes* const output_map_attributes);

/*
 * GEM printer
 *   If the template is paired generates 1 PairedEnd MAP line
 *   If the template is unpaired generates 2 SingleEnd MAP lines
 */
GT_INLINE gt_status gt_output_map_gprint_gem_template(gt_generic_printer* const gprinter,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_fprint_gem_template(FILE* file,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_sprint_gem_template(gt_string* const string,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_bprint_gem_template(gt_output_buffer* const output_buffer,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_ofprint_gem_template(gt_output_file* const output_file,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_bofprint_gem_template(gt_buffered_output_file* const buffered_output_file,gt_template* const template,gt_output_map_attributes* const output_map_attributes);

/*
 * Misc. Handy printers
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_mismatch_summary,gt_map* const map);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_map_block_pretty,gt_map* const map,
    char* const pattern,const uint64_t pattern_length,char* const sequence,const uint64_t sequence_length);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_map,print_map_pretty_sa,gt_map* const map,
    gt_string* const pattern,gt_sequence_archive* const sequence_archive);

#endif /* GT_OUTPUT_MAP_H_ */
