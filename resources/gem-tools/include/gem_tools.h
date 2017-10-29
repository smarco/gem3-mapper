/*
 * PROJECT: GEM-Tools library
 * FILE: gem_tools.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GEM_TOOLS_H_
#define GEM_TOOLS_H_

// Essentials
#include "gt_essentials.h"
#include "gt_dna_string.h"

// Input handlers
#include "gt_input_file.h"
#include "gt_buffered_input_file.h"
// Input parsers/utils
#include "gt_input_parser.h"
#include "gt_input_map_parser.h"
#include "gt_input_map_utils.h"
#include "gt_input_sam_parser.h"
#include "gt_input_fasta_parser.h"
#include "gt_input_generic_parser.h"

// Output handlers
#include "gt_output_buffer.h"
#include "gt_buffered_output_file.h"
// Output printers (MAP,SAM,BAM,...)
#include "gt_output_fasta.h"
#include "gt_output_map.h"
#include "gt_output_sam.h"
#include "gt_output_generic_printer.h"

// GEM-Tools basic data structures: Template/Alignment/Maps/...
#include "gt_misms.h"
#include "gt_map.h"
#include "gt_dna_read.h"
#include "gt_attributes.h"
#include "gt_alignment.h"
#include "gt_alignment_utils.h"
#include "gt_template.h"
#include "gt_template_utils.h"
#include "gt_counters_utils.h"
#include "gt_compact_dna_string.h"
#include "gt_sequence_archive.h"

// HighLevel Modules
#include "gt_stats.h"
#include "gt_gtf.h"

// Utilities
#include "gt_options_menu.h"
#include "gt_json.h"

// GEM Idx Loader
#include "gt_gemIdx_loader.h"

/*
 * General generic I/O loop (overkilling, but useful)
 */
#define GT_BEGIN_READING_WRITING_LOOP(input_file,output_file,paired_end,buffered_output,template) \
  /* Prepare IN/OUT buffers & printers */ \
  gt_status __error_code; \
  gt_buffered_input_file* __buffered_input = gt_buffered_input_file_new(input_file); \
  gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file); \
  gt_buffered_input_file_attach_buffered_output(__buffered_input,buffered_output); \
  /* Prepare Attributes for generic I/O */ \
  gt_generic_parser_attributes* __gparser_attr = gt_input_generic_parser_attributes_new(paired_end); \
  /* I/O Loop */ \
  gt_template* template = gt_template_new(); \
  while ((__error_code=gt_input_generic_parser_get_template(__buffered_input,template,__gparser_attr))) { \
    if (__error_code!=GT_INPUT_STATUS_OK) { \
      gt_error_msg("Fatal error parsing file '%s', line %"PRIu64"\n", \
          gt_input_file_get_file_name(input_file),__buffered_input->current_line_num-1); \
      continue; \
    }
#define GT_END_READING_WRITING_LOOP(input_file,output_file,template) \
  } \
  /* Clean */ \
  gt_buffered_input_file_close(__buffered_input); \
  gt_buffered_output_file_close(buffered_output); \
  gt_input_generic_parser_attributes_delete(__gparser_attr); \
  gt_template_delete(template)

/*
 * File basic I/O
 */
typedef enum { GT_COMPRESSION_NONE, GT_COMPRESSION_BZIP, GT_COMPRESSION_GZIP } gt_compression_t;

GT_INLINE gt_input_file* gt_tools_open_input_file(char* const name_input_file,const gt_compression_t compression_type);
GT_INLINE gt_output_file* gt_tools_open_output_file(char* const name_output_file,const gt_compression_t compression_type);

/*
 * Opening an archive (GEMIndex/MULTIFastaFile)
 */
GT_INLINE gt_sequence_archive* gt_tools_open_sequence_archive(
    char* const name_gem_index_file,char* const name_reference_file,const bool load_sequences);

/*
 * Argument Parsing
 */
GT_INLINE void gt_tools_parse_cvs_uint64(char* const parameters_list,const uint64_t num_params,...);
GT_INLINE uint64_t gt_tools_parse_cvs_float(char* const parameters_list,const uint64_t num_params,...);

#endif /* GEM_TOOLS_H_ */
