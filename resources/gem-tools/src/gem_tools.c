/*
 * PROJECT: GEM-Tools library
 * FILE: gem_tools.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gem_tools.h"

/*
 * File I/O basics
 */
gt_input_file* gt_tools_open_input_file(char* const name_input_file,const gt_compression_t compression_type) {
  switch (compression_type) {
    case GT_COMPRESSION_NONE:
      return ((name_input_file==NULL) ?
          gt_input_stream_open(stdin) :
          gt_input_file_open(name_input_file,false));
      break;
    case GT_COMPRESSION_GZIP:
      return ((name_input_file==NULL) ?
          gt_input_gzip_stream_open(stdin) :
          gt_input_gzip_stream_open(gt_open_FILE(name_input_file,"r")));
      break;
    case GT_COMPRESSION_BZIP:
      return ((name_input_file==NULL) ?
          gt_input_bzip_stream_open(stdin) :
          gt_input_bzip_stream_open(gt_open_FILE(name_input_file,"r")));
      break;
    default:
      GT_INVALID_CASE();
      break;
  }
}
gt_output_file* gt_tools_open_output_file(char* const name_output_file,const gt_compression_t compression_type) {
  switch (compression_type) {
    case GT_COMPRESSION_NONE:
      return ((name_output_file==NULL) ?
          gt_output_stream_new(stdout,SORTED_FILE) :
          gt_output_file_new(name_output_file,SORTED_FILE));
      break;
    case GT_COMPRESSION_GZIP:
      return ((name_output_file==NULL) ?
          gt_output_gzip_stream_new(stdout,SORTED_FILE) :
          gt_output_gzip_stream_new(gt_open_FILE(name_output_file,"w"),SORTED_FILE));
      break;
    case GT_COMPRESSION_BZIP:
      return ((name_output_file==NULL) ?
          gt_output_bzip_stream_new(stdout,SORTED_FILE) :
          gt_output_bzip_stream_new(gt_open_FILE(name_output_file,"w"),SORTED_FILE));
      break;
    default:
      GT_INVALID_CASE();
      break;
  }
}
/*
 * Opening an archive (GEMIndex/MULTIFastaFile)
 */
GT_INLINE gt_sequence_archive* gt_tools_open_sequence_archive(
    char* const name_gem_index_file,char* const name_reference_file,const bool load_sequences) {
  gt_sequence_archive* sequence_archive = NULL;
  gt_log("Loading reference file ...");
  if (name_gem_index_file!=NULL) { // Load GEM-IDX
    sequence_archive = gt_sequence_archive_new(GT_BED_ARCHIVE);
    gt_gemIdx_load_archive(name_gem_index_file,sequence_archive,load_sequences);
  } else {
    gt_input_file* const reference_file = gt_input_file_open(name_reference_file,false);
    sequence_archive = gt_sequence_archive_new(GT_CDNA_ARCHIVE);
    if (gt_input_multifasta_parser_get_archive(reference_file,sequence_archive)!=GT_IFP_OK) {
      gt_fatal_error_msg("Error parsing reference file '%s'\n",name_reference_file);
    }
    gt_input_file_close(reference_file);
  }
  gt_log("Done.");
  return sequence_archive;
}
/*
 * Argument Parsing
 */
GT_INLINE void gt_tools_parse_cvs_uint64(char* const parameters_list,const uint64_t num_params,...) {
  uint64_t num_params_parsed = 0;
  // Start va_args
  va_list v_args;
  va_start(v_args,num_params);
  // Start parsing
  char *opt = strtok(parameters_list,",");
  while (opt!=NULL && num_params_parsed<num_params) {
    uint64_t* const uint64_arg = va_arg(v_args,uint64_t*);
    *uint64_arg = atoll(opt);
    opt = strtok(NULL,",");
  }
  // End va_args
  va_end(v_args);
}
GT_INLINE uint64_t gt_tools_parse_cvs_float(char* const parameters_list,const uint64_t num_params,...) {
  uint64_t num_params_parsed = 0;
  // Start va_args
  va_list v_args;
  va_start(v_args,num_params);
  // Start parsing
  char *opt = strtok(parameters_list,",");
  while (opt!=NULL && num_params_parsed<num_params) {
    float* const float_arg = va_arg(v_args,float*);
    *float_arg = atof(opt);
    opt = strtok(NULL,",");
    ++num_params_parsed;
  }
  // End va_args
  va_end(v_args);
  return num_params_parsed;
}

