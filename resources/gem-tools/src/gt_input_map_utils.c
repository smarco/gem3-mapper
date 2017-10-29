/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_map_utils.c
 * DATE: 01/03/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_input_map_utils.h"

#include "gt_input_file.h"
#include "gt_buffered_input_file.h"
#include "gt_input_parser.h"
#include "gt_input_map_parser.h"

#include "gt_output_buffer.h"
#include "gt_buffered_output_file.h"
#include "gt_output_map.h"

#include "gt_alignment.h"
#include "gt_alignment_utils.h"
#include "gt_template.h"
#include "gt_template_utils.h"

/*
 * Parsers Helpers
 */
GT_INLINE gt_status gt_merge_map_read_template_sync(
    pthread_mutex_t* const input_mutex,gt_buffered_input_file* const buffered_input_master,gt_buffered_input_file* const buffered_input_slave,
    gt_map_parser_attributes* const map_parser_attr,gt_template* const template_master,gt_template* const template_slave,
    gt_buffered_output_file* buffered_output) {
  gt_status error_code_master, error_code_slave;
  gt_output_map_attributes output_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();
  do {
    // Read Synch blocks
    error_code_master=gt_input_map_parser_synch_blocks_by_subset(
        input_mutex,map_parser_attr,buffered_input_master,buffered_input_slave);
    if (error_code_master==GT_INPUT_STATUS_EOF) return GT_INPUT_STATUS_EOF;
    if (error_code_master==GT_INPUT_STATUS_FAIL) gt_fatal_error_msg("Fatal error synchronizing files");
    // Read master (always guaranteed)
    if ((error_code_master=gt_input_map_parser_get_template(buffered_input_master,template_master,NULL))==GT_INPUT_STATUS_FAIL) {
      gt_fatal_error_msg("Fatal error parsing file <<Master>>");
    }
    // Check slave
    if (gt_buffered_input_file_eob(buffered_input_slave)) { // Slave exhausted. Dump master & return EOF
      gt_output_map_bofprint_template(buffered_output,template_master,&output_attributes);
      while ((error_code_master=gt_input_map_parser_get_template(buffered_input_master,template_master,NULL))) {
        if (error_code_master==GT_INPUT_STATUS_FAIL) gt_fatal_error_msg("Fatal error parsing file <<Master>>");
        gt_output_map_bofprint_template(buffered_output,template_master,&output_attributes);
      }
    } else {
      // Read slave
      if ((error_code_slave=gt_input_map_parser_get_template(buffered_input_slave,template_slave,NULL))==GT_INPUT_STATUS_FAIL) {
        gt_fatal_error_msg("Fatal error parsing file <<Slave>>");
      }
      // Synch loop
      while (!gt_streq(gt_template_get_tag(template_master),gt_template_get_tag(template_slave))) {
        // Print non correlative master's template
        gt_output_map_bofprint_template(buffered_output,template_master,&output_attributes);
        // Fetch next master's template
        if (gt_buffered_input_file_eob(buffered_input_master)) gt_fatal_error_msg("<<Slave>> contains more/different reads from <<Master>>");
        if ((error_code_master=gt_input_map_parser_get_template(buffered_input_master,template_master,NULL))!=GT_INPUT_STATUS_OK) {
          gt_fatal_error_msg("Fatal error parsing file <<Master>>");
        }
      }
      return GT_INPUT_STATUS_OK;
    }
  } while (true);
}

GT_INLINE void gt_merge_synch_map_files_(
    pthread_mutex_t* const input_mutex,const bool paired_end,gt_buffered_output_file* const buffered_output_file,
    gt_buffered_input_file** const buffered_input_file,const uint64_t num_files) {
  GT_NULL_CHECK(input_mutex);
  GT_BUFFERED_OUTPUT_FILE_CHECK(buffered_output_file);
  GT_NULL_CHECK(buffered_input_file);
  GT_ZERO_CHECK(num_files);
  // Attach the writing to the first buffered_input_file
  gt_buffered_input_file_attach_buffered_output(buffered_input_file[0],buffered_output_file);
  // Init templates
  uint64_t i;
  gt_template** template = gt_calloc(num_files,gt_template*,false);
  for (i=0;i<num_files;++i) {
    template[i] = gt_template_new(); // Allocate template
  }
  // Merge loop
  gt_map_parser_attributes map_parser_attr = GT_MAP_PARSER_ATTR_DEFAULT(paired_end);
  gt_output_map_attributes output_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();
  gt_status error_code_master, error_code_slave;
  while (gt_input_map_parser_synch_blocks_a(input_mutex,buffered_input_file,num_files,&map_parser_attr)) {
    // Read master (always guaranteed)
    if ((error_code_master=gt_input_map_parser_get_template(buffered_input_file[0],template[0],NULL))==GT_INPUT_STATUS_FAIL) {
      gt_fatal_error_msg("Fatal error parsing file Master::%s",
          gt_input_file_get_file_name(buffered_input_file[0]->input_file));
    }
    for (i=1;i<num_files;++i) {
      // Read slave
      if ((error_code_slave=gt_input_map_parser_get_template(buffered_input_file[i],template[i],NULL))==GT_INPUT_STATUS_FAIL) {
        gt_fatal_error_msg("Fatal error parsing file Slave::%s",
            gt_input_file_get_file_name(buffered_input_file[i]->input_file));
      }
      if (error_code_master!=error_code_slave) {
        gt_fatal_error_msg("<<Slave>> contains more/different reads from <<Master>>, ('%s','%s')",
            gt_input_file_get_file_name(buffered_input_file[0]->input_file),
            gt_input_file_get_file_name(buffered_input_file[i]->input_file));
      }
      // Check tags
      if (!gt_string_equals(template[0]->tag,template[i]->tag)) {
        gt_fatal_error_msg("Files are not synchronized ('%s','%s')\n"
            "\tDifferent TAGs found '"PRIgts"' '"PRIgts"' ",
            gt_input_file_get_file_name(buffered_input_file[0]->input_file),
            gt_input_file_get_file_name(buffered_input_file[i]->input_file),
            PRIgts_content(template[0]->tag),PRIgts_content(template[i]->tag));
      }
    }
    // Merge maps
    gt_template *ptemplate = gt_template_union_template_mmaps_a(template,num_files);
    gt_output_map_bofprint_template(buffered_output_file,ptemplate,&output_attributes); // Print template
    gt_template_delete(ptemplate); // Delete template
  }
  // Free
  for (i=0;i<num_files;++i) {
    gt_template_delete(template[i]);
  }
  gt_free(template);
}

GT_INLINE void gt_merge_synch_map_files_a(
    pthread_mutex_t* const input_mutex,const bool paired_end,gt_output_file* const output_file,
    gt_input_file** const input_map_files,const uint64_t num_files) {
  GT_NULL_CHECK(input_mutex);
  GT_OUTPUT_FILE_CHECK(output_file);
  GT_NULL_CHECK(input_map_files);
  GT_ZERO_CHECK(num_files);
  // Init Buffered Output File
  gt_buffered_output_file* const buffered_output_file = gt_buffered_output_file_new(output_file);
  // Init Buffered Input Files
  uint64_t i;
  gt_buffered_input_file** buffered_input_file = gt_calloc(num_files,gt_buffered_input_file*,false);
  for (i=0;i<num_files;++i) {
    GT_INPUT_FILE_CHECK(input_map_files[i]);
    // Open buffered input file
    buffered_input_file[i] = gt_buffered_input_file_new(input_map_files[i]);
  }
  // Buffered Reading+process
  gt_merge_synch_map_files_(input_mutex,paired_end,buffered_output_file,buffered_input_file,num_files);
  // Clean & free
  for (i=0;i<num_files;++i) {
    gt_buffered_input_file_close(buffered_input_file[i]);
  }
  gt_free(buffered_input_file);
  gt_buffered_output_file_close(buffered_output_file);
}

GT_INLINE void gt_merge_synch_map_files_v(
    pthread_mutex_t* const input_mutex,const bool paired_end,gt_output_file* const output_file,
    gt_input_file* const input_map_master,const uint64_t num_slaves,va_list v_args) {
  GT_NULL_CHECK(input_mutex);
  GT_OUTPUT_FILE_CHECK(output_file);
  GT_INPUT_FILE_CHECK(input_map_master);
  // Setup master
  gt_buffered_output_file* const buffered_output_file = gt_buffered_output_file_new(output_file);
  // Setup slaves
  uint64_t i;
  const uint64_t num_files = num_slaves+1;
  gt_buffered_input_file** buffered_input_file = gt_calloc(num_files,gt_buffered_input_file*,false);
  // Init master (i=0)
  GT_INPUT_FILE_CHECK(input_map_master);
  buffered_input_file[0] = gt_buffered_input_file_new(input_map_master);
  // Init slaves
  for (i=1;i<num_files;++i) {
    // Get input file slave
    gt_input_file* input_map_file = va_arg(v_args,gt_input_file*);
    GT_INPUT_FILE_CHECK(input_map_file);
    // Open buffered input slave
    buffered_input_file[i] = gt_buffered_input_file_new(input_map_file);
  }
  // Buffered Reading+process
  gt_merge_synch_map_files_(input_mutex,paired_end,buffered_output_file,buffered_input_file,num_files);
  // Clean
  for (i=0;i<num_files;++i) {
    gt_buffered_input_file_close(buffered_input_file[i]);
  }
  gt_free(buffered_input_file);
  gt_buffered_output_file_close(buffered_output_file);
}

GT_INLINE void gt_merge_synch_map_files_va(
    pthread_mutex_t* const input_mutex,const bool paired_end,gt_output_file* const output_file,
    gt_input_file* const input_map_master,const uint64_t num_slaves,...) {
  GT_NULL_CHECK(input_mutex);
  GT_OUTPUT_FILE_CHECK(output_file);
  GT_INPUT_FILE_CHECK(input_map_master);
  va_list v_args;
  va_start(v_args,num_slaves);
  gt_merge_synch_map_files_v(input_mutex,paired_end,output_file,input_map_master,num_slaves,v_args);
  va_end(v_args);
}

GT_INLINE void gt_merge_unsynch_map_files(
    pthread_mutex_t* const input_mutex,gt_input_file* const input_map_master,gt_input_file* const input_map_slave,
    const bool paired_end,gt_output_file* const output_file) {
  // Reading+process
  gt_buffered_input_file* buffered_input_master = gt_buffered_input_file_new(input_map_master);
  gt_buffered_input_file* buffered_input_slave  = gt_buffered_input_file_new(input_map_slave);
  gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
  gt_buffered_input_file_attach_buffered_output(buffered_input_master,buffered_output);

  gt_template *template_master = gt_template_new();
  gt_template *template_slave = gt_template_new();

  gt_map_parser_attributes map_parser_attr = GT_MAP_PARSER_ATTR_DEFAULT(paired_end);
  gt_output_map_attributes output_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();

  while (gt_merge_map_read_template_sync(input_mutex,buffered_input_master,buffered_input_slave,
      &map_parser_attr,template_master,template_slave,buffered_output)) {
    // Merge maps
    // gt_template *ptemplate = gt_template_union_template_mmaps_fx_va(gt_mmap_cmp,gt_map_cmp,2,template_master,template_slave);
    gt_template *ptemplate = gt_template_union_template_mmaps(template_master,template_slave);
    // Print template
    gt_output_map_bofprint_template(buffered_output,ptemplate,&output_attributes);
    // Delete template
    gt_template_delete(ptemplate);
  }

  // Clean
  gt_template_delete(template_master);
  gt_template_delete(template_slave);
  gt_buffered_input_file_close(buffered_input_master);
  gt_buffered_input_file_close(buffered_input_slave);
  gt_buffered_output_file_close(buffered_output);
}
