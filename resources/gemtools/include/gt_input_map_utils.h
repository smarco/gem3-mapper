/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_map_utils.h
 * DATE: 01/03/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_INPUT_MAP_UTILS_H_
#define GT_INPUT_MAP_UTILS_H_

#include "gt_essentials.h"

#include "gt_input_file.h"
#include "gt_output_file.h"

// Merge functions (synch files)
#define gt_merge_synch_map_files(input_mutex,paired_end,output_file,input_map_master,input_map_slave) \
  gt_merge_synch_map_files_va(input_mutex,paired_end,output_file,input_map_master,1,input_map_slave)
GT_INLINE void gt_merge_synch_map_files_a(
    pthread_mutex_t* const input_mutex,const bool paired_end,gt_output_file* const output_file,
    gt_input_file** const input_map_files,const uint64_t num_input_map_files);
GT_INLINE void gt_merge_synch_map_files_v(
    pthread_mutex_t* const input_mutex,const bool paired_end,gt_output_file* const output_file,
    gt_input_file* const input_map_master,const uint64_t num_slaves,va_list v_args);
GT_INLINE void gt_merge_synch_map_files_va(
    pthread_mutex_t* const input_mutex,const bool paired_end,gt_output_file* const output_file,
    gt_input_file* const input_map_master,const uint64_t num_slaves,...);

// Merge functions (unsynch files)
GT_INLINE void gt_merge_unsynch_map_files(
    pthread_mutex_t* const input_mutex,gt_input_file* const input_map_master,gt_input_file* const input_map_slave,
    const bool paired_end,gt_output_file* const output_file);

#endif /* GT_INPUT_MAP_UTILS_H_ */
