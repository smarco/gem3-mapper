/*
 * PROJECT: GEMMapper
 * FILE: search_stage_state.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef SEARCH_STAGE_STATE_H_
#define SEARCH_STAGE_STATE_H_

#include "utils/essentials.h"
#include "archive/archive_search.h"
#include "archive/archive_search_se_stepwise.h"

/*
 * Buffer State
 */
typedef enum {
  search_group_buffer_phase_sending,
  search_group_buffer_phase_retrieving
} search_stage_mode_t;
typedef struct {
  uint64_t num_buffers;           // Total buffers
  uint64_t current_buffer_idx;    // Current buffer index
  uint64_t num_searches;          // Total archive-searches (in current buffer)
  uint64_t current_search_idx;    // Current search index (in current buffer)
} search_stage_iterator_t;

#endif /* SEARCH_STAGE_STATE_H_ */
