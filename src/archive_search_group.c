/*
 * PROJECT: GEMMapper
 * FILE: archive_search_group.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_search_group.h"
#include "bpm_align_gpu.h"

/*
 * Archive-search group
 */
GEM_INLINE void archive_search_group_add(
    archive_search_group_t* const archive_search_group,
    archive_search_t* const archive_search,const uint64_t results_buffer_offset) {
  archive_search_member_t* archive_search_member;
  vector_alloc_new(archive_search_group->archive_searches,archive_search_member_t,archive_search_member);
  archive_search_member->archive_search = archive_search;
  archive_search_member->results_buffer_offset = results_buffer_offset;
}
/*
 * Dispatcher
 */
GEM_INLINE archive_search_group_dispatcher_t* archive_search_group_dispatcher_new() {
  GEM_NOT_IMPLEMENTED(); // TODO
}
GEM_INLINE void archive_search_group_dispatcher_delete(archive_search_group_dispatcher_t* const dispatcher) {
  GEM_NOT_IMPLEMENTED(); // TODO
}

GEM_INLINE archive_search_group_t* archive_search_group_dispatcher_request_generating(
    archive_search_group_dispatcher_t* const dispatcher) {
  GEM_NOT_IMPLEMENTED(); // TODO
}
GEM_INLINE void archive_search_group_dispatcher_return_generating(
    archive_search_group_dispatcher_t* const dispatcher,archive_search_group_t* const search_group) {
  GEM_NOT_IMPLEMENTED(); // TODO
}
GEM_INLINE archive_search_group_t* archive_search_group_dispatcher_request_selecting(
    archive_search_group_dispatcher_t* const dispatcher) {
  GEM_NOT_IMPLEMENTED(); // TODO
}
GEM_INLINE void archive_search_group_dispatcher_return_selecting(
    archive_search_group_dispatcher_t* const dispatcher,archive_search_group_t* const search_group) {
  GEM_NOT_IMPLEMENTED(); // TODO
}
/*
 * Step-wise SE-Search
 */
GEM_INLINE void archive_search_generate_candidates(archive_search_t* const archive_search) {
  GEM_NOT_IMPLEMENTED(); // TODO
}
GEM_INLINE void archive_search_copy_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  GEM_NOT_IMPLEMENTED(); // TODO
}
GEM_INLINE void archive_search_select_candidates(
    archive_search_t* const archive_search,
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const uint64_t results_buffer_offset) {
  GEM_NOT_IMPLEMENTED(); // TODO
}
