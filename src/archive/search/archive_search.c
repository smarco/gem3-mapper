/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Archive-Search main module & support data structures
 */

#include "archive/search/archive_search.h"
#include "archive/search/archive_select.h"
#include "matches/classify/matches_classify.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * Archive Search State
 */
const char* archive_search_pe_state_label[] =
{
    [0] = "begin",
    [1] = "search-end1",
    [2] = "search-end2",
    [3] = "find-pairs",
    [4] = "end",
};

/*
 * Setup
 */
void archive_search_configure(
    archive_search_t* const archive_search,
    search_parameters_t* const search_parameters,
    const bool buffered_search) {
  // Configure search-parameters
  memcpy(&archive_search->search_parameters,search_parameters,sizeof(search_parameters_t));
  archive_search->buffered_search = buffered_search;
  // Sequence
  archive_search->sequence = NULL;
}
void archive_search_se_new(
    search_parameters_t* const search_parameters,
    const bool buffered_search,
    archive_search_t** const archive_search) {
  // Prepare Search
  *archive_search = mm_alloc(archive_search_t); // Allocate handler
  archive_search_configure(*archive_search,search_parameters,buffered_search);
  archive_select_configure_se(*archive_search);
}
void archive_search_pe_new(
    search_parameters_t* const search_parameters,
    const bool buffered_search,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2) {
  // Allocate Search
  *archive_search_end1 = mm_alloc(archive_search_t); // Allocate handler
  archive_search_configure(*archive_search_end1,search_parameters,buffered_search);
  *archive_search_end2 = mm_alloc(archive_search_t); // Allocate handler
  archive_search_configure(*archive_search_end2,search_parameters,buffered_search);
  // Configure Select align
  archive_select_configure_pe(*archive_search_end1);
  archive_select_configure_pe(*archive_search_end2);
}
void archive_search_destroy(
    archive_search_t* const archive_search) {
  // Clear Approximate Search
  approximate_search_destroy(&archive_search->approximate_search);
  // Clear Sequence
  sequence_destroy(archive_search->sequence);
  mm_allocator_free(archive_search->mm_allocator,archive_search->sequence);
  archive_search->sequence = NULL;
}
void archive_search_delete(
    archive_search_t* const archive_search) {
  // Free handler
  mm_free(archive_search);
}
/*
 * Prepare Search
 */
void archive_search_prepare_sequence(
    archive_search_t* const archive_search,
    sequence_t* const sequence) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_INIT,PROFILE_LEVEL);
  // Instantiate parameters actual-values
  const uint64_t sequence_length = sequence_get_length(sequence);
  search_instantiate_values(&archive_search->search_parameters,sequence_length);
  archive_search->sequence = sequence;
  // Prepare the pattern
  const bool run_length_pattern = archive_search->archive->text->run_length;
  approximate_search_prepare(
      &archive_search->approximate_search,
      run_length_pattern,sequence);
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_INIT,PROFILE_LEVEL);
}
void archive_search_inject_handlers(
    archive_search_t* const archive_search,
    archive_t* const archive,
    filtering_candidates_t* const filtering_candidates,
    nsearch_schedule_t* const nsearch_schedule,
    mapper_stats_t* const mapper_stats,
    mm_allocator_t* const mm_allocator) {
  // Handlers Injection (Support Data Structures)
  archive_search->archive = archive; // Archive
  approximate_search_init(
      &archive_search->approximate_search,
      archive,&archive_search->search_parameters,
      filtering_candidates,nsearch_schedule,mm_allocator); // Approximate Search
  // Stats
  archive_search->mapper_stats = mapper_stats;
  // MM
  archive_search->mm_allocator = mm_allocator;
}
/*
 * Accessors
 */
bool archive_search_finished(const archive_search_t* const archive_search) {
  return archive_search->approximate_search.search_stage == asearch_stage_end;
}
uint64_t archive_search_get_num_regions_profile(const archive_search_t* const archive_search) {
  return approximate_search_get_num_regions_profile(&archive_search->approximate_search);
}
uint64_t archive_search_get_num_decode_candidates(const archive_search_t* const archive_search) {
  return approximate_search_get_num_decode_candidates(&archive_search->approximate_search);
}
uint64_t archive_search_get_num_kmer_filter_candidates(const archive_search_t* const archive_search) {
  return approximate_search_get_num_filtering_candidates(&archive_search->approximate_search);
}
uint64_t archive_search_get_num_bpm_distance_candidates(const archive_search_t* const archive_search) {
  return approximate_search_get_num_filtering_candidates(&archive_search->approximate_search);
}
uint64_t archive_search_get_num_bpm_align_candidates(const archive_search_t* const archive_search) {
  return approximate_search_get_num_filtering_candidates_buffered(&archive_search->approximate_search);
}
uint64_t archive_search_get_num_bpm_align_canonical_candidates(const archive_search_t* const archive_search) {
  return approximate_search_get_num_filtering_canonical_candidates_buffered(&archive_search->approximate_search);
}
uint64_t archive_search_get_num_bpm_align_candidate_tiles_length(const archive_search_t* const archive_search) {
  return approximate_search_get_num_filtering_candidate_buffered_tiles_length(&archive_search->approximate_search);
}
