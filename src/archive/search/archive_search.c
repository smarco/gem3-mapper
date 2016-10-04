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
    [3] = "recovery",
    [4] = "find-pairs",
    [5] = "end",
};

/*
 * Constants
 */
#define BISULFITE_SEQUENCE_INITIAL_LENGTH 200

/*
 * Setup
 */
void archive_search_init(
    archive_search_t* const archive_search,
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    const bool buffered_search,
    mm_stack_t* const mm_stack) {
  // Archive
  archive_search->archive = archive;
  // Approximate Search
  memcpy(&archive_search->search_parameters,search_parameters,sizeof(search_parameters_t));
  if (mm_stack==NULL) {
    // Sequence
    sequence_init(&archive_search->sequence);
    // BS-Sequence
    if (archive->type == archive_dna_bisulfite) {
      archive_search->bs_original_sequence = mm_alloc(string_t);
      string_init(archive_search->bs_original_sequence,BISULFITE_SEQUENCE_INITIAL_LENGTH);
    } else {
      archive_search->bs_original_sequence = NULL;
    }
  } else {
    // Sequence
    sequence_init_mm(&archive_search->sequence,mm_stack);
    // BS-Sequence
    if (archive->type == archive_dna_bisulfite) {
      archive_search->bs_original_sequence = mm_alloc(string_t);
      string_init_mm(archive_search->bs_original_sequence,BISULFITE_SEQUENCE_INITIAL_LENGTH,mm_stack);
    } else {
      archive_search->bs_original_sequence = NULL;
    }
  }
  // BS Init
  const bisulfite_read_t bisulfite_read_mode = search_parameters->bisulfite_read;
  archive_search->bs_sequence_end = (bisulfite_read_mode==bisulfite_read_2) ? paired_end2 : paired_end1;
  // Approximate Search Init
  approximate_search_init(&archive_search->approximate_search,archive,&archive_search->search_parameters);
  archive_search->buffered_search = buffered_search;
}
void archive_search_se_new(
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    const bool buffered_search,
    mm_stack_t* const mm_stack,
    archive_search_t** const archive_search) {
  // Prepare Search
  *archive_search = mm_alloc(archive_search_t); // Allocate handler
  archive_search_init(*archive_search,
      archive,search_parameters,buffered_search,mm_stack);
  // Select align
  archive_select_configure_se(*archive_search);
}
void archive_search_pe_new(
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    const bool buffered_search,
    mm_stack_t* const mm_stack,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2) {
  // Allocate Search
  *archive_search_end1 = mm_alloc(archive_search_t); // Allocate handler
  archive_search_init(*archive_search_end1,
      archive,search_parameters,buffered_search,mm_stack);
  *archive_search_end2 = mm_alloc(archive_search_t); // Allocate handler
  archive_search_init(*archive_search_end2,
      archive,search_parameters,buffered_search,mm_stack);
  // Select align
  archive_select_configure_pe(*archive_search_end1);
  archive_select_configure_pe(*archive_search_end2);
}
void archive_search_prepare_sequence(archive_search_t* const archive_search) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_PREPARE_SEQUENCE,PROFILE_LEVEL);
  // Generate the pattern(s)
  const bool run_length_pattern = archive_search->archive->text->run_length;
  const bool kmer_filter_compile = !archive_search->buffered_search;
  approximate_search_t* const approximate_search = &archive_search->approximate_search;
  pattern_init(&approximate_search->pattern,&archive_search->sequence,
      &approximate_search->do_quality_search,approximate_search->search_parameters,
      run_length_pattern,kmer_filter_compile,archive_search->mm_stack);
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_PREPARE_SEQUENCE,PROFILE_LEVEL);
}
void archive_search_reset(archive_search_t* const archive_search) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_INIT,PROFILE_LEVEL);
  // Instantiate parameters actual-values
  const uint64_t sequence_length = sequence_get_length(&archive_search->sequence);
  search_instantiate_values(&archive_search->search_parameters,sequence_length);
  // Prepare for sequence
  archive_search_prepare_sequence(archive_search);
  // Clear F/R search states
  approximate_search_reset(&archive_search->approximate_search);
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_INIT,PROFILE_LEVEL);
}
void archive_search_destroy(archive_search_t* const archive_search) {
  // Destroy Sequence
  sequence_destroy(&archive_search->sequence);
  if (archive_search->bs_original_sequence!=NULL) {
    string_destroy(archive_search->bs_original_sequence);
    mm_free(archive_search->bs_original_sequence);
  }
}
void archive_search_delete(archive_search_t* const archive_search) {
  // Destroy archive-search
  archive_search_destroy(archive_search);
  // Free handler
  mm_free(archive_search);
}
/*
 * Accessors
 */
sequence_t* archive_search_get_sequence(const archive_search_t* const archive_search) {
  return (sequence_t*)&archive_search->sequence;
}
bool archive_search_finished(const archive_search_t* const archive_search) {
  return archive_search->approximate_search.search_stage == asearch_stage_end;
}
uint64_t archive_search_get_num_regions_profile(const archive_search_t* const archive_search) {
  return approximate_search_get_num_regions_profile(&archive_search->approximate_search);
}
uint64_t archive_search_get_num_decode_candidates(const archive_search_t* const archive_search) {
  return approximate_search_get_num_decode_candidates(&archive_search->approximate_search);
}
uint64_t archive_search_get_num_verify_candidates(const archive_search_t* const archive_search) {
  return approximate_search_get_num_verify_candidates(&archive_search->approximate_search);
}
