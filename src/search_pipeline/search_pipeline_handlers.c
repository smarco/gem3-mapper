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
 */

#include "search_pipeline/search_pipeline_handlers.h"

/*
 * Setup
 */
search_pipeline_handlers_t* search_pipeline_handlers_new(archive_t* const archive) {
  // Allocate
  search_pipeline_handlers_t* const search_pipeline_handlers = mm_alloc(search_pipeline_handlers_t);
  // Archive
  search_pipeline_handlers->archive = archive;
  // Filtering Candidates
  filtering_candidates_init(&search_pipeline_handlers->fc_decode_end1);
  filtering_candidates_init(&search_pipeline_handlers->fc_decode_end2);
  filtering_candidates_init(&search_pipeline_handlers->fc_kmer_filter_end1);
  filtering_candidates_init(&search_pipeline_handlers->fc_kmer_filter_end2);
  filtering_candidates_init(&search_pipeline_handlers->fc_bpm_distance_end1);
  filtering_candidates_init(&search_pipeline_handlers->fc_bpm_distance_end2);
  filtering_candidates_init(&search_pipeline_handlers->fc_bpm_align_end1);
  filtering_candidates_init(&search_pipeline_handlers->fc_bpm_align_end2);
  // Stats
  search_pipeline_handlers->mapper_stats = mapper_stats_new();
  // MM
  search_pipeline_handlers->mm_slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_16M,MM_UNLIMITED_MEM);
  search_pipeline_handlers->mm_allocator = mm_allocator_new(search_pipeline_handlers->mm_slab);
  // Return
  return search_pipeline_handlers;
}
void search_pipeline_handlers_clear(search_pipeline_handlers_t* const search_pipeline_handlers) {
  filtering_candidates_clear(&search_pipeline_handlers->fc_decode_end1,true);
  filtering_candidates_clear(&search_pipeline_handlers->fc_decode_end2,true);
  filtering_candidates_clear(&search_pipeline_handlers->fc_kmer_filter_end1,true);
  filtering_candidates_clear(&search_pipeline_handlers->fc_kmer_filter_end2,true);
  filtering_candidates_clear(&search_pipeline_handlers->fc_bpm_distance_end1,true);
  filtering_candidates_clear(&search_pipeline_handlers->fc_bpm_distance_end2,true);
  filtering_candidates_clear(&search_pipeline_handlers->fc_bpm_align_end1,true);
  filtering_candidates_clear(&search_pipeline_handlers->fc_bpm_align_end2,true);
  mm_allocator_clear(search_pipeline_handlers->mm_allocator);
}
void search_pipeline_handlers_delete(search_pipeline_handlers_t* const search_pipeline_handlers) {
  filtering_candidates_destroy(&search_pipeline_handlers->fc_decode_end1,true);
  filtering_candidates_destroy(&search_pipeline_handlers->fc_decode_end2,true);
  filtering_candidates_destroy(&search_pipeline_handlers->fc_kmer_filter_end1,true);
  filtering_candidates_destroy(&search_pipeline_handlers->fc_kmer_filter_end2,true);
  filtering_candidates_destroy(&search_pipeline_handlers->fc_bpm_distance_end1,true);
  filtering_candidates_destroy(&search_pipeline_handlers->fc_bpm_distance_end2,true);
  filtering_candidates_destroy(&search_pipeline_handlers->fc_bpm_align_end1,true);
  filtering_candidates_destroy(&search_pipeline_handlers->fc_bpm_align_end2,true);
  mapper_stats_delete(search_pipeline_handlers->mapper_stats);
  mm_allocator_delete(search_pipeline_handlers->mm_allocator);
  mm_slab_delete(search_pipeline_handlers->mm_slab);
  mm_free(search_pipeline_handlers);
}
/*
 * Injection (Support Data Structures)
 */
void search_pipeline_handlers_prepare_se(
    archive_search_t* const archive_search,
    sequence_t* const sequence,
	 bisulfite_conversion_t const bisulfite_conversion,
    search_pipeline_handlers_t* const search_pipeline_handlers) {
  // Inject Handlers
  archive_search_inject_handlers(
      archive_search,search_pipeline_handlers->archive,
      NULL,&search_pipeline_handlers->nsearch_schedule,
      search_pipeline_handlers->mapper_stats,
      search_pipeline_handlers->mm_allocator);
  // Set bisulfite conversion
  archive_search->approximate_search.bisulfite_conversion=bisulfite_conversion;
  // Prepare sequence
  archive_search_prepare_sequence(archive_search,sequence);
}
void search_pipeline_handlers_prepare_pe(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    sequence_t* const sequence_end1,
    sequence_t* const sequence_end2,
	 bisulfite_conversion_t const bisulfite_conversion_end1,
    bisulfite_conversion_t const bisulfite_conversion_end2,
    search_pipeline_handlers_t* const search_pipeline_handlers) {
  // Inject Handlers
  archive_search_inject_handlers(
      archive_search_end1,search_pipeline_handlers->archive,
      NULL,&search_pipeline_handlers->nsearch_schedule,
      search_pipeline_handlers->mapper_stats,
      search_pipeline_handlers->mm_allocator);
  archive_search_inject_handlers(
      archive_search_end2,search_pipeline_handlers->archive,
      NULL,&search_pipeline_handlers->nsearch_schedule,
      NULL,search_pipeline_handlers->mm_allocator);
  // Set bisulfite conversion
  archive_search_end1->approximate_search.bisulfite_conversion=bisulfite_conversion_end1;
  archive_search_end2->approximate_search.bisulfite_conversion=bisulfite_conversion_end2;

  // Prepare sequences
  archive_search_prepare_sequence(archive_search_end1,sequence_end1);
  archive_search_prepare_sequence(archive_search_end2,sequence_end2);
}
