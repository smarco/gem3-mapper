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

#ifndef SEARCH_PIPELINE_HANDLERS_H_
#define SEARCH_PIPELINE_HANDLERS_H_

#include "utils/essentials.h"
#include "archive/archive.h"
#include "archive/search/archive_search.h"
#include "archive/search/archive_search_se_parameters.h"
#include "mapper/mapper_stats.h"

/*
 * Search Pipeline Handlers
 */
typedef struct {
  /* Archive */
  archive_t* archive;                   // Archive
  /* Mapping Statistics */
  mapper_stats_t* mapper_stats;         // Mapping Statistics
  /* Filtering Candidates */
  filtering_candidates_t fc_decode_end1;
  filtering_candidates_t fc_decode_end2;
  filtering_candidates_t fc_kmer_filter_end1;
  filtering_candidates_t fc_kmer_filter_end2;
  filtering_candidates_t fc_bpm_distance_end1;
  filtering_candidates_t fc_bpm_distance_end2;
  filtering_candidates_t fc_bpm_align_end1;
  filtering_candidates_t fc_bpm_align_end2;
  /* Neighborhood Search */
  nsearch_schedule_t nsearch_schedule;
  /* MM */
  mm_slab_t* mm_slab;                   // MM-Slab
  mm_allocator_t* mm_allocator;         // MM-Allocator
} search_pipeline_handlers_t;

/*
 * Setup
 */
search_pipeline_handlers_t* search_pipeline_handlers_new(archive_t* const archive);
void search_pipeline_handlers_clear(search_pipeline_handlers_t* const search_pipeline_handlers);
void search_pipeline_handlers_delete(search_pipeline_handlers_t* const search_pipeline_handlers);

/*
 * Injection (Support Data Structures)
 */
void search_pipeline_handlers_prepare_se(
    archive_search_t* const archive_search,
    sequence_t* const sequence,
    search_pipeline_handlers_t* const search_pipeline_handlers);
void search_pipeline_handlers_prepare_pe(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    sequence_t* const sequence_end1,
    sequence_t* const sequence_end2,
    search_pipeline_handlers_t* const search_pipeline_handlers);

#endif /* SEARCH_PIPELINE_HANDLERS_H_ */
