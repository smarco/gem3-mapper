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

#ifndef SEARCH_PIPELINE_H_
#define SEARCH_PIPELINE_H_

#include "utils/essentials.h"
#include "search_pipeline/search_stage_region_profile.h"
#include "search_pipeline/search_stage_decode.h"
#include "search_pipeline/search_stage_kmer_filter.h"
#include "search_pipeline/search_stage_bpm_distance.h"
#include "search_pipeline/search_stage_bpm_align.h"
#include "archive/search/archive_search.h"
#include "archive/search/archive_search_cache.h"
#include "mapper/mapper.h"

/*
 * Archive-Search groups
 */
typedef struct {
  /* Search-Stages buffer to verify candidates */
  search_stage_region_profile_t* stage_region_profile;
  search_stage_decode_t* stage_decode;
  search_stage_kmer_filter_t* stage_kmer_filter;
  search_stage_bpm_distance_t* stage_bpm_distance;
  search_stage_bpm_align_t* stage_bpm_align;
  /* Archive-search cache */
  archive_search_cache_t* archive_search_cache;
  /* Support Data Structures */
  search_pipeline_handlers_t* search_pipeline_handlers;
} search_pipeline_t;

/*
 * Setup
 */
search_pipeline_t* search_pipeline_new(
    mapper_parameters_t* const mapper_parameters,
    gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffers_offset,
    const bool paired_end);
void search_pipeline_clear(search_pipeline_t* const search_pipeline);
void search_pipeline_delete(search_pipeline_t* const search_pipeline);

/*
 * Archive-Search allocation
 */
void search_pipeline_allocate_se(
    search_pipeline_t* const search_pipeline,
    archive_search_t** const archive_search);
void search_pipeline_allocate_pe(
    search_pipeline_t* const search_pipeline,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2);
void search_pipeline_free(
    search_pipeline_t* const search_pipeline,
    archive_search_t* const archive_search_end);

#endif /* SEARCH_PIPELINE_H_ */
