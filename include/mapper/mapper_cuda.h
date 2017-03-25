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
 *   Mapper-CUDA module implements the high-level mapping workflow using GPU(s)
 */

#ifndef MAPPER_CUDA_H_
#define MAPPER_CUDA_H_

#include "mapper/mapper.h"
#include "search_pipeline/search_pipeline.h"
#include "stats/report_stats_mstats.h"

/*
 * Mapper-CUDA Search
 */
typedef struct {
  /* Thread Info */
  uint64_t thread_id;
  pthread_t* thread_data;
  /* Parameters */
  mapper_parameters_t* mapper_parameters;
  /* I/O */
  mapper_io_handler_t* mapper_io_handler;
  /* GPU Buffers */
  gpu_buffer_collection_t* gpu_buffer_collection;
  uint64_t gpu_buffers_offset;
  /* Search Pipeline (and auxiliary variables) */
  search_pipeline_t* search_pipeline;
  archive_search_t* pending_region_profile_end1;
  archive_search_t* pending_region_profile_end2;
  archive_search_t* pending_decode_end1;
  archive_search_t* pending_decode_end2;
  archive_search_t* pending_kmer_filter_end1;
  archive_search_t* pending_kmer_filter_end2;
  archive_search_t* pending_bpm_distance_end1;
  archive_search_t* pending_bpm_distance_end2;
  archive_search_t* pending_bpm_align_end1;
  archive_search_t* pending_bpm_align_end2;
  /* Stats */
  mapping_stats_t* mapping_stats; // Per thread stats report structures
  /* Progress */
  ticker_t* ticker;
  uint64_t reads_processed;
} mapper_cuda_search_t;

/*
 * SE-CUDA Mapper
 */
void mapper_cuda_se_run(mapper_parameters_t* const mapper_parameters);

/*
 * PE-CUDA Mapper
 */
void mapper_cuda_pe_run(mapper_parameters_t* const mapper_parameters);

#endif /* MAPPER_CUDA_H_ */
