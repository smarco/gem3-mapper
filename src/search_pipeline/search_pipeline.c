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

#include "search_pipeline/search_pipeline.h"
#include "archive/search/archive_search_se.h"
#include "archive/search/archive_search_pe.h"

/*
 * Setup
 */
search_pipeline_t* search_pipeline_new(
    mapper_parameters_t* const mapper_parameters,
    gpu_buffer_collection_t* const gpu_buffer_collection,
    const uint64_t buffers_offset,
    const bool paired_end) {
  // Parameters
  const mapper_parameters_cuda_t* const cuda = &mapper_parameters->cuda;
  search_parameters_t* const search_parameters = &mapper_parameters->search_parameters;
  archive_t* const archive = mapper_parameters->archive;
  fm_index_t* const fm_index = archive->fm_index;
  const bool cpu_emulation = cuda->cpu_emulation;
  // Alloc
  search_pipeline_t* const search_pipeline = mm_alloc(search_pipeline_t);
  // Allocate archive-search cache
  search_pipeline->archive_search_cache = archive_search_cache_new(search_parameters);
  // Allocate Support Data Structures
  search_pipeline_handlers_t* const search_pipeline_handlers = search_pipeline_handlers_new(archive);
  search_pipeline->search_pipeline_handlers = search_pipeline_handlers;
  // Allocate pipeline-stage region-profile
  uint64_t acc_buffers_offset = buffers_offset;
  region_profile_model_t* const region_profile_model = &search_parameters->region_profile_model;
  const bool region_profile_enabled = gpu_buffer_collection->gpu_region_profile_available && !cpu_emulation;
  search_parameters->gpu_stage_region_profile_enabled = region_profile_enabled;
  search_pipeline->stage_region_profile = search_stage_region_profile_new(gpu_buffer_collection,
      acc_buffers_offset,cuda->num_fmi_bsearch_buffers,region_profile_enabled,
      (uint32_t)region_profile_model->region_th,(uint32_t)region_profile_model->max_steps,
      (uint32_t)region_profile_model->dec_factor);
  acc_buffers_offset += cuda->num_fmi_bsearch_buffers;
  // Allocate pipeline-stage decode-candidates
  sampled_sa_t* const sampled_sa = fm_index->sampled_sa;
  const uint32_t sampling_rate = sampled_sa_get_sa_sampling_rate(sampled_sa);
  const bool decode_sa_enabled = gpu_buffer_collection->gpu_decode_sa_available && !cpu_emulation;
  const bool decode_text_enabled = gpu_buffer_collection->gpu_decode_text_available && !cpu_emulation;
  search_parameters->gpu_stage_decode_enabled = decode_sa_enabled && decode_text_enabled;
  search_pipeline->stage_decode = search_stage_decode_new(
      gpu_buffer_collection,acc_buffers_offset,cuda->num_fmi_decode_buffers,
      sampling_rate,decode_sa_enabled,decode_text_enabled,search_pipeline_handlers);
  acc_buffers_offset += cuda->num_fmi_decode_buffers;
  // Allocate pipeline-stage kmer-filter
  const bool kmer_filter_enabled = gpu_buffer_collection->gpu_kmer_filter_available && !cpu_emulation;
  search_parameters->gpu_stage_kmer_filter_enabled = kmer_filter_enabled;
  search_pipeline->stage_kmer_filter = search_stage_kmer_filter_new(
      gpu_buffer_collection,acc_buffers_offset,cuda->num_kmer_filter_buffers,
      kmer_filter_enabled,search_pipeline_handlers);
  acc_buffers_offset += cuda->num_kmer_filter_buffers;
  // Allocate pipeline-stage BPM-Distance
  const bool bpm_distance_enabled = gpu_buffer_collection->gpu_bpm_distance_available && !cpu_emulation;
  search_parameters->gpu_stage_bpm_distance_enabled = bpm_distance_enabled;
  search_pipeline->stage_bpm_distance = search_stage_bpm_distance_new(
      gpu_buffer_collection,acc_buffers_offset,cuda->num_bpm_distance_buffers,
      paired_end,bpm_distance_enabled,search_pipeline_handlers);
  acc_buffers_offset += cuda->num_bpm_distance_buffers;
  // Allocate pipeline-stage BPM-Align
  const bool bpm_align_enabled = gpu_buffer_collection->gpu_bpm_align_available && !cpu_emulation;
  search_parameters->gpu_stage_bpm_align_enabled = bpm_align_enabled;
  search_pipeline->stage_bpm_align = search_stage_bpm_align_new(
      gpu_buffer_collection,acc_buffers_offset,cuda->num_bpm_align_buffers,
      paired_end,bpm_align_enabled,search_pipeline_handlers);
  // Return
  return search_pipeline;
}
void search_pipeline_clear(search_pipeline_t* const search_pipeline) {
  search_pipeline_handlers_clear(search_pipeline->search_pipeline_handlers);
}
void search_pipeline_delete(search_pipeline_t* const search_pipeline) {
  search_stage_region_profile_delete(search_pipeline->stage_region_profile);
  search_stage_decode_delete(search_pipeline->stage_decode);
  search_stage_kmer_filter_delete(search_pipeline->stage_kmer_filter);
  search_stage_bpm_distance_delete(search_pipeline->stage_bpm_distance);
  search_stage_bpm_align_delete(search_pipeline->stage_bpm_align);
  archive_search_cache_delete(search_pipeline->archive_search_cache);
  search_pipeline_handlers_delete(search_pipeline->search_pipeline_handlers);
  mm_free(search_pipeline);
}
/*
 * Archive-Search allocation
 */
void search_pipeline_allocate_se(
    search_pipeline_t* const search_pipeline,
    archive_search_t** const archive_search) {
  // Alloc
  archive_search_cache_se_alloc(
      search_pipeline->archive_search_cache,archive_search);
}
void search_pipeline_allocate_pe(
    search_pipeline_t* const search_pipeline,
    archive_search_t** const archive_search_end1,
    archive_search_t** const archive_search_end2) {
  // Alloc
  archive_search_cache_pe_alloc(
      search_pipeline->archive_search_cache,
      archive_search_end1,archive_search_end2);
}
void search_pipeline_free(
    search_pipeline_t* const search_pipeline,
    archive_search_t* const archive_search_end) {
  // Destroy
  archive_search_destroy(archive_search_end);
  // Free
  archive_search_cache_free(
      search_pipeline->archive_search_cache,archive_search_end);
}

