/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2013-2017 by Alejandro Chacon <alejandro.chacond@gmail.com>
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
 *            Alejandro Chacon <alejandro.chacond@gmail.com>
 * DESCRIPTION:
 *   Mapper-CUDA module implements the high-level mapping workflow
 *   for PE-mapping using GPU(s)
 */

#include "mapper/mapper_cuda_pe.h"
#include "mapper/mapper.h"
#include "mapper/mapper_io.h"
#include "archive/search/archive_search_pe.h"
#include "archive/search/archive_search_pe_stepwise.h"
#include "text/sequence_bisulfite.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * Check Occupancy
 */
bool mapper_pe_cuda_stage_read_input_sequences_exhausted(mapper_cuda_search_t* const mapper_search) {
  // Check pending search
  if (mapper_search->pending_region_profile_end1!=NULL) return false;
  // Check end_of_block
  if (!buffered_input_file_eob(mapper_search->mapper_io_handler->buffered_fasta_input_end1)) return false;
  // Reload buffer
  const uint64_t error_code = mapper_pe_reload_buffer(mapper_search->mapper_io_handler);
  if (error_code==INPUT_STATUS_EOF) return true;
  // Clear pipeline (release intermediate memory & start pipeline fresh)
  search_pipeline_clear(mapper_search->search_pipeline);
  return false;
}
bool mapper_pe_cuda_stage_region_profile_output_exhausted(mapper_cuda_search_t* const mapper_search) {
  // Check input-stage (Region-Profile)
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  if (!search_stage_region_profile_retrieve_finished(search_pipeline->stage_region_profile)) return false;
  // Check pending searches at current stage
  if (mapper_search->pending_decode_end1!=NULL) return false;
  // Exhausted
  return true;
}
bool mapper_pe_cuda_stage_decode_output_exhausted(mapper_cuda_search_t* const mapper_search) {
  // Check input-stage (Decode-Candidates)
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  if (!search_stage_decode_retrieve_finished(search_pipeline->stage_decode)) return false;
  /// Check pending searches at current stage
  if (mapper_search->pending_kmer_filter_end1!=NULL) return false;
  // Exhausted
  return true;
}
bool mapper_pe_cuda_stage_kmer_filter_output_exhausted(mapper_cuda_search_t* const mapper_search) {
  // Check input-stage (Kmer-filter)
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  if (!search_stage_kmer_filter_retrieve_finished(search_pipeline->stage_kmer_filter)) return false;
  /// Check pending searches at current stage
  if (mapper_search->pending_bpm_distance_end1!=NULL) return false;
  // Exhausted
  return true;
}
bool mapper_pe_cuda_stage_bpm_distance_output_exhausted(mapper_cuda_search_t* const mapper_search) {
  // Check input-stage (BPM-Distance)
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  if (!search_stage_bpm_distance_retrieve_finished(search_pipeline->stage_bpm_distance)) return false;
  /// Check pending searches at current stage
  if (mapper_search->pending_bpm_align_end1!=NULL) return false;
  // Exhausted
  return true;
}
/*
 * Region Profile
 */
void mapper_pe_cuda_region_profile(mapper_cuda_search_t* const mapper_search) {
  PROFILE_START(GP_MAPPER_CUDA_PE_REGION_PROFILE,PROFILE_LEVEL);
  // Parameters
  mapper_io_handler_t* const mapper_io_handler = mapper_search->mapper_io_handler;
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  search_stage_region_profile_t* const stage_region_profile = search_pipeline->stage_region_profile;
  archive_search_t *archive_search_end1 = NULL, *archive_search_end2 = NULL;
  sequence_t *sequence_end1, *sequence_end2;
  // Reschedule search (that couldn't fit into the buffer)
  if (mapper_search->pending_region_profile_end1!=NULL) {
    search_stage_region_profile_send_pe_search(stage_region_profile,
        mapper_search->pending_region_profile_end1,
        mapper_search->pending_region_profile_end2);
    mapper_search->pending_region_profile_end1 = NULL;
  }
  // Generation. Keep processing the current input-block
  while (mapper_read_paired_sequence(mapper_io_handler,false,&sequence_end1,&sequence_end2)) {
    // Prepare archive-search
    search_pipeline_allocate_pe(search_pipeline,&archive_search_end1,&archive_search_end2);
    // Inject Support Data Structures
    search_pipeline_handlers_prepare_pe(
        archive_search_end1,archive_search_end2,
        sequence_end1,sequence_end2,
		  no_conversion,no_conversion,
        search_pipeline->search_pipeline_handlers);
    // Generate Candidates (Search into the archive)
    archive_search_pe_stepwise_init_search(archive_search_end1,archive_search_end2);
#ifdef GPU_REGION_PROFILE_ADAPTIVE
    archive_search_pe_stepwise_region_profile_adaptive_generate(archive_search_end1,archive_search_end2);
#else
    archive_search_pe_stepwise_region_profile_static_generate(archive_search_end1,archive_search_end2);
#endif
    // Send to stage region-profile
    const bool search_sent = search_stage_region_profile_send_pe_search(
        stage_region_profile,archive_search_end1,archive_search_end2);
    if (!search_sent) {
      mapper_search->pending_region_profile_end1 = archive_search_end1; // Pending Search (end/1)
      mapper_search->pending_region_profile_end2 = archive_search_end2; // Pending Search (end/2)
      break;
    }
  }
  PROFILE_STOP(GP_MAPPER_CUDA_PE_REGION_PROFILE,PROFILE_LEVEL);
}
/*
 * Decode Candidates
 */
void mapper_pe_cuda_decode(mapper_cuda_search_t* const mapper_search) {
  PROFILE_START(GP_MAPPER_CUDA_PE_DECODE,PROFILE_LEVEL);
  // Parameters
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  search_stage_region_profile_t* const stage_region_profile = search_pipeline->stage_region_profile;
  search_stage_decode_t* const stage_decode_candidates = search_pipeline->stage_decode;
  archive_search_t *archive_search_end1 = NULL, *archive_search_end2 = NULL;
  // Reschedule search (that couldn't fit into the buffer)
  if (mapper_search->pending_decode_end1!=NULL) {
    search_stage_decode_send_pe_search(stage_decode_candidates,
        mapper_search->pending_decode_end1,mapper_search->pending_decode_end2);
    mapper_search->pending_decode_end1 = NULL;
  }
  // Read from stage region-profile
  bool pending_searches;
  while ((pending_searches=search_stage_region_profile_retrieve_pe_search(
      stage_region_profile,&archive_search_end1,&archive_search_end2))) {
    // Send to stage decode-candidates
    const bool search_sent = search_stage_decode_send_pe_search(
        stage_decode_candidates,archive_search_end1,archive_search_end2);
    if (!search_sent) {
      mapper_search->pending_decode_end1 = archive_search_end1; // Pending Search (end/1)
      mapper_search->pending_decode_end2 = archive_search_end2; // Pending Search (end/2)
      break;
    }
  }
  // Clean
  if (!pending_searches) search_stage_region_profile_clear(stage_region_profile);
  PROFILE_STOP(GP_MAPPER_CUDA_PE_DECODE,PROFILE_LEVEL);
}
/*
 * Kmer-filter Candidates
 */
void mapper_pe_cuda_kmer_filter(mapper_cuda_search_t* const mapper_search) {
  PROFILE_START(GP_MAPPER_CUDA_PE_KMER_FILTER,PROFILE_LEVEL);
  // Parameters
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  search_stage_decode_t* const stage_decode = search_pipeline->stage_decode;
  search_stage_kmer_filter_t* const stage_kmer_filter = search_pipeline->stage_kmer_filter;
  archive_search_t *archive_search_end1 = NULL, *archive_search_end2 = NULL;
  // Reschedule search (that couldn't fit into the buffer)
  if (mapper_search->pending_kmer_filter_end1!=NULL) {
    search_stage_kmer_filter_send_pe_search(stage_kmer_filter,
        mapper_search->pending_kmer_filter_end1,mapper_search->pending_kmer_filter_end2);
    mapper_search->pending_kmer_filter_end1 = NULL;
  }
  // Read from stage decode-candidates
  bool pending_searches;
  while ((pending_searches=search_stage_decode_retrieve_pe_search(
      stage_decode,&archive_search_end1,&archive_search_end2))) {
    // Send to stage kmer-filter
    const bool search_sent = search_stage_kmer_filter_send_pe_search(
        stage_kmer_filter,archive_search_end1,archive_search_end2);
    if (!search_sent) {
      mapper_search->pending_kmer_filter_end1 = archive_search_end1; // Pending Search (end/1)
      mapper_search->pending_kmer_filter_end2 = archive_search_end2; // Pending Search (end/2)
      break;
    }
  }
  // Clean
  if (!pending_searches) search_stage_decode_clear(stage_decode);
  PROFILE_STOP(GP_MAPPER_CUDA_PE_KMER_FILTER,PROFILE_LEVEL);
}
/*
 * BPM-Distance Candidates
 */
void mapper_pe_cuda_bpm_distance(mapper_cuda_search_t* const mapper_search) {
  PROFILE_START(GP_MAPPER_CUDA_PE_BPM_DISTANCE,PROFILE_LEVEL);
  // Parameters
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  search_stage_kmer_filter_t* const stage_kmer_filter = search_pipeline->stage_kmer_filter;
  search_stage_bpm_distance_t* const stage_bpm_distance = search_pipeline->stage_bpm_distance;
  archive_search_t *archive_search_end1 = NULL, *archive_search_end2 = NULL;
  // Reschedule search (that couldn't fit into the buffer)
  if (mapper_search->pending_bpm_distance_end1!=NULL) {
    search_stage_bpm_distance_send_pe_search(stage_bpm_distance,
        mapper_search->pending_bpm_distance_end1,
        mapper_search->pending_bpm_distance_end2);
    mapper_search->pending_bpm_distance_end1 = NULL;
  }
  // Read from stage kmer-filter
  bool pending_searches;
  while ((pending_searches=search_stage_kmer_filter_retrieve_pe_search(
      stage_kmer_filter,&archive_search_end1,&archive_search_end2))) {
    // Send to stage BPM-Distance
    const bool search_sent = search_stage_bpm_distance_send_pe_search(
        stage_bpm_distance,archive_search_end1,archive_search_end2);
    if (!search_sent) {
      mapper_search->pending_bpm_distance_end1 = archive_search_end1; // Pending Search (end/1)
      mapper_search->pending_bpm_distance_end2 = archive_search_end2; // Pending Search (end/2)
      break;
    }
  }
  // Clean
  if (!pending_searches) search_stage_kmer_filter_clear(stage_kmer_filter);
  PROFILE_STOP(GP_MAPPER_CUDA_PE_BPM_DISTANCE,PROFILE_LEVEL);
}
/*
 * BPM-Align Candidates
 */
void mapper_pe_cuda_bpm_align(mapper_cuda_search_t* const mapper_search) {
  PROFILE_START(GP_MAPPER_CUDA_PE_BPM_ALIGN,PROFILE_LEVEL);
  // Parameters
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  search_stage_bpm_distance_t* const stage_bpm_distance = search_pipeline->stage_bpm_distance;
  search_stage_bpm_align_t* const stage_bpm_align = search_pipeline->stage_bpm_align;
  archive_search_t *archive_search_end1 = NULL, *archive_search_end2 = NULL;
  // Reschedule search (that couldn't fit into the buffer)
  if (mapper_search->pending_bpm_align_end1!=NULL) {
    search_stage_bpm_align_send_pe_search(stage_bpm_align,
        mapper_search->pending_bpm_align_end1,
        mapper_search->pending_bpm_align_end2);
    mapper_search->pending_bpm_align_end1 = NULL;
  }
  // Read from stage BPM-Distance
  bool pending_searches;
  while ((pending_searches=search_stage_bpm_distance_retrieve_pe_search(
      stage_bpm_distance,&archive_search_end1,&archive_search_end2))) {
    // Send to stage BPM-Align
    const bool search_sent = search_stage_bpm_align_send_pe_search(
        stage_bpm_align,archive_search_end1,archive_search_end2);
    if (!search_sent) {
      mapper_search->pending_bpm_align_end1 = archive_search_end1; // Pending Search (end/1)
      mapper_search->pending_bpm_align_end2 = archive_search_end2; // Pending Search (end/2)
      break;
    }
  }
  // Clean
  if (!pending_searches) search_stage_bpm_distance_clear(stage_bpm_distance);
  PROFILE_STOP(GP_MAPPER_CUDA_PE_BPM_ALIGN,PROFILE_LEVEL);
}
/*
 * Finish Search
 */
void mapper_pe_cuda_finish_search(mapper_cuda_search_t* const mapper_search) {
  PROFILE_START(GP_MAPPER_CUDA_PE_FINISH_SEARCH,PROFILE_LEVEL);
  // Parameters
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  search_pipeline_t* const search_pipeline = mapper_search->search_pipeline;
  search_stage_bpm_align_t* const stage_bpm_align = search_pipeline->stage_bpm_align;
  archive_search_t *archive_search_end1 = NULL, *archive_search_end2 = NULL;
  // Read from stage BPM-Align
  while (search_stage_bpm_align_retrieve_pe_search(
      stage_bpm_align,&archive_search_end1,&archive_search_end2)) {
    // Finish Search
    archive_search_pe_stepwise_finish_search(archive_search_end1,
        archive_search_end2,stage_bpm_align->paired_matches); // Finish search
    // Calculate Stats
    const uint64_t total_candidates_analized_gpu_end1 = search_pipeline->search_pipeline_handlers->fc_bpm_align_end1.total_candidates_canonical_gpu;
    const uint64_t total_candidates_analized_gpu_end2 = search_pipeline->search_pipeline_handlers->fc_bpm_align_end2.total_candidates_canonical_gpu;
    const uint64_t id_histogram_end1 = MIN(total_candidates_analized_gpu_end1, GEM_HIST_CAND_ALIGNED-1);
    const uint64_t id_histogram_end2 = MIN(total_candidates_analized_gpu_end2, GEM_HIST_CAND_ALIGNED-1);
    COUNTER_ADD(&search_pipeline->search_pipeline_handlers->fc_bpm_distance_end1.candidates_aligned_histo[id_histogram_end1],
                 search_pipeline->search_pipeline_handlers->fc_bpm_align_end1.total_candidates_realigned_swg);
    COUNTER_ADD(&search_pipeline->search_pipeline_handlers->fc_bpm_distance_end2.candidates_aligned_histo[id_histogram_end2],
                 search_pipeline->search_pipeline_handlers->fc_bpm_align_end2.total_candidates_realigned_swg);
    // Output Matches
    mapper_io_handler_output_paired_matches(
        mapper_search->mapper_io_handler,archive_search_end1,archive_search_end2,
        stage_bpm_align->paired_matches,mapper_search->mapping_stats);
    // Update processed
    if (++mapper_search->reads_processed == parameters->io.mapper_ticker_step) {
      ticker_update_mutex(mapper_search->ticker,mapper_search->reads_processed);
      mapper_search->reads_processed=0;
    }
    // Free
    search_pipeline_free(search_pipeline,archive_search_end1);
    search_pipeline_free(search_pipeline,archive_search_end2);
  }
  // Clean
  search_stage_bpm_align_clear(stage_bpm_align);
  PROFILE_STOP(GP_MAPPER_CUDA_PE_FINISH_SEARCH,PROFILE_LEVEL);
}
/*
 * Mapper PE-CUDA
 */
void mapper_cuda_pe_thread_pipeline(mapper_cuda_search_t* const mapper_search) {
  while (!mapper_pe_cuda_stage_read_input_sequences_exhausted(mapper_search)) {
    // Region Profile
    mapper_pe_cuda_region_profile(mapper_search);
    do {
      // Decode
      mapper_pe_cuda_decode(mapper_search);
      do {
        // Kmer-filter
        mapper_pe_cuda_kmer_filter(mapper_search);
        do {
          // BPM-Distance
          mapper_pe_cuda_bpm_distance(mapper_search);
          do {
            // BPM-Align
            mapper_pe_cuda_bpm_align(mapper_search);
            // Finish Search
            mapper_pe_cuda_finish_search(mapper_search);
          } while (!mapper_pe_cuda_stage_bpm_distance_output_exhausted(mapper_search));
        } while (!mapper_pe_cuda_stage_kmer_filter_output_exhausted(mapper_search));
      } while (!mapper_pe_cuda_stage_decode_output_exhausted(mapper_search));
    } while (!mapper_pe_cuda_stage_region_profile_output_exhausted(mapper_search));
  }
}
void* mapper_cuda_pe_thread(mapper_cuda_search_t* const mapper_search) {
  // GEM-thread error handler
  gem_thread_register_id(mapper_search->thread_id+1);
  PROFILE_START(GP_MAPPER_CUDA_PE,PROFILE_LEVEL);
  // Parameters
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  const mapper_parameters_cuda_t* const cuda_parameters = &parameters->cuda;
  // Create search-pipeline & initialize matches
  mapper_search->search_pipeline = search_pipeline_new(parameters,
      mapper_search->gpu_buffer_collection,mapper_search->gpu_buffers_offset,true);
  mapper_search->pending_region_profile_end1 = NULL;
  mapper_search->pending_decode_end1 = NULL;
  mapper_search->pending_kmer_filter_end1 = NULL;
  mapper_search->pending_bpm_distance_end1 = NULL;
  mapper_search->pending_bpm_align_end1 = NULL;
  // Create new I/O handler
  mapper_search->mapper_io_handler = mapper_io_handler_new_pe(
      parameters,cuda_parameters->input_buffer_size,
      mapper_search->search_pipeline->search_pipeline_handlers->mapper_stats,
      mapper_search->search_pipeline->search_pipeline_handlers->mm_allocator);
  // FASTA/FASTQ reading loop
  mapper_search->reads_processed = 0;
  mapper_cuda_pe_thread_pipeline(mapper_search);
  // Clean up
  ticker_update_mutex(mapper_search->ticker,mapper_search->reads_processed);
  search_pipeline_delete(mapper_search->search_pipeline);
  mapper_io_handler_delete(mapper_search->mapper_io_handler);
  PROFILE_STOP(GP_MAPPER_CUDA_PE,PROFILE_LEVEL);
  pthread_exit(0);
}
