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
 *   Archive-Search Single-End module encapsulating
 *   the basic SE-search stages in a step-wise approach. Mainly
 *   used to horizontally process a batch of searches (doing
 *   each stage for all searches at the same time before
 *   progressing to the next search stage)
 */

#include "approximate_search/approximate_search_stages.h"
#include "archive/search/archive_search_se_stepwise.h"
#include "archive/search/archive_search_se.h"
#include "archive/search/archive_select.h"
#include "archive/search/archive_check.h"
#include "archive/score/archive_score_se.h"
#include "approximate_search/approximate_search_stepwise.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * Debug
 */
#define DEBUG_ARCHIVE_SEARCH_SE_STEPWISE GEM_DEEP_DEBUG

void archive_search_se_stepwise_debug_prologue(
    archive_search_t* const archive_search,
    const char* const label) {
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: %s (stage=%s,%s) (Tag=%s)\n",
        label,asearch_stage_label[archive_search->approximate_search.search_stage],
        asearch_processing_state_label[archive_search->approximate_search.processing_state],
        archive_search->sequence->tag.buffer);
    tab_global_inc();
  }
}
void archive_search_se_stepwise_debug_epilogue() {
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) { tab_global_dec(); }
}
/*
 * Stepwise: Region-Profile
 */
void archive_search_se_stepwise_region_profile_static_generate(archive_search_t* const archive_search) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_prologue(archive_search,"Region.Profile.Static.Generate.Static");
  // Region-Profile Generate
  approximate_search_stepwise_region_profile_static_generate(&archive_search->approximate_search);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
}
void archive_search_se_stepwise_region_profile_static_copy(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_COPY,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_prologue(archive_search,"Region.Profile.Static.Copy");
  // Region-Profile Copy
  approximate_search_stepwise_region_profile_static_copy(
      &archive_search->approximate_search,gpu_buffer_fmi_ssearch);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_COPY,PROFILE_LEVEL);
}
void archive_search_se_stepwise_region_profile_static_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_prologue(archive_search,"Region.Profile.Static.Retrieve");
  // Region-Profile Retrieve
  approximate_search_stepwise_region_profile_static_retrieve(
      &archive_search->approximate_search,gpu_buffer_fmi_ssearch);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
}
void archive_search_se_stepwise_region_profile_adaptive_generate(archive_search_t* const archive_search) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_prologue(archive_search,"Region.Profile.Adaptive.Generate.Adaptive");
  // Region-Profile Generate
  approximate_search_stepwise_region_profile_adaptive_generate(&archive_search->approximate_search);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
}
void archive_search_se_stepwise_region_profile_adaptive_copy(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_COPY,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_prologue(archive_search,"Region.Profile.Adaptive.Copy");
  // Region-Profile Copy
  approximate_search_stepwise_region_profile_adaptive_copy(
      &archive_search->approximate_search,gpu_buffer_fmi_asearch);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_COPY,PROFILE_LEVEL);
}
void archive_search_se_stepwise_region_profile_adaptive_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_prologue(archive_search,"Region.Profile.Adaptive.Retrieve");
  // Region-Profile Retrieve
  approximate_search_stepwise_region_profile_adaptive_retrieve(
      &archive_search->approximate_search,gpu_buffer_fmi_asearch);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: Decode-Candidates
 */
void archive_search_se_stepwise_decode_copy(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_DECODE_COPY,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_prologue(archive_search,"Decode.Candidates.Copy");
  // Decode-Candidates Copy
  approximate_search_stepwise_decode_copy(
      &archive_search->approximate_search,gpu_buffer_fmi_decode);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_DECODE_COPY,PROFILE_LEVEL);
}
void archive_search_se_stepwise_decode_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_DECODE_RETRIEVE,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_prologue(archive_search,"Decode.Candidates.Retrieve");
  // Decode-Candidates Retrieve
  approximate_search_stepwise_decode_retrieve(
      &archive_search->approximate_search,gpu_buffer_fmi_decode);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_DECODE_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: Kmer-filter
 */
void archive_search_se_stepwise_kmer_filter_copy(
    archive_search_t* const archive_search,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_KMER_FILTER_COPY,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_prologue(archive_search,"Kmer.Filter.Copy");
  // Kmer-filter Copy
  approximate_search_stepwise_kmer_filter_copy(
      &archive_search->approximate_search,gpu_buffer_kmer_filter);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_KMER_FILTER_COPY,PROFILE_LEVEL);
}
void archive_search_se_stepwise_kmer_filter_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_KMER_FILTER_RETRIEVE,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_prologue(archive_search,"Kmer.Filter.Retrieve");
  // Kmer-filter Retrieve
  approximate_search_stepwise_kmer_filter_retrieve(
      &archive_search->approximate_search,gpu_buffer_kmer_filter);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_KMER_FILTER_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: BPM-Distance
 */
void archive_search_se_stepwise_bpm_distance_copy(
    archive_search_t* const archive_search,
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_BPM_DISTANCE_COPY,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_prologue(archive_search,"BPM.Distance.Copy");
  // Verify-Candidates Copy
  approximate_search_stepwise_bpm_distance_copy(
      &archive_search->approximate_search,gpu_buffer_bpm_distance);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_BPM_DISTANCE_COPY,PROFILE_LEVEL);
}
void archive_search_se_stepwise_bpm_distance_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_BPM_DISTANCE_RETRIEVE,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_prologue(archive_search,"BPM.Distance.Retrieve");
  // Verify-Candidates Retrieve
  approximate_search_stepwise_bpm_distance_retrieve(
      &archive_search->approximate_search,gpu_buffer_bpm_distance);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_BPM_DISTANCE_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: BPM-Align
 */
void archive_search_se_stepwise_bpm_align_update(
    archive_search_t* const archive_search) {
  // Verify-Candidates Copy
  approximate_search_stepwise_bpm_align_update(
      &archive_search->approximate_search);
}
void archive_search_se_stepwise_bpm_align_copy(
    archive_search_t* const archive_search,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_BPM_ALIGN_COPY,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_prologue(archive_search,"BPM.Align.Copy");
  // Verify-Candidates Copy
  approximate_search_stepwise_bpm_align_copy(
      &archive_search->approximate_search,gpu_buffer_bpm_align);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_BPM_ALIGN_COPY,PROFILE_LEVEL);
}
void archive_search_se_stepwise_bpm_align_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_BPM_ALIGN_RETRIEVE,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_prologue(archive_search,"BPM.Align.Retrieve");
  // Verify-Candidates Retrieve
  approximate_search_stepwise_bpm_align_retrieve(
      &archive_search->approximate_search,gpu_buffer_bpm_align,matches);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_BPM_ALIGN_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: Finish Search
 */
void archive_search_se_stepwise_finish_search(
    archive_search_t* const archive_search,
    matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH,PROFILE_LEVEL);
  // DEBUG
  archive_search_se_stepwise_debug_prologue(archive_search,"Finish");
  // Finish Search
  approximate_search(&archive_search->approximate_search,matches);
  // Select Matches
  search_parameters_t* const search_parameters = &archive_search->search_parameters;
  select_parameters_t* const select_parameters = &search_parameters->select_parameters;
  if (!search_parameters->search_paired_parameters.paired_end_search) {
    archive_select_se_matches(select_parameters,matches);
  }
  // Score Matches (Select alignment-Model and process accordingly)
  archive_score_matches_se(archive_search,matches);
  // Check matches
  if (search_parameters->check_type!=archive_check_nothing) {
    archive_check_se_matches(
        archive_search->archive,search_parameters->match_alignment_model,
        &search_parameters->swg_penalties,archive_search->sequence,
        matches,search_parameters->check_type,archive_search->mm_allocator);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_global_inc();
    archive_search_se_print(gem_log_get_stream(),archive_search,matches);
    tab_global_dec();
    tab_global_dec();
  }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH,PROFILE_LEVEL);
}
