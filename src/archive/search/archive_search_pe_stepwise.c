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
 *   Archive-Search Paired-End module encapsulating
 *   the basic PE-search stages in a step-wise approach. Mainly
 *   used to horizontally process a batch of searches (doing
 *   each stage for all searches at the same time before
 *   progressing to the next search stage)
 */

#include "archive/search/archive_search_se_parameters.h"
#include "archive/search/archive_search_pe_stepwise.h"
#include "archive/search/archive_search_se_stepwise.h"
#include "archive/search/archive_search.h"
#include "archive/search/archive_search_pe.h"
#include "archive/search/archive_select.h"
#include "archive/search/archive_check.h"
#include "archive/score/archive_score_se.h"
#include "archive/score/archive_score_pe.h"
#include "approximate_search/approximate_search_stepwise.h"

/*
 * Debug
 */
#define DEBUG_ARCHIVE_SEARCH_PE_STEPWISE GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * Debug
 */
void archive_search_pe_stepwise_debug_preface(
    archive_search_t* const archive_search,
    const char* const label) {
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) {
    tab_fprintf(stderr,
        "[GEM]>ArchiveSearch.STEPWISE.PE :: %s (Tag=%s)\n",
        label,archive_search->sequence->tag.buffer);
    tab_global_inc();
  }
}
void archive_search_pe_stepwise_debug_epilogue() {
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) { tab_global_dec(); }
}

/*
 * Stepwise: Init Search
 */
void archive_search_pe_stepwise_init_search(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_INIT,PROFILE_LEVEL);
  archive_search_end1->searched = false;
  archive_search_end2->searched = false;
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.PE :: Region.Profile.Generate\n");
    tab_fprintf(gem_log_get_stream(),"  => Tag %s\n",archive_search_end1->sequence->tag.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence.End1 %s\n",archive_search_end1->sequence->read.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence.End2 %s\n",archive_search_end2->sequence->read.buffer);
  }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_INIT,PROFILE_LEVEL);
}
/*
 * Stepwise: Region-Profile
 */
void archive_search_pe_stepwise_region_profile_static_generate(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
  // DEBUG
  archive_search_pe_stepwise_debug_preface(archive_search_end1,"Region.Profile.Generate");
  // Region-Profile Generate
  approximate_search_stepwise_region_profile_static_generate(&archive_search_end1->approximate_search);
  approximate_search_stepwise_region_profile_static_generate(&archive_search_end2->approximate_search);
  // DEBUG
  archive_search_pe_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
}
void archive_search_pe_stepwise_region_profile_static_copy(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_COPY,PROFILE_LEVEL);
  // DEBUG
  archive_search_pe_stepwise_debug_preface(archive_search_end1,"Region.Profile.Copy");
  // Region-Profile Copy
  approximate_search_stepwise_region_profile_static_copy(
      &archive_search_end1->approximate_search,gpu_buffer_fmi_ssearch);
  approximate_search_stepwise_region_profile_static_copy(
      &archive_search_end2->approximate_search,gpu_buffer_fmi_ssearch);
  // DEBUG
  archive_search_pe_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_COPY,PROFILE_LEVEL);
}
void archive_search_pe_stepwise_region_profile_static_retrieve(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
  // DEBUG
  archive_search_pe_stepwise_debug_preface(archive_search_end1,"Region.Profile.Retrieve");
  // Region-Profile Retrieve
  approximate_search_stepwise_region_profile_static_retrieve(
      &archive_search_end1->approximate_search,gpu_buffer_fmi_ssearch);
  approximate_search_stepwise_region_profile_static_retrieve(
      &archive_search_end2->approximate_search,gpu_buffer_fmi_ssearch);
  // DEBUG
  archive_search_pe_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
}
void archive_search_pe_stepwise_region_profile_adaptive_generate(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
  // DEBUG
  archive_search_pe_stepwise_debug_preface(archive_search_end1,"Region.Profile.Generate");
  // Region-Profile Generate
  approximate_search_stepwise_region_profile_adaptive_generate(&archive_search_end1->approximate_search);
  approximate_search_stepwise_region_profile_adaptive_generate(&archive_search_end2->approximate_search);
  // DEBUG
  archive_search_pe_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
}
void archive_search_pe_stepwise_region_profile_adaptive_copy(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_COPY,PROFILE_LEVEL);
  // DEBUG
  archive_search_pe_stepwise_debug_preface(archive_search_end1,"Region.Profile.Copy");
  // Region-Profile Copy
  approximate_search_stepwise_region_profile_adaptive_copy(
      &archive_search_end1->approximate_search,gpu_buffer_fmi_asearch);
  approximate_search_stepwise_region_profile_adaptive_copy(
      &archive_search_end2->approximate_search,gpu_buffer_fmi_asearch);
  // DEBUG
  archive_search_pe_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_COPY,PROFILE_LEVEL);
}
void archive_search_pe_stepwise_region_profile_adaptive_retrieve(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
  // DEBUG
  archive_search_pe_stepwise_debug_preface(archive_search_end1,"Region.Profile.Retrieve");
  // Region-Profile Retrieve
  approximate_search_stepwise_region_profile_adaptive_retrieve(
      &archive_search_end1->approximate_search,gpu_buffer_fmi_asearch);
  approximate_search_stepwise_region_profile_adaptive_retrieve(
      &archive_search_end2->approximate_search,gpu_buffer_fmi_asearch);
  // DEBUG
  archive_search_pe_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: Decode-Candidates
 */
void archive_search_pe_stepwise_decode_copy(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_DECODE_COPY,PROFILE_LEVEL);
  // DEBUG
  archive_search_pe_stepwise_debug_preface(archive_search_end1,"Decode.Candidates.Copy");
  // Decode-Candidates Copy
  approximate_search_stepwise_decode_copy(
      &archive_search_end1->approximate_search,gpu_buffer_fmi_decode);
  approximate_search_stepwise_decode_copy(
      &archive_search_end2->approximate_search,gpu_buffer_fmi_decode);
  // DEBUG
  archive_search_pe_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_DECODE_COPY,PROFILE_LEVEL);
}
void archive_search_pe_stepwise_decode_retrieve(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_DECODE_RETRIEVE,PROFILE_LEVEL);
  // DEBUG
  archive_search_pe_stepwise_debug_preface(archive_search_end1,"Decode.Candidates.Retrieve");
  // Decode-Candidates Retrieve
  approximate_search_stepwise_decode_retrieve(
      &archive_search_end1->approximate_search,gpu_buffer_fmi_decode);
  approximate_search_stepwise_decode_retrieve(
      &archive_search_end2->approximate_search,gpu_buffer_fmi_decode);
  // DEBUG
  archive_search_pe_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_DECODE_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: Kmer-filter
 */
void archive_search_pe_stepwise_kmer_filter_copy(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_KMER_FILTER_COPY,PROFILE_LEVEL);
  // DEBUG
  archive_search_pe_stepwise_debug_preface(archive_search_end1,"Kmer.Filter.Copy");
  // Kmer-filter Copy
  approximate_search_stepwise_kmer_filter_copy(
      &archive_search_end1->approximate_search,gpu_buffer_kmer_filter);
  approximate_search_stepwise_kmer_filter_copy(
      &archive_search_end2->approximate_search,gpu_buffer_kmer_filter);
  // DEBUG
  archive_search_pe_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_KMER_FILTER_COPY,PROFILE_LEVEL);
}
void archive_search_pe_stepwise_kmer_filter_retrieve(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_KMER_FILTER_RETRIEVE,PROFILE_LEVEL);
  // DEBUG
  archive_search_pe_stepwise_debug_preface(archive_search_end1,"Kmer.Filter.Retrieve");
  // Kmer-filter Retrieve
  approximate_search_stepwise_kmer_filter_retrieve(
      &archive_search_end1->approximate_search,gpu_buffer_kmer_filter);
  approximate_search_stepwise_kmer_filter_retrieve(
      &archive_search_end2->approximate_search,gpu_buffer_kmer_filter);
  // DEBUG
  archive_search_pe_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_KMER_FILTER_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: BPM-Distance
 */
void archive_search_pe_stepwise_bpm_distance_copy(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_BPM_DISTANCE_COPY,PROFILE_LEVEL);
  // DEBUG
  archive_search_pe_stepwise_debug_preface(archive_search_end1,"BPM.Distance.Copy");
  // Verify-Candidates Copy
  approximate_search_stepwise_bpm_distance_copy(
      &archive_search_end1->approximate_search,gpu_buffer_bpm_distance);
  approximate_search_stepwise_bpm_distance_copy(
      &archive_search_end2->approximate_search,gpu_buffer_bpm_distance);
  // DEBUG
  archive_search_pe_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_BPM_DISTANCE_COPY,PROFILE_LEVEL);
}
void archive_search_pe_stepwise_bpm_distance_retrieve(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_BPM_DISTANCE_RETRIEVE,PROFILE_LEVEL);
  // DEBUG
  archive_search_pe_stepwise_debug_preface(archive_search_end1,"BPM.Distance.Retrieve");
  // Verify-Candidates Retrieve
  approximate_search_stepwise_bpm_distance_retrieve(
      &archive_search_end1->approximate_search,gpu_buffer_bpm_distance);
  approximate_search_stepwise_bpm_distance_retrieve(
      &archive_search_end2->approximate_search,gpu_buffer_bpm_distance);
  // DEBUG
  archive_search_pe_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_BPM_DISTANCE_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: BPM-Align
 */
void archive_search_pe_stepwise_bpm_align_update(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  // Verify-Candidates Copy
  approximate_search_stepwise_bpm_align_update(
      &archive_search_end1->approximate_search);
  approximate_search_stepwise_bpm_align_update(
      &archive_search_end2->approximate_search);
}
void archive_search_pe_stepwise_bpm_align_copy(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_BPM_ALIGN_COPY,PROFILE_LEVEL);
  // DEBUG
  archive_search_pe_stepwise_debug_preface(archive_search_end1,"BPM.Align.Copy");
  // Verify-Candidates Copy
  approximate_search_stepwise_bpm_align_copy(
      &archive_search_end1->approximate_search,gpu_buffer_bpm_align);
  approximate_search_stepwise_bpm_align_copy(
      &archive_search_end2->approximate_search,gpu_buffer_bpm_align);
  // DEBUG
  archive_search_pe_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_BPM_ALIGN_COPY,PROFILE_LEVEL);
}
void archive_search_pe_stepwise_bpm_align_retrieve(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_BPM_ALIGN_RETRIEVE,PROFILE_LEVEL);
  // DEBUG
  archive_search_pe_stepwise_debug_preface(archive_search_end1,"BPM.Align.Retrieve");
  // Verify-Candidates Retrieve
  approximate_search_stepwise_bpm_align_retrieve(
      &archive_search_end1->approximate_search,gpu_buffer_bpm_align,matches);
  approximate_search_stepwise_bpm_align_retrieve(
      &archive_search_end2->approximate_search,gpu_buffer_bpm_align,matches);
  // DEBUG
  archive_search_pe_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_BPM_ALIGN_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: Finish Search
 */
void archive_search_pe_stepwise_finish_search(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH,PROFILE_LEVEL);
  // DEBUG
  archive_search_pe_stepwise_debug_preface(archive_search_end1,"Finish");
  // Finish SE-Search
  search_parameters_t* const search_parameters = &archive_search_end1->search_parameters;
  archive_search_se_stepwise_finish_search(archive_search_end1,paired_matches->matches_end1);
  archive_search_se_stepwise_finish_search(archive_search_end2,paired_matches->matches_end2);
  archive_search_end1->pe_search_state = archive_search_pe_state_find_pairs;
  archive_search_end1->searched = true;
  archive_search_end2->searched = true;
  archive_search_pe_continue(archive_search_end1,archive_search_end2,paired_matches);
  // PE Score (Select alignment-Model and process accordingly)
  archive_score_matches_pe(archive_search_end1,archive_search_end2,paired_matches);
  // PE Check matches
  if (search_parameters->check_type!=archive_check_nothing) {
    archive_check_pe_matches(
        archive_search_end1->archive,search_parameters->match_alignment_model,
        &search_parameters->swg_penalties,archive_search_end1->sequence,
        archive_search_end2->sequence,paired_matches,
        search_parameters->check_type,archive_search_end1->mm_allocator);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) {
    tab_global_inc();
    archive_search_pe_print(gem_log_get_stream(),archive_search_end1,archive_search_end2,paired_matches);
    tab_global_dec();
    tab_global_dec();
  }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH,PROFILE_LEVEL);
}

