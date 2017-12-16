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

#ifndef ARCHIVE_SEARCH_PE_STEPWISE_H_
#define ARCHIVE_SEARCH_PE_STEPWISE_H_

#include "utils/essentials.h"
#include "archive/search/archive_search.h"
#include "gpu/gpu_buffer_fmi_ssearch.h"
#include "gpu/gpu_buffer_fmi_asearch.h"
#include "gpu/gpu_buffer_fmi_decode.h"
#include "gpu/gpu_buffer_kmer_filter.h"
#include "gpu/gpu_buffer_bpm_distance.h"
#include "gpu/gpu_buffer_bpm_align.h"
#include "matches/matches.h"
#include "matches/paired_matches.h"

/*
 * Stepwise: Init Search
 */
void archive_search_pe_stepwise_init_search(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);

/*
 * Stepwise: Region-Profile
 */
void archive_search_pe_stepwise_region_profile_static_generate(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);
void archive_search_pe_stepwise_region_profile_static_copy(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);
void archive_search_pe_stepwise_region_profile_static_retrieve(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);

void archive_search_pe_stepwise_region_profile_adaptive_generate(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);
void archive_search_pe_stepwise_region_profile_adaptive_copy(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);
void archive_search_pe_stepwise_region_profile_adaptive_retrieve(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);

/*
 * Stepwise: Decode-Candidates
 */
void archive_search_pe_stepwise_decode_copy(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
void archive_search_pe_stepwise_decode_retrieve(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);

/*
 * Stepwise: Kmer-filter
 */
void archive_search_pe_stepwise_kmer_filter_copy(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter);
void archive_search_pe_stepwise_kmer_filter_retrieve(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter);

/*
 * Stepwise: BPM-Distance
 */
void archive_search_pe_stepwise_bpm_distance_copy(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance);
void archive_search_pe_stepwise_bpm_distance_retrieve(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance,
    matches_t* const matches);

/*
 * Stepwise: BPM-Align
 */
void archive_search_pe_stepwise_bpm_align_update(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);
void archive_search_pe_stepwise_bpm_align_copy(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);
void archive_search_pe_stepwise_bpm_align_retrieve(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    matches_t* const matches);

/*
 * Stepwise: Finish Search
 */
void archive_search_pe_stepwise_finish_search(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

#endif /* ARCHIVE_SEARCH_PE_STEPWISE_H_ */
