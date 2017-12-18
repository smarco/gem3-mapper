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
 *   Approximate-String-Matching (ASM) module encapsulating
 *   the basic stages in a step-wise approach. Mainly used to
 *   horizontally process a batch of searches (doing each stage
 *   for all searches at the same time before progressing to
 *   the next search stage)
 */

#ifndef APPROXIMATE_SEARCH_STEPWISE_H_
#define APPROXIMATE_SEARCH_STEPWISE_H_

#include "utils/essentials.h"
#include "approximate_search/approximate_search.h"
#include "gpu/gpu_buffer_fmi_ssearch.h"
#include "gpu/gpu_buffer_fmi_asearch.h"
#include "gpu/gpu_buffer_fmi_decode.h"
#include "gpu/gpu_buffer_kmer_filter.h"
#include "gpu/gpu_buffer_bpm_distance.h"
#include "gpu/gpu_buffer_bpm_align.h"

/*
 * AM Stepwise :: Region Profile
 */
void approximate_search_stepwise_region_profile_static_generate(
    approximate_search_t* const search);
void approximate_search_stepwise_region_profile_static_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);
void approximate_search_stepwise_region_profile_static_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);

void approximate_search_stepwise_region_profile_adaptive_generate(
    approximate_search_t* const search);
void approximate_search_stepwise_region_profile_adaptive_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);
void approximate_search_stepwise_region_profile_adaptive_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);

/*
 * AM Stepwise :: Decode Candidates
 */
void approximate_search_stepwise_decode_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
void approximate_search_stepwise_decode_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);

/*
 * AM Stepwise :: Kmer-filter
 */
void approximate_search_stepwise_kmer_filter_copy(
    approximate_search_t* const search,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter);
void approximate_search_stepwise_kmer_filter_retrieve(
    approximate_search_t* const search,
    gpu_buffer_kmer_filter_t* const gpu_buffer_kmer_filter);

/*
 * AM Stepwise :: BPM-Distance
 */
void approximate_search_stepwise_bpm_distance_copy(
    approximate_search_t* const search,
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance);
void approximate_search_stepwise_bpm_distance_retrieve(
    approximate_search_t* const search,
    gpu_buffer_bpm_distance_t* const gpu_buffer_bpm_distance);

/*
 * AM Stepwise :: BPM-Align
 */
void approximate_search_stepwise_bpm_align_update(
    approximate_search_t* const search);
void approximate_search_stepwise_bpm_align_copy(
    approximate_search_t* const search,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align);
void approximate_search_stepwise_bpm_align_retrieve(
    approximate_search_t* const search,
    gpu_buffer_bpm_align_t* const gpu_buffer_bpm_align,
    matches_t* const matches);
#endif /* APPROXIMATE_SEARCH_STEPWISE_H_ */
