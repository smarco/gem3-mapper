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
 *   Approximate-String-Matching (ASM) module to produce
 *   region-profiles (i.e. key-partition with few times candidates)
 */

#ifndef APPROXIMATE_SEARCH_REGION_PROFILE_H_
#define APPROXIMATE_SEARCH_REGION_PROFILE_H_

#include "approximate_search/approximate_search.h"
#include "gpu/gpu_buffer_fmi_ssearch.h"
#include "gpu/gpu_buffer_fmi_asearch.h"

/*
 * Region Profile Adaptive
 */
void approximate_search_region_profile(approximate_search_t* const search);

/*
 * Region Partition Fixed
 */
void approximate_search_region_profile_static_partition(approximate_search_t* const search);
void approximate_search_region_profile_static_compute(
    approximate_search_t* const search);

/*
 * Static Buffered Copy/Retrieve
 */
void approximate_search_region_profile_static_buffered_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);
void approximate_search_region_profile_static_buffered_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch);
void approximate_search_region_profile_static_buffered_recompute(
    approximate_search_t* const search);

/*
 * Adaptive Buffered Copy/Retrieve
 */
void approximate_search_region_profile_adaptive_buffered_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);
void approximate_search_region_profile_adaptive_buffered_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch);

#endif /* APPROXIMATE_SEARCH_REGION_PROFILE_H_ */
