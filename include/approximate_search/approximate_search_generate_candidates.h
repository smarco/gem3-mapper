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
 *   Approximate-String-Matching (ASM) module to generate filtering-candidates
 *   from a region-profile/key-partition
 */


#ifndef APPROXIMATE_SEARCH_GENERATE_CANDIDATES_H_
#define APPROXIMATE_SEARCH_GENERATE_CANDIDATES_H_

#include "approximate_search/approximate_search.h"
#include "gpu/gpu_buffer_fmi_decode.h"

/*
 * Generate Candidates from Region-Profile
 */
void approximate_search_generate_candidates_limit_exact_matches(
    approximate_search_t* const search);
void approximate_search_generate_candidates(
    approximate_search_t* const search);

/*
 * Buffered Copy/Retrieve
 */
void approximate_search_generate_candidates_buffered_copy(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);
void approximate_search_generate_candidates_buffered_retrieve(
    approximate_search_t* const search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode);

#endif /* APPROXIMATE_SEARCH_GENERATE_CANDIDATES_H_ */
