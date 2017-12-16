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
 *   Filtering candidates module provides functions to process all
 *   bwt-encoded positions (originated from a candidate generation process)
 *   and compose them into filtering-regions (using decoded text coordinates)
 *   This "buffered" module operates in batches of filtering-candidates and
 *   makes use of GPU-buffers to offloads the decoding of positions to a GPU
 */

#ifndef FILTERING_CANDIDATES_PROCESS_BUFFERED_H_
#define FILTERING_CANDIDATES_PROCESS_BUFFERED_H_

#include "filtering/candidates/filtering_candidates.h"
#include "filtering/candidates/filtering_candidates_buffered.h"
#include "align/pattern/pattern.h"
#include "gpu/gpu_buffer_fmi_decode.h"

/*
 * Decode SA-Positions Buffered (from GPU-Buffer)
 */
void filtering_candidates_buffered_decode_sa_positions(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern);

/*
 * Decode sampled SA-Positions Buffered (from GPU-Buffer)
 */
void filtering_candidates_buffered_decode_sampled_sa_positions(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_offset_begin);

/*
 * Decode Text-Positions Buffered (from GPU-Buffer)
 */
void filtering_candidates_buffered_decode_text_positions(
    filtering_candidates_t* const filtering_candidates,
    filtering_candidates_buffered_t* const filtering_candidates_buffered,
    pattern_t* const pattern,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode,
    const uint64_t buffer_offset_begin);

/*
 * Process Candidates Buffered (from GPU-Buffer)
 */
void filtering_candidates_buffered_process_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const bool compose_alignment_regions);

#endif /* FILTERING_CANDIDATES_PROCESS_BUFFERED_H_ */
