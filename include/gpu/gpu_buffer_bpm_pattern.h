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
 */

#ifndef GPU_BUFFER_BPM_PATTERN_H_
#define GPU_BUFFER_BPM_PATTERN_H_

#include "utils/essentials.h"
#include "align/pattern/pattern.h"
#include "gpu/gpu_config.h"

/*
 * Constants
 */
#define GPU_BPM_ALPHABET_LENGTH  GPU_BPM_FILTER_PEQ_ALPHABET_SIZE
#define GPU_BPM_ENTRY_LENGTH     GPU_BPM_FILTER_PEQ_ENTRY_LENGTH
#define GPU_BPM_ENTRY_SIZE       (GPU_BPM_FILTER_PEQ_ENTRY_LENGTH / UINT8_LENGTH)
#define GPU_BPM_NUM_SUB_ENTRIES  GPU_BPM_FILTER_PEQ_SUBENTRIES
#define GPU_BPM_SUBENTRY_LENGTH  GPU_BPM_FILTER_PEQ_SUBENTRY_LENGTH

/*
 * Aliased query-entry
 */
typedef struct {
  uint32_t bitmap[GPU_BPM_ALPHABET_LENGTH][GPU_BPM_NUM_SUB_ENTRIES];
} gpu_bpm_peq_entry_t;

/*
 * Compile/Decompile Query Entries [DTO]
 */
void gpu_buffer_bpm_pattern_compile(
    gpu_bpm_peq_entry_t* const pattern_entry,
    const bpm_pattern_t* const bpm_pattern);
void gpu_buffer_bpm_pattern_decompile(
    gpu_bpm_peq_entry_t* const pattern_entry,
    const uint32_t pattern_length,
    bpm_pattern_t* const bpm_pattern,
    mm_allocator_t* const mm_allocator);

#endif /* GPU_BUFFER_BPM_PATTERN_H_ */

