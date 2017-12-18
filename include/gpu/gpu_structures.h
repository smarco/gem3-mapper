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
 *   GPU-adaptor module provides functions to write a GEM3-GPU index
 */

#ifndef GPU_STRUCTURES_H_
#define GPU_STRUCTURES_H_

#include "utils/essentials.h"
#include "fm_index/bwt/bwt.h"
#include "fm_index/rank_mtable.h"
#include "text/dna_text.h"

/*
 * GPU Structures Write
 */
void gpu_structures_write(
    const char* const index_file_name_prefix,
    dna_text_t* const enc_text,
    const uint64_t forward_text_length,
    bwt_builder_t* const bwt_builder,
    rank_mtable_t* const rank_mtable,
    uint64_t* const sa_gem,
    const uint64_t sa_sampling);

#endif /* GPU_STRUCTURES_H_ */
