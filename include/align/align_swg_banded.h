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
 *   Smith-Waterman-Gotoh (SWG) alignment module
 */

#ifndef ALIGN_SWG_BANDED_H_
#define ALIGN_SWG_BANDED_H_

#include "utils/essentials.h"
#include "matches/align/match_align_dto.h"

/*
 * SWG Banded
 */
void align_swg_banded(
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    const bool begin_free,
    const bool end_free,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    mm_stack_t* const mm_stack);

#endif /* ALIGN_SWG_BANDED_H_ */
