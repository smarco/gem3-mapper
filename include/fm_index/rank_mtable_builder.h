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
 *   FM-Index module provides a builder for a rank-memoization table
 *   that stores pre-computed bwt-intevals for all the k-mer (k=11)
 *   words from {ACGT}
 */

#ifndef RANK_MTABLE_BUILDER_H_
#define RANK_MTABLE_BUILDER_H_

#include "utils/essentials.h"
#include "fm_index/bwt/bwt.h"
#include "fm_index/rank_mtable.h"

/*
 * Write mtable
 */
void rank_mtable_builder_write(
    fm_t* const file_manager,
    rank_mtable_t* const rank_mtable);

/*
 * Generation
 */
rank_mtable_t* rank_mtable_builder_new(
    const bwt_builder_t* const bwt_builder,
    const bool verbose);

/*
 * Delete
 */
void rank_mtable_builder_delete(rank_mtable_t* const rank_mtable);

#endif /* RANK_MTABLE_BUILDER_H_ */
