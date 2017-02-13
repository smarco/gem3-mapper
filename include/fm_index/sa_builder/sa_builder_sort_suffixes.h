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
 */

#ifndef SA_BUILDER_SORT_SUFFIXES_H_
#define SA_BUILDER_SORT_SUFFIXES_H_

#include "utils/essentials.h"
#include "fm_index/sa_builder/sa_builder.h"

/*
 * Sorting Suffixes
 *   1.- Count all suffixes
 *   2.- Store all suffixes
 *   3.- Sort all suffixes
 */
void sa_builder_sort_suffixes(
    sa_builder_t* const sa_builder,
    dna_text_t* const enc_bwt,
    sampled_sa_builder_t* const sampled_sa,
    const bool verbose);

#endif /* SA_BUILDER_SORT_SUFFIXES_H_ */
