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
 *   Archive module to build the Suffix-Array (SA) and generate
 *   the Burrows-Wheeler-transform (BWT) & Ferragina-Manzini Index (FM-Index)
 */

#ifndef ARCHIVE_BUILDER_INDEX_H_
#define ARCHIVE_BUILDER_INDEX_H_

#include "utils/essentials.h"
#include "archive/builder/archive_builder.h"

/*
 * Build BWT (SA)
 */
void archive_builder_index_build_bwt(
    archive_builder_t* const archive_builder,
    const bool gpu_index,
    const bool dump_bwt,
    const bool dump_explicit_sa,
    const bool verbose);

/*
 * Display
 */
void archive_builder_index_print_explicit_sa(
    archive_builder_t* const archive_builder,
    const char* const extension);
void archive_builder_index_print_bwt(
    archive_builder_t* const archive_builder,
    const char* const extension,
    const bool verbose);

#endif /* ARCHIVE_BUILDER_INDEX_H_ */
