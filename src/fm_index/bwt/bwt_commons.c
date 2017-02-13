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
 *   BWT common constants & tables
 */

#include "fm_index/bwt/bwt_commons.h"

/*
 * Profile
 */
uint64_t _bwt_ranks = 0; // Bwt rank counter

/*
 * XOR table to mask bitmap depending based on the character (enc)
 */
const int64_t xor_table_1[] = {-1ll, -1ll, -1ll, -1ll, -0ll, -0ll, -0ll, -0ll};
const int64_t xor_table_2[] = {-1ll, -1ll, -0ll, -0ll, -1ll, -1ll, -0ll, -0ll};
const int64_t xor_table_3[] = {-1ll, -0ll, -1ll, -0ll, -1ll, -0ll, -1ll, -0ll};
