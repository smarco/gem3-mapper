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
 *   Mapper-profile module provides functions to display the
 *   GEM-profiling counters involved in the CUDA-mapping workflow
 */

#ifndef MAPPER_PROFILE_CUDA_H_
#define MAPPER_PROFILE_CUDA_H_

#include "utils/essentials.h"
#include "mapper/mapper_profile_counters.h"

/*
 * Mapper SE
 */
void mapper_profile_print_mapper_se_cuda(
    FILE* const stream,const bool map_output,const uint64_t num_threads);

/*
 * Mapper PE
 */
void mapper_profile_print_mapper_pe_cuda(
    FILE* const stream,const bool map_output,const uint64_t num_threads);

#endif /* MAPPER_PROFILE_CUDA_H_ */
