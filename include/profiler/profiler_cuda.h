/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2011-2017 by Alejandro Chacon  <alejandro.chacon@uab.es>
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
 *            Alejandro Chacon  <alejandro.chacon@uab.es>
 */

#ifndef PROFILE_CUDA_H_
#define PROFILE_CUDA_H_

/*
 * Include
 */
#include "system/commons.h"

/*
 * Tags Colors
 */
#define PROFILER_CUDA_TAGS_NUM_COLORS 7
extern const uint32_t profiler_cuda_tags_colors[];

/*
 * Profile Start/Stop
 */
void PROFILE_CUDA_START(char* const name,const uint64_t cid);
void PROFILE_CUDA_STOP(void);

#endif /* PROFILE_CUDA_H_ */
