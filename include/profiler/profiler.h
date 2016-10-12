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

#ifndef PROFILER_H_
#define PROFILER_H_

#include "system/commons.h"
#include "profiler/profiler_gem.h"
#include "profiler/profiler_cuda.h"
#include "profiler/profiler_vtune.h"

#ifdef GEM_PROFILE /* GEM_PROFILE ENABLED */

// Profile Levels
#define PHIGH true
#define PMED  true
#define PLOW  true

// Profile Aggregated
#define PROFILE_START(label,level) \
  if (level) { \
    PROF_START(label); \
    PROFILE_CUDA_START(#label,label); \
  }
#define PROFILE_STOP(label,level) \
  if (level) { \
    PROF_STOP(label); \
    PROFILE_CUDA_STOP(); \
  }
#define PROFILE_PAUSE(label,level) \
  if (level) { \
    PROF_PAUSE(label); \
    PROFILE_CUDA_STOP(); \
  }
#define PROFILE_CONTINUE(label,level) \
  if (level) { \
    PROF_CONTINUE(label); \
    PROFILE_CUDA_START(#label,label); \
  }
#else /* GEM_PROFILE DISABLED */
  #define PROFILE_START(label,level)
  #define PROFILE_STOP(label,level)
  #define PROFILE_PAUSE(label,level)
  #define PROFILE_CONTINUE(label,level)
#endif /* GEM_PROFILE */
#endif /* PROFILER_H_ */
