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
 * DESCRIPTION: Basic includes for GEM
 */

#ifndef ESSENTIALS_H_
#define ESSENTIALS_H_

// Common constants/SysIncludes
#include "system/commons.h"

// Error/Msg handling module
#include "system/report.h"
#include "system/errors.h"

// GThreads
#include "system/gthread.h"

// File Manager
#include "system/fm.h"

// Memory Manager module
#include "system/mm.h"
#include "system/mm_slab.h"
#include "system/mm_stack.h"

// Basic Runtime
#include "system/gruntime.h"

// Basic Data Structures
#include "utils/vector.h"
#include "utils/priority_queue.h"
#include "utils/segmented_vector.h"
#include "utils/hash.h"
#include "utils/string_buffer.h"

// Basic Profiling
#include "profiler/profiler.h"

// Mapper Profiling Counters
#include "mapper/mapper_profile_counters.h"

#endif /* ESSENTIALS_H_ */
