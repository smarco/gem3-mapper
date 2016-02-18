/*
 * PROJECT: GEMMapper
 * FILE: essentials.h
 * DATE: 06/06/2012
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
#include "system/mm_pool.h"

// Basic Profiling
#include "system/profiler.h"

// Basic Runtime
#include "system/gruntime.h"

// Basic Data Structures
#include "utils/vector.h"
#include "utils/priority_queue.h"
#include "utils/segmented_vector.h"
#include "utils/hash.h"
#include "utils/string_buffer.h"

// Mapper Profiling Counters
#include "mapper/mapper_profile_counters.h"

#endif /* ESSENTIALS_H_ */
