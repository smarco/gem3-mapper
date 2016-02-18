/*
 * PROJECT: GEMMapper
 * FILE: profiler.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple time/functional profiler module
 */

#ifndef PROFILER_H_
#define PROFILER_H_

#include "system/commons.h"
#include "system/profiler_gem.h"
#include "system/profiler_cuda.h"

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
