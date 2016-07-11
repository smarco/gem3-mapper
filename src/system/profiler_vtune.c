/*
 * PROJECT: GEMMapper
 * FILE: profiler_cuda.c
 * DATE: 06/06/2012
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: CUDA profile module
 */

#include "system/profiler_vtune.h"

/*
 * CUDA Support
 */
#ifdef GEM_VTUNE

#include "ittnotify.h"

void PROFILE_VTUNE_START() {
  __itt_resume();
}
void PROFILE_VTUNE_STOP() {
  __itt_pause();
}
#else
void PROFILE_VTUNE_START() {}
void PROFILE_VTUNE_STOP() {}
#endif
