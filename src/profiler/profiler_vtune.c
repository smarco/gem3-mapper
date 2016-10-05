/*
 * PROJECT: GEMMapper
 * FILE: profiler_cuda.c
 * DATE: 06/06/2012
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: CUDA profile module
 */

#include "profiler/profiler_vtune.h"

/*
 * CUDA Support
 */
#ifdef GEM_VTUNE

#include "ittnotify.h"

void PROFILE_VTUNE_START(void) {
  __itt_resume();
}
void PROFILE_VTUNE_STOP(void) {
  __itt_pause();
}
#else
void PROFILE_VTUNE_START(void) {}
void PROFILE_VTUNE_STOP(void) {}
#endif
