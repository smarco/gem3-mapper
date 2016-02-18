/*
 * PROJECT: GEMMapper
 * FILE: gpu_config.c
 * DATE: 04/09/2014
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "gpu/gpu_config.h"
#include "resources/gpu_modules/gpu_interface.h"

/*
 * CUDA Supported
 */
#ifdef HAVE_CUDA

bool gpu_supported() {
  return (gpu_get_num_supported_devices_(GPU_ARCH_SUPPORTED) > 0);
}

/*
 * CUDA NOT-Supported
 */
#else

bool gpu_supported() { return false; }

#endif /* HAVE_CUDA */
