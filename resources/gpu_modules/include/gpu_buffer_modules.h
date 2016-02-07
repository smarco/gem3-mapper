/* Include all the supported modules */
#include "gpu_fmi_primitives.h"
#include "gpu_sa_primitives.h"
#include "gpu_bpm_primitives.h"

#ifndef GPU_BUFFER_MODULES_H_
#define GPU_BUFFER_MODULES_H_

typedef union{
  gpu_bpm_buffer_t        bpm;
  gpu_fmi_search_buffer_t search;
  gpu_fmi_decode_buffer_t decode;
} gpu_buffer_modules_t;

#endif /* GPU_BUFFER_MODULES_H_ */
