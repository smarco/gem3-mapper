/* include all the index modules */
#include "gpu_fmi_index.h"
#include "gpu_sa_index.h"

#ifndef GPU_INDEX_MODULES_H_
#define GPU_INDEX_MODULES_H_

/*****************************
Global Objects (General)
*****************************/

typedef struct {
  gpu_fmi_buffer_t fmi;
  gpu_sa_buffer_t  sa;
  gpu_module_t     activeModules;
} gpu_index_buffer_t;

#endif /* GPU_INDEX_MODULES_H_ */
