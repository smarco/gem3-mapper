/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

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
