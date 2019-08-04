/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_SA_PRIMITIVES_C_
#define GPU_SA_PRIMITIVES_C_

#include "../include/gpu_sa_primitives.h"


gpu_sa_entry_t* gpu_sa_buffer_get_index_(const void* const saBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) saBuffer;
  return(mBuff->index->sa.h_sa);
}

gpu_sa_decode_text_pos_t* gpu_sa_decode_buffer_get_ref_pos_(const void* const saBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) saBuffer;
  return(mBuff->data.decode.textPositions.h_textPos);
}

#endif /* GPU_SA_PRIMITIVES_C_ */
