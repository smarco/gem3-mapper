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

//TODO implementar un buffer solo para el SA.
//     poder hacer la conversion de SA por a REF pos
//     con una llamada independiente

#endif /* GPU_SA_PRIMITIVES_C_ */
