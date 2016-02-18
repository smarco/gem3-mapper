/*
 * PROJECT: GEMMapper
 * FILE: profiler_cuda.c
 * DATE: 06/06/2012
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: CUDA profile module
 */

#include "system/profiler_cuda.h"

/*
 * CUDA Support
 */
#ifdef HAVE_CUDA
#include "nvToolsExt.h"

const uint32_t profiler_cuda_tags_colors[] = {
    0x0000ff00,
    0x000000ff,
    0x00ffff00,
    0x00ff00ff,
    0x0000ffff,
    0x00ff0000,
    0x00ffffff
};

/*
 * GEM_PROFILE
 */
#ifdef GEM_PROFILE
void PROFILE_CUDA_START(char* const name,const uint64_t cid) {
  nvtxEventAttributes_t event_attr = {0};
  event_attr.version = NVTX_VERSION;
  event_attr.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  event_attr.colorType = NVTX_COLOR_ARGB;
  event_attr.color = profiler_cuda_tags_colors[cid%PROFILER_CUDA_TAGS_NUM_COLORS];
  event_attr.messageType = NVTX_MESSAGE_TYPE_ASCII;
  event_attr.message.ascii = name;
  nvtxRangePushEx(&event_attr);
}
void PROFILE_CUDA_STOP() {
  nvtxRangePop();
}
#else
  void PROFILE_CUDA_START(char* const name,const uint64_t cid) {}
  void PROFILE_CUDA_STOP() {}
#endif /* GEM_PROFILE */
/*
 * CUDA NOT-Supported
 */
#else
  void PROFILE_CUDA_START(char* const name,const uint64_t cid) {}
  void PROFILE_CUDA_STOP() {}
#endif
