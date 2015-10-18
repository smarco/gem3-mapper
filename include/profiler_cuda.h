/*
 *
 */

#ifndef PROFILE_CUDA_H_
#define PROFILE_CUDA_H_

#include "commons.h"

#ifdef GEM_NOPROFILE

#define PROFILE_CUDA_START(name,cid)
#define PROFILE_CUDA_STOP

#else
 
#include "nvToolsExt.h"

extern uint32_t colors[];
extern int num_colors;

#define PROFILE_CUDA_START(name,cid) { \
	int color_id = cid; \
	color_id = color_id%num_colors;\
	nvtxEventAttributes_t eventAttrib = {0}; \
	eventAttrib.version = NVTX_VERSION; \
	eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
	eventAttrib.colorType = NVTX_COLOR_ARGB; \
	eventAttrib.color = colors[color_id]; \
	eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII; \
	eventAttrib.message.ascii = name; \
	nvtxRangePushEx(&eventAttrib); \
}
#define PROFILE_CUDA_STOP nvtxRangePop();

#endif

#endif /* PROFILE_CUDA_H_ */
