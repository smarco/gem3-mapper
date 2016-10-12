/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2011-2017 by Alejandro Chacon  <alejandro.chacon@uab.es>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 *            Alejandro Chacon  <alejandro.chacon@uab.es>
 */

#include "profiler/profiler_cuda.h"

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
void PROFILE_CUDA_STOP(void) {
  nvtxRangePop();
}
#else
  void PROFILE_CUDA_START(char* const name,const uint64_t cid) {}
  void PROFILE_CUDA_STOP(void) {}
#endif /* GEM_PROFILE */
/*
 * CUDA NOT-Supported
 */
#else
  void PROFILE_CUDA_START(char* const name,const uint64_t cid) {}
  void PROFILE_CUDA_STOP() {}
#endif
