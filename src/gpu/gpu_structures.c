/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2013-2017 by Alejandro Chacon <alejandro.chacond@gmail.com>
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
 *            Alejandro Chacon <alejandro.chacond@gmail.com>
 * DESCRIPTION:
 *   GPU-adaptor module provides functions to write a GEM3-GPU index
 */

#include "gpu/gpu_structures.h"
#include "gpu/gpu_config.h"

/*
 * CUDA Supported
 */
#ifdef HAVE_CUDA

/*
 * GPU Index Write
 */
void gpu_structures_write(
    const char* const index_file_name_prefix,
    dna_text_t* const enc_text,
    const uint64_t forward_text_length,
    bwt_builder_t* const bwt_builder,
    rank_mtable_t* const rank_mtable,
    uint64_t* const sa_gem,
    const uint64_t sa_sampling) {
  // Prepare ticket
  ticker_t ticker;
  ticker_percentage_reset(&ticker,true,"Building GPU-Index",0,0,true);
  // Configure index name
  char* const gpu_index_name = gem_strcat(index_file_name_prefix,".gem.gpu");
  // Prepare GPU DTOs (Data Transfer Objects)
  gpu_gem_fmi_dto_t gpu_gem_fmi_dto = {
      .c                        = bwt_builder->bwt.c,
      .C                        = bwt_builder->bwt.C,
      .mayor_counters           = bwt_builder->bwt.mayor_counters,
      .bwt_mem                  = bwt_builder->bwt.bwt_mem,
      .bwt_length               = bwt_builder->bwt.length,
      .index_coding             = GPU_INDEX_GEM_FULL,
      /* (numLevels, Size) => (9, 0.8MB) (10, 3.2MB) (11, 13MB) (12, 52MB) (13, 210MB) */
      .num_levels_fmi_table     = RANK_MTABLE_LEVELS,
      .skip_levels_fmi_table    = rank_mtable->min_matching_depth,
      .occ_threashold_fmi_table = RANK_MTABLE_MMD_THRESHOLD,
  };
  gpu_gem_sa_dto_t gpu_gem_sa_dto = {
      .sa           = sa_gem,
      .sa_sampling  = (uint32_t)sa_sampling,
      .index_coding = GPU_INDEX_GEM_FULL,
      .sa_length    = bwt_builder->bwt.length,
  };
  gpu_gem_ref_dto_t gpu_gem_ref_dto = {
      .reference  = (char*) dna_text_get_text(enc_text),
      .ref_coding = GPU_REF_GEM_FULL,
      .ref_length = forward_text_length,
  };
  // Translate & write GPU index
  gpu_io_save_indexed_structures_GEM_(gpu_index_name,
      &gpu_gem_fmi_dto,&gpu_gem_ref_dto,&gpu_gem_sa_dto,GPU_ALL_MODULES);
  // Free
  free(gpu_index_name);
  // Finish ticker
  ticker_finish(&ticker);
}

/*
 * CUDA NOT-Supported
 */
#else

/*
 * GPU Structures Write
 */
void gpu_structures_write(
    const char* const index_file_name_prefix,
    dna_text_t* const enc_text,
    const uint64_t forward_text_length,
    bwt_builder_t* const bwt_builder,
    rank_mtable_t* const rank_mtable,
    uint64_t* const sa_gem,
    const uint64_t sa_sampling) { GEM_CUDA_NOT_SUPPORTED(); }

#endif /* HAVE_CUDA */
