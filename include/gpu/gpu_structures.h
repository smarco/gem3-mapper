/*
 * PROJECT: GEMMapper
 * FILE: gpu_structures.h
 * DATE: 06/06/2012
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#ifndef GPU_STRUCTURES_H_
#define GPU_STRUCTURES_H_

#include "utils/essentials.h"
#include "fm_index/bwt.h"
#include "data_structures/dna_text.h"

/*
 * GPU Structures Write
 */
void gpu_structures_write(
    const char* const index_file_name_prefix,
    dna_text_t* const enc_text,
    const uint64_t forward_text_length,
    bwt_builder_t* const bwt_builder,
    uint64_t* const sa_gem,
    const uint32_t sa_sampling);

#endif /* GPU_STRUCTURES_H_ */
