/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
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
 * DESCRIPTION:
 *   BWT data structure (basic implementation sampled).
 *   Provides basic routines to encode a DNA text into a
 *   Burrows-Wheeler transform using a compact bit representation
 *   and counter buckets as to enhance Occ/rank queries
 */

#ifndef BWT_H_
#define BWT_H_

/*
 * BWT-Models
 */
#include "fm_index/bwt/bwt_sbasic.h"
#include "fm_index/bwt/bwt_basic.h"

/*
 * Instantiate BWT Models
 */
// Forward
typedef bwt_sbasic_t bwt_t;
typedef bwt_sbasic_builder_t bwt_builder_t;

/*
 * BWT Builder
 */
// Forward
#define bwt_builder_new bwt_sbasic_builder_new
#define bwt_builder_write bwt_sbasic_builder_write
#define bwt_builder_delete bwt_sbasic_builder_delete

/*
 * BWT Loader
 */
// Forward
#define bwt_read_mem bwt_sbasic_read_mem
#define bwt_delete bwt_sbasic_delete

/*
 * BWT General Accessors
 */
// Forward
#define bwt_builder_get_length bwt_sbasic_builder_get_length
#define bwt_builder_get_size bwt_sbasic_builder_get_size
#define bwt_get_length bwt_sbasic_get_length
#define bwt_get_size bwt_sbasic_get_size
#define bwt_is_same_bucket bwt_sbasic_is_same_bucket

/*
 * BWT Character Accessors
 */
// Forward
#define bwt_char bwt_sbasic_char
#define bwt_char_character bwt_sbasic_char_character

/*
 * BWT ERank (Exclusive Rank Function)
 */
// Forward
#define bwt_builder_erank bwt_sbasic_builder_erank
#define bwt_erank bwt_sbasic_erank
#define bwt_erank_character bwt_sbasic_erank_character
#define bwt_erank_interval bwt_sbasic_erank_interval
#define bwt_sampling_erank bwt_sbasic_sampling_erank

/*
 * BWT Prefetched ERank
 */
// Forward
#define bwt_prefetch bwt_sbasic_prefetch
#define bwt_prefetched_erank bwt_sbasic_prefetched_erank
#define bwt_prefetched_erank_interval bwt_sbasic_prefetched_erank_interval

/*
 * BWT Precomputed ERank (Precomputation of the block's elements)
 */
// Forward
#define bwt_precompute bwt_sbasic_precompute
#define bwt_precompute_interval bwt_sbasic_precompute_interval
#define bwt_prefetched_precompute bwt_sbasic_prefetched_precompute
#define bwt_prefetched_precompute_interval bwt_sbasic_prefetched_precompute_interval
#define bwt_precomputed_erank bwt_sbasic_precomputed_erank
#define bwt_precomputed_erank_interval bwt_sbasic_precomputed_erank_interval

/*
 * BWT LF (Last to first)
 */
// Forward
#define bwt_LF_ bwt_sbasic_LF_
#define bwt_LF bwt_sbasic_LF
#define bwt_prefetched_LF bwt_sbasic_prefetched_LF
#define bwt_LF__enc bwt_sbasic_LF__enc
#define bwt_LF__character bwt_sbasic_LF__character
#define bwt_prefetched_LF__enc bwt_sbasic_prefetched_LF__enc

/*
 * Display
 */
// Forward
#define bwt_builder_print bwt_sbasic_builder_print
#define bwt_print bwt_sbasic_print

#endif /* BWT_H_ */
