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
 *   FM-Index module provides a rank-memoization table that stores
 *   pre-computed bwt-intevals for all the k-mer (k=11) words from {ACGT}
 */

#ifndef RANK_MTABLE_H_
#define RANK_MTABLE_H_

#include "utils/essentials.h"

/*
 * Constants
 */
#define RANK_MTABLE_SEARCH_DEPTH  11                          // Number of character that can be searched up in the table
#define RANK_MTABLE_LEVELS       (RANK_MTABLE_SEARCH_DEPTH+1) // One fake level (zero-HI)

#define RANK_MTABLE_MMD_THRESHOLD 20                          // Minimum Matching Depth (MMD)

#define RANK_MTABLE_SIZE(search_depth) ( (5*(POW4(search_depth)-1))/3 + 2 )  /* S := (5*(4^(n)-1))/3 + 2 */
#define RANK_MTABLE_LEVEL_SIZE(level)  (5*POW4(level-1))

typedef struct {
  // Meta-info
  uint64_t table_size;         // Total number of ranks stored
  uint64_t num_levels;         // Total depth of the table
  // Table
  uint64_t* level_skip;        // Skip from levels-to-level (Pre-computed)
  uint64_t** sa_ranks_levels;  // Pointers to the levels
  // Optimization info
  uint64_t min_matching_depth; // Minimum depth to achieve less than MMD_THRESHOLD matches
  /* MM */
  mm_t* mm_sa_ranks;
} rank_mtable_t;

typedef struct {
  uint64_t hi_position; // Effective HI-position in the table
  uint64_t level;       // Level on the table
} rank_mquery_t;

/*
 * Setup
 */
void rank_mtable_init_levels(rank_mtable_t* const rank_mtable);
rank_mtable_t* rank_mtable_read(fm_t* const file_manager);
rank_mtable_t* rank_mtable_read_mem(mm_t* const memory_manager);
void rank_mtable_delete(rank_mtable_t* const rank_mtable);

/*
 * Accessors
 */
uint64_t rank_mtable_get_size(const rank_mtable_t* const rank_mtable);

/*
 * Query
 */
void rank_mquery_new(rank_mquery_t* const query);
void rank_mquery_add_char(
    const rank_mtable_t* const rank_mtable,
    rank_mquery_t* const query,
    uint8_t const enc_char);
uint64_t rank_mquery_get_level(const rank_mquery_t* const query);
uint64_t rank_mquery_is_exhausted(const rank_mquery_t* const query);

/*
 * Fetch rank value
 */
void rank_mtable_prefetch(
    const rank_mtable_t* const rank_mtable,
    const rank_mquery_t* const query);
void rank_mtable_fetch(
    const rank_mtable_t* const rank_mtable,
    const rank_mquery_t* const query,
    uint64_t* const lo,
    uint64_t* const hi);

/*
 * Display
 */
void rank_mtable_print(
    FILE* const stream,
    rank_mtable_t* const rank_mtable,
    const bool print_content);

#endif /* RANK_MTABLE_H_ */
