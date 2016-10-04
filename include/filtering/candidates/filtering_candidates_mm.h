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
 *   Memory management (MM) module for filtering candidates that
 *   centralizes all memory sources (mm_stacks) and specializes
 *   them according to its expected use
 */

#ifndef FILTERING_CANDIDATES_MM_H_
#define FILTERING_CANDIDATES_MM_H_

#include "utils/essentials.h"

/*
 * Filtering Candidates MM
 */
typedef struct {
  // MM
  mm_slab_t* mm_slab;                   // MM-Slab
  mm_stack_t* mm_general;               // MM-Stack
  mm_stack_t* mm_filtering;             // MM-Stack
  // Specialized Stacks
  mm_stack_t* mm_text;                  // MM-Stack
  mm_stack_t* mm_positions;             // MM-Stack
  mm_stack_t* mm_regions;               // MM-Stack
  mm_stack_t* mm_alignment_tiles;       // MM-Stack
  mm_stack_t* mm_alignment_regions;     // MM-Stack
} filtering_candidates_mm_t;

/*
 * Filtering Candidates Buffered MM
 */
typedef struct {
  mm_slab_t* mm_slab;                   // MM-Slab
  mm_stack_t* mm_positions;             // MM-Stack
  mm_stack_t* mm_regions;               // MM-Stack
  mm_stack_t* mm_regions_attr;          // MM-Stack
} filtering_candidates_buffered_mm_t;

/*
 * Filtering Candidates MM
 */
void filtering_candidates_mm_init(
    filtering_candidates_mm_t* const filtering_candidates_mm);
void filtering_candidates_mm_clear(
    filtering_candidates_mm_t* const filtering_candidates_mm);
void filtering_candidates_mm_destroy(
    filtering_candidates_mm_t* const filtering_candidates_mm);
void filtering_candidates_mm_inject_buffered_mm(
    filtering_candidates_mm_t* const filtering_candidates_mm,
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm);

/*
 * Filtering Candidates Buffered MM
 */
void filtering_candidates_buffered_mm_init(
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm);
void filtering_candidates_buffered_mm_clear(
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm);
void filtering_candidates_buffered_mm_clear_positions(
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm);
void filtering_candidates_buffered_mm_clear_regions(
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm);
void filtering_candidates_buffered_mm_clear_regions_attr(
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm);
void filtering_candidates_buffered_mm_destroy(
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm);

#endif /* FILTERING_CANDIDATES_MM_H_ */
