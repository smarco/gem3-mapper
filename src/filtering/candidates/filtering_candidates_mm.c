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

#include "filtering/candidates/filtering_candidates_mm.h"

/*
 * Filtering Candidates MM
 */
void filtering_candidates_mm_init(
    filtering_candidates_mm_t* const filtering_candidates_mm) {
  // MM
  filtering_candidates_mm->mm_slab = mm_slab_new_(
      BUFFER_SIZE_4M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM,"filtering_candidates_mm.4MB/8MB");
  mm_stack_t* const mm_general = mm_stack_new(filtering_candidates_mm->mm_slab);
  mm_stack_t* const mm_filtering = mm_stack_new(filtering_candidates_mm->mm_slab);
  filtering_candidates_mm->mm_general = mm_general;
  filtering_candidates_mm->mm_filtering = mm_filtering;
  // Specialized Stacks
  filtering_candidates_mm->mm_text = mm_filtering;
  filtering_candidates_mm->mm_positions = mm_filtering;
  filtering_candidates_mm->mm_regions = mm_filtering;
  filtering_candidates_mm->mm_alignment_tiles = mm_filtering;
  filtering_candidates_mm->mm_alignment_regions = mm_filtering;
}
void filtering_candidates_mm_clear(
    filtering_candidates_mm_t* const filtering_candidates_mm) {
  mm_stack_clear(filtering_candidates_mm->mm_general);
  mm_stack_clear(filtering_candidates_mm->mm_filtering);
}
void filtering_candidates_mm_destroy(
    filtering_candidates_mm_t* const filtering_candidates_mm) {
  mm_stack_delete(filtering_candidates_mm->mm_general);
  mm_stack_delete(filtering_candidates_mm->mm_filtering);
  mm_slab_delete(filtering_candidates_mm->mm_slab);
}
void filtering_candidates_mm_inject_buffered_mm(
    filtering_candidates_mm_t* const filtering_candidates_mm,
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm) {
  filtering_candidates_mm->mm_alignment_tiles = filtering_candidates_buffered_mm->mm_regions;
  filtering_candidates_mm->mm_alignment_regions = filtering_candidates_buffered_mm->mm_regions;
}
/*
 * Filtering Candidates Buffered MM
 */
void filtering_candidates_buffered_mm_init(
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm) {
  filtering_candidates_buffered_mm->mm_slab = mm_slab_new_(
      BUFFER_SIZE_1M,BUFFER_SIZE_4M,MM_UNLIMITED_MEM,"filtering_candidates_buffered_mm.1MB/2MB");
  filtering_candidates_buffered_mm->mm_positions = mm_stack_new(filtering_candidates_buffered_mm->mm_slab);
  filtering_candidates_buffered_mm->mm_regions = mm_stack_new(filtering_candidates_buffered_mm->mm_slab);
  filtering_candidates_buffered_mm->mm_regions_attr = mm_stack_new(filtering_candidates_buffered_mm->mm_slab);
}
void filtering_candidates_buffered_mm_clear(
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm) {
  mm_stack_clear(filtering_candidates_buffered_mm->mm_positions);
  mm_stack_clear(filtering_candidates_buffered_mm->mm_regions);
  mm_stack_clear(filtering_candidates_buffered_mm->mm_regions_attr);
}
void filtering_candidates_buffered_mm_clear_positions(
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm) {
  mm_stack_clear(filtering_candidates_buffered_mm->mm_positions);
}
void filtering_candidates_buffered_mm_clear_regions(
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm) {
  mm_stack_clear(filtering_candidates_buffered_mm->mm_regions);
}
void filtering_candidates_buffered_mm_clear_regions_attr(
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm) {
  mm_stack_clear(filtering_candidates_buffered_mm->mm_regions_attr);
}
void filtering_candidates_buffered_mm_destroy(
    filtering_candidates_buffered_mm_t* const filtering_candidates_buffered_mm) {
  mm_stack_delete(filtering_candidates_buffered_mm->mm_positions);
  mm_stack_delete(filtering_candidates_buffered_mm->mm_regions);
  mm_stack_delete(filtering_candidates_buffered_mm->mm_regions_attr);
  mm_slab_delete(filtering_candidates_buffered_mm->mm_slab);
}
