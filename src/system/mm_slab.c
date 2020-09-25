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
 *   SlabMemory. Relative big amounts of contiguous memory allocated
 *   all at once (like the LINUX slab allocator). Objects can be
 *   allocated inside the slab reducing the overhead of malloc/free
 *   cycles along the program
 */

#include "system/mm_slab.h"
#include "system/gthread.h"

/*
 * Debug
 */
#define MM_SLAB_DEBUG_REQUESTED_SEGMENTS  false
#define MM_SLAB_LOG                       false

/*
 * Constants
 */
#define MM_SLAB_NUM_INITIAL_SEGMENTS 100

/*
 * Setup
 */
void mm_slab_add_new_segment(mm_slab_t* const mm_slab) {
  // Allocate new slab-segment
  mm_slab_segment_t* const mm_slab_segment = mm_alloc(mm_slab_segment_t);
  mm_slab_segment->segment_id = (mm_slab->segment_id_generator)++; // Set ID
  // Calculate number of slab units
  const uint64_t slab_unit_size = mm_slab->slab_unit_size;
  const uint64_t total_slabs_units = mm_slab->slab_segment_size/slab_unit_size;
  mm_slab_segment->total_slabs_units = total_slabs_units;
  mm_slab_segment->busy_slabs_units = 0;
  // Allocate mm-segment
  mm_t* mem_manager;
//  const uint64_t available_memory = (mm_slab->max_memory - mm_slab->requested_memory);
//  if ( (available_memory > mm_slab->slab_segment_size) &&
//       (mm_get_available_mem() > mm_slab->slab_segment_size+BUFFER_SIZE_128M) ) {
    mem_manager = mm_bulk_malloc(mm_slab->slab_segment_size,false);
    mm_slab->requested_memory += mm_slab->slab_segment_size;
//  } else {
//    mem_manager = mm_bulk_mmalloc_temp(mm_slab->slab_segment_size);
//  }
  mm_slab_segment->mm = mem_manager;
  // Add the segment
  vector_insert(mm_slab->slabs_segments,mm_slab_segment,mm_slab_segment_t*);
  // Allocate the slab units
  mm_slab_segment->slab_units = mm_calloc(total_slabs_units,mm_slab_unit_t,false);
  // Add the units as free
  int64_t i;
  void* memory = mm_get_mem(mem_manager)+((total_slabs_units-1)*slab_unit_size);
  for (i=total_slabs_units-1;i>=0;--i,memory-=slab_unit_size) {
    mm_slab_unit_t* const mm_slab_unit = mm_slab_segment->slab_units+i;
    mm_slab_unit->slab_segment = mm_slab_segment;
    mm_slab_unit->memory = memory;
    vector_insert(mm_slab->slabs_units_free,mm_slab_unit,mm_slab_unit_t*);
  }
  gem_cond_log(MM_SLAB_DEBUG_REQUESTED_SEGMENTS,
      "[Slab:%"PRIu64"]Has allocated %"PRIu64" MB",mm_slab->slab_id,
      (vector_get_used(mm_slab->slabs_segments)*mm_slab->slab_segment_size)/1024/1024);
  gem_cond_log(MM_SLAB_LOG,"[GEM]> mm_slab(%"PRIu64").addSegment(%"PRIu64" MB)",
      mm_slab->slab_id,CONVERT_B_TO_MB(mm_slab->slab_segment_size));
}
mm_slab_t* mm_slab_new_(
    const uint64_t slab_unit_size,
    const uint64_t slab_group_size,
    const uint64_t max_allocatable_memory) {
  gem_cond_fatal_error(slab_group_size < slab_unit_size,MM_SLAB_WRONG_DIMENSIONS,slab_group_size,slab_unit_size);
  mm_slab_t* const mm_slab = mm_alloc(mm_slab_t);
  mm_slab->slab_id = gem_rand_IID(0,UINT16_MAX);
  gem_cond_log(MM_SLAB_LOG,"[GEM]> mm_slab(%"PRIu64").new()",mm_slab->slab_id);
  mm_slab->max_memory = max_allocatable_memory;
  mm_slab->requested_memory = 0;
  mm_slab->segment_id_generator = 0; // Init ids
  // Set page size
  const int64_t sz = sysconf(_SC_PAGESIZE);
  gem_cond_fatal_error(sz==-1,SYS_SYSCONF);
  // Set sizes (always multiples of the page size)
  const uint64_t segment_pages = ((slab_group_size+(sz-1))/sz); // SysPages
  mm_slab->slab_segment_size = segment_pages*sz;
  mm_slab->slab_unit_size = slab_unit_size;
  // gem_cond_error_msg(mm_slab->slab_segment_size%mm_slab->slab_unit_size!=0,
  //   SLAB_WASTED_MEM,mm_slab->slab_segment_size,mm_slab->slab_unit_size);
  // Allocate vectors
  mm_slab->slabs_segments = vector_new(MM_SLAB_NUM_INITIAL_SEGMENTS,mm_slab_segment_t*);
  mm_slab->slabs_units_free = vector_new(MM_SLAB_NUM_INITIAL_SEGMENTS,mm_slab_unit_t*);
  // Allocate init segment
  mm_slab_add_new_segment(mm_slab);
  // Return
  return mm_slab;
}
void mm_slab_delete(
    mm_slab_t* const mm_slab) {
  // Free all slab segments
  VECTOR_ITERATE(mm_slab->slabs_segments,slabs_segment,ss_i,mm_slab_segment_t*) {
    mm_bulk_free((*slabs_segment)->mm); // Free memory
    mm_free((*slabs_segment)->slab_units); // Free all slab units
    mm_free(*slabs_segment);
  }
  vector_delete(mm_slab->slabs_segments);
  vector_delete(mm_slab->slabs_units_free);
  mm_free(mm_slab);
}
/*
 *  Accessors
 */
mm_slab_unit_t* mm_slab_request(
    mm_slab_t* const mm_slab) {
  // Check slabs available (Add new one if required)
  if (gem_expect_false(vector_get_used(mm_slab->slabs_units_free)==0)) {
    mm_slab_add_new_segment(mm_slab);
  }
  // Serve free slab
  mm_slab_unit_t* const mm_slab_unit = *vector_get_last_elm(mm_slab->slabs_units_free,mm_slab_unit_t*);
  vector_dec_used(mm_slab->slabs_units_free); // Remove from slabs available
  ++(mm_slab_unit->slab_segment->busy_slabs_units); // Decrement free slabs
  return mm_slab_unit;
}
void mm_slab_return(
    mm_slab_t* const mm_slab,
    mm_slab_unit_t* const mm_slab_unit) {
  // Restore slab as free
  --(mm_slab_unit->slab_segment->busy_slabs_units);
  vector_insert(mm_slab->slabs_units_free,mm_slab_unit,mm_slab_unit_t*);
}
uint64_t mm_slab_get_slab_size(
    mm_slab_t* const mm_slab) {
  return mm_slab->slab_unit_size;
}
/*
 * Display/Profile
 */
void mm_slab_print(
    FILE* const stream,
    mm_slab_t* const mm_slab,
    const bool show_internals) {
  GEM_CHECK_NULL(stream);
  const uint64_t num_segments = vector_get_used(mm_slab->slabs_segments);
  const uint64_t slab_units_per_segment = mm_slab->slab_segment_size/mm_slab->slab_unit_size;
  const uint64_t slab_units_total = num_segments*slab_units_per_segment;
  const uint64_t slab_units_free = vector_get_used(mm_slab->slabs_units_free);
  // General Header
  fprintf(stream,"Slab.totalSize %"PRIu64" bytes\n",CONVERT_B_TO_MB(num_segments*mm_slab->slab_segment_size));
  // Display dimensions
  fprintf(stream,"  --> Slab.segment.total %"PRIu64"\n",num_segments);
  fprintf(stream,"  --> Slab.segment.size %"PRIu64"\n",mm_slab->slab_segment_size);
  fprintf(stream,"    --> Slab.unit.size %"PRIu64" bytes\n",mm_slab->slab_unit_size);
  fprintf(stream,"    --> Slab.unit.perSegment %"PRIu64" bytes\n",slab_units_per_segment);
  fprintf(stream,"    --> Slab.unit.total %"PRIu64"\n",slab_units_total);
  fprintf(stream,"      --> Slab.unit.free %"PRIu64"  (%02.4f%%)\n",
      slab_units_free,PERCENTAGE(slab_units_free,slab_units_total));
  if (show_internals) {
    fprintf(stream,">>> Slab.segments.internals\n");
    fprintf(stream,"ID \t Busy(%%busy)/Total \t SlabUnits \t MemoryAddrs\n");
    VECTOR_ITERATE(mm_slab->slabs_segments,slabs_segment,ss_i,mm_slab_segment_t*) {
      fprintf(stream,"%"PRIu64"\t%"PRIu64"(%02.4f%%)/%"PRIu64"\t%"PRIu64"\t%p\n",
          (*slabs_segment)->segment_id,(*slabs_segment)->busy_slabs_units,
          PERCENTAGE((*slabs_segment)->busy_slabs_units,(*slabs_segment)->total_slabs_units),
          (*slabs_segment)->total_slabs_units,mm_get_allocated((*slabs_segment)->mm),mm_get_mem((*slabs_segment)->mm));
    }
    fprintf(stream,">>> Slab.units.internals\n");
    fprintf(stream,"ID \t MemoryAddrs\n");
    VECTOR_ITERATE(mm_slab->slabs_units_free,slabs_unit,su_i,mm_slab_unit_t*) {
      fprintf(stream,"%"PRIu64" \t %p",(*slabs_unit)->slab_segment->segment_id,(*slabs_unit)->memory);
    }
  }
}
