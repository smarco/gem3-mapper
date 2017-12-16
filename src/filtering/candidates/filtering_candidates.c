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
 *   Filtering candidates module provides functions to store and handle
 *   all positions/regions during a search based on filtering, that is,
 *   generation of candidates & verification of candidates
 */

#include "align/alignment.h"
#include "filtering/candidates/filtering_candidates.h"
#include "filtering/region/filtering_region.h"
#include "filtering/region/filtering_region_verify.h"
#include "filtering/region/filtering_region_align.h"
#include "matches/classify/matches_classify.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES  GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Constants
 */
#define REGIONS_BUFFER_INIT                      100
#define CANDIDATE_POSITIONS_INIT                 1000

/*
 * Setup
 */
void filtering_candidates_init(filtering_candidates_t* const filtering_candidates) {
  // Candidates
  filtering_candidates->filtering_positions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_position_t*);
  filtering_candidates->filtering_regions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_region_t*);
  filtering_candidates->discarded_regions = vector_new(CANDIDATE_POSITIONS_INIT,filtering_region_t*);
  // Cache
  filtering_region_cache_init(&filtering_candidates->filtering_region_cache);
  // Stats Histogram
  for(uint64_t id_cand_histo = 0; id_cand_histo < GEM_HIST_CAND_ALIGNED + 1; id_cand_histo++){
    COUNTER_RESET(&filtering_candidates->candidates_aligned_histo[id_cand_histo]);
  }
}
void filtering_candidates_init_alignment(
    filtering_candidates_t* const filtering_candidates,
    alignment_t* const alignment,
    pattern_t* const pattern,
    const uint64_t text_length,
    const uint64_t max_error) {
  // Check alignment
  if (alignment->alignment_tiles!=NULL) return; // Already initialized
  // Allocate & initialize alignment
  const uint64_t num_tiles = pattern->pattern_tiled.num_tiles;
  alignment->alignment_tiles =
      filtering_candidates_allocate_alignment_tiles(
          filtering_candidates,num_tiles);
  alignment_init(
      alignment,pattern->key_length,0,text_length,max_error,
      num_tiles,pattern->pattern_tiled.tile_length);
}
void filtering_candidates_clear(
    filtering_candidates_t* const filtering_candidates,
    const bool free_memory) {
  // Candidate Positions
  if (free_memory) {
    VECTOR_ITERATE(filtering_candidates->filtering_positions,position,p,filtering_position_t*) {
      filtering_candidates_free_position(filtering_candidates,*position);
    }
  }
  vector_clear(filtering_candidates->filtering_positions);
  // Candidate Regions
  if (free_memory) {
    VECTOR_ITERATE(filtering_candidates->filtering_regions,region,r,filtering_region_t*) {
      filtering_candidates_free_region(filtering_candidates,*region);
    }
  }
  vector_clear(filtering_candidates->filtering_regions);
  // Discarded Regions
  if (free_memory) {
    VECTOR_ITERATE(filtering_candidates->discarded_regions,dis_region,dr,filtering_region_t*) {
      filtering_candidates_free_region(filtering_candidates,*dis_region);
    }
  }
  vector_clear(filtering_candidates->discarded_regions);
}
void filtering_candidates_destroy(
    filtering_candidates_t* const filtering_candidates,
    const bool free_memory) {
  // Clear structures
  filtering_candidates_clear(filtering_candidates,free_memory);
  // Candidates
  vector_delete(filtering_candidates->filtering_positions);
  vector_delete(filtering_candidates->filtering_regions);
  vector_delete(filtering_candidates->discarded_regions);
}
/*
 * Handlers Injection (Support Data Structures)
 */
void filtering_candidates_inject_handlers(
    filtering_candidates_t* const filtering_candidates,
    archive_t* const archive,
    search_parameters_t* const search_parameters,
    mm_allocator_t* const mm_allocator) {
  // Inject Hanlders
  filtering_candidates->archive = archive;
  filtering_candidates->search_parameters = search_parameters;
  filtering_candidates->mm_allocator = mm_allocator;
  // Clear structures
  filtering_candidates_clear(filtering_candidates,false);
}
/*
 * Allocators
 */
filtering_position_t* filtering_candidates_allocate_position(
    filtering_candidates_t* const filtering_candidates) {
  filtering_position_t* const position = mm_allocator_malloc(
      filtering_candidates->mm_allocator,sizeof(filtering_position_t));
  vector_insert(filtering_candidates->filtering_positions,position,filtering_position_t*);
  return position;
}
void filtering_candidates_free_position(
    const filtering_candidates_t* const filtering_candidates,
    filtering_position_t* const filtering_position) {
  mm_allocator_free(filtering_candidates->mm_allocator,filtering_position);
}
filtering_region_t* filtering_candidates_allocate_region(
    filtering_candidates_t* const filtering_candidates) {
  filtering_region_t* const region = mm_allocator_malloc(
      filtering_candidates->mm_allocator,sizeof(filtering_region_t));
  vector_insert(filtering_candidates->filtering_regions,region,filtering_region_t*);
  return region;
}
void filtering_candidates_free_region(
    const filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region) {
  // Free text trace
  text_trace_destroy(&filtering_region->text_trace,filtering_candidates->mm_allocator);
  // Free alignment tiles
  filtering_candidates_free_alignment_tiles(
      filtering_candidates,filtering_region->alignment.alignment_tiles);
  // Free scaffold
  match_scaffold_destroy(
      &filtering_region->match_scaffold,
      filtering_candidates->mm_allocator);
  // Free handler
  mm_allocator_free(filtering_candidates->mm_allocator,filtering_region);
}
alignment_tile_t* filtering_candidates_allocate_alignment_tiles(
    filtering_candidates_t* const filtering_candidates,
    const uint64_t num_alignment_tiles) {
  return mm_allocator_malloc(
      filtering_candidates->mm_allocator,
      num_alignment_tiles*sizeof(alignment_tile_t));
}
void filtering_candidates_free_alignment_tiles(
    const filtering_candidates_t* const filtering_candidates,
    alignment_tile_t* const alignment_tile) {
  if (alignment_tile!=NULL) {
    mm_allocator_free(filtering_candidates->mm_allocator,alignment_tile);
  }
}
/*
 * Filtering Positions
 */
uint64_t filtering_candidates_get_num_positions(
    const filtering_candidates_t* const filtering_candidates) {
  return vector_get_used(filtering_candidates->filtering_positions);
}
void filtering_candidates_set_num_positions(
    const filtering_candidates_t* const filtering_candidates,
    const uint64_t num_positions) {
  vector_set_used(filtering_candidates->filtering_positions,num_positions);
}
filtering_position_t** filtering_candidates_get_positions(
    const filtering_candidates_t* const filtering_candidates) {
  return vector_get_mem(filtering_candidates->filtering_positions,filtering_position_t*);
}
void filtering_candidates_clear_positions(
    const filtering_candidates_t* const filtering_candidates,
    const bool free_positions) {
  if (free_positions) {
    VECTOR_ITERATE(filtering_candidates->filtering_positions,position,p,filtering_position_t*) {
      filtering_candidates_free_position(filtering_candidates,*position);
    }
  }
  vector_clear(filtering_candidates->filtering_positions);
}
/*
 * Filtering Regions
 */
uint64_t filtering_candidates_get_num_regions(
    const filtering_candidates_t* const filtering_candidates) {
  return vector_get_used(filtering_candidates->filtering_regions);
}
void filtering_candidates_set_num_regions(
    const filtering_candidates_t* const filtering_candidates,
    const uint64_t num_regions) {
  vector_set_used(filtering_candidates->filtering_regions,num_regions);
}
filtering_region_t** filtering_candidates_get_regions(
    const filtering_candidates_t* const filtering_candidates) {
  return vector_get_mem(filtering_candidates->filtering_regions,filtering_region_t*);
}
void filtering_candidates_clear_regions(
    const filtering_candidates_t* const filtering_candidates,
    const bool free_regions) {
  if (free_regions) {
    VECTOR_ITERATE(filtering_candidates->filtering_regions,region,r,filtering_region_t*) {
      filtering_candidates_free_region(filtering_candidates,*region);
    }
  }
  vector_clear(filtering_candidates->filtering_regions);
}
/*
 * Discarded Regions
 */
uint64_t filtering_candidates_get_num_discarded_regions(
    const filtering_candidates_t* const filtering_candidates) {
  return vector_get_used(filtering_candidates->discarded_regions);
}
void filtering_candidates_set_num_discarded_regions(
    const filtering_candidates_t* const filtering_candidates,
    const uint64_t num_discarded_regions) {
  vector_set_used(filtering_candidates->discarded_regions,num_discarded_regions);
}
void filtering_candidates_add_num_discarded_regions(
    const filtering_candidates_t* const filtering_candidates,
    const uint64_t num_discarded_regions) {
  vector_add_used(filtering_candidates->discarded_regions,num_discarded_regions);
}
filtering_region_t** filtering_candidates_get_discarded_regions(
    const filtering_candidates_t* const filtering_candidates) {
  return vector_get_mem(filtering_candidates->discarded_regions,filtering_region_t*);
}
filtering_region_t** filtering_candidates_reserve_discarded_regions(
    const filtering_candidates_t* const filtering_candidates,
    const uint64_t num_regions) {
  vector_reserve_additional(filtering_candidates->discarded_regions,num_regions);
  return vector_get_free_elm(filtering_candidates->discarded_regions,filtering_region_t*);
}
void filtering_candidates_clear_discarded_regions(
    const filtering_candidates_t* const filtering_candidates,
    const bool free_regions) {
  if (free_regions) {
    VECTOR_ITERATE(filtering_candidates->discarded_regions,dis_region,dr,filtering_region_t*) {
      filtering_candidates_free_region(filtering_candidates,*dis_region);
    }
  }
  vector_clear(filtering_candidates->discarded_regions);
}
/*
 * Sorting
 */
int64_t filtering_position_cmp_position(filtering_position_t** a,filtering_position_t** b) {
  return (int64_t)(*a)->text_end_position - (int64_t)(*b)->text_end_position;
}
int64_t filtering_region_cmp_align_distance(filtering_region_t** a,filtering_region_t** b) {
  const int64_t cmp = (int64_t)(*a)->alignment.distance_min_bound - (int64_t)(*b)->alignment.distance_min_bound;
  if (cmp) return cmp;
  // Compare position (Helps stability)
  return (int64_t)(*a)->text_begin_position - (int64_t)(*b)->text_begin_position;
}
int64_t filtering_region_cmp_scaffold_coverage(filtering_region_t** a,filtering_region_t** b) {
  return (int64_t)(*b)->match_scaffold.scaffolding_coverage - (int64_t)(*a)->match_scaffold.scaffolding_coverage;
}
int64_t filtering_region_cmp_align_rank(filtering_region_t** a,filtering_region_t** b) {
  return (int64_t)(*b)->alignment.distance_rank - (int64_t)(*a)->alignment.distance_rank;
}
#define VECTOR_SORT_NAME                 filtering_positions
#define VECTOR_SORT_TYPE                 filtering_position_t*
#define VECTOR_SORT_CMP(a,b)             filtering_position_cmp_position(a,b)
#include "utils/vector_sort.h"
void filtering_candidates_sort_positions(filtering_candidates_t* const filtering_candidates) {
  PROF_ADD_COUNTER(GP_FC_SORT_BY_POSITION,vector_get_used(filtering_candidates->filtering_positions));
  vector_sort_filtering_positions(filtering_candidates->filtering_positions);
}
#define VECTOR_SORT_NAME                 align_distance
#define VECTOR_SORT_TYPE                 filtering_region_t*
#define VECTOR_SORT_CMP(a,b)             filtering_region_cmp_align_distance(a,b)
#include "utils/vector_sort.h"
void filtering_candidates_sort_regions_by_align_distance(filtering_candidates_t* const filtering_candidates) {
  PROF_ADD_COUNTER(GP_FC_SORT_BY_ALIGN_DIST,vector_get_used(filtering_candidates->filtering_regions));
  vector_sort_align_distance(filtering_candidates->filtering_regions);
}
#define VECTOR_SORT_NAME                 scaffold_coverage
#define VECTOR_SORT_TYPE                 filtering_region_t*
#define VECTOR_SORT_CMP(a,b)             filtering_region_cmp_scaffold_coverage(a,b)
#include "utils/vector_sort.h"
void filtering_candidates_sort_regions_by_scaffold_coverage(filtering_candidates_t* const filtering_candidates) {
  PROF_ADD_COUNTER(GP_FC_SORT_BY_COVERAGE,vector_get_used(filtering_candidates->filtering_regions));
  vector_sort_scaffold_coverage(filtering_candidates->filtering_regions);
}
void filtering_candidates_sort_discarded_by_scaffold_coverage(filtering_candidates_t* const filtering_candidates) {
  PROF_ADD_COUNTER(GP_FC_SORT_BY_COVERAGE,vector_get_used(filtering_candidates->discarded_regions));
  vector_sort_scaffold_coverage(filtering_candidates->discarded_regions);
}
#define VECTOR_SORT_NAME                 align_rank
#define VECTOR_SORT_TYPE                 filtering_region_t*
#define VECTOR_SORT_CMP(a,b)             filtering_region_cmp_align_rank(a,b)
#include "utils/vector_sort.h"
void filtering_candidates_sort_discarded_by_rank(filtering_candidates_t* const filtering_candidates) {
  PROF_ADD_COUNTER(GP_FC_SORT_BY_COVERAGE,vector_get_used(filtering_candidates->discarded_regions));
  vector_sort_align_rank(filtering_candidates->discarded_regions);
}
/*
 * Display
 */
void filtering_candidates_print_regions_by_status(
    FILE* const stream,
    vector_t* const filtering_regions,
    const filtering_region_status_t status,
    const bool print_alignment_regions) {
  uint64_t i, total_printed = 0;
  const uint64_t num_regions = vector_get_used(filtering_regions);
  filtering_region_t** const fregion = vector_get_mem(filtering_regions,filtering_region_t*);
  // Count
  for (i=0;i<num_regions;++i) {
    if (fregion[i]->status!=status) continue;
    ++total_printed;
  }
  if (total_printed == 0) return;
  tab_fprintf(stream,"  => Regions.%s  (%"PRIu64")\n",filtering_region_status_label[status],total_printed);
  // Print
  tab_global_inc();
  for (i=0;i<num_regions;++i) {
    if (fregion[i]->status!=status) continue;
    filtering_region_print(stream,fregion[i],false,print_alignment_regions,true);
  }
  tab_global_dec();
}
void filtering_candidates_print_regions(
    FILE* const stream,
    filtering_candidates_t* const filtering_candidates,
    const bool print_alignment_regions) {
  tab_fprintf(stream,"[GEM]>Filtering.Regions\n");
  vector_t* const filtering_regions = filtering_candidates->filtering_regions;
  vector_t* const discarded_regions = filtering_candidates->discarded_regions;
  filtering_candidates_print_regions_by_status(
      stream,filtering_regions,filtering_region_pending,print_alignment_regions);
  filtering_candidates_print_regions_by_status(
      stream,filtering_regions,filtering_region_unverified,print_alignment_regions);
  filtering_candidates_print_regions_by_status(
      stream,discarded_regions,filtering_region_verified_discarded,print_alignment_regions);
  filtering_candidates_print_regions_by_status(
      stream,filtering_regions,filtering_region_accepted,print_alignment_regions);
  filtering_candidates_print_regions_by_status(
      stream,discarded_regions,filtering_region_accepted_subdominant,print_alignment_regions);
  const uint64_t total_regions =
      filtering_candidates_get_num_regions(filtering_candidates) +
      vector_get_used(filtering_candidates->discarded_regions);
  if (total_regions > 0) tab_fprintf(stream,"  => Total.Regions %"PRIu64"\n",total_regions);
}
void filtering_candidates_print_positions(
    FILE* const stream,
    filtering_candidates_t* const filtering_candidates) {
  const uint64_t num_candidate_positions = filtering_candidates_get_num_positions(filtering_candidates);
  filtering_position_t** const candidate_positions = filtering_candidates_get_positions(filtering_candidates);
  uint64_t i;
  tab_fprintf(stream,"[GEM]>Filtering.Positions\n");
  for (i=0;i<num_candidate_positions;++i) {
    tab_fprintf(stream,"  => Key[%lu,%lu) ~> Text[%lu,%lu) Distance=%lu\n",
        candidate_positions[i]->source_region_begin,
        candidate_positions[i]->source_region_end,
        candidate_positions[i]->text_begin_position,
        candidate_positions[i]->text_end_position,
        candidate_positions[i]->align_distance);
  }
}
