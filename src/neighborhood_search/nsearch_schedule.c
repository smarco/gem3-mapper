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
 */

#include "neighborhood_search/nsearch_schedule.h"
#include "neighborhood_search/nsearch_operation.h"
#include "neighborhood_search/nsearch_partition.h"
#include "neighborhood_search/nsearch_hamming.h"
#include "neighborhood_search/nsearch_levenshtein.h"
#include "neighborhood_search/nsearch_levenshtein_scheduled.h"
#include "neighborhood_search/nsearch_filtering.h"

/*
 * Setup
 */
void nsearch_schedule_init(
    nsearch_schedule_t* const nsearch_schedule,
    const nsearch_model_t nsearch_model,
    const uint64_t max_complete_error,
    archive_t* const archive,
    pattern_t* const pattern,
    region_profile_t* const region_profile,
    search_parameters_t* const search_parameters,
    filtering_candidates_t* const filtering_candidates,
    matches_t* const matches) {
  // Search Structures
  nsearch_schedule->search_id = 0;
  nsearch_schedule->archive = archive;
  nsearch_schedule->pattern = pattern;
  nsearch_schedule->region_profile = region_profile;
  nsearch_schedule->filtering_candidates = filtering_candidates;
  nsearch_schedule->matches = matches;
  // Search Parameters
  nsearch_schedule->search_parameters = search_parameters;
  nsearch_schedule->nsearch_model = nsearch_model;
  nsearch_schedule->max_error = max_complete_error;
  nsearch_schedule->quick_abandon = false;
  // Search Operations
  const uint64_t key_length = nsearch_schedule->pattern->key_length;
  const uint64_t max_error = nsearch_schedule->max_error;
  const uint64_t max_text_length = key_length + max_error;
  const uint64_t max_pending_ops = (uint64_t)ceil(gem_log2((float)(max_error+1))) + 1;
  nsearch_operation_t* const pending_searches =
      mm_allocator_calloc(nsearch_schedule->mm_allocator,max_pending_ops,nsearch_operation_t,false);
  nsearch_schedule->pending_searches = pending_searches;
  nsearch_schedule->num_pending_searches = 0;
  uint64_t i;
  for (i=0;i<max_pending_ops;++i) {
    nsearch_operation_init(pending_searches+i,
        key_length,max_text_length,nsearch_schedule->mm_allocator);
  }
  // Profiler
  nsearch_schedule->profile.ns_nodes = 0;
  nsearch_schedule->profile.ns_nodes_success = 0;
  nsearch_schedule->profile.ns_nodes_fail = 0;
}
void nsearch_schedule_inject_mm(
    nsearch_schedule_t* const nsearch_schedule,
    mm_allocator_t* const mm_allocator) {
  nsearch_schedule->mm_allocator = mm_allocator;
}
/*
 * Add pending search
 */
bool nsearch_schedule_add_pending_search(
    nsearch_schedule_t* const nsearch_schedule,
    const search_direction_t search_direction,
    const uint64_t local_key_begin,
    const uint64_t local_key_length,
    const uint64_t global_key_begin,
    const uint64_t global_key_length,
    const uint64_t min_local_error,
    const uint64_t min_global_error,
    const uint64_t max_global_error) {
  // Check impossible search configuration
  if (min_local_error > local_key_length) return false;
  if (min_global_error > max_global_error) return false;
  // Queue search
  const uint64_t num_op = nsearch_schedule->num_pending_searches;
  nsearch_operation_t* const nsearch_operation = nsearch_schedule->pending_searches + num_op;
  nsearch_operation->search_direction = search_direction;
  nsearch_operation->local_key_begin = local_key_begin;
  nsearch_operation->local_key_end = local_key_begin+local_key_length;
  nsearch_operation->global_key_begin = global_key_begin;
  nsearch_operation->global_key_end = global_key_begin+global_key_length;
  nsearch_operation->min_local_error = min_local_error;
  nsearch_operation->min_global_error = min_global_error;
  nsearch_operation->max_global_error = max_global_error;
  ++(nsearch_schedule->num_pending_searches);
  // Prepare DP
  const bool search_forward = (search_direction==direction_forward);
  const bool supercondensed = (search_forward) ?
          (nsearch_operation->global_key_begin == 0) :
          (nsearch_operation->global_key_end == nsearch_schedule->pattern->key_length);
  nsearch_levenshtein_state_prepare(&nsearch_operation->nsearch_state,supercondensed);
  // Return
  return true;
}
/*
 * Schedule the search
 */
void nsearch_schedule_search_step(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t chunk_offset,
    const uint64_t chunk_length,
    const uint64_t chunk_min_error,
    const uint64_t chunk_max_error) {
  // Adjust error
  bool feasible_search;
  if (chunk_max_error==0 || chunk_length<=chunk_max_error) {
    // Exact search of the chunk
    feasible_search = nsearch_schedule_add_pending_search(
        nsearch_schedule,direction_forward,
        chunk_offset,chunk_length,chunk_offset,chunk_length,
        chunk_min_error,chunk_min_error,chunk_max_error);
    if (!feasible_search) return; // Impossible search
    // Update current mcs
#ifndef NSEARCH_ENUMERATE
    const uint64_t current_mcs = nsearch_schedule_compute_min_error(nsearch_schedule);
    matches_update_mcs(nsearch_schedule->matches,current_mcs);
#endif
    // Select nsearch-alignment model & Perform the search (Solve pending extensions)
    switch (nsearch_schedule->nsearch_model) {
      case nsearch_model_hamming: {
        nsearch_hamming_scheduled_search(nsearch_schedule);
        break;
      }
      case nsearch_model_levenshtein: {
        nsearch_levenshtein_scheduled_search(nsearch_schedule);
        break;
      }
      default:
        GEM_INVALID_CASE();
        break;
    }
  } else {
    const uint64_t num_pending_searches = nsearch_schedule->num_pending_searches;
    // Compute error partition
    nsearch_partition_t epartition;
    nsearch_partition_compute(&epartition,chunk_offset,chunk_length);
    nsearch_partition_compute_error(&epartition,chunk_min_error,chunk_max_error);
    /*
     * First Partition (forward)
     */
    feasible_search = nsearch_schedule_add_pending_search(
        nsearch_schedule,direction_forward,
        epartition.offset_1,epartition.length_1,chunk_offset,chunk_length,
        epartition.extend_1_local_min_error,epartition.extend_1_global_min_error,
        epartition.extend_1_global_max_error);
    if (!feasible_search) return; // Impossible search
    nsearch_schedule_search_step(
        nsearch_schedule,epartition.offset_0,epartition.length_0,
        epartition.search_0_min_error,epartition.search_0_max_error);
    if (nsearch_schedule->quick_abandon) return; // Quick abandon
    /*
     * Second partition (backward)
     */
    nsearch_schedule->num_pending_searches = num_pending_searches; // Restore pending searches point
    feasible_search = nsearch_schedule_add_pending_search(
        nsearch_schedule,direction_backward,
        epartition.offset_0,epartition.length_0,chunk_offset,chunk_length,
        epartition.extend_0_local_min_error,epartition.extend_0_global_min_error,
        epartition.extend_0_global_max_error);
    if (!feasible_search) return; // Impossible search
    nsearch_schedule_search_step(
        nsearch_schedule,epartition.offset_1,epartition.length_1,
        epartition.search_1_min_error,epartition.search_1_max_error);
    // Unnecessary but => if (nsearch_schedule->quick_abandon) return; // Quick abandon
  }
}
void nsearch_schedule_search(nsearch_schedule_t* const nsearch_schedule) {
  PROF_START(GP_NS_GENERATION);
  // Search
  nsearch_schedule_search_step(
      nsearch_schedule,0,nsearch_schedule->pattern->key_length,
      0,nsearch_schedule->max_error);
  // Adjust MCS
#ifndef NSEARCH_ENUMERATE
  if (!nsearch_schedule->quick_abandon) {
    matches_update_mcs(nsearch_schedule->matches,nsearch_schedule->max_error+1);
  }
#endif
  // Profile
  PROF_ADD_COUNTER(GP_NS_NODES,nsearch_schedule->profile.ns_nodes);
  PROF_ADD_COUNTER(GP_NS_NODES_SUCCESS,nsearch_schedule->profile.ns_nodes_success);
  PROF_ADD_COUNTER(GP_NS_NODES_FAIL,nsearch_schedule->profile.ns_nodes_fail);
  PROF_STOP(GP_NS_GENERATION);
}
/*
 * Schedule the search (preconditioned by region profile)
 */
void nsearch_schedule_search_preconditioned_step(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t region_offset,
    const uint64_t num_regions,
    const uint64_t chunk_min_error,
    const uint64_t chunk_max_error) {
  // PRE: (num_regions > 1)
  region_profile_t* const region_profile = nsearch_schedule->region_profile;
  region_search_t* const regions = region_profile->filtering_region;
  const uint64_t global_region_offset = regions[region_offset].begin;
  const uint64_t global_region_length = regions[region_offset+num_regions-1].end - regions[region_offset].begin;
  if (chunk_max_error==0 || global_region_length<=chunk_max_error) {
    nsearch_schedule_search_step(nsearch_schedule,
        global_region_offset,global_region_length,chunk_min_error,chunk_max_error);
    if (nsearch_schedule->quick_abandon) return; // Quick abandon
  } else {
    bool feasible_search;
    const uint64_t num_pending_searches = nsearch_schedule->num_pending_searches;
    // Error Partition
    nsearch_partition_t epartition;
    nsearch_partition_preconditioned_compute(&epartition,regions,region_offset,num_regions);
    nsearch_partition_compute_error(&epartition,chunk_min_error,chunk_max_error);
    /*
     * Fist Partition
     */
    // Extend first partition (forward)
    feasible_search = nsearch_schedule_add_pending_search(
        nsearch_schedule,direction_forward,
        epartition.offset_1,epartition.length_1,global_region_offset,global_region_length,
        epartition.extend_1_local_min_error,epartition.extend_1_global_min_error,
        epartition.extend_1_global_max_error);
    if (!feasible_search) return; // Impossible search
    // Search first partition (forward)
    if (epartition.region_0!=NULL) {
      nsearch_schedule_search_step(
          nsearch_schedule,epartition.offset_0,epartition.length_0,
          epartition.search_0_min_error,epartition.search_0_max_error);
    } else {
      nsearch_schedule_search_preconditioned_step(
          nsearch_schedule,epartition.region_offset_0,epartition.num_regions_0,
          epartition.search_0_min_error,epartition.search_0_max_error);
    }
    if (nsearch_schedule->quick_abandon) return; // Quick abandon
    /*
     * Second Partition
     */
    // Extend second partition (backward)
    nsearch_schedule->num_pending_searches = num_pending_searches; // Restore pending searches point
    feasible_search = nsearch_schedule_add_pending_search(
        nsearch_schedule,direction_backward,
        epartition.offset_0,epartition.length_0,global_region_offset,global_region_length,
        epartition.extend_0_local_min_error,epartition.extend_0_global_min_error,
        epartition.extend_0_global_max_error);
    if (!feasible_search) return; // Impossible search
    // Search second partition (backward)
    if (epartition.region_1!=NULL) {
      nsearch_schedule_search_step(
          nsearch_schedule,epartition.offset_1,epartition.length_1,
          epartition.search_1_min_error,epartition.search_1_max_error);
    } else {
      nsearch_schedule_search_preconditioned_step(
          nsearch_schedule,epartition.region_offset_1,epartition.num_regions_1,
          epartition.search_1_min_error,epartition.search_1_max_error);
    }
    // Unnecessary but => if (nsearch_schedule->quick_abandon) return; // Quick abandon
  }
}
void nsearch_schedule_search_preconditioned(nsearch_schedule_t* const nsearch_schedule) {
  PROF_START(GP_NS_GENERATION);
  // Search
  region_profile_t* const region_profile = nsearch_schedule->region_profile;
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
  if (num_filtering_regions <= 1) {
    nsearch_schedule_search_step(nsearch_schedule,
        0,nsearch_schedule->pattern->key_length,0,nsearch_schedule->max_error);
  } else {
    nsearch_schedule_search_preconditioned_step(nsearch_schedule,
        0,num_filtering_regions,0,nsearch_schedule->max_error);
  }
  // Adjust MCS
#ifndef NSEARCH_ENUMERATE
  if (!nsearch_schedule->quick_abandon) {
    matches_update_mcs(nsearch_schedule->matches,nsearch_schedule->max_error+1);
  }
#endif
  // Profile
  PROF_ADD_COUNTER(GP_NS_NODES,nsearch_schedule->profile.ns_nodes);
  PROF_ADD_COUNTER(GP_NS_NODES_SUCCESS,nsearch_schedule->profile.ns_nodes_success);
  PROF_ADD_COUNTER(GP_NS_NODES_FAIL,nsearch_schedule->profile.ns_nodes_fail);
  PROF_STOP(GP_NS_GENERATION);
}
/*
 * Utils
 */
uint64_t nsearch_schedule_compute_min_error(
    nsearch_schedule_t* const nsearch_schedule) {
  uint64_t i, total_min_error = 0;
  for (i=0;i<nsearch_schedule->num_pending_searches;++i) {
    nsearch_operation_t* const pending_search = nsearch_schedule->pending_searches + i;
    total_min_error += pending_search->min_local_error;
  }
  return total_min_error;
}
/*
 * Display
 */
void nsearch_schedule_print(
    FILE* const stream,
    nsearch_schedule_t* const nsearch_schedule) {
  const uint64_t key_length = nsearch_schedule->pattern->key_length;
  uint64_t offset = 0;
  while (offset != key_length) {
    uint64_t i;
    for (i=0;i<nsearch_schedule->num_pending_searches;++i) {
      nsearch_operation_t* const pending_search = nsearch_schedule->pending_searches + i;
      if (pending_search->local_key_begin == offset) {
        nsearch_operation_print(stream,pending_search);
        offset = pending_search->local_key_end;
        break;
      }
    }
  }
}
void nsearch_schedule_print_region_segment(
    FILE* const stream,
    const uint64_t length,
    const char begin_c,
    const char middle_c,
    const char end_c) {
  uint64_t j;
  if (length < 5) {
    for (j=0;j<length;++j) fprintf(stream,"*");
  } else {
    fprintf(stream,"%c",begin_c);
    for (j=0;j<length-2;++j) fprintf(stream,"%c",middle_c);
    fprintf(stream,"%c",end_c);
  }
}
void nsearch_schedule_print_region_error(
    FILE* const stream,
    const uint64_t length,
    const uint64_t min_local,
    const uint64_t min_global,
    const uint64_t max_global) {
  uint64_t j;
  if (length < 8) {
    for (j=0;j<length;++j) fprintf(stream,"*");
  } else {
    const uint64_t left = (length-8)/2;
    const uint64_t right = (length-8)-left;
    for (j=0;j<left;++j) fprintf(stream," ");
    fprintf(stream,"{%"PRIu64"}{%"PRIu64",%"PRIu64"}",min_local,min_global,max_global);
    for (j=0;j<right;++j) fprintf(stream," ");
  }
}
void nsearch_schedule_print_region_limits(
    FILE* const stream,
    const uint64_t length,
    const uint64_t min,
    const uint64_t max) {
  uint64_t j;
  if (length < 5) {
    for (j=0;j<length;++j) fprintf(stream,"*");
  } else {
    const uint64_t left = (length-5)/2;
    const uint64_t right = (length-5)-left;
    for (j=0;j<left;++j) fprintf(stream," ");
    fprintf(stream,"{%"PRIu64",%"PRIu64"}",min,max);
    for (j=0;j<right;++j) fprintf(stream," ");
  }
}
typedef struct {
  uint64_t nsearch_schedule_pos;
  uint64_t plength;
} nsearch_schedule_print_data_t;
void nsearch_schedule_print_pretty(
    FILE* const stream,
    nsearch_schedule_t* const nsearch_schedule) {
  // Save allocator state & allocate mem
  mm_allocator_t* const mm_allocator = nsearch_schedule->mm_allocator;
  mm_allocator_push_state(mm_allocator);
  const uint64_t num_pending_searches = nsearch_schedule->num_pending_searches;
  nsearch_schedule_print_data_t* const print_data =
      mm_allocator_calloc(mm_allocator,num_pending_searches,nsearch_schedule_print_data_t,true);
  // Set proper amplification factor
  const uint64_t key_length = nsearch_schedule->pattern->key_length;
  uint64_t amplification = 1, i;
  if (key_length < 100) {
    amplification = 100 / key_length;
    if (amplification == 0) amplification = 1;
  }
  // Compute print info
  nsearch_operation_t* const pending_searches = nsearch_schedule->pending_searches;
  int64_t offset = 0, j=0;
  while (offset != key_length) {
    for (i=0;i<num_pending_searches;++i) {
      nsearch_operation_t* const pending_search = pending_searches + i;
      if (pending_search->local_key_begin == offset) {
        print_data[j].nsearch_schedule_pos = i;
        print_data[j].plength = amplification*(pending_search->local_key_end-pending_search->local_key_begin);
        offset = pending_search->local_key_end;
        ++j;
        break;
      }
    }
  }
  // Print Header
  fprintf(stream,"[GEM]>NSearch[%"PRIu64"]\n",nsearch_schedule->search_id++);
  //  // Print Scheduled Operations
  //  fprintf(stream,"  => Scheduled.Operations\n");
  //  for (i=0;i<nsearch_schedule->num_pending_searches;++i) {
  //    nsearch_operation_t* const pending_search = pending_searches + (nsearch_schedule->num_pending_searches-i-1);
  //    fprintf(stream,"    => Key.Local[%lu,%lu) Key.Global[%lu,%lu)\n",
  //        pending_search->local_key_begin,pending_search->local_key_end,
  //        pending_search->global_key_begin,pending_search->global_key_end);
  //  }
  //  fprintf(stream,"\n");
  // Print local min/max limits
  fprintf(stream,"  => Error          ");
  for (i=0;i<nsearch_schedule->num_pending_searches;++i) {
    nsearch_schedule_print_data_t* const operation_data = print_data + i;
    nsearch_operation_t* const pending_search = pending_searches + operation_data->nsearch_schedule_pos;
    nsearch_schedule_print_region_error(stream,operation_data->plength,
        pending_search->min_local_error,pending_search->min_global_error,
        pending_search->max_global_error);
  }
  fprintf(stream,"\n");
  // Print intervals
  fprintf(stream,"  => Partition      ");
  for (i=0;i<num_pending_searches;++i) {
    nsearch_schedule_print_data_t* const operation_data = print_data + i;
    nsearch_operation_t* const pending_search = pending_searches + operation_data->nsearch_schedule_pos;
    const char direcction_char = (pending_search->search_direction==direction_forward) ? '>' : '<';
    nsearch_schedule_print_region_segment(stream,operation_data->plength,'{',direcction_char,'}');
  }
  fprintf(stream,"\n");
  // Print regions
  region_profile_t* const region_profile = nsearch_schedule->region_profile;
  if (region_profile->num_filtering_regions > 0) {
    const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
    fprintf(stream,"  => Regions        ");
    for (i=0;i<num_filtering_regions;++i) {
      region_search_t* const region = region_profile->filtering_region + i;
      const uint64_t plength = amplification*(region->end-region->begin);
      nsearch_schedule_print_region_segment(stream,plength,'|','-','|');
    }
    fprintf(stream,"\n");
    fprintf(stream,"  => Regions.Error  ");
    for (i=0;i<num_filtering_regions;++i) {
      region_search_t* const region = region_profile->filtering_region + i;
      const uint64_t plength = amplification*(region->end-region->begin);
      nsearch_schedule_print_region_limits(stream,plength,region->min,region->max);
    }
    fprintf(stream,"\n");
  }
  // Succint print
  // nsearch_schedule_print(stream,nsearch_schedule);
  mm_allocator_pop_state(mm_allocator);
}
void nsearch_schedule_print_profile(
    FILE* const stream,
    nsearch_schedule_t* const nsearch_schedule) {
  fprintf(stream,"%"PRIu64"\t%"PRIu64"\n",
      nsearch_schedule->profile.ns_nodes_success,
      nsearch_schedule->profile.ns_nodes);
}
