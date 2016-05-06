/*
 * PROJECT: GEMMapper
 * FILE: nsearch_schedule.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "neighborhood_search/nsearch_schedule.h"
#include "neighborhood_search/nsearch_partition.h"
#include "neighborhood_search/nsearch_hamming.h"
#include "neighborhood_search/nsearch_levenshtein.h"
#include "neighborhood_search/nsearch_levenshtein_state.h"

/*
 * Setup
 */
void nsearch_schedule_init(
    nsearch_schedule_t* const nsearch_schedule,const nsearch_model_t nsearch_model,
    fm_index_t* const fm_index,region_profile_t* const region_profile,
    uint8_t* const key,const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const intervals_result,mm_stack_t* const mm_stack) {
  // Search
  nsearch_schedule->fm_index = fm_index;
  nsearch_schedule->nsearch_model = nsearch_model;
  nsearch_schedule->key = key;
  nsearch_schedule->key_length = key_length;
  nsearch_schedule->region_profile = region_profile;
  nsearch_schedule->max_error = max_error;
  const uint64_t max_text_length = key_length + max_error;
  nsearch_schedule->max_text_length = max_text_length;
  nsearch_schedule->intervals_result = intervals_result;
  nsearch_schedule->mm_stack = mm_stack;
  // Pending Search Operations
  const uint64_t max_pending_ops = max_error+1; // TODO Shrink
  nsearch_operation_t* const pending_searches = mm_stack_calloc(mm_stack,max_pending_ops,nsearch_operation_t,false);
  nsearch_schedule->pending_searches = pending_searches;
  if (nsearch_model==nsearch_model_levenshtein) {
    const uint64_t num_columns = max_text_length + 3;
    const uint64_t column_length = key_length + 1;
    uint64_t i;
    for (i=0;i<max_pending_ops;++i) {
      nsearch_levenshtein_state_init(&pending_searches[i].nsearch_state,num_columns,column_length,mm_stack);
    }
  }
  nsearch_schedule->num_pending_searches = 0;
  nsearch_schedule->search_id = 0;
  nsearch_schedule->search_string = mm_stack_calloc(mm_stack,max_text_length,char,true);
  // Profiler
  nsearch_schedule->ns_nodes = 0;
  nsearch_schedule->ns_nodes_mtable = 0;
  nsearch_schedule->ns_nodes_success = 0;
  nsearch_schedule->ns_nodes_closed = 0;
  nsearch_schedule->ns_nodes_fail_optimize = 0;
  COUNTER_RESET(&nsearch_schedule->ns_nodes_closed_depth);
  TIMER_RESET(&nsearch_schedule->ns_timer);
}
/*
 * Add pending search
 */
bool nsearch_schedule_add_pending_search(
    nsearch_schedule_t* const nsearch_schedule,const search_direction_t search_direction,
    const uint64_t offset,const uint64_t length,
    const uint64_t min_local_error,const uint64_t max_local_error,
    const uint64_t min_global_error,const uint64_t max_global_error) {
  // Check impossible search configuration
  if (min_local_error > length) return false;
  // Queue search
  const uint64_t num_op = nsearch_schedule->num_pending_searches;
  nsearch_schedule->pending_searches[num_op].search_direction = search_direction;
  nsearch_schedule->pending_searches[num_op].begin = offset;
  nsearch_schedule->pending_searches[num_op].end = offset+length;
  nsearch_schedule->pending_searches[num_op].min_local_error = min_local_error;
  nsearch_schedule->pending_searches[num_op].max_local_error = max_local_error;
  nsearch_schedule->pending_searches[num_op].min_global_error = min_global_error;
  nsearch_schedule->pending_searches[num_op].max_global_error = max_global_error;
  nsearch_schedule->pending_searches[num_op].max_text_length = length+max_local_error;
  ++(nsearch_schedule->num_pending_searches);
  // Return
  return true;
}
/*
 * Schedule the search
 */
void nsearch_schedule_search_step(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t chunk_offset,const uint64_t chunk_length,
    const uint64_t chunk_min_error,const uint64_t chunk_max_error) {
  // Adjust error
  bool feasible_search;
  if (chunk_max_error==0 || chunk_length<=chunk_max_error) {
    // Exact search of the chunk
    feasible_search = nsearch_schedule_add_pending_search(
        nsearch_schedule,direction_forward,chunk_offset,chunk_length,
        chunk_min_error,chunk_max_error,chunk_min_error,chunk_max_error);
    if (!feasible_search) return; // Impossible search
    // Perform the search (Solve pending extensions)
    nsearch_schedule_print_pretty(stderr,nsearch_schedule); // DEBUG
    const uint64_t num_pending_searches = nsearch_schedule->num_pending_searches - 1;
    nsearch_operation_t* const nsearch_operation = nsearch_schedule->pending_searches + num_pending_searches;
    // Select nsearch-alignment model
    switch (nsearch_schedule->nsearch_model) {
      case nsearch_model_hamming:
        nsearch_hamming_perform_scheduled_search(nsearch_schedule,
            num_pending_searches,nsearch_operation,nsearch_operation->begin,0,0);
        break;
      case nsearch_model_levenshtein: {
        nsearch_levenshtein_state_t* const nsearch_state = &nsearch_operation->nsearch_state;
        nsearch_state->key_begin = nsearch_operation->begin;
        nsearch_state->key_end = nsearch_operation->end;
        nsearch_state->global_text_length = 0;
        nsearch_state->local_text_length = 0;
        nsearch_levenshtein_state_prepare_supercondensed(nsearch_state,nsearch_operation->max_local_error);
        nsearch_levenshtein_perform_scheduled_search(nsearch_schedule,num_pending_searches,nsearch_operation,0);
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
        nsearch_schedule,direction_forward,epartition.offset_1,epartition.length_1,
        epartition.extend_1_local_min_error,epartition.extend_1_local_max_error,
        epartition.extend_1_global_min_error,epartition.extend_1_global_max_error);
    if (!feasible_search) return; // Impossible search
    nsearch_schedule_search_step(
        nsearch_schedule,epartition.offset_0,epartition.length_0,
        epartition.search_0_min_error,epartition.search_0_max_error);
    /*
     * Second partition (backward)
     */
    nsearch_schedule->num_pending_searches = num_pending_searches; // Restore pending searches point
    feasible_search = nsearch_schedule_add_pending_search(
        nsearch_schedule,direction_backward,epartition.offset_0,epartition.length_0,
        epartition.extend_0_local_min_error,epartition.extend_0_local_max_error,
        epartition.extend_0_global_min_error,epartition.extend_0_global_max_error);
    if (!feasible_search) return; // Impossible search
    nsearch_schedule_search_step(
        nsearch_schedule,epartition.offset_1,epartition.length_1,
        epartition.search_1_min_error,epartition.search_1_max_error);
  }
}
void nsearch_schedule_search(nsearch_schedule_t* const nsearch_schedule) {
  TIMER_START(&nsearch_schedule->ns_timer);
  nsearch_schedule_search_step(nsearch_schedule,0,nsearch_schedule->key_length,0,nsearch_schedule->max_error);
  TIMER_STOP(&nsearch_schedule->ns_timer);
}
/*
 * Schedule the search (preconditioned by region profile)
 */
void nsearch_schedule_search_preconditioned_step(
    nsearch_schedule_t* const nsearch_schedule,
    const uint64_t region_offset,const uint64_t num_regions,
    const uint64_t chunk_min_error,const uint64_t chunk_max_error) {
  // PRE: (num_regions > 1)
  region_search_t* const regions = nsearch_schedule->region_profile->filtering_region;
  const uint64_t global_region_length = regions[region_offset+num_regions-1].end - regions[region_offset].begin;
  if (chunk_max_error==0 || global_region_length<=chunk_max_error) {
    nsearch_schedule_search_step(
        nsearch_schedule,regions[region_offset].begin,
        global_region_length,chunk_min_error,chunk_max_error);
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
        nsearch_schedule,direction_forward,epartition.offset_1,epartition.length_1,
        epartition.extend_1_local_min_error,epartition.extend_1_local_max_error,
        epartition.extend_1_global_min_error,epartition.extend_1_global_max_error);
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
    /*
     * Second Partition
     */
    // Extend second partition (backward)
    nsearch_schedule->num_pending_searches = num_pending_searches; // Restore pending searches point
    feasible_search = nsearch_schedule_add_pending_search(
        nsearch_schedule,direction_backward,epartition.offset_0,epartition.length_0,
        epartition.extend_0_local_min_error,epartition.extend_0_local_max_error,
        epartition.extend_0_global_min_error,epartition.extend_0_global_max_error);
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
  }
}
void nsearch_schedule_search_preconditioned(nsearch_schedule_t* const nsearch_schedule) {
  TIMER_START(&nsearch_schedule->ns_timer);
  const uint64_t num_filtering_regions = nsearch_schedule->region_profile->num_filtering_regions;
  if (num_filtering_regions <= 1) {
    nsearch_schedule_search_step(nsearch_schedule,
        0,nsearch_schedule->key_length,0,nsearch_schedule->max_error);
  } else {
    nsearch_schedule_search_preconditioned_step(nsearch_schedule,
        0,num_filtering_regions,0,nsearch_schedule->max_error);
  }
  TIMER_STOP(&nsearch_schedule->ns_timer);
}
/*
 * Display
 */
void nsearch_schedule_print(FILE* const stream,nsearch_schedule_t* const nsearch_schedule) {
  uint64_t offset = 0;
  while (offset != nsearch_schedule->key_length) {
    uint64_t i;
    for (i=0;i<nsearch_schedule->num_pending_searches;++i) {
      nsearch_operation_t* const pending_search = nsearch_schedule->pending_searches + i;
      if (pending_search->begin == offset) {
        fprintf(stream,">>[%lu,%lu) {%lu,%lu}/(%lu,%lu)\n",
            pending_search->begin,pending_search->end,
            pending_search->min_local_error,pending_search->max_local_error,
            pending_search->min_global_error,pending_search->max_global_error);
        offset = pending_search->end;
        break;
      }
    }
  }
}
typedef struct {
  uint64_t nsearch_schedule_pos;
  uint64_t plength;
} nsearch_schedule_print_data_t;
void nsearch_schedule_print_region_segment(
    FILE* const stream,const uint64_t length,
    const char begin_c,const char middle_c,const char end_c) {
  uint64_t j;
  if (length < 5) {
    for (j=0;j<length;++j) fprintf(stream,"*");
  } else {
    fprintf(stream,"%c",begin_c);
    for (j=0;j<length-2;++j) fprintf(stream,"%c",middle_c);
    fprintf(stream,"%c",end_c);
  }
}
void nsearch_schedule_print_region_limits(
    FILE* const stream,const uint64_t length,
    const uint64_t min,const uint64_t max) {
  uint64_t j;
  if (length < 5) {
    for (j=0;j<length;++j) fprintf(stream,"*");
  } else {
    const uint64_t left = (length-5)/2;
    const uint64_t right = (length-5)-left;
    for (j=0;j<left;++j) fprintf(stream," ");
    fprintf(stream,"{%lu,%lu}",min,max);
    for (j=0;j<right;++j) fprintf(stream," ");
  }
}
void nsearch_schedule_print_pretty(FILE* const stream,nsearch_schedule_t* const nsearch_schedule) {
  // Save stack state & allocate mem
  mm_stack_t* const mm_stack = nsearch_schedule->mm_stack;
  mm_stack_push_state(mm_stack);
  const uint64_t num_pending_searches = nsearch_schedule->num_pending_searches;
  nsearch_schedule_print_data_t* const print_data =
      mm_stack_calloc(mm_stack,num_pending_searches,nsearch_schedule_print_data_t,true);
  // Set proper amplification factor
  uint64_t amplification = 1;
  if (nsearch_schedule->key_length < 100) {
    amplification = 100 / nsearch_schedule->key_length;
    if (amplification == 0) amplification = 1;
  }
  // Compute print info
  nsearch_operation_t* const pending_searches = nsearch_schedule->pending_searches;
  int64_t offset = 0, j=0, i;
  while (offset != nsearch_schedule->key_length) {
    for (i=0;i<num_pending_searches;++i) {
      nsearch_operation_t* const pending_search = pending_searches + i;
      if (pending_search->begin == offset) {
        print_data[j].nsearch_schedule_pos = i;
        print_data[j].plength = amplification*(pending_search->end-pending_search->begin);
        offset = pending_search->end;
        ++j;
        break;
      }
    }
  }
  // Print Header & local min/max limits
  fprintf(stream,"[GEM]>NSearch[%lu]\n",nsearch_schedule->search_id++);
  for (i=0;i<nsearch_schedule->num_pending_searches;++i) {
    nsearch_schedule_print_data_t* const operation_data = print_data + i;
    nsearch_operation_t* const pending_search = pending_searches + operation_data->nsearch_schedule_pos;
    nsearch_schedule_print_region_limits(stream,operation_data->plength,
        pending_search->min_local_error,pending_search->max_local_error);
  }
  fprintf(stream,"\n");
  // Print intervals
  for (i=0;i<num_pending_searches;++i) {
    nsearch_schedule_print_data_t* const operation_data = print_data + i;
    nsearch_operation_t* const pending_search = pending_searches + operation_data->nsearch_schedule_pos;
    const char direcction_char = (pending_search->search_direction==direction_forward) ? '>' : '<';
    nsearch_schedule_print_region_segment(stream,operation_data->plength,'{',direcction_char,'}');
  }
  fprintf(stream,"\n");
  // Print regions
  region_profile_t* const region_profile = nsearch_schedule->region_profile;
  if (region_profile!=NULL) {
    const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
    for (i=0;i<num_filtering_regions;++i) {
      region_search_t* const region = region_profile->filtering_region + i;
      const uint64_t plength = amplification*(region->end-region->begin);
      nsearch_schedule_print_region_segment(stream,plength,'|','-','|');
    }
    fprintf(stream,"\n");
    for (i=0;i<num_filtering_regions;++i) {
      region_search_t* const region = region_profile->filtering_region + i;
      const uint64_t plength = amplification*(region->end-region->begin);
      nsearch_schedule_print_region_limits(stream,plength,region->min,region->max);
    }
    fprintf(stream,"\n");
  }
  // Succint print
  // nsearch_schedule_print(stream,nsearch_schedule);
  mm_stack_pop_state(mm_stack);
}
void nsearch_schedule_print_profile(FILE* const stream,nsearch_schedule_t* const nsearch_schedule) {
  fprintf(stderr,"%lu\t%lu\t%2.3f\n",nsearch_schedule->ns_nodes_success,
      nsearch_schedule->ns_nodes,TIMER_GET_TOTAL_S(&nsearch_schedule->ns_timer));
}
void nsearch_schedule_print_search_string(FILE* const stream,nsearch_schedule_t* const nsearch_schedule) {
  fprintf(stream,"%s",nsearch_schedule->search_string);
}
