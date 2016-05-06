/*
 * PROJECT: GEMMapper
 * FILE: neighborhood_search.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "neighborhood_search/nsearch.h"
#include "neighborhood_search/dp_matrix.h"
#include "data_structures/dna_text.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * DEBUG
 */
#define NS_PRINT_STRING false
#define NS_PRINT_MATRIX false

char *ns_string = NULL;

/*
 * Constants
 */
#define NSC_PROPER_CELL_SET_MASK      UINT64_ONE_LAST_MASK
#define NSC_PROPER_CELL_EXTRACT_MASK  (~UINT64_ONE_LAST_MASK)

/*
 * Profiler/Stats
 */
uint64_t _ns_nodes = 0;
uint64_t _ns_nodes_mtable = 0;
uint64_t _ns_nodes_success = 0;
uint64_t _ns_nodes_closed = 0;
uint64_t _ns_nodes_fail_optimize = 0;
gem_counter_t _ns_nodes_closed_depth;

/*
 * DTO
 */
typedef struct {
  // FM-Index
  fm_index_t* fm_index;
  // Key
  uint8_t* key;
  uint64_t key_length;
  uint64_t column_length;
  // Search Parameters
  uint64_t max_error;
  uint64_t max_text_length;
  // Results
  interval_set_t* intervals_result;
  // MM
  mm_stack_t* mm_stack;
} neighborhood_search_t;

/*
 * Debug
 */
void neighborhood_dp_column_print(
    neighborhood_search_t* const neighborhood_search,
    uint64_t* const dp_column,const uint64_t column_position,
    const uint64_t lo,const uint64_t hi) {
  const uint64_t key_length = neighborhood_search->key_length;
  uint64_t i, min = UINT64_MAX;
  for (i=0;i<=key_length;++i) min = MIN(min,dp_column[i]);
  gem_slog(">> [%"PRIu64"](#%"PRIu64"){min=%"PRIu64",last=%"PRIu64"}\n",column_position,hi-lo,min,dp_column[key_length]);
}
void neighborhood_dp_matrix_print(
    neighborhood_search_t* const neighborhood_search,
    dp_column_t* const dp_columns,const uint64_t num_rows,const uint64_t num_columns) {
  int64_t h, v;
  for (v=0;v<num_rows;++v) {
    for (h=0;h<num_columns;++h) {
      const uint64_t cell = dp_columns[h].cells[v];
      gem_slog("%c%02llu ",
          cell & NSC_PROPER_CELL_SET_MASK ? '*' : ' ',
          (uint64_t) cell & NSC_PROPER_CELL_EXTRACT_MASK);
    }
    gem_slog("\n");
  }
}
void neighborhood_dp_matrix_traceback(
    neighborhood_search_t* const neighborhood_search,
    dp_column_t* const dp_columns,const uint64_t column_position) {
  const uint8_t* const key = neighborhood_search->key;
  const uint64_t key_length = neighborhood_search->key_length;
  int64_t h = column_position;
  int64_t v = key_length;
  while (h > 0 && v > 0) {
    const bool match = dna_encode(ns_string[h-1])==key[key_length-v];
    const uint64_t current = (dp_columns[h].cells[v] & NSC_PROPER_CELL_EXTRACT_MASK);
    const uint64_t del = (dp_columns[h].cells[v-1] & NSC_PROPER_CELL_EXTRACT_MASK) + 1;
    const uint64_t sub = (dp_columns[h-1].cells[v-1] & NSC_PROPER_CELL_EXTRACT_MASK) + (match ? 0 : 1);
    // const uint64_t ins = (dp_columns[h-1][v] & NSC_PROPER_CELL_EXTRACT_MASK) + 1;
    if (current == del) {
      gem_slog("D"); --v;
    } else if (current == sub) {
      if (!match) {
        gem_slog("%c",ns_string[h]);
      } else {
        gem_slog(".");
      }
      --v; --h;
    } else {
      gem_slog("I"); --h;
    }
  }
  while (v-- > 0) gem_slog("D");
  while (h-- > 0) gem_slog("I");
  gem_slog("\n");
}
void neighborhood_search_debug_match(
    neighborhood_search_t* const neighborhood_search,
    dp_column_t* const dp_columns,const uint64_t column_position,
    const uint64_t min_val,const uint64_t num_matches_found) {
  const uint64_t key_length = neighborhood_search->key_length;
  dp_column_t* const next_column = dp_columns + column_position;
  uint64_t i;
  for (i=0;i<=key_length;++i) gem_slog("%"PRId64" ",next_column->cells[i]>1000 ? -1 : (int64_t)next_column->cells[i]);
  gem_slog("\n");
  neighborhood_dp_matrix_traceback(neighborhood_search,dp_columns,column_position);
  gem_cond_debug_block(NS_PRINT_MATRIX) { neighborhood_dp_matrix_print(neighborhood_search,dp_columns,10,10); }
  gem_slog("\n[%02"PRIu64"](%03"PRIu64")> %.*s\n",min_val,num_matches_found,(int)column_position,ns_string);
}
/*
 * Condensed Neighborhood
 */
uint64_t neighborhood_compute_dp_column(
    neighborhood_search_t* const neighborhood_search,
    uint64_t* const base_column,uint64_t* const next_column,
    const uint64_t column_position,const uint8_t text_enc) {
  const uint8_t* const key = neighborhood_search->key;
  const uint64_t key_length = neighborhood_search->key_length;
  // Fill column
  uint64_t i, min = UINT64_MAX;
  next_column[0] = column_position + 1; // Loop-peeling
  for (i=1;i<=key_length;++i) {
    const uint64_t sub = base_column[i-1] + (text_enc==key[key_length-i] ? 0 : 1);
    const uint64_t ins = base_column[i] + 1;
    const uint64_t del = next_column[i-1] + 1;
    next_column[i] = MIN(sub,MIN(ins,del));
    min = MIN(min,next_column[i]);
  }
  return min;
}
void neighborhood_condensed_search(
    neighborhood_search_t* const neighborhood_search,
    uint64_t** const dp_columns,const uint64_t column_position,
    const uint64_t lo,const uint64_t hi) {
  // Parameters
  bwt_t* const bwt = neighborhood_search->fm_index->bwt;
  const uint64_t key_length = neighborhood_search->key_length;
  const uint64_t max_error = neighborhood_search->max_error;
  // DP-Columns
  uint64_t* const base_column = dp_columns[column_position];
  uint64_t* const next_column = dp_columns[column_position+1];
  // Advance for all characters
  uint64_t min_val;
  uint8_t enc;
  for (enc=0;enc<DNA_RANGE;++enc) {
    // Compute DP-next
    min_val = neighborhood_compute_dp_column(neighborhood_search,base_column,next_column,column_position,enc);
    if (min_val <= max_error) {
      // Compute erank
      const uint64_t next_lo = bwt_erank(bwt,enc,lo);
      const uint64_t next_hi = bwt_erank(bwt,enc,hi);
      if (next_lo >= next_hi) continue;
      // Search (Condensed)
      if (next_column[key_length]<=max_error) {
        interval_set_add(neighborhood_search->intervals_result,
            next_lo,next_hi,next_column[key_length],column_position+1);
      } else if (min_val <= max_error && column_position+1<neighborhood_search->max_text_length) {
        neighborhood_condensed_search(neighborhood_search,dp_columns,column_position+1,next_lo,next_hi);
      }
    }
  }
}
/*
 * SuperCondensed Neighborhood
 */
uint64_t neighborhood_supercondensed_compute_dp_column(
    neighborhood_search_t* const neighborhood_search,
    uint64_t* const base_column,uint64_t* const next_column,
    const uint64_t column_position,const uint8_t text_enc) {
  const uint8_t* const key = neighborhood_search->key;
  const uint64_t key_length = neighborhood_search->key_length;
  // Fill column
  uint64_t i, column_min = UINT64_MAX;
  next_column[0] = (column_position + 1) | NSC_PROPER_CELL_SET_MASK; // Loop-peeling
  for (i=1;i<=key_length;++i) {
    const uint64_t del = next_column[i-1] + 1;
    const uint64_t sub = base_column[i-1] + (text_enc==key[key_length-i] ? 0 : 1);
    const uint64_t ins = base_column[i] + 1;
    const uint64_t m_del = del & NSC_PROPER_CELL_EXTRACT_MASK;
    const uint64_t m_sub = sub & NSC_PROPER_CELL_EXTRACT_MASK;
    const uint64_t m_ins = ins & NSC_PROPER_CELL_EXTRACT_MASK;
    const uint64_t min2 = MIN(m_del,m_sub);
    const uint64_t min3 = MIN(min2,m_ins);
    if (min3 == m_del) {
      next_column[i] = del;
    } else if (min3 == m_sub) {
      next_column[i] = sub;
    } else { // min3 == m_ins
      next_column[i] = ins;
    }
    column_min = MIN(column_min,next_column[i]);
  }
  return column_min;
}
void neighborhood_supercondensed_search(
    neighborhood_search_t* const neighborhood_search,
    uint64_t** const dp_columns,const uint64_t column_position,
    const uint64_t lo,const uint64_t hi) {
  // Parameters
  bwt_t* const bwt = neighborhood_search->fm_index->bwt;
  const uint64_t key_length = neighborhood_search->key_length;
  const uint64_t max_error = neighborhood_search->max_error;
  // DP-Columns
  uint64_t* const base_column = dp_columns[column_position];
  uint64_t* const next_column = dp_columns[column_position+1];
  // Advance for all characters
  uint64_t min_val;
  uint8_t enc;
  for (enc=0;enc<DNA_RANGE;++enc) {
    gem_cond_debug_block(NS_PRINT_STRING) { ns_string[column_position] = dna_decode(enc); }
    // Compute DP-next
    min_val = neighborhood_supercondensed_compute_dp_column(neighborhood_search,base_column,next_column,column_position,enc);
    if (min_val <= max_error) {
      // Compute erank
      const uint64_t next_lo = bwt_erank(bwt,enc,lo);
      const uint64_t next_hi = bwt_erank(bwt,enc,hi);
      if (next_lo >= next_hi) continue;
      // Search (Supercondensed)
      if (next_column[key_length]<=max_error) {
        interval_set_add(neighborhood_search->intervals_result,
            next_lo,next_hi,next_column[key_length],column_position+1);
      } else if (min_val <= max_error && column_position+1<neighborhood_search->max_text_length) {
        neighborhood_supercondensed_search(neighborhood_search,dp_columns,column_position+1,next_lo,next_hi);
      }
    }
  }
}
/*
 * Best-Match
 */
uint64_t neighborhood_best_matches_search_bwt(
    neighborhood_search_t* const neighborhood_search,
    dp_column_t* const dp_columns,const uint64_t current_column,
    const uint64_t lo,const uint64_t hi);
uint64_t neighborhood_best_matches_search_mtable(
    neighborhood_search_t* const neighborhood_search,
    dp_column_t* const dp_columns,const uint64_t current_column,
    rank_mquery_t* const query);
dp_column_t* neighborhood_best_matches_search_init(neighborhood_search_t* const neighborhood_search) {
  // Allocate columns
  mm_stack_t* const mm_stack = neighborhood_search->mm_stack;
  const uint64_t max_text_length = neighborhood_search->max_text_length;
  dp_column_t* const dp_columns = mm_stack_calloc(mm_stack,max_text_length+2,dp_column_t,false);
  uint64_t i;
  for (i=0;i<=max_text_length;++i) {
    dp_columns[i].cells = mm_stack_calloc(mm_stack,neighborhood_search->column_length,uint64_t,false);
  }
  // Init columns
  const uint64_t key_length = neighborhood_search->key_length;
  const uint64_t max_error = neighborhood_search->max_error;
  dp_columns[0].cells[0] = 1;
  for (i=1;i<=key_length;++i) dp_columns[0].cells[i] = NS_ENCODE_DISTANCE(i);
  const uint64_t column_start_band = max_error+1;
  for (i=0;i<max_text_length;++i) {
    dp_column_t* base_column = dp_columns + i;
    dp_column_t* next_column = base_column + 1;
    uint64_t _band_high_offset = max_error + i + 1;
    next_column->band_high_offset = MIN(_band_high_offset,key_length);
    if (i < column_start_band) {
      next_column->band_low_offset = 1;
      next_column->cells[0] = NS_ENCODE_DISTANCE(i + 1); // Loop-peeling
    } else {
      next_column->band_low_offset = i - column_start_band + 1;
      next_column->cells[next_column->band_low_offset-1] = NS_INF; // Loop-peeling
    }
    base_column->cells[next_column->band_high_offset] = NS_INF;
    next_column->cells[key_length] = NS_INF;
  }
  // Return
  return dp_columns;
}
void ns_best_matches_compute_dp_column_banded(
    neighborhood_search_t* const neighborhood_search,
    dp_column_t* const dp_columns,const uint64_t current_column,
    const uint64_t max_error,const uint8_t text_enc,
    uint64_t* const min_val,uint64_t* const align_distance) {
  // Parameters
  const uint8_t* const key = neighborhood_search->key;
  const uint64_t key_length = neighborhood_search->key_length;
  dp_column_t* const base_column = dp_columns + (current_column-1);
  dp_column_t* const next_column = dp_columns + current_column;
  // Fill columns
  const uint64_t band_low_offset = next_column->band_low_offset;
  const uint64_t band_high_offset = next_column->band_high_offset;
  uint64_t i, column_min = NS_INF;
  for (i=band_low_offset;i<=band_high_offset;++i) {
    // Compute cells
    const uint64_t del = next_column->cells[i-1] + 2;
    const uint64_t sub = base_column->cells[i-1] + (text_enc==key[key_length-i] ? 0 : 2);
    const uint64_t ins = base_column->cells[i] + 2;
    // Compute min
    const uint64_t min2 = MIN(del,sub);
    const uint64_t min3 = MIN(min2,ins);
    next_column->cells[i] = min3;
    if (NS_HAS_PRIORITY(min3,0)) column_min = MIN(column_min,min3);
  }
  *min_val = NS_HAS_PRIORITY(column_min,0) ? NS_DECODE_DISTANCE(column_min) : NS_INF;
  *align_distance = NS_DECODE_DISTANCE(next_column->cells[key_length]);
}
uint64_t neighborhood_best_matches_search_continue(
    neighborhood_search_t* const neighborhood_search,
    dp_column_t* const dp_columns,const uint64_t current_column,
    const uint64_t min_val,const uint64_t align_distance,
    const uint64_t max_error,rank_mquery_t* const next_query,
    const uint64_t lo,const uint64_t hi) {
  // Search
  uint64_t matches_found = 0;
  if (min_val < align_distance) {
    if (current_column < neighborhood_search->max_text_length) {
      matches_found = (next_query==NULL || rank_mquery_is_exhausted(next_query)) ?
          neighborhood_best_matches_search_bwt(neighborhood_search,dp_columns,current_column+1,lo,hi) :
          neighborhood_best_matches_search_mtable(neighborhood_search,dp_columns,current_column+1,next_query);
      PROF_BLOCK() { if (matches_found == 0) { ++_ns_nodes_fail_optimize; } }
    }
    if (matches_found==0 && align_distance<=max_error) {
      matches_found = hi - lo;
      gem_cond_debug_block(NS_PRINT_STRING) {
        neighborhood_search_debug_match(neighborhood_search,dp_columns,current_column,min_val,matches_found);
      }
      PROF_BLOCK() {++_ns_nodes_success;}
      interval_set_add(neighborhood_search->intervals_result,lo,hi,align_distance,current_column);
    }
  } else { // min_val == next_column[key_length]
    matches_found = hi - lo;
    gem_cond_debug_block(NS_PRINT_STRING) {
      neighborhood_search_debug_match(neighborhood_search,dp_columns,current_column,min_val,matches_found);
    }
    PROF_BLOCK() {++_ns_nodes_success;}
    interval_set_add(neighborhood_search->intervals_result,lo,hi,align_distance,current_column);
  }
  return matches_found;
}
uint64_t neighborhood_best_matches_search_bwt(
    neighborhood_search_t* const neighborhood_search,
    dp_column_t* const dp_columns,const uint64_t current_column,
    const uint64_t lo,const uint64_t hi) {
  // Profiler
  PROF_BLOCK() { ++_ns_nodes; }
  // Parameters
  bwt_t* const bwt = neighborhood_search->fm_index->bwt;
  const uint64_t max_error = neighborhood_search->max_error;
  // Advance for all characters
  uint64_t total_matches_found = 0;
  uint8_t enc;
  bwt_block_locator_t bwt_block_locator_lo;
  bwt_block_locator_t bwt_block_locator_hi;
  bwt_block_elms_t bwt_block_elms_lo, bwt_block_elms_hi;
  bwt_precompute(bwt,lo,&bwt_block_locator_lo,&bwt_block_elms_lo);
  bwt_precompute(bwt,hi,&bwt_block_locator_hi,&bwt_block_elms_hi);
  for (enc=0;enc<DNA_RANGE;++enc) {
    gem_cond_debug_block(NS_PRINT_STRING) { ns_string[current_column-1] = dna_decode(enc); }
    // Compute DP-next
    uint64_t min_val, align_distance;
    ns_best_matches_compute_dp_column_banded(
        neighborhood_search,dp_columns,current_column,max_error,enc,&min_val,&align_distance);
    if (min_val <= max_error) {
      // Compute erank
      const uint64_t next_lo = bwt_precomputed_erank(bwt,enc,&bwt_block_locator_lo,&bwt_block_elms_lo);
      const uint64_t next_hi = bwt_precomputed_erank(bwt,enc,&bwt_block_locator_hi,&bwt_block_elms_hi);
      if (next_lo >= next_hi) continue;
      // Keep searching
      total_matches_found += neighborhood_best_matches_search_continue(
          neighborhood_search,dp_columns,current_column,
          min_val,align_distance,max_error,NULL,next_lo,next_hi);
    }
  }
  PROF_BLOCK() {
    if (total_matches_found==0) { ++_ns_nodes_closed; COUNTER_ADD(&_ns_nodes_closed_depth,current_column); }
  }
  return total_matches_found;
}
uint64_t neighborhood_best_matches_search_mtable(
    neighborhood_search_t* const neighborhood_search,dp_column_t* const dp_columns,
    const uint64_t current_column,rank_mquery_t* const query) {
  // Profiler
  PROF_BLOCK() { ++_ns_nodes; ++_ns_nodes_mtable; }
  // Parameters
  const rank_mtable_t* const rank_mtable = neighborhood_search->fm_index->rank_table;
  const uint64_t max_error = neighborhood_search->max_error;
  // Advance for all characters
  uint64_t total_matches_found = 0;
  uint8_t enc;
  for (enc=0;enc<DNA_RANGE;++enc) {
    gem_cond_debug_block(NS_PRINT_STRING) { ns_string[current_column-1] = dna_decode(enc); }
    // Compute DP-next
    uint64_t min_val, align_distance;
    ns_best_matches_compute_dp_column_banded(
        neighborhood_search,dp_columns,current_column,max_error,enc,&min_val,&align_distance);
    if (min_val <= max_error) {
      // Compute erank using mtable
      rank_mquery_t next_query;
      uint64_t next_lo, next_hi;
      next_query = *query;
      rank_mquery_add_char(rank_mtable,&next_query,enc);
      rank_mtable_fetch(rank_mtable,&next_query,&next_lo,&next_hi);
      if (next_lo >= next_hi) continue;
      // Keep searching
      total_matches_found += neighborhood_best_matches_search_continue(
          neighborhood_search,dp_columns,current_column,min_val,
          align_distance,max_error,&next_query,next_lo,next_hi);
    }
  }
  PROF_BLOCK() {
    if (total_matches_found==0) { ++_ns_nodes_closed; COUNTER_ADD(&_ns_nodes_closed_depth,current_column); }
  }
  return total_matches_found;
}
/*
 * Main
 */
void neighborhood_search(
    fm_index_t* const fm_index,uint8_t* const key,
    const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const intervals_result,mm_stack_t* const mm_stack) {
  PROFILE_START(GP_NSEARCH,PROFILE_LEVEL);
  // DEBUG/PROF
  gem_cond_debug_block(NS_PRINT_STRING) { if (ns_string==NULL) ns_string = malloc(1000); }
#ifdef GEM_PROFILE
  const uint64_t _ns_nodes_start = _ns_nodes;
  const uint64_t _ns_nodes_mtable_start = _ns_nodes_mtable;
  const uint64_t _ns_nodes_success_start = _ns_nodes_success;
  const uint64_t _ns_nodes_closed_start = _ns_nodes_closed;
  const uint64_t _ns_nodes_fail_optimize_start = _ns_nodes_fail_optimize;
#endif
  // Create neighborhood search DTO
  neighborhood_search_t neighborhood_search;
  neighborhood_search.fm_index = fm_index;
  neighborhood_search.key = key;
  neighborhood_search.key_length = key_length;
  neighborhood_search.column_length = key_length+1;
  neighborhood_search.max_error = max_error;
  neighborhood_search.max_text_length = key_length + max_error;
  neighborhood_search.intervals_result = intervals_result;
  neighborhood_search.mm_stack = mm_stack;

//  // Allocate columns
//  mm_stack_push_state(mm_stack);
//  uint64_t** const dp_columns = mm_stack_calloc(mm_stack,neighborhood_search.max_text_length+2,uint64_t*,false);
//  uint64_t i;
//  for (i=0;i<=neighborhood_search.max_text_length;++i) {
//    dp_columns[i] = mm_stack_calloc(mm_stack,neighborhood_search.column_length,uint64_t,false);
//  }

//  // Condensed
//  for (i=0;i<=key_length;++i) dp_columns[0][i] = i;
//  neighborhood_condensed_search(&neighborhood_search,dp_columns,0,0,fm_index_get_length(fm_index));

//  // Supercondensed
//  for (i=0;i<=key_length;++i) dp_columns[0][i] = i;
//  neighborhood_supercondensed_search(&neighborhood_search,dp_columns,0,0,fm_index_get_length(fm_index));

  // Best-match
  PROF_START(GP_NS_BEST_MATCH);
  mm_stack_push_state(mm_stack);
  rank_mquery_t query;
  dp_column_t* const dp_columns = neighborhood_best_matches_search_init(&neighborhood_search);
  rank_mquery_new(&query);
  neighborhood_best_matches_search_mtable(&neighborhood_search,dp_columns,1,&query);
  PROF_STOP(GP_NS_BEST_MATCH);
  // PROF
  PROF_ADD_COUNTER(GP_NS_NODES_EXPLORED,(_ns_nodes-_ns_nodes_start));
  PROF_ADD_COUNTER(GP_NS_NODES_EXPLORED_MTABLE,(_ns_nodes_mtable-_ns_nodes_mtable_start));
  PROF_ADD_COUNTER(GP_NS_NODE_SUCCESS,(_ns_nodes_success-_ns_nodes_success_start));
  PROF_ADD_COUNTER(GP_NS_FAILED_OPT,(_ns_nodes_fail_optimize-_ns_nodes_fail_optimize_start));
  PROF_ADD_COUNTER(GP_NS_NODE_CLOSED,(_ns_nodes_closed-_ns_nodes_closed_start));
  // Free
  mm_stack_pop_state(mm_stack);
  PROFILE_STOP(GP_NSEARCH,PROFILE_LEVEL);
}



