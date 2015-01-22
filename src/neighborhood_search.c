/*
 * PROJECT: GEMMapper
 * FILE: neighborhood_search.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "neighborhood_search.h"

//typedef struct {
//  fm_index_t* fm_index;
//  uint8_t* key;
//  uint64_t length;
//  uint64_t column_length;
//  uint64_t max_error;
//  interval_set_t* intervals_result;
//  mm_stack_t* mm_stack;
//} neighborhood_search_t;

//GEM_INLINE void neighborhood_search_step(
//    neighborhood_search_t* const neighborhood_search,
//    const uint64_t position,uint16_t* const base_column,const bool expand_node) {
//  const mm_stack_t* mm_stack = neighborhood_search->mm_stack;
//  const uint64_t column_length = neighborhood_search->column_length;
//  // Select case
//  if (expand_node) {
//    uint8_t enc_char;
//    for (enc_char=0;enc_char<DNA_RANGE;++enc_char) {
//      // Allowed char
//      if (/* TODO */) continue;
//      // Allocate column
//      uint16_t* const column = mm_stack_calloc(mm_stack,column_length,uint16_t,false);
//      // Update expand column
//      int64_t delta;
//      if ((delta=neighborhood_search_update_column(neighborhood_search,position,base_column,column,enc_char)) >= 0) {
//        if (position < neighborhood_search->length) {
//          // Eligible path
//          neighborhood_search_step(&neighborhood_search,position+1,column,delta>0);
//        } else {
//          // End of the search (Add it to intervals_result)
//          // TODO
//        }
//      }
//    }
//  } else {
//    // Not expansive, exact search
//    uint16_t* const column = mm_stack_calloc(mm_stack,column_length,uint16_t,false); // Allocate column
//    neighborhood_search_update_column(neighborhood_search,
//        position+1,base_column,column,neighborhood_search->key[]);
//  }
//}
GEM_INLINE void neighborhood_search(
    const fm_index_t* const fm_index,const uint8_t* const key,const uint64_t length,
    const uint64_t max_error,interval_set_t* const intervals_result,mm_stack_t* const mm_stack) {
//  // Create neighborhood search DTO
//  neighborhood_search_t neighborhood_search;
//  neighborhood_search.fm_index = fm_index;
//  neighborhood_search.key = key;
//  neighborhood_search.column_length = length+1;
//  neighborhood_search.max_error = max_error;
//  neighborhood_search.intervals_result = intervals_result;
//  neighborhood_search.mm_stack = mm_stack;
//  // Allocate initial column
//  mm_stack_push_state(mm_stack);
//  uint16_t* const base_column = mm_stack_calloc(mm_stack,length+1,uint16_t,false);
//  uint64_t i;
//  for (i=0;i<=length;++i) base_column[i] = i;
//  // Start search
//  neighborhood_search_step(&neighborhood_search,1,base_column,max_error>0);
//  // Free
//  mm_stack_pop_state(mm_stack);
}




