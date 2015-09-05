/*
 * PROJECT: GEMMapper
 * FILE: neighborhood_search.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "neighborhood_search_hamming.h"

/*
 * Profile
 */
uint64_t _ns_nodes = 0;
uint64_t _ns_nodes_success = 0;
uint64_t _ns_nodes_closed = 0;

/*
 * Pending Search Operation
 */
typedef enum { direction_forward, direction_backward } direction_t;
typedef struct {
  // Chunk
  uint64_t offset;
  uint64_t length;
  // Direction
  direction_t direction;
  // Error
  uint64_t max_error;
} pending_search_operation_t;

/*
 * Bidirectional Search
 */
typedef struct {
  // FM-Index
  fm_index_t* fm_index;
  // Key
  uint8_t* key;
  uint64_t key_length;
  // Search Parameters
  uint64_t max_error;
  uint8_t* search_error_mask_min; // Search min-error mask
  // Pending Search Operations
  pending_search_operation_t* pending_searches;
  uint64_t num_pending_searches;
  // Results
  interval_set_t* intervals_result;
  // MM
  mm_stack_t* mm_stack;
  // Debug
  char* search_string;
} bidirectional_search_t;

/*
 * Setup
 */
GEM_INLINE void bidirectional_search_init(
    bidirectional_search_t* const bidirectional_search,fm_index_t* const fm_index,
    uint8_t* const key,const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const intervals_result,mm_stack_t* const mm_stack) {
  bidirectional_search->fm_index = fm_index;
  bidirectional_search->key = key;
  bidirectional_search->key_length = key_length;
  bidirectional_search->max_error = max_error;
  bidirectional_search->intervals_result = intervals_result;
  bidirectional_search->mm_stack = mm_stack;
  // Search min-error mask
  bidirectional_search->search_error_mask_min = mm_stack_calloc(mm_stack,key_length,uint8_t,true);
  // Pending Search Operations
  const uint64_t max_pending_ops = max_error; // FIXME
  bidirectional_search->pending_searches = mm_stack_calloc(mm_stack,max_pending_ops,pending_search_operation_t,false);
  bidirectional_search->num_pending_searches = 0;
  // Debug
  bidirectional_search->search_string = mm_stack_calloc(mm_stack,key_length,char,true);
}
/*
 * Utils
 */
GEM_INLINE void bidirectional_search_push_operation(
    bidirectional_search_t* const bidirectional_search,const direction_t direction,
    const uint64_t offset,const uint64_t length,const uint64_t max_error) {
  const uint64_t num_op = bidirectional_search->num_pending_searches;
  bidirectional_search->pending_searches[num_op].direction = direction;
  bidirectional_search->pending_searches[num_op].offset = offset;
  bidirectional_search->pending_searches[num_op].length = length;
  bidirectional_search->pending_searches[num_op].max_error = max_error;
  ++(bidirectional_search->num_pending_searches);
}
GEM_INLINE void bidirectional_search_error_mask_update_min(
    bidirectional_search_t* const bidirectional_search,const uint64_t offset,
    const uint64_t length,const uint8_t search_error) {
  uint64_t i;
  const uint64_t limit = offset+length;
  for (i=offset;i<limit;++i) {
    bidirectional_search->search_error_mask_min[i] = search_error;
  }
//  // DEBUG
//  for (i=0;i<bidirectional_search->key_length;++i) {
//    fprintf(stdout,"%d",bidirectional_search->search_error_mask_min[i]);
//  }
//  fprintf(stdout,"\n");
}
/*
 *
 */
GEM_INLINE void neighborhood_search_hamming_extend(
    bidirectional_search_t* const bidirectional_search,const uint64_t pending_searches,
    const direction_t direction,const uint64_t position,const uint64_t limit,
    const uint64_t base_chunk_error,const uint64_t current_error,const uint64_t max_error) {
  // Expand node
  uint8_t enc;
  for (enc=0;enc<DNA_RANGE;++enc) {
    // Compute distance function (DP or whatever)
    const uint64_t next_error = current_error + ((enc!=bidirectional_search->key[position]) ? 1 : 0);
    if (next_error <= max_error) {
      // Search character
      bidirectional_search->search_string[position] = dna_decode(enc);
      ++_ns_nodes; // PROFILING
      // Keep searching
      const uint64_t next_position = position + 1;
      if (next_position < limit) {
        neighborhood_search_hamming_extend(bidirectional_search,pending_searches,
            direction,next_position,limit,base_chunk_error,next_error,max_error);
      } else { // (next_position == limit) End of the chunk-search
        // Check min-error mask
        const uint64_t error_in_chunk = next_error - base_chunk_error;
        if (error_in_chunk < bidirectional_search->search_error_mask_min[position]) continue;
        // Check pending searches
        if (pending_searches==0) {
          // End of the search (Print search-string)
          ++_ns_nodes_success;
          fprintf(stdout,"%s\n",bidirectional_search->search_string);
        } else {
          // Search next chunk
          pending_search_operation_t* const pending_search = bidirectional_search->pending_searches + (pending_searches-1);
          neighborhood_search_hamming_extend(bidirectional_search,pending_searches-1,pending_search->direction,
              pending_search->offset,pending_search->offset+pending_search->length,
              next_error,next_error,pending_search->max_error);
        }
      }
    }
  }
}
GEM_INLINE void neighborhood_search_hamming_bidirectional_schedule__search(
    bidirectional_search_t* const bidirectional_search,
    const uint64_t chunk_offset,const uint64_t chunk_length,const uint64_t chunk_max_error) {
  if (chunk_max_error==0) {
    // Exact search of the chunk
    int64_t i;
    for (i=0;i<chunk_length;++i) {
      bidirectional_search->search_string[chunk_offset+i] = dna_decode(bidirectional_search->key[chunk_offset+i]);
      ++_ns_nodes;
    }
    // Perform the search (Solve pending extensions)
    const uint64_t num_pending_searches = bidirectional_search->num_pending_searches;
    pending_search_operation_t* const pending_search = bidirectional_search->pending_searches + (num_pending_searches-1);
    neighborhood_search_hamming_extend(bidirectional_search,num_pending_searches-1,pending_search->direction,
        pending_search->offset,pending_search->offset+pending_search->length,0,0,pending_search->max_error);
  } else {
    // Compute error
    const uint64_t search_max_error = chunk_max_error/2;
    const uint64_t extend_max_error = chunk_max_error;
    // Search first partition (forward)
    const uint64_t first_half_offset = chunk_offset;
    const uint64_t first_half_length = chunk_length/2;
    const uint64_t second_half_offset = chunk_offset + first_half_length;
    const uint64_t second_half_length = chunk_length - first_half_length;
    // Pending searches
    const uint64_t num_pending_searches = bidirectional_search->num_pending_searches;
    // Schedule search
    bidirectional_search_push_operation(bidirectional_search,
        direction_forward,second_half_offset,second_half_length,extend_max_error);
    neighborhood_search_hamming_bidirectional_schedule__search(
        bidirectional_search,first_half_offset,first_half_length,search_max_error);
    bidirectional_search_error_mask_update_min(bidirectional_search,first_half_offset,first_half_length,search_max_error+1);
    // Search second partition (backward)
    bidirectional_search->num_pending_searches = num_pending_searches; // Restore pending searches point
    bidirectional_search_push_operation(bidirectional_search,
        direction_backward,first_half_offset,first_half_length,extend_max_error);
    neighborhood_search_hamming_bidirectional_schedule__search(
        bidirectional_search,second_half_offset,second_half_length,search_max_error);
  }
}
/*
 * Main
 */
GEM_INLINE void neighborhood_search_hamming_bidirectional(
    fm_index_t* const fm_index,uint8_t* const key,
    const uint64_t key_length,const uint64_t max_error,
    interval_set_t* const intervals_result,mm_stack_t* const mm_stack) {
  // Push stack state
  mm_stack_push_state(mm_stack);
  // Create bidirectional search
  bidirectional_search_t bidirectional_search;
  bidirectional_search_init(&bidirectional_search,fm_index,key,key_length,max_error,intervals_result,mm_stack);
  // Best-match
  neighborhood_search_hamming_bidirectional_schedule__search(&bidirectional_search,0,key_length,max_error);
  // Free used stack
  mm_stack_pop_state(mm_stack,false);
  // PROFILE
  fprintf(stderr," => Total nodes     %lu\n",_ns_nodes);
  fprintf(stderr," => Total solutions %lu\n",_ns_nodes_success);
}
/*
 * Brute force
 */
void neighborhood_search_hamming_brute_force_search(
    char* const search_string,uint8_t* const key,
    const uint64_t key_length,const uint64_t current_position,
    const uint64_t current_error,const uint64_t max_error) {
  // Expand node
  uint8_t enc;
  for (enc=0;enc<DNA_RANGE;++enc) {
    const uint64_t next_error = current_error + (enc!=(key[current_position]) ? 1 : 0);
    if (next_error <= max_error) {
      search_string[current_position] = dna_decode(enc);
      ++_ns_nodes; // PROFILING
      const uint64_t next_position = current_position + 1;
      if (next_position < key_length) {
        neighborhood_search_hamming_brute_force_search(
            search_string,key,key_length,next_position,next_error,max_error);
      } else {
        ++_ns_nodes_success;
        fprintf(stdout,"%s\n",search_string);
      }
    }
  }
}
void neighborhood_search_hamming_brute_force(uint8_t* const key,const uint64_t key_length,const uint64_t max_error) {
  // Init string
  char* const search_string = malloc(key_length+1);
  memset(search_string,0,key_length);
  // Search
  neighborhood_search_hamming_brute_force_search(search_string,key,key_length,0,0,max_error);
  // PROFILE
  fprintf(stderr," => Total nodes     %lu\n",_ns_nodes);
  fprintf(stderr," => Total solutions %lu\n",_ns_nodes_success);
}




