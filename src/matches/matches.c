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

#include "matches/matches.h"
#include "matches/matches_cigar.h"
#include "matches/classify/matches_classify.h"
#include "text/restriction_text.h"

/*
 * Constants
 */
#define MATCHES_INIT_INTERVAL_MATCHES  100
#define MATCHES_INIT_GLOBAL_MATCHES    100
#define MATCHES_INIT_CIGAR_OPS        5000

/*
 * Setup
 */
matches_t* matches_new(void) {
  // Allocate handler
  matches_t* const matches = mm_alloc(matches_t);
  // Matches Counters
  matches->counters = matches_counters_new();
  matches->max_complete_stratum = 0;
  matches->limited_exact_matches = false;
  matches->matches_extended = false;
  // MM
  matches->mm_slab = mm_slab_new_(BUFFER_SIZE_1M,BUFFER_SIZE_1M,MM_UNLIMITED_MEM);
  matches->mm_allocator = mm_allocator_new(matches->mm_slab);
  // Match-Traces
  matches->match_traces = vector_new(MATCHES_INIT_GLOBAL_MATCHES,match_trace_t*);
  matches->match_traces_local = vector_new(MATCHES_INIT_GLOBAL_MATCHES,match_trace_t*);
  matches->match_traces_extended = vector_new(MATCHES_INIT_GLOBAL_MATCHES,match_trace_t*);
  matches->match_traces_begin = ihash_new(matches->mm_allocator);
  matches->match_traces_end = ihash_new(matches->mm_allocator);
  matches->match_replaced = false;
  // Match-Blocks
  matches->match_blocks = vector_new(1, match_block_t);
  matches->match_total_covered_bases = 0;
  // Potential split sites
  matches->match_potential_split_sites = vector_new(1, uint64_t);
  // CIGAR buffer
  matches->cigar_vector = vector_new(MATCHES_INIT_CIGAR_OPS,cigar_element_t);
  // Init metrics
  matches_metrics_init(&matches->metrics);
  // Return
  return matches;
}
void matches_clear(matches_t* const matches) {
  matches->max_complete_stratum = 0;
  matches->limited_exact_matches = false;
  matches->matches_extended = false;
  matches_counters_clear(matches->counters);
  matches_metrics_init(&matches->metrics);
  vector_clear(matches->match_traces);
  vector_clear(matches->match_traces_local);
  vector_clear(matches->match_traces_extended);
  ihash_clear(matches->match_traces_begin);
  ihash_clear(matches->match_traces_end);
  vector_clear(matches->match_blocks);
  vector_clear(matches->match_potential_split_sites);
  matches->match_total_covered_bases = 0;
  matches->match_replaced = false;
  vector_clear(matches->cigar_vector);
  mm_allocator_clear(matches->mm_allocator);
}
void matches_delete(matches_t* const matches) {
  // Delete all
  matches_counters_delete(matches->counters);
  vector_delete(matches->match_traces);
  vector_delete(matches->match_traces_local);
  vector_delete(matches->match_traces_extended);
  ihash_delete(matches->match_traces_begin);
  ihash_delete(matches->match_traces_end);
  vector_delete(matches->cigar_vector);
  vector_delete(matches->match_blocks);
  vector_delete(matches->match_potential_split_sites);
  mm_allocator_delete(matches->mm_allocator);
  mm_slab_delete(matches->mm_slab);
  mm_free(matches);
}
/*
 * Accessors
 */
bool matches_is_mapped(const matches_t* const matches) {
  return matches_get_num_match_traces(matches) > 0;
}
void matches_recompute_metrics(matches_t* const matches) {
  // Parameters
  matches_metrics_t* const metrics = &matches->metrics;
  // ReCompute distance metrics
  match_trace_t** const match_traces = matches_get_match_traces(matches);
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  matches_metrics_init(metrics);
  uint64_t i;
  for (i=0;i<num_matches;++i) {
    const match_trace_t* const match = match_traces[i];
    matches_metrics_update(
        metrics,match->event_distance,
        match->edit_distance,match->swg_score);
  }
}
uint64_t matches_get_first_stratum_matches(matches_t* const matches) {
  const uint64_t min_edit_distance = matches->metrics.min1_edit_distance;
  return (min_edit_distance==UINT32_MAX) ? 0 : matches_counters_get_count(matches->counters,min_edit_distance);
}
uint64_t matches_get_subdominant_stratum_matches(matches_t* const matches) {
  const uint64_t first_stratum_matches = matches_get_first_stratum_matches(matches);
  return matches_counters_get_total_count(matches->counters) - first_stratum_matches;
}
uint8_t matches_get_primary_mapq(matches_t* const matches) {
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  if (num_matches == 0) return 0;
  return matches_get_primary_match(matches)->mapq_score;
}
void matches_update_mcs(
    matches_t* const matches,
    const uint64_t current_mcs) {
  matches->max_complete_stratum = MAX(matches->max_complete_stratum,current_mcs);
}
void matches_update_limited_exact_matches(
    matches_t* const matches,
    const uint64_t num_exact_matches_limited) {
  if (num_exact_matches_limited > 0) {
    matches->limited_exact_matches = true;
    matches_counters_add(matches->counters,0,num_exact_matches_limited);
    matches->max_complete_stratum = 0;
  }
}
/*
 * Index
 */
void matches_index_match(
    matches_t* const matches,
    match_trace_t* const match_trace) {
  // Compute dimensions
  const uint64_t begin_position = match_trace->match_alignment.match_position;
  const uint64_t effective_length = match_trace->match_alignment.effective_length;
  // Store begin/end positions of the match (as to fast index matches)
  ihash_insert(matches->match_traces_begin,begin_position,match_trace);
  ihash_insert(matches->match_traces_end,begin_position+effective_length,match_trace);
}
match_trace_t* matches_lookup_match(
    matches_t* const matches,
    const uint64_t begin_position,
    const uint64_t effective_length) {
  match_trace_t* const match_trace = ihash_get(matches->match_traces_begin,begin_position,match_trace_t);
  if (match_trace!=NULL) return match_trace;
  return ihash_get(matches->match_traces_end,begin_position+effective_length,match_trace_t);
}
/*
 * Matches
 */
void matches_clear_match_traces(const matches_t* const matches) {
  vector_clear(matches->match_traces);
}
uint64_t matches_get_num_match_traces(const matches_t* const matches) {
  return vector_get_used(matches->match_traces);
}
uint64_t matches_get_num_match_traces_extended(const matches_t* const matches) {
  return vector_get_used(matches->match_traces_extended);
}
match_trace_t* matches_get_primary_match(const matches_t* const matches) {
  return matches_get_match_traces(matches)[0];
}
match_trace_t** matches_get_match_traces(const matches_t* const matches) {
  return vector_get_mem(matches->match_traces,match_trace_t*);
}
/*
 * Sorting Matches
 */
void matches_traces_sort_by_genomic_position(
    match_trace_t** const match_traces,
    const uint64_t num_match_traces) {
  // Sort global matches (match_trace_t) wrt distance
  qsort(match_traces,num_match_traces,sizeof(match_trace_t*),
        (int (*)(const void *,const void *))match_trace_cmp_genomic_position);
}
void matches_traces_sort_by_primary_and_genomic_position(
    match_trace_t** const match_traces,
    const uint64_t num_match_traces) {
  // Sort local matches (match_trace_t) wrt primality and distance
  qsort(match_traces,num_match_traces,sizeof(match_trace_t*),
        (int (*)(const void *,const void *))match_trace_cmp_primary_and_genomic_position);
}
void matches_add_match_trace_insert_sorted(
    matches_t* const matches,
    match_trace_t* const match_trace) {
  // Reserve
  const uint64_t num_match_traces = matches_get_num_match_traces(matches);
  vector_reserve(matches->match_traces,num_match_traces+1,false);
  // Find insertion position
  match_trace_t** const match_traces = matches_get_match_traces(matches);
  int64_t i, j;
  for (i=0;i<num_match_traces;++i) {
    if (match_trace_cmp_swg_score(
        (const match_trace_t** const)&match_trace,
        (const match_trace_t** const)match_traces+i) < 0) break;
  }
  // Shift match-traces from i-th position
  for (j=num_match_traces;j>i;--j) {
    match_traces[j] = match_traces[j-1];
  }
  // Insert at i-th position
  match_traces[i] = match_trace;
  vector_add_used(matches->match_traces,1);
}
void matches_add_match_trace_preserve_sorted(matches_t* const matches) {
  // Pseudo-bubble sort
  const uint64_t num_match_traces = matches_get_num_match_traces(matches) - 1;
  match_trace_t** const match_traces = matches_get_match_traces(matches);
  uint64_t i, j;
  for (i=0;i<num_match_traces;++i) {
    if (match_trace_cmp_swg_score(
        (const match_trace_t** const)match_traces+i,
        (const match_trace_t** const)match_traces+(i+1)) > 0) {
      j = i+1;
      do {
        SWAP(match_traces[j-1],match_traces[j]);
        --j;
      } while (j>0 && match_trace_cmp_swg_score(
          (const match_trace_t** const)match_traces+(j-1),
          (const match_trace_t** const)match_traces+j) > 0);
    }
  }
}
/*
 * Adding Match-traces
 */
void match_trace_locate(
    match_trace_t* const match_trace,
    const locator_t* const locator) {
  GEM_INTERNAL_CHECK(
		  match_trace->match_alignment.effective_length >= 0,
		  "Match effective length must be positive");
  location_t location;
  locator_map(locator,match_trace->match_alignment.match_position,&location);
  match_trace->text_position = location.position;
  match_trace->sequence_name = location.tag;
  match_trace->tag_id = location.tag_id;
  if (location.strand == Reverse) { // Adjust position by the effective length
    match_trace->text_position -= (uint64_t) match_trace->match_alignment.effective_length;
    match_trace->strand = Reverse;
  } else {
    match_trace->strand = Forward;
  }
  match_trace->bs_strand = location.bs_strand;
}
void match_trace_copy(
    match_trace_t* const match_trace_dst,
    match_trace_t* const match_trace_src) {
  /* Text (Reference) */
  match_trace_dst->type = match_trace_src->type;
  match_trace_dst->text_trace = match_trace_src->text_trace;
  match_trace_dst->text = match_trace_src->text;
  match_trace_dst->text_length = match_trace_src->text_length;
  /* Match */
  match_trace_dst->sequence_name = match_trace_src->sequence_name;
  match_trace_dst->strand = match_trace_src->strand;
  match_trace_dst->bs_strand = match_trace_src->bs_strand;
  match_trace_dst->text_position = match_trace_src->text_position;
  /* Score */
  match_trace_dst->edit_distance = match_trace_src->edit_distance;
  /* This is recalculated for local matches, otherwise use global score */
  match_trace_dst->edit_distance_local = match_trace_src->edit_distance;

  match_trace_dst->event_distance = match_trace_src->event_distance;
  match_trace_dst->swg_score = match_trace_src->swg_score;
  match_trace_dst->error_quality = match_trace_src->error_quality;
  match_trace_dst->mapq_score = match_trace_src->mapq_score;
  /* Alignment */
  match_trace_dst->match_alignment = match_trace_src->match_alignment;
  match_trace_dst->match_scaffold = match_trace_src->match_scaffold;
  /* Block */
  match_trace_dst->match_block_index = match_trace_src->match_block_index;
  /* Restriction site indices */
  match_trace_dst->first_restriction_site = match_trace_src->first_restriction_site;
  match_trace_dst->last_restriction_site = match_trace_src->last_restriction_site;
}
void match_trace_get_block(
		matches_t* const matches,
		match_trace_t* const match_trace,
		uint64_t *bounds) {
  const uint64_t cigar_length = match_trace->match_alignment.cigar_length;
  const cigar_element_t* const cigar_array =
  		vector_get_mem(matches->cigar_vector,cigar_element_t) + match_trace->match_alignment.cigar_offset;
  uint64_t lo = 0, hi = 0, i = 0, right_skip = 0;
  if(cigar_length > 0) {
  	if(cigar_array[0].type==cigar_del && cigar_array[0].attributes == cigar_attr_trim) {
  		lo += cigar_array[0].length;
  		i++;
  	}
  	hi = lo;
  	for(;i<cigar_length;i++) {
  		switch(cigar_array[i].type) {
  		case cigar_match:
  			hi += cigar_array[i].length;
  			break;
  		case cigar_mismatch:
  			hi ++;
  			break;
  		case cigar_del:
  			if(cigar_array[i].attributes == cigar_attr_trim) right_skip = cigar_array[i].length;
  			else {
  				hi += cigar_array[i].length;
  			}
  			break;
  		case cigar_ins:
  			break;
  		default:
  		  GEM_INVALID_CASE();
  		  break;
  		}
  	}
  }
	bounds[0] = lo;
	bounds[1] = hi;
	bounds[2] = hi + right_skip;
}

void copy_match_block(match_block_t * const dst, match_block_t * const src) {
	dst->hi = src->hi;
	dst->lo = src->lo;
	dst->overlaps = src->overlaps;
	dst->clip_right = src->clip_right;
}

void matches_add_match_trace_to_blocks(
		search_parameters_t* const search,
    matches_t* const matches,
    match_trace_t* const match_trace) {

	uint64_t max_overlap = search->split_map_max_overlap;
	uint64_t bounds[3];
	match_trace_get_block(matches, match_trace, bounds);
	uint64_t num_match_blocks = vector_get_used(matches->match_blocks);
	match_block_t *blocks = vector_get_mem(matches->match_blocks, match_block_t);
	uint64_t insert_pos = 0;
//	fprintf(stderr,"bounds = [%"PRIu64", %"PRIu64"], num_match_blocks = %"PRIu64"\n", bounds[0], bounds[1], num_match_blocks);
	uint64_t i;
	bool overlaps = false;
	for(i = 0; i < num_match_blocks; i++) {
		if((bounds[0] >= blocks[i].hi || blocks[i].hi - bounds[0] < max_overlap)
			 || (bounds[1] <= blocks[i].lo || bounds[1] - blocks[i].lo < max_overlap)) { // Almost non-overlapping blocks
			blocks[i].overlaps = false;
		} else {
			if(bounds[0] >= blocks[i].lo) {
				bounds[0] = blocks[i].lo;
			}
			if(bounds[1] <= blocks[i].hi) {
				bounds[1] = blocks[i].hi;
			}
			blocks[i].overlaps = true;
			overlaps = true;
		}
	}
//	fprintf(stderr,"...bounds = [%"PRIu64", %"PRIu64"] overlaps:%c\n", bounds[0], bounds[1], overlaps ? 'T' : 'F');
	int64_t* new_index = num_match_blocks ? mm_allocator_calloc(matches->mm_allocator,num_match_blocks, int64_t, false) : NULL;
	if(!overlaps) {
		vector_reserve(matches->match_blocks, num_match_blocks + 1, false);
		blocks = vector_get_mem(matches->match_blocks, match_block_t);
	}
	uint64_t j;
	for(i = j = 0; i < num_match_blocks; i++) {
		if(!blocks[i].overlaps) {
			if(blocks[i].lo > bounds[0]) break;
			new_index[i] = j++;
		} else new_index[i] = -1;
	}
	insert_pos = j++;
	match_block_t block, block1;

	if(i < num_match_blocks) copy_match_block(&block, blocks + i);
	for(;i < num_match_blocks; i++) {
		if(i < num_match_blocks - 1) copy_match_block(&block1, blocks + i + 1);
		if(!block.overlaps) {
			new_index[i] = j;
			copy_match_block(blocks + (j++),&block);
		} else new_index[i] = -1;
		copy_match_block(&block, &block1);
	}
	vector_set_used(matches->match_blocks,j);
	blocks[insert_pos].lo = bounds[0];
	blocks[insert_pos].hi = bounds[1];
	blocks[insert_pos].clip_right = bounds[2] - bounds[1];

//	matches->match_total_covered_bases = 0;
//	for(i = 0; i < j; i++) {
//		if(!i || blocks[i-1].hi <= blocks[i].lo) matches->match_total_covered_bases += blocks[i].hi - blocks[i].lo;
//		else matches->match_total_covered_bases += blocks[i].hi - blocks[i-1].hi;
//	}
//	matches->match_total_covered_bases = bounds[2];
	matches->match_total_covered_bases = matches->key_length;

//  fprintf(stderr,"covered bases = %"PRIu64", insert_pos = %"PRIu64"\n", matches->match_total_covered_bases, insert_pos);

	uint64_t num_match_traces = matches_get_num_match_traces(matches);
	match_trace_t **traces = matches_get_match_traces(matches);
	for(uint64_t i = 0; i < num_match_traces; i++) {
		uint64_t old_block_index = traces[i]->match_block_index;
		traces[i]->match_block_index = new_index[old_block_index] < 0 ? insert_pos : (uint64_t)new_index[old_block_index];
	}
	if(new_index != NULL) mm_allocator_free(matches->mm_allocator, new_index);
	match_trace->match_block_index=insert_pos;
}

match_trace_t* matches_add_match_trace(
		search_parameters_t* const search,
    matches_t* const matches,
    const locator_t* const locator,
    match_trace_t* const match_trace,
    bool* const match_replaced) {

// Init/Clean
match_trace->mapq_score = 0;
  // Check duplicates
  match_trace_t* const match_trace_dup = matches_lookup_match(
      matches,match_trace->match_alignment.match_position,
      match_trace->match_alignment.effective_length);
  if (match_trace_dup==NULL) {
    // Process match-trace (Correct CIGAR & locate)
    match_trace_locate(match_trace,locator);
    restriction_match_trace_locate(match_trace, &search->restriction_site_locator);
//  	fprintf(stderr,"Adding map trace: %s %lu\n", match_trace->sequence_name, match_trace->text_position);
    if(search->split_map_score) {
    	matches_add_match_trace_to_blocks(search, matches, match_trace);
    }
    // Find insertion position
    // Add the match-trace
    PROF_INC_COUNTER(GP_MATCHES_MAPS_ADDED);
    match_trace_t* const new_match_trace = mm_allocator_alloc(matches->mm_allocator,match_trace_t);
    match_trace_copy(new_match_trace,match_trace);

    // Index
    matches_index_match(matches,new_match_trace);

    // Update counters
    matches_counters_add(matches->counters,match_trace->edit_distance,1);
    // Update metrics
    matches_metrics_update(
        &matches->metrics,match_trace->event_distance,
        match_trace->edit_distance,match_trace->swg_score);
    // Insert preserving sorting
    matches_add_match_trace_insert_sorted(matches,new_match_trace);
    *match_replaced = false;
		uint64_t num_match_traces = matches_get_num_match_traces(matches);
		match_trace_t **traces = matches_get_match_traces(matches);
		match_block_t *blocks = vector_get_mem(matches->match_blocks, match_block_t);
		if(search->split_map) for(uint64_t i = 0; i < num_match_traces; i++) {
    	uint64_t bounds[3];
    	match_trace_get_block(matches, traces[i], bounds);
			uint64_t block_index = traces[i]->match_block_index;
//			fprintf(stderr,"..trace %"PRIu64" (%"PRIu64"-%"PRIu64") : block %"PRIu64" (%"PRIu64"-%"PRIu64") %p\n", i, bounds[0], bounds[1], block_index, blocks[block_index].lo, blocks[block_index].hi, blocks + block_index);
			gem_cond_error_msg(bounds[0] < blocks[block_index].lo || bounds[1] > blocks[block_index].hi,
					"Trace not contained in block (trace: %"PRIu64"-%"PRIu64", (block: %"PRIu64"-%"PRIu64")", bounds[0], bounds[1], blocks[block_index].lo, blocks[block_index].hi);
		}

    return new_match_trace;
  } else {
    // Pick up the least distant match
    PROF_INC_COUNTER(GP_MATCHES_MAPS_DUP);
    if (match_trace_dup->swg_score < match_trace->swg_score) {
//			fprintf(stderr,"Duplicates\n");
      // Process match-trace (Correct CIGAR & locate)
      match_trace_locate(match_trace,locator);
      restriction_match_trace_locate(match_trace, &search->restriction_site_locator);
			match_trace->match_block_index = match_trace_dup->match_block_index;
      // Update counters
      matches_counters_sub(matches->counters,match_trace_dup->edit_distance,1);
      matches_counters_add(matches->counters,match_trace->edit_distance,1);
      // Replace the match-trace
      match_trace_copy(match_trace_dup,match_trace);
      matches_add_match_trace_preserve_sorted(matches);
      // Recompute metrics
      matches_recompute_metrics(matches);
			// Recompute match_blocks
			if(search->split_map_score) {
				matches_add_match_trace_to_blocks(search, matches, match_trace_dup);
				uint64_t num_match_traces = matches_get_num_match_traces(matches);
				match_trace_t **traces = matches_get_match_traces(matches);
				match_block_t *blocks = vector_get_mem(matches->match_blocks, match_block_t);
				for(uint64_t i = 0; i < num_match_traces; i++) {
					uint64_t bounds[3];
					match_trace_get_block(matches, traces[i], bounds);
					uint64_t block_index = traces[i]->match_block_index;
				//	fprintf(stderr,"..trace %"PRIu64" (%"PRIu64"-%"PRIu64") : block %"PRIu64" (%"PRIu64"-%"PRIu64") %p\n", i, bounds[0], bounds[1], block_index, blocks[block_index].lo, blocks[block_index].hi, blocks + block_index);
					gem_cond_error_msg(bounds[0] < blocks[block_index].lo || bounds[1] > blocks[block_index].hi,
										"Trace not contained in block (trace: %"PRIu64"-%"PRIu64", (block: %"PRIu64"-%"PRIu64")", bounds[0], bounds[1], blocks[block_index].lo, blocks[block_index].hi);
				} 
			}
      // Return
      *match_replaced = true;
      return match_trace_dup;
    } else {
      // Return
      *match_replaced = false;
      return NULL;
    }
  }
}
void matches_add_match_trace_extended(
		search_parameters_t* const search,
    matches_t* const matches,
    const locator_t* const locator,
    match_trace_t* const match_trace) {
  bool match_replaced;
  match_trace_t* const match_trace_added =
      matches_add_match_trace(search,matches,locator,match_trace,&match_replaced);
  if (match_trace_added!=NULL) {
    match_trace_added->type = match_type_extended;
    if (match_replaced) {
      matches->match_replaced = true;
    } else {
      vector_insert(matches->match_traces_extended,match_trace_added,match_trace_t*);
    }
  }
}
void matches_add_match_trace_local_pending(
    matches_t* const matches,
    match_trace_t* const match_trace) {
  match_trace_t* const new_match_trace = mm_allocator_alloc(matches->mm_allocator,match_trace_t);
  match_trace_copy(new_match_trace,match_trace);
  vector_insert(matches->match_traces_local,new_match_trace,match_trace_t*);
}
void matches_local_pending_add_to_regular_matches(
		search_parameters_t* const search,
    matches_t* const matches,
    const locator_t* const locator) {
  const uint64_t num_match_traces_local = vector_get_used(matches->match_traces_local);
  match_trace_t** const match_traces_local = vector_get_mem(matches->match_traces_local,match_trace_t*);
  bool match_replaced;
  uint64_t i;
  for (i=0;i<num_match_traces_local;++i) {
    match_trace_t* const match_traces = match_traces_local[i];
    matches_add_match_trace(search,matches,locator,match_traces,&match_replaced);
  }
  vector_clear(matches->match_traces_local);
}
void matches_local_pending_add_to_extended_matches(
    matches_t* const matches,
    const locator_t* const locator) {
  const uint64_t num_match_traces_local = vector_get_used(matches->match_traces_local);
  match_trace_t** const match_traces_local = vector_get_mem(matches->match_traces_local,match_trace_t*);
  uint64_t i;
  for (i=0;i<num_match_traces_local;++i) {
    match_trace_t* const match_traces = match_traces_local[i];
    match_traces->type = match_type_extended;
    vector_insert(matches->match_traces_extended,match_traces,match_trace_t*);
  }
  vector_clear(matches->match_traces_local);
}
/*
 * Display
 */
void matches_print(
    FILE* const stream,
    matches_t* const matches) {
  tab_fprintf(stream,"[GEM]>Matches\n");
  tab_global_inc();
  tab_fprintf(stream,"=> Counters\t");
  matches_counters_print(stream,matches->counters,matches->max_complete_stratum);
  fprintf(stream,"\n");
  tab_fprintf(stream,"=> Position.Matches %lu\n",matches_get_num_match_traces(matches));
  tab_fprintf(stream,"  => Positions.Hashed.Begin %lu\n",ihash_get_num_elements(matches->match_traces_begin));
  tab_fprintf(stream,"  => Positions.Hashed.End %lu\n",ihash_get_num_elements(matches->match_traces_end));
  tab_fprintf(stream,"=> Metrics.SE\n");
  tab_global_inc();
  matches_metrics_print(stream,&matches->metrics);
  tab_global_dec();
  tab_global_dec();
}


