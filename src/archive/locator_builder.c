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
 *   Locator builder provides functions to create a locator for the GEM3 index
 */

#include "archive/locator_builder.h"
#include "archive/locator.h"

/*
 * Setup
 */
locator_builder_t* locator_builder_new(mm_slab_t* const mm_slab) {
  // Allocate
  locator_builder_t* const locator_builder = mm_alloc(locator_builder_t);
  // Intervals
  locator_builder->current_interval = NULL;
  locator_builder->intervals = svector_new(mm_slab,locator_interval_t);
  svector_iterator_new(&(locator_builder->intervals_iterator),locator_builder->intervals,SVECTOR_WRITE_ITERATOR,0);
  // Tag Locator
  locator_builder->tags_dictionary = shash_new(NULL);
  locator_builder->tag_locator = svector_new(mm_slab,locator_tag_t);
  svector_iterator_new(&(locator_builder->tag_locator_iterator),locator_builder->tag_locator,SVECTOR_WRITE_ITERATOR,0);
  // Tags Buffer
  locator_builder->tags_buffer = svector_new(mm_slab,char);
  locator_builder->tags_buffer_size = 0;
  // Builder Position
  locator_builder->index_position = 0;
  locator_builder->sequence_offset = 0;
  // Builder
  locator_builder->tag_id_generator = 0;
  // Return
  return locator_builder;
}
void locator_builder_delete(locator_builder_t* const locator_builder) {
  // Free
  svector_delete(locator_builder->intervals);
  svector_delete(locator_builder->tag_locator);
  svector_delete(locator_builder->tags_buffer);
  shash_delete(locator_builder->tags_dictionary);
  mm_free(locator_builder);
}
void locator_builder_write(
    fm_t* const file_manager,
    locator_builder_t* const locator_builder) {
  // Write Locator header
  fm_write_uint64(file_manager,LOCATOR_MODEL_NO);
  fm_write_uint64(file_manager,svector_get_used(locator_builder->intervals)); // num_intervals
  fm_write_uint64(file_manager,svector_get_used(locator_builder->tag_locator)); // num_tags
  fm_write_uint64(file_manager,locator_builder->tags_buffer_size); // tags_buffer_size
  // Write tags
  svector_iterator_t tag_locator_iterator;
  svector_iterator_new(&tag_locator_iterator,locator_builder->tag_locator,SVECTOR_READ_ITERATOR,0);
  uint64_t tags_length_written = 0;
  while (!svector_read_iterator_eoi(&tag_locator_iterator)) {
    // Get next locator tag
    locator_tag_t* const locator_tag = svector_iterator_get_element(&tag_locator_iterator,locator_tag_t);
    svector_read_iterator_next(&tag_locator_iterator);
    // Write the tag
    const char* const tag = locator_builder_get_tag(locator_builder,locator_tag);
    fm_write_mem(file_manager,tag,locator_tag->length+1);
    // Adjust the offset
    locator_tag->offset = tags_length_written;
    tags_length_written += locator_tag->length+1;
  }
  // Write tag-locator
  svector_write(file_manager,locator_builder->tag_locator);
  // Write intervals
  svector_write(file_manager,locator_builder->intervals);
}
/*
 * Locator Builder Accessors
 */
uint64_t locator_builder_get_num_intervals(locator_builder_t* const locator_builder) {
  return svector_get_used(locator_builder->intervals);
}
uint64_t locator_builder_get_num_tags(locator_builder_t* const locator_builder) {
  return locator_builder->tag_id_generator;
}
/*
 * Skip text/index
 */
void locator_builder_skip_index(
    locator_builder_t* const locator_builder,
    const uint64_t length) {
  locator_builder->index_position += length;
}
/*
 * Locator Builder Interval Accessors
 */
locator_interval_t* locator_builder_get_interval(
    locator_builder_t* const locator_builder,
    const uint64_t interval_index) {
  return svector_get_element(locator_builder->intervals,interval_index,locator_interval_t);
}
/*
 * Locator Builder Tag Accessors
 */
char* locator_builder_get_tag(
    locator_builder_t* const locator_builder,
    locator_tag_t* const locator_tag) {
  return svector_get_element(locator_builder->tags_buffer,locator_tag->offset,char);
}
locator_tag_t* locator_builder_get_locator_tag_by_id(
    locator_builder_t* const locator_builder,
    const int64_t tag_id) {
  return svector_get_element(locator_builder->tag_locator,tag_id,locator_tag_t);
}
locator_tag_t* locator_builder_get_locator_tag_by_tag(
    locator_builder_t* const locator_builder,
    char* const tag) {
  return shash_get(locator_builder->tags_dictionary,tag,locator_tag_t);
}
/*
 * Intervals handling
 */
void locator_builder_open_interval(
    locator_builder_t* const locator_builder,
    const int64_t tag_id) {
  // Open a new interval => [begin_position,end_position)
  svector_iterator_t* const intervals_iterator = &(locator_builder->intervals_iterator);
  locator_builder->current_interval = svector_iterator_get_element(intervals_iterator,locator_interval_t);
  locator_builder->current_interval->tag_id = tag_id;
  locator_builder->current_interval->begin_position = locator_builder->index_position;
}
void locator_builder_close_interval(
    locator_builder_t* const locator_builder,
    const uint64_t interval_text_length,
    const uint64_t interval_index_length,
    const locator_interval_type type,
    const strand_t strand,
    const bs_strand_t bs_strand) {
  // Close interval & write it => [begin_position,end_position)
  locator_builder->current_interval->type = type;
  locator_builder->current_interval->strand = strand;
  locator_builder->current_interval->bs_strand = bs_strand;
  locator_builder->current_interval->end_position =
      locator_builder->current_interval->begin_position+interval_index_length;
  locator_builder->current_interval->rl_begin_position = 0;
  locator_builder->current_interval->rl_end_position = 0;
  locator_builder->current_interval->sequence_offset = locator_builder->sequence_offset;
  locator_builder->current_interval->sequence_length = interval_text_length;
  svector_write_iterator_next(&(locator_builder->intervals_iterator)); // Next
  // Update builder state
  locator_builder->index_position += interval_index_length;
  locator_builder->sequence_offset += interval_text_length;
}
void locator_builder_add_interval(
    locator_builder_t* const locator_builder,
    const int64_t tag_id,
    const uint64_t sequence_offset,
    const uint64_t interval_text_length,
    const uint64_t interval_index_length,
    const locator_interval_type type,
    const strand_t strand,
    const bs_strand_t bs_strand) {
  // Add interval
  locator_builder->current_interval =
      svector_iterator_get_element(&(locator_builder->intervals_iterator),locator_interval_t);
  locator_builder->current_interval->type = type;
  locator_builder->current_interval->strand = strand;
  locator_builder->current_interval->bs_strand = bs_strand;
  locator_builder->current_interval->begin_position = locator_builder->index_position;
  locator_builder->current_interval->end_position = locator_builder->index_position+interval_index_length;
  locator_builder->current_interval->rl_begin_position = 0;
  locator_builder->current_interval->rl_end_position = 0;
  locator_builder->current_interval->sequence_offset = sequence_offset;
  locator_builder->current_interval->sequence_length = interval_text_length;
  locator_builder->current_interval->tag_id = (tag_id==-1) ? locator_builder->tag_id_generator : tag_id;
  svector_write_iterator_next(&(locator_builder->intervals_iterator));
  // Update builder state
  locator_builder->index_position += interval_index_length;
}
void locator_builder_add_rc_interval(
    locator_builder_t* const locator_builder,
    locator_interval_t* const locator_interval) {
  const uint64_t interval_length = locator_interval_get_index_length(locator_interval);
  const strand_t rc_strand = dna_text_strand_get_complement(locator_interval->strand);
  locator_builder_add_interval(locator_builder,locator_interval->tag_id,locator_interval->sequence_offset,
      locator_interval->sequence_length,interval_length,locator_interval->type,rc_strand,locator_interval->bs_strand);
}
int64_t locator_builder_add_sequence(
    locator_builder_t* const locator_builder,
    char* const tag,
    const uint64_t tag_length) {
  // Reset the sequence offset
  locator_builder->sequence_offset = 0;
  // Check if its already in the dictionary of tags
  locator_tag_t* current_locator_tag = locator_builder_get_locator_tag_by_tag(locator_builder,tag);
  if (current_locator_tag!=NULL) {
    return current_locator_tag->tag_id; // Return found TAG
  } else {
    // Setup the current locator_tag
    current_locator_tag = svector_iterator_get_element(&(locator_builder->tag_locator_iterator),locator_tag_t);
    current_locator_tag->tag_id = (locator_builder->tag_id_generator)++;
    // Add it to the tag_buffer
    char* const tag_in_buffer = svector_insert_char_buffer(
        locator_builder->tags_buffer,&current_locator_tag->offset,tag,tag_length);
    locator_builder->tags_buffer_size += tag_length+1;
    current_locator_tag->length = tag_length;
    svector_write_iterator_next(&(locator_builder->tag_locator_iterator)); // Next
    // Add it to the dictionary
    shash_insert(locator_builder->tags_dictionary,tag_in_buffer,tag_length,current_locator_tag);
    // Return added TAG
    return current_locator_tag->tag_id;
  }
}
/*
 * Display
 */
void locator_builder_print(
    FILE* const stream,
    locator_builder_t* const locator_builder,
    const bool display_intervals) {
  GEM_CHECK_NULL(stream);
  // Print locator info
  tab_fprintf(stream,"[GEM]>Locator\n");
  locator_print_summary(stream,
      svector_get_used(locator_builder->intervals),
      svector_get_used(locator_builder->tag_locator),
      locator_builder->tags_buffer_size);
  // Print intervals
  if (display_intervals) {
    tab_fprintf(stream,"  => Locator.intervals\n");
    uint64_t i = 0;
    svector_iterator_t* const intervals_iterator = &(locator_builder->intervals_iterator);
    svector_iterator_new(intervals_iterator,locator_builder->intervals,SVECTOR_READ_ITERATOR,0);
    while (!svector_read_iterator_eoi(intervals_iterator)) {
      // Get interval & tag
      locator_interval_t* const interval = svector_iterator_get_element(intervals_iterator,locator_interval_t);
      locator_tag_t* const locator_tag = locator_builder_get_locator_tag_by_id(locator_builder,interval->tag_id);
      const char* const interval_tag = locator_builder_get_tag(locator_builder,locator_tag);
      // Print interval information
      tab_fprintf(stream,"       [%"PRIu64"] ",i++);
      locator_interval_print(stream,interval,interval_tag);
      // Next!
      svector_read_iterator_next(intervals_iterator);
    }
  }
  // Flush
  fflush(stream);
}
