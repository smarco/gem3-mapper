/*
 * PROJECT: GEMMapper
 * FILE: locator_builder.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef LOCATOR_BUILDER_H_
#define LOCATOR_BUILDER_H_

#include "utils/essentials.h"
#include "archive/locator.h"

/*
 * Locator builder (To create a new locator structure)
 */
typedef struct {
  /* Intervals */
  svector_t* intervals;                    // Stores all intervals (locator_interval_t)
  svector_iterator_t intervals_iterator;   // Interval write iterator
  /* Tags Locator */
  shash_t* tags_dictionary;                // Tag's Dictionary (locator_tag_t*)
  svector_t* tag_locator;                  // Tag_Locator Array (locator_tag_t)
  svector_iterator_t tag_locator_iterator; // Tag_Locator write iterator
  /* Tags Buffer */
  svector_t* tags_buffer;                  // Stores all the tags (char)
  uint64_t tags_buffer_size;               // Current sum of all the length of the tags stored (Note. tags_buffer->size != tags_buffer_size)
  /* Builder */
  locator_interval_t* current_interval;    // Current Interval
  int64_t tag_id_generator;
  uint64_t index_position;                 // Current position in the index
  uint64_t sequence_offset;                // Current position in the text (Offset wrt the current sequence)
} locator_builder_t;

/*
 * Setup
 */
locator_builder_t* locator_builder_new(mm_slab_t* const mm_slab);
void locator_builder_delete(locator_builder_t* const locator_builder);
void locator_builder_write(
    fm_t* const file_manager,
    locator_builder_t* const locator_builder);

/*
 * Locator Builder Accessors
 */
uint64_t locator_builder_get_num_intervals(locator_builder_t* const locator_builder);
uint64_t locator_builder_get_num_tags(locator_builder_t* const locator_builder);

/*
 * Skip text/index
 */
void locator_builder_skip_index(
    locator_builder_t* const locator_builder,
    const uint64_t length);

/*
 * Locator Builder Interval Accessors
 */
locator_interval_t* locator_builder_get_interval(
    locator_builder_t* const locator_builder,
    const uint64_t interval_index);

/*
 * Locator Builder Tag Accessors
 */
char* locator_builder_get_tag(
    locator_builder_t* const locator_builder,
    locator_tag_t* const locator_tag);
locator_tag_t* locator_builder_get_locator_tag_by_id(
    locator_builder_t* const locator_builder,
    const int64_t tag_id);
locator_tag_t* locator_builder_get_locator_tag_by_tag(
    locator_builder_t* const locator_builder,
    char* const tag);

/*
 * Interval handling
 */
void locator_builder_open_interval(
    locator_builder_t* const locator_builder,
    const int64_t tag_id);
void locator_builder_close_interval(
    locator_builder_t* const locator_builder,
    const uint64_t interval_text_length,
    const uint64_t interval_index_length,
    const locator_interval_type type,
    const strand_t strand,
    const bs_strand_t bs_strand);
void locator_builder_add_interval(
    locator_builder_t* const locator_builder,
    const int64_t tag_id,
    const uint64_t sequence_offset,
    const uint64_t interval_text_length,
    const uint64_t interval_index_length,
    const locator_interval_type type,
    const strand_t strand,
    const bs_strand_t bs_strand);
void locator_builder_add_rc_interval(
    locator_builder_t* const locator_builder,
    locator_interval_t* const locator_interval);
int64_t locator_builder_add_sequence(
    locator_builder_t* const locator_builder,
    char* const tag,
    const uint64_t tag_length);

/*
 * Display
 */
void locator_builder_print(
    FILE* const stream,
    locator_builder_t* const locator_builder,
    const bool display_intervals);


#endif /* LOCATOR_BUILDER_H_ */
