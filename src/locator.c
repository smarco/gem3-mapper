/*
 * PROJECT: GEMMapper
 * FILE: locator.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides functionality as to recover sequence location from an archive of sequences
 */

#include "locator.h"

/*
 * FM-Inde Model & Version
 */
#define LOCATOR_MODEL_NO  3000ul

/*
 * Constants
 */
#define LOCATOR_NUM_INITIAL_TAGS 128

/*
 * Setup
 */
#define LOCATOR_CALCULATE_SIZES(locator) \
  const uint64_t intervals_locator_size = locator->num_intervals*sizeof(locator_interval_t); \
  const uint64_t tag_locator_size = locator->num_tags*sizeof(locator_tag_t); \
  const uint64_t tag_buffer_size = locator->tags_buffer_size;
GEM_INLINE locator_t* locator_read_mem(mm_t* const memory_manager) {
  MM_CHECK(memory_manager);
  // Allocate locator
  locator_t* const locator = mm_alloc(locator_t);
  // Read locator
  const uint64_t locator_model_no = mm_read_uint64(memory_manager);
  gem_cond_fatal_error(locator_model_no!=LOCATOR_MODEL_NO,LOCATOR_WRONG_MODEL_NO,locator_model_no,LOCATOR_MODEL_NO);
  locator->num_intervals = mm_read_uint64(memory_manager);
  locator->num_tags = mm_read_uint64(memory_manager);
  locator->tags_buffer_size = mm_read_uint64(memory_manager);
  // Calculate sizes
  LOCATOR_CALCULATE_SIZES(locator);
  // Null Locator memory
  locator->mm = NULL;
  // Read intervals
  locator->intervals = mm_read_mem(memory_manager,intervals_locator_size);
  // Read tags
  locator->tags_dictionary = NULL;
  locator->tag_locator = mm_read_mem(memory_manager,tag_locator_size);
  locator->tags_buffer = mm_read_mem(memory_manager,tag_buffer_size);
  // Return
  return locator;
}
GEM_INLINE void locator_setup_inverse_locator(locator_t* const locator) {
  LOCATOR_CHECK(locator);
  // TODO
}
GEM_INLINE void locator_delete(locator_t* const locator) {
  LOCATOR_CHECK(locator);
  // Free memory managers
  if (locator->mm) mm_bulk_free(locator->mm);
  // Free handler
  mm_free(locator);
}
/*
 * Accessors/Utils [Locator]
 */
GEM_INLINE uint64_t locator_get_size(locator_t* const locator) {
  LOCATOR_CHECK(locator);
  const uint64_t intervals_locator_size = locator->num_intervals*sizeof(locator_interval_t);
  const uint64_t tag_locator_size = locator->num_tags*sizeof(locator_tag_t);
  const uint64_t tag_buffer_size = locator->tags_buffer_size;
  return intervals_locator_size + tag_locator_size + tag_buffer_size;
}
GEM_INLINE locator_interval_t* locator_get_interval(const locator_t* const locator,const uint64_t interval_index) {
  LOCATOR_CHECK(locator);
  gem_fatal_check(interval_index >= locator->num_intervals,
      LOCATOR_INTERVAL_INDEX_OOB,interval_index,0ul,locator->num_intervals);
  return locator->intervals + interval_index;
}
/*
 * Accessors/Utils [Locator Interval]
 */
GEM_INLINE char* locator_interval_get_tag(const locator_t* const locator,const locator_interval_t* const interval) {
  LOCATOR_CHECK(locator);
  LOCATOR_INTERVAL_CHECK(interval);
  // return locator->tag_locator[ABS(interval->tag_id)-1].tag;
  // TODO: improve this incoherence
  //       OR enhance and load at start
  return locator->tags_buffer+locator->tag_locator[ABS(interval->tag_id)-1].offset;
}
GEM_INLINE uint64_t locator_interval_get_index_length(const locator_interval_t* const interval) {
  LOCATOR_INTERVAL_CHECK(interval);
  return (interval->end_position - interval->begin_position);
}
GEM_INLINE uint64_t locator_interval_get_text_length(const locator_interval_t* const interval) {
  LOCATOR_INTERVAL_CHECK(interval);
  return interval->sequence_length;
}
GEM_INLINE bool locator_interval_is_forward(const locator_interval_t* const interval) {
  LOCATOR_INTERVAL_CHECK(interval);
  return (interval->tag_id > 0);
}
/*
 * Builder
 */
GEM_INLINE locator_builder_t* locator_builder_new(mm_slab_t* const mm_slab) {
  MM_SLAB_CHECK(mm_slab);
  // Allocate
  locator_builder_t* const locator_builder = mm_alloc(locator_builder_t);
  // Intervals
  locator_builder->current_interval = NULL;
  locator_builder->intervals = svector_new(mm_slab,locator_interval_t);
  svector_iterator_new(&(locator_builder->intervals_iterator),locator_builder->intervals,SVECTOR_WRITE_ITERATOR,0);
  // Tag Locator
  locator_builder->tags_dictionary = shash_new();
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
GEM_INLINE void locator_builder_delete(locator_builder_t* const locator_builder) {
  LOCATOR_BUILDER_CHECK(locator_builder);
  // Free
  svector_delete(locator_builder->intervals);
  svector_delete(locator_builder->tag_locator);
  svector_delete(locator_builder->tags_buffer);
  shash_delete(locator_builder->tags_dictionary);
  mm_free(locator_builder);
}
GEM_INLINE void locator_builder_write(fm_t* const file_manager,locator_builder_t* const locator_builder) {
  FM_CHECK(file_manager);
  LOCATOR_BUILDER_CHECK(locator_builder);
  // Write Locator header
  fm_write_uint64(file_manager,LOCATOR_MODEL_NO);
  fm_write_uint64(file_manager,svector_get_used(locator_builder->intervals)); // num_intervals
  fm_write_uint64(file_manager,svector_get_used(locator_builder->tag_locator)); // num_tags
  fm_write_uint64(file_manager,locator_builder->tags_buffer_size); // tags_buffer_size
  // Write Intervals
  svector_write(file_manager,locator_builder->intervals);
  // Translation table
  svector_write(file_manager,locator_builder->tag_locator);
  // Write Tags
  svector_iterator_t tag_locator_iterator;
  svector_iterator_new(&tag_locator_iterator,locator_builder->tag_locator,SVECTOR_READ_ITERATOR,0);
  while (!svector_read_iterator_eoi(&tag_locator_iterator)) {
    locator_tag_t* const locator_tag = svector_iterator_get_element(&tag_locator_iterator,locator_tag_t);
    fm_write_mem(file_manager,locator_tag->tag,locator_tag->length+1);
    svector_read_iterator_next(&tag_locator_iterator);
  }
}
/*
 * Builder accessors
 */
GEM_INLINE uint64_t locator_builder_get_num_intervals(locator_builder_t* const locator_builder) {
  LOCATOR_BUILDER_CHECK(locator_builder);
  return svector_get_used(locator_builder->intervals);
}
GEM_INLINE uint64_t locator_builder_get_num_tags(locator_builder_t* const locator_builder) {
  LOCATOR_BUILDER_CHECK(locator_builder);
  return locator_builder->tag_id_generator;
}
GEM_INLINE locator_interval_t* locator_builder_get_interval(locator_builder_t* const locator_builder,const uint64_t interval_index) {
  LOCATOR_BUILDER_CHECK(locator_builder);
  return svector_get_element(locator_builder->intervals,interval_index,locator_interval_t);
}
GEM_INLINE locator_tag_t* locator_builder_get_tag(locator_builder_t* const locator_builder,const int64_t tag_id) {
  return svector_get_element(locator_builder->tag_locator,ABS(tag_id)-1,locator_tag_t);
}
GEM_INLINE locator_tag_t* locator_builder_locate_tag(locator_builder_t* const locator_builder,char* const tag) {
  LOCATOR_BUILDER_CHECK(locator_builder);
  GEM_CHECK_NULL(tag);
  return shash_get(locator_builder->tags_dictionary,tag,locator_tag_t);
}
GEM_INLINE locator_tag_t* locator_builder_interval_get_tag(
    locator_builder_t* const locator_builder,locator_interval_t* const interval) {
  LOCATOR_BUILDER_CHECK(locator_builder);
  return locator_builder_get_tag(locator_builder,interval->tag_id);
}
/*
 * Builder add intervals
 */
GEM_INLINE void locator_builder_open_interval(locator_builder_t* const locator_builder,const int64_t tag_id) {
  LOCATOR_BUILDER_CHECK(locator_builder);
  // Open a new interval => [begin_position,end_position)
  locator_builder->current_interval = svector_iterator_get_element(&(locator_builder->intervals_iterator),locator_interval_t);
  locator_builder->current_interval->tag_id = (tag_id!=0) ? tag_id : locator_builder->tag_id_generator;
  locator_builder->current_interval->begin_position = locator_builder->index_position;
}
GEM_INLINE void locator_builder_close_interval(
    locator_builder_t* const locator_builder,
    const uint64_t interval_text_length,const uint64_t interval_index_length,const locator_interval_type type) {
  LOCATOR_BUILDER_CHECK(locator_builder);
  // Close interval & write it => [begin_position,end_position)
  locator_builder->current_interval->type = type;
  locator_builder->current_interval->end_position = locator_builder->current_interval->begin_position+interval_index_length;
  locator_builder->current_interval->sequence_offset = locator_builder->sequence_offset;
  locator_builder->current_interval->sequence_length = interval_text_length;
  svector_write_iterator_next(&(locator_builder->intervals_iterator)); // Next
  // Update builder state
  locator_builder->index_position += interval_index_length;
  locator_builder->sequence_offset += interval_text_length;
}
GEM_INLINE void locator_builder_add_interval(
    locator_builder_t* const locator_builder,const int64_t tag_id,const uint64_t sequence_offset,
    const uint64_t interval_text_length,const uint64_t interval_index_length,const locator_interval_type type) {
  LOCATOR_BUILDER_CHECK(locator_builder);
  // Add interval
  locator_builder->current_interval = svector_iterator_get_element(&(locator_builder->intervals_iterator),locator_interval_t);
  locator_builder->current_interval->type = type;
  locator_builder->current_interval->begin_position = locator_builder->index_position;
  locator_builder->current_interval->end_position = locator_builder->index_position+interval_index_length;
  locator_builder->current_interval->sequence_offset = sequence_offset;
  locator_builder->current_interval->sequence_length = interval_text_length;
  locator_builder->current_interval->tag_id = (tag_id!=0) ? tag_id : locator_builder->tag_id_generator;
  svector_write_iterator_next(&(locator_builder->intervals_iterator));
  // Update builder state
  locator_builder->index_position += interval_index_length;
}
GEM_INLINE void locator_builder_add_rc_interval(
    locator_builder_t* const locator_builder,locator_interval_t* const locator_interval) {
  LOCATOR_BUILDER_CHECK(locator_builder);
  const uint64_t interval_length = locator_interval_get_index_length(locator_interval);
  const int64_t rc_tag_id = -(locator_interval->tag_id);
  locator_builder_add_interval(locator_builder,
      rc_tag_id,locator_interval->sequence_offset,
      interval_length,interval_length,locator_interval->type);
}
/*
 * Builder add intervals
 */
GEM_INLINE int64_t locator_builder_add_sequence(
    locator_builder_t* const locator_builder,char* const tag,const uint64_t tag_length) {
  LOCATOR_BUILDER_CHECK(locator_builder);
  // Reset the sequence offset
  locator_builder->sequence_offset = 0;
  // Check if its already in the dictionary of tags
  locator_tag_t* current_tag = locator_builder_locate_tag(locator_builder,tag);
  if (current_tag!=NULL) {
    // Return found TAG
    return current_tag->tag_id;
  } else {
    // Setup the current locator_tag & add it to the tag_buffer
    current_tag = svector_iterator_get_element(&(locator_builder->tag_locator_iterator),locator_tag_t);
    ++(locator_builder->tag_id_generator);
    current_tag->tag_id = locator_builder->tag_id_generator;
    current_tag->offset = locator_builder->tags_buffer_size;
    current_tag->tag = svector_insert_char_buffer(locator_builder->tags_buffer,tag,tag_length+1); // Insert EOS
    locator_builder->tags_buffer_size += tag_length+1;
    current_tag->length = tag_length;
    svector_write_iterator_next(&(locator_builder->tag_locator_iterator)); // Next
    // Add it to the dictionary
    shash_insert(locator_builder->tags_dictionary,current_tag->tag,current_tag->length,current_tag);
    // Return added TAG
    return current_tag->tag_id;
  }
}
GEM_INLINE void locator_builder_skip_index(locator_builder_t* const locator_builder,const uint64_t length) {
  LOCATOR_BUILDER_CHECK(locator_builder);
  locator_builder->index_position += length;
}
/*
 * Locating functions
 */
GEM_INLINE uint64_t locator_lookup_interval_index(const locator_t* const locator,const uint64_t index_position) {
  LOCATOR_CHECK(locator);
  // Binary search of the interval
  const locator_interval_t* const intervals = locator->intervals;
  uint64_t lo = 0;
  uint64_t hi = locator->num_intervals-1;
  do {
    const uint64_t half = (hi+lo+1)/2;
    if (intervals[half].begin_position > index_position) {
      hi = half-1;
    } else {
      lo = half;
    }
  } while (hi > lo);
  // Return Interval
  if (!(intervals[lo].begin_position <= index_position && index_position < intervals[lo].end_position)) {
    GEM_INTERNAL_CHECK(intervals[lo].begin_position <= index_position &&
        index_position < intervals[lo].end_position,"Locator-Interval Binary Search. Wrong Boundaries");
  }
  return lo;
}
GEM_INLINE locator_interval_t* locator_lookup_interval(const locator_t* const locator,const uint64_t index_position) {
  LOCATOR_CHECK(locator);
  return locator->intervals + locator_lookup_interval_index(locator,index_position);
}
/*
 * Map functions (High level mapping)
 */
// Direct Locator (Position-to-location mapping)
GEM_INLINE void locator_map(const locator_t* const locator,const uint64_t index_position,location_t* const location) {
  LOCATOR_CHECK(locator);
  GEM_CHECK_NULL(location);
  // Find interval
  const locator_interval_t* const intervals = locator_lookup_interval(locator,index_position);
  // Fill @location
  if (locator_interval_is_forward(intervals)) {
    location->position = intervals->sequence_offset + (index_position-intervals->begin_position);
    location->direction = Forward;
    location->tag = locator_interval_get_tag(locator,intervals);
  } else { // Reverse
    location->position = intervals->sequence_offset +
        (intervals->end_position-intervals->begin_position) - (index_position-intervals->begin_position);
      // TODO Memoizate +(intervals->end_position-intervals->begin_position) in the index
    location->direction = Reverse;
    location->tag = locator_interval_get_tag(locator,intervals);
  }
}
// Inverse Locator (Location-to-position mapping)
GEM_INLINE uint64_t inverse_locator_map(
    locator_t* const locator,const uint8_t* const tag,const strand_t strand,const int64_t text_position) {
  locator_interval_t* const intervals = locator->intervals;
  const uint64_t num_intervals = locator->num_intervals;
  uint64_t i;
//  if (strand==Forward) {
//    for (i=0;i<num_intervals;++i) {
//      if (intervals[i].tag_id > 0 && intervals[i].sequence_offset <= text_position &&
//          text_position < intervals[i].sequence_offset+intervals[i].sequence_length) {
//        if (gem_streq((char*)tag,locator_interval_get_tag(locator,intervals+i))) {
//          return intervals[i].begin_position + (text_position-intervals[i].sequence_offset);
//        }
//      }
//    }
//  } else {
//    for (i=0;i<num_intervals;++i) {
//      if (intervals[i].tag_id < 0 && intervals[i].sequence_offset <= text_position &&
//          text_position < intervals[i].sequence_offset+intervals[i].sequence_length) {
//        if (gem_streq((char*)tag,locator_interval_get_tag(locator,intervals+i))) {
//          return intervals[i].end_position - (text_position-intervals[i].sequence_offset);
//        }
//      }
//    }
//  }
  // TODO Improve by doing binary search
  // TODO Improve by doing stranded search
  for (i=0;i<num_intervals;++i) {
    if (intervals[i].sequence_offset <= text_position &&
        text_position < intervals[i].sequence_offset+intervals[i].sequence_length) {
      if (gem_streq((char*)tag,locator_interval_get_tag(locator,intervals+i))) {
        return intervals[i].begin_position + (text_position-intervals[i].sequence_offset);
      }
    }
  }
  gem_fatal_error(LOCATOR_INVERSE_NOT_FOUND);
  return 0;
}
/*
 * Display
 */
GEM_INLINE void locator_interval_print(FILE* const stream,const char* const tag,locator_interval_t* const interval) {
  // Print interval information
  switch (interval->type) {
    case locator_interval_regular:
      fprintf(stream,"Interval\tIdx=[%lu,%lu)\tText=%s:%c:[%lu,%lu)\n",
          interval->begin_position,interval->end_position,
          tag,locator_interval_is_forward(interval) ? '+' : '-',
          interval->sequence_offset,interval->sequence_offset+locator_interval_get_text_length(interval));
      break;
    case locator_interval_unknown:
      fprintf(stream,"Unknown\tIdx=[%lu,%lu)\tText=%s:%c:[%lu,%lu)\n",
          interval->begin_position,interval->end_position,
          tag,locator_interval_is_forward(interval) ? '+' : '-',
          interval->sequence_offset,interval->sequence_offset+locator_interval_get_text_length(interval));
      break;
    case locator_interval_variant:
      fprintf(stream,"Variant\tIdx=[%lu,%lu)\tText=%s:%c:[%lu,%lu)\n",
          interval->begin_position,interval->end_position,
          tag,locator_interval_is_forward(interval) ? '+' : '-',
          interval->sequence_offset,interval->sequence_offset+locator_interval_get_text_length(interval));
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
GEM_INLINE void locator_print_(FILE* const stream,const uint64_t num_intervals,const uint64_t num_tags,const uint64_t tags_buffer_size) {
  // Print locator info
  tab_fprintf(stream,"  => Intervals.Size %lu MB (%lu intervals x %luB per interval)\n",
      CONVERT_B_TO_MB(sizeof(locator_interval_t)*num_intervals),num_intervals,sizeof(locator_interval_t));
  tab_fprintf(stream,"  => Tags.Locator.Size %lu B (%lu locators x %luB per locator)\n",
      CONVERT_B_TO_MB(sizeof(locator_tag_t)*num_tags),num_tags,sizeof(locator_tag_t));
  tab_fprintf(stream,"  => Tags.Buffer.Size %lu B\n",tags_buffer_size);
}
GEM_INLINE void locator_print(FILE* const stream,const locator_t* const locator,const bool display_intervals) {
  GEM_CHECK_NULL(stream);
  LOCATOR_CHECK(locator);
  // Print locator info
  tab_fprintf(stream,"[GEM]>Locator\n");
  locator_print_(stream,locator->num_intervals,locator->num_tags,locator->tags_buffer_size);
  // Print intervals
  if (display_intervals) {
    tab_fprintf(stream,"  => Locator.intervals\n");
    uint64_t i;
    for (i=0;i<locator->num_intervals;++i) {
      // Get interval & tag
      locator_interval_t* const interval = locator->intervals + i;
      char* const tag = locator_interval_get_tag(locator,interval);
      // Print interval information
      tab_fprintf(stream,"       [%lu] ",i);
      locator_interval_print(stream,tag,interval);
    }
  }
  // Flush
  fflush(stream);
}
GEM_INLINE void locator_builder_print(FILE* const stream,locator_builder_t* const locator_builder,const bool display_intervals) {
  GEM_CHECK_NULL(stream);
  LOCATOR_BUILDER_CHECK(locator_builder);
  // Print locator info
  tab_fprintf(stream,"[GEM]>Locator\n");
  locator_print_(stream,
      svector_get_used(locator_builder->intervals),
      svector_get_used(locator_builder->tag_locator),
      locator_builder->tags_buffer_size);
  // Print intervals
  if (display_intervals) {
    tab_fprintf(stream,"  => Locator.intervals\n");
    uint64_t i = 0;
    svector_iterator_new(&(locator_builder->intervals_iterator),locator_builder->intervals,SVECTOR_READ_ITERATOR,0);
    while (!svector_read_iterator_eoi(&(locator_builder->intervals_iterator))) {
      // Get interval & tag
      locator_interval_t* const interval = svector_iterator_get_element(&(locator_builder->intervals_iterator),locator_interval_t);
      locator_tag_t* const tag = locator_builder_interval_get_tag(locator_builder,interval);
      // Print interval information
      tab_fprintf(stream,"       [%lu] ",i);
      locator_interval_print(stream,tag->tag,interval);
      // Next!
      ++i;
      svector_read_iterator_next(&(locator_builder->intervals_iterator));
    }
  }
  // Flush
  fflush(stream);
}
