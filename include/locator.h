/*
 * PROJECT: GEMMapper
 * FILE: locator.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides functionality as to recover sequence location from an archive of sequences
 */

#ifndef LOCATOR_H_
#define LOCATOR_H_

#include "essentials.h"
#include "segmented_vector.h"

/*
 * Checkers
 */
#define LOCATOR_CHECK(locator) \
  GEM_CHECK_NULL(locator)
#define LOCATOR_BUILDER_CHECK(locator_builder) \
  GEM_CHECK_NULL(locator_builder)
#define LOCATOR_INTERVAL_CHECK(interval) \
  GEM_CHECK_NULL(interval)

// Tag-Locator
typedef struct {
  uint64_t tag_id;
  char* tag;
  uint64_t length;
  uint64_t offset;
} locator_tag_t;
// Interval-Locator (Stored into the index)
typedef enum {
  locator_interval_regular = 0,
  locator_interval_unknown = 1,
  locator_interval_variant = UINT64_MAX
} locator_interval_type;
typedef struct {
  /*
   * Interval type
   */
  locator_interval_type type;
  /*
   * Positions over INDEX (One single sequence composed by pasting all reference sequences)
   *   Interval  = [begin_position,end_position) (Zero Based)
   *  |Interval| = end_position-begin_position
   */
  uint64_t begin_position;  // Global bottom location (over all indexed text)
  uint64_t end_position;    // Global top location (over all indexed text)
  /*
   * Position over TEXT => Eg: 'Chunk of Chr1 starting from 30'
   */
  uint64_t sequence_offset; // Offset relative to the sequence/chromosome (wrt all sequences/chromosome)
  uint64_t sequence_length; // Sequence total length
  /*
   * TAG ID/info (Sequence/Chromosome name)
   *   - Additionally, the strand is encoded as the sign of the tag_id
   *       tag_id if offset by 1. Range (-inf,1] U [1,+inf)
   *       (tag_id > 0) ? "Forward" : "Reverse"
   *
   */
  int64_t tag_id;
} locator_interval_t;
typedef struct {
  // Forward
  uint64_t forward_low_interval_idx;
  uint64_t forward_high_interval_idx;
  // Reverse
  uint64_t reverse_low_interval_idx;
  uint64_t reverse_high_interval_idx;
  // Tag locator
  locator_tag_t* tag_locator;
} locator_inverse_t;
typedef struct {
  /* Intervals */
  locator_interval_t* intervals;    // Intervals
  uint64_t num_intervals;
  /* Tag locator */
  shash_t* tags_dictionary;         // Tag's Dictionary use for inverse locator (locator_inverse_t*)
  locator_tag_t* tag_locator;       // Tag Locator // FIXME: Doesn't make sense to have a pointer and an ID
  uint64_t num_tags;
  /* Tags Buffer */
  char* tags_buffer;                // Tag Buffer
  uint64_t tags_buffer_size;
  /* MM */
  mm_t* mm;
} locator_t;
// Sequence Location (Container to return results from query locations)
typedef struct {
  int64_t position;
  strand_t direction;
  char* tag;
} location_t;
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
  uint64_t tag_id_generator;
  uint64_t index_position;                 // Current position in the index
  uint64_t sequence_offset;                // Current position in the text (Offset wrt the current sequence)
} locator_builder_t;

/*
 * Loader/Setup
 */
GEM_INLINE locator_t* locator_read(fm_t* const file_manager);
GEM_INLINE locator_t* locator_read_mem(mm_t* const memory_manager);
GEM_INLINE void locator_setup_inverse_locator(locator_t* const locator);
GEM_INLINE void locator_delete(locator_t* const locator);

/*
 * Accessors
 */
// [Locator]
GEM_INLINE uint64_t locator_get_size(locator_t* const locator);
GEM_INLINE locator_interval_t* locator_get_interval(const locator_t* const locator,const uint64_t interval_index);
// [Locator Interval]
GEM_INLINE char* locator_interval_get_tag(const locator_t* const locator,const locator_interval_t* const interval);
GEM_INLINE uint64_t locator_interval_get_index_length(const locator_interval_t* const interval);
GEM_INLINE uint64_t locator_interval_get_text_length(const locator_interval_t* const interval);
GEM_INLINE bool locator_interval_is_forward(const locator_interval_t* const interval);

/*
 * Builder
 */
GEM_INLINE locator_builder_t* locator_builder_new(mm_slab_t* const mm_slab);
GEM_INLINE void locator_builder_delete(locator_builder_t* const locator_builder);
GEM_INLINE void locator_builder_write(fm_t* const file_manager,locator_builder_t* const locator_builder);
// [Builder accessors]
GEM_INLINE uint64_t locator_builder_get_num_intervals(locator_builder_t* const locator_builder);
GEM_INLINE uint64_t locator_builder_get_num_tags(locator_builder_t* const locator_builder);
GEM_INLINE locator_interval_t* locator_builder_get_interval(locator_builder_t* const locator_builder,const uint64_t interval_index);
GEM_INLINE locator_tag_t* locator_builder_get_tag(locator_builder_t* const locator_builder,const int64_t tag_id);
GEM_INLINE locator_tag_t* locator_builder_locate_tag(locator_builder_t* const locator_builder,char* const tag);
GEM_INLINE locator_tag_t* locator_builder_interval_get_tag(
    locator_builder_t* const locator_builder,locator_interval_t* const interval);
// [Builder intervals]
GEM_INLINE void locator_builder_add_interval(
    locator_builder_t* const locator_builder,const int64_t tag_id,const uint64_t sequence_offset,
    const uint64_t interval_text_length,const uint64_t interval_index_length,const locator_interval_type type);
GEM_INLINE void locator_builder_add_rc_interval(
    locator_builder_t* const locator_builder,locator_interval_t* const locator_interval);
GEM_INLINE int64_t locator_builder_add_sequence(
    locator_builder_t* const locator_builder,char* const tag,const uint64_t tag_length);
GEM_INLINE void locator_builder_open_interval(locator_builder_t* const locator_builder,const int64_t tag_id);
GEM_INLINE void locator_builder_close_interval(
    locator_builder_t* const locator_builder,
    const uint64_t interval_text_length,const uint64_t interval_index_length,const locator_interval_type type);
// [Builder skip text/index]
GEM_INLINE void locator_builder_skip_index(locator_builder_t* const locator_builder,const uint64_t length);

/*
 * Locating functions
 */
GEM_INLINE uint64_t locator_lookup_interval_index(const locator_t* const locator,const uint64_t index_position);
GEM_INLINE locator_interval_t* locator_lookup_interval(const locator_t* const locator,const uint64_t index_position);

/*
 * Map functions (High level mapping)
 */
// [Direct Locator] (Position-to-location mapping)
GEM_INLINE void locator_map(const locator_t* const locator,const uint64_t index_position,location_t* const location);
// [Inverse Locator] (Location-to-position mapping)
GEM_INLINE uint64_t inverse_locator_map(
    locator_t* const locator,const uint8_t* const tag,const strand_t strand,const int64_t text_position);

/*
 * Display
 */
GEM_INLINE void locator_print(FILE* const stream,const locator_t* const locator,const bool display_intervals);
GEM_INLINE void locator_builder_print(FILE* const stream,locator_builder_t* const locator_builder,const bool display_intervals);

/*
 * Error Messages
 */
#define GEM_ERROR_LOCATOR_INTERVAL_INDEX_OOB "Locator. Requested locator-interval index (%lu) out-of-bounds [%lu,%lu)]"
#define GEM_ERROR_LOCATOR_INVERSE_NOT_FOUND "Locator. Inverse location not found"

#endif /* LOCATOR_H_ */
