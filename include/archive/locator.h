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
 *   Locator data structure enables translation between multiple sequences
 *   coordinates to one single sequence coordinates an viceversa (i.e. GEM3
 *   index pastes all individual sequences into a monolithic sequence-index)
 */

#ifndef LOCATOR_H_
#define LOCATOR_H_

#include "utils/essentials.h"
#include "utils/segmented_vector.h"
#include "text/dna_text.h"

/*
 * Locator Model & Version
 */
#define LOCATOR_MODEL_NO  3006ull

/*
 * Locator
 */
typedef struct {
  uint64_t tag_id;
  uint64_t length;
  uint64_t offset;
} locator_tag_t;
typedef enum {
  locator_interval_chromosomal_assembly = 0,  // Standard chromosome
  locator_interval_unlocalized_contigs = 1,   // Chromosome known but location unknown
  locator_interval_unplaced_contig = 2,       // Chromosome unknown
  locator_interval_alt_contig  = 3,           // Long clustered variations
  locator_interval_uncalled = 4,              // Uncalled/Unknown region
  locator_interval_variant = UINT64_MAX,      // Variant
} locator_interval_type;
typedef struct {
  // Interval type
  locator_interval_type type; // Sequence type (Regular, Variant, ...)
  strand_t strand;            // DNA Strand
  bs_strand_t bs_strand;      // Bisulfite strand
  // Positions over INDEX (One single sequence composed by pasting all reference sequences)
  //    Interval  = [begin_position,end_position) (Zero Based)
  //   |Interval| = end_position-begin_position
  uint64_t begin_position;    // Global bottom location (over all text)
  uint64_t end_position;      // Global top location (over all text)
  // Positions over RL-INDEX (One single sequence composed by pasting all reference sequences & RL)
  uint64_t rl_begin_position; // Global bottom location (over all rl-text)
  uint64_t rl_end_position;   // Global top location (over all rl-text)
  // Position over TEXT => Eg: 'Chunk of Chr1 starting from 30'
  uint64_t sequence_offset;   // Offset relative to the sequence/chromosome (wrt all sequences/chromosome)
  uint64_t sequence_length;   // Sequence total length
  // TAG ID
  int64_t tag_id;             // ID (Sequence/Chromosome name ID)
} locator_interval_t;
typedef struct {
  /* Intervals */
  locator_interval_t* intervals;      // Intervals
  uint64_t* intervals_lookahead;      // LUT Look-ahead search
  uint64_t* intervals_begin_position; // LUT Intervals[i].begin_position
  uint64_t num_intervals;
  /* Tag locator */
  locator_tag_t* tag_locator;       // Tag Locator
  uint64_t num_tags;
  /* Tags Buffer*/
  char* tags_buffer;                // Tag Buffer
  uint64_t tags_buffer_size;
  /* MM */
  mm_t* mm;
} locator_t;

/*
 * Locator query (Sequence Location: Container to return results from query locations)
 */
typedef struct {
  char* tag;
  strand_t strand;
  bs_strand_t bs_strand;
  int64_t position;
} location_t;

/*
 * Loader/Setup
 */
locator_t* locator_read_mem(mm_t* const memory_manager);
void locator_setup_inverse_locator(locator_t* const locator);
void locator_delete(locator_t* const locator);

/*
 * Locator Accessors
 */
uint64_t locator_get_size(locator_t* const locator);
locator_interval_t* locator_get_interval(
    const locator_t* const locator,
    const uint64_t interval_index);

/*
 * Interval Accessors
 */
char* locator_interval_get_tag(
    const locator_t* const locator,
    const locator_interval_t* const interval);
uint64_t locator_interval_get_index_length(const locator_interval_t* const interval);
uint64_t locator_interval_get_text_length(const locator_interval_t* const interval);

/*
 * Text-Locating functions
 */
locator_interval_t* locator_lookup_interval(
    const locator_t* const locator,
    const uint64_t index_position);

/*
 * RL-Locating functions
 */
locator_interval_t* locator_lookup_rl_interval(
    const locator_t* const locator,
    const uint64_t rl_index_position);

/*
 * Map functions (High level mapping)
 */
// [Direct Locator] (Position-to-location mapping)
void locator_map(
    const locator_t* const locator,
    const uint64_t index_position,
    location_t* const location);
// [Inverse Locator] (Location-to-position mapping)
locator_interval_t* locator_inverse_map(
    locator_t* const locator,
    const uint8_t* const tag,
    const strand_t strand,
    const bs_strand_t bs_strand,
    const uint64_t text_position);
uint64_t locator_inverse_map_position(
    locator_t* const locator,
    const uint8_t* const tag,
    const strand_t strand,
    const bs_strand_t bs_strand,
    const uint64_t text_position);


/*
 * Display
 */
void locator_interval_print(
    FILE* const stream,
    locator_interval_t* const interval,
    const char* const interval_tag);
void locator_print_summary(
    FILE* const stream,
    const uint64_t num_intervals,
    const uint64_t num_tags,
    const uint64_t tags_buffer_size);
void locator_print(
    FILE* const stream,
    const locator_t* const locator,
    const bool display_intervals);

#endif /* LOCATOR_H_ */
