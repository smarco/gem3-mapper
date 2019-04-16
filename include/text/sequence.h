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
 * DESCRIPTION: Simple data structure to store genomic reads
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include "utils/essentials.h"

/*
 * Checkers
 */
#define SEQUENCE_QUALITY_IS_VALID(character) (33 <= (character))

/*
 * Constants
 */
#define SEQUENCE_TAG_INITIAL_LENGTH 80
#define SEQUENCE_TAG_ATTRIBUTE_INITIAL_LENGTH 80
#define SEQUENCE_INITIAL_LENGTH 200
#define SEQUENCE_BISULFITE_INITIAL_LENGTH 200

/*
 * Bisulfite Sequence
 */
typedef enum {
  bisulfite_disabled,
  bisulfite_inferred_C2T_G2A, // First read C2T converted, second (if present) G2A
  bisulfite_inferred_G2A_C2T, // First read G2A converted, second (if present) C2T
  bisulfite_C2T,  // Perform C2T conversion
  bisulfite_G2A,  // Perform G2A conversion
	bisulfite_non_stranded,  // Perform both conversions and merge
} bisulfite_read_t;

typedef enum {
    no_conversion,
    C2T_conversion,
    G2A_conversion
} bisulfite_conversion_t;

/*
 * Sequence
 */
typedef enum { single_end=0, paired_end1=1, paired_end2=2 } sequence_end_t;
typedef struct {
  /* Sequence */
  string_t tag;                    // Sequence Tag
  string_t read;                   // Sequence Read
  string_t qualities;              // Sequence Qualities
  /* Attributes */
  bool has_qualities;
  sequence_end_t end_info;
  string_t casava_tag;
  string_t extra_tag;
} sequence_t;

/*
 * Constructor
 */
void sequence_init(
    sequence_t* const sequence,
    const bisulfite_read_t bisulfite_read_mode,
    mm_allocator_t* const mm_allocator);
void sequence_clear(sequence_t* const sequence);
void sequence_destroy(sequence_t* const sequence);

/*
 * Accessors
 */
uint64_t sequence_get_length(sequence_t* const sequence);

char* sequence_get_tag(sequence_t* const sequence);
void sequence_set_tag(
    sequence_t* const sequence,
    char* const text,
    const uint64_t length);

char* sequence_get_read(sequence_t* const sequence);
void sequence_set_read(
    sequence_t* const sequence,
    char* const text,
    const uint64_t length);

char* sequence_get_qualities(sequence_t* const sequence);
void sequence_set_qualities(
    sequence_t* const sequence,
    char* const text,
    const uint64_t length);

bool sequence_has_qualities(const sequence_t* const sequence);
bool sequence_has_casava_tag(const sequence_t* const sequence);
bool sequence_has_extra_tag(const sequence_t* const sequence);

sequence_end_t sequence_get_end_info(const sequence_t* const sequence);

/*
 * Utils
 */
bool sequence_equals(sequence_t* const sequence_a,sequence_t* const sequence_b);
void sequence_generate_reverse(sequence_t* const sequence,sequence_t* const rev_sequence);
void sequence_generate_reverse_complement(sequence_t* const sequence,sequence_t* const rc_sequence);

/*
 * Display
 */
void sequence_print(FILE* const stream,sequence_t* const sequence);

#endif /* SEQUENCE_H_ */
