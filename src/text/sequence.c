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

#include "text/dna_text.h"
#include "text/sequence.h"

/*
 * Constructor
 */
void sequence_init(
    sequence_t* const sequence,
    const bisulfite_read_t bisulfite_read_mode,
    mm_allocator_t* const mm_allocator) {
  // Sequence fields
  string_init(&sequence->tag,SEQUENCE_TAG_INITIAL_LENGTH,mm_allocator);
  string_init(&sequence->read,SEQUENCE_INITIAL_LENGTH,mm_allocator);
  string_init(&sequence->qualities,SEQUENCE_INITIAL_LENGTH,mm_allocator);
  sequence->end_info = single_end;
  // Extra fields
  string_init(&sequence->casava_tag,SEQUENCE_TAG_ATTRIBUTE_INITIAL_LENGTH,mm_allocator);
  string_init(&sequence->extra_tag,SEQUENCE_TAG_ATTRIBUTE_INITIAL_LENGTH,mm_allocator);
  // Bisulfite
  if (bisulfite_read_mode!=bisulfite_disabled) {
    sequence->bs_sequence_end = (bisulfite_read_mode==bisulfite_read_2) ? paired_end2 : paired_end1;
    string_init(&sequence->bs_original_sequence,SEQUENCE_BISULFITE_INITIAL_LENGTH,mm_allocator);
  }
}
void sequence_clear(sequence_t* const sequence) {
  string_clear(&sequence->tag);
  string_clear(&sequence->read);
  string_clear(&sequence->qualities);
  sequence->end_info = single_end;
  string_clear(&sequence->casava_tag);
  string_clear(&sequence->extra_tag);
}
void sequence_destroy(sequence_t* const sequence) {
  string_destroy(&sequence->tag);
  string_destroy(&sequence->read);
  string_destroy(&sequence->qualities);
  string_destroy(&sequence->casava_tag);
  string_destroy(&sequence->extra_tag);
}
/*
 * Accessors
 */
uint64_t sequence_get_length(sequence_t* const sequence) {
  return string_get_length(&sequence->read);
}
char* sequence_get_tag(sequence_t* const sequence) {
  return string_get_buffer(&sequence->tag);
}
void sequence_set_tag(
    sequence_t* sequence,
    char* const text,
    const uint64_t length) {
  string_set_buffer(&sequence->tag,text,length);
}
char* sequence_get_read(sequence_t* const sequence) {
  return string_get_buffer(&sequence->read);
}
void sequence_set_read(
    sequence_t* sequence,
    char* const text,
    const uint64_t length) {
  string_set_buffer(&sequence->read,text,length);
}
char* sequence_get_qualities(sequence_t* const sequence) {
  return string_get_buffer(&sequence->qualities);
}
void sequence_set_qualities(
    sequence_t* sequence,
    char* const text,
    const uint64_t length) {
  string_set_buffer(&sequence->qualities,text,length);
}
bool sequence_has_qualities(const sequence_t* const sequence) {
  return !string_is_null(&sequence->qualities);
}
bool sequence_has_casava_tag(const sequence_t* const sequence) {
  return !string_is_null(&sequence->casava_tag);
}
bool sequence_has_extra_tag(const sequence_t* const sequence) {
  return !string_is_null(&sequence->extra_tag);
}
sequence_end_t sequence_get_end_info(const sequence_t* const sequence) {
  return sequence->end_info;
}
/*
 * Utils
 */
bool sequence_equals(sequence_t* const sequence_a,sequence_t* const sequence_b) {
  return string_equals(&sequence_a->read,&sequence_b->read);
}
void sequence_generate_reverse(sequence_t* const sequence,sequence_t* const rev_sequence) {
  // Prepare rc_string (Read)
  const uint64_t seq_buffer_length = string_get_length(&sequence->read);
  string_resize(&rev_sequence->read,seq_buffer_length,false);
  string_clear(&rev_sequence->read);
  // Reverse Read
  int64_t pos;
  const char* const seq_buffer = string_get_buffer(&sequence->read);
  for (pos=seq_buffer_length-1;pos>=0;--pos) {
    string_append_char(&rev_sequence->read,seq_buffer[pos]);
  }
  string_append_eos(&rev_sequence->read);
  // Reverse Qualities
  if (sequence_has_qualities(sequence)) {
    string_copy_reverse(&rev_sequence->qualities,&sequence->qualities);
  }
}
void sequence_generate_reverse_complement(sequence_t* const sequence,sequence_t* const rc_sequence) {
  // Prepare rc_string (Read)
  const uint64_t seq_buffer_length = string_get_length(&sequence->read);
  string_resize(&rc_sequence->read,seq_buffer_length,false);
  string_clear(&rc_sequence->read);
  // Reverse-Complement Read
  int64_t pos;
  const char* const seq_buffer = string_get_buffer(&sequence->read);
  for (pos=seq_buffer_length-1;pos>=0;--pos) {
    string_append_char(&rc_sequence->read,dna_complement(seq_buffer[pos]));
  }
  string_append_eos(&rc_sequence->read);
  // Reverse Qualities
  if (sequence_has_qualities(sequence)) {
    string_copy_reverse(&rc_sequence->qualities,&sequence->qualities);
  }
}
/*
 * Display
 */
void sequence_print(FILE* const stream,sequence_t* const sequence) {
  // tab_fprintf(stream,"[GEM]>Sequence\n");
  if (sequence_has_qualities(sequence)) {
    fprintf(stream,"@%s\n",sequence_get_tag(sequence));
    fprintf(stream,"%s+\n",sequence_get_read(sequence));
    fprintf(stream,"%s",sequence_get_qualities(sequence));
  } else {
    fprintf(stream,">%s\n",sequence_get_tag(sequence));
    fprintf(stream,"%s",sequence_get_read(sequence));
  }
}

