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

#include "matches/align/match_alignment_region.h"
#include "matches/matches_cigar.h"
#include "align/alignment.h"

/*
 * Accessors
 */
void match_alignment_region_init(
    match_alignment_region_t* const match_alignment_region,
    const match_alignment_region_type type,
    const uint64_t error,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length,
    const uint64_t key_begin,
    const uint64_t key_end,
    const uint64_t text_begin,
    const uint64_t text_end) {
  match_alignment_region_set_type(match_alignment_region,type);
  match_alignment_region->error = error;
  match_alignment_region->cigar_buffer_offset = cigar_buffer_offset;
  match_alignment_region->cigar_length = cigar_length;
  match_alignment_region->key_begin = key_begin;
  match_alignment_region->key_end = key_end;
  match_alignment_region->text_begin = text_begin;
  match_alignment_region->text_end = text_end;
}
match_alignment_region_type match_alignment_region_get_type(
    const match_alignment_region_t* const match_alignment_region) {
  return match_alignment_region->type;
}
void match_alignment_region_set_type(
    match_alignment_region_t* const match_alignment_region,
    const match_alignment_region_type type) {
  if (type==0) {
    match_alignment_region->type = match_alignment_region_exact;
  } else {
    match_alignment_region->type = match_alignment_region_approximate;
  }
}
uint64_t match_alignment_region_get_error(
    const match_alignment_region_t* const match_alignment_region) {
  return match_alignment_region->error;
}
void match_alignment_region_set_error(
    match_alignment_region_t* const match_alignment_region,
    const uint64_t error) {
  match_alignment_region->error = error;
}
uint64_t match_alignment_region_get_cigar_buffer_offset(
    const match_alignment_region_t* const match_alignment_region) {
  return match_alignment_region->cigar_buffer_offset;
}
uint64_t match_alignment_region_get_cigar_length(
    const match_alignment_region_t* const match_alignment_region) {
  return match_alignment_region->cigar_length;
}
uint64_t match_alignment_region_get_key_begin(
    const match_alignment_region_t* const match_alignment_region) {
  return match_alignment_region->key_begin;
}
uint64_t match_alignment_region_get_key_end(
    const match_alignment_region_t* const match_alignment_region) {
  return match_alignment_region->key_end;
}
uint64_t match_alignment_region_get_text_begin(
    const match_alignment_region_t* const match_alignment_region) {
  return match_alignment_region->text_begin;
}
uint64_t match_alignment_region_get_text_end(
    const match_alignment_region_t* const match_alignment_region) {
  return match_alignment_region->text_end;
}
void match_alignment_region_set_cigar_buffer_offset(
    match_alignment_region_t* const match_alignment_region,
    const uint64_t cigar_buffer_offset) {
  match_alignment_region->cigar_buffer_offset = cigar_buffer_offset;
}
void match_alignment_region_set_cigar_length(
    match_alignment_region_t* const match_alignment_region,
    const uint64_t cigar_length) {
  match_alignment_region->cigar_length = cigar_length;
}
void match_alignment_region_set_key_begin(
    match_alignment_region_t* const match_alignment_region,
    const uint64_t key_begin) {
  match_alignment_region->key_begin = key_begin;
}
void match_alignment_region_set_key_end(
    match_alignment_region_t* const match_alignment_region,
    const uint64_t key_end) {
  match_alignment_region->key_end = key_end;
}
void match_alignment_region_set_text_begin(
    match_alignment_region_t* const match_alignment_region,
    const uint64_t text_begin) {
  match_alignment_region->text_begin = text_begin;
}
void match_alignment_region_set_text_end(
    match_alignment_region_t* const match_alignment_region,
    const uint64_t text_end) {
  match_alignment_region->text_end = text_end;
}
/*
 * Key/Text Operators
 */
int match_alignment_region_text_coverage(
    match_alignment_region_t* const match_alignment_region) {
  return (int)match_alignment_region->text_end - (int)match_alignment_region->text_begin;
}
int match_alignment_region_text_distance(
    match_alignment_region_t* const match_alignment_region_a,
    match_alignment_region_t* const match_alignment_region_b) {
  return (int)match_alignment_region_b->text_begin - (int)match_alignment_region_a->text_end;
}
bool match_alignment_region_text_overlap(
    match_alignment_region_t* const match_alignment_region_a,
    match_alignment_region_t* const match_alignment_region_b) {
  return !(match_alignment_region_a->text_end <= match_alignment_region_b->text_begin ||
           match_alignment_region_b->text_end <= match_alignment_region_a->text_begin);
}
/*
 * Compare
 */
int64_t match_alignment_region_cmp_key_position(
    const match_alignment_region_t* const a,
    const match_alignment_region_t* const b) {
  return (int64_t)a->key_begin - (int64_t)b->key_begin;
}
int64_t match_alignment_region_cmp_text_position(
    const match_alignment_region_t* const a,
    const match_alignment_region_t* const b) {
  return (int64_t)a->text_begin - (int64_t)b->text_begin;
}
/*
 * Check
 */
void match_alignment_region_check(
    match_alignment_region_t* const match_alignment_region,
    uint8_t* const key,
    uint8_t* const text,
    vector_t* const cigar_vector) {
  // Parameters
  const uint64_t key_region_length = match_alignment_region->key_end - match_alignment_region->key_begin;
  uint8_t* const key_region = key + match_alignment_region->key_begin;
  const uint64_t text_region_length = match_alignment_region->text_end - match_alignment_region->text_begin;
  uint8_t* const text_region = text + match_alignment_region->text_begin;
  // Check
  bool alignment_correct;
  if (match_alignment_region->error == 0) {
    // Check all matching characters
    uint64_t i;
    for (i=0;i<text_region_length;++i) {
      if (key[i] != text[i]) {
        alignment_correct = false;
        break;
      }
    }
    alignment_correct = true;
  } else {
    alignment_correct = alignment_check(
        stderr,key_region,key_region_length,
        text_region,text_region_length,
        cigar_vector,match_alignment_region->cigar_buffer_offset,
        match_alignment_region->cigar_length,false,false);
  }
  if (!alignment_correct) {
    gem_fatal_error_msg("Match alignment-region check. Invalid alignment (region)");
  }
}
/*
 * Display
 */
void match_alignment_region_print(
    FILE* const stream,
    match_alignment_region_t* const match_alignment_region,
    const uint64_t match_alignment_region_id,
    matches_t* const matches) {
  // Print alignment-region
  switch (match_alignment_region_get_type(match_alignment_region)) {
    case match_alignment_region_exact:
      tab_fprintf(stream,"    %"PRIu64"[exact]\t",match_alignment_region_id);
      break;
    case match_alignment_region_approximate:
      tab_fprintf(stream,"    %"PRIu64"[approximate]\t",match_alignment_region_id);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  tab_fprintf(stream,"-> [%"PRIu64",%"PRIu64") ~> [+%"PRIu64",+%"PRIu64")",
      match_alignment_region->key_begin,match_alignment_region->key_end,
      match_alignment_region->text_begin,match_alignment_region->text_end);
  // Print CIGAR
  if (matches!=NULL &&
      match_alignment_region_get_type(match_alignment_region)==match_alignment_region_approximate &&
      match_alignment_region->cigar_length>0) {
    tab_fprintf(stream,"\tCIGAR=");
    match_cigar_print(stream,matches->cigar_vector,
        match_alignment_region->cigar_buffer_offset,
        match_alignment_region->cigar_length);
  }
  tab_fprintf(stream,"\n");
}
void match_alignment_region_print_pretty(
    FILE* const stream,
    match_alignment_region_t* const match_alignment_region,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text) {
  // Compute offset
  uint64_t i, offset_text = 0, offset_key = 0;
  if (match_alignment_region->key_begin > match_alignment_region->text_begin) {
    offset_key = match_alignment_region->key_begin-match_alignment_region->text_begin;
  } else {
    offset_text = match_alignment_region->text_begin-match_alignment_region->key_begin;
    for (i=0;i<offset_text;++i) fprintf(stream," ");
  }
  // Print Key
  for (i=offset_key;i<key_length;++i) {
    if (match_alignment_region->key_begin <= i && i < match_alignment_region->key_end) {
      if (text[offset_text+i]==key[i]) {
        fprintf(stream,"%c",dna_decode(key[i]));
      } else {
        fprintf(stream,"%c",tolower(dna_decode(key[i])));
      }
    } else {
      fprintf(stream,"-");
    }
  }
  fprintf(stream,"\n");
}
