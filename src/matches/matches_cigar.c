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

#include "matches/matches_cigar.h"
#include "matches/classify/matches_classify.h"

/*
 * Error Messages
 */
#define GEM_ERROR_MATCHES_CIGAR_ZERO_LENGTH "Matches. CIGAR length cannot be zero"

/*
 * CIGAR Buffer Handling
 */
void matches_cigar_buffer_add_cigar_element__attr(
    cigar_element_t** const cigar_buffer_sentinel,
    const cigar_t cigar_element_type,
    const uint64_t element_length,
    const cigar_attr_t attributes) {
  if ((*cigar_buffer_sentinel)->type == cigar_element_type) {
    (*cigar_buffer_sentinel)->length += element_length;
    (*cigar_buffer_sentinel)->attributes = cigar_attr_none;
  } else {
    if ((*cigar_buffer_sentinel)->type != cigar_null) ++(*cigar_buffer_sentinel);
    (*cigar_buffer_sentinel)->type = cigar_element_type;
    (*cigar_buffer_sentinel)->length = element_length;
    (*cigar_buffer_sentinel)->attributes = attributes;
  }
}
void matches_cigar_buffer_add_cigar_element(
    cigar_element_t** const cigar_buffer_sentinel,
    const cigar_t cigar_element_type,
    const uint64_t element_length) {
  matches_cigar_buffer_add_cigar_element__attr(cigar_buffer_sentinel,
      cigar_element_type,element_length,cigar_attr_none);
}
void matches_cigar_buffer_add_mismatch(
    cigar_element_t** const cigar_buffer_sentinel,
    const uint8_t mismatch) {
  if ((*cigar_buffer_sentinel)->type!=cigar_null) ++(*cigar_buffer_sentinel);
  (*cigar_buffer_sentinel)->type = cigar_mismatch;
  (*cigar_buffer_sentinel)->mismatch = mismatch; // Mismatch
  (*cigar_buffer_sentinel)->attributes = cigar_attr_none;
}
/*
 * CIGAR Vector Handling
 */
void matches_cigar_vector_append_insertion(
    vector_t* const cigar_vector,
    uint64_t* const current_cigar_length,
    const uint64_t indel_length,
    const cigar_attr_t attributes) {
  if (*current_cigar_length > 0) {
    cigar_element_t* cigar_element = vector_get_last_elm(cigar_vector,cigar_element_t);
    if (cigar_element->type==cigar_ins) {
      cigar_element->length += indel_length;
      cigar_element->attributes = cigar_attr_none;
      return;
    }
  }
  // Append a new one
  vector_reserve_additional(cigar_vector,1); // Reserve
  cigar_element_t* const cigar_element = vector_get_free_elm(cigar_vector,cigar_element_t);// Add CIGAR element
  cigar_element->type = cigar_ins;
  cigar_element->length = indel_length;
  cigar_element->attributes = attributes;
  vector_inc_used(cigar_vector); // Increment used
  *current_cigar_length += 1;
}
void matches_cigar_vector_append_deletion(
    vector_t* const cigar_vector,
    uint64_t* const current_cigar_length,
    const uint64_t indel_length,
    const cigar_attr_t attributes) {
  if (*current_cigar_length > 0) {
    cigar_element_t* cigar_element = vector_get_last_elm(cigar_vector,cigar_element_t);
    if (cigar_element->type==cigar_del) {
      cigar_element->length += indel_length;
      cigar_element->attributes = cigar_attr_none;
      return;
    }
  }
  // Append a new one
  vector_reserve_additional(cigar_vector,1); // Reserve
  cigar_element_t* const cigar_element = vector_get_free_elm(cigar_vector,cigar_element_t);// Add CIGAR element
  cigar_element->type = cigar_del;
  cigar_element->length = indel_length;
  cigar_element->attributes = attributes;
  vector_inc_used(cigar_vector); // Increment used
  *current_cigar_length += 1;
}
void matches_cigar_vector_append_match(
    vector_t* const cigar_vector,
    uint64_t* const current_cigar_length,
    const uint64_t match_length,
    const cigar_attr_t attributes) {
  // Check previous cigar-element (for merging)
  if (*current_cigar_length > 0) {
    cigar_element_t* cigar_element = vector_get_last_elm(cigar_vector,cigar_element_t);
    if (cigar_element->type==cigar_match) {
      cigar_element->length += match_length;
      cigar_element->attributes = cigar_attr_none;
      return;
    }
  }
  // Append a new one
  vector_reserve_additional(cigar_vector,1); // Reserve
  cigar_element_t* const cigar_element = vector_get_free_elm(cigar_vector,cigar_element_t);// Add CIGAR element
  cigar_element->type = cigar_match;
  cigar_element->length = match_length;
  cigar_element->attributes = attributes;
  vector_inc_used(cigar_vector); // Increment used
  *current_cigar_length += 1;
}
void matches_cigar_vector_append_mismatch(
    vector_t* const cigar_vector,
    uint64_t* const current_cigar_length,
    const uint8_t mismatch,
    const cigar_attr_t attributes) {
  vector_reserve_additional(cigar_vector,1); // Reserve
  cigar_element_t* const cigar_element = vector_get_free_elm(cigar_vector,cigar_element_t);// Add CIGAR element
  cigar_element->type = cigar_mismatch;
  cigar_element->mismatch = mismatch;
  cigar_element->attributes = attributes;
  vector_inc_used(cigar_vector); // Increment used
  *current_cigar_length += 1;
}
void matches_cigar_vector_append_cigar_element(
    vector_t* const cigar_vector,
    uint64_t* const cigar_length,
    cigar_element_t* const cigar_element) {
  switch (cigar_element->type) {
    case cigar_match:
      matches_cigar_vector_append_match(cigar_vector,cigar_length,cigar_element->length,cigar_element->attributes);
      break;
    case cigar_mismatch:
      matches_cigar_vector_append_mismatch(cigar_vector,cigar_length,cigar_element->mismatch,cigar_element->attributes);
      break;
    case cigar_ins:
      matches_cigar_vector_append_insertion(cigar_vector,cigar_length,cigar_element->length,cigar_element->attributes);
      break;
    case cigar_del:
      matches_cigar_vector_append_deletion(cigar_vector,cigar_length,cigar_element->length,cigar_element->attributes);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
void matches_cigar_vector_insert_cigar_element(
    vector_t* const cigar_vector,
    const int64_t position,
    cigar_element_t* const cigar_element) {
  // Reserve element & shift elements
  vector_reserve_additional(cigar_vector,1);
  cigar_element_t* const cigar_element_buffer = vector_get_mem(cigar_vector,cigar_element_t);
  const uint64_t cigar_vector_used = vector_get_used(cigar_vector);
  int64_t i;
  for (i=cigar_vector_used-1;i>=position;--i) {
    cigar_element_buffer[i+1] = cigar_element_buffer[i];
  }
  // Insert new element
  cigar_element_buffer[position] = *cigar_element;
  // Set used
  vector_set_used(cigar_vector,cigar_vector_used+1);
}
/*
 * CIGAR Vector Utils
 */
void matches_cigar_reverse(
    vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length) {
  // Reverse CIGAR
  cigar_element_t* const cigar_buffer = vector_get_elm(cigar_vector,cigar_buffer_offset,cigar_element_t);
  const uint64_t middle_point = cigar_length/2;
  uint64_t i;
  for (i=0;i<middle_point;++i) {
    cigar_element_t* const origin = cigar_buffer + i;
    cigar_element_t* const flipped = cigar_buffer + (cigar_length-1-i);
    SWAP(*origin,*flipped);
    if (origin->type == cigar_mismatch) origin->mismatch = dna_encoded_complement(origin->mismatch);
    if (flipped->type == cigar_mismatch) flipped->mismatch = dna_encoded_complement(flipped->mismatch);
  }
  if (cigar_length%2) {
    cigar_element_t* const middle = cigar_buffer + middle_point;
    if (middle->type == cigar_mismatch) middle->mismatch = dna_encoded_complement(middle->mismatch);
  }
}
void matches_cigar_reverse_colorspace(
    vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length) {
  // Reverse CIGAR
  cigar_element_t* const cigar_buffer = vector_get_elm(cigar_vector,cigar_buffer_offset,cigar_element_t);
  const uint64_t middle_point = cigar_length/2;
  uint64_t i;
  for (i=0;i<middle_point;++i) {
    cigar_element_t* const origin = cigar_buffer + i;
    cigar_element_t* const flipped = cigar_buffer + (cigar_length-1-i);
    SWAP(*origin,*flipped);
  }
}
uint64_t matches_cigar_compute_event_distance(
    vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length) {
  // Sum up all cigar elements
  const cigar_element_t* const cigar_buffer = vector_get_elm(cigar_vector,cigar_buffer_offset,cigar_element_t);
  uint64_t i, event_distance = 0;
  for (i=0;i<cigar_length;++i) {
    switch (cigar_buffer[i].type) {
      case cigar_match:
        break;
      case cigar_mismatch:
      case cigar_ins:
      case cigar_del:
        ++event_distance;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  return event_distance;
}
uint64_t matches_cigar_compute_edit_distance(
    vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length) {
  // Sum up all cigar elements
  const cigar_element_t* const cigar_buffer = vector_get_elm(cigar_vector,cigar_buffer_offset,cigar_element_t);
  uint64_t i, edit_distance = 0;
  for (i=0;i<cigar_length;++i) {
    switch (cigar_buffer[i].type) {
      case cigar_match: break;
      case cigar_mismatch:
        ++edit_distance;
        break;
      case cigar_ins:
      case cigar_del:
        edit_distance += cigar_buffer[i].length;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  return edit_distance;
}
uint64_t matches_cigar_compute_edit_distance__excluding_clipping(
    vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length) {
  // Sum up all cigar elements
  const cigar_element_t* const cigar_buffer = vector_get_elm(cigar_vector,cigar_buffer_offset,cigar_element_t);
  uint64_t i, edit_distance = 0;
  for (i=0;i<cigar_length;++i) {
    switch (cigar_buffer[i].type) {
      case cigar_match: break;
      case cigar_mismatch:
        edit_distance++;
        break;
      case cigar_del:
        if (i==0 || i==cigar_length-1) break;
      // No break
      case cigar_ins:
        edit_distance += cigar_buffer[i].length;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  return edit_distance;
}
uint64_t matches_cigar_compute_matching_bases(
    vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length) {
  // Sum up all cigar elements
  const cigar_element_t* const cigar_buffer = vector_get_elm(cigar_vector,cigar_buffer_offset,cigar_element_t);
  uint64_t i, matching_bases = 0;
  for (i=0;i<cigar_length;++i) {
    switch (cigar_buffer[i].type) {
      case cigar_match:
        matching_bases += cigar_buffer[i].length;
        break;
      case cigar_mismatch:
      case cigar_ins:
      case cigar_del:
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  return matching_bases;
}
int64_t matches_cigar_element_effective_length(
    const cigar_element_t* const cigar_element) {
  switch (cigar_element->type) {
    case cigar_match:
      return cigar_element->length;
      break;
    case cigar_mismatch:
      return 1;
      break;
    case cigar_del:
      break;
    case cigar_ins:
      return cigar_element->length;
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  return 0;
}
int64_t matches_cigar_compute_effective_length(
    vector_t* const cigar_vector,
    const uint64_t cigar_offset,
    const uint64_t cigar_length) {
  uint64_t effective_length = 0;
  // Traverse all CIGAR elements
  const cigar_element_t* const cigar_buffer = vector_get_elm(cigar_vector,cigar_offset,cigar_element_t);
  uint64_t i;
  for (i=0;i<cigar_length;++i) {
    effective_length += matches_cigar_element_effective_length(cigar_buffer+i);
  }
  // Return effective length
  return effective_length;
}
float matches_cigar_compute_error_quality(
    vector_t* const cigar_vector,
    const uint64_t cigar_offset,
    const uint64_t cigar_length,
    uint8_t* const quality_mask,
    const uint64_t quality_mask_length){
  // No quality mask or Exact Match
  if (quality_mask==NULL || cigar_length==0) return (float)SEQUENCE_QUALITIES_MAX;
  // Traverse CIGAR
  const cigar_element_t* cigar_element = vector_get_elm(cigar_vector,cigar_offset,cigar_element_t);
  uint64_t error_quality = 0, error_count = 0;
  uint64_t i, j, read_pos = 0;
  for (i=0;i<cigar_length;++i,++cigar_element) {
    switch (cigar_element->type) {
      case cigar_match:
        read_pos += cigar_element->length;
        break;
      case cigar_mismatch:
        error_quality += quality_mask[read_pos];
        ++error_count;
        ++read_pos;
        break;
      case cigar_ins:
        if (i==0 || i==cigar_length-1) break; // Skip
        error_quality += quality_mask[read_pos];
        ++error_count;
        break;
      case cigar_del:
        if (i>0 && i<cigar_length-1) { // Skip trims
          for (j=0;j<cigar_element->length;++j,++error_count) {
            error_quality += quality_mask[read_pos+j];
          }
        }
        read_pos += cigar_element->length;
        break;
      default:
        break;
    }
  }
  // Return average
  return (error_count!=0) ? (float)error_quality/(float)error_count : (float)SEQUENCE_QUALITIES_MAX;
}
/*
 * CIGAR Vector Compare
 */
int matches_cigar_cmp(
    vector_t* const cigar0_vector,
    const uint64_t cigar0_offset,
    const uint64_t cigar0_length,
    vector_t* const cigar1_vector,
    const uint64_t cigar1_offset,
    const uint64_t cigar1_length) {
  if (cigar0_length != cigar1_length) return -1;
  // Locate CIGARs
  cigar_element_t* const match0_cigar = vector_get_elm(cigar0_vector,cigar0_offset,cigar_element_t);
  cigar_element_t* const match1_cigar = vector_get_elm(cigar1_vector,cigar1_offset,cigar_element_t);
  // Compare
  uint64_t i;
  for (i=0;i<cigar0_length;++i) {
    if (match0_cigar[i].type != match1_cigar[i].type) return -1;
    if (match0_cigar[i].attributes != match1_cigar[i].attributes) return -1;
    switch (match0_cigar[i].type) {
      case cigar_match:
      case cigar_ins:
      case cigar_del:
        if (match0_cigar[i].length != match1_cigar[i].length) return -1;
        break;
      case cigar_mismatch:
        if (match0_cigar[i].mismatch != match1_cigar[i].mismatch) return -1;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  return 0;
}
/*
 * Display
 */
void match_cigar_print(
    FILE* const stream,
    vector_t* const cigar_vector,
    const uint64_t cigar_buffer_offset,
    const uint64_t cigar_length) {
  uint64_t j;
  for (j=0;j<cigar_length;++j) {
    cigar_element_t* const cigar_element = vector_get_elm(cigar_vector,cigar_buffer_offset+j,cigar_element_t);
    switch (cigar_element->type) {
      case cigar_match: fprintf(stream,"%dM",cigar_element->length); break;
      case cigar_mismatch: fprintf(stream,"1X"); break;
      case cigar_ins:
        if (cigar_element->attributes == cigar_attr_homopolymer) {
          fprintf(stream,"%di",cigar_element->length);
        } else {
          fprintf(stream,"%dI",cigar_element->length);
        }
        break;
      case cigar_del:
        if (cigar_element->attributes == cigar_attr_homopolymer) {
          fprintf(stream,"%dd",cigar_element->length);
        } else if (cigar_element->attributes == cigar_attr_trim) {
          fprintf(stream,"%dS",cigar_element->length);
        } else {
          fprintf(stream,"%dD",cigar_element->length);
        }
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}

