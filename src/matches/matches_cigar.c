/*
 * PROJECT: GEMMapper
 * FILE: matches_cigar.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Data structure to store alignment matches {sequence,position,strand,CIGAR}
 */

#include "matches/matches_cigar.h"
#include "matches/matches_classify.h"

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
    // FIXME In case of indel, flip @origin->indel.indel_text (only SAM.MD field is using it)
  }
  if (cigar_length%2) {
    cigar_element_t* const middle = cigar_buffer + middle_point;
    if (middle->type == cigar_mismatch) middle->mismatch = dna_encoded_complement(middle->mismatch);
    // FIXME In case of indel, flip @origin->indel.indel_text (only SAM.MD field is using it)
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
        edit_distance += cigar_buffer[i].type;
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
        edit_distance += cigar_buffer[i].type;
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
int64_t matches_cigar_element_effective_length(const cigar_element_t* const cigar_element) {
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
int64_t matches_cigar_effective_length(
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
/*
 * CIGAR Vector Compare
 */
int matches_cigar_cmp(
    vector_t* const cigar_vector_match0,
    match_trace_t* const match0,
    vector_t* const cigar_vector_match1,
    match_trace_t* const match1) {
  const uint64_t match0_cigar_length = match0->match_alignment.cigar_length;
  const uint64_t match1_cigar_length = match1->match_alignment.cigar_length;
  if (match0_cigar_length != match1_cigar_length) return -1;
  // Locate CIGARs
  cigar_element_t* const match0_cigar = vector_get_elm(cigar_vector_match0,match0->match_alignment.cigar_offset,cigar_element_t);
  cigar_element_t* const match1_cigar = vector_get_elm(cigar_vector_match1,match1->match_alignment.cigar_offset,cigar_element_t);
  // Compare
  uint64_t i;
  for (i=0;i<match0_cigar_length;++i) {
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

