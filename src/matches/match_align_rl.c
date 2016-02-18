/*
 * PROJECT: GEMMapper
 * FILE: match_align_rl.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "matches/match_align_rl.h"
#include "matches/matches_cigar.h"
#include "archive/archive_text_rl.h"

/*
 * RL-Translate CIGAR elements
 */
void match_align_rl_translate_cigar_match(
    cigar_element_t** const cigar_buffer,
    const cigar_element_t* const rl_cigar_element,
    const uint8_t* const rl_key_runs,
    const uint8_t* const rl_text_runs,
    uint64_t rl_text_pos,
    uint64_t rl_key_pos,
    const bool left_gap_alignment) {
  // Traverse all matching characters
  uint64_t i;
  for (i=0;i<rl_cigar_element->length;++i) {
    // Check Run Length
    const uint64_t rl_key = archive_text_rl_get_run_length(rl_key_runs,rl_key_pos);
    const uint64_t rl_text = archive_text_rl_get_run_length(rl_text_runs,rl_text_pos);
    // Compute the CIGAR-translated
    cigar_t cigar_op;
    uint64_t match_length, indel_length;
    if (rl_key > rl_text) {
      cigar_op = cigar_del;
      indel_length = rl_key - rl_text;
      match_length = rl_text;
    } else {
      cigar_op = cigar_ins;
      indel_length = rl_text - rl_key;
      match_length = rl_key;
    }
    // Add CIGAR-translated
    if (left_gap_alignment) {
      if (indel_length > 0) {
        matches_cigar_buffer_add_cigar_element__attr(cigar_buffer,cigar_op,indel_length,cigar_attr_homopolymer);
      }
      matches_cigar_buffer_add_cigar_element(cigar_buffer,cigar_match,match_length);
    } else {
      matches_cigar_buffer_add_cigar_element(cigar_buffer,cigar_match,match_length);
      if (indel_length > 0) {
        matches_cigar_buffer_add_cigar_element__attr(cigar_buffer,cigar_op,indel_length,cigar_attr_homopolymer);
      }
    }
    // Next
    ++rl_key_pos;
    ++rl_text_pos;
  }
}
void match_align_rl_translate_cigar_mismatch(
    cigar_element_t** const cigar_buffer,
    const cigar_element_t* const rl_cigar_element,
    const uint8_t* const rl_key_runs,
    const uint8_t* const rl_text_runs,
    const uint64_t rl_text_pos,
    const uint64_t rl_key_pos,
    const bool left_gap_alignment) {
  // Check Run Length
  const uint64_t rl_key = archive_text_rl_get_run_length(rl_key_runs,rl_key_pos);
  const uint64_t rl_text = archive_text_rl_get_run_length(rl_text_runs,rl_text_pos);
  // Compute the CIGAR-translated
  cigar_t cigar_op;
  uint64_t mismatch_length, indel_length, i;
  if (rl_key > rl_text) {
    cigar_op = cigar_del;
    mismatch_length = rl_text;
    indel_length = rl_key-rl_text;
  } else {
    cigar_op = cigar_ins;
    mismatch_length = rl_key;
    indel_length = rl_text-rl_key;
  }
  // Add CIGAR-translated
  if (left_gap_alignment) {
    if (indel_length > 0)  {
      matches_cigar_buffer_add_cigar_element__attr(cigar_buffer,cigar_op,indel_length,cigar_attr_homopolymer);
    }
    for (i=0;i<mismatch_length;++i) matches_cigar_buffer_add_mismatch(cigar_buffer,rl_cigar_element->mismatch);
  } else {
    for (i=0;i<mismatch_length;++i) matches_cigar_buffer_add_mismatch(cigar_buffer,rl_cigar_element->mismatch);
    if (indel_length > 0)  {
      matches_cigar_buffer_add_cigar_element__attr(cigar_buffer,cigar_op,indel_length,cigar_attr_homopolymer);
    }
  }
}
void match_align_rl_translate_cigar_ins(
    cigar_element_t** const cigar_buffer,
    const cigar_element_t* const rl_cigar_element,
    const uint8_t* const rl_text_runs,
    const uint64_t rl_text_pos) {
  // Compute the CIGAR-translated
  const uint64_t indel_length = archive_text_rl_get_decoded_length(rl_text_runs,rl_text_pos,rl_cigar_element->length);
  // Add CIGAR-translated
  matches_cigar_buffer_add_cigar_element(cigar_buffer,cigar_ins,indel_length);
}
void match_align_rl_translate_cigar_del(
    cigar_element_t** const cigar_buffer,
    const cigar_element_t* const rl_cigar_element,
    const uint8_t* const rl_key_runs,
    const uint64_t rl_key_pos) {
  // Compute the CIGAR-translated
  const uint64_t indel_length = archive_text_rl_get_decoded_length(rl_key_runs,rl_key_pos,rl_cigar_element->length);
  // Add CIGAR-translated
  matches_cigar_buffer_add_cigar_element(cigar_buffer,cigar_del,indel_length);
}
/*
 * RL-Translate CIGAR
 */
void match_align_rl_translate_cigar(
    match_trace_t* const match_trace,
    vector_t* const cigar_vector,
    match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters) {
  // Parameters
  const bool left_gap_alignment = align_parameters->left_gap_alignment;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  const uint64_t key_length = align_input->key_length;
  const uint8_t* const rl_key_runs = align_input->rl_key_runs;
  const uint8_t* const rl_text_runs = align_input->rl_text_runs;
  // RL-CIGAR
  const uint64_t rl_cigar_length = match_alignment->cigar_length;
  cigar_element_t* const rl_cigar_buffer = vector_get_elm(cigar_vector,match_alignment->cigar_offset,cigar_element_t);
  // Allocate Translated CIGAR
  const uint64_t cigar_offset = vector_get_used(cigar_vector);
  vector_reserve_additional(cigar_vector,rl_key_runs[key_length-1]);
  cigar_element_t* cigar_buffer = vector_get_free_elm(cigar_vector,cigar_element_t); // Sentinel
  cigar_element_t* const cigar_buffer_base = cigar_buffer;
  cigar_buffer->type = cigar_null; // Trick
  // Translate all CIGAR elements
  uint64_t text_pos = match_trace->match_alignment.match_text_offset, key_pos = 0;
  uint64_t i;
  for (i=0;i<rl_cigar_length;++i) {
    cigar_element_t* rl_cigar_element = rl_cigar_buffer + i;
    switch (rl_cigar_element->type) {
      case cigar_match:
        match_align_rl_translate_cigar_match(&cigar_buffer,rl_cigar_element,
            rl_key_runs,rl_text_runs,text_pos,key_pos,left_gap_alignment);
        text_pos += rl_cigar_element->length;
        key_pos += rl_cigar_element->length;
        break;
      case cigar_mismatch:
        match_align_rl_translate_cigar_mismatch(&cigar_buffer,rl_cigar_element,
            rl_key_runs,rl_text_runs,text_pos,key_pos,left_gap_alignment);
        ++text_pos;
        ++key_pos;
        break;
      case cigar_ins:
        match_align_rl_translate_cigar_ins(&cigar_buffer,rl_cigar_element,rl_text_runs,text_pos);
        text_pos += rl_cigar_element->length;
        break;
      case cigar_del:
        match_align_rl_translate_cigar_del(&cigar_buffer,rl_cigar_element,rl_key_runs,key_pos);
        key_pos += rl_cigar_element->length;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  // Set CIGAR buffer used
  if (cigar_buffer->type!=cigar_null) ++(cigar_buffer);
  const uint64_t num_cigar_elements = cigar_buffer - cigar_buffer_base;
  vector_add_used(cigar_vector,num_cigar_elements);
  // Setup translated CIGAR
  match_alignment->cigar_length = num_cigar_elements;
  match_alignment->cigar_offset = cigar_offset;
}
