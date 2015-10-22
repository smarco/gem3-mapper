/*
 * PROJECT: GEMMapper
 * FILE: align_swg_score.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */


#include "align_swg_score.h"

/*
 * SWG Score
 */
GEM_INLINE int32_t align_swg_score_deletion(const swg_penalties_t* const swg_penalties,const int32_t length) {
  const int32_t gap_open_score = swg_penalties->gap_open_score;
  const int32_t gap_extension = swg_penalties->gap_extension_score;
  return gap_open_score + gap_extension*length;
}
GEM_INLINE int32_t align_swg_score_insertion(const swg_penalties_t* const swg_penalties,const int32_t length) {
  const int32_t gap_open_score = swg_penalties->gap_open_score;
  const int32_t gap_extension = swg_penalties->gap_extension_score;
  return gap_open_score + gap_extension*length;
}
GEM_INLINE int32_t align_swg_score_mismatch(const swg_penalties_t* const swg_penalties) {
  return swg_penalties->generic_mismatch_score;
}
GEM_INLINE int32_t align_swg_score_match(const swg_penalties_t* const swg_penalties,const int32_t match_length) {
  return swg_penalties->generic_match_score * match_length;
}
GEM_INLINE int32_t align_swg_score_cigar_element(
    const swg_penalties_t* const swg_penalties,const cigar_element_t* const cigar_element) {
  switch (cigar_element->type) {
    case cigar_match:
      return align_swg_score_match(swg_penalties,cigar_element->length);
      break;
    case cigar_mismatch:
      return align_swg_score_mismatch(swg_penalties);
      break;
    case cigar_ins:
      return align_swg_score_insertion(swg_penalties,cigar_element->length);
      break;
    case cigar_del:
      return align_swg_score_deletion(swg_penalties,cigar_element->length);
      break;
    case cigar_null:
    default:
      GEM_INVALID_CASE();
      break;
  }
  return 0;
}
GEM_INLINE int32_t align_swg_score_cigar(
    const swg_penalties_t* const swg_penalties,vector_t* const cigar_vector,
    const uint64_t cigar_offset,const uint64_t cigar_length) {
  // Parameters
  const cigar_element_t* const cigar_buffer = vector_get_elm(cigar_vector,cigar_offset,cigar_element_t);
  // Traverse all CIGAR elements
  int32_t score = 0, i;
  for (i=0;i<cigar_length;++i) {
    score += align_swg_score_cigar_element(swg_penalties,cigar_buffer+i);
  }
  // Return score
  return score;
}
GEM_INLINE int32_t align_swg_score_cigar__excluding_clipping(
    const swg_penalties_t* const swg_penalties,vector_t* const cigar_vector,
    const uint64_t cigar_offset,const uint64_t cigar_length) {
  // Parameters
  const cigar_element_t* const cigar_buffer = vector_get_elm(cigar_vector,cigar_offset,cigar_element_t);
  // Ignore trims
  if (cigar_length==0) return 0;
  const int64_t last_cigar_element = cigar_length-1;
  const int64_t begin_element = (cigar_buffer[0].type==cigar_del) ? 1 : 0;
  const int64_t last_element = (cigar_buffer[last_cigar_element].type==cigar_del) ? last_cigar_element-1 : last_cigar_element;
  // Traverse all CIGAR elements
  int32_t score = 0, i;
  for (i=begin_element;i<=last_element;++i) {
    score += align_swg_score_cigar_element(swg_penalties,cigar_buffer+i);
  }
  // Return score
  return score;
}
