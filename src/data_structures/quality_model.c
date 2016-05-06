/*
 * PROJECT: GEMMapper
 * FILE: quality_model.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 *            Paolo Ribeca <paolo.ribeca@gmail.com>
 * DESCRIPTION:
 */

#include "data_structures/quality_model.h"

/*
 * Error Msg
 */
#define GEM_ERROR_QUALITY_NEGATIVE "Negative quality ('%c') at position %"PRIu64""

/*
 * Quality Model (Models quality scores from sequence into qm_values)
 */
void quality_model(
    sequence_t* const sequence,
    const quality_model_t quality_model,
    const quality_format_t quality_format,
    const uint64_t quality_threshold,
    uint8_t* const quality_mask) {
  switch (quality_model) {
    case quality_model_type_flat:
      quality_model_flat(sequence,quality_format,quality_threshold,quality_mask);
      break;
    case quality_model_type_gem:
      quality_model_gem(sequence,quality_format,quality_threshold,quality_mask);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
void quality_model_flat(
    sequence_t* const sequence,
    const quality_format_t quality_format,
    const uint64_t quality_threshold,
    uint8_t* const quality_mask) {
  const uint64_t key_length = sequence_get_length(sequence);
  uint64_t i;
  for (i=0;i<key_length;++i) quality_mask[i] = qm_real;
}
void quality_model_gem(
    sequence_t* const sequence,
    const quality_format_t quality_format,
    const uint64_t quality_threshold,
    uint8_t* const quality_mask) {
  uint64_t enc_diff;
  switch (quality_format) {
    case qualities_offset_33:
      enc_diff = 33;
      break;
    case qualities_offset_64:
      enc_diff = 64;
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  const uint64_t key_length = sequence_get_length(sequence);
  uint64_t i;
  for (i=0;i<key_length;++i) {
    const int64_t quality = *string_char_at(&sequence->qualities,i) - enc_diff;
    gem_cond_fatal_error(quality < 0,QUALITY_NEGATIVE,string_get_buffer(&sequence->qualities)[i],i);
    quality_mask[i] = (quality >= quality_threshold) ? qm_real : qm_pseudo;
  }
}
