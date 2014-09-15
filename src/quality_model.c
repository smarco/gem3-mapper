/*
 * PROJECT: GEMMapper
 * FILE: quality_model.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 *            Paolo Ribeca <paolo.ribeca@gmail.com>
 * DESCRIPTION:
 */

#include "quality_model.h"

/*
 * Quality Model (Models quality scores from sequence into qm_values)
 */
GEM_INLINE void quality_model(
    sequence_t* const sequence,const quality_model_t quality_model,
    const quality_format_t quality_format,const uint64_t quality_threshold,uint8_t* const quality_mask) {
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
GEM_INLINE void quality_model_flat(
    sequence_t* const sequence,
    const quality_format_t quality_format,const uint64_t quality_threshold,uint8_t* const quality_mask) {
  // TODO CHECK sequence not null
  const uint64_t key_length = sequence_get_length(sequence);
  uint64_t i;
  for (i=0;i<key_length;++i) quality_mask[i] = qm_real;
}
GEM_INLINE void quality_model_gem(
    sequence_t* const sequence,
    const quality_format_t quality_format,const uint64_t quality_threshold,uint8_t* const quality_mask) {
  // TODO CHECK sequence not null
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







//GEM_INLINE void quality_score_flat(
//    const sequence_t* const sequence,const quality_format_t format,const uint64_t quality_threshold,
//    const mismatch* const misms,const uint64_t mism_num) {
//  return mism_num;
//}

//GEM_INLINE void quality_score_gem(const sequence_t* const sequence,
// const quality_format_t format,const uint64_t quality_threshold,
// const mismatch* const misms,const uint64_t mism_num) {
//
//  register uint64_t enc_diff, res = 0, i;
//  switch (format) {
//  case Offset_33:
//    enc_diff = 33;
//    break;
//  case Offset_64:
//    enc_diff = 64;
//    break;
//  default:
//    assert(0);
//    break;
//  }
//  for (i = 0; i < mism_num; ++i)
//    if (misms[i].mismatch < 256) {
//      register const int quality=qualities[misms[i].position]-enc_diff;
//      /* This should have been already checked before by the model */
//      gem_cond_fatal_error((quality < 0),QUAL_NEG,qualities[misms[i].position]);
//      res += quality;
//    }
//  return res;
//}
