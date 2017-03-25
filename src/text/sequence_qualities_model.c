/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2011-2017 by Paolo Ribeca  <paolo.ribeca@gmail.com>
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
 *            Paolo Ribeca <paolo.ribeca@gmail.com>
 */

#include "text/sequence_qualities_model.h"

/*
 * Error Msg
 */
#define GEM_ERROR_QUALITY_NEGATIVE "Negative quality ('%c') at position %"PRIu64""

/*
 * Quality Model (from sequence into qm_values)
 */
void sequence_qualities_model_process(
    sequence_t* const sequence,
    const sequence_qualities_format_t quality_format,
    const uint64_t quality_threshold,
    uint8_t* const quality_mask) {
  uint64_t enc_diff;
  switch (quality_format) {
    case sequence_qualities_offset_33:
      enc_diff = 33;
      break;
    case sequence_qualities_offset_64:
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
    gem_cond_error(quality < 0,QUALITY_NEGATIVE,string_get_buffer(&sequence->qualities)[i],i);
    quality_mask[i] = (quality >= quality_threshold) ? qm_real : qm_pseudo;
  }
}

