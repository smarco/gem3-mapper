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

#ifndef SEQUENCE_QUALITIES_MODEL_H_
#define SEQUENCE_QUALITIES_MODEL_H_

#include "utils/essentials.h"
#include "text/sequence.h"

/*
 * Quality Data Structures
 */
typedef enum {
  sequence_qualities_model_flat,
  sequence_qualities_model_gem
} sequence_qualities_model_t;
typedef enum {
  sequence_qualities_ignore,
  sequence_qualities_offset_33,
  sequence_qualities_offset_64
} sequence_qualities_format_t;
typedef enum {
  qm_real=0,
  qm_pseudo=1
} quality_type_t;

/*
 * Quality Models
 */
void sequence_qualities_model_process(
    sequence_t* const sequence,
    const sequence_qualities_model_t qualities_model,
    const sequence_qualities_format_t quality_format,
    const uint64_t quality_threshold,
    uint8_t* const quality_mask);

#endif /* SEQUENCE_QUALITIES_MODEL_H_ */
