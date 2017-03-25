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
 * DESCRIPTION:
 *   Archive text module provides data structures and functions to
 *   manipulate and query a genomic-text (ACGTN) encoded run-length (RL)
 */

#include "archive/archive_text_rl.h"
#include "archive/sampled_rl.h"

/*
 * Encode RL-Text
 */
void archive_text_rl_encode(
    const uint8_t* const text,
    const uint64_t text_length,
    uint8_t* const rl_text,
    uint64_t* const rl_text_length,
    uint32_t* const rl_runs_acc) {
  // Init
  rl_text[0] = text[0];
  rl_runs_acc[0] = 1;
  // RL-Encode all text
  uint64_t rl_pos = 0, text_pos, acc = 0;
  for (text_pos=1;text_pos<text_length;++text_pos) {
    if (rl_text[rl_pos]==text[text_pos] &&
        rl_runs_acc[rl_pos] < TEXT_RL_MAX_RUN_LENGTH) {
      // Increase run
      ++rl_runs_acc[rl_pos];
    } else {
      // Store accumulated
      rl_runs_acc[rl_pos] += acc;
      acc = rl_runs_acc[rl_pos];
      // Next run
      ++rl_pos;
      rl_text[rl_pos] = text[text_pos];
      rl_runs_acc[rl_pos] = 1;
    }
  }
  // Store accumulated
  rl_runs_acc[rl_pos] += acc;
  // Set rl-length
  *rl_text_length = rl_pos+1;
}
/*
 * Translate position
 */
uint64_t archive_text_rl_position_translate(
    archive_text_t* const archive_text,
    const uint64_t position_rl,
    mm_allocator_t* const mm_allocator) {
  // Retrieve first sampled position
  const uint64_t sampled_position_rl = position_rl / SAMPLED_RL_SAMPLING_RATE;
  const uint64_t remainder_rl = position_rl % SAMPLED_RL_SAMPLING_RATE;
  // Translate sampled position
  const uint64_t sampled_position = sampled_rl_get_sample(archive_text->sampled_rl,sampled_position_rl);
  if (remainder_rl==0) {
    return sampled_position;
  } else {
    // Retrieve text-RL
    mm_allocator_push_state(mm_allocator);
    const uint64_t remainder_max_expanded_length = remainder_rl * TEXT_RL_MAX_RUN_LENGTH;
    text_trace_t text_trace;
    archive_text_retrieve(archive_text,sampled_position,
        remainder_max_expanded_length,false,true,&text_trace,mm_allocator);
    // Compute offset from @sampled_position_rl to @position_rl
    const uint64_t remainder = text_trace.rl_runs_acc[remainder_rl-1];
    // Return
    mm_allocator_pop_state(mm_allocator);
    return sampled_position + remainder;
  }
}
/*
 * Utils
 */
uint64_t archive_text_rl_get_run_length(
    const uint32_t* const rl_runs_acc,
    const uint64_t rl_position) {
  return (rl_position>0) ? rl_runs_acc[rl_position]-rl_runs_acc[rl_position-1] : rl_runs_acc[rl_position];
}
uint64_t archive_text_rl_get_decoded_offset_inc(
    const uint32_t* const rl_runs_acc,
    const uint64_t rl_position) {
  return rl_runs_acc[rl_position];
}
uint64_t archive_text_rl_get_decoded_offset_exl(
    const uint32_t* const rl_runs_acc,
    const uint64_t rl_position) {
  return (rl_position>0) ? rl_runs_acc[rl_position-1] : 0;
}
uint64_t archive_text_rl_get_decoded_length(
    const uint32_t* const rl_runs_acc,
    const uint64_t rl_position,
    const uint64_t length) {
  return rl_runs_acc[rl_position+length-1] - archive_text_rl_get_decoded_offset_exl(rl_runs_acc,rl_position);
}
uint64_t archive_text_rl_get_encoded_offset(
    const uint32_t* const rl_runs_acc,
    const uint64_t rl_text_length,
    const uint64_t text_position) { // TODO Impl log2 solution
  uint64_t i;
  for (i=0;i<rl_text_length;++i) {
    if (rl_runs_acc[i] > text_position) return i;
  }
  return rl_text_length;
}


