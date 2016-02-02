/*
 * PROJECT: GEMMapper
 * FILE: archive_text_rl.c
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "archive_text_rl.h"
#include "sampled_rl.h"

/*
 * Encode RL-Text
 */
void archive_text_rl_encode(
    const uint8_t* const text,const uint64_t text_length,
    uint8_t* const rl_text,uint8_t* const rl_runs,
    uint64_t* const rl_text_length) {
  // Init
  rl_text[0] = text[0];
  rl_runs[0] = 1;
  // RL-Encode all text
  uint64_t rl_pos = 0, text_pos, acc = 0;
  for (text_pos=1;text_pos<text_length;++text_pos) {
    if (rl_text[rl_pos]==text[text_pos] &&
        rl_runs[rl_pos] < TEXT_RL_MAX_RUN_LENGTH) {
      // Increase run
      ++rl_runs[rl_pos];
    } else {
      // Store accumulated
      rl_runs[rl_pos] += acc;
      acc = rl_runs[rl_pos];
      // Next run
      ++rl_pos;
      rl_text[rl_pos] = text[text_pos];
      rl_runs[rl_pos] = 1;
    }
  }
  // Store accumulated
  rl_runs[rl_pos] += acc;
  // Set rl-length
  *rl_text_length = rl_pos+1;
}
/*
 * Translate position
 */
uint64_t archive_text_rl_translate(
    archive_text_t* const archive_text,const uint64_t position_rl,
    mm_stack_t* const mm_stack) {
  // Retrieve first sampled position
  const uint64_t sampled_position_rl = position_rl / SAMPLED_RL_SAMPLING_RATE;
  const uint64_t remainder_rl = position_rl % SAMPLED_RL_SAMPLING_RATE;
  const uint64_t remainder_max_expanded_length = remainder_rl * TEXT_RL_MAX_RUN_LENGTH;
  // Translate sampled position
  const uint64_t sampled_position = sampled_rl_get_sample(archive_text->sampled_rl,sampled_position_rl);
  // Retrieve text-RL
  text_trace_t text_trace;
  mm_stack_push_state(mm_stack);
  archive_text_retrieve(archive_text,sampled_position,
      remainder_max_expanded_length,false,true,&text_trace,mm_stack);
  // Compute offset from @sampled_position_rl to @position_rl
  const uint64_t remainder = (remainder_rl>0) ?
      text_trace.rl_runs[remainder_rl-1] : text_trace.rl_runs[remainder_rl];
  // Return
  mm_stack_pop_state(mm_stack);
  return sampled_position + remainder;
}
/*
 * Utils
 */
uint64_t archive_text_rl_decode_run_length(uint8_t* const rl_runs,const uint64_t rl_position) {
  return (rl_position>0) ? rl_runs[rl_position]-rl_runs[rl_position-1] : rl_runs[rl_position];
}
uint64_t archive_text_rl_position_decode_inc(uint8_t* const rl_runs,const uint64_t rl_position) {
  return rl_runs[rl_position];
}
uint64_t archive_text_rl_position_decode_exl(uint8_t* const rl_runs,const uint64_t rl_position) {
  return (rl_position>0) ? rl_runs[rl_position-1] : 0;
}
uint64_t archive_text_rl_position_encode(
    uint8_t* const rl_runs,const uint64_t rl_text_length,
    const uint64_t text_position) { // TODO Impl log2 solution
  uint64_t i;
  for (i=0;i<rl_text_length;++i) {
    if (rl_runs[i] > text_position) return i;
  }
  return rl_text_length;
}


