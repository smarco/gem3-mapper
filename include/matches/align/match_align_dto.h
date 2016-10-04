/*
 * PROJECT: GEMMapper
 * FILE: match_align_dto.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */
#ifndef MATCH_ALIGN_DTO_H_
#define MATCH_ALIGN_DTO_H_

typedef struct match_align_input_t match_align_input_t;
typedef struct match_align_parameters_t match_align_parameters_t;

#include "utils/essentials.h"
#include "align/alignment.h"
#include "align/align_bpm_pattern.h"
#include "align/align_swg_score.h"
#include "align/align_swg_simd.h"
#include "text/text_collection.h"

struct match_align_input_t {
  /* Sequence */
  uint64_t sequence_clip_left;             // Input sequence clipped (masked)
  uint64_t sequence_clip_right;            // Input sequence clipped (masked)
  /* Key */
  uint8_t* key;
  uint64_t key_length;
  uint64_t key_trim_left;                  // Alignment trim
  uint64_t key_trim_right;                 // Alignment trim
  bpm_pattern_t* bpm_pattern;
  bpm_pattern_t* bpm_pattern_tiles;
  swg_query_profile_t* swg_query_profile;
  /* Text */
  uint64_t text_trace_offset;              // Text-Trace Offset
  uint64_t text_position;                  // Text position
  uint8_t* text;                           // Text (candidate)
  uint64_t text_length;                    // Text full length (whole decoded text)
  alignment_t* alignment;                  // Text alignment
  /* RL-Input */
  bool run_length;
  uint32_t* rl_key_runs_acc;
  uint32_t* rl_text_runs_acc;
};
struct match_align_parameters_t {
  bool* allowed_enc;
  swg_penalties_t* swg_penalties;
  uint64_t max_error;
  uint64_t max_bandwidth;
  uint64_t global_min_identity;
  int64_t global_min_swg_threshold;
  bool local_alignment;
  uint64_t max_aligned_gap_length;
  uint64_t local_min_identity;
  int64_t local_min_swg_threshold;
  bool left_gap_alignment;
  bool force_full_swg;
  uint64_t scaffolding_matching_min_length;
  uint64_t scaffolding_min_coverage;
  bool cigar_curation;
  uint64_t cigar_curation_min_end_context;
};

#endif /* MATCH_ALIGN_DTO_H_ */
