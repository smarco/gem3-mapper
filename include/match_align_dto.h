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

#include "essentials.h"
#include "align_bpm_pattern.h"
#include "align_swg_score.h"
#include "align_swg_simd.h"
#include "text_collection.h"

struct match_align_input_t {
  /* Key */
  uint8_t* key;
  uint64_t key_length;
  uint64_t key_trim_left;
  uint64_t key_trim_right;
  bpm_pattern_t* bpm_pattern;
  bpm_pattern_t* bpm_pattern_tiles;
  swg_query_profile_t* swg_query_profile;
  /* Text */
  uint64_t text_trace_offset;              // Text-Trace Offset
  uint64_t text_position_translated;       // Translated text-position from RL-Space (In case of RL-text)
  uint64_t text_position;                  // Text position
  uint8_t* text;                           // Text (candidate)
  uint64_t text_length;                    // Text full length (whole decoded text)
  uint64_t text_offset_begin;              // Text effective-begin offset (eg. bounded by verification)
  uint64_t text_offset_end;                // Text effective-end offset   (eg. bounded by verification)
  uint64_t text_offset_base_begin;         // Text base begin offset (range projected by key)
  uint64_t text_offset_base_end;           // Text base end offset   (range projected by key)
  uint64_t align_distance_bound;
  /* RL-Input */
  bool run_length;
  uint8_t* rl_key_runs;
  uint8_t* rl_text_runs;
};
struct match_align_parameters_t {
  bool emulated_rc_search;
  bool* allowed_enc;
  swg_penalties_t* swg_penalties;
  uint64_t max_error;
  uint64_t max_bandwidth;
  uint64_t global_min_identity;
  int64_t global_min_swg_threshold;
  uint64_t local_min_identity;
  int64_t local_min_swg_threshold;
  bool left_gap_alignment;
  bool scaffolding;
  uint64_t scaffolding_matching_min_length;
  uint64_t scaffolding_min_coverage;
  bool cigar_curation;
  uint64_t cigar_curation_min_end_context;
};

#endif /* MATCH_ALIGN_DTO_H_ */
