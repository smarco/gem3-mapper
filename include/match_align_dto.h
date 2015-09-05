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
#include "bpm_align.h"
#include "swg_align.h"

struct match_align_input_t {
  /* Key */
  uint8_t* key;
  uint64_t key_length;
  bpm_pattern_t* bpm_pattern;
  swg_query_profile_t* swg_query_profile;
  /* Text */
  uint64_t text_trace_offset;
  uint64_t text_position;
  uint8_t* text;
  uint64_t text_length;
  uint64_t text_offset_begin;
  uint64_t text_offset_end;
};
struct match_align_parameters_t {
  bool emulated_rc_search;
  bool* allowed_enc;
  swg_penalties_t* swg_penalties;
  int64_t swg_threshold;
  uint64_t max_error;
  uint64_t min_identity;
  uint64_t max_bandwidth;
  bool left_gap_alignment;
  bool scaffolding;
  uint64_t scaffolding_matching_min_length;
  uint64_t scaffolding_homopolymer_min_context;
  uint64_t scaffolding_min_coverage;
  bool cigar_curation;
  uint64_t cigar_curation_min_end_context;
};

#endif /* MATCH_ALIGN_DTO_H_ */
