/*
 * PROJECT: GEMMapper
 * FILE: archive_search.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_search.h"
#include "archive_select.h"
#include "archive_score.h"
#include "matches_classify.h"

/*
 * Debug
 */
#define DEBUG_ARCHIVE_SEARCH_READ_NAME GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * Archive Search State
 */
const char* archive_search_pe_state_label[] =
{
    [0] = "begin",
    [1] = "search-end1",
    [2] = "search-end2",
    [3] = "recovery",
    [4] = "find-pairs",
    [5] = "end",
};

/*
 * Setup
 */
GEM_INLINE archive_search_t* archive_search_new(
    archive_t* const archive,search_parameters_t* const search_parameters,
    select_parameters_t* const select_parameters) {
  // Allocate handler
  archive_search_t* const archive_search = mm_alloc(archive_search_t);
  // Archive
  archive_search->archive = archive;
  // Sequence
  sequence_init(&archive_search->sequence);
  sequence_init(&archive_search->rc_sequence);
  // Approximate Search
  archive_search->as_parameters.search_parameters = search_parameters;
  archive_search->select_parameters = select_parameters;
  approximate_search_init(
      &archive_search->forward_search_state,archive,
      &archive_search->as_parameters,false);
  approximate_search_init(
      &archive_search->reverse_search_state,archive,
      &archive_search->as_parameters,true);
  // Archive search control (Flow control) [DEFAULTS]
  archive_search->probe_strand = true;
  archive_search->emulate_rc_search = !archive->indexed_complement;
  // Return
  return archive_search;
}
GEM_INLINE void archive_search_prepare_sequence(archive_search_t* const archive_search) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_PREPARE_SEQUENCE,PROFILE_LEVEL);
  // Check the index characteristics & generate reverse-complement (if needed)
  if (archive_search->archive->indexed_complement) {
    archive_search->emulate_rc_search = false;
  } else {
    sequence_generate_reverse_complement(&archive_search->sequence,&archive_search->rc_sequence);
    archive_search->emulate_rc_search = !sequence_equals(&archive_search->sequence,&archive_search->rc_sequence);
  }
  // Generate the pattern(s)
  const bool prepare_rl_pattern = archive_search->archive->text->run_length;
  approximate_search_t* const forward_search_state = &archive_search->forward_search_state;
  pattern_prepare(&archive_search->sequence,&forward_search_state->pattern,
      &forward_search_state->region_profile,forward_search_state->as_parameters,
      prepare_rl_pattern,&forward_search_state->do_quality_search,archive_search->mm_stack);
  if (archive_search->emulate_rc_search) {
    approximate_search_t* const reverse_search_state = &archive_search->reverse_search_state;
    pattern_prepare(&archive_search->rc_sequence,&reverse_search_state->pattern,
        &reverse_search_state->region_profile,reverse_search_state->as_parameters,
        prepare_rl_pattern,&reverse_search_state->do_quality_search,archive_search->mm_stack);
  } else {
    pattern_clear(&archive_search->reverse_search_state.pattern);
  }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_PREPARE_SEQUENCE,PROFILE_LEVEL);
}
GEM_INLINE void archive_search_reset(archive_search_t* const archive_search) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_INIT,PROFILE_LEVEL);
  // Instantiate parameters actual-values
  const uint64_t sequence_length = sequence_get_length(&archive_search->sequence);
  search_instantiate_values(&archive_search->as_parameters,sequence_length);
  // Prepare for sequence
  archive_search_prepare_sequence(archive_search);
  // Clear F/R search states
  approximate_search_reset(&archive_search->forward_search_state);
  if (archive_search->emulate_rc_search) {
    approximate_search_reset(&archive_search->reverse_search_state);
  }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_INIT,PROFILE_LEVEL);
}
GEM_INLINE void archive_search_delete(archive_search_t* const archive_search) {
  // Delete Sequence
  sequence_destroy(&archive_search->sequence);
  sequence_destroy(&archive_search->rc_sequence);
  // Destroy search states
  approximate_search_destroy(&archive_search->forward_search_state);
  approximate_search_destroy(&archive_search->reverse_search_state);
  // Free handler
  mm_free(archive_search);
}
/*
 * Memory Injection (Support Data Structures)
 */
void archive_search_inject_mm_stack(
    archive_search_t* const archive_search,mm_stack_t* const mm_stack) {
  archive_search->mm_stack = mm_stack;
  approximate_search_inject_mm_stack(&archive_search->forward_search_state,mm_stack);
  approximate_search_inject_mm_stack(&archive_search->reverse_search_state,mm_stack);
}
void archive_search_inject_mapper_stats(
    archive_search_t* const archive_search,mapper_stats_t* mapper_stats) {
  archive_search->mapper_stats = mapper_stats;
}
void archive_search_inject_interval_set(
    archive_search_t* const archive_search,interval_set_t* const interval_set) {
  approximate_search_inject_interval_set(&archive_search->forward_search_state,interval_set);
  approximate_search_inject_interval_set(&archive_search->reverse_search_state,interval_set);
}
void archive_search_inject_text_collection(
    archive_search_t* const archive_search,text_collection_t* const text_collection) {
  archive_search->text_collection = text_collection;
  approximate_search_inject_text_collection(&archive_search->forward_search_state,text_collection);
  approximate_search_inject_text_collection(&archive_search->reverse_search_state,text_collection);
}
void archive_search_inject_filtering_candidates(
    archive_search_t* const archive_search,
    filtering_candidates_t* const filtering_candidates_forward,
    filtering_candidates_t* const filtering_candidates_reverse) {
  approximate_search_inject_filtering_candidates(
      &archive_search->forward_search_state,filtering_candidates_forward);
  approximate_search_inject_filtering_candidates(
      &archive_search->reverse_search_state,filtering_candidates_reverse);
}
/*
 * Accessors
 */
GEM_INLINE sequence_t* archive_search_get_sequence(const archive_search_t* const archive_search) {
  return (sequence_t*)&archive_search->sequence;
}
GEM_INLINE bool archive_search_finished(const archive_search_t* const archive_search) {
  if (archive_search->archive->indexed_complement) {
    return archive_search->forward_search_state.search_stage == asearch_stage_end;
  } else {
    return archive_search->forward_search_state.search_stage == asearch_stage_end &&
           archive_search->reverse_search_state.search_stage == asearch_stage_end;
  }
}
GEM_INLINE uint64_t archive_search_get_search_exact_matches(const archive_search_t* const archive_search) {
  if (archive_search->archive->indexed_complement) {
    return approximate_search_get_num_exact_filtering_candidates(&archive_search->forward_search_state);
  } else {
    return approximate_search_get_num_exact_filtering_candidates(&archive_search->forward_search_state) +
           approximate_search_get_num_exact_filtering_candidates(&archive_search->reverse_search_state);
  }
}
GEM_INLINE uint64_t archive_search_get_max_region_length(const archive_search_t* const archive_search) {
  if (archive_search->archive->indexed_complement) {
    return archive_search->forward_search_state.region_profile.max_region_length;
  } else {
    return MAX(archive_search->forward_search_state.region_profile.max_region_length,
               archive_search->reverse_search_state.region_profile.max_region_length);
  }
}
GEM_INLINE uint64_t archive_search_get_num_zero_regions(const archive_search_t* const archive_search) {
  if (archive_search->archive->indexed_complement) {
    return archive_search->forward_search_state.region_profile.num_zero_regions;
  } else {
    return MAX(archive_search->forward_search_state.region_profile.num_zero_regions,
               archive_search->reverse_search_state.region_profile.num_zero_regions);
  }
}
GEM_INLINE uint64_t archive_search_get_num_regions_profile(const archive_search_t* const archive_search) {
  uint64_t num_regions_profile = approximate_search_get_num_regions_profile(&archive_search->forward_search_state);
  if (!archive_search->archive->indexed_complement) {
    num_regions_profile += approximate_search_get_num_regions_profile(&archive_search->reverse_search_state);
  }
  return num_regions_profile;
}
GEM_INLINE uint64_t archive_search_get_num_decode_candidates(const archive_search_t* const archive_search) {
  uint64_t num_decode_candidates = approximate_search_get_num_decode_candidates(&archive_search->forward_search_state);
  if (!archive_search->archive->indexed_complement) {
    num_decode_candidates += approximate_search_get_num_decode_candidates(&archive_search->reverse_search_state);
  }
  return num_decode_candidates;
}
GEM_INLINE uint64_t archive_search_get_num_verify_candidates(const archive_search_t* const archive_search) {
  uint64_t num_verify_candidates = approximate_search_get_num_verify_candidates(&archive_search->forward_search_state);
  if (!archive_search->archive->indexed_complement) {
    num_verify_candidates += approximate_search_get_num_verify_candidates(&archive_search->reverse_search_state);
  }
  return num_verify_candidates;
}
/*
 * Utils
 */
GEM_INLINE void archive_search_hold_verification_candidates(archive_search_t* const archive_search) {
  approximate_search_hold_verification_candidates(&archive_search->forward_search_state);
  if (archive_search->emulate_rc_search) {
    approximate_search_hold_verification_candidates(&archive_search->reverse_search_state);
  }
}
GEM_INLINE void archive_search_release_verification_candidates(archive_search_t* const archive_search) {
  approximate_search_release_verification_candidates(&archive_search->forward_search_state);
  if (archive_search->emulate_rc_search) {
    approximate_search_release_verification_candidates(&archive_search->reverse_search_state);
  }
}
