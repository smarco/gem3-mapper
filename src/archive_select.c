/*
 * PROJECT: GEMMapper
 * FILE: archive_select.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_select.h"
#include "archive_score.h"
#include "archive_check.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Realigning Matches
 */
void archive_select_realign_match_interval(
    archive_search_t* const archive_search,matches_t* const matches,
    match_interval_t* const match_interval,match_trace_t* const match_trace,mm_stack_t* const mm_stack) {
  as_parameters_t* const as_parameters = &archive_search->as_parameters;
  search_parameters_t* const search_parameters = archive_search->as_parameters.search_parameters;
  const alignment_model_t alignment_model = search_parameters->alignment_model;
  match_align_input_t align_input;
  match_align_parameters_t align_parameters;
  match_scaffold_t match_scaffold = { .scaffold_regions = 0, .num_scaffold_regions = 0 };
  if (match_interval->distance==0 || alignment_model==alignment_model_none) {
    align_input.text_trace_offset = match_trace->trace_offset;
    align_input.text = match_interval->text;
    align_input.text_offset_begin = 0;
    align_input.text_offset_end = match_interval->text_length;
    align_input.key_length = sequence_get_length(&archive_search->sequence);
    align_input.text_position = match_trace->match_alignment.match_position;
    align_parameters.emulated_rc_search = match_interval->emulated_rc_search;
    align_parameters.max_error = match_interval->distance;
    align_parameters.swg_penalties = &search_parameters->swg_penalties;
    align_parameters.swg_threshold = as_parameters->swg_threshold_nominal;
    match_align_exact(matches,match_trace,&align_input,&align_parameters);
  } else {
    // Search-state (strand-based)
    approximate_search_t* const search_state = (match_interval->emulated_rc_search) ?
        &archive_search->reverse_search_state : &archive_search->forward_search_state;
    const bool left_gap_alignment = !match_interval->emulated_rc_search;
    // Align
    switch (alignment_model) {
      case alignment_model_hamming:
        align_input.key = search_state->pattern.key;
        align_input.key_length = search_state->pattern.key_length;
        align_input.text_trace_offset = match_trace->trace_offset;
        align_input.text_position = match_trace->match_alignment.match_position;
        align_input.text = match_interval->text;
        align_input.text_offset_begin = 0;
        align_input.text_offset_end = search_state->pattern.key_length;
        align_parameters.emulated_rc_search = match_interval->emulated_rc_search;
        align_parameters.allowed_enc = search_parameters->allowed_enc;
        match_align_hamming(matches,match_trace,&align_input,&align_parameters);
        break;
      case alignment_model_levenshtein:
        align_input.key = search_state->pattern.key;
        align_input.key_length = search_state->pattern.key_length,
        align_input.bpm_pattern =  &search_state->pattern.bpm_pattern;
        align_input.text_trace_offset = match_trace->trace_offset;
        align_input.text_position = match_trace->match_alignment.match_position;
        align_input.text = match_interval->text;
        align_input.text_offset_begin = 0;
        align_input.text_offset_end = match_interval->text_length;
        align_input.text_length = match_interval->text_length;
        align_parameters.emulated_rc_search = match_interval->emulated_rc_search;
        align_parameters.max_error = match_interval->distance;
        align_parameters.left_gap_alignment = left_gap_alignment;
        align_parameters.swg_penalties = &search_parameters->swg_penalties;
        match_align_levenshtein(matches,match_trace,&align_input,&align_parameters,mm_stack);
        break;
      case alignment_model_gap_affine:
        align_input.key = search_state->pattern.key;
        align_input.key_length = search_state->pattern.key_length;
        align_input.bpm_pattern =  &search_state->pattern.bpm_pattern;
        align_input.text_trace_offset = match_trace->trace_offset;
        align_input.text_position = match_trace->match_alignment.match_position;
        align_input.text = match_interval->text;
        align_input.text_length = match_interval->text_length;
        align_input.text_offset_begin = 0;
        align_input.text_offset_end = match_interval->text_length;
        align_input.align_distance_bound = match_interval->distance;
        align_input.align_match_begin_column = 0;
        align_input.align_match_end_column = match_interval->text_length;
        align_parameters.emulated_rc_search = match_interval->emulated_rc_search;
        align_parameters.max_error = match_interval->distance;
        align_parameters.max_bandwidth = match_interval->distance;
        align_parameters.left_gap_alignment = left_gap_alignment;
        align_parameters.scaffolding = search_parameters->alignment_scaffolding,
        align_parameters.scaffolding_min_coverage = as_parameters->alignment_scaffolding_min_coverage_nominal;
        align_parameters.scaffolding_matching_min_length = as_parameters->alignment_scaffolding_min_matching_length_nominal;
        align_parameters.scaffolding_homopolymer_min_context = as_parameters->alignment_scaffolding_homopolymer_min_context_nominal;
        align_parameters.allowed_enc = search_parameters->allowed_enc;
        align_parameters.swg_penalties = &search_parameters->swg_penalties;
        align_parameters.cigar_curation = search_parameters->cigar_curation;
        align_parameters.cigar_curation_min_end_context = as_parameters->cigar_curation_min_end_context_nominal;
        // Smith-Waterman-Gotoh Alignment (Gap-affine)
        match_align_smith_waterman_gotoh(matches,match_trace,&align_input,&align_parameters,&match_scaffold,mm_stack);
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
/*
 * Decoding Matches (Retrieving & Processing matches)
 */
void archive_select_process_trace_matches(
    archive_search_t* const archive_search,matches_t* const matches,
    const uint64_t reported_strata,uint64_t* const last_stratum_reported_matches) {
  const uint64_t last_stratum_distance = reported_strata-1;
  // Count already decoded matches & discard unwanted matches
  uint64_t num_matches_last_stratum = 0;
  match_trace_t* position_matches_it = matches_get_match_trace_buffer(matches);
  VECTOR_ITERATE(matches->position_matches,match_trace,match_trace_num,match_trace_t) {
    if (match_trace->distance <= last_stratum_distance) {
      // Count matches last stratum
      if (match_trace->distance == last_stratum_distance) {
        if (num_matches_last_stratum >= *last_stratum_reported_matches) {
          continue; // Too many matches decoded in the last stratum
        }
        ++num_matches_last_stratum;
      }
      // Add the match (Store it in the vector, removing unwanted ones)
      if (position_matches_it != match_trace) *position_matches_it = *match_trace;
      ++position_matches_it;
    }
  }
  const uint64_t num_initial_matches = vector_get_used(matches->position_matches);
  vector_update_used(matches->position_matches,position_matches_it); // Update used
  if (num_initial_matches != vector_get_used(matches->position_matches)) {
    // Because the matches had been reallocated, the indexed-positions are no longer valid
    matches_index_rebuild(matches,archive_search->mm_stack);
    matches_recompute_metrics(matches);
  }
  // Return matches to decode in the last stratum
  *last_stratum_reported_matches -= num_matches_last_stratum;
}
void archive_select_process_interval_matches(
    archive_search_t* const archive_search,matches_t* const matches,
    const uint64_t reported_strata,const uint64_t last_stratum_reported_matches) {
  // Parameters
  const archive_t* const archive = archive_search->archive;
  text_collection_t* const text_collection = archive_search->text_collection;
  const uint64_t last_stratum_distance = reported_strata-1;
  // Traverse all interval-matches
  uint64_t num_matches_last_stratum = 0;
  VECTOR_ITERATE(matches->interval_matches,match_interval,match_interval_num,match_interval_t) {
    if (num_matches_last_stratum >= last_stratum_reported_matches) break;
    // Get next match-interval
    if ((match_interval->lo >= match_interval->hi) || (match_interval->distance > last_stratum_distance)) continue;
    // Decode position (Use the first position of the interval, lo)
    match_trace_t match_trace;
    match_trace.match_alignment.match_position = fm_index_decode(archive->fm_index,match_interval->lo);
    match_trace.match_alignment.score = match_interval->distance;
    // (Re)Align match (retrieved with the seq-read(pattern))
    if (match_interval->distance > 0) {
      match_trace.trace_offset = archive_text_retrieve(archive->text,text_collection,
          match_trace.match_alignment.match_position,match_interval->text_length, /* + delta TODO */
          false,archive_search->mm_stack);
      // Set interval text
      const text_trace_t* const text_trace = text_collection_get_trace(text_collection,match_trace.trace_offset);
      match_interval->text = text_trace->text;
    }
    archive_select_realign_match_interval(archive_search,matches,match_interval,&match_trace,archive_search->mm_stack);
    // Add the match
    matches_add_match_trace(matches,&match_trace,false,archive->locator,archive_search->mm_stack);
    const bool last_stratum_match = (match_interval->distance == last_stratum_distance);
    if (last_stratum_match) ++num_matches_last_stratum;
    // Build the rest of the interval
    uint64_t sa_position;
    for (sa_position=match_interval->lo+1;sa_position<match_interval->hi;++sa_position) {
      // Check number of matches (only last stratum)
      if (last_stratum_match && (num_matches_last_stratum++ >= last_stratum_reported_matches)) break;
      // Decode position
      match_trace.match_alignment.match_position = fm_index_decode(archive->fm_index,sa_position);
      // (Re)Align the match [Already DONE] & Add the match
      matches_add_match_trace(matches,&match_trace,false,archive->locator,archive_search->mm_stack);
    }
  }
}
void archive_select_decode_matches(
    archive_search_t* const archive_search,matches_t* const matches,
    const uint64_t reported_strata,uint64_t last_stratum_reported_matches) {
  // Decode trace-matches. Count already decoded matches & discard unwanted matches
  archive_select_process_trace_matches(
      archive_search,matches,reported_strata,&last_stratum_reported_matches);
  // Decode interval-matches
  archive_select_process_interval_matches(
      archive_search,matches,reported_strata,last_stratum_reported_matches);
}
/*
 * Select Paired-Matches
 */
void archive_select_se_matches(
    archive_search_t* const archive_search,
    const bool paired_mapping,matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SELECT_SE_MATCHES,PROFILE_LEVEL);
  // Instantiate Search Parameters Values
  select_parameters_t* const select_parameters = archive_search->select_parameters;
  select_parameters_instantiate_values(select_parameters,sequence_get_length(&archive_search->sequence));
  // Calculate the number of matches to decode wrt input parameters
  uint64_t reported_strata = 0, last_stratum_reported_matches = 0;
  matches_counters_compute_matches_to_decode(
      matches->counters,select_parameters->min_reported_strata_nominal,select_parameters->min_reported_matches,
      select_parameters->max_reported_matches,&reported_strata,&last_stratum_reported_matches);
  if (reported_strata==0) {
    // Remove all matches
    matches_get_clear_match_traces(matches);
    PROFILE_STOP(GP_ARCHIVE_SELECT_SE_MATCHES,PROFILE_LEVEL);
    return;
  }
  // Decode matches
  archive_select_decode_matches(archive_search,matches,reported_strata,last_stratum_reported_matches);
  PROFILE_STOP(GP_ARCHIVE_SELECT_SE_MATCHES,PROFILE_LEVEL);
}
void archive_select_pe_matches(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Update stats (Check number of paired-matches)
  const uint64_t num_matches = paired_matches_get_num_maps(paired_matches);
  if (num_matches==0) return;
  PROFILE_START(GP_ARCHIVE_SELECT_PE_MATCHES,PROFILE_LEVEL);
  // Sample unique
  if (num_matches==1) {
    const paired_map_t* const paired_map = paired_matches_get_maps(paired_matches);
    if (paired_map->pair_relation==pair_relation_concordant) {
      mapper_stats_template_length_sample(archive_search_end1->mapper_stats,paired_map->template_length);
    }
  }
  // Discard surplus
  select_parameters_t* const select_parameters = archive_search_end1->select_parameters;
  const uint64_t num_paired_matches = vector_get_used(paired_matches->paired_maps);
  if (num_paired_matches > select_parameters->max_reported_matches) {
    vector_set_used(paired_matches->paired_maps,select_parameters->max_reported_matches);
  }
  PROFILE_STOP(GP_ARCHIVE_SELECT_PE_MATCHES,PROFILE_LEVEL);
}
