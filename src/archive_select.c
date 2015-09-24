/*
 * PROJECT: GEMMapper
 * FILE: archive_select.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_select.h"
#include "archive_score.h"
#include "archive_text.h"
#include "match_align.h"
#include "match_align_dto.h"
#include "matches_classify.h"
#include "output_map.h"

/*
 * Realigning Matches
 */
GEM_INLINE void archive_select_realign_match_interval(
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
GEM_INLINE void archive_select_process_trace_matches(
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
GEM_INLINE void archive_select_process_interval_matches(
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
    match_trace.match_alignment.match_position = fm_index_lookup(archive->fm_index,match_interval->lo);
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
      match_trace.match_alignment.match_position = fm_index_lookup(archive->fm_index,sa_position);
      // (Re)Align the match [Already DONE] & Add the match
      matches_add_match_trace(matches,&match_trace,false,archive->locator,archive_search->mm_stack);
    }
  }
}
GEM_INLINE void archive_select_decode_matches(
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
GEM_INLINE void archive_select_matches(
    archive_search_t* const archive_search,
    const bool paired_mapping,matches_t* const matches) {
  PROF_START(GP_ARCHIVE_SELECT_SE_MATCHES);
  // Instantiate Search Parameters Values
  select_parameters_t* const select_parameters = archive_search->select_parameters;
  select_instantiate_values(select_parameters,sequence_get_length(&archive_search->sequence));
  // Calculate the number of matches to decode wrt input parameters
  uint64_t reported_strata = 0, last_stratum_reported_matches = 0;
  matches_counters_compute_matches_to_decode(
      matches->counters,select_parameters->min_reported_strata_nominal,select_parameters->min_reported_matches,
      select_parameters->max_reported_matches,&reported_strata,&last_stratum_reported_matches);
  // Decode matches
  if (reported_strata > 0) {
    archive_select_decode_matches(archive_search,matches,reported_strata,last_stratum_reported_matches);
  } else {
    // Remove all matches
    matches_get_clear_match_traces(matches);
    PROF_STOP(GP_ARCHIVE_SELECT_SE_MATCHES);
    return;
  }
  // Select alignment-Model and process accordingly
  archive_score_matches_se(archive_search,paired_mapping,matches);
  // Check matches
  if (select_parameters->check_correct || select_parameters->check_optimum || select_parameters->check_complete) {
    search_parameters_t* const search_parameters = archive_search->as_parameters.search_parameters;
    archive_check_matches(
        archive_search->archive,search_parameters->alignment_model,
        &search_parameters->swg_penalties,&archive_search->sequence,
        matches,select_parameters->check_optimum,select_parameters->check_complete,
        archive_search->mm_stack);
  }
  PROF_STOP(GP_ARCHIVE_SELECT_SE_MATCHES);
}
GEM_INLINE void archive_select_paired_matches(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Update stats (Check number of paired-matches)
  const uint64_t num_matches = paired_matches_get_num_maps(paired_matches);
  if (num_matches==0) return;
  PROF_START(GP_ARCHIVE_SELECT_PE_MATCHES);
  // Sample unique
  if (num_matches==1) {
    const paired_map_t* const paired_map = paired_matches_get_maps(paired_matches);
    if (paired_map->pair_relation==pair_relation_concordant) {
      mapper_stats_template_length_sample(archive_search_end1->mapper_stats,paired_map->template_length);
    }
  }
  // Select alignment-Model and process accordingly
  archive_score_matches_pe(archive_search_end1,archive_search_end2,paired_matches);
  // Discard surplus
  select_parameters_t* const select_parameters = archive_search_end1->select_parameters;
  const uint64_t num_paired_matches = vector_get_used(paired_matches->paired_maps);
  if (num_paired_matches > select_parameters->max_reported_matches) {
    vector_set_used(paired_matches->paired_maps,select_parameters->max_reported_matches);
  }
  // Check matches
  if (select_parameters->check_correct || select_parameters->check_optimum || select_parameters->check_complete) {
    search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
    archive_check_paired_matches(
        archive_search_end1->archive,search_parameters->alignment_model,
        &search_parameters->swg_penalties,&archive_search_end1->sequence,
        &archive_search_end2->sequence,paired_matches,select_parameters->check_optimum,
        select_parameters->check_complete,archive_search_end1->mm_stack);
  }
  PROF_STOP(GP_ARCHIVE_SELECT_PE_MATCHES);
}
/*
 * Check Matches
 */
GEM_INLINE void archive_check_matches(
    archive_t* const archive,const alignment_model_t alignment_model,
    swg_penalties_t* swg_penalties,sequence_t* const sequence,
    matches_t* const matches,const bool check_optimum_alignment,
    const bool check_complete,mm_stack_t* const mm_stack) {
  PROF_INC_COUNTER(GP_CHECK_NUM_READS);
  // Prepare key(s)
  mm_stack_push_state(mm_stack);
  FILE* const stream = gem_error_get_stream();
  const char* const sequence_buffer = string_get_buffer(&sequence->read);
  const uint64_t key_length = sequence_get_length(sequence);
  uint8_t* const key = mm_stack_calloc(mm_stack,key_length,uint8_t,false);
  uint64_t i;
  for (i=0;i<key_length;++i) key[i] = dna_encode(sequence_buffer[i]);
  // Traverse all matches and check their CIGAR
  uint8_t* text;
  match_trace_t opt_match_trace;
  uint64_t index_begin_position, index_end_position, text_length, text_offset, match_effective_length;
  bool correct_alignment = true, best_alignment = true;
  PROF_ADD_COUNTER(GP_CHECK_NUM_MAPS,matches_get_num_match_traces(matches));
  VECTOR_ITERATE(matches->position_matches,match_trace,match_trace_num,match_trace_t) {
    // Retrieve location
    locator_t* const locator = archive->locator;
    const uint8_t* const tag = (uint8_t*)match_trace->sequence_name;
    const uint64_t index_position = locator_inverse_map(locator,tag,match_trace->strand,match_trace->text_position);
    // Expand text-range
    locator_interval_t* const locator_interval = locator_lookup_interval(archive->locator,index_position);
    match_effective_length = match_trace->match_alignment.effective_length;
    const uint64_t boundary_error = (uint64_t)((double)match_effective_length*0.20);
    const uint64_t extra_length = match_effective_length+boundary_error;
    index_begin_position = BOUNDED_SUBTRACTION(index_position,boundary_error,locator_interval->begin_position);
    index_end_position = BOUNDED_ADDITION(index_position,extra_length,locator_interval->end_position);
    text_length = index_end_position-index_begin_position;
    // Retrieve text
    uint8_t* reverse_text;
    uint8_t* forward_text = dna_text_retrieve_sequence(
        archive->text->enc_text,index_begin_position,text_length,mm_stack);
    if (match_trace->strand==Forward) {
      text = forward_text;
      text_offset = index_position-index_begin_position;
    } else { // Reverse
      reverse_text = mm_stack_calloc(mm_stack,text_length,uint8_t,false);
      for (i=0;i<text_length;++i) reverse_text[text_length-i-1] = dna_encoded_complement(forward_text[i]);
      text = reverse_text;
      text_offset = BOUNDED_SUBTRACTION(index_end_position,match_effective_length+index_position,0);
    }
    // Check correctness
    correct_alignment = align_check_match(stream,key,key_length,text+text_offset,match_effective_length,
        matches->cigar_vector,match_trace->match_alignment.cigar_offset,match_trace->match_alignment.cigar_length,true);
    if (!correct_alignment) {
      PROF_INC_COUNTER(GP_CHECK_INCORRECT);
      break; // FAIL
    }
    // Calculate optimum alignment
    vector_t* const cigar_vector = matches->cigar_vector;
    opt_match_trace.sequence_name = match_trace->sequence_name;
    opt_match_trace.strand = match_trace->strand;
    opt_match_trace.match_alignment.match_position = index_begin_position;
    opt_match_trace.match_alignment.cigar_offset = vector_get_used(cigar_vector);
    if (check_optimum_alignment) {
      switch (alignment_model) {
        case alignment_model_none: break;
        case alignment_model_hamming: break;
        case alignment_model_levenshtein:
          GEM_NOT_IMPLEMENTED();
          break;
        case alignment_model_gap_affine: {
          match_align_input_t match_align_input = {
              .key = key,
              .key_length = key_length,
              .text = text,
              .text_length = text_length,
          };
          match_align_parameters_t match_align_parameters = {
              .swg_penalties = swg_penalties,
          };
          swg_align_match_base(&match_align_input,&match_align_parameters,
              &opt_match_trace.match_alignment,cigar_vector,mm_stack);
          opt_match_trace.swg_score = opt_match_trace.match_alignment.score;
          best_alignment = (opt_match_trace.swg_score == match_trace->swg_score);
          }
          break;
        default:
          break;
      }
      // Check result
      if(!best_alignment) {
        PROF_INC_COUNTER(GP_CHECK_SUBOPTIMAL);
        if (match_trace_num > 0) PROF_INC_COUNTER(GP_CHECK_SUBOPTIMAL_SUBDOMINANT);
        PROF_ADD_COUNTER(GP_CHECK_SUBOPTIMAL_DIFF,ABS(match_trace->swg_score-opt_match_trace.swg_score));
        PROF_ADD_COUNTER(GP_CHECK_SUBOPTIMAL_SCORE,match_trace->swg_score);
        PROF_ADD_COUNTER(GP_CHECK_SUBOPTIMAL_DISTANCE,match_trace->distance);
        break;
      }
    }
  }
  // Print Read
  if (!correct_alignment || !best_alignment) {
    fprintf(stream,"Match check failed. Global Match alignment => ");
    output_map_alignment_pretty(stream,match_trace,matches,key,key_length,text+text_offset,match_effective_length,mm_stack);
#ifdef GEM_DEBUG
    match_scaffold_t* const match_scaffold = (match_scaffold_t*) match_trace->match_scaffold;
    if (match_scaffold->match_alignment.cigar_length > 0) {
      fprintf(stream,"|> Supporting.Edit.Alignment = ");
      match_cigar_print(stream,matches->cigar_vector,match_scaffold->match_alignment.cigar_offset,
          match_scaffold->match_alignment.cigar_length);
      fprintf(stream,"\n");
    }
    fprintf(stream,"|> Supporting.Scaffold\n");
    match_scaffold_print(stream,matches,match_trace->match_scaffold);
#endif
    if (check_optimum_alignment && !best_alignment) {
      fprintf(stream,"Found better alignment Score(match-found=%d,best-alg=%d) => ",
          match_trace->swg_score,opt_match_trace.swg_score);
      const uint64_t opt_alignment_offset = opt_match_trace.match_alignment.match_position - index_begin_position;
      if (match_trace->strand==Reverse) {
        opt_match_trace.match_alignment.match_position =
            index_begin_position + (text_length - opt_alignment_offset - opt_match_trace.match_alignment.effective_length);
      }
      location_t location;
      locator_map(archive->locator,opt_match_trace.match_alignment.match_position,&location);
      opt_match_trace.text_position = location.position;
      opt_match_trace.sequence_name = location.tag;
      output_map_alignment_pretty(stream,&opt_match_trace,matches,key,key_length,
          text+opt_alignment_offset,opt_match_trace.match_alignment.effective_length,mm_stack);
    }
    fprintf(stream,"|> Read\n");
    if (sequence_has_qualities(sequence)) {
      fprintf(stream,"@%s\n",sequence_get_tag(sequence));
      fprintf(stream,"%s\n+\n",sequence_get_read(sequence));
      fprintf(stream,"%s\n",sequence_get_qualities(sequence));
    } else {
      fprintf(stream,">%s\n",sequence_get_tag(sequence));
      fprintf(stream,"%s\n",sequence_get_read(sequence));
    }
    fflush(stream); // exit(1);
  }
  mm_stack_pop_state(mm_stack,false);
}
GEM_INLINE void archive_check_paired_matches(
    archive_t* const archive,const alignment_model_t alignment_model,
    swg_penalties_t* swg_penalties,sequence_t* const sequence_end1,
    sequence_t* const sequence_end2,paired_matches_t* const paired_matches,
    const bool check_optimum_alignment,const bool check_complete,mm_stack_t* const mm_stack) {
  // Check individually end1
  archive_check_matches(archive,alignment_model,swg_penalties,sequence_end1,
      paired_matches->matches_end1,check_optimum_alignment,check_complete,mm_stack);
  // Check individually end2
  archive_check_matches(archive,alignment_model,swg_penalties,sequence_end2,
      paired_matches->matches_end2,check_optimum_alignment,check_complete,mm_stack);
  // TODO More checks related with template-size orientation, etc (But this might be stats better)
}

