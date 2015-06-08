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
    align_input.key_length = sequence_get_length(&archive_search->sequence);
    align_input.text_trace_offset = match_trace->trace_offset;
    align_input.text_position = match_trace->match_alignment.match_position;
    align_input.text = match_interval->text;
    align_input.text_offset_begin = 0;
    align_input.text_offset_end = match_interval->length;
    align_parameters.emulated_rc_search = match_interval->emulated_rc_search;
    align_parameters.max_error = match_interval->distance;
    align_parameters.swg_penalties = &search_parameters->swg_penalties;
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
        align_parameters.emulated_rc_search = match_interval->emulated_rc_search;
        align_parameters.allowed_enc = search_parameters->allowed_enc;
        match_align_hamming(matches,match_trace,&align_input,&align_parameters);
        break;
      case alignment_model_levenshtein:
        align_input.key = search_state->pattern.key;
        align_input.bpm_pattern =  &search_state->pattern.bpm_pattern;
        align_input.text_trace_offset = match_trace->trace_offset;
        align_input.text_position = match_trace->match_alignment.match_position;
        align_input.text = match_interval->text;
        align_input.text_offset_begin = 0;
        align_input.text_length = match_interval->length;
        align_parameters.emulated_rc_search = match_interval->emulated_rc_search;
        align_parameters.max_error = match_interval->distance;
        align_parameters.left_gap_alignment = left_gap_alignment;
        match_align_levenshtein(matches,match_trace,&align_input,&align_parameters,mm_stack);
        break;
      case alignment_model_gap_affine:
        align_input.key = search_state->pattern.key;
        align_input.key_length = search_state->pattern.key_length;
        align_input.bpm_pattern =  &search_state->pattern.bpm_pattern;
        align_input.text_trace_offset = match_trace->trace_offset;
        align_input.text_position = match_trace->match_alignment.match_position;
        align_input.text = match_interval->text;
        align_input.text_length = match_interval->length;
        align_input.text_offset_begin = 0;
        align_input.text_offset_end = match_interval->length;
        align_parameters.emulated_rc_search = match_interval->emulated_rc_search;
        align_parameters.max_error = match_interval->distance;
        align_parameters.max_bandwidth = match_interval->distance;
        align_parameters.left_gap_alignment = left_gap_alignment;
        align_parameters.min_coverage = as_parameters->region_scaffolding_coverage_threshold_nominal,
        align_parameters.min_matching_length = as_parameters->region_scaffolding_min_length_nominal,
        align_parameters.min_context_length = as_parameters->region_scaffolding_min_context_length_nominal,
        align_parameters.allowed_enc = search_parameters->allowed_enc;
        align_parameters.swg_penalties = &search_parameters->swg_penalties;
        // Scaffold the alignment
        if (search_parameters->allow_region_chaining) {
          match_scaffold_alignment(matches,&align_input,&align_parameters,&match_scaffold,mm_stack);
        }
        // Smith-Waterman-Gotoh Alignment (Gap-affine)
        match_align_smith_waterman_gotoh(matches,match_trace,&align_input,&align_parameters,&match_scaffold,mm_stack);
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
GEM_INLINE void archive_select_correct_cigar(
    const archive_t* const archive,matches_t* const matches,
    const uint64_t cigar_buffer_offset,const uint64_t cigar_length) {
  if (archive->filter_type==Iupac_dna) {
    matches_cigar_reverse(matches,cigar_buffer_offset,cigar_length);
  } else { // Iupac_colorspace_dna
    matches_cigar_reverse_colorspace(matches,cigar_buffer_offset,cigar_length);
//    /* In this case, we also have to re-encode the mismatches */
//    if (in_colorspace) { // TODO
//      for (i=0;i<mism_num;++i) {
//        misms[i].mismatch=misms_colorspace[misms[i].mismatch];
//          // const ch_t misms_colorspace[CH_T_RANGE]= { ['A']='a', ['C']='b', ['G']='c', ['T']='d', ['N']='e' };
//      }
//    }
  }
}
/*
 * Decoding Matches (Retrieving & Processing matches)
 */
GEM_INLINE void archive_select_locate_match_trace(
    const archive_t* const archive,matches_t* const matches,match_trace_t* const match_trace,
    const uint64_t seq_length,const int64_t match_length,const bool emulated_rc_search) {
  GEM_INTERNAL_CHECK(match_length >= 0,"Match effective length must be positive");
  location_t location;
  locator_map(archive->locator,match_trace->match_alignment.match_position,&location);
  match_trace->text_position = location.position;
  match_trace->sequence_name = location.tag;
  if (location.direction == Reverse) { // Adjust position by the effective length
    match_trace->text_position -= (uint64_t) match_length;
    match_trace->strand = Reverse;
    GEM_INTERNAL_CHECK(!emulated_rc_search,
        "Archive-Select. Locating match-trace. "
        "Impossible combination (search_strand==Reverse while emulated searching FR-index)");
  } else {
    match_trace->strand = emulated_rc_search ? Reverse : Forward;
  }
}
GEM_INLINE uint64_t archive_select_process_trace_matches(
    archive_search_t* const archive_search,matches_t* const matches,
    const uint64_t strata_to_decode,const uint64_t matches_to_decode_last_stratum) {
  const archive_t* const archive = archive_search->archive;
  const uint64_t seq_length = sequence_get_length(&archive_search->sequence);
  const uint64_t last_stratum_distance = strata_to_decode-1;
  /*
   * Match-Trace
   *   Count already decoded matches & discard unwanted matches
   */
  uint64_t num_matches_last_stratum = 0;
  match_trace_t* position_matches_it = matches_get_match_traces(matches);
  VECTOR_ITERATE(matches->position_matches,match_trace,match_trace_num,match_trace_t) {
    if (match_trace->distance <= last_stratum_distance) {
      // Count matches last stratum
      if (match_trace->distance == last_stratum_distance) {
        if (num_matches_last_stratum < matches_to_decode_last_stratum) {
          ++num_matches_last_stratum;
        } else {
          continue; // Too many matches decoded in the last stratum
        }
      }
      if (match_trace->text_position == UINT64_MAX) { // Not decoded yet
        // 1.- (Re)Align match [Already DONE]
        // 2.- Correct CIGAR (Reverse it if the search was performed in the reverse strand, emulated)
        if (match_trace->emulated_rc_search) {
          archive_select_correct_cigar(archive,matches,
              match_trace->match_alignment.cigar_offset,match_trace->match_alignment.cigar_length);
        }
        // 3.- Locate-map the match
        archive_select_locate_match_trace(archive,matches,match_trace,seq_length,
            match_trace->match_alignment.effective_length,match_trace->emulated_rc_search);
      }
      // 4.- Add the match (Store it in the vector, removing unwanted ones)
      if (position_matches_it != match_trace) *position_matches_it = *match_trace;
      ++position_matches_it; // TODO regenerate matches-index on the fly
    }
  }
  vector_update_used(matches->position_matches,position_matches_it); // Update used
  matches_clear_index(matches);
  // Return matches in the last stratum
  return num_matches_last_stratum;
}
GEM_INLINE void archive_select_process_interval_matches(
    archive_search_t* const archive_search,matches_t* const matches,
    const uint64_t strata_to_decode,const uint64_t matches_to_decode_last_stratum,
    uint64_t num_matches_last_stratum) {
  // Parameters
  const archive_t* const archive = archive_search->archive;
  text_collection_t* const text_collection = archive_search->text_collection;
  const uint64_t seq_length = sequence_get_length(&archive_search->sequence);
  const uint64_t last_stratum_distance = strata_to_decode-1;
  // Traverse all interval-matches
  VECTOR_ITERATE(matches->interval_matches,match_interval,match_interval_num,match_interval_t) {
    if (num_matches_last_stratum >= matches_to_decode_last_stratum) break;
    // Get next match-interval
    if ((match_interval->lo >= match_interval->hi) || (match_interval->distance > last_stratum_distance)) continue;
    // 0.- Decode position (Use the first position of the interval, lo)
    match_trace_t match_trace;
    match_trace.match_alignment.match_position = fm_index_lookup(archive->fm_index,match_interval->lo);
    match_trace.match_alignment.score = match_interval->distance;
    // 1.- (Re)Align match (retrieved with the seq-read(pattern))
    if (match_interval->distance > 0) {
      match_trace.trace_offset = archive_text_retrieve(archive->text,text_collection,
          match_trace.match_alignment.match_position,match_interval->length, /* + delta TODO */
          false,archive_search->mm_stack);
      // Set interval text
      const text_trace_t* const text_trace = text_collection_get_trace(text_collection,match_trace.trace_offset);
      match_interval->text = text_trace->text;
    }
    archive_select_realign_match_interval(archive_search,matches,match_interval,&match_trace,archive_search->mm_stack);
    // 2.- Correct CIGAR (Reverse it if the search was performed in the reverse strand, emulated)
    if (match_interval->emulated_rc_search) {
      archive_select_correct_cigar(archive,matches,
          match_trace.match_alignment.cigar_offset,match_trace.match_alignment.cigar_length);
    }
    // 3.- Locate-map the match
    archive_select_locate_match_trace(archive,matches,&match_trace,
        seq_length,match_trace.match_alignment.effective_length,match_interval->emulated_rc_search);
    // 4.- Add the match
    matches_add_match_trace_t(matches,&match_trace,false,archive_search->mm_stack);
    const bool last_stratum_match = (match_interval->distance == last_stratum_distance);
    if (last_stratum_match) ++num_matches_last_stratum;
    // 5.- Build the rest of the interval
    uint64_t sa_position;
    for (sa_position=match_interval->lo+1;sa_position<match_interval->hi;++sa_position) {
      // Check number of matches (only last stratum)
      if (last_stratum_match && (++num_matches_last_stratum >= matches_to_decode_last_stratum)) break;
      // 0.- Decode position
      match_trace.match_alignment.match_position = fm_index_lookup(archive->fm_index,sa_position);
      // 1.- (Re)Align the match [Already DONE]
      // 2.- Correct CIGAR [Already DONE]
      // 3.- Locate-map the match
      archive_select_locate_match_trace(archive,matches,&match_trace,
          seq_length,match_trace.match_alignment.effective_length,match_interval->emulated_rc_search);
      // 4.- Add the match
      matches_add_match_trace_t(matches,&match_trace,false,archive_search->mm_stack);
    }
  }
}
GEM_INLINE void archive_select_decode_matches(
    archive_search_t* const archive_search,matches_t* const matches,
    const uint64_t strata_to_decode,const uint64_t matches_to_decode_last_stratum) {
  // Decode trace-matches. Count already decoded matches & discard unwanted matches
  const uint64_t num_matches_last_stratum =
      archive_select_process_trace_matches(archive_search,matches,strata_to_decode,matches_to_decode_last_stratum);
  // Decode interval.matches (until we reach @max_decoded_matches)
  archive_select_process_interval_matches(
      archive_search,matches,strata_to_decode,matches_to_decode_last_stratum,num_matches_last_stratum);
}
GEM_INLINE void archive_select_decode_trace_matches(
    archive_search_t* const archive_search,matches_t* const matches,
    match_trace_t* const match_trace,const uint64_t num_match_traces,
    const bool force_strand,const strand_t strand) {
  uint64_t mate_pos;
  for (mate_pos=0;mate_pos<num_match_traces;++mate_pos) {
    match_trace_t* const mate_trace = match_trace + mate_pos;
    if (mate_trace->text_position != UINT64_MAX) continue; // Already decoded
//    if (correct_cigar) {
//      // Correct CIGAR
//      archive_select_correct_cigar(archive_search->archive,matches,
//          match_trace->strand,match_trace->cigar_buffer_offset,match_trace->cigar_length);
//    }
    // Locate-map the match
    archive_select_locate_match_trace(archive_search->archive,matches,
        mate_trace,archive_search->forward_search_state.pattern.key_length,
        mate_trace->match_alignment.effective_length,false);
    if (force_strand) mate_trace->strand = strand; // Restore original
  }
}
GEM_INLINE void archive_select_decode_trace_matches_all(
    archive_search_t* const archive_search,matches_t* const matches,
    const bool force_strand,const strand_t strand) {
  match_trace_t* const match_trace = matches_get_match_traces(matches);
  const uint64_t num_match_traces = matches_get_num_match_traces(matches);
  archive_select_decode_trace_matches(archive_search,matches,
      match_trace,num_match_traces,force_strand,strand);
}
/*
 * Filters
 */
GEM_INLINE void archive_select_filter_matches_mapq(
    matches_t* const matches,const uint8_t mapq_threshold) {
  match_trace_t* match_trace_out = matches_get_match_traces(matches);
  VECTOR_ITERATE(matches->position_matches,match_trace_in,n,match_trace_t) {
    if (match_trace_in->mapq_score >= mapq_threshold) {
      *match_trace_out = *match_trace_in;
      ++match_trace_out;
    }
  }
  vector_update_used(matches->position_matches,match_trace_out);
}
GEM_INLINE void archive_select_filter_paired_matches_mapq(
    paired_matches_t* const paired_matches,const uint8_t mapq_threshold) {
  paired_match_t* paired_match_out = vector_get_mem(paired_matches->matches,paired_match_t);
  VECTOR_ITERATE(paired_matches->matches,paired_match_in,n,paired_match_t) {
    if (paired_match_in->mapq_score >= mapq_threshold) {
      *paired_match_out = *paired_match_in;
      ++paired_match_out;
    }
  }
  vector_update_used(paired_matches->matches,paired_match_out);
}
/*
 * Select Paired-Matches
 */
GEM_INLINE void archive_select_calculate_matches_to_decode(
    matches_t* const matches,const uint64_t max_decoded_matches,const uint64_t min_decoded_strata,
    const uint64_t min_reported_matches,const uint64_t max_reported_matches,
    uint64_t* const total_strata_to_decode,uint64_t* const matches_to_decode_last_stratum_out) {
  // Compact counters (Shrink the counters to the last non-zero stratum)
  const uint64_t max_nz_stratum = matches_counters_compact(matches); // Maximum reachable stratum w.r.t counters
  if (max_nz_stratum==0) return;
  // Maximum stratum to decode (increased by @max_decoded_matches)
  const uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
  uint64_t strata_to_decode, total_matches;
  for (strata_to_decode=0,total_matches=0;strata_to_decode<max_nz_stratum;++strata_to_decode) {
    total_matches += counters[strata_to_decode];
    if (total_matches > max_decoded_matches) {
      total_matches -= counters[strata_to_decode];
      break;
    }
  }
  // Maximum stratum to decode (increased by @min_decoded_strata)
  if (min_decoded_strata > 0) {
    const uint64_t min_nz_stratum = matches_counters_get_min_distance(matches);
    const uint64_t mandatory_strata = min_nz_stratum + min_decoded_strata;
    for (;strata_to_decode<max_nz_stratum && strata_to_decode<mandatory_strata;++strata_to_decode) {
      total_matches += counters[strata_to_decode];
    }
  }
  // Maximum stratum to decode (increased by @min_reported_matches)
  for (;strata_to_decode<max_nz_stratum && total_matches<min_reported_matches;++strata_to_decode) {
    total_matches += counters[strata_to_decode];
  }
  // Maximum stratum to decode (lowered by @max_reported_matches & @min_reported_matches)
  for (;strata_to_decode>0;--strata_to_decode) {
    const uint64_t prev_acc = total_matches - counters[strata_to_decode-1];
    if (total_matches <= max_reported_matches || prev_acc < min_reported_matches) break;
    total_matches = prev_acc;
  }
  // Decode matches
  if (total_matches!=0) { // => (strata_to_decode > 0)
    *total_strata_to_decode = strata_to_decode;
    *matches_to_decode_last_stratum_out = UINT64_MAX;
    if (total_matches > max_reported_matches) {
      const uint64_t prev_acc = total_matches - counters[strata_to_decode-1];
      const uint64_t max_matches_from_last_stratum = max_reported_matches - prev_acc;
      *matches_to_decode_last_stratum_out = max_matches_from_last_stratum;
    }
  } else {
    *total_strata_to_decode = 0;
    *matches_to_decode_last_stratum_out = 0;
  }
}
GEM_INLINE void archive_select_matches(
    archive_search_t* const archive_search,const bool curate_matches,
    const bool score_matches,const bool sort_matches,matches_t* const matches) {
  PROF_START(GP_ARCHIVE_SELECT_SE_MATCHES);
  // Instantiate Search Parameters Values
  select_parameters_t* const select_parameters = archive_search->select_parameters;
  archive_select_instantiate_values(select_parameters,sequence_get_length(&archive_search->sequence));
  // Check if we need to decode sth
  if (select_parameters->max_decoded_matches==0 &&
      select_parameters->min_decoded_strata==0 &&
      select_parameters->min_reported_matches==0) return;
  if (select_parameters->min_reported_matches==0 && select_parameters->max_reported_matches==0) return;
  // Calculate the number of matches to decode wrt input parameters
  uint64_t strata_to_decode = 0, matches_to_decode_last_stratum = 0;
  archive_select_calculate_matches_to_decode(
      matches,select_parameters->max_decoded_matches,
      select_parameters->min_decoded_strata,select_parameters->min_reported_matches,
      select_parameters->max_reported_matches,&strata_to_decode,&matches_to_decode_last_stratum);
  if (strata_to_decode > 0) {
    // Decode matches
    archive_select_decode_matches(archive_search,matches,strata_to_decode,matches_to_decode_last_stratum);
    // Score matches
    if (score_matches) archive_score_matches_se(archive_search,matches);
    // Curate matches
    if (curate_matches) matches_curate(matches,0.20);
    // Filter by score
    if (select_parameters->mapq_threshold > 0) {
      archive_select_filter_matches_mapq(matches,select_parameters->mapq_threshold);
    }
    // Sort matches
    if (sort_matches) {
      switch (archive_search->select_parameters->sorting) {
        case matches_sorting_mapq: matches_sort_by_mapq_score(matches); break;
        case matches_sorting_distance: matches_sort_by_distance(matches); break;
        default: GEM_INVALID_CASE(); break;
      }
    }
  } else {
    // Remove all matches
    matches_get_clear_match_traces(matches);
  }
  PROF_STOP(GP_ARCHIVE_SELECT_SE_MATCHES);
}
GEM_INLINE void archive_select_paired_matches(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Update stats (Check number of paired-matches)
  PROF_START(GP_ARCHIVE_SELECT_PE_MATCHES);
  const uint64_t num_matches = vector_get_used(paired_matches->matches);
  if (num_matches > 0) {
    // Sample if unique
    if (num_matches==1) {
      const paired_match_t* const paired_match = vector_get_mem(paired_matches->matches,paired_match_t);
      if (paired_match->pair_orientation==pair_orientation_concordant) {
        mapper_stats_template_length_sample(archive_search_end1->mapper_stats,paired_match->template_length);
      }
    }
    // Score matches
    archive_score_matches_pe(archive_search_end1,archive_search_end2,paired_matches);
    // Filter by score
    const select_parameters_t* const select_parameters = archive_search_end1->select_parameters;
    if (select_parameters->mapq_threshold > 0)  {
      archive_select_filter_paired_matches_mapq(paired_matches,select_parameters->mapq_threshold);
    }
    // Sort paired-matches
    switch (select_parameters->sorting) {
      case matches_sorting_mapq: paired_matches_sort_by_mapq_score(paired_matches); break;
      case matches_sorting_distance: paired_matches_sort_by_distance(paired_matches); break;
      default: GEM_INVALID_CASE(); break;
    }
    // Discard surplus
    const uint64_t num_paired_matches = vector_get_used(paired_matches->matches);
    if (num_paired_matches > select_parameters->max_reported_matches) {
      vector_set_used(paired_matches->matches,select_parameters->max_reported_matches);
    }
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
    const uint64_t index_position = inverse_locator_map(locator,tag,match_trace->strand,match_trace->text_position);
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

