/*
 * PROJECT: GEMMapper
 * FILE: archive_check.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive/archive_check.h"
#include "align/align.h"
#include "matches/matches_cigar.h"

/*
 * Check Matches
 */
void archive_check_se_match_retrieve_text(
    archive_t* const archive,
    match_trace_t* const match_trace,
    match_align_input_t* const match_align_input,
    mm_stack_t* const mm_stack) {
  // Retrieve location
  locator_t* const locator = archive->locator;
  const uint8_t* const tag = (uint8_t*)match_trace->sequence_name;
  const uint64_t index_position = locator_inverse_map_position(
      locator,tag,Forward,match_trace->bs_strand,match_trace->text_position);
  // Expand text-range
  locator_interval_t* const locator_interval = locator_lookup_interval(archive->locator,index_position);
  const uint64_t effective_length = match_trace->match_alignment.effective_length;
  const uint64_t boundary_error = (uint64_t)((double)effective_length*0.20);
  const uint64_t extra_length = effective_length+boundary_error;
  const uint64_t index_begin_pos = BOUNDED_SUBTRACTION(index_position,boundary_error,locator_interval->begin_position);
  const uint64_t index_end_pos = BOUNDED_ADDITION(index_position,extra_length,locator_interval->end_position);
  const uint64_t text_length = index_end_pos-index_begin_pos;
  match_align_input->text_length = text_length;
  match_align_input->text_position = index_begin_pos;
  // Retrieve text
  uint8_t* reverse_text;
  uint8_t* forward_text = dna_text_retrieve_sequence(archive->text->enc_text,index_begin_pos,text_length,mm_stack);
  if (match_trace->strand==Forward) {
    const uint64_t text_offset_begin = index_position-index_begin_pos;
    match_align_input->text = forward_text+text_offset_begin;
  } else { // Reverse
    uint64_t i;
    reverse_text = mm_stack_calloc(mm_stack,text_length,uint8_t,false);
    for (i=0;i<text_length;++i) {
      reverse_text[text_length-i-1] = dna_encoded_complement(forward_text[i]);
    }
    const uint64_t text_offset_begin = BOUNDED_SUBTRACTION(index_end_pos,effective_length+index_position,0);
    match_align_input->text = reverse_text+text_offset_begin;
  }
}
bool archive_check_se_match_check_optimum(
    match_trace_t* const match_trace,
    const alignment_model_t alignment_model,
    swg_penalties_t* const swg_penalties,
    match_align_input_t* const match_align_input,
    match_trace_t* const optimum_match_trace,
    vector_t* const cigar_vector,
    mm_stack_t* const mm_stack) {
  optimum_match_trace->sequence_name = match_trace->sequence_name;
  optimum_match_trace->strand = match_trace->strand;
  optimum_match_trace->match_alignment.match_position = match_align_input->text_position;
  optimum_match_trace->match_alignment.cigar_offset = vector_get_used(cigar_vector);
  switch (alignment_model) {
    case alignment_model_none: break;
    case alignment_model_hamming: break;
    case alignment_model_levenshtein:
      GEM_NOT_IMPLEMENTED();
      break;
    case alignment_model_gap_affine: {
      match_align_parameters_t match_align_parameters = {
          .swg_penalties = swg_penalties,
      };
      align_swg_base(match_align_input,&match_align_parameters,
          &optimum_match_trace->match_alignment,cigar_vector,mm_stack);
      optimum_match_trace->swg_score = optimum_match_trace->match_alignment.score;
      return (optimum_match_trace->swg_score == match_trace->swg_score);
      }
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  return false;
}
void archive_check_se_match_print(
    FILE* const stream,
    const char* const label,
    match_trace_t* const match_trace,
    matches_t* const matches,
    match_align_input_t* const match_align_input,
    sequence_t* const sequence,
    mm_stack_t* const mm_stack) {
  tab_fprintf(stream,"[GEM]>Check.SE\n");
  tab_global_inc();
  tab_fprintf(stream,"=> Match check failed (%s)\n",label);
  tab_fprintf(stream,"=> Sequence\n");
  sequence_print(stream,sequence);
  tab_fprintf(stream,"=> Alignment\n");
  tab_global_inc();
// output_map_alignment_pretty // TODO
//  (%s:%"PRIu64":%c:)\n",match_trace->sequence_name,
//      match_trace->text_position,(match_trace->strand==Forward)?'+':'-');
// match_trace_print(stream,matches->cigar_vector,);
  match_alignment_print_pretty(stream,&match_trace->match_alignment,
      matches->cigar_vector,match_align_input->key,match_align_input->key_length,
      match_align_input->text,match_trace->match_alignment.effective_length,mm_stack);
  tab_global_dec();
  match_scaffold_t* const match_scaffold = (match_scaffold_t*) match_trace->match_scaffold;
  tab_fprintf(stream,"=> Supporting.Edit.Alignment ");
  match_cigar_print(stream,matches->cigar_vector,
      match_scaffold->match_alignment.cigar_offset,
      match_scaffold->match_alignment.cigar_length);
  fprintf(stream,"\n");
  tab_fprintf(stream,"=> Supporting.Scaffold\n");
  tab_global_inc();
  match_scaffold_print(stream,matches,match_trace->match_scaffold);
  tab_global_dec();
  tab_global_dec();
  fflush(stream); // exit(1);
}
void archive_check_se_match_print_incorrect(
    FILE* const stream,
    match_trace_t* const match_trace,
    matches_t* const matches,
    match_align_input_t* const match_align_input,
    sequence_t* const sequence,
    mm_stack_t* const mm_stack) {
  archive_check_se_match_print(stream,"INCORRECT",
      match_trace,matches,match_align_input,sequence,mm_stack);
}
void archive_check_se_match_print_suboptimum(
    FILE* const stream,
    match_trace_t* const match_trace,
    matches_t* const matches,
    match_align_input_t* const match_align_input,
    sequence_t* const sequence,
    match_trace_t* const optimum_match_trace,
    locator_t* const locator,
    mm_stack_t* const mm_stack) {
  archive_check_se_match_print(stream,"SUBOPTIMUM",
      match_trace,matches,match_align_input,sequence,mm_stack);
  tab_global_inc();
  tab_fprintf(stream,"=> Suboptimum.Alignment.Score (match-found=%d,best-alg=%d)\n",
      match_trace->swg_score,optimum_match_trace->swg_score);
  const uint64_t opt_alignment_offset =
      optimum_match_trace->match_alignment.match_position - match_align_input->text_position;
  if (match_trace->strand==Reverse) {
    optimum_match_trace->match_alignment.match_position = match_align_input->text_position +
        (match_align_input->text_length -
         opt_alignment_offset -
         optimum_match_trace->match_alignment.effective_length);
  }
  location_t location;
  locator_map(locator,optimum_match_trace->match_alignment.match_position,&location);
  optimum_match_trace->text_position = location.position;
  optimum_match_trace->sequence_name = location.tag;
  tab_fprintf(stream,"=> Optimum.Alignment");
  tab_global_inc();
  match_alignment_print_pretty(
      stream,&optimum_match_trace->match_alignment,
      matches->cigar_vector,match_align_input->key,
      match_align_input->key_length,match_align_input->text+opt_alignment_offset,
      optimum_match_trace->match_alignment.effective_length,mm_stack);
  tab_global_dec();
  tab_global_dec();
}
void archive_check_se_matches(
    archive_t* const archive,
    const alignment_model_t alignment_model,
    swg_penalties_t* swg_penalties,
    sequence_t* const sequence,
    matches_t* const matches,
    const archive_check_type check_type,
    mm_stack_t* const mm_stack) {
  // Parameters
  PROF_INC_COUNTER(GP_CHECK_NUM_READS);
  mm_stack_push_state(mm_stack);
  FILE* const stream = gem_error_get_stream();
  const char* const sequence_buffer = string_get_buffer(&sequence->read);
  const uint64_t key_length = sequence_get_length(sequence);
  // Prepare
  uint8_t* const key = mm_stack_calloc(mm_stack,key_length,uint8_t,false);
  uint64_t i;
  for (i=0;i<key_length;++i) key[i] = dna_encode(sequence_buffer[i]);
  match_align_input_t match_align_input;
  match_align_input.key = key;
  match_align_input.key_length = key_length;
  // Traverse all matches and check their CIGAR
  PROF_ADD_COUNTER(GP_CHECK_NUM_MAPS,matches_get_num_match_traces(matches));
  VECTOR_ITERATE(matches->position_matches,match_trace,match_trace_num,match_trace_t) {
    archive_check_se_match_retrieve_text(archive,match_trace,&match_align_input,mm_stack);
    // Check correctness
    match_alignment_t* const match_alignment = &match_trace->match_alignment;
    const bool correct_alignment = align_check(stream,key,key_length,match_align_input.text,
        match_alignment->effective_length,matches->cigar_vector,match_alignment->cigar_offset,
        match_alignment->cigar_length,true);
    if (!correct_alignment) {
      PROF_INC_COUNTER(GP_CHECK_INCORRECT);
      archive_check_se_match_print_incorrect(stream,match_trace,matches,&match_align_input,sequence,mm_stack);
      break; // FAIL
    }
    if (check_type==archive_check_correct) continue;
    if (check_type==archive_check_correct__first_optimum) {
      if (match_trace_num > 0) continue;
    }
    // Check optimum alignment
    if (match_trace_num==0) {
      PROF_INC_COUNTER(GP_CHECK_PRIMARY_SUBOPTIMAL);
    } else {
      PROF_INC_COUNTER(GP_CHECK_SUBDOMINANT_SUBOPTIMAL);
    }
    // Compute optimal & check
    match_trace_t optimum_match_trace;
    const bool optimal_alignment =
        archive_check_se_match_check_optimum(match_trace,alignment_model,swg_penalties,
        &match_align_input,&optimum_match_trace,matches->cigar_vector,mm_stack);
    // Check result
    if (!optimal_alignment) {
      if (match_trace_num==0) {
        PROF_INC_COUNTER(GP_CHECK_SUBDOMINANT_SUBOPTIMAL_FAIL);
        PROF_ADD_COUNTER(GP_CHECK_SUBDOMINANT_SUBOPTIMAL_DIFF,ABS(match_trace->swg_score-optimum_match_trace.swg_score));
        PROF_ADD_COUNTER(GP_CHECK_SUBDOMINANT_SUBOPTIMAL_SCORE,match_trace->swg_score);
        PROF_ADD_COUNTER(GP_CHECK_SUBDOMINANT_SUBOPTIMAL_DISTANCE,match_trace->distance);
      } else {
        PROF_INC_COUNTER(GP_CHECK_PRIMARY_SUBOPTIMAL_FAIL);
        PROF_ADD_COUNTER(GP_CHECK_PRIMARY_SUBOPTIMAL_SCORE,match_trace->swg_score);
        PROF_ADD_COUNTER(GP_CHECK_PRIMARY_SUBOPTIMAL_DIFF,ABS(match_trace->swg_score-optimum_match_trace.swg_score));
        PROF_ADD_COUNTER(GP_CHECK_PRIMARY_SUBOPTIMAL_DISTANCE,match_trace->distance);
      }
      archive_check_se_match_print_suboptimum(stream,match_trace,matches,
          &match_align_input,sequence,&optimum_match_trace,archive->locator,mm_stack);
      break; // FAIL
    }
  }
  mm_stack_pop_state(mm_stack);
}
void archive_check_pe_matches(
    archive_t* const archive,
    const alignment_model_t alignment_model,
    swg_penalties_t* swg_penalties,
    sequence_t* const sequence_end1,
    sequence_t* const sequence_end2,
    paired_matches_t* const paired_matches,
    const archive_check_type check_type,
    mm_stack_t* const mm_stack) {
  // Check individually end1
  archive_check_se_matches(archive,alignment_model,swg_penalties,sequence_end1,
      paired_matches->matches_end1,check_type,mm_stack);
  // Check individually end2
  archive_check_se_matches(archive,alignment_model,swg_penalties,sequence_end2,
      paired_matches->matches_end2,check_type,mm_stack);
  // TODO More checks related with template-size orientation, etc (But this might be stats better)
}


