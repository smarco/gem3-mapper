/*
 * PROJECT: GEMMapper
 * FILE: match_align.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "match_align.h"
#include "match_align_dto.h"

/*
 * Curate Alignment
 */
GEM_INLINE bool match_align_curate_cigar_trim(
    match_align_parameters_t* const align_parameters,const cigar_element_t* const cigar_element,
    uint64_t* const trim_length,uint64_t* const match_position) {
  switch (cigar_element->type) {
    case cigar_match:
      // Small Match
      if (cigar_element->match_length <= align_parameters->min_context_length) {
        *trim_length += cigar_element->match_length;
        if (match_position!=NULL) *match_position += cigar_element->match_length;
        return true;
      }
    break;
    case cigar_mismatch:
      ++(*trim_length);
      if (match_position!=NULL) ++(*match_position);
      return true;
      break;
    case cigar_soft_trim:
    case cigar_del:
      *trim_length += cigar_element->indel.indel_length;
      return true;
      break;
    case cigar_ins:
      if (match_position!=NULL) *match_position += cigar_element->indel.indel_length;
      return true;
      break;
    default:
      return false;
      break;
  }
  return false;
}
GEM_INLINE void match_align_curate_cigar(
    match_trace_t* const match_trace,vector_t* const cigar_vector,
    match_align_parameters_t* const align_parameters) {
  // Parameters
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  // Prepare CIGAR-buffer to curate
  cigar_element_t* const cigar_buffer = vector_get_elm(cigar_vector,match_alignment->cigar_offset,cigar_element_t);
  cigar_element_t* last_curated_element = NULL;
  const uint64_t cigar_length = match_alignment->cigar_length;
  uint64_t curated_cigar_length = 0, i = 0;
  // Handle first CIGAR element(s)
  while (i<cigar_length) {
    cigar_element_t* const cigar_element = cigar_buffer + i;
    uint64_t indel_length = 0;
    if (match_align_curate_cigar_trim(align_parameters,cigar_element,&indel_length,&match_alignment->match_position)) {
      if (curated_cigar_length > 0) {
        cigar_buffer->indel.indel_length += indel_length;
      } else {
        cigar_buffer->type = cigar_soft_trim;
        cigar_buffer->indel.indel_length = indel_length;
        ++curated_cigar_length;
      }
      last_curated_element = cigar_buffer;
      ++i;
    } else {
      break;
    }
  }
  // Traverse all CIGAR elements
  for (;i<cigar_length;) {
    cigar_element_t* const cigar_element = cigar_buffer + i;
    cigar_element_t* const curated_element = cigar_buffer + curated_cigar_length;
    switch (cigar_element->type) {
      case cigar_match:
        last_curated_element = NULL;
        break;
      case cigar_mismatch:
        if (last_curated_element!=NULL && (last_curated_element->type==cigar_del || last_curated_element->type==cigar_soft_trim)) {
          if (i+1<cigar_length && (cigar_buffer[i+1].type==cigar_del || cigar_buffer[i+1].type==cigar_soft_trim)) {
            ++last_curated_element->indel.indel_length;
            ++i;
            continue;
          }
        }
        last_curated_element = NULL;
        break;
      case cigar_del:
      case cigar_soft_trim:
        // Check previous element
        if (last_curated_element!=NULL && (last_curated_element->type==cigar_del || last_curated_element->type==cigar_soft_trim)) {
          last_curated_element->indel.indel_length += cigar_element->indel.indel_length;
          ++i;
          continue;
        }
        last_curated_element = curated_element;
        break;
      case cigar_ins:
        // Check Insertion/Deletion proper order
        if (i+1<cigar_length && (cigar_buffer[i+1].type==cigar_del || cigar_buffer[i+1].type==cigar_soft_trim)) {
          SWAP(cigar_buffer[i],cigar_buffer[i+1]); continue;
        }
        // Check previous element
        if (last_curated_element!=NULL && last_curated_element->type==cigar_ins) {
          last_curated_element->indel.indel_length += cigar_element->indel.indel_length; // FIXME Combine indel text !!
          ++i;
          continue;
        }
        last_curated_element = curated_element;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
    if (i != curated_cigar_length) *curated_element = *cigar_element;
    ++curated_cigar_length;
    ++i;
  }
  // Handle last CIGAR element(s)
  if (curated_cigar_length > 0) {
    // Check match or deletion at the end
    last_curated_element = cigar_buffer + (curated_cigar_length-1);
    uint64_t indel_length = 0;
    if (match_align_curate_cigar_trim(align_parameters,last_curated_element,&indel_length,NULL)) {
      // Chain all the mismatches & deletions at the end
      while (last_curated_element > cigar_buffer) {
        --last_curated_element;
        if (!match_align_curate_cigar_trim(align_parameters,last_curated_element,&indel_length,NULL)) break;
        --curated_cigar_length;
      }
      // Merge all of them
      last_curated_element = cigar_buffer + (curated_cigar_length-1);
      last_curated_element->type = cigar_soft_trim;
      last_curated_element->indel.indel_length = indel_length;
    }
  }
  // Set curated-CIGAR length
  match_alignment->cigar_length = curated_cigar_length;
}
/*
 * Exact Alignment
 *   @align_input->key_length
 *   @align_input->text_position
 *   @align_input->text_trace_offset
 *   @align_input->text_offset_end
 *   @align_parameters->emulated_rc_search
 *   @align_parameters->swg_penalties
 *   @match_trace->match_alignment.score
 */
GEM_INLINE void match_align_exact(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters) {
  PROF_START(GP_MATCHES_ALIGN_EXACT);
  // Parameters
  const uint64_t key_length = align_input->key_length;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  // Configure match-trace
  match_trace->trace_offset = align_input->text_trace_offset;
  match_trace->sequence_name = NULL;
  match_trace->text_position = UINT64_MAX;
  match_trace->emulated_rc_search = align_parameters->emulated_rc_search;
  match_trace->distance = match_alignment->score;
  match_trace->swg_score = swg_score_match(align_parameters->swg_penalties,(int32_t)key_length);
  // Insert exact-match CIGAR
  match_alignment->match_position = (match_alignment->score==0) ? // Adjust position (if needed)
      align_input->text_position + (align_input->text_offset_end - key_length) : align_input->text_position;
  match_alignment->cigar_offset = vector_get_used(matches->cigar_vector);
  match_alignment->cigar_length = 1;
  match_alignment->effective_length = key_length;
  cigar_element_t cigar_element;
  cigar_element.type = cigar_match;
  cigar_element.match_length = key_length;
  vector_insert(matches->cigar_vector,cigar_element,cigar_element_t);
  PROF_STOP(GP_MATCHES_ALIGN_EXACT);
}
/*
 * Hamming Alignment (Only mismatches)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text_trace_offset
 *   @align_input->text_position
 *   @align_input->text
 *   @align_input->text_offset_begin
 *   @align_parameters->emulated_rc_search
 *   @align_parameters->allowed_enc
 *   @align_parameters->swg_penalties
 */
GEM_INLINE void match_align_hamming(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters) {
  PROF_START(GP_MATCHES_ALIGN_HAMMING);
  // Parameters
  const uint8_t* const key = align_input->key;
  const uint64_t key_length = align_input->key_length;
  uint8_t* const text = align_input->text;
  const uint64_t text_offset_begin = align_input->text_offset_begin;
  const uint64_t text_end_offset = text_offset_begin + key_length;
  const bool* const allowed_enc = align_parameters->allowed_enc;
  // Configure match-trace
  match_trace->trace_offset = align_input->text_trace_offset;
  match_trace->sequence_name = NULL;
  match_trace->text_position = UINT64_MAX;
  match_trace->emulated_rc_search = align_parameters->emulated_rc_search;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_alignment->match_position = align_input->text_position;
  match_alignment->cigar_offset = vector_get_used(matches->cigar_vector);
  // Hamming Check
  vector_reserve_additional(matches->cigar_vector,key_length);
  cigar_element_t* cigar_element = vector_get_free_elm(matches->cigar_vector,cigar_element_t);
  uint64_t i, mismatches;
  for (i=text_offset_begin,mismatches=0;i<text_end_offset;++i) {
    const uint8_t candidate_enc = text[i];
    // Check Mismatch
    if (!allowed_enc[candidate_enc] || candidate_enc != key[i]) {
      ++mismatches;
      cigar_element->type = cigar_mismatch;
      cigar_element->mismatch = candidate_enc;
    }
  }
  vector_update_used(matches->cigar_vector,cigar_element);
  match_alignment->cigar_length = mismatches;
  match_alignment->effective_length = key_length;
  match_trace->distance = mismatches;
  match_trace->swg_score = swg_score_cigar(align_parameters->swg_penalties,
      matches->cigar_vector,match_alignment->cigar_offset,match_alignment->cigar_length);
  PROF_STOP(GP_MATCHES_ALIGN_HAMMING);
}
/*
 * Levenshtein
 *   @align_input->key
 *   @align_input->bpm_pattern
 *   @align_input->text_trace_offset
 *   @align_input->text_position
 *   @align_input->text
 *   @align_input->text_offset_begin
 *   @align_input->text_length
 *   @align_parameters->emulated_rc_search
 *   @align_parameters->max_error
 *   @align_parameters->left_gap_alignment
 */
GEM_INLINE void match_align_levenshtein(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    mm_stack_t* const mm_stack) {
  PROF_START(GP_MATCHES_ALIGN_LEVENSHTEIN);
  // Configure match-trace
  match_trace->trace_offset = align_input->text_trace_offset;
  match_trace->sequence_name = NULL;
  match_trace->text_position = UINT64_MAX;
  match_trace->emulated_rc_search = align_parameters->emulated_rc_search;
  // Levenshtein Align
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  match_alignment->match_position = align_input->text_position + align_input->text_offset_begin;
  bpm_align_match(align_input,align_parameters->max_error,
      align_parameters->left_gap_alignment,match_alignment,matches->cigar_vector,mm_stack);
  match_trace->distance = match_alignment->score;
  match_trace->swg_score = swg_score_cigar(align_parameters->swg_penalties,
      matches->cigar_vector,match_alignment->cigar_offset,match_alignment->cigar_length);
  PROF_STOP(GP_MATCHES_ALIGN_LEVENSHTEIN);
}
/*
 * Smith-Waterman-Gotoh
 */
/*
 * Align Matching Region (@region_matching)
 *   @align_input->key
 *   @align_input->text
 *   @align_parameters->swg_penalties
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->allowed_enc
 *   @align_parameters->max_bandwidth
 *   @match_alignment->cigar_length (Cumulative)
 */
GEM_INLINE void match_align_matching_region(
    const region_matching_t* const region_matching,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    match_alignment_t* const match_alignment,vector_t* const cigar_vector,mm_stack_t* const mm_stack) {
  // Parameters
  uint8_t* const key = align_input->key;
  const uint64_t key_matching_length = region_matching->key_end-region_matching->key_begin;
  uint8_t* const text = align_input->text;
  const uint64_t text_matching_length = region_matching->text_end-region_matching->text_begin;
  // Select matching-region type
  switch (region_matching->matching_type) {
    case region_matching_exact:
      matches_cigar_vector_append_match(cigar_vector,&match_alignment->cigar_length,key_matching_length);
      match_alignment->score = swg_score_match(align_parameters->swg_penalties,key_matching_length);
      break;
    case region_matching_approximate:
      if (region_matching->cigar_length > 0) {
        // Copy the CIGAR from the matching region
        const uint64_t cigar_buffer_offset = region_matching->cigar_buffer_offset;
        uint64_t i;
        match_alignment->score = 0;
        for (i=0;i<region_matching->cigar_length;++i) {
          cigar_element_t* const scaffolding_elm = vector_get_elm(cigar_vector,cigar_buffer_offset+i,cigar_element_t);
          match_alignment->score += swg_score_cigar_element(align_parameters->swg_penalties,scaffolding_elm);
          matches_cigar_vector_append_cigar_element(cigar_vector,&match_alignment->cigar_length,scaffolding_elm);
        }
      } else {
        // Compute the CIGAR from the matching region
        match_align_input_t align_chunk_input = {
            .key = key+region_matching->key_begin,
            .key_length = key_matching_length,
            .text = text+region_matching->text_begin,
            .text_length = text_matching_length,
        };
        match_align_parameters_t align_chunk_parameters = {
            .max_bandwidth = region_matching->error+1, // FIXME Include CIGAR when generating approximate regions
            .left_gap_alignment = align_parameters->left_gap_alignment,
            .allowed_enc = align_parameters->allowed_enc,
            .swg_penalties = align_parameters->swg_penalties,
        };
        swg_align_match(&align_chunk_input,&align_chunk_parameters,false,false,match_alignment,cigar_vector,mm_stack);
      }
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
/*
 * Chained SWG Alignment
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text
 *   @align_input->text_length
 *   @align_parameters->swg_penalties
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->allowed_enc
 *   @align_parameters->max_bandwidth
 *   @match_scaffold->scaffold_regions
 *   @match_scaffold->num_scaffold_regions
 *   @match_alignment->cigar_length (Cumulative)
 *   @match_alignment->effective_length (Cumulative)
 *   @match_alignment->score (Cumulative)
 */
GEM_INLINE void match_align_smith_waterman_gotoh_chained(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,vector_t* const cigar_vector,mm_stack_t* const mm_stack) {
  // Parameters
  uint8_t* const key = align_input->key;
  const uint64_t key_length = align_input->key_length;
  uint8_t* const text = align_input->text;
  const uint64_t text_length = align_input->text_length;
  const uint64_t max_bandwidth = align_parameters->max_bandwidth;
  match_alignment_t* const match_alignment = &match_trace->match_alignment;
  region_matching_t* const scaffold_regions = match_scaffold->scaffold_regions;
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  // Auxiliary Variables
  match_align_input_t align_chunk_input;
  uint64_t key_chunk_begin_offset, key_chunk_length, i, cigar_length, cigar_vector_used;
  uint64_t text_chunk_end_offset, text_chunk_begin_offset, text_chunk_length, match_begin_position;
  int32_t total_score = 0;
  // Chain to FIRST matching region
  const region_matching_t* const first_region_matching = scaffold_regions;
  key_chunk_length = first_region_matching->key_begin;
  text_chunk_length = BOUNDED_ADDITION(key_chunk_length,max_bandwidth,first_region_matching->text_begin);
  text_chunk_begin_offset = first_region_matching->text_begin - text_chunk_length;
  align_chunk_input.key = key;
  align_chunk_input.key_length = key_chunk_length;
  align_chunk_input.text = text + text_chunk_begin_offset;
  align_chunk_input.text_length = text_chunk_length;
  match_alignment->match_position += text_chunk_begin_offset; // Offset match position
  // Keep CIGAR state + Position & align
  cigar_length = match_alignment->cigar_length;
  cigar_vector_used = vector_get_used(matches->cigar_vector);
  match_begin_position = match_alignment->match_position;
  swg_align_match(&align_chunk_input,align_parameters,true,false,match_alignment,matches->cigar_vector,mm_stack);
  // Check alignment (Rollback & trim if needed)
  if (match_alignment->score == ALIGN_DISTANCE_INF || match_alignment->score < 0) {
    // Restore the CIGAR
    match_alignment->cigar_length = cigar_length;
    vector_set_used(matches->cigar_vector,cigar_vector_used);
    // Trim the alignment
    if (key_chunk_length > 0) {
      matches_cigar_vector_append_deletion(cigar_vector,&match_alignment->cigar_length,key_chunk_length);
      total_score += swg_score_deletion(align_parameters->swg_penalties,key_chunk_length); // Update total score
    }
    // Save match position
    match_begin_position = match_begin_position + text_chunk_length;
  } else {
    total_score += match_alignment->score; // Update total score
    match_begin_position = match_alignment->match_position; // Save match position
  }
  // Add matching region to CIGAR
  match_align_matching_region(first_region_matching,align_input,align_parameters,match_alignment,cigar_vector,mm_stack);
  total_score += match_alignment->score;
  // Chain matching regions
  for (i=1;i<num_scaffold_regions;++i) {
    const region_matching_t* const prev_region_matching = scaffold_regions + (i-1);
    const region_matching_t* const current_region_matching = scaffold_regions + i;
    // Align gap between regions matching
    const uint64_t key_chunk_begin_offset = prev_region_matching->key_end;
    const uint64_t key_chunk_length = current_region_matching->key_begin - key_chunk_begin_offset;
    const uint64_t text_chunk_begin_offset = prev_region_matching->text_end;
    const uint64_t text_chunk_length = current_region_matching->text_begin - text_chunk_begin_offset;
    align_chunk_input.key = key+key_chunk_begin_offset;
    align_chunk_input.key_length = key_chunk_length;
    align_chunk_input.text = text+text_chunk_begin_offset;
    align_chunk_input.text_length = text_chunk_length;
    // Keep CIGAR state & align
    cigar_length = match_alignment->cigar_length;
    cigar_vector_used = vector_get_used(matches->cigar_vector);
    swg_align_match(&align_chunk_input,align_parameters,false,false,match_alignment,matches->cigar_vector,mm_stack);
    // Check alignment result
    if (match_alignment->score == ALIGN_DISTANCE_INF /* || TODO */ ) {
      // Restore the CIGAR
      match_alignment->cigar_length = cigar_length;
      vector_set_used(matches->cigar_vector,cigar_vector_used);
      // Delete the read chunk
      if (key_chunk_length > 0) {
        matches_cigar_vector_append_deletion(cigar_vector,&match_alignment->cigar_length,key_chunk_length);
        total_score += swg_score_deletion(align_parameters->swg_penalties,key_chunk_length);
      }
      // Insert the text chunk
      if (text_chunk_length > 0) {
        matches_cigar_vector_append_insertion(cigar_vector,&match_alignment->cigar_length,text_chunk_length,align_chunk_input.text);
        total_score += swg_score_insertion(align_parameters->swg_penalties,text_chunk_length);
      }
    } else {
      total_score += match_alignment->score;
    }
    // Add matching region to CIGAR
    match_align_matching_region(current_region_matching,align_input,align_parameters,match_alignment,cigar_vector,mm_stack);
    total_score += match_alignment->score;
  }
  // Chain from LAST matching region
  const region_matching_t* const last_region_matching = scaffold_regions + (num_scaffold_regions-1);
  key_chunk_begin_offset = last_region_matching->key_end;
  key_chunk_length = key_length - last_region_matching->key_end;
  text_chunk_end_offset = BOUNDED_ADDITION(last_region_matching->text_end,key_chunk_length+max_bandwidth,text_length);
  text_chunk_begin_offset = last_region_matching->text_end;
  text_chunk_length = text_chunk_end_offset-last_region_matching->text_end;
  align_chunk_input.key = key+key_chunk_begin_offset;
  align_chunk_input.key_length = key_chunk_length;
  align_chunk_input.text = text+text_chunk_begin_offset;
  align_chunk_input.text_length = text_chunk_length;
  // Keep CIGAR state & align
  cigar_length = match_alignment->cigar_length;
  cigar_vector_used = vector_get_used(matches->cigar_vector);
  swg_align_match(&align_chunk_input,align_parameters,false,true,match_alignment,matches->cigar_vector,mm_stack);
  // Check alignment (Rollback & trim if needed)
  if (match_alignment->score == ALIGN_DISTANCE_INF || match_alignment->score < 0) {
    // Restore the CIGAR
    match_alignment->cigar_length = cigar_length;
    vector_set_used(matches->cigar_vector,cigar_vector_used);
    // Trim the alignment
    if (key_chunk_length > 0) {
      matches_cigar_vector_append_deletion(cigar_vector,&match_alignment->cigar_length,key_chunk_length);
      total_score += swg_score_deletion(align_parameters->swg_penalties,key_chunk_length); // Update total score
    }
  }
  // Post-processing
  match_alignment->match_position = match_begin_position; // Restore match position
  match_alignment->score = total_score;
}
/*
 * Smith-Waterman-Gotoh Alignment (Gap-affine)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->text_trace_offset
 *   @align_input->text_position
 *   @align_input->text
 *   @align_input->text_length
 *   @align_input->text_offset_begin
 *   @align_input->text_offset_end
 *   @align_parameters->emulated_rc_search
 *   @align_parameters->swg_penalties
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->allowed_enc
 *   @align_parameters->max_bandwidth
 *   @align_parameters->cigar_curation
 *   @match_scaffold->scaffold_regions
 *   @match_scaffold->num_scaffold_regions
 */
GEM_INLINE void match_align_smith_waterman_gotoh(
    matches_t* const matches,match_trace_t* const match_trace,
    match_align_input_t* const align_input,match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
  PROF_START(GP_MATCHES_ALIGN_SWG);
  // Configure match-trace
  match_trace->trace_offset = align_input->text_trace_offset;
  match_trace->sequence_name = NULL;
  match_trace->text_position = UINT64_MAX;
  match_trace->emulated_rc_search = align_parameters->emulated_rc_search;
  match_trace->match_alignment.match_position = align_input->text_position;
  match_trace->match_alignment.cigar_offset = vector_get_used(matches->cigar_vector);
  match_trace->match_alignment.cigar_length = 0;
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
#ifdef GEM_DEBUG
  match_trace->match_scaffold = (match_scaffold!=NULL && num_scaffold_regions > 0) ? match_scaffold : NULL;
#endif
  // Check the number of matching regions
  if (match_scaffold!=NULL && num_scaffold_regions > 0) {
    // Chain matching regions and align gaps (SWG)
    match_align_smith_waterman_gotoh_chained(matches,match_trace,
        align_input,align_parameters,match_scaffold,matches->cigar_vector,mm_stack);
  } else {
    // Full SWG
    align_input->text = align_input->text + align_input->text_offset_begin;
    align_input->text_length = align_input->text_offset_end - align_input->text_offset_begin;
    match_trace->match_alignment.match_position += align_input->text_offset_begin;
    swg_align_match(align_input,align_parameters,true,true,
        &match_trace->match_alignment,matches->cigar_vector,mm_stack);
  }
  // Check for Bad alignment (discarded)
  if (match_trace->match_alignment.score == ALIGN_DISTANCE_INF) {
    match_trace->swg_score = ALIGN_DISTANCE_INF;
    match_trace->distance = ALIGN_DISTANCE_INF;
  } else {
    // Curate alignment
    if (align_parameters->cigar_curation) {
      match_align_curate_cigar(match_trace,matches->cigar_vector,align_parameters);
    }
    // Compute Score + Effective-Length
    match_alignment_t* const match_alignment = &match_trace->match_alignment;
    match_alignment->effective_length = matches_cigar_effective_length(
        matches->cigar_vector,match_alignment->cigar_offset,match_alignment->cigar_length);
    match_trace->swg_score = swg_score_cigar(align_parameters->swg_penalties,
        matches->cigar_vector,match_alignment->cigar_offset,match_alignment->cigar_length);
    match_trace->distance = matches_cigar_compute_edit_distance(matches,
        match_trace->match_alignment.cigar_offset,match_trace->match_alignment.cigar_length);
  }
  PROF_STOP(GP_MATCHES_ALIGN_SWG);
}
