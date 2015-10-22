/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "match_scaffold.h"
#include "match_align_dto.h"
#include "matches_cigar.h"
#include "align.h"
#include "align_bpm.h"
#include "align_ond.h"

/*
 * Debug
 */
#define DEBUG_MATCH_SCAFFOLD       GEM_DEEP_DEBUG
#define DEBUG_MATCH_SCAFFOLD_EDIT  GEM_DEEP_DEBUG
#define DEBUG_MATCH_SCAFFOLD_SWG   GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

/*
 * Setup
 */
GEM_INLINE void match_scaffold_init(match_scaffold_t* const match_scaffold) {
  // Scaffold matching region
  match_scaffold->scaffold_regions = NULL;
  match_scaffold->num_scaffold_regions = 0;
  match_scaffold->scaffolding_coverage = 0;
  // Underlying Alignment
  match_scaffold->match_alignment.cigar_length = 0;
}
/*
 * Scaffold Match (based on levenshtein alignment)
 */
GEM_INLINE bool match_scaffold_match_is_homopolymer(
    const uint8_t* const text,const uint64_t text_pos,
    const uint64_t homopolymer_min_context,const uint8_t* indel_text,
    const uint64_t indel_length) {
  if (text_pos < homopolymer_min_context) return false;
  // Search for minimum context support for an homopolymer event
  const uint64_t context_init_position = text_pos - homopolymer_min_context;
  const uint64_t context_end_position = text_pos;
  const uint8_t monomer = text[context_init_position];
  // Check context
  uint64_t j;
  for (j=context_init_position+1;j<=context_end_position;j++) {
    if (monomer!=text[j]) return false;
  }
  // Check indel-text
  for (j=0;j<indel_length;++j) {
    if (monomer!=indel_text[j]) return false;
  }
  return true;
}
GEM_INLINE void match_scaffold_match_add_match_matching_approximate(
    const uint64_t key_begin,const uint64_t key_end,
    const uint64_t text_begin,const uint64_t text_end,
    const uint64_t cigar_offset,const uint64_t cigar_length,
    const uint64_t score,match_scaffold_t* const match_scaffold) {
  // Add matching region
  region_matching_t* const region_matching = match_scaffold->scaffold_regions + match_scaffold->num_scaffold_regions;
  ++match_scaffold->num_scaffold_regions;
  // Setup matching-region
  region_matching->matching_type = region_matching_approximate;
  region_matching->error = score;
  region_matching->key_begin = key_begin;
  region_matching->key_end = key_end;
  region_matching->text_begin = text_begin;
  region_matching->text_end = text_end;
  region_matching->cigar_buffer_offset = cigar_offset;
  region_matching->cigar_length = cigar_length;
}
GEM_INLINE region_matching_t* match_scaffold_match_add_match(
    uint8_t* const text,uint64_t* const key_pos,uint64_t* const text_pos,
    const uint64_t cigar_offset,const uint64_t match_length,
    match_scaffold_t* const match_scaffold) {
  // Add matching region
  region_matching_t* const region_matching = match_scaffold->scaffold_regions + match_scaffold->num_scaffold_regions;
  ++match_scaffold->num_scaffold_regions;
  // Set-up matching region
  region_matching->matching_type = region_matching_exact;
  region_matching->key_begin = *key_pos;
  region_matching->text_begin = *text_pos;
  region_matching->cigar_buffer_offset = cigar_offset;
  region_matching->cigar_length = 1;
  *key_pos += match_length;
  *text_pos += match_length;
  match_scaffold->scaffolding_coverage += match_length;
  region_matching->key_end = *key_pos;
  region_matching->text_end = *text_pos;
  // Return matching region
  return region_matching;
}
GEM_INLINE void match_scaffold_match_add_mismatch(
    uint64_t* const key_pos,uint64_t* const text_pos,const uint64_t match_length,
    match_scaffold_t* const match_scaffold,region_matching_t* const last_scaffold_region) {
  // Extend matching-region
  last_scaffold_region->matching_type = region_matching_approximate;
  last_scaffold_region->cigar_length += 2; // Add the mismatch + matching
  *key_pos += match_length;
  *text_pos += match_length;
  match_scaffold->scaffolding_coverage += match_length;
  last_scaffold_region->key_end = *key_pos;
  last_scaffold_region->text_end = *text_pos;
}
GEM_INLINE void match_scaffold_match_add_homopolymer(
    uint8_t* const text,uint64_t* const key_pos,uint64_t* const text_pos,
    const uint64_t context_cigar_offset,const uint64_t context_match_length,
    const cigar_t indel_type,const uint64_t indel_length,
    match_scaffold_t* const match_scaffold,region_matching_t* const last_scaffold_region) {
  // Check previous matching-region
  if (last_scaffold_region!=NULL) {
    region_matching_t* const region_matching = last_scaffold_region;
    region_matching->matching_type = region_matching_approximate;
    ++region_matching->cigar_length;
    if (indel_type==cigar_ins) {
      *text_pos += indel_length;
      region_matching->text_end = *text_pos;
    } else { // cigar_del
      *key_pos += indel_length;
      region_matching->key_end = *key_pos;
      match_scaffold->scaffolding_coverage += indel_length;
    }
  } else {
    // Add matching region
    region_matching_t* const region_matching = match_scaffold->scaffold_regions + match_scaffold->num_scaffold_regions;
    ++match_scaffold->num_scaffold_regions;
    // Set-up homopolymer matching region
    region_matching->matching_type = region_matching_approximate;
    region_matching->cigar_buffer_offset = context_cigar_offset;
    region_matching->cigar_length = 2;
    region_matching->key_begin = *key_pos - context_match_length;
    region_matching->text_begin = *text_pos - context_match_length;
    if (indel_type==cigar_ins) {
      *text_pos += indel_length;
      match_scaffold->scaffolding_coverage += context_match_length;
    } else { // cigar_del
      *key_pos += indel_length;
      match_scaffold->scaffolding_coverage += context_match_length + indel_length;
    }
    region_matching->key_end = *key_pos;
    region_matching->text_end = *text_pos;
  }
}
/*
 * Compose the scaffolding (determine the matching regions based on the levenshtein alignment)
 *   @align_input->key
 *   @align_input->text
 *   @align_parameters->min_matching_length
 *   @align_parameters->min_context_length
 */
GEM_INLINE void match_scaffold_levenshtein_compose_matching_regions(
    matches_t* const matches,match_align_input_t* const align_input,
    uint64_t key_pos,uint64_t text_pos,const uint64_t cigar_offset,
    const uint64_t cigar_length,match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
  // Parameters
  uint8_t* const text = align_input->text;
  const uint64_t matching_min_length = align_parameters->scaffolding_matching_min_length;
  cigar_element_t* const cigar_buffer = vector_get_elm(matches->cigar_vector,cigar_offset,cigar_element_t);
  // Init scaffolding
  match_scaffold->scaffold_regions = mm_stack_calloc(mm_stack,cigar_length,region_matching_t,false);
  match_scaffold->num_scaffold_regions = 0;
  match_scaffold->scaffolding_coverage = 0;
  // Traverse CIGAR elements and compose scaffolding
  region_matching_t* last_scaffold_region = NULL;
  uint64_t i = 0;
  for (i=0;i<cigar_length;) {
    switch (cigar_buffer[i].type) {
      case cigar_match: {
        // Check matching length
        const uint64_t match_length = cigar_buffer[i].length;
        if (match_length < matching_min_length) {
          text_pos += match_length;
          key_pos += match_length;
          last_scaffold_region = NULL; // Nullify last region-matching
        } else {
          last_scaffold_region = match_scaffold_match_add_match(text,&key_pos,&text_pos,
              cigar_offset + i,match_length,match_scaffold); // Add Match
        }
        ++i;
        break;
      }
      case cigar_mismatch:
        // Tolerate mismatches if well delimited between exact matching regions
        if (last_scaffold_region!=NULL && i+1<cigar_length &&
            cigar_buffer[i+1].type==cigar_match && cigar_buffer[i+1].length >= matching_min_length) {
          match_scaffold_match_add_mismatch(&key_pos,&text_pos,
              cigar_buffer[i+1].length+1,match_scaffold,last_scaffold_region);
          i += 2;
        } else {
          last_scaffold_region = NULL; // Nullify last region-matching
          ++text_pos; ++key_pos;
          ++i;
        }
        break;
      case cigar_ins:
        text_pos += cigar_buffer[i].length;
        last_scaffold_region = NULL; // Nullify last region-matching
        ++i;
        break;
      case cigar_del:
        key_pos += cigar_buffer[i].length;
        last_scaffold_region = NULL; // Nullify last region-matching
        ++i;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
/*
 * Compose the scaffolding (determine the matching regions based on the SWG score)
 */
GEM_INLINE void match_scaffold_smith_waterman_gotoh_compose_matching_regions(
    matches_t* const matches,match_align_input_t* const align_input,
    uint64_t key_pos,uint64_t text_pos,const uint64_t cigar_offset,
    const uint64_t cigar_length,match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
  // Parameters
  const swg_penalties_t* const swg_penalties = align_parameters->swg_penalties;
  cigar_element_t* const cigar_buffer = vector_get_elm(matches->cigar_vector,cigar_offset,cigar_element_t);
  // Init scaffolding
  match_scaffold->scaffold_regions = mm_stack_calloc(mm_stack,cigar_length,region_matching_t,false);
  match_scaffold->num_scaffold_regions = 0;
  match_scaffold->scaffolding_coverage = 0;
  // Traverse CIGAR elements and compose scaffolding
  int64_t i, local_score = 0, max_local_score = 0;
  uint64_t max_local_begin_key_pos = key_pos, max_local_end_key_pos = key_pos;
  uint64_t max_local_begin_text_pos = text_pos, max_local_end_text_pos = text_pos;
  uint64_t max_local_begin_cigar = 0, max_local_end_cigar = 0;
  for (i=0;i<cigar_length;) {
    // Add CIGAR operation
    switch (cigar_buffer[i].type) {
      case cigar_match: {
        const int32_t match_length = cigar_buffer[i].length;
        local_score += align_swg_score_match(swg_penalties,match_length);
        text_pos += match_length;
        key_pos += match_length;
        break;
      }
      case cigar_mismatch:
        local_score += align_swg_score_mismatch(swg_penalties);
        ++text_pos;
        ++key_pos;
        break;
      case cigar_ins: {
        const int32_t indel_length = cigar_buffer[i].length;
        local_score += align_swg_score_insertion(swg_penalties,indel_length);
        text_pos += indel_length;
        break;
      }
      case cigar_del: {
        const int32_t indel_length = cigar_buffer[i].length;
        local_score += align_swg_score_deletion(swg_penalties,indel_length);
        key_pos += indel_length;
        break;
      }
      default:
        GEM_INVALID_CASE();
        break;
    }
    ++i;
    // Check local score
    if (local_score > max_local_score) {
      // Max local score
      max_local_score = local_score;
      max_local_end_key_pos = key_pos;
      max_local_end_text_pos = text_pos;
      max_local_end_cigar = i;
    } else if (local_score < 0) {
      // Store local-max alignment chunk
      if (max_local_score >= align_parameters->swg_threshold) {
        match_scaffold_match_add_match_matching_approximate(
            max_local_begin_key_pos,max_local_end_key_pos,
            max_local_begin_text_pos,max_local_end_text_pos,
            cigar_offset+max_local_begin_cigar,max_local_end_cigar-max_local_begin_cigar,
            max_local_score,match_scaffold); // Add approximate matching region
      }
      // Reset
      local_score = 0;
      max_local_score = 0;
      max_local_begin_key_pos = key_pos;
      max_local_begin_text_pos = text_pos;
      max_local_begin_cigar = i;
    }
  }
  // Check local score
  if (max_local_score >= align_parameters->swg_threshold) {
    match_scaffold_match_add_match_matching_approximate(
        max_local_begin_key_pos,max_local_end_key_pos,
        max_local_begin_text_pos,max_local_end_text_pos,
        cigar_offset+max_local_begin_cigar,max_local_end_cigar-max_local_begin_cigar,
        max_local_score,match_scaffold); // Add approximate matching region
  }
}
/*
 * Scaffold the alignment (based on levenshtein alignment)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->bpm_pattern
 *   @align_input->text_position
 *   @align_input->text
 *   @align_input->text_length
 *   @align_input->text_offset_begin
 *   @align_input->text_offset_end
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->max_error
 */
GEM_INLINE void match_scaffold_levenshtein_alignment(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
  // Fill Matrix (Pv,Mv)
  mm_stack_push_state(mm_stack); // Save stack state
  match_align_input_t bpm_align_input = {
      .key = align_input->key,
      .key_length = align_input->key_length,
      .bpm_pattern = align_input->bpm_pattern,
      .text = align_input->text + align_input->text_offset_begin,
      .text_length = align_input->text_offset_end - align_input->text_offset_begin,
  };
  bpm_align_matrix_t bpm_align_matrix;
  align_bpm_compute_matrix(&bpm_align_input,align_parameters->max_error,&bpm_align_matrix,mm_stack);
  // Set distance
  match_scaffold->match_alignment.score = bpm_align_matrix.min_score;
  if (bpm_align_matrix.min_score == ALIGN_DISTANCE_INF) {
    match_scaffold->match_alignment.cigar_length = 0;
    mm_stack_pop_state(mm_stack,false); // Free
    return;
  }
  // Backtrace and generate CIGAR
  match_alignment_t* const match_alignment = &match_scaffold->match_alignment;
  const uint64_t match_position = align_input->text_position + align_input->text_offset_begin;
  match_alignment->match_position = match_position; // Sentinel
  align_bpm_backtrace_matrix(&bpm_align_input,align_parameters->left_gap_alignment,
      &bpm_align_matrix,match_alignment,matches->cigar_vector);
  mm_stack_pop_state(mm_stack,false); // Free
  // Store the offset
  match_alignment->match_position = (match_alignment->match_position - match_position) + align_input->text_offset_begin;
}
/*
 * Scaffold the alignment (based on levenshtein alignment)
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->bpm_pattern
 *   @align_input->text_position
 *   @align_input->text
 *   @align_input->text_length
 *   @align_input->text_offset_begin
 *   @align_input->text_offset_end
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->min_matching_length
 *   @align_parameters->min_context_length
 */
GEM_INLINE bool match_scaffold_levenshtein(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
  PROF_INC_COUNTER(GP_MATCH_SCAFFOLD_EDIT_SCAFFOLDS);
  PROFILE_START(GP_MATCH_SCAFFOLD_EDIT,PROFILE_LEVEL);
  // Compute the levenshtein alignment (CIGAR)
  match_scaffold_levenshtein_alignment(matches,align_input,align_parameters,match_scaffold,mm_stack);
  match_alignment_t* const match_alignment = &match_scaffold->match_alignment;
  if (match_alignment->score==ALIGN_DISTANCE_INF) {
    match_scaffold->num_scaffold_regions = 0;
    match_scaffold->scaffolding_coverage = 0;
    PROFILE_STOP(GP_MATCH_SCAFFOLD_EDIT,PROFILE_LEVEL);
    return false;
  }
  // Compose matching region of the scaffolding
  match_scaffold_levenshtein_compose_matching_regions(matches,
      align_input,0,match_alignment->match_position,match_alignment->cigar_offset,
      match_alignment->cigar_length,align_parameters,match_scaffold,mm_stack);
  PROFILE_STOP(GP_MATCH_SCAFFOLD_EDIT,PROFILE_LEVEL);
  PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_EDIT_COVERAGE,
      (100*match_scaffold->scaffolding_coverage)/align_input->key_length);
  // DEBUG
  gem_cond_debug_block(DEBUG_MATCH_SCAFFOLD_EDIT) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Match.Scaffold.Edit\n");
    tab_global_inc();
    match_scaffold_print(gem_log_get_stream(),matches,match_scaffold);
    tab_global_dec();
    /*
     * Debug
     */
    /*
    match_alignment_t ond_match_alignment;
    ond_match_alignment.match_position = 0;
    match_align_input_t ond_align_input;
    ond_align_input.key = align_input->key;
    ond_align_input.key_length = align_input->key_length;
    ond_align_input.text = align_input->text; //  + align_input->align_match_begin_column;
    ond_align_input.text_length = align_input->text_length; // align_input->align_match_end_column - align_input->align_match_begin_column;
    align_ond_match(&ond_align_input,ond_align_input.text_length,
        &ond_match_alignment,matches->cigar_vector,mm_stack);
    // Display
    match_alignment_print_pretty(stderr,&ond_match_alignment,
        matches->cigar_vector,ond_align_input.key,ond_align_input.key_length,
        ond_align_input.text,ond_align_input.text_length,mm_stack);
    return match_scaffold->num_scaffold_regions > 0;
    */
    /*
     * Debug
     */
  }
  return match_scaffold->num_scaffold_regions > 0;
}
/*
 * Scaffold the alignment (based on SWG-distance)
 */
GEM_INLINE bool match_scaffold_smith_waterman_gotoh(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
  PROF_INC_COUNTER(GP_MATCH_SCAFFOLD_SWG_SCAFFOLDS);
  PROFILE_START(GP_MATCH_SCAFFOLD_SWG,PROFILE_LEVEL);
  // Compute the levenshtein alignment (CIGAR)
  match_scaffold_levenshtein_alignment(matches,align_input,align_parameters,match_scaffold,mm_stack);
  match_alignment_t* const match_alignment = &match_scaffold->match_alignment;
  if (match_alignment->score==ALIGN_DISTANCE_INF) {
    match_scaffold->num_scaffold_regions = 0;
    match_scaffold->scaffolding_coverage = 0;
    PROFILE_STOP(GP_MATCH_SCAFFOLD_SWG,PROFILE_LEVEL);
    return false;
  }
  // Compose matching region of the scaffolding
  match_scaffold_smith_waterman_gotoh_compose_matching_regions(matches,
      align_input,0,match_alignment->match_position,match_alignment->cigar_offset,
      match_alignment->cigar_length,align_parameters,match_scaffold,mm_stack);
  PROFILE_STOP(GP_MATCH_SCAFFOLD_SWG,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_MATCH_SCAFFOLD_SWG) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Match.Scaffold.SWG\n");
    match_scaffold_print(gem_log_get_stream(),matches,match_scaffold);
  }
  return match_scaffold->num_scaffold_regions > 0;
}
/*
 * Exact extend of matching regions
 */
GEM_INLINE void match_scaffold_exact_extend(
    matches_t* const matches,match_scaffold_t* const match_scaffold,
    const uint8_t* const key,const uint64_t key_length,
    const uint8_t* const text,const uint64_t text_length,
    const bool* const allowed_enc) {
  PROFILE_START(GP_MATCH_SCAFFOLD_EXTEND_REGIONS,PROFILE_LEVEL);
  // Extend all matching regions
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  const uint64_t last_region = num_scaffold_regions-1;
  uint64_t i, inc_coverage = 0;
  for (i=0;i<num_scaffold_regions;++i) {
    // Try to left extend
    region_matching_t* const region_matching = match_scaffold->scaffold_regions + i;
    const int64_t left_key_max = (i==0) ? 0 : match_scaffold->scaffold_regions[i-1].key_end;
    const int64_t left_text_max = (i==0) ? 0 : match_scaffold->scaffold_regions[i-1].text_end;
    int64_t left_key = region_matching->key_begin-1;
    int64_t left_text = region_matching->text_begin-1;
    while (left_key_max<=left_key && left_text_max<=left_text) {
      // Check match
      const uint8_t candidate_enc = text[left_text];
      if (!allowed_enc[candidate_enc] || candidate_enc != key[left_key]) break;
      --left_key;
      --left_text;
      ++inc_coverage;
    }
    region_matching->key_begin = left_key+1;
    region_matching->text_begin = left_text+1;
    // Try to right extend
    const int64_t right_key_max = (i==last_region) ? key_length-1 : match_scaffold->scaffold_regions[i+1].key_begin-1;
    const int64_t right_text_max = (i==last_region) ? text_length-1 : match_scaffold->scaffold_regions[i+1].text_begin-1;
    int64_t right_key = region_matching->key_end;
    int64_t right_text = region_matching->text_end;
    while (right_key_max>=right_key && right_text_max>=right_text) {
      // Check match
      const uint8_t candidate_enc = text[right_text];
      if (!allowed_enc[candidate_enc] || candidate_enc != key[right_key]) break;
      ++right_key;
      ++right_text;
      ++inc_coverage;
    }
    region_matching->key_end = right_key;
    region_matching->text_end = right_text;
  }
  match_scaffold->scaffolding_coverage += inc_coverage;
  PROFILE_STOP(GP_MATCH_SCAFFOLD_EXTEND_REGIONS,PROFILE_LEVEL);
}
/*
 * Chain matching regions
 *   @match_scaffold->num_scaffold_regions
 *   @match_scaffold->scaffold_regions
 */
GEM_INLINE void match_scaffold_chain_matching_regions(
    matches_t* const matches,match_scaffold_t* const match_scaffold,
    const uint8_t* const key,const uint64_t key_length,
    const uint8_t* const text,const bool* const allowed_enc,
    const uint64_t max_error,mm_stack_t* const stack) {
  PROFILE_START(GP_MATCH_SCAFFOLD_CHAIN_REGIONS,PROFILE_LEVEL);
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  // Sort matching regions
  match_scaffold_sort_regions_matching(match_scaffold);
  // Check overlapping
  bool overlapping_text = false, unsorted_read = false, max_diff = false;
  region_matching_t* last_region_matching = NULL;
  uint64_t i, coverage = 0;
  for (i=0;i<num_scaffold_regions;++i) {
    region_matching_t* const region_matching = match_scaffold->scaffold_regions + i;
    if (last_region_matching!=NULL) {
      if (last_region_matching->text_end > region_matching->text_begin) {
        overlapping_text = true; break; // Might be an deletion in the reference
      }
      if (last_region_matching->key_end > region_matching->key_begin) {
        unsorted_read = true; break;  // Might be a deletion in the read
      }
      const int64_t text_gap = region_matching->text_begin - last_region_matching->text_end;
      const int64_t key_gap = region_matching->key_begin - last_region_matching->key_end;
      const int64_t diff = text_gap-key_gap;
      if (ABS(diff) > max_error) {
        max_diff = true; break;
      }
    }
    coverage += region_matching->key_end - region_matching->key_begin;
    last_region_matching = region_matching;
  }
  if (overlapping_text || unsorted_read || max_diff) {
    match_scaffold->num_scaffold_regions = 0; // Disable region chaining
    match_scaffold->scaffolding_coverage = 0;
  } else {
    match_scaffold->scaffolding_coverage = coverage;
  }
  PROFILE_STOP(GP_MATCH_SCAFFOLD_CHAIN_REGIONS,PROFILE_LEVEL);
}
/*
 * Compute an scaffold for the alignment
 *   @align_input->key
 *   @align_input->key_length
 *   @align_input->bpm_pattern
 *   @align_input->text_position
 *   @align_input->text
 *   @align_input->text_length
 *   @align_input->text_offset_begin
 *   @align_input->text_offset_end
 *   @align_parameters->max_error
 *   @align_parameters->left_gap_alignment
 *   @align_parameters->scaffolding_min_coverage
 *   @align_parameters->scaffolding_matching_min_length
 *   @align_parameters->scaffolding_homopolymer_min_context
 *   @match_scaffold->num_scaffold_regions
 *   @match_scaffold->scaffold_regions
 */
GEM_INLINE void match_scaffold_alignment(
    matches_t* const matches,match_align_input_t* const align_input,
    match_align_parameters_t* const align_parameters,
    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
  PROFILE_START(GP_MATCH_SCAFFOLD_ALIGNMENT,PROFILE_LEVEL);
  // Parameters
  const uint8_t* const key = align_input->key;
  const uint64_t key_length = align_input->key_length;
  const uint8_t* const text = align_input->text;
  const uint64_t text_length = align_input->text_length;
  const uint64_t max_error = align_parameters->max_error;
  // Chain matching regions
  if (match_scaffold->num_scaffold_regions > 0) {
    const bool* const allowed_enc = align_parameters->allowed_enc;
    // Find a compatible chain of matching-regions
    match_scaffold_chain_matching_regions(matches,match_scaffold,key,key_length,text,allowed_enc,max_error,mm_stack);
    // Extend matching-regions as to maximize coverage
    if (match_scaffold->num_scaffold_regions > 0) {
      PROF_INC_COUNTER(GP_MATCH_SCAFFOLD_CHAIN_REGIONS_SUCCESS);
      PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_CHAIN_REGIONS_COVERAGE,(100*match_scaffold->scaffolding_coverage)/key_length);
      if (match_scaffold->scaffolding_coverage < key_length) {
        match_scaffold_exact_extend(matches,match_scaffold,key,key_length,text,text_length,allowed_enc);
      }
      PROF_ADD_COUNTER(GP_MATCH_SCAFFOLD_EXTEND_REGIONS_COVERAGE,(100*match_scaffold->scaffolding_coverage)/key_length);
    }
    // Set score as matching bases
    match_scaffold->match_alignment.score = key_length - match_scaffold->scaffolding_coverage;
  }
  // Scaffold from Levenshtein-alignment
  if (match_scaffold->scaffolding_coverage < align_parameters->scaffolding_min_coverage) {
    match_scaffold_levenshtein(matches,align_input,align_parameters,match_scaffold,mm_stack);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_MATCH_SCAFFOLD) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Match.Scaffold (scaffold alignment)\n");
    tab_global_inc();
    match_scaffold_print(stderr,matches,match_scaffold);
    tab_global_dec();
  }
  PROFILE_STOP(GP_MATCH_SCAFFOLD_ALIGNMENT,PROFILE_LEVEL);
}
/*
 * Sorting
 */
int region_matching_cmp_text_position(const region_matching_t* const a,const region_matching_t* const b) {
  return a->text_begin - b->text_begin;
}
GEM_INLINE void match_scaffold_sort_regions_matching(match_scaffold_t* const match_scaffold) {
  // Sort Scaffold regions (region_matching_t) wrt their starting position in the text
  void* array = match_scaffold->scaffold_regions;
  const size_t count = match_scaffold->num_scaffold_regions;
  qsort(array,count,sizeof(region_matching_t),(int (*)(const void *,const void *))region_matching_cmp_text_position);
}
/*
 * Display
 */
GEM_INLINE void match_scaffold_print(
    FILE* const stream,matches_t* const matches,match_scaffold_t* const match_scaffold) {
  const uint64_t num_scaffold_regions = match_scaffold->num_scaffold_regions;
  tab_fprintf(stream,"[GEM]>Matching.Scaffolded.Regions\n");
  tab_fprintf(stream,"  => Num.scaffold.regions %"PRIu64"\n",num_scaffold_regions);
  tab_fprintf(stream,"  => Scaffold.coverage %"PRIu64"\n",match_scaffold->scaffolding_coverage);
  uint64_t i;
  for (i=0;i<num_scaffold_regions;++i) {
    region_matching_t* const region_matching = match_scaffold->scaffold_regions + i;
    // Print matching region
    switch (region_matching->matching_type) {
      case region_matching_exact: tab_fprintf(stream,"    %"PRIu64"[exact]\t",i); break;
      case region_matching_approximate: tab_fprintf(stream,"    %"PRIu64"[approximate]\t",i); break;
      default: GEM_INVALID_CASE(); break;
    }
    tab_fprintf(stream,"-> [%"PRIu64",%"PRIu64") ~> [+%"PRIu64",+%"PRIu64")",
        region_matching->key_begin,region_matching->key_end,
        region_matching->text_begin,region_matching->text_end);
    // Print CIGAR
    if (matches!=NULL && region_matching->matching_type==region_matching_approximate && region_matching->cigar_length>0) {
      tab_fprintf(stream,"\tCIGAR=");
      match_cigar_print(stream,matches->cigar_vector,region_matching->cigar_buffer_offset,region_matching->cigar_length);
    }
    tab_fprintf(stream,"\n");
  }
}

