/*
 * PROJECT: GEMMapper
 * FILE: archive_search.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_search.h"

/*
 * Archive Search Setup
 */
GEM_INLINE archive_search_t* archive_search_new(
    archive_t* const archive,search_parameters_t* const search_parameters,
    select_parameters_t* const select_parameters) {
  ARCHIVE_CHECK(archive);
  // Allocate handler
  archive_search_t* const archive_search = mm_alloc(archive_search_t);
  // Archive
  archive_search->archive = archive;
  // Sequence
  archive_search->sequence = sequence_new(); // Input
  archive_search->rc_sequence = sequence_new(); // Generated
  // Approximate Search
  archive_search->search_actual_parameters.search_parameters = search_parameters;
  archive_search->select_parameters = select_parameters;
  archive_search->forward_search_state = approximate_search_new(
      archive->locator,archive->graph,archive->enc_text,archive->fm_index,
      &archive_search->search_actual_parameters);
  archive_search->reverse_search_state = approximate_search_new(
      archive->locator,archive->graph,archive->enc_text,archive->fm_index,
      &archive_search->search_actual_parameters);
  // Archive search control (Flow control) [DEFAULTS]
  archive_search->probe_strand = true;
  archive_search->search_reverse = !archive->indexed_complement;
  // Return
  return archive_search;
}
GEM_INLINE void archive_search_clear(archive_search_t* const archive_search) {
  // Clear F/R search state
  approximate_search_clear(archive_search->forward_search_state);
  approximate_search_clear(archive_search->reverse_search_state);
}
GEM_INLINE void archive_search_delete(archive_search_t* const archive_search) {
  // Delete Sequence
  sequence_delete(archive_search->sequence);
  sequence_delete(archive_search->rc_sequence);
  // Delete Approximate Search
  approximate_search_delete(archive_search->forward_search_state);
  approximate_search_delete(archive_search->reverse_search_state);
  // Free handler
  mm_free(archive_search);
}
/*
 * Archive Search [Accessors]
 */
GEM_INLINE sequence_t* archive_search_get_sequence(const archive_search_t* const archive_search) {
  return archive_search->sequence;
}
GEM_INLINE uint64_t archive_search_get_num_potential_canditates(const archive_search_t* const archive_search) {
  if (archive_search->archive->indexed_complement) {
    return archive_search->forward_search_state->num_potential_candidates;
  } else {
    return archive_search->forward_search_state->num_potential_candidates +
        archive_search->reverse_search_state->num_potential_candidates;
  }
}
/*
 * SingleEnd Indexed Search (SE Online Approximate String Search)
 */
// [Initialize]
GEM_INLINE void archive_search_prepare_sequence(archive_search_t* const archive_search,mm_stack_t* const mm_stack) {
  // Check the index characteristics & generate reverse-complement (if needed)
  if (archive_search->archive->indexed_complement) {
    archive_search->search_reverse = false;
  } else {
    if (archive_search->archive->filter_type == Iupac_colorspace_dna) {
      sequence_generate_reverse_complement(archive_search->sequence,archive_search->rc_sequence);
    } else {
      sequence_generate_reverse(archive_search->sequence,archive_search->rc_sequence);
    }
    archive_search->search_reverse =
        !(string_equals(archive_search->sequence->read,archive_search->rc_sequence->read));
  }
  // Generate the pattern(s)
  approximate_search_prepare_pattern(archive_search->forward_search_state,archive_search->sequence,mm_stack);
  if (archive_search->search_reverse) {
    approximate_search_prepare_pattern(archive_search->reverse_search_state,archive_search->rc_sequence,mm_stack);
  }
}
GEM_INLINE void archive_search_single_end(
    archive_search_t* const archive_search,matches_t* const matches,mm_search_t* const mm_search) {
  ARCHIVE_SEARCH_CHECK(archive_search);
  // Prepare pattern(s) & instantiate Search-Parameters values
  approximate_search_instantiate_values(
      &archive_search->search_actual_parameters,sequence_get_length(archive_search->sequence));
  archive_search_prepare_sequence(archive_search,mm_search->mm_stack);
  // Clean Matches
  archive_search_clear(archive_search);
  // Search the pattern(s)
  approximate_search_t* const forward_asearch = archive_search->forward_search_state;
  if (!archive_search->search_reverse) {
    // Compute the full search
    forward_asearch->stop_search_stage = asearch_end; // Don't stop until search is done
    forward_asearch->search_strand = Forward;
    approximate_search(forward_asearch,matches,mm_search);
  } else {
    // Configure search stage to stop at
    forward_asearch->stop_search_stage =
        (archive_search->search_actual_parameters.complete_strata_after_best_nominal < forward_asearch->max_differences
            && archive_search->probe_strand) ? asearch_neighborhood : asearch_end;
    // Run the search (FORWARD)
    forward_asearch->search_strand = Forward; // Configure forward search
    approximate_search(forward_asearch,matches,mm_search);
    // Check the number of matches & keep searching
    if (!forward_asearch->max_matches_reached) {
      // Keep on searching
      approximate_search_t* const reverse_asearch = archive_search->reverse_search_state;
      reverse_asearch->stop_search_stage = asearch_end; // Force a full search
      // Run the search (REVERSE)
      reverse_asearch->search_strand = Reverse; // Configure reverse search
      approximate_search(reverse_asearch,matches,mm_search);
      // Resume forward search (if not completed before)
      if (forward_asearch->current_search_stage != asearch_end && !forward_asearch->max_matches_reached) {
        forward_asearch->stop_search_stage = asearch_end;
        forward_asearch->search_strand = Forward; // Configure forward search
        approximate_search(forward_asearch,matches,mm_search);
      }
    }
  }
//  #ifdef GEM_MAPPER_DEBUG
//   const uint64_t check_level = search_params->internal_parameters.check_alignments;
//  extern bool keepon;
//  if (keepon && (check_level&CHECK_ALG_P_CORRECTNESS)) { // Check p-correctness
//     const uint64_t old_misms = search_params->max_mismatches;
//    search_params->max_mismatches = caller_max_misms; // Full check
//    fmi_check_pos_correctness(archive->index,search_params,matches,false,mpool);
//    search_params->max_mismatches = old_misms;
//  }
//  if (keepon && (check_level&CHECK_ALG_I_CORRECTNESS)) { // Check i-correctness
//     const uint64_t old_misms = search_params->max_mismatches;
//    search_params->max_mismatches = caller_max_misms; // Full check
//    fmi_check_int_correctness(archive->index,search_params,matches,false,mpool);
//    search_params->max_mismatches = old_misms;
//  }
//  if (keepon && (check_level&CHECK_ALG_COMPLETENESS) && search_params->num_wildcards==0) { // Check completeness
//     const uint64_t old_misms = search_params->max_mismatches;
//    search_params->max_mismatches = (search_params->fast_mapping_degree>0) ?
//        matches->max_complete_stratum-1 : old_misms; /*caller_max_misms*/ // Full check
//    fmi_check_completeness(archive->index,search_params,matches,false,mpool);
//    search_params->max_mismatches = old_misms;
//  }
//  #endif
}
