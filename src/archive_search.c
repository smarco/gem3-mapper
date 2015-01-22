/*
 * PROJECT: GEMMapper
 * FILE: archive_search.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_search.h"
#include "archive_select.h"

/*
 * Constants
 */
#define ARCHIVE_SEARCH_END1 0
#define ARCHIVE_SEARCH_END2 1

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
  sequence_init(&archive_search->sequence);
  sequence_init(&archive_search->rc_sequence);
  // Approximate Search
  archive_search->search_actual_parameters.search_parameters = search_parameters;
  archive_search->select_parameters = select_parameters;
  approximate_search_init(
      &archive_search->forward_search_state,
      archive->locator,archive->graph,archive->enc_text,archive->fm_index,
      &archive_search->search_actual_parameters);
  approximate_search_init(
      &archive_search->reverse_search_state,
      archive->locator,archive->graph,archive->enc_text,archive->fm_index,
      &archive_search->search_actual_parameters);
  // Archive search control (Flow control) [DEFAULTS]
  archive_search->probe_strand = true;
  archive_search->search_reverse = !archive->indexed_complement;
  // Return
  return archive_search;
}
GEM_INLINE void archive_search_configure(
    archive_search_t* const archive_search,mm_search_t* const mm_search) {
  // Text-Collection
  archive_search->text_collection = &mm_search->text_collection;
  // Clear F/R search states
  approximate_search_configure(
      &archive_search->forward_search_state,&mm_search->text_collection,
      &mm_search->interval_set,mm_search->mm_stack);
  approximate_search_configure(
      &archive_search->reverse_search_state,&mm_search->text_collection,
      &mm_search->interval_set,mm_search->mm_stack);
  // MM
  archive_search->mm_stack = mm_search->mm_stack;
}
GEM_INLINE void archive_search_prepare_sequence(archive_search_t* const archive_search) {
  // Check the index characteristics & generate reverse-complement (if needed)
  if (archive_search->archive->indexed_complement) {
    archive_search->search_reverse = false;
  } else {
    if (archive_search->archive->filter_type == Iupac_colorspace_dna) {
      sequence_generate_reverse_complement(&archive_search->sequence,&archive_search->rc_sequence);
    } else {
      sequence_generate_reverse(&archive_search->sequence,&archive_search->rc_sequence);
    }
    archive_search->search_reverse = !sequence_equals(&archive_search->sequence,&archive_search->rc_sequence);
  }
  // Generate the pattern(s)
  approximate_search_prepare_pattern(
      &archive_search->forward_search_state,&archive_search->sequence);
  if (archive_search->search_reverse) {
    approximate_search_prepare_pattern(
        &archive_search->reverse_search_state,&archive_search->rc_sequence);
  }
}
GEM_INLINE void archive_search_reset(archive_search_t* const archive_search,const uint64_t sequence_length) {
  // Instantiate parameters actual-values
  approximate_search_instantiate_values(&archive_search->search_actual_parameters,sequence_length);
  // Clear F/R search states
  approximate_search_reset(&archive_search->forward_search_state);
  approximate_search_reset(&archive_search->reverse_search_state);
  // Prepare for sequence
  archive_search_prepare_sequence(archive_search);
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
 * Archive Search [Accessors]
 */
GEM_INLINE sequence_t* archive_search_get_sequence(const archive_search_t* const archive_search) {
  return (sequence_t*)&archive_search->sequence;
}
GEM_INLINE uint64_t archive_search_get_search_canditates(const archive_search_t* const archive_search) {
  if (archive_search->archive->indexed_complement) {
    return approximate_search_get_num_potential_candidates(&archive_search->forward_search_state);
  } else {
    return approximate_search_get_num_potential_candidates(&archive_search->forward_search_state) +
           approximate_search_get_num_potential_candidates(&archive_search->reverse_search_state);
  }
}
/*
 * Archive Search (Step-wise Search building-blocks)
 */
GEM_INLINE void archive_search_continue(
    archive_search_t* const archive_search,
    const approximate_search_state_t stop_before_state,matches_t* const matches) {
  // Run the search (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  forward_asearch->stop_before_state = stop_before_state; // Stop before state
  forward_asearch->search_strand = Forward; // Configure forward search
  approximate_search(forward_asearch,matches);
  if (archive_search->search_reverse) {
    // Run the search (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    reverse_asearch->stop_before_state = stop_before_state; // Stop before state
    reverse_asearch->search_strand = Reverse; // Configure reverse search
    approximate_search(reverse_asearch,matches);
  }
}
GEM_INLINE void archive_search_generate_candidates(archive_search_t* const archive_search) {
  ARCHIVE_SEARCH_CHECK(archive_search);
  PROF_START(GP_ARCHIVE_SEARCH_GENERATE_CANDIDATES);
  // Reset initial values (Prepare pattern(s), instantiate parameters values, ...)
  archive_search_reset(archive_search,sequence_get_length(&archive_search->sequence));
  // Run the search (stop before filtering)
  archive_search_continue(archive_search,asearch_verify_candidates,NULL);
  PROF_STOP(GP_ARCHIVE_SEARCH_GENERATE_CANDIDATES);
}
GEM_INLINE void archive_search_copy_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  PROF_START(GP_ARCHIVE_SEARCH_COPY_CANDIDATES);
  const archive_t* const archive = archive_search->archive;
  // Add candidates (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  forward_asearch->bpm_buffer_offset = bpm_gpu_buffer_get_num_candidates(bpm_gpu_buffer);
  forward_asearch->bpm_buffer_candidates = filtering_candidates_bpm_buffer_add(
      &forward_asearch->filtering_candidates,archive->locator,archive->fm_index,archive->enc_text,
      &forward_asearch->pattern,forward_asearch->search_strand,
      forward_asearch->search_actual_parameters,bpm_gpu_buffer,archive_search->mm_stack);
  if (archive_search->search_reverse) {
    // Add candidates (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    reverse_asearch->bpm_buffer_offset = bpm_gpu_buffer_get_num_candidates(bpm_gpu_buffer);
    reverse_asearch->bpm_buffer_candidates = filtering_candidates_bpm_buffer_add(
        &reverse_asearch->filtering_candidates,archive->locator,archive->fm_index,archive->enc_text,
        &reverse_asearch->pattern,reverse_asearch->search_strand,
        reverse_asearch->search_actual_parameters,bpm_gpu_buffer,archive_search->mm_stack);
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_COPY_CANDIDATES);
}
GEM_INLINE void archive_search_retrieve_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer,matches_t* const matches) {
  PROF_START(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES);
  // Verified candidates (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  approximate_search_verify_using_bpm_buffer(forward_asearch,matches,bpm_gpu_buffer,
      forward_asearch->bpm_buffer_offset,forward_asearch->bpm_buffer_offset+forward_asearch->bpm_buffer_candidates);
  if (archive_search->search_reverse) {
    // Verified candidates (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    approximate_search_verify_using_bpm_buffer(reverse_asearch,matches,bpm_gpu_buffer,
        reverse_asearch->bpm_buffer_offset,reverse_asearch->bpm_buffer_offset+reverse_asearch->bpm_buffer_candidates);
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES);
}
/*
 * SingleEnd Indexed Search (SE Online Approximate String Search)
 */
GEM_INLINE void archive_search_single_end(archive_search_t* const archive_search,matches_t* const matches) {
  ARCHIVE_SEARCH_CHECK(archive_search);
  PROF_START(GP_ARCHIVE_SEARCH_SE);
  // Reset initial values (Prepare pattern(s), instantiate parameters values, ...)
  archive_search_reset(archive_search,sequence_get_length(&archive_search->sequence));
  // Search the pattern(s)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  if (!archive_search->search_reverse) {
    // Compute the full search
    forward_asearch->search_strand = Forward;
    approximate_search(forward_asearch,matches);
  } else {
    // Configure search stage to stop at
    const search_actual_parameters_t* const actual_parameters = &archive_search->search_actual_parameters;
    const bool lower_max_difference = actual_parameters->complete_strata_after_best_nominal < forward_asearch->max_differences;
    if (lower_max_difference && archive_search->probe_strand) {
      forward_asearch->stop_before_state = asearch_neighborhood;
    }
    // Run the search (FORWARD)
    forward_asearch->search_strand = Forward; // Configure forward search
    approximate_search(forward_asearch,matches);
    // Check the number of matches & keep searching
    if (!forward_asearch->max_matches_reached) {
      // Keep on searching
      approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
      // Run the search (REVERSE)
      reverse_asearch->search_strand = Reverse; // Configure reverse search
      approximate_search(reverse_asearch,matches);
      // Resume forward search (if not completed before)
      if (forward_asearch->search_state != asearch_end && !forward_asearch->max_matches_reached) {
        forward_asearch->stop_before_state = asearch_end;
        approximate_search(forward_asearch,matches);
      }
    }
  }
  // Select matches
  archive_select_matches(archive_search,matches);
  // Check matches // TODO
  const uint64_t check_matches_mask = archive_search->select_parameters->check_matches_mask;
  if (check_matches_mask & check_correctness) {
    archive_check_matches_correctness(
        archive_search->archive,&archive_search->sequence,matches,archive_search->mm_stack);
  } else if (check_matches_mask & check_optimum) {
    archive_check_matches_optimum(
        archive_search->archive,&archive_search->sequence,matches,archive_search->mm_stack);
  } else if (check_matches_mask & check_completness) {
    archive_check_matches_completness(
        archive_search->archive,&archive_search->sequence,matches,archive_search->mm_stack);
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_SE);
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
/*
 * PairedEnd Indexed Search (PE Online Approximate String Search)
 */
GEM_INLINE void archive_search_paired_end_extend_candidates(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,const uint64_t extended_end) {
  /* TODO */
  /* TODO */
  /* TODO */
  /* TODO */
  /* TODO */
  /* TODO */
  /* TODO */
  /* TODO */
  /* TODO */
  /* TODO */
  /* TODO */
//  PROF_START(GP_ARCHIVE_SEARCH_COPY_CANDIDATES);
//  const archive_t* const archive = archive_search_end1->archive;
//  const mm_stack_t* const mm_stack = archive_search_end1->mm_stack;
//  // Add candidates (FORWARD)
//  approximate_search_t* asearch_candidate;
//  approximate_search_t* asearch_extend;
//  if (extended_end==ARCHIVE_SEARCH_END1) {
//    asearch_extend = &archive_search_end1->forward_search_state;
//    asearch_candidate = &archive_search_end2->forward_search_state;
//  } else {
//    asearch_candidate = &archive_search_end1->forward_search_state;
//    asearch_extend = &archive_search_end2->forward_search_state;
//  }
//  filtering_candidates_extend_matches(
//      asearch_candidate->filtering_candidates,archive,Forward,extended_end,
//      &asearch_candidate->pattern,asearch_candidate->search_actual_parameters,
//      &asearch_extend->pattern,asearch_extend->search_actual_parameters,paired_matches,mm_stack);
//  if (archive_search_end1->search_reverse) {
//    // Add candidates (REVERSE)
//    if (extended_end==ARCHIVE_SEARCH_END1) {
//      asearch_extend = &archive_search_end1->reverse_search_state;
//      asearch_candidate = &archive_search_end2->reverse_search_state;
//    } else {
//      asearch_candidate = &archive_search_end1->reverse_search_state;
//      asearch_extend = &archive_search_end2->reverse_search_state;
//    }
//    filtering_candidates_extend_matches(
//        asearch_candidate->filtering_candidates,archive,Reverse,extended_end,
//        &asearch_candidate->pattern,asearch_candidate->search_actual_parameters,
//        &asearch_extend->pattern,asearch_extend->search_actual_parameters,paired_matches,mm_stack);
//  }
//  PROF_STOP(GP_ARCHIVE_SEARCH_COPY_CANDIDATES);
}
GEM_INLINE void archive_search_paired_process_candidates(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2) {
  PROF_START(GP_PAIRED_DISCARD_FILTERING_REGIONS);
  // Check the index type
  if (!archive_search_end1->archive->indexed_complement) { // TODO For RC-indexes
    // Reduce the number of search candidates (paired selection of regions)
    approximate_search_t* const forward_search_end1 = &archive_search_end1->forward_search_state;
    approximate_search_t* const reverse_search_end1 = &archive_search_end1->reverse_search_state;
    approximate_search_t* const forward_search_end2 = &archive_search_end2->forward_search_state;
    approximate_search_t* const reverse_search_end2 = &archive_search_end2->reverse_search_state;
    if (forward_search_end1->search_state == asearch_verify_candidates &&
        reverse_search_end1->search_state == asearch_verify_candidates &&
        forward_search_end2->search_state == asearch_verify_candidates &&
        reverse_search_end2->search_state == asearch_verify_candidates) {
      // Initialize (Invalidates all candidate-regions as lacking mate)
      approximate_search_init_paired_filtering(forward_search_end1);
      approximate_search_init_paired_filtering(reverse_search_end1);
      approximate_search_init_paired_filtering(forward_search_end2);
      approximate_search_init_paired_filtering(reverse_search_end2);
      // Search for compatible pairs (Validates candidate-regions)
      const search_parameters_t* const parameters = archive_search_end1->search_actual_parameters.search_parameters;
      if (parameters->pair_orientation_FR != pair_orientation_invalid) {
        approximate_search_verify_paired_process_candidates(forward_search_end1,reverse_search_end2,
            parameters->min_template_length,parameters->max_template_length,false);
      }
      if (parameters->pair_orientation_RF != pair_orientation_invalid) {
        approximate_search_verify_paired_process_candidates(reverse_search_end1,forward_search_end2,
            parameters->min_template_length,parameters->max_template_length,false);
      }
      if (parameters->pair_orientation_FF != pair_orientation_invalid) {
        approximate_search_verify_paired_process_candidates(forward_search_end1,forward_search_end2,
            parameters->min_template_length,parameters->max_template_length,true);
      }
      if (parameters->pair_orientation_RR != pair_orientation_invalid) {
        approximate_search_verify_paired_process_candidates(reverse_search_end1,reverse_search_end2,
            parameters->min_template_length,parameters->max_template_length,true);
      }
      // Profile (Count the number of discarded filtering regions)
      PROF_BLOCK() {
        uint64_t num_regions = 0, num_regions_discarded = 0;
        num_regions_discarded += filtering_candidates_count_candidate_regions(
            &forward_search_end1->filtering_candidates,filtering_region_unpaired);
        num_regions_discarded += filtering_candidates_count_candidate_regions(
            &reverse_search_end1->filtering_candidates,filtering_region_unpaired);
        num_regions_discarded += filtering_candidates_count_candidate_regions(
            &forward_search_end2->filtering_candidates,filtering_region_unpaired);
        num_regions_discarded += filtering_candidates_count_candidate_regions(
            &reverse_search_end2->filtering_candidates,filtering_region_unpaired);
        num_regions += filtering_candidates_get_num_candidate_regions(&forward_search_end1->filtering_candidates);
        num_regions += filtering_candidates_get_num_candidate_regions(&reverse_search_end1->filtering_candidates);
        num_regions += filtering_candidates_get_num_candidate_regions(&forward_search_end2->filtering_candidates);
        num_regions += filtering_candidates_get_num_candidate_regions(&reverse_search_end2->filtering_candidates);
        PROF_ADD_COUNTER(GP_PAIRED_NOT_CONCORDANT_FILTERING_REGIONS,num_regions_discarded);
        PROF_ADD_COUNTER(GP_PAIRED_ALL_FILTERING_REGIONS,num_regions);
      }
    }
  }
  PROF_STOP(GP_PAIRED_DISCARD_FILTERING_REGIONS);
}
GEM_INLINE bool archive_search_paired_end_fulfilled(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  /* TODO */
  return false;
}
GEM_INLINE bool archive_search_paired_end_use_extension(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,const uint64_t extended_end) {
  // search_parameters_t* const search_parameters = archive_search_end1->search_actual_parameters.search_parameters;
  // TODO dynamic choice (learning & balancing)
  return false;
//  if (extended_end==ARCHIVE_SEARCH_END1) {
//    const uint64_t candidates_end1 = archive_search_get_search_canditates(archive_search_end1);
//    return candidates_end1>0 && candidates_end1<=search_parameters->max_extendable_candidates;
//  } else { // extended_end==ARCHIVE_SEARCH_END2
//    const uint64_t candidates_end2 = archive_search_get_search_canditates(archive_search_end2);
//    return candidates_end2>0 && candidates_end2<=search_parameters->max_extendable_candidates;
//  }
}
GEM_INLINE void archive_search_paired_end(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  search_parameters_t* const search_parameters = archive_search_end1->search_actual_parameters.search_parameters;
  const bool map_both_ends = search_parameters->map_both_ends;
  /*
   * Init
   */
  archive_search_reset(archive_search_end1,sequence_get_length(&archive_search_end1->sequence));
  archive_search_reset(archive_search_end2,sequence_get_length(&archive_search_end2->sequence));
  /*
   * Try combined filtering
   */
  // Generate candidates (End/1)
  archive_search_continue(archive_search_end1,asearch_verify_candidates,paired_matches_end1(paired_matches));
  // Try extension (End/1)
  if (!map_both_ends && archive_search_paired_end_use_extension(archive_search_end1,archive_search_end2,ARCHIVE_SEARCH_END1)) {
    archive_search_paired_end_extend_candidates(archive_search_end1,archive_search_end2,paired_matches,ARCHIVE_SEARCH_END1);
    if (archive_search_paired_end_fulfilled(archive_search_end1,archive_search_end2,paired_matches)) return; // Enough
  } // TODO Update MCS
  // Generate candidates (End/2)
  archive_search_continue(archive_search_end2,asearch_verify_candidates,paired_matches_end2(paired_matches));
  // Try extension (End/2)
  if (!map_both_ends && archive_search_paired_end_use_extension(archive_search_end1,archive_search_end2,ARCHIVE_SEARCH_END2)) {
    archive_search_paired_end_extend_candidates(archive_search_end1,archive_search_end2,paired_matches,ARCHIVE_SEARCH_END2);
    if (archive_search_paired_end_fulfilled(archive_search_end1,archive_search_end2,paired_matches)) return; // Enough
  } // TODO Update MCS
  /*
   * Verify candidates
   */
//  // Combined filtering
//  archive_search_paired_process_candidates(archive_search_end1,archive_search_end2); // FIXME UncommentME
  /*
   * Resume the search of both ends & pair
   */
  // Finish search (End/1)
  archive_search_continue(archive_search_end1,asearch_end,paired_matches_end1(paired_matches));
  archive_select_matches(archive_search_end1,paired_matches_end1(paired_matches)); // Select Matches (End/1)
  // Finish search (End/2)
  archive_search_continue(archive_search_end2,asearch_end,paired_matches_end2(paired_matches));
  archive_select_matches(archive_search_end2,paired_matches_end2(paired_matches)); // Select Matches (End/2)
  // Pair matches (Cross-link matches from both ends)
  paired_matches_find_pairs(paired_matches,search_parameters);
  // Find discordant (if required)
  // paired_matches_find_discordant_pairs(paired_matches,search_parameters); // TODO Activate
  /*
   * Post-processing
   */
  // Select matches
  archive_select_paired_matches(archive_search_end1,archive_search_end2,paired_matches);
}

