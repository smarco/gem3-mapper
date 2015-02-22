/*
 * PROJECT: GEMMapper
 * FILE: archive_search.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_search.h"
#include "archive_select.h"

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
      &archive_search->forward_search_state,archive,
      &archive_search->search_actual_parameters);
  approximate_search_init(
      &archive_search->reverse_search_state,archive,
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
  PROF_START(GP_ARCHIVE_SEARCH_PREPARE_SEQUENCE);
  // Check the index characteristics & generate reverse-complement (if needed)
  if (archive_search->archive->indexed_complement) {
    archive_search->search_reverse = false;
  } else {
    if (archive_search->archive->filter_type == Iupac_dna) {
      sequence_generate_reverse_complement(&archive_search->sequence,&archive_search->rc_sequence);
    } else {
      sequence_generate_reverse(&archive_search->sequence,&archive_search->rc_sequence);
    }
    archive_search->search_reverse = !sequence_equals(&archive_search->sequence,&archive_search->rc_sequence);
  }
  // Generate the pattern(s)
  approximate_search_pattern_prepare(&archive_search->forward_search_state,&archive_search->sequence);
  if (archive_search->search_reverse) {
    approximate_search_pattern_prepare(&archive_search->reverse_search_state,&archive_search->rc_sequence);
  } else {
    approximate_search_pattern_clear(&archive_search->reverse_search_state);
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_PREPARE_SEQUENCE);
}
GEM_INLINE void archive_search_prepare_sequence_reverse(archive_search_t* const archive_search) {
  // Check if it was generated already or it's null // FIXME No possible use (dead-code)
  if (approximate_search_pattern_is_null(&archive_search->reverse_search_state)) {
    if (archive_search->archive->filter_type == Iupac_colorspace_dna) {
      sequence_generate_reverse_complement(&archive_search->sequence,&archive_search->rc_sequence);
    } else {
      sequence_generate_reverse(&archive_search->sequence,&archive_search->rc_sequence);
    }
    // Generate the pattern
    approximate_search_pattern_prepare(&archive_search->reverse_search_state,&archive_search->rc_sequence);
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
    archive_search_t* const archive_search,const bool verify_candidates,matches_t* const matches) {
  // Run the search (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  forward_asearch->verify_candidates = verify_candidates;
  forward_asearch->search_strand = Forward; // Configure forward search
  approximate_search(forward_asearch,matches);
  if (archive_search->search_reverse) {
    // Run the search (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    reverse_asearch->verify_candidates = verify_candidates;
    reverse_asearch->search_strand = Reverse; // Configure reverse search
    approximate_search(reverse_asearch,matches);
  }
}
GEM_INLINE void archive_search_generate_candidates(archive_search_t* const archive_search,matches_t* const matches) {
  ARCHIVE_SEARCH_CHECK(archive_search);
  PROF_START(GP_ARCHIVE_SEARCH_GENERATE_CANDIDATES);
  archive_search_continue(archive_search,false,NULL); // Run the search (stop before filtering)
  PROF_STOP(GP_ARCHIVE_SEARCH_GENERATE_CANDIDATES);
}
GEM_INLINE void archive_search_copy_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  PROF_START(GP_ARCHIVE_SEARCH_COPY_CANDIDATES);
  archive_t* const archive = archive_search->archive;
  // Add candidates (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  forward_asearch->bpm_buffer_offset = bpm_gpu_buffer_get_num_candidates(bpm_gpu_buffer);
  forward_asearch->bpm_buffer_candidates = filtering_candidates_bpm_buffer_add(
      &forward_asearch->filtering_candidates,archive,&forward_asearch->pattern,
      forward_asearch->search_strand,forward_asearch->search_actual_parameters,
      bpm_gpu_buffer,archive_search->mm_stack);
  if (archive_search->search_reverse) {
    // Add candidates (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    reverse_asearch->bpm_buffer_offset = bpm_gpu_buffer_get_num_candidates(bpm_gpu_buffer);
    reverse_asearch->bpm_buffer_candidates = filtering_candidates_bpm_buffer_add(
        &reverse_asearch->filtering_candidates,archive,&reverse_asearch->pattern,
        reverse_asearch->search_strand,reverse_asearch->search_actual_parameters,
        bpm_gpu_buffer,archive_search->mm_stack);
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
GEM_INLINE void archive_search_hold_verification_candidates(archive_search_t* const archive_search) {
  filtering_candidates_set_all_regions_pending(&archive_search->forward_search_state.filtering_candidates);
  if (archive_search->search_reverse) {
    filtering_candidates_set_all_regions_pending(&archive_search->reverse_search_state.filtering_candidates);
  }
}
GEM_INLINE void archive_search_release_verification_candidates(archive_search_t* const archive_search) {
  filtering_candidates_set_all_regions_unverified(&archive_search->forward_search_state.filtering_candidates);
  if (archive_search->search_reverse) {
    filtering_candidates_set_all_regions_unverified(&archive_search->reverse_search_state.filtering_candidates);
  }
}
GEM_INLINE void archive_search_finish_search(archive_search_t* const archive_search,matches_t* const matches) {
  // Run the search up to the end
  PROF_START(GP_ARCHIVE_SEARCH_FINISH_SEARCH);
  archive_search_continue(archive_search,true,matches);
  PROF_STOP(GP_ARCHIVE_SEARCH_FINISH_SEARCH);
  // Select matches
  archive_select_matches(archive_search,matches);
}
/*
 * SingleEnd Indexed Search (SE Online Approximate String Search)
 */
GEM_INLINE void archive_search_single_end(archive_search_t* const archive_search,matches_t* const matches) {
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
    if (lower_max_difference && archive_search->probe_strand) forward_asearch->stop_before_neighborhood_search = true;
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
        approximate_search(forward_asearch,matches);
      }
    }
  }
  // Select matches
  archive_select_matches(archive_search,matches);
  // Check matches
  select_parameters_t* const select_parameters = archive_search->select_parameters;
  if (select_parameters->check_correct || select_parameters->check_optimum) {
    search_parameters_t* const search_parameters = archive_search->search_actual_parameters.search_parameters;
    archive_check_matches_correct(
        archive_search->archive,search_parameters->alignment_model,
        &search_parameters->swg_penalties,&archive_search->sequence,
        matches,select_parameters->check_optimum,archive_search->mm_stack);
  }
  if (select_parameters->check_complete) {
    archive_check_matches_completness( // TODO
        archive_search->archive,&archive_search->sequence,matches,archive_search->mm_stack);
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_SE);
}
/*
 * PairedEnd Indexed Search (PE Online Approximate String Search)
 */
GEM_INLINE void archive_search_paired_end_extend_candidates(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,const sequence_end_t candidate_end) {
  PROF_START(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES);
  // Parameters
  search_parameters_t* const search_parameters = archive_search_end1->search_actual_parameters.search_parameters;
  archive_t* const archive = archive_search_end1->archive;
  // Extend in all possible concordant orientations
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  // Configure extension
  archive_search_t* candidate_archive_search;
  matches_t* candidate_matches;
  matches_t* extended_matches;
  search_actual_parameters_t* search_actual_parameters;
  uint64_t matches_found = 0;
  if (candidate_end==paired_end2) {
    // Extend towards end/2
    candidate_archive_search = archive_search_end2;
    search_actual_parameters = &candidate_archive_search->search_actual_parameters;
    extended_matches = matches_end1;
    candidate_matches = matches_end2;
  } else {
    // Extend towards end/1
    candidate_archive_search = archive_search_end1;
    search_actual_parameters = &candidate_archive_search->search_actual_parameters;
    extended_matches = matches_end2;
    candidate_matches = matches_end1;
  }
  filtering_candidates_t* const filtering_candidates = &candidate_archive_search->forward_search_state.filtering_candidates;
  pattern_t* const pattern = &candidate_archive_search->forward_search_state.pattern;
  text_collection_t* const text_collection = candidate_archive_search->text_collection;
  mm_stack_t* const mm_stack = candidate_archive_search->mm_stack;
  // Iterate over all matches of the extended end
  VECTOR_ITERATE(extended_matches->global_matches,extended_match,n,match_trace_t) {
    if (search_parameters->pair_orientation_FR == pair_orientation_concordant) {
      // Extend (filter nearby region)
      if (extended_match->strand==Forward) {
        matches_found += filtering_candidates_extend_match(
            filtering_candidates,archive,text_collection,extended_match,pattern,
            Reverse,true,search_actual_parameters,paired_matches,candidate_end,mm_stack);
      } else { // end2==Reverse
        extended_match->index_position = inverse_locator_map(archive->locator, // FIXME Hack?
            (uint8_t*)extended_match->sequence_name,extended_match->strand,extended_match->text_position);
        matches_found += filtering_candidates_extend_match(
            filtering_candidates,archive,text_collection,extended_match,pattern,
            Forward,false,search_actual_parameters,paired_matches,candidate_end,mm_stack);
      }
      // Decode & Pair found matches
      if (matches_found > 0) {
        PROF_ADD_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_FOUND,matches_found);
        match_trace_t* const mates_array =
            vector_get_mem(candidate_matches->global_matches,match_trace_t) +
            (vector_get_used(candidate_matches->global_matches) - matches_found);
        archive_select_extended_matches(candidate_archive_search,candidate_matches,mates_array,matches_found);
        // Pair with extended matches
        paired_matches_pair_match_with_mates(paired_matches,search_parameters,
            pair_orientation_concordant,extended_match,candidate_end,mates_array,matches_found);
      }
    }
  } // End Iteration
  PROF_STOP(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES);
}
GEM_INLINE void archive_search_paired_end_discard_filtering_regions(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2) {
  PROF_START(GP_ARCHIVE_SEARCH_PE_DISCARD_FILTERING_REGIONS);
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
      filtering_candidates_set_all_regions_pending(&forward_search_end1->filtering_candidates);
      filtering_candidates_set_all_regions_pending(&reverse_search_end1->filtering_candidates);
      filtering_candidates_set_all_regions_pending(&forward_search_end2->filtering_candidates);
      filtering_candidates_set_all_regions_pending(&reverse_search_end2->filtering_candidates);
      // Search for compatible pairs (Validates candidate-regions)
      const search_parameters_t* const parameters = archive_search_end1->search_actual_parameters.search_parameters;
      if (parameters->pair_orientation_FR != pair_orientation_invalid) {
        filtering_candidates_paired_regions_filtering(
            &forward_search_end1->filtering_candidates,&reverse_search_end2->filtering_candidates,
            parameters->min_template_length,parameters->max_template_length,false);
      }
      if (parameters->pair_orientation_RF != pair_orientation_invalid) {
        filtering_candidates_paired_regions_filtering(
            &reverse_search_end1->filtering_candidates,&forward_search_end2->filtering_candidates,
            parameters->min_template_length,parameters->max_template_length,false);
      }
      if (parameters->pair_orientation_FF != pair_orientation_invalid) {
        filtering_candidates_paired_regions_filtering(
            &forward_search_end1->filtering_candidates,&forward_search_end2->filtering_candidates,
            parameters->min_template_length,parameters->max_template_length,true);
      }
      if (parameters->pair_orientation_RR != pair_orientation_invalid) {
        filtering_candidates_paired_regions_filtering(
            &reverse_search_end1->filtering_candidates,&reverse_search_end2->filtering_candidates,
            parameters->min_template_length,parameters->max_template_length,true);
      }
      // Profile (Count the number of discarded filtering regions)
      PROF_BLOCK() {
        uint64_t num_regions = 0, num_regions_discarded = 0;
        num_regions_discarded += filtering_candidates_count_candidate_regions(
            &forward_search_end1->filtering_candidates,filtering_region_pending);
        num_regions_discarded += filtering_candidates_count_candidate_regions(
            &reverse_search_end1->filtering_candidates,filtering_region_pending);
        num_regions_discarded += filtering_candidates_count_candidate_regions(
            &forward_search_end2->filtering_candidates,filtering_region_pending);
        num_regions_discarded += filtering_candidates_count_candidate_regions(
            &reverse_search_end2->filtering_candidates,filtering_region_pending);
        num_regions += filtering_candidates_get_num_candidate_regions(&forward_search_end1->filtering_candidates);
        num_regions += filtering_candidates_get_num_candidate_regions(&reverse_search_end1->filtering_candidates);
        num_regions += filtering_candidates_get_num_candidate_regions(&forward_search_end2->filtering_candidates);
        num_regions += filtering_candidates_get_num_candidate_regions(&reverse_search_end2->filtering_candidates);
        PROF_ADD_COUNTER(GP_ARCHIVE_SEARCH_PE_DISCARD_FILTERING_REGIONS_NOT_CONCORDANT,num_regions_discarded);
        PROF_ADD_COUNTER(GP_ARCHIVE_SEARCH_PE_DISCARD_FILTERING_REGIONS_TOTAL,num_regions);
      }
    }
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_PE_DISCARD_FILTERING_REGIONS);
}
GEM_INLINE bool archive_search_paired_end_fulfilled(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  return vector_get_used(paired_matches->concordant_matches) > 0;
}
GEM_INLINE bool archive_search_paired_end_use_extension(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    const sequence_end_t candidate_end,const paired_matches_t* const paired_matches) {
  search_parameters_t* const search_parameters = archive_search_end1->search_actual_parameters.search_parameters;
  // Check mapping parameters
  if (search_parameters->paired_mapping_mode == paired_mapping_map_both_ends ||
      search_parameters->paired_mapping_mode == paired_mapping_paired_filtering) return false;
  // Check the number of samples to derive the expected insert size
  const uint64_t num_samples = COUNTER_GET_NUM_SAMPLES(&paired_matches->unique_template_size);
  if (num_samples < search_parameters->min_unique_pair_samples) return false;
  // Check the number of candidates
  if (candidate_end==paired_end2) { // Extend end/1 => TODO dynamic choice (learning & balancing)
    // if (archive_search_end2->sequence.read.length > TODO) return false;
    const uint64_t candidates_end1 = archive_search_get_search_canditates(archive_search_end1);
    // TODO Substitute @max_extendable_candidates with a @COUNTER of number of candidates to check per read
    return candidates_end1>0 && candidates_end1<=search_parameters->max_extendable_candidates;
  } else { // extended_end==ARCHIVE_SEARCH_END2
    const uint64_t candidates_end2 = archive_search_get_search_canditates(archive_search_end2);
    return candidates_end2>0 && candidates_end2<=search_parameters->max_extendable_candidates;
  }
}
/*
 * PairedEnd Indexed Search (PE Online Approximate String Search)
 */
GEM_INLINE void archive_search_paired_end_continue(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  PROF_START(GP_ARCHIVE_SEARCH_PE);
  // Matches
  matches_t* const matches_end1 = (paired_matches!=NULL) ? paired_matches->matches_end1 : NULL;
  matches_t* const matches_end2 = (paired_matches!=NULL) ? paired_matches->matches_end2 : NULL;
  // Callback (switch to proper search stage)
  bool try_extending_end1 = false, try_extending_end2 = false;
  archive_search_paired_end_callback:
  switch (archive_search_end1->pe_search_state) {
    case archive_search_pe_begin: // Beginning of the search (Init)
      archive_search_reset(archive_search_end1,sequence_get_length(&archive_search_end1->sequence));
      archive_search_reset(archive_search_end2,sequence_get_length(&archive_search_end2->sequence));
    // No break
    case archive_search_pe_search_end1:
      // Generate candidates (End/1)
      archive_search_generate_candidates(archive_search_end1,matches_end1);
      // Try extension (End/1)
      try_extending_end1 = archive_search_paired_end_use_extension(
          archive_search_end1,archive_search_end2,paired_end2,paired_matches);
      if (try_extending_end1) {
        // Finish search (End/1)
        archive_search_finish_search(archive_search_end1,matches_end1);
        // Try extending (End/1)
        PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END1);
        archive_search_end1->pe_search_state = archive_search_pe_extend_end1;
        // Callback
        goto archive_search_paired_end_callback;
      }
    // No break
    case archive_search_pe_search_end2:
      // Generate candidates (End/2)
      archive_search_generate_candidates(archive_search_end2,matches_end2);
      // Verify candidates. Combined filtering
      //   // TODO Activate => archive_search_paired_end_discard_filtering_regions
      //   // TODO Finish search considering all candidates, even the discarded ones
      // Try extension (End/2)
      try_extending_end2 = archive_search_paired_end_use_extension(
          archive_search_end1,archive_search_end2,paired_end1,paired_matches);
      if (try_extending_end2) {
        // Put on-hold search (End/1)
        archive_search_hold_verification_candidates(archive_search_end1);
        // Finish search (End/2)
        archive_search_finish_search(archive_search_end2,matches_end2);
        // Try extending (End/2)
        PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END2);
        archive_search_end1->pe_search_state = archive_search_pe_extend_end2;
        // Callback
        goto archive_search_paired_end_callback;
      }
    // No break
    case archive_search_pe_verify_both_ends: { // Verify candidates for both ends
      search_parameters_t* const search_parameters = archive_search_end1->search_actual_parameters.search_parameters;
      if (search_parameters->paired_mapping_mode == paired_mapping_paired_filtering) {
        archive_search_paired_end_discard_filtering_regions(archive_search_end1,archive_search_end2);
      }
      archive_search_finish_search(archive_search_end1,matches_end1); // Finish search (End/1)
      archive_search_finish_search(archive_search_end2,matches_end2); // Finish search (End/2)
      archive_search_end1->pe_search_state = archive_search_pe_both_ends_verified;
      goto archive_search_paired_end_callback;
    }
    // No break
    case archive_search_pe_extend_end1: // Extend (End/1)
      archive_search_paired_end_extend_candidates(archive_search_end1,archive_search_end2,paired_matches,paired_end2);
    // No break
    case archive_search_pe_extended_end1: // End1 extended
      if (archive_search_paired_end_fulfilled(archive_search_end1,archive_search_end2,paired_matches)) {
        paired_matches->max_complete_stratum = 0; // No algorithmic guarantee
        PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END1_SUCCESS);
        archive_search_end1->pe_search_state = archive_search_pe_end;
        PROF_STOP(GP_ARCHIVE_SEARCH_PE); return; // Enough
      }
      /*
       * Extension failed because
       *  (1) Insert size is beyond the expected distribution (Any insert size filtering/search must be discarded)
       *  (2) Sensitivity of the candidates/end1 search is not enough
       */
      // Archive search (End/2)
      archive_search_single_end(archive_search_end2,matches_end2);
      archive_search_end1->pe_search_state = archive_search_pe_both_ends_verified;
      goto archive_search_paired_end_callback; // Callback
      break;
    case archive_search_pe_extend_end2: // Extend (End/2)
      archive_search_paired_end_extend_candidates(archive_search_end1,archive_search_end2,paired_matches,paired_end1);
    // No break
    case archive_search_pe_extended_end2: // End2 extended
      if (archive_search_paired_end_fulfilled(archive_search_end1,archive_search_end2,paired_matches)) {
        PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END2_SUCCESS);
        paired_matches->max_complete_stratum = 0; // No algorithmic guarantee
        PROF_STOP(GP_ARCHIVE_SEARCH_PE); return; // Enough
      }
      /*
       * Extension failed because
       *  (1) Insert size is beyond the expected distribution (Any insert size filtering/search must be discarded)
       *  (2) Sensitivity of the candidates/end2 search is not enough
       */
      // Finish search (End/1)
      archive_search_release_verification_candidates(archive_search_end1);
      archive_search_finish_search(archive_search_end1,matches_end1);
    // No break
    case archive_search_pe_both_ends_verified: { // Candidates of both ends had been verified
      // Check number of matches found
      const uint64_t num_matches_end1 = vector_get_used(paired_matches->matches_end1->global_matches);
      const uint64_t num_matches_end2 = vector_get_used(paired_matches->matches_end2->global_matches);
//      if (!try_extending_end1 && num_matches_end1==0) {
//        // Rescue last resort // TODO
//      }
//      if (!try_extending_end2 && num_matches_end2==0) {
//        // Rescue last resort // TODO
//      }
      if (num_matches_end1 > 0 || num_matches_end2 > 0) {
        // Pair matches (Cross-link matches from both ends)
        paired_matches_find_pairs(paired_matches,archive_search_end1->search_actual_parameters.search_parameters);
        // Find discordant (if required)
        // paired_matches_find_discordant_pairs(paired_matches,search_parameters); // TODO Activate
        // Select matches
        archive_select_paired_matches(archive_search_end1,archive_search_end2,paired_matches);
      }
      archive_search_end1->pe_search_state = archive_search_pe_end;
    }
    // No break
    case archive_search_pe_end: break; // End of the current workflow
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_PE);
}
GEM_INLINE void archive_search_paired_end(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  archive_search_end1->pe_search_state = archive_search_pe_begin;
  archive_search_end2->pe_search_state = archive_search_pe_begin;
  archive_search_paired_end_continue(archive_search_end1,archive_search_end2,paired_matches);
}

