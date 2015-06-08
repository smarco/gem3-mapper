/*
 * PROJECT: GEMMapper
 * FILE: archive_search.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_search.h"
#include "archive_select.h"
#include "archive_score.h"

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
GEM_INLINE void archive_search_configure_(
    archive_search_t* const archive_search,const sequence_end_t sequence_end,
    mm_search_t* const mm_search) {
  // Configure F/R search states
  switch (sequence_end) {
    case single_end:
    case paired_end1:
      approximate_search_configure(
          &archive_search->forward_search_state,&mm_search->filtering_candidates_forward_end1,
          &mm_search->text_collection,&mm_search->interval_set,mm_search->mm_stack);
      approximate_search_configure(
          &archive_search->reverse_search_state,&mm_search->filtering_candidates_reverse_end1,
          &mm_search->text_collection,&mm_search->interval_set,mm_search->mm_stack);
      break;
    case paired_end2:
      approximate_search_configure(
          &archive_search->forward_search_state,&mm_search->filtering_candidates_forward_end2,
          &mm_search->text_collection,&mm_search->interval_set,mm_search->mm_stack);
      approximate_search_configure(
          &archive_search->reverse_search_state,&mm_search->filtering_candidates_reverse_end2,
          &mm_search->text_collection,&mm_search->interval_set,mm_search->mm_stack);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Text-Collection
  archive_search->text_collection = &mm_search->text_collection;
  // Stats
  archive_search->mapper_stats = mm_search->mapper_stats;
  // MM
  archive_search->mm_stack = mm_search->mm_stack;
}
GEM_INLINE void archive_search_configure(archive_search_t* const archive_search,mm_search_t* const mm_search) {
  archive_search_configure_(archive_search,single_end,mm_search);
}
GEM_INLINE void archive_search_pe_configure(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    mm_search_t* const mm_search) {
  archive_search_configure_(archive_search_end1,paired_end1,mm_search);
  archive_search_configure_(archive_search_end2,paired_end2,mm_search);
}
GEM_INLINE void archive_search_prepare_sequence(archive_search_t* const archive_search) {
  PROF_START(GP_ARCHIVE_SEARCH_PREPARE_SEQUENCE);
  // Check the index characteristics & generate reverse-complement (if needed)
  if (archive_search->archive->indexed_complement) {
    archive_search->emulate_rc_search = false;
  } else {
    if (archive_search->archive->filter_type == Iupac_dna) {
      sequence_generate_reverse_complement(&archive_search->sequence,&archive_search->rc_sequence);
    } else {
      sequence_generate_reverse(&archive_search->sequence,&archive_search->rc_sequence);
    }
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
  PROF_STOP(GP_ARCHIVE_SEARCH_PREPARE_SEQUENCE);
}
GEM_INLINE void archive_search_reset(archive_search_t* const archive_search) {
  PROF_START(GP_ARCHIVE_SEARCH_SE_INIT);
  // Instantiate parameters actual-values
  const uint64_t sequence_length = sequence_get_length(&archive_search->sequence);
  approximate_search_instantiate_values(&archive_search->as_parameters,sequence_length);
  // Prepare for sequence
  archive_search_prepare_sequence(archive_search);
  // Clear F/R search states
  approximate_search_reset(&archive_search->forward_search_state);
  if (archive_search->emulate_rc_search) {
    approximate_search_reset(&archive_search->reverse_search_state);
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_SE_INIT);
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
    return approximate_search_get_num_filtering_candidates(&archive_search->forward_search_state);
  } else {
    return approximate_search_get_num_filtering_candidates(&archive_search->forward_search_state) +
           approximate_search_get_num_filtering_candidates(&archive_search->reverse_search_state);
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
GEM_INLINE bool archive_search_finished(const archive_search_t* const archive_search) {
  if (archive_search->archive->indexed_complement) {
    return archive_search->forward_search_state.search_state == asearch_end;
  } else {
    return archive_search->forward_search_state.search_state == asearch_end &&
           archive_search->reverse_search_state.search_state == asearch_end;
  }
}
/*
 * SingleEnd Indexed Search (Step-wise)
 */
GEM_INLINE void archive_search_continue(
    archive_search_t* const archive_search,const bool verify_candidates,matches_t* const matches) {
  // Run the search (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  forward_asearch->verify_candidates = verify_candidates;
  approximate_search(forward_asearch,matches); // Forward search
  if (archive_search->emulate_rc_search) {
    // Run the search (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    reverse_asearch->verify_candidates = verify_candidates;
    approximate_search(reverse_asearch,matches); // Reverse emulated-search
  }
}
GEM_INLINE void archive_search_generate_candidates(archive_search_t* const archive_search) {
  ARCHIVE_SEARCH_CHECK(archive_search);
  // Reset initial values (Prepare pattern(s), instantiate parameters values, ...)
  archive_search_reset(archive_search);
  // Run the search (stop before filtering)
  PROF_START(GP_ARCHIVE_SEARCH_SE_GENERATE_CANDIDATES);
  archive_search_continue(archive_search,false,NULL);
  PROF_STOP(GP_ARCHIVE_SEARCH_SE_GENERATE_CANDIDATES);
}
GEM_INLINE void archive_search_verify_candidates(archive_search_t* const archive_search,matches_t* const matches) {
  PROF_START(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES);
  // Verify candidates (FORWARD)
  approximate_search_verify(&archive_search->forward_search_state,matches);
  if (archive_search->emulate_rc_search) {
    // Verify candidates (REVERSE)
    approximate_search_verify(&archive_search->reverse_search_state,matches);
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES);
}
GEM_INLINE void archive_search_finish_search(archive_search_t* const archive_search,matches_t* const matches) {
  // Run the search up to the end
  PROF_START(GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH);
  archive_search_continue(archive_search,true,matches);
  PROF_STOP(GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH);
  // Select matches
  archive_select_matches(archive_search,true,true,true,matches);
}
GEM_INLINE void archive_search_copy_candidates(
    archive_search_t* const archive_search,bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  PROF_START(GP_ARCHIVE_SEARCH_COPY_CANDIDATES);
  // Add candidates (FORWARD)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  forward_asearch->bpm_buffer_offset = bpm_gpu_buffer_get_num_candidates(bpm_gpu_buffer);
  forward_asearch->bpm_buffer_candidates = filtering_candidates_bpm_buffer_add(
      forward_asearch->filtering_candidates,&forward_asearch->pattern,bpm_gpu_buffer);
  if (archive_search->emulate_rc_search) {
    // Add candidates (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    reverse_asearch->bpm_buffer_offset = bpm_gpu_buffer_get_num_candidates(bpm_gpu_buffer);
    reverse_asearch->bpm_buffer_candidates = filtering_candidates_bpm_buffer_add(
        reverse_asearch->filtering_candidates,&reverse_asearch->pattern,bpm_gpu_buffer);
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
  if (archive_search->emulate_rc_search) {
    // Verified candidates (REVERSE)
    approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
    approximate_search_verify_using_bpm_buffer(reverse_asearch,matches,bpm_gpu_buffer,
        reverse_asearch->bpm_buffer_offset,reverse_asearch->bpm_buffer_offset+reverse_asearch->bpm_buffer_candidates);
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES);
}
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
/*
 * SingleEnd Indexed Search (SE Online Approximate String Search)
 */
GEM_INLINE void archive_search_single_end(archive_search_t* const archive_search,matches_t* const matches) {
  PROF_START(GP_ARCHIVE_SEARCH_SE);
  // Reset initial values (Prepare pattern(s), instantiate parameters values, ...)
  archive_search_reset(archive_search);
  // Search the pattern(s)
  approximate_search_t* const forward_asearch = &archive_search->forward_search_state;
  if (!archive_search->emulate_rc_search) {
    // Compute the full search
    approximate_search(forward_asearch,matches); // Forward search
  } else {
    // Configure search stage to stop at
    const as_parameters_t* const actual_parameters = &archive_search->as_parameters;
    const bool lower_max_difference = actual_parameters->complete_strata_after_best_nominal < forward_asearch->max_differences;
    if (lower_max_difference && archive_search->probe_strand) forward_asearch->stop_before_neighborhood_search = true;
    // Run the search (FORWARD)
    approximate_search(forward_asearch,matches); // Forward search
    // Check the number of matches & keep searching
    if (!forward_asearch->max_matches_reached) {
      // Keep on searching
      approximate_search_t* const reverse_asearch = &archive_search->reverse_search_state;
      // Run the search (REVERSE)
      approximate_search(reverse_asearch,matches); // Reverse emulated-search
      // Resume forward search (if not completed before)
      if (forward_asearch->search_state != asearch_end && !forward_asearch->max_matches_reached) {
        approximate_search(forward_asearch,matches);
      }
    }
  }
  // Select matches
  archive_select_matches(archive_search,true,true,true,matches);
  // Check matches
  select_parameters_t* const select_parameters = archive_search->select_parameters;
  if (select_parameters->check_correct || select_parameters->check_optimum || select_parameters->check_complete) {
    search_parameters_t* const search_parameters = archive_search->as_parameters.search_parameters;
    archive_check_matches(
        archive_search->archive,search_parameters->alignment_model,
        &search_parameters->swg_penalties,&archive_search->sequence,
        matches,select_parameters->check_optimum,select_parameters->check_complete,
        archive_search->mm_stack);
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_SE);
}
/*
 * PairedEnd Indexed Search (PE Online Approximate String Search)
 */
GEM_INLINE void archive_search_paired_end_extend_matches(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,const sequence_end_t candidate_end) {
  PROF_START(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES);
  // Parameters
  search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
  mapper_stats_t* const mapper_stats = archive_search_end1->mapper_stats;
  archive_t* const archive = archive_search_end1->archive;
  // Extend in all possible concordant orientations
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  /*
   * Configure extension
   *   All extensions are done against the forward strand. If the candidate is in the reverse,
   *   the reverse-candidate region is generated reverse-complementing the forward. Thus, all
   *   the extension is done using @forward_search_state (matches contain all the strand needed info)
   */
  archive_search_t* candidate_archive_search;
  matches_t* candidate_matches;
  matches_t* extended_matches;
  as_parameters_t* as_parameters;
  if (candidate_end==paired_end2) {
    // Extend towards end/2
    candidate_archive_search = archive_search_end2;
    as_parameters = &candidate_archive_search->as_parameters;
    extended_matches = matches_end1;
    candidate_matches = matches_end2;
  } else {
    // Extend towards end/1
    candidate_archive_search = archive_search_end1;
    as_parameters = &candidate_archive_search->as_parameters;
    extended_matches = matches_end2;
    candidate_matches = matches_end1;
  }
  filtering_candidates_t* const filtering_candidates = candidate_archive_search->forward_search_state.filtering_candidates;
  pattern_t* const pattern = &candidate_archive_search->forward_search_state.pattern;
  text_collection_t* const text_collection = candidate_archive_search->text_collection;
  mm_stack_t* const mm_stack = candidate_archive_search->mm_stack;
  const uint64_t num_base_candidate_matches = matches_get_num_match_traces(candidate_matches);
  uint64_t matches_found = 0, total_matches_found = 0;
  // Iterate over all matches of the extended end
  VECTOR_ITERATE(extended_matches->position_matches,extended_match,en,match_trace_t) {
    if (search_parameters->pair_orientation_FR == pair_orientation_concordant) {
      // Extend (filter nearby region)
      strand_t candidate_strand;
      if (extended_match->strand==Forward) {
        candidate_strand = Reverse;
        matches_found = filtering_candidates_extend_match(
            filtering_candidates,archive,text_collection,extended_match,pattern,true,true,
            as_parameters,mapper_stats,paired_matches,candidate_end,mm_stack);
      } else { // end2==Reverse
        candidate_strand = Forward;
        extended_match->match_alignment.match_position =  // Locates against forward strand // TODO improve using Text-Dimensions
            inverse_locator_map(archive->locator,(uint8_t*)extended_match->sequence_name,
            extended_match->strand,extended_match->text_position);
        matches_found = filtering_candidates_extend_match(
            filtering_candidates,archive,text_collection,extended_match,pattern,false,false,
            as_parameters,mapper_stats,paired_matches,candidate_end,mm_stack);
      }
      // Decode & Pair found matches
      total_matches_found += matches_found;
	  match_trace_t* const mates_array = matches_get_match_traces(candidate_matches) + num_base_candidate_matches;
      archive_select_decode_trace_matches(candidate_archive_search,
          candidate_matches,mates_array,total_matches_found,true,candidate_strand);
      if (matches_found > 0) {
        // Pair with extended matches
        paired_matches_pair_match_with_mates(paired_matches,search_parameters,mapper_stats,
            pair_orientation_concordant,extended_match,candidate_end,mates_array,total_matches_found);
      }
    }
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES);
}
GEM_INLINE void archive_search_paired_end_generate_extension_candidates(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,const sequence_end_t candidate_end) {
  // Check
  archive_t* const archive = archive_search_end1->archive;
  gem_check(!archive->indexed_complement || archive_search_end1->emulate_rc_search,ARCHIVE_SEARCH_INDEX_COMPLEMENT_REQUIRED);
  // Generate extension candidates
  search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
  mapper_stats_t* const mapper_stats = archive_search_end1->mapper_stats;
  approximate_search_t* const forward_asearch_end1 = &archive_search_end1->forward_search_state;
  approximate_search_t* const forward_asearch_end2 = &archive_search_end2->forward_search_state;
  if (candidate_end==paired_end2) {
    filtering_candidates_process_extension_candidates(
        forward_asearch_end1->filtering_candidates,forward_asearch_end2->filtering_candidates,
        archive,archive_search_end1->text_collection,&forward_asearch_end1->pattern,&forward_asearch_end2->pattern,
        search_parameters,mapper_stats,paired_matches,archive_search_end1->mm_stack);
  } else {
    filtering_candidates_process_extension_candidates(
        forward_asearch_end2->filtering_candidates,forward_asearch_end1->filtering_candidates,
        archive,archive_search_end1->text_collection,&forward_asearch_end2->pattern,&forward_asearch_end1->pattern,
        search_parameters,mapper_stats,paired_matches,archive_search_end1->mm_stack);
  }
}
GEM_INLINE void archive_search_paired_end_discard_filtering_regions(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_PAIRED_FILTERED);
  PROF_START(GP_ARCHIVE_SEARCH_PE_DISCARD_FILTERING_REGIONS);
  // Filtering candidates
  archive_t* const archive = archive_search_end1->archive;
  mapper_stats_t* const mapper_stats = archive_search_end1->mapper_stats;
  approximate_search_t* const forward_search_end1 = &archive_search_end1->forward_search_state;
  approximate_search_t* const forward_search_end2 = &archive_search_end2->forward_search_state;
  approximate_search_t* reverse_search_end1, *reverse_search_end2;
  const bool indexed_complement = archive_search_end1->archive->indexed_complement;
  if (indexed_complement) {
    const bool verify_candidates_both_ends =
        (forward_search_end1->search_state == asearch_verify_candidates &&
         forward_search_end2->search_state == asearch_verify_candidates);
    if (verify_candidates_both_ends) {
      // Initialize (Invalidates all candidate-regions as lacking mate)
      filtering_candidates_set_all_regions_pending(forward_search_end1->filtering_candidates);
      filtering_candidates_set_all_regions_pending(forward_search_end2->filtering_candidates);
      // Search for compatible pairs (Validates candidate-regions)
      search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
      filtering_candidates_paired_regions_filtering(
          forward_search_end1->filtering_candidates,forward_search_end2->filtering_candidates,
          archive->text,search_parameters,mapper_stats,paired_matches);
    }
  } else {
    reverse_search_end1 = &archive_search_end1->reverse_search_state;
    reverse_search_end2 = &archive_search_end2->reverse_search_state;
    const bool verify_candidates_both_ends =
        (forward_search_end1->search_state == asearch_verify_candidates &&
         reverse_search_end1->search_state == asearch_verify_candidates &&
         forward_search_end2->search_state == asearch_verify_candidates &&
         reverse_search_end2->search_state == asearch_verify_candidates);
    if (verify_candidates_both_ends) {
      // Initialize (Invalidates all candidate-regions as lacking mate)
      filtering_candidates_set_all_regions_pending(forward_search_end1->filtering_candidates);
      filtering_candidates_set_all_regions_pending(forward_search_end2->filtering_candidates);
      filtering_candidates_set_all_regions_pending(reverse_search_end1->filtering_candidates);
      filtering_candidates_set_all_regions_pending(reverse_search_end2->filtering_candidates);
      // Search for compatible pairs (Validates candidate-regions)
      search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
      filtering_candidates_paired_regions_filtering(
          forward_search_end1->filtering_candidates,reverse_search_end2->filtering_candidates,
          archive->text,search_parameters,mapper_stats,paired_matches);
      filtering_candidates_paired_regions_filtering(
          reverse_search_end1->filtering_candidates,forward_search_end2->filtering_candidates,
          archive->text,search_parameters,mapper_stats,paired_matches);
    }
  }
  // Profile (Count the number of discarded filtering regions)
  PROF_BLOCK() {
    uint64_t num_regions = 0, num_regions_discarded = 0;
    num_regions_discarded += filtering_candidates_count_candidate_regions(
        forward_search_end1->filtering_candidates,filtering_region_pending);
    num_regions_discarded += filtering_candidates_count_candidate_regions(
        forward_search_end2->filtering_candidates,filtering_region_pending);
    if (!indexed_complement) {
      num_regions_discarded += filtering_candidates_count_candidate_regions(
          reverse_search_end1->filtering_candidates,filtering_region_pending);
      num_regions_discarded += filtering_candidates_count_candidate_regions(
          reverse_search_end2->filtering_candidates,filtering_region_pending);
    }
    num_regions += filtering_candidates_get_num_candidate_regions(forward_search_end1->filtering_candidates);
    num_regions += filtering_candidates_get_num_candidate_regions(forward_search_end2->filtering_candidates);
    if (!indexed_complement) {
      num_regions += filtering_candidates_get_num_candidate_regions(reverse_search_end1->filtering_candidates);
      num_regions += filtering_candidates_get_num_candidate_regions(reverse_search_end2->filtering_candidates);
    }
    PROF_ADD_COUNTER(GP_ARCHIVE_SEARCH_PE_DISCARD_FILTERING_REGIONS_NOT_CONCORDANT,num_regions_discarded);
    PROF_ADD_COUNTER(GP_ARCHIVE_SEARCH_PE_DISCARD_FILTERING_REGIONS_TOTAL,num_regions);
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_PE_DISCARD_FILTERING_REGIONS);
}
GEM_INLINE bool archive_search_paired_end_fulfilled(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  return vector_get_used(paired_matches->matches) > 0;
}
GEM_INLINE bool archive_search_paired_end_use_paired_filtering(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    const paired_matches_t* const paired_matches) {
  search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
  // Check mapping parameters
  if (search_parameters->paired_mapping_mode != paired_mapping_paired_filtering) {
    return false;
  } else {
    const uint64_t max_error =
        archive_search_end1->as_parameters.max_filtering_error_nominal +
        archive_search_end2->as_parameters.max_filtering_error_nominal;
    return mapper_stats_template_length_estimation_within_ci(archive_search_end1->mapper_stats,max_error);
  }
}
GEM_INLINE bool archive_search_paired_end_feasible_extension(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    const sequence_end_t candidate_end,const paired_matches_t* const paired_matches) {
  search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
  // Check mapping parameters
  if (search_parameters->paired_mapping_mode == paired_mapping_map_both_ends) return false;
  if (search_parameters->paired_mapping_mode == paired_mapping_paired_filtering) return false;
  // Check the number of samples to derive the expected insert size
  const uint64_t length_end1 = sequence_get_length(&archive_search_end1->sequence);
  const uint64_t length_end2 = sequence_get_length(&archive_search_end2->sequence);
  const uint64_t max_error = (uint64_t)(search_parameters->max_search_error*(double)(length_end1+length_end2));
  if (!mapper_stats_template_length_estimation_within_ci(archive_search_end1->mapper_stats,max_error)) return false;
  // Check the number of candidates
  if (candidate_end==paired_end2) { // Extend end/1 => TODO dynamic choice (learning & balancing)
    if (length_end2 >= 1000) return false;
    if (archive_search_get_search_exact_matches(archive_search_end1) > 0) return false;
    return true;
//    const uint64_t candidates_end1 = archive_search_get_search_canditates(archive_search_end1);
//    return candidates_end1>0 && candidates_end1<=search_parameters->max_extendable_candidates;
  } else { // Extend end/2
    if (length_end1 >= 1000) return false;
    if (archive_search_get_search_exact_matches(archive_search_end2) > 0) return false;
    return true;
//    const uint64_t candidates_end2 = archive_search_get_search_canditates(archive_search_end2);
//    return candidates_end2>0 && candidates_end2<=search_parameters->max_extendable_candidates;
  }
}
GEM_INLINE bool archive_search_paired_end_use_extension(matches_t* const matches) {
  return matches_classify_unique(matches,matches->max_complete_stratum) > 0.98;
}
GEM_INLINE bool archive_search_paired_end_use_recovery_by_extension(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    const sequence_end_t candidate_end,const paired_matches_t* const paired_matches) {
  search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
  // Check the number of samples to derive the expected insert size
  const uint64_t max_error =
      archive_search_end1->as_parameters.max_filtering_error_nominal +
      archive_search_end2->as_parameters.max_filtering_error_nominal;
  if (!mapper_stats_template_length_estimation_within_ci(archive_search_end1->mapper_stats,max_error)) return false;
  // Check the number of candidates
  if (candidate_end==paired_end2) { // Extend end/1
    const uint64_t matches_end1 = matches_get_num_match_traces(paired_matches->matches_end1);
    return matches_end1>0 && matches_end1<=search_parameters->max_extendable_candidates;
  } else { // extended_end==ARCHIVE_SEARCH_END2
    const uint64_t matches_end2 = matches_get_num_match_traces(paired_matches->matches_end2);
    return matches_end2>0 && matches_end2<=search_parameters->max_extendable_candidates;
  }
}
/*
 * PE Online Approximate String Search
 */
GEM_INLINE void archive_search_paired_end_continue(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  PROF_START(GP_ARCHIVE_SEARCH_PE);
  // Parameters
  search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
  mapper_stats_t* const mapper_stats = archive_search_end1->mapper_stats;
  matches_t* const matches_end1 = (paired_matches!=NULL) ? paired_matches->matches_end1 : NULL;
  matches_t* const matches_end2 = (paired_matches!=NULL) ? paired_matches->matches_end2 : NULL;
  // Callback (switch to proper search stage)
  archive_search_paired_end_callback:
  switch (archive_search_end1->pe_search_state) {
    case archive_search_pe_begin: // Beginning of the search (Init)
      archive_search_end1->paired_extending = false;
      archive_search_end1->paired_filtering = false;
      archive_search_end2->paired_extending = false;
      archive_search_end2->paired_filtering = false;
    // No break
    case archive_search_pe_search_end1:
      // Full Search (End/1)
      archive_search_reset(archive_search_end1);
      archive_search_finish_search(archive_search_end1,matches_end1);
      // Test if extension is feasible (End/1)
      const bool feasible_extension_end1 = archive_search_paired_end_feasible_extension(
          archive_search_end1,archive_search_end2,paired_end2,paired_matches);
      if (feasible_extension_end1) {
        archive_search_end1->paired_extending = archive_search_paired_end_use_extension(matches_end1);
        if (archive_search_end1->paired_extending) {
          // Go for extending (End/1)
          PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END1);
          archive_search_reset(archive_search_end2); // Init (End/2)
          archive_search_end1->pe_search_state = archive_search_pe_extend_end1;
          goto archive_search_paired_end_callback; // Callback
        }
      }
    // No break
    case archive_search_pe_search_end2:
      // Full Search (End/2)
      archive_search_reset(archive_search_end2);
      archive_search_finish_search(archive_search_end2,matches_end2);
      // Callback
      archive_search_end1->pe_search_state = archive_search_pe_paired_filtering_verified;
      goto archive_search_paired_end_callback;
    // No break
//    case archive_search_pe_paired_filtering_discard: // Discard paired-filtering candidates
//      archive_search_end1->paired_filtering =
//          archive_search_paired_end_use_paired_filtering(archive_search_end1,archive_search_end2,paired_matches);
//      if (archive_search_end1->paired_filtering) {
//        // Discard non feasible matches using paired-filtering
//        archive_search_paired_end_discard_filtering_regions(archive_search_end1,archive_search_end2,paired_matches);
//      } else {
//        // Callback
//        archive_search_end1->pe_search_state = archive_search_pe_both_ends_verify;
//        goto archive_search_paired_end_callback;
//      }
//    // No break
//    case archive_search_pe_paired_filtering_verify: // Verify paired-filtering candidates
//      // Paired-filtering verify
//      archive_search_verify_candidates(archive_search_end1,paired_matches->matches_end1);
//      archive_search_verify_candidates(archive_search_end2,paired_matches->matches_end2);
//      archive_select_decode_trace_matches_all(archive_search_end1,paired_matches->matches_end1,false,Forward);
//      archive_select_decode_trace_matches_all(archive_search_end2,paired_matches->matches_end2,false,Forward);
//      // Callback
//      archive_search_end1->pe_search_state = archive_search_pe_paired_filtering_verified;
//      goto archive_search_paired_end_callback;
//      break;
    case archive_search_pe_extend_end1: // Extend (End/1)
      PROF_START(GP_ARCHIVE_SEARCH_PE_EXTENSION);
      archive_search_paired_end_extend_matches(archive_search_end1,archive_search_end2,paired_matches,paired_end2);
      PROF_STOP(GP_ARCHIVE_SEARCH_PE_EXTENSION);
    // No break
    case archive_search_pe_extended_end1: // End1 extended
      if (archive_search_paired_end_fulfilled(archive_search_end1,archive_search_end2,paired_matches)) {
        PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END1_SUCCESS);
        archive_search_end1->pe_search_state = archive_search_pe_end;
        goto archive_search_paired_end_callback; // Callback
      }
      /*
       * Extension failed because
       *  (1) Insert size is beyond the expected distribution (Any insert size filtering/search must be discarded)
       *  (2) Sensitivity of the candidates/end1 search is not enough
       */
      // Archive search (End/2)
      archive_search_finish_search(archive_search_end2,matches_end2); // Full Search (End/2)
      archive_search_end1->pe_search_state = archive_search_pe_both_ends_verified;
      goto archive_search_paired_end_callback; // Callback
      break;
//    case archive_search_pe_extend_end2: // Extend (End/2)
//      PROF_START(GP_ARCHIVE_SEARCH_PE_EXTENSION);
//      archive_search_paired_end_extend_matches(archive_search_end1,archive_search_end2,paired_matches,paired_end1);
//      PROF_STOP(GP_ARCHIVE_SEARCH_PE_EXTENSION);
//    // No break
//    case archive_search_pe_extended_end2: // End2 extended
//      if (archive_search_paired_end_fulfilled(archive_search_end1,archive_search_end2,paired_matches)) {
//        PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END2_SUCCESS);
//        archive_search_end1->pe_search_state = archive_search_pe_end;
//        goto archive_search_paired_end_callback; // Callback
//      }
//      /*
//       * Extension failed because
//       *  (1) Insert size is beyond the expected distribution (Any insert size filtering/search must be discarded)
//       *  (2) Sensitivity of the candidates/end2 search is not enough
//       */
//      // Finish search (End/1)
//      archive_search_release_verification_candidates(archive_search_end1);
//      if (!archive_search_finished(archive_search_end1)) {
//        archive_search_finish_search(archive_search_end1,matches_end1); // Finish Search (End/1)
//      }
//      archive_search_end1->pe_search_state = archive_search_pe_both_ends_verified;
//      goto archive_search_paired_end_callback; // Callback
//      break;
    case archive_search_pe_paired_filtering_verified: { // Paired filtering applied
      // Pair matches (Cross-link matches from both ends)
      const uint64_t num_matches_end1 = matches_get_num_match_traces(paired_matches->matches_end1);
      const uint64_t num_matches_end2 = matches_get_num_match_traces(paired_matches->matches_end2);
      if (num_matches_end1 > 0 && num_matches_end2 > 0) {
        paired_matches_find_pairs(paired_matches,search_parameters,mapper_stats);
      }
      // Check number of paired-matches
      if (archive_search_paired_end_fulfilled(archive_search_end1,archive_search_end2,paired_matches)) {
        PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_PAIRED_FILTERING_SUCCESS);
        // Callback
        archive_search_end1->pe_search_state = archive_search_pe_end;
        goto archive_search_paired_end_callback;
      } else {
        // Clear previous paired_matches (restart)
        paired_matches_clear(paired_matches);
        // Release candidates
        archive_search_release_verification_candidates(archive_search_end1);
        archive_search_release_verification_candidates(archive_search_end2);
      }
    }
    // No break
    case archive_search_pe_both_ends_verify: // Verify candidates for both ends
      // Finish search both ends
      if (!archive_search_finished(archive_search_end1)) {
        archive_search_finish_search(archive_search_end1,matches_end1); // Finish search (End/1)
      }
      if (!archive_search_finished(archive_search_end2)) {
        archive_search_finish_search(archive_search_end2,matches_end2); // Finish search (End/2)
      }
    // No break
    case archive_search_pe_both_ends_verified: { // Candidates of both ends had been verified
      // Pair matches (Cross-link matches from both ends)
      const uint64_t num_matches_end1 = matches_get_num_match_traces(paired_matches->matches_end1);
      const uint64_t num_matches_end2 = matches_get_num_match_traces(paired_matches->matches_end2);
      if (num_matches_end1 > 0 && num_matches_end2 > 0) {
        paired_matches_find_pairs(paired_matches,search_parameters,mapper_stats);
        paired_matches_find_discordant_pairs(paired_matches,search_parameters); // Find discordant (if required)
      }
      const bool search_fulfilled = archive_search_paired_end_fulfilled(archive_search_end1,archive_search_end2,paired_matches);
      if (search_fulfilled || !search_parameters->recovery_by_extension) {
        archive_search_end1->pe_search_state = archive_search_pe_end;
        goto archive_search_paired_end_callback;
      }
    }
    // No break
    case archive_search_pe_recovery: // Paired-end recovery by extension
      // Recover Matches (by extension if not done before)
      if (!archive_search_end1->paired_extending) {
        if (archive_search_paired_end_use_recovery_by_extension(archive_search_end1,archive_search_end2,paired_end2,paired_matches)) {
          PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_RECOVER_BY_EXTENSION_END1);
          PROF_START(GP_ARCHIVE_SEARCH_PE_RECOVER_BY_EXTENSION);
          archive_search_paired_end_extend_matches(archive_search_end1,archive_search_end2,paired_matches,paired_end2);
          PROF_STOP(GP_ARCHIVE_SEARCH_PE_RECOVER_BY_EXTENSION);
        }
      }
      if (!archive_search_end2->paired_extending) {
        if (archive_search_paired_end_use_recovery_by_extension(archive_search_end1,archive_search_end2,paired_end1,paired_matches)) {
          PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_RECOVER_BY_EXTENSION_END2);
          PROF_START(GP_ARCHIVE_SEARCH_PE_RECOVER_BY_EXTENSION);
          archive_search_paired_end_extend_matches(archive_search_end1,archive_search_end2,paired_matches,paired_end1);
          PROF_STOP(GP_ARCHIVE_SEARCH_PE_RECOVER_BY_EXTENSION);
        }
      }
      PROF_BLOCK() { if (vector_get_used(paired_matches->matches) > 0) {
        PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_RECOVER_BY_EXTENSION_HIT);
      }}
    // No break
    case archive_search_pe_end: // End of the current workflow
      // Select matches
      if (vector_get_used(paired_matches->matches) > 0) {
        paired_matches->max_complete_stratum =
            ((matches_end1->max_complete_stratum!=ALL) ? matches_end1->max_complete_stratum : 0) +
            ((matches_end2->max_complete_stratum!=ALL) ? matches_end2->max_complete_stratum : 0);
        archive_select_paired_matches(archive_search_end1,archive_search_end2,paired_matches);
        // Check matches
        select_parameters_t* const select_parameters = archive_search_end1->select_parameters;
        if (select_parameters->check_correct || select_parameters->check_optimum || select_parameters->check_complete) {
          search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
          archive_check_paired_matches(
              archive_search_end1->archive,search_parameters->alignment_model,
              &search_parameters->swg_penalties,&archive_search_end1->sequence,
              &archive_search_end2->sequence,paired_matches,select_parameters->check_optimum,
              select_parameters->check_complete,archive_search_end1->mm_stack);
        }
      }
      archive_search_end1->pe_search_state = archive_search_pe_end;
      break; // End of the workflow
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_PE);
}
GEM_INLINE void archive_search_pe_generate_candidates(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  PROF_START(GP_ARCHIVE_SEARCH_PE);
  PROF_START(GP_ARCHIVE_SEARCH_PE_GENERATE_CANDIDATES);
  // Beginning of the search (Init)
  archive_search_end1->pe_search_state = archive_search_pe_begin;
  archive_search_end1->paired_extending = false;
  archive_search_end1->paired_filtering = false;
  archive_search_end2->pe_search_state = archive_search_pe_begin;
  archive_search_end2->paired_extending = false;
  archive_search_end2->paired_filtering = false;
  // Select PE-mapping-mode
  search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
  switch (search_parameters->paired_mapping_mode) {
    case paired_mapping_map_both_ends:
    case paired_mapping_paired_filtering:
      // Generate candidates (End/1 & End/2)
      archive_search_generate_candidates(archive_search_end1);
      archive_search_generate_candidates(archive_search_end2);
      // Use Paired-filtering (if applies)
      if (archive_search_paired_end_use_paired_filtering(archive_search_end1,archive_search_end2,paired_matches)) {
        archive_search_paired_end_discard_filtering_regions(archive_search_end1,archive_search_end2,paired_matches);
        archive_search_end1->pe_search_state = archive_search_pe_paired_filtering_verify;
      } else {
        archive_search_end1->pe_search_state = archive_search_pe_both_ends_verify;
      }
      break;
    case paired_mapping_map_extension:
      archive_search_generate_candidates(archive_search_end1);
      archive_search_end1->paired_extending = archive_search_paired_end_feasible_extension(
          archive_search_end1,archive_search_end2,paired_end2,paired_matches);
      if (archive_search_end1->paired_extending) {
        PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END1);
        archive_search_reset(archive_search_end2); // Init end/2
        archive_search_paired_end_generate_extension_candidates(
            archive_search_end1,archive_search_end2,paired_matches,paired_end2); // Generate candidates (from extension)
        archive_search_end1->pe_search_state = archive_search_pe_extended_end1;
      } else {
        archive_search_generate_candidates(archive_search_end2);
        // Use Paired-filtering (if applies)
        archive_search_end1->paired_filtering =
            archive_search_paired_end_use_paired_filtering(archive_search_end1,archive_search_end2,paired_matches);
        if (archive_search_end1->paired_filtering) {
          archive_search_paired_end_discard_filtering_regions(archive_search_end1,archive_search_end2,paired_matches);
          archive_search_end1->pe_search_state = archive_search_pe_paired_filtering_verify;
        } else {
          archive_search_end1->pe_search_state = archive_search_pe_both_ends_verify;
        }
      }
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROF_STOP(GP_ARCHIVE_SEARCH_PE_GENERATE_CANDIDATES);
  PROF_STOP(GP_ARCHIVE_SEARCH_PE);
}
GEM_INLINE void archive_search_pe_finish_search(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  PROF_START(GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH);
  // Parameters
  search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
  mapper_stats_t* const mapper_stats = archive_search_end1->mapper_stats;
  // Handle special cases & restart regular workflow
  if (archive_search_end1->pe_search_state==archive_search_pe_extended_end1) {
    // Pair matches (Cross-link matches from both ends)
    const uint64_t num_matches_end1 = matches_get_num_match_traces(paired_matches->matches_end1);
    const uint64_t num_matches_end2 = matches_get_num_match_traces(paired_matches->matches_end2);
    if (num_matches_end1 > 0 && num_matches_end2 > 0) {
      archive_select_decode_trace_matches_all(archive_search_end1,paired_matches->matches_end1,false,Forward);
      archive_select_decode_trace_matches_all(archive_search_end2,paired_matches->matches_end2,false,Forward);
      paired_matches_find_pairs(paired_matches,search_parameters,mapper_stats);
    }
    // Check number of matches
    if (archive_search_paired_end_fulfilled(archive_search_end1,archive_search_end2,paired_matches)) {
      PROF_INC_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END1_SUCCESS);
      archive_search_end1->pe_search_state = archive_search_pe_end; // Callback
    } else {
      // Clear previous paired_matches (restart)
      paired_matches_clear(paired_matches);
      // Map-both-ends
      archive_search_finish_search(archive_search_end1,paired_matches->matches_end1); // Finish search (End/1)
      archive_search_finish_search(archive_search_end2,paired_matches->matches_end2); // Full search (End/2)
      archive_search_end1->pe_search_state = archive_search_pe_both_ends_verified;    // Callback
    }
  } else if (archive_search_end1->pe_search_state==archive_search_pe_paired_filtering_verify) {
    archive_select_decode_trace_matches_all(archive_search_end1,paired_matches->matches_end1,false,Forward);
    archive_select_decode_trace_matches_all(archive_search_end2,paired_matches->matches_end2,false,Forward);
    // Callback
    archive_search_end1->pe_search_state = archive_search_pe_paired_filtering_verified;
  }
  // PE search (continue search until the end)
  archive_search_paired_end_continue(archive_search_end1,archive_search_end2,paired_matches);
  PROF_STOP(GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH);
}
GEM_INLINE void archive_search_paired_end(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Init
  archive_search_end1->pe_search_state = archive_search_pe_begin;
  archive_search_end2->pe_search_state = archive_search_pe_begin;
  // PE search
  archive_search_paired_end_continue(archive_search_end1,archive_search_end2,paired_matches);
}

