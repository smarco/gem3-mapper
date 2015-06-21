/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_filtering.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search_filtering.h"
#include "approximate_search_filtering_base.h"
#include "approximate_search_filtering_control.h"
#include "approximate_search_neighborhood.h"

#include "mapper_profile.h"
#include "region_profile_schedule.h"

/*
 * Read Recovery
 */
GEM_INLINE void approximate_search_read_recovery(
    approximate_search_t* const search,const bool verify_candidates,matches_t* const matches) {
  PROF_START(GP_AS_READ_RECOVERY);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,">ASM::Read Recovery\n");
  }
  /*
   * If the number of wildcards (or errors required) is greater than the maximum number of
   * differences allowed, we try to recover as many matches as possible.
   * We extract feasible regions from the read and filter them trying to recover anything
   * out of bad quality reads.
   */
  // Parameters
  as_parameters_t* const actual_parameters = search->as_parameters;
  search_parameters_t* const parameters = actual_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  // Select proper key
  const archive_t* archive = search->archive;
  const uint8_t* key;
  uint64_t key_length;
  if (!archive->text->run_length) {
    key = pattern->key;
    key_length = pattern->key_length;
  } else {
    key = pattern->rl_key;
    key_length = pattern->rl_key_length;
  }
  // Compute the region profile
  region_profile_generate_adaptive(region_profile,search->archive->fm_index,key,key_length,
      parameters->allowed_enc,&parameters->rp_recovery,search->max_differences+1,false);
  // Append the results to filter
  region_profile_schedule_filtering_fixed(region_profile,ALL,REGION_FILTER_DEGREE_ZERO,parameters->filtering_threshold);
  approximate_search_generate_exact_candidates(search,matches);
  approximate_search_update_mcs(search,region_profile->errors_allowed + pattern->num_wildcards);
  // Process candidates
  filtering_candidates_process_candidates(filtering_candidates,
      search->archive,pattern,actual_parameters,verify_candidates,search->mm_stack);
  // Verify/Process candidates
  if (verify_candidates) {
    approximate_search_verify_candidates(search,matches);
  } else {
    // Process candidates (just prepare)
    filtering_candidates_process_candidates(filtering_candidates,search->archive,pattern,actual_parameters,false,search->mm_stack);
    search->search_state = asearch_verify_candidates;
  }
  PROF_STOP(GP_AS_READ_RECOVERY);
}
GEM_INLINE void approximate_search_exact_filtering_adaptive(
    approximate_search_t* const search,const bool verify_candidates,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_EXACT);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,">ASM::Adaptive Filtering (Exact)\n");
  }
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  // Region-Minimal Profile (Reduce the number of candidates per region and maximize number of regions)
  approximate_search_generate_region_profile(search,
      &parameters->rp_minimal,region_profile_adaptive_lightweight,0,false);
  // Check corner cases
  if (search->search_state==asearch_no_regions || search->search_state==asearch_exact_matches) {
    asearch_control_next_state(search,asearch_exact_filtering_adaptive,matches);
  } else {
    // Region-Minimal Profile Stats
    approximate_search_generate_region_profile_minimal_stats(search);
    // Generate exact-candidates
    const bool verify_ahead = verify_candidates && (parameters->mapping_mode==mapping_adaptive_filtering_fast);
    region_profile_schedule_filtering_fixed(region_profile,
        search->max_differences+1,REGION_FILTER_DEGREE_ZERO,parameters->filtering_threshold);
    if (!verify_ahead) {
      approximate_search_generate_exact_candidates(search,matches);
    } else {
      approximate_search_generate_inexact_candidates(search,false,verify_ahead,matches);
    }
    // Update MCS (maximum complete stratum) [Hint]
    approximate_search_update_mcs(search,region_profile->errors_allowed + pattern->num_wildcards);
    // Process candidates (just prepare to verification)
    filtering_candidates_process_candidates(search->filtering_candidates,
        search->archive,pattern,actual_parameters,verify_candidates,search->mm_stack);
    // Verify Candidates (if needed)
    if (verify_candidates) {
      approximate_search_verify_candidates(search,matches);
    } else {
      search->search_state = asearch_verify_candidates;
    }
  }
  PROF_STOP(GP_AS_FILTERING_EXACT);
}
GEM_INLINE void approximate_search_exact_filtering_boost(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_EXACT_BOOST);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,">ASM::Adaptive Filtering (Boost)\n");
  }
//  // Parameters
//  const as_parameters_t* const actual_parameters = search->as_parameters;
//  const search_parameters_t* const parameters = actual_parameters->search_parameters;
//  pattern_t* const pattern = &search->pattern;
//  // Region-Boost Profile
////  fm_index_t* const fm_index = search->archive->fm_index;
////  const uint64_t base_num_regions = search->region_profile.num_filtering_regions;
////  const uint64_t potential_num_regions = pattern->key_length/(fm_index->proper_length);
////  if (base_num_regions+1 >= potential_num_regions) {
////    asearch_control_next_state(search,asearch_exact_filtering_boost,matches);
////    PROF_STOP(GP_AS_FILTERING_EXACT_BOOST);
////  }
//
//
//  approximate_search_generate_region_profile(search,
//      &parameters->rp_boost,region_profile_adaptive_extensive,0,false);
//  if (//base_num_regions <= search->region_profile.num_filtering_regions ||
//      search->search_state == asearch_no_regions || search->search_state == asearch_exact_matches) {
//    asearch_control_next_state(search,asearch_exact_filtering_boost,matches);
//    PROF_STOP(GP_AS_FILTERING_EXACT_BOOST);
//    return;
//  }
//
//  approximate_search_generate_region_profile_boost_stats(search); // Stats
//  // Generate exact-candidates (Dynamic filtering incorporated)
//  approximate_search_generate_candidates(search,matches,region_filter_adaptive_exact,0,true,true);
//  // Process candidates (just prepare to verification)
//  filtering_candidates_process_candidates(search->filtering_candidates,
//      search->archive,pattern,actual_parameters,true,search->mm_stack);
//  // Verify candidates
//  approximate_search_verify_candidates(search,matches);
//  search->search_state = asearch_exact_filtering_boost; // Correct State

  // Next State
  asearch_control_next_state(search,asearch_exact_filtering_boost,matches);
  PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_BOOST_MAPPED,matches_is_mapped(matches)?1:0);
  PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_BOOST_MCS,search->max_complete_stratum);
  PROF_STOP(GP_AS_FILTERING_EXACT_BOOST);
}
GEM_INLINE void approximate_search_inexact_filtering(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_INEXACT);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,">ASM::Adaptive Filtering (Inexact)\n");
  }
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  // Region-Delimit Profile (Maximize number of unique-regions and try to isolate standard/repetitive regions)
  approximate_search_generate_region_profile(search,
      &parameters->rp_delimit,region_profile_adaptive_lightweight,0,true);
  if (search->region_profile.num_filtering_regions <= 1) {
    asearch_control_next_state(search,asearch_inexact_filtering,matches);
    PROF_STOP(GP_AS_FILTERING_INEXACT);
    return;
  }
  if (search->search_state == asearch_no_regions || search->search_state == asearch_exact_matches) {
    PROF_STOP(GP_AS_FILTERING_INEXACT);
    return;
  }
  approximate_search_generate_region_profile_delimit_stats(search); // Stats
  // Generate exact-candidates (Dynamic filtering incorporated)
  const uint64_t proper_length = fm_index_get_proper_length(search->archive->fm_index);
  const uint64_t sensibility_misms_length = parameters->filtering_region_factor*proper_length;
  region_profile_schedule_filtering_adaptive(region_profile,search->max_differences,sensibility_misms_length);
  approximate_search_generate_inexact_candidates(search,true,true,matches);
  // Update MCS (maximum complete stratum) [Hint]
  approximate_search_update_mcs(search,region_profile->errors_allowed + pattern->num_wildcards);
  // Process candidates (just prepare to verification)
  filtering_candidates_process_candidates(search->filtering_candidates,
      search->archive,pattern,actual_parameters,true,search->mm_stack);
  // Verify candidates
  approximate_search_verify_candidates(search,matches);
  search->search_state = asearch_inexact_filtering; // Correct State
  // Next State
  asearch_control_next_state(search,asearch_inexact_filtering,matches);
  PROF_ADD_COUNTER(GP_AS_FILTERING_INEXACT_MAPPED,matches_is_mapped(matches)?1:0);
  PROF_ADD_COUNTER(GP_AS_FILTERING_INEXACT_MCS,search->max_complete_stratum);
  PROF_STOP(GP_AS_FILTERING_INEXACT);
}
GEM_INLINE void approximate_search_local_alignments_filtering(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_LOCAL_ALIGNMENTS);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,">ASM::Local-Alignment Filtering\n");
  }
  // Local-Align discarded candidates
  filtering_candidates_local_alignment(
      search->filtering_candidates,search->archive,search->text_collection,
      &search->pattern,search->emulated_rc_search,search->as_parameters,
      matches,search->mm_stack);
  // Next State
  asearch_control_next_state(search,asearch_local_alignment,matches);
  PROF_ADD_COUNTER(GP_AS_FILTERING_LOCAL_ALIGNMENTS_MAPPED,matches_is_mapped(matches)?1:0);
  PROF_ADD_COUNTER(GP_AS_FILTERING_LOCAL_ALIGNMENTS_MCS,search->max_complete_stratum);
  PROF_STOP(GP_AS_FILTERING_LOCAL_ALIGNMENTS);
}
/*
 * Basic Cases
 */
GEM_INLINE void approximate_search_adaptive_mapping_basic_cases(
    approximate_search_t* const search,const bool verify_candidates,matches_t* const matches) {
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,">ASM::Basic Cases\n");
  }
  // Parameters
  pattern_t* const pattern = &search->pattern;
  const uint64_t key_length = pattern->key_length;
  const uint64_t num_wildcards = pattern->num_wildcards;
  // Check if all characters are wildcards
  if (key_length==num_wildcards) {
    search->hi_exact_matches = 0;
    search->lo_exact_matches = 0;
    approximate_search_update_mcs(search,key_length);
    search->search_state = asearch_end;
    return;
  }
  // Check if recovery is needed
  if (num_wildcards > 0 && num_wildcards >= search->max_differences) {
    search->search_state = asearch_read_recovery;
    return;
  }
  // Exact search
  if (search->max_differences==0) {
    search->search_state = asearch_neighborhood;
    return;
  }
  // Very short reads (Neighborhood search)
  if (key_length <= RANK_MTABLE_SEARCH_DEPTH || key_length < search->archive->fm_index->proper_length) {
    search->search_state = asearch_neighborhood;
    return;
  }
  // Otherwise, go to standard exact filtering
  search->search_state = asearch_exact_filtering_adaptive;
  return;
}
/*
 * Approximate Search based on Adaptive filtering
 */
GEM_INLINE void approximate_search_filtering_adaptive(approximate_search_t* const search,matches_t* const matches) {
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  pattern_t* const pattern = &search->pattern;
  const bool verify_candidates = search->verify_candidates;
  // Process proper search-stage
  while (search->search_state != asearch_end) {
    switch (search->search_state) {
      case asearch_begin:
        // Search Start. Check basic cases
        approximate_search_adaptive_mapping_basic_cases(search,verify_candidates,matches);
        if (search->search_state==asearch_verify_candidates && !verify_candidates) return; // Return if no-filtering
        break;
      case asearch_read_recovery:
        // Read recovery
        approximate_search_read_recovery(search,verify_candidates,matches);
        if (!verify_candidates) return; // Return if no-filtering
        break;
      case asearch_exact_filtering_adaptive:
        // Exact-Filtering (Adaptive)
        approximate_search_exact_filtering_adaptive(search,verify_candidates,matches);
        if (!verify_candidates) return; // Return if no-filtering
        break;
      case asearch_verify_candidates:
        // Verify Candidates
        approximate_search_verify_candidates(search,matches);
        break;
      case asearch_candidates_verified:
        // Next State
        PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MAPPED,matches_is_mapped(matches)?1:0);
        PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_MCS,search->max_complete_stratum);
        asearch_control_next_state(search,asearch_exact_filtering_adaptive,matches);
        break;
      case asearch_exact_filtering_boost:
        // Exact-Filtering (Fixed)
        approximate_search_exact_filtering_boost(search,matches);
        break;
      case asearch_inexact_filtering:
        // Inexact-Filtering
        approximate_search_inexact_filtering(search,matches);
        break;
      case asearch_neighborhood:
        approximate_search_neighborhood_search(search,matches);
        break;
      case asearch_exact_matches: // Exact Matches
        // Add interval
        matches_add_interval_match(matches,search->lo_exact_matches,
            search->hi_exact_matches,pattern->key_length,0,search->emulated_rc_search);
        search->max_differences = actual_parameters->complete_strata_after_best_nominal; // Adjust max-differences
        approximate_search_update_mcs(search,1);
        search->search_state = asearch_end;
        break;
      case asearch_local_alignment: // Local alignments
        approximate_search_local_alignments_filtering(search,matches);
        break;
      case asearch_no_regions: // No regions found
        approximate_search_update_mcs(search,pattern->num_wildcards);
        search->search_state = asearch_end;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  // Done!!
  PROF_ADD_COUNTER(GP_AS_ADAPTIVE_MCS,search->max_complete_stratum);
}
/*
 * Approximate Search based on Filtering Complete Search
 */
GEM_INLINE void approximate_search_filtering_complete(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_EXACT);
  if (search->max_differences==0) {
    approximate_search_neighborhood_exact_search(search,matches); // Exact Search
  } else {
    // Parameters
    const as_parameters_t* const actual_parameters = search->as_parameters;
    search_parameters_t* const parameters = actual_parameters->search_parameters;
    fm_index_t* const fm_index = search->archive->fm_index;
    pattern_t* const pattern = &search->pattern;
    region_profile_t* const region_profile = &search->region_profile;
    // Select proper key
    const archive_t* archive = search->archive;
    const uint8_t* key;
    uint64_t key_length;
    if (!archive->text->run_length) {
      key = pattern->key;
      key_length = pattern->key_length;
    } else {
      key = pattern->rl_key;
      key_length = pattern->rl_key_length;
    }
    // Compute the region profile
    region_profile_generate_adaptive_limited(region_profile,fm_index,key,key_length,
        parameters->allowed_enc,&parameters->rp_minimal,search->max_differences+1);
    if (region_profile_has_exact_matches(region_profile)) {
      const region_search_t* const first_region = region_profile->filtering_region;
      matches_add_interval_match(matches,first_region->lo,
          first_region->hi,pattern->key_length,0,search->emulated_rc_search);
      approximate_search_update_mcs(search,1); // Update MCS
    } else {
      approximate_search_generate_region_profile_minimal_stats(search); // Stats
      // Generate exact-candidates
      parameters->filtering_threshold = ALL; // No restriction
      region_profile_schedule_filtering_fixed(region_profile,ALL,REGION_FILTER_DEGREE_ZERO,parameters->filtering_threshold);
      approximate_search_generate_inexact_candidates(search,false,false,matches);
      // Process candidates (just prepare to verification)
      filtering_candidates_process_candidates(search->filtering_candidates,
          search->archive,pattern,actual_parameters,true,search->mm_stack);
      // Verify candidates
      filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
      filtering_candidates_verify_candidates(
          filtering_candidates,search->archive,search->text_collection,
          pattern,actual_parameters,matches,search->mm_stack);
      // Align candidates
      filtering_candidates_align_candidates(filtering_candidates,search->archive->text,search->text_collection,
          pattern,search->emulated_rc_search,actual_parameters,false,matches,search->mm_stack);
      // Update MCS (maximum complete stratum) [Hint]
      approximate_search_update_mcs(search,region_profile->errors_allowed + pattern->num_wildcards);
    }
  }
  PROF_STOP(GP_AS_FILTERING_EXACT);
}

