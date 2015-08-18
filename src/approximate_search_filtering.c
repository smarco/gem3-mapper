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
 * Basic Cases
 */
GEM_INLINE void approximate_search_adaptive_mapping_basic_cases(
    approximate_search_t* const search,const bool verify_candidates,matches_t* const matches) {
  gem_cond_debug_block(DEBUG_SEARCH_STATE) {
    tab_fprintf(stderr,">ASM::Basic Cases\n");
  }
  // Parameters
  search_parameters_t* const search_parameters = search->as_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  const uint64_t key_length = pattern->key_length;
  const uint64_t num_wildcards = pattern->num_wildcards;
  // All characters are wildcards
  if (key_length==num_wildcards) {
    search->hi_exact_matches = 0;
    search->lo_exact_matches = 0;
    approximate_search_update_mcs(search,key_length);
    search->search_state = asearch_end;
    return;
  }
  /*
   * Recovery
   *   If the number of wildcards (or errors required) is too high we try to recover as
   *   many matches as possible. We extract feasible regions from the read and filter
   *   them trying to recover anything out of bad quality reads.
   */
  if (num_wildcards > 0) {
    switch (search_parameters->mapping_mode) {
      case mapping_adaptive_filtering_fast:
      case mapping_adaptive_filtering_thorough:
        if (num_wildcards >= search->as_parameters->alignment_max_error_nominal) {
          search->search_state = asearch_read_recovery;
          return;
        }
        break;
      default:
        if (num_wildcards >= search->max_complete_error) {
          search->search_state = asearch_read_recovery;
          return;
        }
        break;
    }
  }
  // Exact search
  if (search->max_complete_error==0) {
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
 * Exact Filtering Adaptive
 */
GEM_INLINE void approximate_search_exact_filtering_adaptive(
    approximate_search_t* const search,const region_profiling_strategy_t profiling_strategy,
    const bool verify_candidates,matches_t* const matches) {
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_fprintf(stderr,">ASM::Adaptive Filtering (Exact)\n"); }
  // Region-Minimal Profile (Reduce the number of candidates per region and maximize number of regions)
  approximate_search_generate_region_profile(search,profiling_strategy,search->mm_stack);
  // Check corner cases
  if (search->search_state==asearch_no_regions || search->search_state==asearch_exact_matches) {
    asearch_control_next_state(search,asearch_exact_filtering_adaptive,matches);
  } else {
    // Parameters
    const as_parameters_t* const actual_parameters = search->as_parameters;
    const search_parameters_t* const parameters = actual_parameters->search_parameters;
    // Generate exact-candidates
    region_profile_schedule_filtering_fixed(&search->region_profile,ALL,
        REGION_FILTER_DEGREE_ZERO,parameters->filtering_threshold);
    const bool verify_ahead = verify_candidates && (parameters->mapping_mode==mapping_adaptive_filtering_fast);
    if (!verify_ahead) {
      approximate_search_generate_exact_candidates(search,matches);
    } else {
      approximate_search_generate_inexact_candidates(search,false,verify_ahead,matches);
    }
    // Process candidates (just prepare to verification)
    filtering_candidates_process_candidates(search->filtering_candidates,
        search->archive,&search->pattern,actual_parameters,verify_candidates,search->mm_stack);
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
#ifdef GEM_PROFILE
  const bool already_mapped = matches_is_mapped(matches);
#endif
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_fprintf(stderr,">ASM::Adaptive Filtering (Boost)\n"); }
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  const uint64_t old_mcs = search->max_complete_stratum;
  // Region-Boost Profile
  approximate_search_generate_region_profile(search,region_profile_adaptive_boost,search->mm_stack);
  if (search->search_state == asearch_no_regions || search->search_state == asearch_exact_matches) {
    asearch_control_next_state(search,asearch_exact_filtering_boost,matches);
    PROF_STOP(GP_AS_FILTERING_EXACT_BOOST);
    return;
  }
  // Generate exact-candidates
  region_profile_schedule_filtering_fixed(&search->region_profile,
      ALL,REGION_FILTER_DEGREE_ZERO,parameters->filtering_threshold);
  approximate_search_generate_exact_candidates(search,matches);
  // Process candidates (just prepare to verification)
  filtering_candidates_process_candidates(search->filtering_candidates,
      search->archive,&search->pattern,actual_parameters,true,search->mm_stack);
  // Verify candidates
  approximate_search_verify_candidates(search,matches);
  if (!search->max_matches_reached && old_mcs > search->max_complete_stratum) {
    approximate_search_update_mcs(search,old_mcs);
  }
  search->search_state = asearch_exact_filtering_boost; // Correct State
  // Next State
  asearch_control_next_state(search,asearch_exact_filtering_boost,matches);
#ifdef GEM_PROFILE
  if (!already_mapped) PROF_ADD_COUNTER(GP_AS_FILTERING_EXACT_BOOST_MAPPED,matches_is_mapped(matches)?1:0);
#endif
  PROF_STOP(GP_AS_FILTERING_EXACT_BOOST);
}
GEM_INLINE void approximate_search_inexact_filtering(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_INEXACT);
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_fprintf(stderr,">ASM::Adaptive Filtering (Inexact)\n"); }
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  const search_parameters_t* const parameters = actual_parameters->search_parameters;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  // Region-Delimit Profile (Maximize number of unique-regions and try to isolate standard/repetitive regions)
  approximate_search_generate_region_profile(search,region_profile_adaptive_lightweight,search->mm_stack);
  if (search->region_profile.num_filtering_regions <= 1) {
    asearch_control_next_state(search,asearch_inexact_filtering,matches);
    PROF_STOP(GP_AS_FILTERING_INEXACT);
    return;
  }
  if (search->search_state == asearch_no_regions || search->search_state == asearch_exact_matches) {
    PROF_STOP(GP_AS_FILTERING_INEXACT);
    return;
  }
  // Generate exact-candidates (Dynamic filtering incorporated)
  const uint64_t proper_length = fm_index_get_proper_length(search->archive->fm_index);
  const uint64_t sensibility_misms_length = parameters->filtering_region_factor*proper_length;
  region_profile_schedule_filtering_adaptive(region_profile,search->max_complete_error,sensibility_misms_length);
  approximate_search_generate_inexact_candidates(search,true,true,matches);
  approximate_search_update_mcs(search,region_profile->errors_allowed + pattern->num_wildcards); // Update MCS
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
GEM_INLINE void approximate_search_unbounded_align(approximate_search_t* const search,matches_t* const matches) {
  PROF_START(GP_AS_FILTERING_UNBOUNDED_ALIGN);
#ifdef GEM_PROFILE
  const bool already_mapped = matches_is_mapped(matches);
#endif
  gem_cond_debug_block(DEBUG_SEARCH_STATE) { tab_fprintf(stderr,">ASM::Local-Alignment Filtering\n"); }
  // Unbounded-Align discarded candidates
  filtering_candidates_align_unbounded(search->filtering_candidates,
      search->archive->text,search->archive->locator,search->text_collection,
      &search->pattern,search->emulated_rc_search,search->as_parameters,
      matches,search->mm_stack);
  // Next State
  asearch_control_next_state(search,asearch_unbounded_alignment,matches);
#ifdef GEM_PROFILE
  if (!already_mapped) PROF_ADD_COUNTER(GP_AS_FILTERING_UNBOUNDED_ALIGN_MAPPED,matches_is_mapped(matches)?1:0);
#endif
  PROF_STOP(GP_AS_FILTERING_UNBOUNDED_ALIGN);
}
/*
 * Approximate Search based on Adaptive filtering
 */
GEM_INLINE void approximate_search_filtering_adaptive(approximate_search_t* const search,matches_t* const matches) {
  // Parameters
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
        PROF_START(GP_AS_READ_RECOVERY);
        approximate_search_exact_filtering_adaptive(search,region_profile_adaptive_recovery,verify_candidates,matches);
        PROF_STOP(GP_AS_READ_RECOVERY);
        if (!verify_candidates) return; // Return if no-filtering
        break;
      case asearch_exact_filtering_adaptive:
        // Exact-Filtering (Adaptive)
        PROF_START(GP_AS_FILTERING_EXACT);
        approximate_search_exact_filtering_adaptive(search,region_profile_adaptive_lightweight,verify_candidates,matches);
        PROF_STOP(GP_AS_FILTERING_EXACT);
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
        approximate_search_update_mcs(search,1);
        search->search_state = asearch_end;
        break;
      case asearch_unbounded_alignment: // Unbounded alignments
        approximate_search_unbounded_align(search,matches);
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
  // Parameters
  const as_parameters_t* const actual_parameters = search->as_parameters;
  search_parameters_t* const parameters = actual_parameters->search_parameters;
  archive_t* const archive = search->archive;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  // Exact Search
  if (search->max_complete_error==0) {
    approximate_search_neighborhood_exact_search(search,matches);
    return;
  }
  // Compute the region profile
  approximate_search_generate_region_profile(search,region_profile_adaptive_limited,search->mm_stack);
  if (search->search_state==asearch_no_regions) return;
  if (search->search_state==asearch_exact_matches) {
    matches_add_interval_match(matches,search->lo_exact_matches,
        search->hi_exact_matches,pattern->key_length,0,search->emulated_rc_search);
    return;
  }
  // Generate exact-candidates
  parameters->filtering_threshold = ALL; // No restriction
  region_profile_schedule_filtering_fixed(region_profile,ALL,REGION_FILTER_DEGREE_ZERO,parameters->filtering_threshold);
  approximate_search_generate_exact_candidates(search,matches);
  // Verify candidates
  filtering_candidates_t* const filtering_candidates = search->filtering_candidates;
  filtering_candidates_process_candidates(filtering_candidates,archive,
      pattern,actual_parameters,true,search->mm_stack);
  filtering_candidates_verify_candidates(
      filtering_candidates,search->archive,search->text_collection,
      pattern,actual_parameters,matches,search->mm_stack);
  // Align candidates
  filtering_candidates_align_candidates(filtering_candidates,
      search->archive->text,search->archive->locator,search->text_collection,
      pattern,search->emulated_rc_search,actual_parameters,false,matches,search->mm_stack);
  // Update MCS (maximum complete stratum)
  approximate_search_update_mcs(search,region_profile->errors_allowed + pattern->num_wildcards);
  PROF_STOP(GP_AS_FILTERING_EXACT);
}

