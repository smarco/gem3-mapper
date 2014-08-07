/*
 * PROJECT: GEMMapper
 * FILE: approximate_search.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search.h"

#include "quality_model.h"
#include "rank_mtable.h"
#include "neighborhood_search.h"

/*
 * Explicit search values
 */
// Unique
#define ALL_MAPPINGS UINT64_MAX
#define UNIQUE 1
// Filtering policies
#define DYNAMIC_SCHEDULING true
#define STATIC_SCHEDULING false
#define DYNAMIC_FILTERING true
#define STATIC_FILTERING false
// Mismatches to perform a Neighborhood Search
#define ZERO_ERRORS 0
#define ONE_ERROR   1
#define TWO_ERRORS  2
// Fast-mapping modes
#define FM_NO_FAST_MODE 0
#define FM_EXPLORING UINT64_MAX
#define FM_ADAPTIVE (UINT64_MAX-1)

/*
 * Region profile default parameters
 */
// Probing Scheme = (500,1,8,2)
#define PRP_REGION_THRESHOLD 1000
#define PRP_MAX_STEPS 3
#define PRP_DEC_FACTOR 3
#define PRP_REGION_TYPE_THRESHOLD 2
// Loose Scheme = (20,3,3,2)
#define SRP_REGION_THRESHOLD 20
#define SRP_MAX_STEPS 4
#define SRP_DEC_FACTOR 2
#define SRP_REGION_TYPE_THRESHOLD 2
// Tight Scheme = (50,7,3,2)
#define HRP_REGION_THRESHOLD 50
#define HRP_MAX_STEPS 10
#define HRP_DEC_FACTOR 4
#define HRP_REGION_TYPE_THRESHOLD 2
// Recovery Scheme = (20,4,2,2)
#define RRP_REGION_THRESHOLD 200
#define RRP_MAX_STEPS 1
#define RRP_DEC_FACTOR 8
#define RRP_REGION_TYPE_THRESHOLD 2
// Filtering thresholds
#define FILTERING_THRESHOLD 350
#define PA_FILTERING_THRESHOLD 2500
#define FILTERING_REGION_FACTOR ((double)1.0)

///////////////////////////////////////////////////////////////////////////////
// [Approximate Search Parameters]
///////////////////////////////////////////////////////////////////////////////
GEM_INLINE void approximate_search_initialize_replacements(approximate_search_parameters_t* const search_parameters) {
  // Reset
  memset(search_parameters->allowed_chars,0,256*sizeof(bool));
  memset(search_parameters->allowed_enc,0,DNA_EXT_RANGE*sizeof(bool));
  search_parameters->mismatch_alphabet[0] = DNA_CHAR_A;
  search_parameters->mismatch_alphabet[1] = DNA_CHAR_C;
  search_parameters->mismatch_alphabet[2] = DNA_CHAR_G;
  search_parameters->mismatch_alphabet[3] = DNA_CHAR_T;
  search_parameters->mismatch_alphabet_length = 4;
  search_parameters->allowed_chars[DNA_CHAR_A] = true;
  search_parameters->allowed_chars[DNA_CHAR_C] = true;
  search_parameters->allowed_chars[DNA_CHAR_G] = true;
  search_parameters->allowed_chars[DNA_CHAR_T] = true;
  search_parameters->allowed_chars[DNA_CHAR_N] = false;
  search_parameters->allowed_enc[ENC_DNA_CHAR_A] = true;
  search_parameters->allowed_enc[ENC_DNA_CHAR_C] = true;
  search_parameters->allowed_enc[ENC_DNA_CHAR_G] = true;
  search_parameters->allowed_enc[ENC_DNA_CHAR_T] = true;
  search_parameters->allowed_enc[ENC_DNA_CHAR_N] = false;
}
GEM_INLINE void approximate_search_parameters_init(approximate_search_parameters_t* const search_parameters) {
  /*
   * Initialize (DEFAULTS)
   */
  // Mapping strategy
  search_parameters->mapping_mode = mapping_adaptive_filtering;
  search_parameters->fast_mapping_degree = 0;
  // Qualities
  search_parameters->quality_model = quality_model_type_gem;
  search_parameters->quality_format = qualities_ignore;
  search_parameters->quality_threshold = 26;
  // Mismatch/Indels Parameters
  search_parameters->max_search_error = 0.04;
  search_parameters->max_filtering_error = 0.2;
  search_parameters->complete_strata_after_best = 0.0;
  search_parameters->min_matching_length = 0.2;
  // Matches search
  search_parameters->max_matches = ALL;
  // Replacements
  approximate_search_initialize_replacements(search_parameters);
  // Soft RP
  search_parameters->srp_region_th = SRP_REGION_THRESHOLD;
  search_parameters->srp_max_steps = SRP_MAX_STEPS;
  search_parameters->srp_dec_factor = SRP_DEC_FACTOR;
  search_parameters->srp_region_type_th = SRP_REGION_TYPE_THRESHOLD;
  // Hard RP
  search_parameters->hrp_region_th = HRP_REGION_THRESHOLD;
  search_parameters->hrp_max_steps = HRP_MAX_STEPS;
  search_parameters->hrp_dec_factor = HRP_DEC_FACTOR;
  search_parameters->hrp_region_type_th = HRP_REGION_TYPE_THRESHOLD;
  // Recover Read
  search_parameters->rrp_region_th = RRP_REGION_THRESHOLD;
  search_parameters->rrp_max_steps = RRP_MAX_STEPS;
  search_parameters->rrp_dec_factor = RRP_DEC_FACTOR;
  search_parameters->rrp_region_type_th = RRP_REGION_TYPE_THRESHOLD;
  // Filtering Thresholds
  search_parameters->filtering_region_factor = FILTERING_REGION_FACTOR;
  search_parameters->filtering_threshold = FILTERING_THRESHOLD;
  search_parameters->pa_filtering_threshold = PA_FILTERING_THRESHOLD;
  // Check alignments
  search_parameters->check_matches = check_none;
}
GEM_INLINE void approximate_search_configure_mapping_strategy(
    approximate_search_parameters_t* const search_parameters,
    const mapping_mode_t mapping_mode,const float mapping_degree) {
  search_parameters->mapping_mode = mapping_mode;
  search_parameters->fast_mapping_degree = mapping_degree;
}
GEM_INLINE void approximate_search_configure_quality_model(
    approximate_search_parameters_t* const search_parameters,
    const quality_model_t quality_model,const quality_format_t quality_format,const uint64_t quality_threshold) {
  search_parameters->quality_model = quality_model;
  search_parameters->quality_format = quality_format;
  search_parameters->quality_threshold = quality_threshold;
}
GEM_INLINE void approximate_search_configure_error_model(
    approximate_search_parameters_t* const search_parameters,
    float max_search_error,float max_filtering_error,
    float complete_strata_after_best,float min_matching_length) {
  search_parameters->max_search_error = max_search_error;
  search_parameters->max_filtering_error = max_filtering_error;
  search_parameters->complete_strata_after_best = complete_strata_after_best;
  search_parameters->min_matching_length = min_matching_length;
}
GEM_INLINE void approximate_search_configure_replacements(
    approximate_search_parameters_t* const search_parameters,
    char* const mismatch_alphabet,const uint64_t mismatch_alphabet_length) {
  // Reset
  approximate_search_initialize_replacements(search_parameters);
  // Filter replacements
  uint64_t i, count;
  for (i=0,count=0;i<mismatch_alphabet_length;i++) {
    if (is_dna(mismatch_alphabet[i])) {
      const char c = dna_normalized(mismatch_alphabet[i]);
      search_parameters->mismatch_alphabet[count] = c;
      search_parameters->allowed_chars[(uint8_t)c] = true;
      search_parameters->allowed_enc[dna_encode(c)] = true;
      ++count;
    }
  }
  gem_cond_fatal_error(count==0,ASM_REPLACEMENT_EMPTY);
  search_parameters->mismatch_alphabet_length = count;
}
GEM_INLINE void approximate_search_configure_matches(
    approximate_search_parameters_t* const search_parameters,const uint64_t max_matches) {
  search_parameters->max_matches = max_matches;
}
GEM_INLINE void approximate_search_instantiate_values(
    approximate_search_parameters_t* const search_parameters,const uint64_t pattern_length) {
  search_parameters->fast_mapping_degree_nominal = integer_proportion(search_parameters->fast_mapping_degree,pattern_length);
  search_parameters->max_search_error_nominal = integer_proportion(search_parameters->max_search_error,pattern_length);
  search_parameters->max_filtering_error_nominal = integer_proportion(search_parameters->max_filtering_error,pattern_length);
  search_parameters->complete_strata_after_best_nominal = integer_proportion(search_parameters->complete_strata_after_best,pattern_length);
  search_parameters->min_matching_length_nominal = integer_proportion(search_parameters->min_matching_length,pattern_length);
}
///////////////////////////////////////////////////////////////////////////////
// [Approximate Search]
///////////////////////////////////////////////////////////////////////////////
GEM_INLINE approximate_search_t* approximate_search_new(
    locator_t* const locator,graph_text_t* const graph,dna_text_t* const enc_text,fm_index_t* const fm_index,
    approximate_search_parameters_t* const search_parameters,mm_stack_t* const mm_stack) {
  // Allocate handler
  approximate_search_t* const approximate_search = mm_alloc(approximate_search_t);
  // Index Structures & Parameters
  approximate_search->locator = locator;
  approximate_search->graph = graph;
  approximate_search->enc_text = enc_text;
  approximate_search->fm_index = fm_index;
  approximate_search->search_parameters = search_parameters;
  // Search Auxiliary Structures
  filtering_candidates_new(&approximate_search->filtering_candidates); // Filtering Candidates
  // Interval Set TODO


  // MM
  approximate_search->mm_stack = mm_stack; // Memory stack
  // Return
  return approximate_search;
}
GEM_INLINE void approximate_search_clear(approximate_search_t* const approximate_search) {
  // Reset Approximate Search State
  approximate_search->current_search_stage = asearch_init; // Current Stage of the search
  approximate_search->stop_search_stage  = asearch_end; // Search stage to stop at
  approximate_search->max_differences = approximate_search->search_parameters->max_search_error_nominal;
  approximate_search->max_complete_stratum = ALL;
  approximate_search->max_matches_reached = false;
  // Filtering Candidates
  filtering_candidates_clear(&approximate_search->filtering_candidates);
  // Interval Set TODO



}
GEM_INLINE void approximate_search_delete(approximate_search_t* const approximate_search) {
  // Filtering Candidates
  filtering_candidates_delete(&approximate_search->filtering_candidates);
  // Free handler
  mm_free(approximate_search);
}
///////////////////////////////////////////////////////////////////////////////
// [Approximate Search Pattern]
///////////////////////////////////////////////////////////////////////////////
GEM_INLINE void approximate_search_prepare_pattern(
    approximate_search_t* const approximate_search,
    const approximate_search_parameters_t* const search_parameters,sequence_t* const sequence) {
  // Set quality search
  approximate_search->do_quality_search =
      (search_parameters->quality_format!=qualities_ignore) && sequence_has_qualities(sequence);
  // Get StackMem-allocator
  mm_stack_t* const mm_stack = approximate_search->mm_stack;
  // Allocate pattern memory
  const uint64_t read_length = sequence_get_length(sequence);
  pattern_t* const pattern = &approximate_search->pattern;
  pattern->key_length = read_length;
  pattern->key = mm_stack_calloc(mm_stack,read_length,uint8_t,false);
  // Build quality model & mask
  if (approximate_search->do_quality_search) {
    pattern->quality_mask =  mm_stack_calloc(mm_stack,read_length,uint8_t,false);
    quality_model(sequence,search_parameters->quality_model,
        search_parameters->quality_format,search_parameters->quality_threshold,pattern->quality_mask);
  } else {
    pattern->quality_mask = NULL;
  }
  /*
   * Check all characters in the key & encode key
   * Counts the number of wildcards(characters not allowed as replacements) & low_quality_bases
   */
  uint64_t i, num_wildcards=0, num_low_quality_bases=0;
  const char* const read = sequence_get_read(sequence);
  if (pattern->quality_mask == NULL) {
    for (i=0;i<read_length;++i) {
      const char character = read[i];
      if (!search_parameters->allowed_chars[(uint8_t)character]) ++num_wildcards;
      pattern->key[i] = dna_encode(character);
    }
  } else {
    for (i=0;i<read_length;++i) {
      const char character = read[i];
      if (!search_parameters->allowed_chars[(uint8_t)character]) {
        ++num_low_quality_bases; ++num_wildcards;
      } else if (pattern->quality_mask[i]!=qm_real) {
        ++num_low_quality_bases;
      }
      pattern->key[i] = dna_encode(character);
    }
  }
  pattern->num_wildcards = num_wildcards;
  pattern->num_low_quality_bases = num_low_quality_bases;
  // Calculate the effective number of differences
  const int64_t max_allowed_error = (int64_t)read_length - (int64_t) search_parameters->min_matching_length_nominal;
  uint64_t effective_filtering_max_error;
  if (gem_expect_false(max_allowed_error<=0)) { // Constrained by min_matching_length_nominal
    effective_filtering_max_error = 0;
  } else {
    // Constrained by num_low_quality_bases
    effective_filtering_max_error =
        search_parameters->max_filtering_error_nominal + approximate_search->pattern.num_low_quality_bases;
    if (effective_filtering_max_error > max_allowed_error) {
      effective_filtering_max_error = max_allowed_error;
    }
  }
  pattern->max_effective_filtering_error = effective_filtering_max_error;
  // Prepare region profile
  region_profile_new(&approximate_search->region_profile,read_length,mm_stack);
  // Prepare BPM pattern
  bpm_pattern_compile(&pattern->bpm_pattern,UINT128_SIZE,pattern->key,read_length,mm_stack);
}
///////////////////////////////////////////////////////////////////////////////
// [Approximate Search Control]
///////////////////////////////////////////////////////////////////////////////
/*
 * If delta parameter is set (and is below the maximum number of mismatches),
 * finds the minimum non zero stratum (mnzs) and adjusts
 * the maximum number of mismatches to (mnzs+delta)
 */
GEM_INLINE void approximate_search_adjust_max_differences_using_strata(
    approximate_search_t* const approximate_search,matches_t* const matches) {
  const uint64_t max_differences = approximate_search->max_differences;
  const uint64_t delta = approximate_search->search_parameters->complete_strata_after_best_nominal;
  if (delta < max_differences) {
    const int64_t fms = matches_counters_get_min_matching_stratum(matches);
    if (fms>=0 && fms+delta < max_differences) {
      approximate_search->max_differences = fms+delta;
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
// [Region Scheduling]
///////////////////////////////////////////////////////////////////////////////
GEM_INLINE void approximate_search_schedule_fixed_filtering_degree(
    region_profile_t* const region_profile,const uint64_t filtering_degree) {
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    region_locator_t* const loc = region_profile->loc+position;
    loc->id = position;
    loc->value = 0;
    region->min = filtering_degree;
  }
}
/*
 * PRE: (region_profile->num_filtering_regions <= max_mismatches)
 * Tries to assign the best possible filtering degree distribution among
 * all the regions to fulfill the requirements of a search up to @max_mismatches
 */
GEM_INLINE void approximate_search_schedule_filtering_degree(
    region_profile_t* const region_profile,
    const uint64_t max_differences,const uint64_t sensibility_misms_length) {
  const uint64_t num_regions = region_profile->num_filtering_regions;
  const uint64_t num_zregions = region_profile->num_standard_regions; // FIXME: Name
  const uint64_t num_nzregions = num_regions - num_zregions;

  // Pre-Sort regions
  region_profile_sort_by_estimated_mappability(region_profile);

  // Try to schedule a distribution of the errors over the regions
  uint64_t misms_nzregions = max_differences + 1 - num_zregions;
  uint64_t nzregions_left = num_nzregions;

  REGION_LOCATOR_ITERATE(region_profile,region,position) {
    if (region->type == region_unique) { // region_unique; can be filtered allowing errors
      const uint64_t length = region->start-region->end;
      if (nzregions_left >= misms_nzregions || length < sensibility_misms_length || max_differences==1) {
        region->min=REGION_FILTER_DEGREE_ZERO;
        misms_nzregions-=REGION_FILTER_DEGREE_ZERO;
      } else if (2*nzregions_left >= misms_nzregions || length < 2*sensibility_misms_length || max_differences==2) {
        region->min=REGION_FILTER_DEGREE_ONE;
        misms_nzregions-=REGION_FILTER_DEGREE_ONE;
      } else {
        region->min=REGION_FILTER_DEGREE_TWO;
        misms_nzregions-=REGION_FILTER_DEGREE_TWO;
      }
      --nzregions_left;

//      // TODO: Check this out !!!
//      if (region->hi-region->lo == 0 && nzregions_left>0){ // && region->min<DEGREE_TWO) {
//        ++(region->min);
//        --misms_nzregions;
//      }

    } else { // region_standard; Already has some/several matches
      region->min = REGION_FILTER_DEGREE_ZERO;
    }
  }

}
///////////////////////////////////////////////////////////////////////////////
// [Filtering Regions]
///////////////////////////////////////////////////////////////////////////////
GEM_INLINE void approximate_search_add_regions_to_filter(
    filtering_candidates_t* const filtering_candidates,
    const region_profile_t* const region_profile) {
  // Add all candidates from the hi/lo regions of the profile to filtering_candidates
  REGION_PROFILE_ITERATE(region_profile,region,position) {
    const uint64_t num_candidates = region->hi-region->lo;
    if (gem_expect_false(num_candidates==0)) continue;
    // Add region candidates (zero degree)
    filtering_candidates_add_interval(filtering_candidates,
        region->lo,region->hi,region->start,region->end,ZERO_ERRORS);
  }
}
/*
 * Filters all the regions up to the scheduled degree
 *  - Dynamic Scheduling: Assigns a filtering degree to each region as the search goes on
 *      (as opposed to a static scheduling giving in advanced)
 *  - Dynamic Filtering: Filters candidates per each queried region (thus, it can reduce
 *      the scope of the search by delta)
 *  - Locator guided: The region filtering is conducted by the locator order of region.
 *
 *  // TODO Needs a customed STATS OBJ
 */
GEM_INLINE void approximate_search_filter_regions(
    approximate_search_t* const search,
    const bool dynamic_scheduling,const uint64_t sensibility_error_length,
    const uint64_t filtering_threshold,matches_t* const matches) {
  const approximate_search_parameters_t* const parameters = search->search_parameters;
  region_profile_t* const region_profile = &search->region_profile;
  filtering_candidates_t* const filtering_candidates = &search->filtering_candidates;
  interval_set_t* const intervals_result = &search->intervals_result;
  const uint8_t* const key = search->pattern.key;

  const uint64_t num_regions = region_profile->num_filtering_regions;
  const uint64_t num_standard_regions = region_profile->num_standard_regions; // FIXME
  const uint64_t num_unique_regions = num_regions - num_standard_regions; // FIXME

  PROF_BEGIN(SC_FILTER_REGIONS);
  PROF_ADD_COUNTER(GSC_CAND_FILTER_REGIONS,MIN(search->max_differences+1,num_regions));

  // Dynamic Scheduling (Pre-Sort regions)
  if (dynamic_scheduling) region_profile_sort_by_estimated_mappability(&search->region_profile);

  uint64_t candidates, misms_required = 0;
  int64_t num_unique_regions_left = num_unique_regions;
  REGION_LOCATOR_ITERATE(region_profile,region,position) {
    PROF_INC_COUNTER(GSC_FILTER_REGIONS);
    bool perform_search = true;

    // Dynamic Schedule
    if (dynamic_scheduling) {
      int64_t misms_unique_regions = (int64_t)(search->max_differences+1) - (int64_t)(misms_required+num_standard_regions);
      if (region->type == region_unique) { // region_unique; can be filtered adding errors
        const uint64_t length = region->start-region->end;
        if (num_unique_regions_left >= misms_unique_regions || length < sensibility_error_length || search->max_differences==1) {
          region->min=REGION_FILTER_DEGREE_ZERO;
          misms_unique_regions-=REGION_FILTER_DEGREE_ZERO;
        } else if (2*num_unique_regions_left >= misms_unique_regions || length < 2*sensibility_error_length || search->max_differences==2) {
          region->min=REGION_FILTER_DEGREE_ONE;
          misms_unique_regions-=REGION_FILTER_DEGREE_ONE;
        } else {
          region->min=REGION_FILTER_DEGREE_TWO;
          misms_unique_regions-=REGION_FILTER_DEGREE_TWO;
        }
        --num_unique_regions_left;
      } else { // region_standard (already has some/several matches)
        region->min = REGION_FILTER_DEGREE_ZERO;
      }
    }

    /* Filter up to n errors (n>=2) */
    while (region->min >= REGION_FILTER_DEGREE_TWO) {
      PROF_INC_COUNTER(GSC_ATH_D2);
      const uint64_t max_error = region->min-1;
      if (perform_search) {
        interval_set_clear(intervals_result);
        neighborhood_search(key+region->end,region->start-region->end,max_error,intervals_result);
        perform_search = false;
      }
      candidates = interval_set_count_intervals_length(intervals_result);
      if (candidates <= filtering_threshold) {
        PROF_INC_COUNTER(GSC_ATH_D2_HIT);
        filtering_candidates_add_interval_set(filtering_candidates,intervals_result,region->start,region->end);
        misms_required+=region->min;
        break;
      } else {
        --(region->min);
      }
    }

    /* Filter up to 1 mismatches */
    if (region->min == REGION_FILTER_DEGREE_ONE) {
      PROF_INC_COUNTER(GSC_ATH_D1);
      if (perform_search) {
        interval_set_clear(intervals_result);
        neighborhood_search(key+region->end,region->start-region->end,ONE_ERROR,intervals_result);
//        num_intervals = fmi_base_one_mismatched_search(fmi,
//            key+region->end,region->start-region->end,
//            repls,repls_len,result_vector,mpool); // TODO
        perform_search = false;
      }
      candidates = interval_set_count_intervals_length_thresholded(intervals_result,ONE_ERROR);
      if (candidates <= filtering_threshold) {
        PROF_INC_COUNTER(GSC_ATH_D1_HIT);
        filtering_candidates_add_interval_set_thresholded(
            filtering_candidates,intervals_result,region->start,region->end,ONE_ERROR);
        misms_required+=REGION_FILTER_DEGREE_ONE;
      } else {
        --(region->min);
      }
    }

    /* Filter up to 0 mismatches */
    if (region->min == REGION_FILTER_DEGREE_ZERO) {
      PROF_INC_COUNTER(GSC_ATH_D0);
      filtering_candidates_add_interval(filtering_candidates,
          region->lo,region->hi,region->start,region->end,ZERO_ERRORS);
      misms_required+=REGION_FILTER_DEGREE_ZERO;
    }

    /* Otherwise, ignore the region (intractable) */

    /* Check candidates (Dynamic filtering)*/
    if (parameters->complete_strata_after_best_nominal < search->max_differences) {
      PROF_ADD_COUNTER(GSC_ATH_FILTER_CAND,filtering_candidates_get_pending_candidates(filtering_candidates));
      PROF_END(SC_FILTER_REGIONS);
      PROF_BEGIN(SC_ATH_FILTER);

      filtering_candidates_verify_pending(filtering_candidates,matches,true,true);
      approximate_search_adjust_max_differences_using_strata(search,matches);

      PROF_END(SC_ATH_FILTER);
      PROF_BEGIN(SC_FILTER_REGIONS);
      PROF_DEC_COUNTER(SC_FILTER_REGIONS);
    }

    /* Check misms condition */
    if (misms_required > search->max_differences) {
      PROF_ADD_COUNTER(GSC_SAVED_FILTER_REGIONS,num_regions-(position+1));
      break;
    }
  }
  // Set the minimum number of mismatches required
  region_profile->misms_required = misms_required;
}
///////////////////////////////////////////////////////////////////////////////
// [Mapping workflows 4.0]
///////////////////////////////////////////////////////////////////////////////
/*
 * If the number of wildcards (or errors required) is greater than the maximum number of
 * differences allowed, we try to recover as many matches as possible.
 * We extract feasible regions from the read) and filter them trying to recover anything
 * out of bad quality reads
 */
GEM_INLINE void approximate_search_read_recovery(
    approximate_search_t* const approximate_search,
    const uint64_t max_complete_stratum_value,matches_t* const matches) {
  approximate_search_parameters_t* const parameters = approximate_search->search_parameters;
  // Compute the region profile
  region_profile_generate_adaptive(
      &approximate_search->region_profile,approximate_search->fm_index,
      &approximate_search->pattern,approximate_search->search_parameters->allowed_enc,
      parameters->rrp_region_th,parameters->rrp_max_steps,
      parameters->rrp_dec_factor,parameters->rrp_region_type_th,
      approximate_search->max_differences+1);
  // region_profile_extend_last_region(approximate_search,parameters->rrp_region_type_th); // TODO

  PROF_BEGIN(GP_READ_RECOVERY);
  // Append the results to filter
  approximate_search_add_regions_to_filter(
      &approximate_search->filtering_candidates,&approximate_search->region_profile);

  // Verify !!
  filtering_candidates_verify_pending(&approximate_search->filtering_candidates,matches,true,true);
  PROF_END(GP_READ_RECOVERY);

  // Update MCS
  approximate_search->max_complete_stratum =
      MIN(max_complete_stratum_value,approximate_search->max_complete_stratum);
}
/*
 * Exact search
 */
GEM_INLINE void approximate_search_exact_search(
    approximate_search_t* const approximate_search,matches_t* const matches) {
  pattern_t* const pattern = &approximate_search->pattern;
  PROF_BEGIN(SC_EXACT_SEARCH);
  // FM-Index basic exact search
  uint64_t hi, lo;
  fm_index_bsearch(approximate_search->fm_index,pattern->key,pattern->key_length,&hi,&lo);
  // Add interval to matches
  matches_add_interval_match(matches,hi,lo,pattern->key_length,0,approximate_search->search_strand);
  // Update MCS
  approximate_search->max_complete_stratum = 1;
  PROF_END(SC_EXACT_SEARCH);
}
/*
 * Basic brute force search
 */
GEM_INLINE void approximate_search_basic(
    approximate_search_t* const approximate_search,matches_t* const matches) {
  approximate_search_parameters_t* const parameters = approximate_search->search_parameters;
  pattern_t* const pattern = &approximate_search->pattern;
  interval_set_t* const intervals_result = &approximate_search->intervals_result;
  // Basic search (Brute force mitigated by mrank_table)
  interval_set_clear(&approximate_search->intervals_result); // Clear
  neighborhood_search(pattern->key,pattern->key_length,parameters->max_search_error_nominal,intervals_result);
  // Add results
  matches_add_interval_set(matches,intervals_result);
  // Update MCS
  approximate_search->max_complete_stratum = parameters->max_search_error_nominal+1;
}
/*
 * Scout adaptive filtering
 *   Tries to guide a filtering search to extract regions from the key and filter
 *   each region up to a level where the number of candidates is low or moderate.
 *   It performs several profiles and verifications trying to get as close as possible
 *   to the required error degree (thus, not wasting resources in deeper searches than needed)
 *   [Sketch of the work-flow]
 *     1.- SoftRegionProfile
 *     2.- HardRegionProfile
 */
GEM_INLINE bool approximate_search_scout_adaptive_filtering(
    approximate_search_t* const search,matches_t* const matches) {
  fm_index_t* const fm_index = search->fm_index;
  approximate_search_parameters_t* const parameters = search->search_parameters;
  region_profile_t* const region_profile = &search->region_profile;
  pattern_t* const pattern = &search->pattern;

  //
  // Soft Region Profile
  //
  PROF_BEGIN(GP_REGION_PROFILE);
  region_profile_generate_adaptive(
      region_profile,fm_index,pattern,parameters->allowed_enc,
      parameters->srp_region_th,parameters->srp_max_steps,
      parameters->srp_dec_factor,parameters->srp_region_type_th,
      search->max_differences+1);
  PROF_END(GP_REGION_PROFILE);
  PROF_ADD_COUNTER(GP_NUM_SOFT_REGIONS,region_get_num_regions(region_profile));

  // FMI_SHOW_REGION_PROFILE(region_profile,"\nKEY:%s\n[STH]",key); TODO

  /*
   * Try to lower the search scope (max_differences) using @parameters->max_complete_strata
   */
  if (parameters->complete_strata_after_best_nominal < search->max_differences) {
    // Check if it has exact matches
    bool exact_matches = false;
    if (region_profile_has_exact_matches(region_profile)) {
      exact_matches = true;
    } else {
      region_profile_extend_first_region(region_profile,fm_index,parameters->srp_region_type_th);
      exact_matches = region_profile_has_exact_matches(region_profile);
    }
    // Try to lower @search->max_differences using @parameters->max_complete_strata
    if (exact_matches) {
      search->max_differences = parameters->complete_strata_after_best_nominal;
      if (parameters->complete_strata_after_best_nominal==0) {
        const filtering_region_t* const first_region = region_profile->filtering_region;
        matches_add_interval_match(matches,first_region->hi,first_region->lo,pattern->key_length,0,search->search_strand);
        search->max_differences = parameters->complete_strata_after_best_nominal;
        search->max_complete_stratum = 1;
        return true; // Well done !
      }
      return false;
    } else if (region_profile->num_filtering_regions <= 1) {

      // FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME
      if (pattern->num_wildcards > search->max_differences) {
        return false; // Why?? because it won't lower the delta
      }

//      // Hard probing
//      PROF_BEGIN(SC_REGION_PROFILE);
//      fmi_region_profile(fmi,key,key_len,search_params->allowed_chars,
//          PRP_REGION_THRESHOLD,PRP_MAX_STEPS,
//          PRP_DEC_FACTOR,PRP_REGION_TYPE_THRESHOLD,
//          search_params->max_mismatches+1,region_profile,mpool);
//      PROF_END(SC_REGION_PROFILE);
//      if (region_profile->num_filtering_regions > 0) {
//        PROF_BEGIN(SC_PROBING_DELTA);
//        fmi_region_profile_extend_last_region(fmi,key,key_len,
//            search_params->allowed_repl,PRP_REGION_TYPE_THRESHOLD,region_profile);
//        fmi_probe_matches(fmi,search_params,region_profile,true,matches,mpool);
//        PROF_END(SC_PROBING_DELTA);
//      }
//      FMI_CLEAR_PROFILE__RETURN_FALSE();
    }
  }

  // FIXME
  return false;

//  // Check that we have at least two regions
//  if (region_profile->num_filtering_regions == 0) return false; // Nothing to do here
//
//  // Zero-Filtering
//  if (region_profile->num_filtering_regions > search_params->max_mismatches) {
//    fmi_region_profile_extend_last_region(fmi,key,key_len,search_params->allowed_repl,
//        internals->srp_region_type_th,region_profile);
//    FMI_ZERO_FILTERING(); return true;
//  }
//
//  // Soft Probing the matches to lower max_mismatches using delta
//  if (search_params->delta < search_params->max_mismatches) {
//    PROF_BEGIN(SC_PROBING_DELTA);
//    fmi_probe_matches(fmi,search_params,region_profile,false,matches,mpool);
//    if (region_profile->num_filtering_regions > search_params->max_mismatches) {
//      PROF_INC_COUNTER(GSC_DELTA_PROBE_HIT); PROF_END(SC_PROBING_DELTA);
//      FMI_ZERO_FILTERING();
//      return true; // Well done !
//    }
//    PROF_END(SC_PROBING_DELTA);
//  }
//
//  //
//  // Hard Region Profile
//  //
//  PROF_BEGIN(SC_REGION_PROFILE);
//  fmi_region_profile(fmi,key,key_len,search_params->allowed_chars,
//      internals->hrp_region_th,internals->hrp_max_steps,
//      internals->hrp_dec_factor,internals->hrp_region_type_th,
//      search_params->max_mismatches+1,region_profile,mpool);
//  PROF_END(SC_REGION_PROFILE);
//  FMI_SHOW_REGION_PROFILE(region_profile,"[ATH]");
//  PROF_ADD_COUNTER(GSC_NUM_HARD_REGIONS,region_profile->num_filtering_regions);
//  if (region_profile->num_filtering_regions == 0) return false; // Nothing to do here
//  else if (region_profile->num_filtering_regions > 1) { // FIXME: OR max==1
//    fmi_region_profile_extend_last_region(fmi,key,key_len,
//        search_params->allowed_repl,internals->hrp_region_type_th,region_profile);
//  }
//
//  // Filtering regions (Dynamic)
//  PROF_BEGIN(SC_ATH_QUERY);
//  fmi_filter_regions(fmi,key,key_len,region_profile,DYNAMIC_SCHEDULING,
//      internals->filtering_region_factor*fmi->proper_length,
//      internals->filtering_threshold,search_params,matches,mpool);
//  FMI_SHOW_REGION_PROFILE(region_profile,"  --> [Schedule]");
//  PROF_END(SC_ATH_QUERY);
//
//  // Check candidates
//  PROF_ADD_COUNTER(GSC_ATH_FILTER_CAND,vector_get_used(mpool->fbuf1));
//  PROF_BEGIN(SC_ATH_FILTER);
//  fmi_matches_filter__append_decoded(fmi,matches,search_params,
//      mpool->fbuf1,mpool->ibuf1,mpool->ibuf2,true,true,mpool);
//  GEM_SWAP(mpool->ibuf1,mpool->ibuf2);
//  approximate_search_adjust_max_differences_using_strata(search_params,matches);
//  PROF_END(SC_ATH_FILTER);
//  vector_clean(mpool->fbuf1);
//
//  if (region_profile->misms_required > search_params->max_mismatches) {
//    PROF_INC_COUNTER(GSC_ATH_HIT);
//    return true;
//  } else {
//    return false;
//  }
}
///////////////////////////////////////////////////////////////////////////////
// [HighLevel Mapping workflows]
///////////////////////////////////////////////////////////////////////////////
/*
 * // A.K.A. brute-force mapping
 */
GEM_INLINE void approximate_search_neighborhood_search(approximate_search_t* const search,matches_t* const matches) {
  if (search->max_differences==0) {
    // Exact Search
    approximate_search_exact_search(search,matches);
    search->current_search_stage = asearch_end; // Update search state
    return; // End
  } else {
    // Basic brute force search
    approximate_search_basic(search,matches);
    search->current_search_stage = asearch_end; // Update search state
    return; // End
  }
}
/*
 * [GEM-workflow 4.0] Adaptive mapping (Formerly known as fast-mapping)
 *
 *   Filtering-only approach indented to adjust the degree of filtering w.r.t
 *   the structure of the read. Thus, in general terms, a read with many regions
 *   will enable this approach to align the read up to more mismatches than a read
 *   with less number of regions.
 *   Fast-mapping (in all its kinds) tries to detect the proper degree of filtering
 *   to achieve a compromise between speed and depth of the search (max_mismatches)
 *   It has different execution modes (@search_params->fast_mapping_degree)
 *     FM_STANDARD :: Using the max_mismatches value, extracts the regions from the read
 *                    and tries a filtering scheduling to adjust the search
 *                    up to that number of mismatches
 *     FM_ADAPTIVE :: Extracts the regions from the read but schedules each one of them
 *                    independently of the rest (and the maximum number of mismatches).
 *                    So each region is filtered up to the number of mismatches for which
 *                    the number of candidates is below the threshold
 *     FM_DEGREE={1..n} :: Every region is scheduled up to the given degree (no matter
 *                         the region structure). EXACT_FILTER=1, ONE_MISM=2, etc.
 *
 *   NOTE: Some parts are similar (or equal) to others in the regular workflow.
 *         However we don't unify the code because as this approach is constantly evolving
 *         we want to keep it separate to keep the main algorithm isolated
 */
GEM_INLINE void approximate_search_adaptive_mapping(approximate_search_t* const search,matches_t* const matches) {
  const approximate_search_parameters_t* const parameters = search->search_parameters;
  fm_index_t* const fm_index = search->fm_index;
  pattern_t* const pattern = &search->pattern;
  region_profile_t* const region_profile = &search->region_profile;
  // Compute the region profile
  PROF_BEGIN(GP_REGION_PROFILE);
  region_profile_generate_adaptive(
      region_profile,fm_index,pattern,parameters->allowed_enc,
      parameters->srp_region_th,parameters->srp_max_steps,
      parameters->srp_dec_factor,parameters->srp_region_type_th,
      search->max_differences+1);
  PROF_END(GP_REGION_PROFILE);
  PROF_ADD_COUNTER(GP_NUM_SOFT_REGIONS,region_get_num_regions(region_profile));
  // Zero-region reads
  const uint64_t num_wildcards = search->pattern.num_wildcards;
  if (region_profile->num_filtering_regions==0) {
    search->max_complete_stratum = num_wildcards;
    return;
  }
  // Extend last region
  region_profile_extend_last_region(region_profile,fm_index,pattern,parameters->allowed_enc,parameters->srp_region_type_th);
  region_profile_print(stderr,region_profile,false,false); // DEBUG
  if (region_profile_has_exact_matches(region_profile)) {
    const filtering_region_t* const first_region = region_profile->filtering_region;
    matches_add_interval_match(matches,first_region->hi,first_region->lo,pattern->key_length,0,search->search_strand);
    search->max_differences = parameters->complete_strata_after_best_nominal;
    search->max_complete_stratum = 1;
    return;
  }
  // Process the proper fast-mapping scheme
  const uint64_t sensibility_error_length =
      parameters->filtering_region_factor*fm_index_get_proper_length(search->fm_index);
  switch (parameters->fast_mapping_degree_nominal) {
    case FM_ADAPTIVE:
      approximate_search_filter_regions(search,DYNAMIC_SCHEDULING,
          sensibility_error_length,parameters->filtering_threshold,matches);
      break;
    case FM_EXPLORING:
      approximate_search_schedule_filtering_degree(region_profile,UINT32_MAX,sensibility_error_length);
      approximate_search_filter_regions(search,STATIC_SCHEDULING,
          sensibility_error_length,parameters->filtering_threshold,matches);
      break;
    default: // FM_DEGREE={0..n}
      approximate_search_schedule_fixed_filtering_degree(region_profile,parameters->fast_mapping_degree+1);
      approximate_search_filter_regions(search,STATIC_SCHEDULING,
          sensibility_error_length,parameters->filtering_threshold,matches);
      break;
  }
  // Filter the candidates
  PROF_ADD_COUNTER(GSC_ATH_FILTER_CAND,filtering_candidates_get_pending_candidates(&search->filtering_candidates));
  PROF_INC_COUNTER(GSC_ATH_HIT); PROF_BEGIN(SC_ATH_FILTER);
  filtering_candidates_verify_pending(&search->filtering_candidates,matches,true,true);
  PROF_END(SC_ATH_FILTER);
  // Set the
  if (!search->max_matches_reached) {
    // Update MCS (maximum complete stratum)
    const int64_t max_complete_stratum = region_profile->misms_required + num_wildcards;
    search->max_complete_stratum = MIN(max_complete_stratum,search->max_complete_stratum);
    // TODO fast-maping-stats
    // if (max_complete_stratum<20) PROF_INC_COUNTER(GSC_FAST_MAPPING_MCS+max_complete_stratum);
  } else {
    search->max_complete_stratum = 0;
  }
}
/*
 * [GEM-workflow 4.0] Incremental mapping (A.K.A. Complete mapping)
 */
GEM_INLINE void approximate_search_incremental_mapping(
    approximate_search_t* const approximate_search,matches_t* const matches) {
  approximate_search_exact_search(approximate_search,matches);
  approximate_search->current_search_stage = asearch_end;
  return; // FIXME FIXME

    // TODO uncomment
//  approximate_search_parameters_t* const parameters = approximate_search->search_parameters;
//  pattern_t* const pattern = &approximate_search->pattern;
//
//  // Check if recovery is needed
//  if (pattern->num_wildcards >= parameters->max_differences) {
//    approximate_search_read_recovery(approximate_search,pattern->num_wildcards,matches);
//    // Update search state
//    approximate_search->current_search_stage = asearch_end;
//    return;
//  }
//
//  // Exact search
//  if (parameters->max_differences==0) {
//    approximate_search_exact_search(approximate_search,matches);
//    // Update search state
//    approximate_search->current_search_stage = asearch_end;
//    return;
//  }
//
//  // Very short reads (Neighborhood search)
//  if (pattern->key_length <= RANK_MTABLE_SEARCH_DEPTH ||
//      pattern->key_length < approximate_search->fm_index->proper_length) {
//    PROF_BEGIN(SC_SMALL_READS);
//    approximate_search_basic(approximate_search,matches);
//    // Update search state
//    approximate_search->current_search_stage = asearch_end;
//    PROF_END(SC_SMALL_READS);
//    return;
//  }

  // TODO Implement
//  // Adaptive Regions Filtering
//  PROF_BEGIN(SC_ATH);
//  if (fmi_adaptive_region_filtering(fmi,key,key_len,search_params,&region_profile,matches,mpool)) {
//    // TODO: This should be adaptive, not fixed wrt region_profile.misms_required
//    fm_matches_update_max_complete_stratum(matches,search_params->max_mismatches+search_params->num_wildcards+1);
//    PROF_END(SC_ATH);
//    return true; // Well done!
//  }
//  PROF_END(SC_ATH);
//
//  /* After filtering, if the required level of complete maximum stratum is not reached,
//   * we fulfill the search by using a progressive scheme up to the maximum number of
//   * mismatches and applying the constraints derived from the first step (Adaptive Filtering) */
//
//  // Check wildcards
//  const uint64_t misms_required = region_profile->misms_required+num_wildcards;
//  if (misms_required > approximate_search->max_differences) {
//    approximate_search_read_recovery(approximate_search,misms_required,matches);
//    // Update search state
//    approximate_search->current_search_stage = asearch_end;
//    return;
//  }

//  // Quick search quit shorcut
//  if (search_params->internal_parameters.quit_if_hard) {
//    fm_matches_recover_state(matches); // FIXME: Not needed, just wipe out everything
//    return false;
//  }

//  // Shortcut
//  fm_matches_recover_state(matches);
//  if (num_wildcards>0) return true;
//  vector* results = fmi_base_mismatched_search_pure_TLS(
//      fmi,key,key_len,search_params->max_mismatches,repls,repls_len,mpool);
//  fm_matches_append_search_results(search_params,matches,results,0);
//  fm_matches_update_max_complete_stratum(matches,search_params->max_mismatches+1);

//  // Progressive Adaptive Regions Filtering
//  PROF_BEGIN(SC_PROGRESSIVE);
//  fm_matches_update_max_complete_stratum(matches,search_params->max_mismatches+1);
//  fmi_progressive_adaptive_region_filtering(fmi,key,
//      key_len,search_params,&region_profile,matches,mpool);
//  PROF_END(SC_PROGRESSIVE);
}
/*
 * Approximate String Matching using the FM-index.
 */
GEM_INLINE void approximate_search(approximate_search_t* const approximate_search,matches_t* const matches) {
  PROF_BEGIN(GP_ASM);
  const approximate_search_parameters_t* const search_parameters = approximate_search->search_parameters;
  // Check if all characters are wildcards
  if (approximate_search->pattern.key_length==approximate_search->pattern.num_wildcards) {
    approximate_search->max_complete_stratum = approximate_search->pattern.key_length;
    approximate_search->current_search_stage = asearch_end;
    return;
  }
  /*
   * Select mapping strategy
   */
  switch (search_parameters->mapping_mode) {
    case mapping_neighborhood_search:
      approximate_search_neighborhood_search(approximate_search,matches); // A.K.A. brute-force mapping
      break;
    case mapping_incremental_mapping:
      approximate_search_incremental_mapping(approximate_search,matches); // A.K.A. complete-mapping
      break;
    case mapping_adaptive_filtering:
      approximate_search_adaptive_mapping(approximate_search,matches);    // A.K.A. fast-mapping
      break;
    case mapping_fast:
      GEM_NOT_IMPLEMENTED();
      break;
    case mapping_fixed_filtering:
      GEM_NOT_IMPLEMENTED();
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROF_END(GP_ASM);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//
///*
// * Filters all the regions in the profile up to zero mismatches
// * From the results tries to reduce the maximum number of mismatches
// * using delta.
// */
//GEM_INLINE void fmi_probe_matches(const _FMI_* const fmi,
//    fmi_search_parameters* const search_params,region_profile* const region_profile,
//    const bool sparse_probing,matches* const matches,vector_pool* const mpool) {
//  // We add all the exact search of the region to the filter query buffer
//  vector* const query_buffer = mpool->fbuf1;
//  vector_init(query_buffer,filter_query);
//
//  // Add candidates to filter
//  if (!sparse_probing) {
//    approximate_search_add_regions_to_filter(region_profile,query_buffer);
//  } else {
//    fmi_append_sampled_region_zero_degree_to_filter(region_profile,
//        search_params->internal_parameters.filtering_threshold,query_buffer);
//  }
//
//  // Filter all the exact regions
//  PROF_ADD_COUNTER(GSC_DELTA_PROBE_CAND,vector_get_used(query_buffer));
//  fmi_matches_filter__append_decoded(fmi,matches,search_params,
//      query_buffer,mpool->ibuf1,mpool->ibuf2,true,false,mpool);
//
////  // Adjust the maximum number of mismatches
////  approximate_search_adjust_max_differences_using_strata(search_params,matches); // TODO; Skip as done internally
//
//  // Clean matches buffers
//  vector_clean(mpool->ibuf1);
//  vector_clean(mpool->ibuf2);
//  vector_clean(mpool->fbuf1);
//}
//
//
//
///*
// * Schedules and filters (max_mismatches+1) regions (degree zero)
// */
//#define FMI_ZERO_FILTERING()
//  fmi_schedule_region_zero_degree_filtering(region_profile,search_params->max_mismatches+1,mpool);
//  fmi_filter_regions(fmi,key,key_len,region_profile,STATIC_SCHEDULING,
//      internals->filtering_region_factor*fmi->proper_length,
//      internals->filtering_threshold,search_params,matches,mpool);
//  PROF_ADD_COUNTER(GSC_ZERO_FILTER_CAND,vector_get_used(mpool->fbuf1)); PROF_BEGIN(SC_ZERO_FILTER);
//  fmi_matches_filter__append_decoded(fmi,matches,search_params,
//      mpool->fbuf1,mpool->ibuf1,mpool->ibuf2,true,true,mpool);
//  PROF_END(SC_ZERO_FILTER)
//#define FMI_CLEAR_PROFILE__RETURN_FALSE()
//region_profile->num_filtering_regions = 0; return false
//
//
//
//
//
//
//
//
//
//#define FMI_ASIGN_DEGREE(i,start,end,remaining_misms) {
//  bool all_wildcards = true;
//  int64_t pos = start-1;
//  all_wildcards = true;
//  while (all_wildcards && pos>=end) all_wildcards = !allowed_chars[key[pos--]];
//  if (!all_wildcards) {
//    misms_region[i].degree = DEGREE_ZERO;
//    --remaining_misms;
//  } else {
//    misms_region[i].degree = DONT_USE;
//  }
//}
//
//GEM_INLINE void fmi_assign_region_degree(
//    const ch_t* const key,const uint64_t key_len,
//    fmi_search_parameters* const search_params,const uint64_t eff_mismatches,
//    region_profile* region_profile,vector_pool* const mpool) {
//  FMI_SHOW_REGION_PROFILE(region_profile, "[PA-RegProf]");
//  if (region_profile->num_filtering_regions==1) PROF_INC_COUNTER(GSC_ONE_REGION);
//  if (region_profile->num_filtering_regions==0) PROF_INC_COUNTER(GSC_ZERO_REGION);
//
//  // Fill regions profile gaps
//  fmi_fill_region_profile_gaps(key,key_len,eff_mismatches,search_params,region_profile,mpool);
//
//  // Adjust degree of regions
//  filtering_region* const filtering_reg = region_profile->filtering_region;
//  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
//  const uint64_t misms_base = eff_mismatches-region_profile->misms_required;
//  uint64_t i;
//  for (i=0; i<num_filtering_regions; ++i) {
//    filtering_reg[i].degree = filtering_reg[i].min;
//    filtering_reg[i].max = misms_base + filtering_reg[i].min;
//  }
//
//  // Distribute the remaining mismatches (Flood)
//  int64_t misms_left = (int64_t)(eff_mismatches+1) - (int64_t)(region_profile->misms_required);
//  uint64_t level = 0;
//  while (misms_left > 0) {
//    int64_t j;
//    ++level;
//    for (j=0; j<num_filtering_regions && misms_left>0; ++j) {
//      if (filtering_reg[j].degree<level &&
//          ((filtering_reg[j].type!=GAP_REGION && filtering_reg[j].start-filtering_reg[j].end >= level) ||
//           (filtering_reg[j].type==GAP_REGION && FMI_REGION_GET_NUMBER_BASES(filtering_reg+j) >= level)) ) {
//        ++(filtering_reg[j].degree);
//        --misms_left;
//      }
//    }
//  }
//  FMI_SHOW_REGION_PROFILE(region_profile, "  -->[PA-RegProf-Degree]");
//
//  // Assign the mismatch regions
//  mismatch_region* misms_region = region_profile->mismatch_region;
//  uint64_t num_misms_region;
//  for (i=0,num_misms_region=0; i<num_filtering_regions; ++i) {
//    const uint64_t reg_degree = filtering_reg[i].degree;
//    if (reg_degree == 0) {
//      misms_region->start = filtering_reg[i].start;
//      misms_region->end = filtering_reg[i].end;
//      misms_region->degree = DONT_USE;
//      ++num_misms_region; ++misms_region;
//    } else {
//      const uint64_t length_region = filtering_reg[i].start-filtering_reg[i].end;
//      const uint64_t size_chunk = length_region/reg_degree;
//      uint64_t last_pos = filtering_reg[i].start, j;
//      for (j=1; j<=reg_degree; ++j, ++misms_region) {
//        misms_region->start = last_pos;
//        last_pos -= size_chunk;
//        misms_region->end = (j<reg_degree) ? last_pos : filtering_reg[i].end;
//        misms_region->degree = DEGREE_ZERO;
//      }
//      num_misms_region += reg_degree;
//    }
//  }
//  region_profile->num_mismatch_region = num_misms_region;
//  FMI_SHOW_MISMS_PROFILE(region_profile, "[PA-MismsProf]");
//}
//
//
//
///*
// *
// */
//GEM_INLINE void fmi_progressive_adaptive_region_filtering(
//    const _FMI_* const fmi,const ch_t* const key,const uint64_t key_len,
//    fmi_search_parameters* const search_params,region_profile* region_profile,
//    matches* const matches,vector_pool* const mpool) {
//  slch_t* repls = search_params->repls;
//  const uint64_t repls_len = search_params->repls_len;
//  const bool* const allowed_repl = search_params->allowed_repl;
//  const uint64_t num_wildcards = search_params->num_wildcards;
//  const uint64_t eff_max_mismatches = search_params->max_mismatches-num_wildcards;
//  const uint64_t unique = search_params->unique_mapping;
//
//  // Assign the mismatches left to the regions.
//  fmi_assign_region_degree(key,key_len,
//      search_params,eff_max_mismatches,region_profile,mpool);
//
//  // Iterate over the regions and perform each progressive search
//  vector* exclusion_set = mpool->buf3;
//  vector* result_set = mpool->rbuf1;
//  uint64_t i, j, candidates, accum;
//
//  PROF_BEGIN(SC_PA_FAIL);
//  PROF_BEGIN(SC_PAF_QUERY);
//  vector_init(exclusion_set,interval_t);
//  vector_init(result_set,interval_t);
//  filtering_region* filtering_reg = region_profile->filtering_region;
//  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
//  mismatch_region* const misms_regions = region_profile->mismatch_region;
//  const uint64_t num_mismatch_region = region_profile->num_mismatch_region;
//  bool keep_searching = true;
//  for (i=0, accum=0; keep_searching && i<num_filtering_regions; ++i) {
//    for (j=0; keep_searching && j<filtering_reg[i].degree; ++j) {
//      const uint64_t reg_num = accum+j;
//      // Assign the PA scheme
//      region_profile->mismatch_region = misms_regions+reg_num;
//      region_profile->num_mismatch_region = num_mismatch_region-reg_num;
//      // Check validity of scheme
//      if (region_profile->mismatch_region->start==filtering_reg[i].start &&
//          filtering_reg[i].degree<=filtering_reg[i].min) continue;
//      // Perform the search
//      if (reg_num==0) { // First scheme saved as exclusion set
//        fmi_base_mismatched_region_search(fmi,key,key_len,eff_max_mismatches,
//            repls,repls_len,allowed_repl,region_profile,exclusion_set,mpool);
//        fm_matches_append_search_results(search_params,matches,exclusion_set,num_wildcards);
//        if (FMI_MATCHES_GET_NUM_MATCHES(matches)>unique) { PROF_END(SC_PAF_QUERY);
//          fm_matches_update_max_complete_stratum(matches,0);
//          return;
//        }
//      } else {
//        // Perform the progressive search // TODO: Predictor of PA-fails
//        vector_clean(result_set);
//        fmi_base_mismatched_region_search(fmi,key,region_profile->mismatch_region->start,
//            eff_max_mismatches,repls,repls_len,allowed_repl,region_profile,result_set,mpool);
//        // Subtract the matches previously added
//        interval_set_subtract(result_set,exclusion_set);
//        COUNT_CANDIDATES(candidates,result_set);
//        if (candidates <= search_params->internal_parameters.pa_filtering_threshold) {
//          RESERVE__ADD_SET_INTERVALS_TO_FILTER_QUERIES(
//              result_set,0,vector_get_used(result_set),
//              candidates,region_profile->mismatch_region->start,0);
//          // Join the new matches to the previously found
//          interval_set_union(exclusion_set,result_set);
//        } else { // Skip to TLS approach using full PA profile
//          PROF_BEGIN(SC_PA_FULL);
//          if (reg_num<10) { PROF_INC_COUNTER(GSC_PA_FAIL_LEVEL+reg_num); START_TIMER(TSC_PA_FAIL_LEVEL+reg_num);}
//
//          fmi_generate_full_progressive_profile(region_profile,misms_regions,reg_num,num_mismatch_region);
//          FMI_SHOW_MISMS_PROFILE(region_profile, "---> Full[i=%lu]",reg_num);
//          vector_clean(result_set);
//          fmi_base_mismatched_region_search(fmi,key,key_len,eff_max_mismatches,
//              repls,repls_len,allowed_repl,region_profile,result_set,mpool);
//          interval_set_subtract(result_set,exclusion_set);
//          fm_matches_append_search_results(search_params,matches,result_set,num_wildcards);
//
//          PROF_ADD_COUNTER(GSC_PA_INTERVALS+reg_num,result_set->used); COUNT_CANDIDATES(candidates,result_set);
//          PROF_END(SC_PA_FULL); PROF_END(SC_PA_FAIL);
//          if (reg_num<10) { STOP_TIMER(TSC_PA_FAIL_LEVEL+reg_num); PROF_ADD_COUNTER(GSC_PA_FAIL_INT+reg_num,candidates); }
//          keep_searching = false; break;
//        }
//      }
//    }
//    accum += filtering_reg[i].degree;
//  }
//  PROF_END(SC_PAF_QUERY);
//
//  // Filter all
//  PROF_ADD_COUNTER(GSC_PAF_FILTER_CAND, vector_get_used(mpool->fbuf1)); PROF_BEGIN(SC_PAF_FILTER);
//  fmi_matches_filter__append_decoded(fmi,matches,search_params,
//      mpool->fbuf1,mpool->ibuf1,mpool->ibuf2,false,true,mpool);
//  PROF_END(SC_PAF_FILTER);
//}
//
//
//
//
//
//GEM_INLINE void fmi_append_sampled_region_zero_degree_to_filter(
//    region_profile* const region_profile,const uint64_t filtering_threshold,
//    vector* const vector_buffer) {
//  const filtering_region* const region = region_profile->filtering_region;
//  const uint64_t num_filtering_regions = region_profile->num_filtering_regions;
//  uint64_t i, total_candidates;
//  for (i=0,total_candidates=0; i<num_filtering_regions; ++i) {
//    total_candidates+=region[i].hi-region[i].lo;
//  }
//  // Append sampled candidates
//  const uint64_t sampling_gap = (total_candidates+filtering_threshold)/filtering_threshold;
//  uint64_t total_used = vector_get_used(vector_buffer);
//  for (i=0; i<num_filtering_regions; ++i) {
//    const uint64_t count = region[i].hi-region[i].lo;
//    if (__builtin_expect(count==0,false)) continue;
//    uint64_t num_sampled = (count+sampling_gap)/sampling_gap;
//    vector_reserve(vector_buffer,total_used+num_sampled);
//    filter_query* queries_to_filter = (filter_query*)vector_get_mem(vector_buffer)+total_used;
//    uint64_t pos;
//    for (num_sampled=0,pos=region[i].lo; pos<region[i].hi; pos+=sampling_gap) {
//      queries_to_filter->start = region[i].start;
//      queries_to_filter->end = region[i].end;
//      queries_to_filter->misms = 0;
//      queries_to_filter->end_reg_idx = pos;
//      ++queries_to_filter; ++num_sampled;
//    }
//    total_used+=num_sampled;
//  }
//  vector_set_used(vector_buffer,total_used);
//}
//
//
////GEM_TOPLEVEL_INLINE uint64_t fmi_mismatched_search_extend(
////    const _FMI_* const fmi,fmi_extend_parameters* const extend_parameters,
////    const idx_t init_position,const idx_t end_position,
////    matches* const matches,vector_pool* const mpool) {
////  const uint64_t key_len = extend_parameters->key_len;
////
////  // Check wildcards and bad-quality bases
////  uint64_t total_wildcards=0, bad_quality_bases=0;
////  if (!fmi_check_key_alphabet__count_wildcards(fmi,
////      extend_parameters->forward_key,key_len,
////      extend_parameters->mismatch_mask,extend_parameters->allowed_chars,
////      &total_wildcards,&bad_quality_bases)) {
////    return 0;
////  }
////
////  // Perform the extended search
////  const uint64_t proposed_max_distance = extend_parameters->max_distance; // +bad_quality_bases;
////  const uint64_t max_align_diff = key_len-extend_parameters->min_anchor_size;
////  const uint64_t max_distance = (proposed_max_distance>max_align_diff)?max_align_diff:proposed_max_distance;
////  PROF_BEGIN(SC_ONLINE_ASM_SEARCH);
////  const uint64_t found_matches =
////      fmi_matches_search(fmi,extend_parameters,init_position,end_position,
////        max_distance,extend_parameters->max_extended_matches,
////        extend_parameters->allowed_repl,matches,mpool);
////  PROF_END(SC_ONLINE_ASM_SEARCH);
////  return found_matches;
////}
//
//
//
//
///*
//
//#define HARDCODED_KEY_MAX_LENGHT_FOR_CHECKING 1000
//#define UPPER_BOUND_MISMS 200
//
//#define FMI_CHECK_GET_NEXT_MISMS()
//  if (misms_num > 0) {
//    misms_pos = misms->position;
//    misms_char = misms->mismatch;
//    --misms_num;
//    ++misms;
//  } else {
//    misms_pos = -1;
//  }
//
//bool fmi_check_alignment_given_mismatches(ch_t* key,ch_t* decoded,uint64_t key_len,
//    uint64_t misms_num,mismatch* misms,uint64_t position,uint64_t max_indel_length,
//    direction_t direction,bool dump_results) {
//  const uint64_t init_misms_num = misms_num;
//  char scratch_key[HARDCODED_KEY_MAX_LENGHT_FOR_CHECKING];
//  char scratch_ref[HARDCODED_KEY_MAX_LENGHT_FOR_CHECKING];
//  char operations[HARDCODED_KEY_MAX_LENGHT_FOR_CHECKING];
//
//  bool bad_result = false; // Be optimistic !
//  uint64_t misms_pos, misms_char;
//  uint64_t total_misms;
//  int j, k, z, t;
//
//  // Init the buffers
//  memset(scratch_key, ' ', HARDCODED_KEY_MAX_LENGHT_FOR_CHECKING);
//  memset(scratch_ref, ' ', HARDCODED_KEY_MAX_LENGHT_FOR_CHECKING);
//
//  // Begin the check
//  FMI_CHECK_GET_NEXT_MISMS();
//  total_misms = 0;
//  for (j=0,k=0,t=0; j<key_len;) {
//    if (j==misms_pos) { // Check misms
//      if (misms_char >= 256) { // Indel (we believe the indel, let's see what happens)
//        uint64_t indl = misms_char / 256;
//        uint64_t size = indl/2;
//        if ((indl % 2) == 0) {
//          for (z=k; z<k+size; z++) { // Insertion in key
//            scratch_key[t] = ' ';
//            scratch_ref[t] = decoded[z];
//            operations[t] = '-'; t++;
//          }
//          k+=size;
//        } else {
//          for (z=j; z<j+indl/2; z++) { // Deletion in key
//            scratch_key[t] = key[z];
//            scratch_ref[t] = ' ';
//            operations[t] = '-'; t++;
//          }
//          j+=size;
//        }
//      } else {
//        if (key[j] == misms_char || misms_char != decoded[k]) {
//          bad_result=true;
//        }
//        scratch_key[t] = key[j];
//        scratch_ref[t] = decoded[k];
//        operations[t] = misms_char; //'+';
//        j++; k++; t++;
//      }
//      total_misms++;
//      FMI_CHECK_GET_NEXT_MISMS(); // Next misms
//    } else { // No misms (j!=misms_pos)
//      if (key[j] != decoded[k]) {
//        bad_result=true;
//        operations[t] = 'X';
//      } else {
//        operations[t] = '|';
//      }
//      scratch_key[t] = key[j];
//      scratch_ref[t] = decoded[k];
//      j++;k++;t++;
//    }
//  }
//  scratch_key[t] = 0;
//  operations[t] = 0;
//  scratch_ref[t] = 0;
//
//  bool error = bad_result || total_misms!=init_misms_num;
//  if (error || dump_results) {
//    fprintf(stderr, "\n%s", error?"[FAIL]":"[OK]");
//    fprintf(stderr, "[%ld-Correctness][F-KEY] %s. \t {Position %lu,%s}  {MismsFound=%lu.MismsExpected=%lu}\n",
//        id_checked, key, position, direction==Forward ? "FORWARD": "REVERSE", total_misms, init_misms_num);
//    fprintf(stderr, "  [KEY] %s \n", scratch_key);
//    fprintf(stderr, "        %s \n", operations);
//    fprintf(stderr, "  [REF] %s \n", scratch_ref);
//  }
//
//  return bad_result;
//}
//
//bool fmi_check_alignment_given_max_mismatches(ch_t* key,ch_t* decoded,uint64_t key_len,
//    uint64_t position,uint64_t max_mismatches,direction_t direction,bool dump_results) {
//  char operations[HARDCODED_KEY_MAX_LENGHT_FOR_CHECKING];
//  uint64_t i, num_misms;
//
//  for (i=0,num_misms=0; i<key_len; ++i) {
//    if (key[i]==decoded[i]) {
//      operations[i]='|';
//    } else {
//      operations[i]=decoded[i];
//      num_misms++;
//    }
//  }
//  operations[i] = 0;
//
//  bool error = num_misms>max_mismatches;
//  if (error || dump_results) {
//    fprintf(stderr, "\n%s", error?"[FAIL]":"[OK]");
//    fprintf(stderr, "[%ld-Correctness][I-KEY] %s. \t {Position %lu,%s} {MismsFound=%lu.MAXMisms=%lu}\n",
//        id_checked, key, position, direction==Forward ? "FORWARD": "REVERSE", num_misms, max_mismatches);
//    fprintf(stderr, "  [KEY] %s \n", key);
//    fprintf(stderr, "        %s \n", operations);
//    fprintf(stderr, "  [REF] %s \n", decoded);
//    return false;
//  } else {
//    return true;
//  }
//}
//
//void fmi_check_pos_correctness(
//    const _FMI_* const fmi,fmi_search_parameters* const search_params,
//    matches* const matches,bool dump_results,vector_pool* const mpool) {
//  ch_t decoded[HARDCODED_KEY_MAX_LENGHT_FOR_CHECKING];
//  keys_info* key_info;
//  ch_t* key;
//
//  // Check the correctness of the mappings decoded
//  uint64_t misms_num, i;
//  mismatch* misms;
//  pos_match_t* positions = (pos_match_t*)vector_get_mem(matches->rbuf_pos);
//  for (i=0; i<vector_get_used(matches->rbuf_pos); ++i, ++positions) {
//    // Decode the aligned key from the index
//    const uint64_t len = search_params->key_len+2*(search_params->max_indel_len+search_params->max_differences);
//    const uint64_t text_length = fmi->text_length;
//    int64_t to_decode = (positions->position+len>text_length?text_length-positions->position:len);
//    fmi_decode(fmi,positions->position,to_decode,decoded);
//
//    // Retrieve the searched key
//    key_info = (keys_info*)vector_get_mem(matches->qbuf_keys_info) + positions->key_id;
//    key = (ch_t*)vector_get_mem(matches->qbuf_keys_buffer) + key_info->displacement;
//    misms_num = positions->mismatches;
//    misms = (mismatch*)vector_get_mem(matches->rbuf_mismatches)+positions->displacement;
//    // Check the correctness of the alignment
//    fmi_check_alignment_given_mismatches(key,decoded,search_params->key_len,
//        misms_num,misms,positions->position,search_params->max_indel_len,
//        key_info->direction,dump_results);
//  }
//}
//void fmi_check_int_correctness(
//    const _FMI_* const fmi,fmi_search_parameters* const search_params,
//    matches* const matches,bool dump_results,vector_pool* const mpool) {
//  ch_t decoded[HARDCODED_KEY_MAX_LENGHT_FOR_CHECKING];
//  keys_info* key_info;
//  ch_t* key;
//
//  // Check the correctness of the mapping not yet decoded (as intervals)
//  uint64_t i, pos;
//  int_match_t* interval_match = vector_get_mem(matches->rbuf_int);
//  for (i=0; i<vector_get_used(matches->rbuf_int); ++i, ++interval_match) {
//    for (pos=interval_match->lo; pos<interval_match->hi; ++pos) {
//      // Decode the aligned key from the index
//      uint64_t position = fmi_lookup(fmi,pos);
//      fmi_decode(fmi,position,search_params->key_len,decoded);
//      // Retrieve the searched key
//      key_info = (keys_info*)vector_get_mem(matches->qbuf_keys_info) + interval_match->key_id;
//      key = (ch_t*)vector_get_mem(matches->qbuf_keys_buffer) + key_info->displacement;
//      // Check the correctness of the alignment
//      fmi_check_alignment_given_max_mismatches(key,decoded,search_params->key_len,
//          position,search_params->max_mismatches,key_info->direction,dump_results);
//    }
//  }
//}
//
//#define SHOW_TLS_COUNTERS(max_misms) {
//  uint64_t l;
//  fprintf(stderr, "TLS[%ld", (int64_t)cntrs_real[0]);
//  for (l=1; l<=max_misms; ++l) {
//    fprintf(stderr, ":%ld", (int64_t)cntrs_real[l]);
//  }
//  fprintf(stderr, "]\n");
//}
//#define SHOW_MATCHES_DIFF_FULL(max_misms) {
//  uint64_t l;
//  fprintf(stderr, "TLS[%ld", (int64_t)cntrs_real[0]);
//  for (l=1; l<=max_misms; ++l) {
//    fprintf(stderr, ":%ld", (int64_t)cntrs_real[l]);
//  }
//  fprintf(stderr, "] ");
//  fprintf(stderr, "- GEM[%ld", (int64_t)cntrs_matches[0]);
//  for (l=1; l<=max_misms; ++l) {
//    fprintf(stderr, ":%ld", (int64_t)cntrs_matches[l]);
//  }
//  fprintf(stderr, "] ");
//  fprintf(stderr, "= [%ld", (int64_t)cntrs_real[0] - (int64_t)cntrs_matches[0]);
//  for (l=1; l<=max_misms; ++l) {
//    fprintf(stderr, ":%ld", (int64_t)cntrs_real[l] - (int64_t)cntrs_matches[l]);
//  }
//  fprintf(stderr, "] \n");
//}
//
//
//void fmi_check_completeness(const _FMI_* const fmi,fmi_search_parameters* const search_params,
//                            matches* const matches,bool dump_results,vector_pool* const mpool) {
//  char scratch_key[HARDCODED_KEY_MAX_LENGHT_FOR_CHECKING];
//  char representative_key[HARDCODED_KEY_MAX_LENGHT_FOR_CHECKING];
//  const uint64_t key_len = search_params->key_len;
//  slch_t* repls = search_params->repls;
//  uint64_t repls_len = search_params->repls_len;
//  uint64_t eff_max_misms;
//  bool weak_passed = true;
//  uint64_t cntrs_real[UPPER_BOUND_MISMS];
//  int64_t i, j;
//  uint64_t cntrs_matches[UPPER_BOUND_MISMS];
//  vector *results = vector_new(1000, sizeof(interval_t)), *aux;
//  interval_t *ints;
//
//  // Set the proper max_misms
//  if (search_params->fast_mapping_degree>0) {
//    if (matches->max_complete_stratum==0) return;
//    eff_max_misms = matches->max_complete_stratum-1;
//    // Fill counters
//    uint64_t *matches_count = vector_get_mem(matches->rbuf_counts);
//    for (i=eff_max_misms; i>=0; --i) {
//      cntrs_real[i] = 0;
//      if (i < vector_get_used(matches->rbuf_counts)) {
//        cntrs_matches[i] = matches_count[i];
//      } else {
//        cntrs_matches[i] = 0;
//      }
//      if (cntrs_matches[i]>0 && i+search_params->delta < search_params->max_mismatches) {
//        eff_max_misms = i+search_params->delta;
//      }
//    }
//  } else {
//    eff_max_misms = search_params->max_mismatches;
//    // Fill counters
//    uint64_t *matches_count = vector_get_mem(matches->rbuf_counts);
//    for (i=search_params->max_mismatches; i>=0; --i) {
//      cntrs_real[i] = 0;
//      if (i < vector_get_used(matches->rbuf_counts)) {
//        cntrs_matches[i] = matches_count[i];
//      } else {
//        cntrs_matches[i] = 0;
//      }
//      if (cntrs_matches[i]>0 && i+search_params->delta < search_params->max_mismatches) {
//        eff_max_misms = i+search_params->delta;
//      }
//    }
//  }
//
//  // Run simple TLS and count the number of matches
//  const uint64_t max_misms_used = (search_params->fast_mapping_degree>0) ?
//      matches->max_complete_stratum-1 : search_params->max_mismatches;
//  vector_clean(results);
//  keys_info* key_info = vector_get_mem(matches->qbuf_keys_info);
//  for (i=0; i<UNSAFE_MIN(2,vector_get_used(matches->qbuf_keys_info)); ++i, ++key_info) {
//    ch_t* key = (ch_t*)vector_get_mem(matches->qbuf_keys_buffer) + key_info->displacement;
//    strncpy(scratch_key,(char*)key,key_len); scratch_key[key_len]=0;
//    if (i==0) {
//      strncpy(representative_key,(char*)key,key_len); representative_key[key_len]=0;
//    }
//    aux = fmi_base_mismatched_search_pure_TLS(fmi,(ch_t*)scratch_key,
//        key_len,max_misms_used,repls,repls_len,mpool);
//    vector_reserve(results, vector_get_used(results)+vector_get_used(aux));
//    ints = vector_get_mem_next_elm(results,interval_t);
//    INTERVAL_ITERATE(aux) {
//      cntrs_real[interval->misms]+=(interval->hi-interval->lo);
//      *ints = *interval;
//      ++ints;
//    } END_INTERVAL_ITERATE;
//    vector_update_used(results,ints);
//   //  SHOW_TLS_COUNTERS(max_misms_used);
//  }
//
//  // Dump counters if requested
//  if (dump_results) {
//    SHOW_MATCHES_DIFF_FULL(eff_max_misms);
//    fprintf(stderr, "[+scope]"); SHOW_MATCHES_DIFF_FULL(eff_max_misms+40);
//  }
//
//  // WEAK CHECK: The number of matches must be the same
//  for (i=0; i<=eff_max_misms; ++i) {
//    // Count the number of matches at each level with indels
//    uint64_t num_indel_matchs = 0, p, q;
//    pos_match_t* positions = (pos_match_t*)vector_get_mem(matches->rbuf_pos);
//    for (p=0; p < vector_get_used(matches->rbuf_pos); ++p, ++positions) {
//      mismatch *misms = (mismatch*)vector_get_mem(matches->rbuf_mismatches)+positions->displacement;
//      uint64_t misms_num = positions->mismatches;
//      if (misms_num!=i) continue;
//      for (q=0; q<misms_num; ++q) {
//        if (misms[q].mismatch >= 256) {
//          ++num_indel_matchs;
//          break;
//        }
//      }
//    }
//
//    if (cntrs_real[i] != cntrs_matches[i]-num_indel_matchs) {
//      fprintf(stderr, "[FAIL][%ld-Completeness][KEY] %s Misms Level %lu. Found %lu vs Real %lu \n",
//          id_checked,representative_key,i,cntrs_matches[i],cntrs_real[i]);
//      SHOW_MATCHES_DIFF_FULL(eff_max_misms);
//      weak_passed = false;
//    } else {
////      fprintf(stderr, "      [%ld-Completeness][KEY] %s. Misms Level %lu. Found %lu vs Real %lu \n",
////          id_checked,representative_key,i,cntrs_matches[i],cntrs_real[i]);
//    }
//  }
//
//  // STRONG CHECK: Check that every result from the check TLS search is in matches
//  if (weak_passed) {
//    INTERVAL_ITERATE(results) {
//      if (interval->misms > eff_max_misms) continue;
//      idx_t hi, lo, pos;
//      hi = interval->hi;
//      lo = interval->lo;
//      for (i=lo; i<hi; i++) {
//        bool found = false;
//
//        // Look it up into the intervals buffer
//        int_match_t* ints_found = vector_get_mem(matches->rbuf_int);
//        uint64_t p;
//        idx_t mismatches;
//        for (p=0;p<vector_get_used(matches->rbuf_int);++ints_found,++p) {
//          if (ints_found->lo <= i && i < ints_found->hi) {
//            found = true; mismatches = ints_found->mismatches; break;
//          }
//        }
//
//        if (!found) {
//          // Look it up into the positions buffer
//          pos = fmi_lookup(fmi, i);
//          // Check the position into the matches
//          pos_match_t* positions = (pos_match_t*)vector_get_mem(matches->rbuf_pos);
//          for (j=0; j<vector_get_used(matches->rbuf_pos); ++j, ++positions) {
//            if (positions->position == pos) {
//              found = true; mismatches = positions->mismatches;
//              if (mismatches == interval->misms) break;
//            }
//          }
//        }
//
//        if (!found) { // Not found... Too bad
//          fprintf(stderr, "[FAIL][%ld-Completeness][KEY] %s. RealMatchSA[%lu] NOT IN FoundMatches. ",id_checked,representative_key,i);
//          SHOW_MATCHES_DIFF_FULL(eff_max_misms);
//        } else if (mismatches != interval->misms) { // Check misms
//          fprintf(stderr, "[FAIL][%ld-Completeness][KEY] %s. RealMatchSA[m=%lu] vs FoundMatches[m=%lu]. ",
//              id_checked,representative_key,interval->misms,mismatches);
//          SHOW_MATCHES_DIFF_FULL(eff_max_misms);
//        }
//      }
//    } END_INTERVAL_ITERATE;
//  } else {
//    fprintf(stderr, "[FAIL] Completeness. Strong check is not carried out because of previous errors \n");
//  }
//
//  // Free vector
//  vector_delete(results);
//
//  // Output final count and done
//  if (id_checked % 1000 == 0) fprintf(stderr, "Check until read %lu... \n", id_checked);
//
//  return;
//}
//
//
//*/
//
//
//
