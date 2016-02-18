/*
 * PROJECT: GEMMapper
 * FILE: region_profile.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef REGION_PROFILE_H_
#define REGION_PROFILE_H_

#include "utils/essentials.h"
#include "fm_index/fm_index.h"

/*
 * Constants
 */
#define REGION_MAX_LENGTH_PL_FACTOR 2.0

// Degree to filter region
#define REGION_FILTER_NONE        0
#define REGION_FILTER_DEGREE_ZERO 1
#define REGION_FILTER_DEGREE_ONE  2
#define REGION_FILTER_DEGREE_TWO  3

/*
 * Region Type:
 *   - region_gap: Constitutes a gap in the region profile (Eg. Contains wildcards).
 *   - region_standard: Region with some candidates. Should only be filtered up to 0-misms
 *   - region_unique: Very few candidates or none. Can be filtered up to n-errors (n>=0)
 */
typedef enum { region_unique, region_standard, region_gap } region_type;
typedef struct {
  uint64_t region_th;      // Max. number of candidates allowed per region
  uint64_t max_steps;      // Max. number of characters to explore to improve the region
  uint64_t dec_factor;     // Decreasing factor per step in region exploration
  uint64_t region_type_th; // Threshold to classify regions {ZERO,NON_ZERO}
} region_profile_model_t;
// Filtering regions
typedef struct {
  // Ranges of the region [begin,end)
  region_type type;
  uint64_t begin;
  uint64_t end;
  // Filtering error
  uint64_t degree; // Degree assigned to this region
  uint64_t max;
  uint64_t min;
  // Region exact-search (candidates)
  uint64_t hi;
  uint64_t lo;
} region_search_t;
// Region Profile
typedef struct {
  uint64_t id;
  uint64_t value;
} region_locator_t;
typedef struct {
  /* Regions */
  region_search_t* filtering_region; // Filtering regions
  uint64_t num_filtering_regions;    // Total number of filtering regions
  uint64_t num_standard_regions;     // Number of Standard Regions
  uint64_t num_unique_regions;       // Number of Unique Regions
  uint64_t num_zero_regions;         // Number of Zero candidate Regions
  /* Pattern Stats */
  uint64_t pattern_length;           // Length of the patter
  uint64_t errors_allowed;           // Total error allowed (the minimum required to get a novel match)
  // Profile Stats
  uint64_t total_candidates;         // Total number of candidates (from exact matching regions)
  uint64_t max_region_length;        // Largest region length
  double mappability_p;
  double mappability_2p;
  /* Locator (region sorting) */
  region_locator_t* loc;
} region_profile_t;

/*
 * Setup
 */
void region_profile_new(
    region_profile_t* const region_profile,
    const uint64_t pattern_length,
    mm_stack_t* const mm_stack);
void region_profile_clear(region_profile_t* const region_profile);

/*
 * Accessors
 */
uint64_t region_get_num_regions(region_profile_t* const region_profile);
bool region_profile_has_exact_matches(region_profile_t* const region_profile);

/*
 * Utils
 */
void region_profile_query_character(
    fm_index_t* const fm_index,
    rank_mquery_t* const rank_mquery,
    uint64_t* const lo,
    uint64_t* const hi,
    const uint8_t enc_char);
void region_profile_extend_last_region(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const bool* const allowed_enc,
    const uint64_t rp_region_type_th);

/*
 * Sort
 */
void region_profile_sort_by_estimated_mappability(region_profile_t* const region_profile);
void region_profile_sort_by_candidates(region_profile_t* const region_profile);

/*
 * Display
 */
void region_profile_print_region(
    FILE* const stream,
    region_search_t* const region,
    const uint64_t position);
void region_profile_print(
    FILE* const stream,
    const region_profile_t* const region_profile,
    const bool sorted);

/*
 * Iterator
 */
#define REGION_PROFILE_ITERATE(region_profile,region,position) \
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions; \
  region_search_t* region = region_profile->filtering_region; \
  uint64_t position; \
  for (position=0;position<num_filtering_regions;++position,++region)
#define REGION_LOCATOR_ITERATE(region_profile,region,position) \
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions; \
  region_search_t* const filtering_region = region_profile->filtering_region; \
  region_locator_t* const loc = region_profile->loc; \
  region_search_t* region; \
  uint64_t position; \
  for (position=0,region=filtering_region+loc[0].id; \
       position<num_filtering_regions; \
       ++position,region=filtering_region+loc[position].id)

#endif /* REGION_PROFILE_H_ */
