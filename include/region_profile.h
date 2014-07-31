/*
 * PROJECT: GEMMapper
 * FILE: region_profile.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef REGION_PROFILE_H_
#define REGION_PROFILE_H_

#include "essentials.h"

#include "fm_index.h"
#include "pattern.h"

#include "stats_vector.h"
#include "stats_matrix.h"

// Degree to filter region
#define REGION_FILTER_NONE 0
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
// Filtering regions
typedef struct {
  // Ranges of the region [end,start) (Yes!, backwards like the FM-search)
  uint64_t start;
  uint64_t end;
  region_type type;
  // Mismatches required to match the region
  uint64_t max;
  uint64_t min;
  // Region filtering results
  uint64_t hi;
  uint64_t lo;
  // Degree assigned to this region
  uint64_t degree;
} filtering_region_t;
// Mismatch regions
typedef struct {
  // Ranges of the region [end,start) (Yes!, backwards like the FM-search)
  uint64_t start;
  uint64_t end;
  // Degree assigned to this region
  uint64_t degree;
} mismatch_region_t;
// Region Profile
typedef struct {
  uint64_t id;
  uint64_t value;
} region_locator_t;
typedef struct {
  /* Pattern */
  uint64_t pattern_length;
  /* Filtering regions */
  filtering_region_t* filtering_region;
  uint64_t num_filtering_regions;  // Total number of regions := region_unique + region_standard + region_gap
  uint64_t num_standard_regions;   // Number of Standard Regions
  /* Mismatch regions */
  mismatch_region_t* mismatch_region;
  uint64_t num_mismatch_region;
  /* Region Partition Properties */
  uint64_t misms_required;              // Total mismatches required to get a new match
  /* Locator for region sorting */
  region_locator_t* loc;
} region_profile_t;
/*
 * Region Profile Stats
 */
typedef struct {

  uint64_t num_profiles;          // Number of samples
  uint64_t total_regions;         // Total number of regions extracted
  uint64_t total_unique_regions;  // Total number of unique regions extracted
  uint64_t total_candidates;      // Total number of candidates determined

  stats_vector_t* regions_per_read;     // Number of regions per read
  stats_vector_t* candidates_per_read;  // Total candidates per read

  stats_matrix_t* region_length_vs_region_candidates; // Number of candidates depending on the region length
} region_profile_aggregated_stats_t;
typedef struct {
  // TODO // TODO
  // TODO // TODO
} region_profile_stats_t;

/*
 * Setup
 */
GEM_INLINE void region_profile_new(
    region_profile_t* const region_profile,
    pattern_t* const pattern,mm_stack_t* const mm_stack);


/*
 * Accessors
 */
GEM_INLINE uint64_t region_get_num_regions(region_profile_t* const region_profile);
GEM_INLINE bool region_profile_has_exact_matches(region_profile_t* const region_profile);

/*
 * Region Profile Generation
 */
GEM_INLINE void region_profile_generate_fixed(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,pattern_t* const pattern,bool* const allowed_enc,
    const uint64_t num_regions);
GEM_INLINE void region_profile_generate_overlapped(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,pattern_t* const pattern,bool* const allowed_enc,
    const uint64_t region_length); // Much more todo
GEM_INLINE void region_profile_generate_adaptive(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,pattern_t* const pattern,const bool* const allowed_enc,
    const uint64_t rp_region_th,const uint64_t rp_max_steps,
    const uint64_t rp_dec_factor,const uint64_t rp_region_type_th,const uint64_t max_regions);
GEM_INLINE void region_profile_generate_full_progressive(
    region_profile_t* const region_profile,
    mismatch_region_t* const base_region,const uint64_t start_region,const uint64_t total_regions);

/*
 * Region Profile Utils
 */
GEM_INLINE void region_profile_sort_by_estimated_mappability(region_profile_t* const region_profile);
GEM_INLINE void region_profile_sort_by_candidates(region_profile_t* const region_profile);
GEM_INLINE void region_profile_fill_gaps(
    region_profile_t* const region_profile,const uint64_t eff_mismatches);
GEM_INLINE void region_profile_extend_first_region(
    region_profile_t* const region_profile,fm_index_t* const fm_index,const uint64_t rp_region_type_th);
GEM_INLINE void region_profile_extend_last_region(
    region_profile_t* const region_profile,fm_index_t* const fm_index,const uint64_t rp_region_type_th);

/*
 * Display
 */
GEM_INLINE void region_profile_print(FILE* const stream,const region_profile_t* const region_profile);

/*
 * Stats
 */
GEM_INLINE region_profile_stats_t* region_profile_stats_new();
GEM_INLINE void region_profile_stats_delete(region_profile_stats_t* const region_profile_stats);
// Record/Sample
GEM_INLINE void region_profile_stats_record(region_profile_stats_t* const region_profile_stats,region_profile_t* const region_profile);

/*
 * Iterator
 */
#define REGION_PROFILE_ITERATE(region_profile,region,position) \
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions; \
  filtering_region_t* region = region_profile->filtering_region; \
  uint64_t position; \
  for (position=0;position<num_filtering_regions;++position,++region)
#define REGION_LOCATOR_ITERATE(region_profile,region,position) \
  const uint64_t num_filtering_regions = region_profile->num_filtering_regions; \
  filtering_region_t* const filtering_region = region_profile->filtering_region; \
  region_locator_t* const loc = region_profile->loc; \
  filtering_region_t* region; \
  uint64_t position; \
  for (position=0,region=filtering_region+loc[0].id; \
       position<num_filtering_regions; \
       ++position,region=filtering_region+loc[position].id)



#endif /* REGION_PROFILE_H_ */
