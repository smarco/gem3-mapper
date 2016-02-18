/*
 * PROJECT: GEMMapper
 * FILE: region_profile_schedule.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef REGION_PROFILE_SCHEDULE_H_
#define REGION_PROFILE_SCHEDULE_H_

#include "utils/essentials.h"
#include "filtering/region_profile.h"

/*
 * Region Scheduling
 */
void region_schedule_filtering_adaptive(
    region_search_t* const region,
    const uint64_t num_standard_regions_left,
    const uint64_t num_unique_regions_left,
    const uint64_t max_complete_error,
    const uint64_t sensibility_error_length,
    const uint64_t errors_allowed);

/*
 * Region Profile Scheduling
 */
void region_profile_schedule_filtering_fixed(
    region_profile_t* const region_profile,
    const uint64_t regions_required,
    const uint64_t filtering_degree,
    const uint64_t filtering_threshold);
void region_profile_schedule_filtering_adaptive(
    region_profile_t* const region_profile,
    const uint64_t max_complete_error,
    const uint64_t sensibility_misms_length);

/*
 * Display
 */
void region_profile_schedule_print(
    region_profile_t* const region_profile,
    const uint64_t max_differences,
    const uint64_t sensibility_error_length);


#endif /* REGION_PROFILE_SCHEDULE_H_ */
