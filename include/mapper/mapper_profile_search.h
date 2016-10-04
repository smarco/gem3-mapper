/*
 * PROJECT: GEMMapper
 * FILE: mapper_profile_search.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:  // TODO
 */

#ifndef MAPPER_PROFILE_SEARCH_H_
#define MAPPER_PROFILE_SEARCH_H_

#include "utils/essentials.h"
#include "mapper/mapper_profile_counters.h"

/*
 * Region Profile
 */
void mapper_profile_print_region_profile(FILE* const stream);

/*
 * Candidates Generation
 */
void mapper_profile_print_candidate_generation(FILE* const stream);

/*
 * Candidate Verification
 */
void mapper_profile_print_candidate_verification(FILE* const stream);

/*
 * Candidate realign
 */
void mapper_profile_print_candidate_realign(FILE* const stream);

/*
 * Neighborhood Search
 */
void mapper_profile_print_neighborhood_search(FILE* const stream);

/*
 * Approximate Search
 */
void mapper_profile_print_approximate_search(FILE* const stream);

/*
 * Approximate Search Profile Summary
 */
void mapper_profile_print_approximate_search_summary(FILE* const stream);

#endif /* MAPPER_PROFILE_SEARCH_H_ */
