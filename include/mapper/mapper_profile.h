/*
 * PROJECT: GEMMapper
 * FILE: mapper_profile.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:  // TODO
 */

#ifndef MAPPER_PROFILE_H_
#define MAPPER_PROFILE_H_

#include "utils/essentials.h"
#include "system/profiler.h"
#include "mapper_profile_counters.h"


/*
 * Mapper SE
 */
void mapper_profile_print_mapper_se(
    FILE* const stream,const bool map_output,const uint64_t num_threads);

/*
 * Mapper PE
 */
void mapper_profile_print_mapper_pe(
    FILE* const stream,const bool map_output,const uint64_t num_threads);

#endif /* MAPPER_PROFILE_H_ */
