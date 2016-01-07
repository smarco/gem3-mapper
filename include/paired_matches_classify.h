/*
 * PROJECT: GEMMapper
 * FILE: paired_matches_classify.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef PAIRED_MATCHES_CLASSIFY_H_
#define PAIRED_MATCHES_CLASSIFY_H_

#include "essentials.h"
#include "matches_classify.h"

/*
 * PE Classify
 */
matches_class_t paired_matches_classify(paired_matches_t* const paired_matches);
void paired_matches_classify_compute_predictors(
    paired_matches_t* const paired_matches,matches_predictors_t* const predictors,
    const swg_penalties_t* const swg_penalties,const uint64_t total_read_length,
    const uint64_t max_region_length,uint64_t const proper_length,
    uint64_t const overriding_mcs,const uint64_t num_zero_regions);

#endif /* PAIRED_MATCHES_CLASSIFY_H_ */
