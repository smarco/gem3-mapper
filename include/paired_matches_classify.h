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
 * Constants
 */
#define PAIRED_MATCHES_MIN_CI          0.5
#define PAIRED_MATCHES_UNIQUE_CI       0.998
#define PAIRED_MATCHES_MMAPS_CI        0.995
#define PAIRED_MATCHES_TIES_CI         0.90

/*
 * PE Classify
 */
matches_class_t paired_matches_classify(paired_matches_t* const paired_matches);
void paired_matches_classify_compute_predictors(
    paired_matches_t* const paired_matches,matches_predictors_t* const predictors,
    const swg_penalties_t* const swg_penalties,const uint64_t total_read_length,
    const uint64_t max_region_length,uint64_t const proper_length,
    uint64_t const overriding_mcs,const uint64_t num_zero_regions);

double paired_matches_classify_unique(matches_predictors_t* const predictors);
double paired_matches_classify_mmaps(matches_predictors_t* const predictors);
double paired_matches_classify_ties(matches_predictors_t* const predictors);

#endif /* PAIRED_MATCHES_CLASSIFY_H_ */
