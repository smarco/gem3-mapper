/*
 * PROJECT: GEMMapper
 * FILE: matches_classify.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCHES_CLASSIFY_H_
#define MATCHES_CLASSIFY_H_

#include "utils/essentials.h"
#include "matches/matches.h"
#include "matches/paired_matches.h"

/*
 * Classify
 */
matches_class_t matches_classify(matches_t* const matches);
paired_matches_class_t paired_matches_classify(paired_matches_t* const paired_matches);

#endif /* MATCHES_CLASSIFY_H_ */
