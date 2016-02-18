/*
 * PROJECT: GEMMapper
 * FILE: matches_classify_logit_models.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MATCHES_CLASSIFY_LOGIT_MODELS_H_
#define MATCHES_CLASSIFY_LOGIT_MODELS_H_

#include "matches/matches_classify_logit.h"

/*
 * SE-Model DEFAULT
 */
extern const matches_classify_logit_model_t logit_model_single_end_default;

#define MATCHES_MIN_CI                 0.90
#define MATCHES_UNIQUE_CI              0.998
#define MATCHES_MMAPS_CI               0.998
#define MATCHES_TIES_CI                0.95

/*
 * PE-Model DEFAULT
 */
extern const matches_classify_logit_model_t logit_model_paired_end_default;

#define PAIRED_MATCHES_MIN_CI          0.90
#define PAIRED_MATCHES_UNIQUE_CI       0.99
#define PAIRED_MATCHES_MMAPS_CI        0.98
#define PAIRED_MATCHES_TIES_CI         0.95

#endif
