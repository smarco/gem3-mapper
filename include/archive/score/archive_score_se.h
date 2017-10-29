/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Archive MAPQ-scoring module produces MAPQ scores for SE-searches
 */


#ifndef ARCHIVE_SCORE_SE_H_
#define ARCHIVE_SCORE_SE_H_

#include "archive/search/archive_search.h"
#include "matches/matches.h"
#include "matches/classify/matches_predictors.h"

/*
 * Archive Scoring Utils
 */
uint8_t archive_score_probability_scale(
    const double probability,
    const double sum_probability,
    const uint8_t floor,
    const uint8_t ceil);
uint64_t archive_score_probability_to_mapq(
    const double probability);

/*
 * Archive Scoring Logit Classes
 */
uint8_t archive_score_matches_se_logit_unique(
    search_parameters_t* const search_parameters,
    matches_predictors_t* const matches_predictors,
    matches_classification_t* const matches_classification);
uint8_t archive_score_matches_se_logit_mmap(
    search_parameters_t* const search_parameters,
    matches_predictors_t* const matches_predictors,
    matches_classification_t* const matches_classification);

/*
 * Archive Scoring SE
 */
void archive_score_matches_se(
    archive_search_t* const archive_search,
    matches_t* const matches);

#endif /* ARCHIVE_SCORE_SE_H_ */
