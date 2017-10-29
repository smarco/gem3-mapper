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
 */

#ifndef MATCHES_CLASSIFY_H_
#define MATCHES_CLASSIFY_H_

#include "utils/essentials.h"
#include "archive/search/archive_search_se_parameters.h"
#include "matches/matches.h"

/*
 * Matches Classify
 */
void matches_classify(matches_t* const matches);

/*
 * Matches Accuracy Tests
 */
bool matches_classify_pattern_viable(pattern_t* const pattern);
bool matches_classify_max_matches_searched(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern);
bool matches_classify_min_depth_searched(
    matches_classification_t* const matches_classification);

/*
 * Alignment algorithm/paradigm fallback
 */
bool matches_classify_neighbourhood_fallback(
    matches_t* const matches,
    search_parameters_t* const search_parameters,
    pattern_t* const pattern);
bool matches_classify_local_alignment_fallback(
    matches_t* const matches,
    const local_alignment_t alignment_local);

/*
 * Adjusting search maximum error
 */
uint64_t matches_classify_adjust_max_error_by_strata_after_best(
    matches_t* const matches,
    const uint64_t max_search_error,
    const uint64_t complete_strata_after_best);
uint64_t matches_classify_compute_max_search_error(
    matches_t* const matches,
    pattern_t* const pattern,
    const uint64_t proper_length);

#endif /* MATCHES_CLASSIFY_H_ */
