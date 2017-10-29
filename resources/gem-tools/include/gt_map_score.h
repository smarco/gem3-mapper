/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_score.h
 * DATE: 03/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_MAP_SCORE_H_
#define GT_MAP_SCORE_H_

#include "gt_essentials.h"

#include "gt_alignment_utils.h"
#include "gt_template_utils.h"

/*
 * Map Score Accessors
 */
GT_INLINE uint64_t gt_map_get_score(gt_map* const map);
GT_INLINE void gt_map_set_score(gt_map* const map,const uint64_t score);
GT_INLINE uint8_t gt_map_get_phred_score(gt_map* const map);
GT_INLINE void gt_map_set_phred_score(gt_map* const map,const uint8_t phred_score);

/*
 * TODO
 * fx(MAP) Scoring functions
 * fx(ALIGNMENT,MAP) Scoring functions
 * fx(TEMPLATE) Scoring functions
 */

#endif /* GT_MAP_SCORE_H_ */
