/*
 * PROJECT: GEM-Tools library
 * FILE: gt_map_metrics.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_MAP_METRICS_H_
#define GT_MAP_METRICS_H_

#include "gt_essentials.h"
#include "gt_map.h"

/*
 * Map Metrics :: Local(over first block) / Global(over all blocks)
 */
// Length of the mapping into the reference sequence (i.e. size over the chromosome)
GT_INLINE uint64_t gt_map_get_length(gt_map* const map);
GT_INLINE uint64_t gt_map_get_segment_length(gt_map* const map);
GT_INLINE uint64_t gt_map_get_global_length(gt_map* const map);
GT_INLINE uint64_t gt_map_get_global_base_length(gt_map* const map);
// GEM Distance (Number of Mismatches/Insert/Delete/Split operations)
GT_INLINE uint64_t gt_map_get_distance(gt_map* const map);
GT_INLINE uint64_t gt_map_get_segment_distance(gt_map* const map);
GT_INLINE uint64_t gt_map_get_global_distance(gt_map* const map);
// Bases aligned (Mismatches not included)
GT_INLINE uint64_t gt_map_get_bases_aligned(gt_map* const map);
GT_INLINE uint64_t gt_map_get_segment_bases_aligned(gt_map* const map);
GT_INLINE uint64_t gt_map_get_global_bases_aligned(gt_map* const map);
// Distance
GT_INLINE uint64_t gt_map_get_levenshtein_distance(gt_map* const map);
GT_INLINE uint64_t gt_map_get_segment_levenshtein_distance(gt_map* const map);
GT_INLINE uint64_t gt_map_get_global_levenshtein_distance(gt_map* const map);
GT_INLINE uint64_t gt_map_get_no_split_distance(gt_map* const map);
// SWG-Distance
GT_INLINE int64_t gt_map_get_swg_score(gt_map* const map);
/*
 * MMap/Vector Metrics
 */
GT_INLINE uint64_t gt_mmap_get_global_distance(gt_map** const mmap,const uint64_t num_blocks);
GT_INLINE uint64_t gt_mmap_get_global_levenshtein_distance(gt_map** const mmap,const uint64_t num_blocks);
GT_INLINE uint64_t gt_map_vector_get_length(gt_vector* const maps);
GT_INLINE uint64_t gt_map_vector_get_distance(gt_vector* const maps);
// Insert/Template Size
GT_INLINE int64_t gt_map_get_observed_template_size(gt_map* const map_a,gt_map* const map_b);

/*
 * Strict Map compare functions
 *   1.- Same sequence name
 *   2.- Same strand
 *   3.- Same number of blocks
 *   4.- for mapBlock in map {
 *         4.1 Same begin position
 *         4.2 Same end position
 *       }
 */
GT_INLINE int64_t gt_map_cmp(gt_map* const map_1,gt_map* const map_2);
/*
 * General Map compare function (differs from the standard one at @gt_map_cmp)
 *   1.- Same sequence name
 *   2.- Same strand
 *   3.- Same number of blocks
 *   4.- If (numBlocks == 1) {
 *           4.1 |begin_a-begin_b|<=range_tolerated || |end_a-end_b|<=range_tolerated
 *       } else {
 *         for mapBlock in map {
 *           5.1 |begin_a-begin_b|<=range_tolerated
 *           5.2 |end_a-end_b|<=range_tolerated
 *         }
 *       }
 */
GT_INLINE int64_t gt_map_range_cmp(gt_map* const map_1,gt_map* const map_2,const uint64_t range_tolerated);
/*
 * CIGAR Map compare function
 */
GT_INLINE int64_t gt_map_cmp_cigar(gt_map* const map_1,gt_map* const map_2);

/*
 * MMap compare functions (based on Map compare functions)
 */
GT_INLINE int64_t gt_mmap_cmp(gt_map** const map_1,gt_map** const map_2,const uint64_t num_maps);
GT_INLINE int64_t gt_mmap_range_cmp(gt_map** const map_1,gt_map** const map_2,const uint64_t num_maps,const uint64_t range_tolerated);
/*
 * As to resolve ties ...
 *   f(map_1,map_2) := (map_1 < map_2)
 */
GT_INLINE bool gt_map_less_than(gt_map* const map_1,gt_map* const map_2);

#endif /* GT_MAP_METRICS_H_ */
