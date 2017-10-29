/*
 * PROJECT: GEM-Tools library
 * FILE: gt_alignment_utils.h
 * DATE: 19/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */


#ifndef GT_ALIGNMENT_UTILS_H_
#define GT_ALIGNMENT_UTILS_H_

#include "gt_essentials.h"
#include "gt_alignment.h"
#include "gt_counters_utils.h"
#include "gt_map_align.h"

#include "gt_output_map.h"
#include "gt_input_parser.h"

/*
 * Alignment basic tools
 */
GT_INLINE uint64_t gt_alignment_get_read_proportion(gt_alignment* const alignment,const float proportion);

/*
 * Alignment's Maps operators (Update global state: counters, ...)
 */
GT_INLINE gt_map* gt_alignment_put_map(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,
    gt_map* const map,const bool check_duplicated,const bool replace_duplicated);

GT_INLINE void gt_alignment_insert_map(gt_alignment* const alignment,gt_map* const map);
GT_INLINE void gt_alignment_insert_map_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,gt_map* const map);
GT_INLINE void gt_alignment_insert_map_gt_vector(gt_alignment* const alignment,gt_vector* const map_vector);
GT_INLINE void gt_alignment_insert_map_fx_gt_vector(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,gt_vector* const map_vector);

GT_INLINE bool gt_alignment_find_map_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,gt_map* const map,
    uint64_t* const found_map_pos,gt_map** const found_map);
GT_INLINE bool gt_alignment_is_map_contained(gt_alignment* const alignment,gt_map* const map);
GT_INLINE bool gt_alignment_is_map_contained_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_alignment* const alignment,gt_map* const map);

GT_INLINE void gt_alignment_reduce_maps(gt_alignment* const alignment,const uint64_t max_num_matches);

/*
 * Alignment's Counters operators
 */
GT_INLINE bool gt_alignment_is_mapped(gt_alignment* const alignment);
GT_INLINE bool gt_alignment_is_thresholded_mapped(gt_alignment* const alignment,const int64_t max_allowed_strata);
GT_INLINE void gt_alignment_recalculate_counters(gt_alignment* const alignment);
GT_INLINE void gt_alignment_recalculate_counters_no_splits(gt_alignment* const alignment);

GT_INLINE int64_t gt_alignment_get_uniq_degree(gt_alignment* const alignment);
GT_INLINE int64_t gt_alignment_get_min_matching_strata(gt_alignment* const alignment);
GT_INLINE bool gt_alignment_get_next_matching_strata(
    gt_alignment* const alignment,const uint64_t begin_strata,
    uint64_t* const next_matching_strata,uint64_t* const num_maps);

/*
 * Alignment's Maps Sorting
 */
GT_INLINE void gt_alignment_sort_by_distance__score(gt_alignment* const alignment);
GT_INLINE void gt_alignment_sort_by_distance__score_no_split(gt_alignment* const alignment);

/*
 * Alignment's Maps Utils
 */
GT_INLINE uint64_t gt_alignment_sum_mismatch_qualities(gt_alignment* const alignment,gt_map* const map);
GT_INLINE uint64_t gt_alignment_get_max_mismatch_quality(gt_alignment* const alignment);

/*
 * Alignment's Maps set-operators
 */
GT_INLINE void gt_alignment_merge_alignment_maps(gt_alignment* const alignment_dst,gt_alignment* const alignment_src);
GT_INLINE void gt_alignment_merge_alignment_maps_fx(
    int64_t (*gt_map_cmp)(gt_map*,gt_map*),
    gt_alignment* const alignment_dst,gt_alignment* const alignment_src);

GT_INLINE gt_alignment* gt_alignment_union_alignment_maps_va(
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...);
#define gt_alignment_union_alignment_maps(alignment_src_A,alignment_src_B) \
        gt_alignment_union_alignment_maps_va(2,alignment_src_A,alignment_src_B)
GT_INLINE gt_alignment* gt_alignment_union_alignment_maps_fx_v(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,va_list v_args);
GT_INLINE gt_alignment* gt_alignment_union_alignment_maps_fx_va(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    const uint64_t num_src_alignments,gt_alignment* const alignment_src,...);
#define gt_alignment_union_alignment_maps_fx(gt_map_cmp_fx,alignment_src_A,alignment_src_B) \
        gt_alignment_union_alignment_maps_fx_va(gt_map_cmp_fx,2,alignment_src_A,alignment_src_B)

GT_INLINE gt_alignment* gt_alignment_subtract_alignment_maps_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_minuend,gt_alignment* const alignment_subtrahend);
GT_INLINE gt_alignment* gt_alignment_subtract_alignment_maps(
    gt_alignment* const alignment_minuend,gt_alignment* const alignment_subtrahend);

GT_INLINE gt_alignment* gt_alignment_intersect_alignment_maps_fx(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_alignment* const alignment_src_A,gt_alignment* const alignment_src_B);
GT_INLINE gt_alignment* gt_alignment_intersect_alignment_maps(
    gt_alignment* const alignment_src_A,gt_alignment* const alignment_src_B);

/*
 * Alignment realignment
 */
GT_INLINE void gt_alignment_recover_mismatches(gt_alignment* const alignment,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_alignment_realign_hamming(gt_alignment* const alignment,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_alignment_realign_levenshtein(gt_alignment* const alignment,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_alignment_realign_weighted(
    gt_alignment* const alignment,gt_sequence_archive* const sequence_archive,int32_t (*gt_weigh_fx)(char*,char*));

/*
 * Alignment trimming
 */
GT_INLINE void gt_alignment_hard_trim(gt_alignment* const alignment,const uint64_t left,const uint64_t right);
GT_INLINE void gt_alignment_hard_trim_min_length(gt_alignment* const alignment,const uint64_t left,const uint64_t right,const uint64_t min_length);
GT_INLINE void gt_alignment_quality_trim(gt_alignment* const alignment,const uint64_t quality_threshold,const uint64_t min_length);

GT_INLINE void gt_alignment_restore_trim(gt_alignment* const alignment);

#endif /* GT_ALIGNMENT_UTILS_H_ */
