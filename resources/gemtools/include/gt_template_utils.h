/*
 * PROJECT: GEM-Tools library
 * FILE: gt_template_utils.h
 * DATE: 19/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_TEMPLATE_UTILS_H_
#define GT_TEMPLATE_UTILS_H_

#include "gt_essentials.h"
#include "gt_alignment_utils.h"
#include "gt_template.h"
#include "gt_map_score.h"

/*
 * Template basic tools
 */
GT_INLINE void gt_template_setup_pair_attributes_to_alignments(gt_template* const template,const bool copy_tags);
GT_INLINE uint64_t gt_template_get_read_proportion(gt_template* const template,const float proportion);

/*
 * Template's MMaps high-level insertion (basic building block)
 */
GT_INLINE gt_map** gt_template_raw_put_mmap(
    int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),gt_template* const template,
    gt_map** const mmap,gt_mmap_attributes* const mmap_attributes);
GT_INLINE gt_map** gt_template_put_mmap(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attr,
    const bool check_duplicated,const bool replace_duplicated);

/*
 * Template's MMaps high-level insertion operators (Update global state: counters, ...)
 */
GT_INLINE void gt_template_insert_mmap(
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attributes);
GT_INLINE void gt_template_insert_mmap_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap,gt_mmap_attributes* const mmap_attributes);
GT_INLINE void gt_template_insert_mmap_gtvector(
    gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attributes);
GT_INLINE void gt_template_insert_mmap_gtvector_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_vector* const mmap,gt_mmap_attributes* const mmap_attributes);

GT_INLINE bool gt_template_find_mmap_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap,
    uint64_t* const found_mmap_pos,gt_map*** const found_mmap,gt_mmap_attributes** const found_mmap_attributes);
GT_INLINE bool gt_template_is_mmap_contained(gt_template* const template,gt_map** const mmap);
GT_INLINE bool gt_template_is_mmap_contained_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),
    gt_template* const template,gt_map** const mmap);

GT_INLINE void gt_template_reduce_mmaps(gt_template* const template,const uint64_t max_num_matches);

/*
 * Template's Insert Size
 */
#define GT_TEMPLATE_INSERT_SIZE_OK 0
#define GT_TEMPLATE_INSERT_SIZE_DIFFERENT_CONTIGS 1
#define GT_TEMPLATE_INSERT_SIZE_SAME_STRAND 2
GT_INLINE int64_t gt_template_get_insert_size(gt_map** const mmap,gt_status* gt_error,uint64_t *start_x,gt_string **ctg);

/*
 * Template's Counters operators
 */
GT_INLINE bool gt_template_is_mapped(gt_template* const template);
GT_INLINE bool gt_template_is_thresholded_mapped(gt_template* const template,const uint64_t max_allowed_strata);
GT_INLINE void gt_template_recalculate_counters(gt_template* const template);
GT_INLINE void gt_template_recalculate_counters_no_splits(gt_template* const template);

GT_INLINE int64_t gt_template_get_min_matching_strata(gt_template* const template);
GT_INLINE int64_t gt_template_get_uniq_degree(gt_template* const template);
GT_INLINE bool gt_template_get_next_matching_strata(
    gt_template* const template,const uint64_t begin_strata,
    uint64_t* const next_matching_strata,uint64_t* const num_maps);

/*
 * Template's Maps Sorting
 */
GT_INLINE void gt_template_sort_by_distance__score(gt_template* const template);
GT_INLINE void gt_template_sort_by_distance__score_no_split(gt_template* const template);

/*
 * Template's MMaps Utils
 */
GT_INLINE uint64_t gt_template_sum_mismatch_qualities(gt_template* const template,gt_map** const mmap);
GT_INLINE uint64_t gt_template_get_max_mismatch_quality(gt_template* const template);

/*
 * Template Set operators
 */
GT_INLINE void gt_template_merge_template_mmaps(gt_template* const template_dst,gt_template* const template_src);
GT_INLINE void gt_template_merge_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template_dst,gt_template* const template_src);

GT_INLINE gt_template* gt_template_union_template_mmaps_v(
    const uint64_t num_src_templates,gt_template* const template_src,va_list v_args);
GT_INLINE gt_template* gt_template_union_template_mmaps_va(
    const uint64_t num_src_templates,gt_template* const template_src,...);
GT_INLINE gt_template* gt_template_union_template_mmaps_a(
    gt_template** const templates,const uint64_t num_src_templates);
#define gt_template_union_template_mmaps(template_src_A,template_src_B) \
        gt_template_union_template_mmaps_va(2,template_src_A,template_src_B)

GT_INLINE gt_template* gt_template_union_template_mmaps_fx_v(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    const uint64_t num_src_templates,gt_template* const template_src,va_list v_args);
GT_INLINE gt_template* gt_template_union_template_mmaps_fx_va(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    const uint64_t num_src_templates,gt_template* const template_src,...);
#define gt_template_union_template_mmaps_fx(gt_mmap_cmp_fx,gt_map_cmp_fx,template_src_A,template_src_B) \
        gt_template_union_template_mmaps_fx_va(gt_mmap_cmp_fx,gt_map_cmp_fx,2,template_src_A,template_src_B)

GT_INLINE gt_template* gt_template_subtract_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template_minuend,gt_template* const template_subtrahend);
GT_INLINE gt_template* gt_template_subtract_template_mmaps(
    gt_template* const template_minuend,gt_template* const template_subtrahend);

GT_INLINE gt_template* gt_template_intersect_template_mmaps_fx(
    int64_t (*gt_mmap_cmp_fx)(gt_map**,gt_map**,uint64_t),int64_t (*gt_map_cmp_fx)(gt_map*,gt_map*),
    gt_template* const template_A,gt_template* const template_B);
GT_INLINE gt_template* gt_template_intersect_template_mmaps(
    gt_template* const template_A,gt_template* const template_B);

/*
 * Template realignment
 */
GT_INLINE void gt_template_recover_mismatches(gt_template* const template,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_template_realign_hamming(gt_template* const template,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_template_realign_levenshtein(gt_template* const template,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_template_realign_weighted(
    gt_template* const template,gt_sequence_archive* const sequence_archive,int32_t (*gt_weigh_fx)(char*,char*));

/*
 * Template trimming
 */
GT_INLINE void gt_template_hard_trim(gt_template* const template,const uint64_t left,const uint64_t right);
GT_INLINE void gt_template_restore_trim(gt_template* const template);

/*
 * Template/Alignment Placeholder (mmap/map)
 */
typedef enum { GT_MMAP_PLACEHOLDER_PAIRED=0, GT_MMAP_PLACEHOLDER_UNPAIRED=1, GT_MAP_PLACEHOLDER=2 } gt_map_placeholder_t;
typedef struct {
  /* Source Template/Alignment (SE) */
  gt_template* template;
  gt_alignment* alignment;
} gt_map_placeholder_se;
typedef struct {
  /* Source Template (PE) */
  gt_template* template;
  uint64_t paired_end_position;
  /* Mate */
  gt_map* mate;
  gt_mmap_attributes* mmap_attributes;
} gt_map_placeholder_pe;
typedef struct {
  gt_map_placeholder_t type;
  /* Template/Alignment info (SE/PE) */
  gt_map* map; // Current Map
  union {
    gt_map_placeholder_se single_end; /* For {GT_MAP_PLACEHOLDER} */
    gt_map_placeholder_pe paired_end; /* For {GT_MMAP_PLACEHOLDER_PAIRED, GT_MMAP_PLACEHOLDER_UNPAIRED} */
  };
  /* SAM info */
  bool secondary_alignment;
  bool not_passing_QC;
  bool PCR_duplicate;
  uint32_t hard_trim_left;
  uint32_t hard_trim_right;
} gt_map_placeholder;

GT_INLINE gt_map_placeholder* gt_map_placeholder_new();
GT_INLINE void gt_map_placeholder_delete(gt_map_placeholder* const map_placeholder);
GT_INLINE void gt_map_placeholder_clear(gt_map_placeholder* const map_placeholder);
GT_INLINE void gt_map_placeholder_set_sam_fields(gt_map_placeholder* const map_placeholder,
    const bool not_passing_QC,const bool PCR_duplicate,const uint32_t hard_trim_left,const uint32_t hard_trim_right);

/*
 * Template/Alignment Placeholder Compare functions
 */
int gt_map_placeholder_cmp_coordinates(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b);

int gt_map_placeholder_cmp_map_distance(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b);
int gt_map_placeholder_cmp_map_phred_scores(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b);
int gt_map_placeholder_cmp_map_gt_scores_ascending(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b);
int gt_map_placeholder_cmp_map_gt_scores_descending(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b);

int gt_map_placeholder_cmp_mmap_phred_scores(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b);
int gt_map_placeholder_cmp_mmap_gt_scores_ascending(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b);
int gt_map_placeholder_cmp_mmap_gt_scores_descending(gt_map_placeholder* const ph_map_a,gt_map_placeholder* const ph_map_b);

/*
 * Fills the placeholder's vector @mmap_placeholder with the map/mmap
 */
GT_INLINE void gt_map_placeholder_add_map(
    gt_map* const map,gt_string* const read,
    gt_vector* const mmap_placeholder,const bool split_segments,
    int (*gt_ph_cmp_fx)(gt_map_placeholder* const,gt_map_placeholder* const),const bool cmp_with_best,
    gt_map_placeholder* const best_mmap_ph,uint64_t* const best_mmap_ph_position,
    gt_map_placeholder* const mmap_ph);
GT_INLINE void gt_map_placeholder_add_mmap(
    gt_map* const map_endA,gt_map* const map_endB,gt_string* const read_endA,const uint64_t paired_end_position,
    gt_vector* const mmap_placeholder,const bool split_segments,
    int (*gt_ph_cmp_fx)(gt_map_placeholder* const,gt_map_placeholder* const),const bool cmp_with_best,
    gt_map_placeholder* const best_mmap_ph,uint64_t* const best_mmap_ph_position,
    gt_map_placeholder* const mmap_ph);

/*
 * Fills the placeholder's vector @mmap_placeholder with all the mmaps from template
 *   - If @include_mate_placeholder is set, inserts a placeholder for the mate
 *       ({ph.paired_end.paired_end_position==1}). i.e. As for SAM output like operations
 *   - If @split_segments is set, individual segments are inserted into @mmap_placeholder (taking into account quimeras/segments)
 *   - If @primary_mmap_end1_pos is not null, sets it to the position in the array of the primary/best mmap(end/1) placeholder
 *       Sorts w.r.t @gt_ph_cmp_fx sorting. If @gt_ph_cmp_fx is null, then is set to the first one added.
 *   - If @primary_mmap_end2_pos is not null, sets it to the position in the array of the primary/best mmap(end/2) placeholder
 *       Sorts w.r.t @gt_ph_cmp_fx sorting. If @gt_ph_cmp_fx is null, then is set to the first one added.
 *       If @include_mate_placeholder is not set, then @primary_mmap_end1_pos==@primary_mmap_end2_pos
 *   - If @placeholder_template is not null, then its values are used as defaults (not_passing_QC,PCR_duplicate,...)
 */
GT_INLINE void gt_map_placeholder_build_from_template(
    gt_template* const template,gt_vector* const mmap_placeholder,
    const bool include_mate_placeholder,const bool split_segments,const uint64_t max_num_maps,
    int (*gt_ph_cmp_fx)(gt_map_placeholder* const,gt_map_placeholder* const),
    uint64_t* const primary_mmap_end1_pos,uint64_t* const primary_mmap_end2_pos,
    gt_map_placeholder* const placeholder_template);
/*
 * Fills the placeholder's vector @mmap_placeholder with all the maps from alignment
 *   - If @split_segments is set, individual segments are inserted into @mmap_placeholder (taking into account quimeras/segments)
 *   - If @primary_map_position is not null, sets it to the position in the array of the primary/best alignment placeholder
 *       Sorts w.r.t @gt_ph_cmp_fx sorting. If @gt_ph_cmp_fx is null, then is set to the first one added.
 *       If @include_mate_placeholder is not set, then @primary_mmap_end1_pos==@primary_mmap_end2_pos
 *   - If @placeholder_template is not null, then its values are used as defaults (not_passing_QC,PCR_duplicate,...)
 */
GT_INLINE void gt_map_placeholder_build_from_alignment(
    gt_alignment* const alignment,gt_vector* const mmap_placeholder,const bool split_segments,const uint64_t max_num_maps,
    int (*gt_ph_cmp_fx)(gt_map_placeholder* const,gt_map_placeholder* const),uint64_t* const primary_map_position,
    gt_map_placeholder* const placeholder_template);

#endif /* GT_TEMPLATE_UTILS_H_ */
