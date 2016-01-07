/*
 * PROJECT: GEMMapper
 * FILE: archive_score_pe.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_score_pe.h"
#include "archive_score_se.h"
#include "archive_search_pe.h"
#include "matches_classify_logit.h"
#include "matches_classify_logit_models.h"
#include "paired_matches_classify.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * PE Score Categories
 */
uint8_t archive_score_matches_pe_default_ties(matches_predictors_t* const predictors) {
  // Classify ties
  const double pr = matches_classify_logit_ties(predictors,&logit_model_paired_end_default);
  if (pr < PAIRED_MATCHES_MIN_CI) return 0;
  if (pr < PAIRED_MATCHES_TIES_CI) return 1;
  return archive_score_probability_scale(pr-PAIRED_MATCHES_TIES_CI,1.-PAIRED_MATCHES_TIES_CI,2,29);
}
uint8_t archive_score_matches_pe_default_mmap(matches_predictors_t* const predictors) {
  // Classify multimaps
  const double pr = matches_classify_logit_mmaps(predictors,&logit_model_paired_end_default);
  if (pr < PAIRED_MATCHES_MIN_CI) return 0;
  if (pr < PAIRED_MATCHES_MMAPS_CI) return 1;
  return archive_score_probability_scale(pr-PAIRED_MATCHES_MMAPS_CI,1.-PAIRED_MATCHES_MMAPS_CI,30,49);
}
uint8_t archive_score_matches_pe_default_unique(matches_predictors_t* const predictors) {
  const double pr = matches_classify_logit_unique(predictors,&logit_model_paired_end_default);
  if (pr < PAIRED_MATCHES_MIN_CI) return 0;
  if (pr < PAIRED_MATCHES_UNIQUE_CI) return 1;
  return archive_score_probability_scale(pr-PAIRED_MATCHES_UNIQUE_CI,1.-PAIRED_MATCHES_UNIQUE_CI,50,60);
}
/*
 * PE Scoring Models
 */
void archive_score_matches_pe_default(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Parameters
  paired_map_t* const paired_map = paired_matches_get_maps(paired_matches);
  const uint64_t num_paired_map = paired_matches_get_num_maps(paired_matches);
  // Score subdominant matches (MAPQ=0)
  uint64_t i;
  for (i=0;i<num_paired_map;++i) paired_map[i].mapq_score = 0;
  if (num_paired_map>0 && paired_map[0].pair_relation==pair_relation_discordant) return;
  // Classify
  matches_predictors_t predictors;
  archive_search_pe_compute_predictors(archive_search_end1,archive_search_end2,paired_matches,&predictors);
  const matches_class_t matches_class = paired_matches_classify(paired_matches);
  switch (matches_class) {
    case matches_class_unmapped:
    case matches_class_tie_indistinguishable:
    case matches_class_tie_swg_score:
    case matches_class_tie_edit_distance:
      break;
    case matches_class_tie_event_distance:
      paired_map[0].mapq_score = archive_score_matches_pe_default_ties(&predictors);
      break;
    case matches_class_mmap:
      paired_map[0].mapq_score = archive_score_matches_pe_default_mmap(&predictors);
      break;
    case matches_class_unique:
      paired_map[0].mapq_score = archive_score_matches_pe_default_unique(&predictors);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
void archive_score_matches_pe_stratify(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Score subdominant matches (MAPQ=0)
  paired_map_t* const paired_map = paired_matches_get_maps(paired_matches);
  const uint64_t num_paired_map = paired_matches_get_num_maps(paired_matches);
  double pr;
  uint64_t i;
  for (i=0;i<num_paired_map;++i) paired_map[i].mapq_score = 0;
  if (num_paired_map>0 && paired_map[0].pair_relation==pair_relation_discordant) return;
  /*
   * Classify
   *   200-250   Unique (high likelihood)
   *   199       Unique (low likelihood)
   *
   *   140-190   MMaps (high likelihood)
   *   139       MMaps (low likelihood)
   *
   *   80-130    Tie (high likelihood)
   *   79        Tie (low likelihood)
   *
   *   4         Tie Edit
   *   3         Tie SWG
   *   2         Ambiguous (not enough resolution)
   *   1         Tie Indistinguishable (low likelihood)
   *   0         Subdominant
   */
  const matches_class_t matches_class = paired_matches_classify(paired_matches);
  const uint64_t mcs_end1 = paired_matches->matches_end1->max_complete_stratum;
  const uint64_t mcs_end2 = paired_matches->matches_end2->max_complete_stratum;
  matches_predictors_t predictors;
  archive_search_pe_compute_predictors(archive_search_end1,archive_search_end2,paired_matches,&predictors);
  if (matches_class==matches_class_unmapped) {
    return;
  } else if (matches_class==matches_class_tie_indistinguishable) {
    paired_map[0].mapq_score = 1;
  } else if ((archive_search_end1->end_class==matches_class_unmapped || mcs_end1<=1) && mcs_end2<=1) {
    paired_map[0].mapq_score = 2;
  } else {
    switch (matches_class) {
      case matches_class_tie_swg_score:
        paired_map[0].mapq_score = 3;
        break;
      case matches_class_tie_edit_distance:
        paired_map[0].mapq_score = 4;
        break;
      case matches_class_tie_event_distance:
        pr = matches_classify_logit_ties(&predictors,&logit_model_paired_end_default);
        paired_map[0].mapq_score = (pr >= 0.80) ? 80 + (uint8_t)((pr-0.80)*250.0) : 79;
        break;
      case matches_class_mmap:
        pr = matches_classify_logit_mmaps(&predictors,&logit_model_paired_end_default);
        paired_map[0].mapq_score = (pr >= 0.80) ? 140 + (uint8_t)((pr-0.80)*250.0) : 139;
        break;
      case matches_class_unique:
        pr = matches_classify_logit_unique(&predictors,&logit_model_paired_end_default);
        paired_map[0].mapq_score = (pr >= 0.80) ? 200 + (uint8_t)((pr-0.80)*250.0) : 199;
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
/*
 * Archive Scoring PE
 */
void archive_score_matches_pe(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Check no-matches
  if (paired_matches_get_num_maps(paired_matches)==0) return;
  // Check no-align model
  search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
  if (search_parameters->alignment_model==alignment_model_none) return;
  // Sort by distance (whichever it's selected)
  paired_matches_sort_by_distance(paired_matches);
  /*
   * Select scoring model
   */
  PROFILE_START(GP_ARCHIVE_SCORE_PE_MATCHES,PROFILE_LEVEL);
  switch (search_parameters->mapq_model) {
    case mapq_model_none: break;
    case mapq_model_classify:
      archive_score_matches_pe_stratify(archive_search_end1,archive_search_end2,paired_matches);
      break;
    case mapq_model_gem:
      archive_score_matches_pe_default(archive_search_end1,archive_search_end2,paired_matches);
      break;
    case mapq_model_dump_predictors: {
      // Score default
      archive_score_matches_pe_default(archive_search_end1,archive_search_end2,paired_matches);
      // Compute predictors
      paired_map_t* const paired_map = paired_matches_get_maps(paired_matches);
      match_trace_t* const match_end1 = paired_map_get_match_end1(paired_matches,paired_map);
      match_trace_t* const match_end2 = paired_map_get_match_end2(paired_matches,paired_map);
      matches_predictors_t predictors;
      archive_search_pe_compute_predictors(archive_search_end1,archive_search_end2,paired_matches,&predictors);
      // Dump predictors
      const uint8_t mapq_score_end1 = (archive_search_end1->pair_searched) ? match_end1->mapq_score : 0;
      const uint8_t mapq_score_end2 = (archive_search_end2->pair_searched) ? match_end2->mapq_score : 0;
      paired_matches_predictors_print(
          stdout,&predictors,sequence_get_tag(&archive_search_end1->sequence),
          paired_map[0].mapq_score,mapq_score_end1,mapq_score_end2);
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROFILE_STOP(GP_ARCHIVE_SCORE_PE_MATCHES,PROFILE_LEVEL);
}
