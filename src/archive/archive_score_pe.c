/*
 * PROJECT: GEMMapper
 * FILE: archive_score_pe.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive/archive_score_pe.h"
#include "archive/archive_score_se.h"
#include "archive/archive_search_pe.h"
#include "matches/matches_classify.h"
#include "matches/matches_classify_logit.h"
#include "matches/matches_classify_logit_models.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * PE Score Categories
 */
uint8_t archive_score_matches_pe_default_ties(
    const matches_predictors_t* const matches_predictors,
    const match_predictors_t* const match_predictors) {
  // Classify ties
  const double pr = matches_classify_logit_ties(
      &logit_model_paired_end_default,matches_predictors,match_predictors);
  if (pr < PAIRED_MATCHES_MIN_CI) return 0;
  if (pr < PAIRED_MATCHES_TIES_CI) return 1;
  return archive_score_probability_scale(pr-PAIRED_MATCHES_TIES_CI,1.-PAIRED_MATCHES_TIES_CI,3,28);
}
uint8_t archive_score_matches_pe_default_mmap(
    const matches_predictors_t* const matches_predictors,
    const match_predictors_t* const match_predictors) {
  // Classify multimaps
  const double pr = matches_classify_logit_mmaps(
      &logit_model_paired_end_default,matches_predictors,match_predictors);
  if (pr < PAIRED_MATCHES_MIN_CI) return 0;
  if (pr < PAIRED_MATCHES_MMAPS_CI) return 2;
  return archive_score_probability_scale(pr-PAIRED_MATCHES_MMAPS_CI,1.-PAIRED_MATCHES_MMAPS_CI,30,48);
}
uint8_t archive_score_matches_pe_default_unique(
    const matches_predictors_t* const matches_predictors,
    const match_predictors_t* const match_predictors) {
  const double pr = matches_classify_logit_unique(
      &logit_model_paired_end_default,matches_predictors,match_predictors);
  if (pr < PAIRED_MATCHES_MIN_CI) return 0;
  if (pr < PAIRED_MATCHES_UNIQUE_CI) return 2;
  return archive_score_probability_scale(pr-PAIRED_MATCHES_UNIQUE_CI,1.-PAIRED_MATCHES_UNIQUE_CI,50,59);
}
/*
 * PE Scoring Models
 */
void archive_score_matches_pe_default(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Unmapped
  if (paired_matches->paired_matches_class==paired_matches_class_unmapped) return;
  // Score subdominant matches (MAPQ=0)
  const uint64_t num_paired_map = paired_matches_get_num_maps(paired_matches);
  paired_map_t* const paired_map = paired_matches_get_maps(paired_matches);
  uint64_t i;
  for (i=0;i<num_paired_map;++i) paired_map[i].mapq_score = 0;
  // Discordant
  if (paired_map[0].pair_relation==pair_relation_discordant) return;
  // Classify
  match_trace_t* const primary_end1 = paired_map_get_match_end1(paired_matches,paired_map);
  match_trace_t* const primary_end2 = paired_map_get_match_end2(paired_matches,paired_map);
  const bool high_quality_ends = (primary_end1->mapq_score >= 30 && primary_end2->mapq_score >= 30);
  matches_predictors_t matches_predictors;
  match_predictors_t match_predictors;
  switch (paired_matches->paired_matches_class) {
    case paired_matches_class_tie_d0:
      break;
    case paired_matches_class_tie_d1:
      matches_predictors_compute_pe(&matches_predictors,paired_matches);
      match_predictors_compute_pe(&match_predictors,paired_matches,paired_map);
      if (high_quality_ends) {
        paired_map[0].mapq_score = 29;
      } else {
        paired_map[0].mapq_score = archive_score_matches_pe_default_ties(&matches_predictors,&match_predictors);
      }
      break;
    case paired_matches_class_mmap:
      matches_predictors_compute_pe(&matches_predictors,paired_matches);
      match_predictors_compute_pe(&match_predictors,paired_matches,paired_map);
      if (high_quality_ends) {
        paired_map[0].mapq_score = 49;
      } else {
        paired_map[0].mapq_score = archive_score_matches_pe_default_mmap(&matches_predictors,&match_predictors);
      }
      break;
    case paired_matches_class_unique:
      matches_predictors_compute_pe(&matches_predictors,paired_matches);
      match_predictors_compute_pe(&match_predictors,paired_matches,paired_map);
      if (high_quality_ends) {
        paired_map[0].mapq_score = 60;
      } else {
        paired_map[0].mapq_score = archive_score_matches_pe_default_unique(&matches_predictors,&match_predictors);
      }
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
void archive_score_matches_pe_stratify(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Score subdominant matches (MAPQ=0)
  paired_map_t* const paired_map = paired_matches_get_maps(paired_matches);
  const uint64_t num_paired_map = paired_matches_get_num_maps(paired_matches);
  uint64_t i;
  for (i=0;i<num_paired_map;++i) paired_map[i].mapq_score = 0;
  // Unmmapped
  if (paired_matches->paired_matches_class==paired_matches_class_unmapped) return;
  // Discordant
  if (paired_map[0].pair_relation==pair_relation_discordant) {
    paired_map[0].mapq_score = 1;
    return;
  }
  // Classify
  match_trace_t* const primary_end1 = paired_map_get_match_end1(paired_matches,paired_map);
  match_trace_t* const primary_end2 = paired_map_get_match_end2(paired_matches,paired_map);
  const bool high_quality_ends = (primary_end1->mapq_score >= 30 && primary_end2->mapq_score >= 30);
  matches_predictors_t matches_predictors;
  match_predictors_t match_predictors;
  matches_predictors_compute_pe(&matches_predictors,paired_matches);
  match_predictors_compute_pe(&match_predictors,paired_matches,paired_map);
  switch (paired_matches->paired_matches_class) {
    case paired_matches_class_tie_d0: {
      if (paired_matches->matches_end1->matches_class == matches_class_tie_d0 &&
          paired_matches->matches_end2->matches_class == matches_class_tie_d0) {
        paired_map[0].mapq_score = 3;
      } else if (paired_matches->matches_end1->matches_class == matches_class_tie_d0 ||
                 paired_matches->matches_end2->matches_class == matches_class_tie_d0) {
        paired_map[0].mapq_score = 4;
      } else if (paired_map[0].template_length == paired_map[1].template_length) {
        paired_map[0].mapq_score = 5;
      } else {
        paired_map[0].mapq_score = 6;
      }
      break;
    }
    case paired_matches_class_tie_d1: {
      if (high_quality_ends) {
        paired_map[0].mapq_score = 8;
      } else {
        const double pr = matches_classify_logit_ties(
            &logit_model_paired_end_default,&matches_predictors,&match_predictors);
        if (pr < 0.90) {
          paired_map[0].mapq_score = 9;
        } else {
          paired_map[0].mapq_score = archive_score_probability_scale(pr-0.90,1.0-0.90,10,50);
        }
      }
      break;
    }
    case paired_matches_class_mmap: {
      if (high_quality_ends) {
        paired_map[0].mapq_score = 58;
      } else {
        const double pr = matches_classify_logit_mmaps(
            &logit_model_paired_end_default,&matches_predictors,&match_predictors);
        if (pr < 0.90) {
          paired_map[0].mapq_score = 59;
        } else {
          paired_map[0].mapq_score = archive_score_probability_scale(pr-0.90,1.0-0.90,60,100);
        }
      }
      break;
    }
    case paired_matches_class_unique: {
      if (high_quality_ends) {
        paired_map[0].mapq_score = 108;
      } else {
        const double pr = matches_classify_logit_unique(
            &logit_model_paired_end_default,&matches_predictors,&match_predictors);
        if (pr < 0.90) {
          paired_map[0].mapq_score = 109;
        } else {
          paired_map[0].mapq_score = archive_score_probability_scale(pr-0.90,1.0-0.90,110,150);
        }
      }
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
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
  search_parameters_t* const search_parameters = &archive_search_end1->search_parameters;
  if (search_parameters->alignment_model==alignment_model_none) return;
  // Classify
  paired_matches_classify(paired_matches);
  /*
   * Select scoring model
   */
  PROFILE_START(GP_ARCHIVE_SCORE_PE_MATCHES,PROFILE_LEVEL);
  switch (search_parameters->mapq_model_pe) {
    case mapq_model_none: break;
    case mapq_model_classify:
      archive_score_matches_pe_stratify(archive_search_end1,archive_search_end2,paired_matches);
      break;
    case mapq_model_gem:
      archive_score_matches_pe_default(archive_search_end1,archive_search_end2,paired_matches);
      break;
    case mapq_model_dump_predictors: {
      // Score default
      archive_score_matches_pe_stratify(archive_search_end1,archive_search_end2,paired_matches);
      // Compute predictors
      matches_predictors_t matches_predictors;
      match_predictors_t match_predictors;
      paired_map_t* const paired_map = paired_matches_is_mapped(paired_matches) ?
          paired_matches_get_maps(paired_matches) : NULL;
      matches_predictors_compute_pe(&matches_predictors,paired_matches);
      match_predictors_compute_pe(&match_predictors,paired_matches,paired_map);
      // Dump predictors
      matches_predictors_pe_print(
          stdout,sequence_get_tag(&archive_search_end1->sequence),
          paired_matches->paired_matches_class,&matches_predictors,&match_predictors);
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROFILE_STOP(GP_ARCHIVE_SCORE_PE_MATCHES,PROFILE_LEVEL);
}
