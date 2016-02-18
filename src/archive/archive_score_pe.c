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
uint8_t archive_score_matches_pe_default_ties(matches_predictors_t* const predictors) {
  // Classify ties
  const double pr = matches_classify_logit_ties(predictors,&logit_model_paired_end_default);
  if (pr < PAIRED_MATCHES_MIN_CI) return 0;
  if (pr < PAIRED_MATCHES_TIES_CI) return 1;
  return archive_score_probability_scale(pr-PAIRED_MATCHES_TIES_CI,1.-PAIRED_MATCHES_TIES_CI,2,29);
}
uint8_t archive_score_matches_pe_default_mmap() {
//  // Classify multimaps
//  const double pr = matches_classify_logit_mmaps(predictors,&logit_model_paired_end_default);
//  if (pr < PAIRED_MATCHES_MIN_CI) return 0;
//  if (pr < PAIRED_MATCHES_MMAPS_CI) return 1;
//  return archive_score_probability_scale(pr-PAIRED_MATCHES_MMAPS_CI,1.-PAIRED_MATCHES_MMAPS_CI,30,49);
  return 30;
}
uint8_t archive_score_matches_pe_default_unique(matches_predictors_t* const predictors) {
  const double pr = matches_classify_logit_unique(predictors,&logit_model_paired_end_default);
  if (pr < PAIRED_MATCHES_MIN_CI) return 0;
  if (pr < PAIRED_MATCHES_UNIQUE_CI) return 1;
  return archive_score_probability_scale(pr-PAIRED_MATCHES_UNIQUE_CI,1.-PAIRED_MATCHES_UNIQUE_CI,31,49);
}
uint8_t archive_score_matches_pe_default_high_quality_ends(
    match_trace_t* const match_end1,
    match_trace_t* const match_end2) {
  const double pr = ((double)match_end1->mapq_score * (double)match_end2->mapq_score);
  return archive_score_probability_scale(pr,3600.0,50,60);
}
/*
 * PE Scoring Models
 */
void archive_score_matches_pe_default(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Parameters
  const uint64_t num_paired_map = paired_matches_get_num_maps(paired_matches);
  if (num_paired_map==0) return;
  // Score subdominant matches (MAPQ=0)
  paired_map_t* const paired_map = paired_matches_get_maps(paired_matches);
  uint64_t i;
  for (i=0;i<num_paired_map;++i) paired_map[i].mapq_score = 0;
  // Classify
  matches_predictors_t predictors;
  if (paired_map[0].pair_relation==pair_relation_discordant) return; // Discordant
  switch (paired_matches->paired_matches_class) {
    case paired_matches_class_unmapped:
    case paired_matches_class_subdominant_end:
    case paired_matches_class_tie_d0:
      break;
    case paired_matches_class_tie_d1:
      archive_search_pe_compute_predictors(archive_search_end1,archive_search_end2,paired_matches,&predictors);
      paired_map[0].mapq_score = archive_score_matches_pe_default_ties(&predictors);
      break;
    case paired_matches_class_mmap: // Lack of samples
      paired_map[0].mapq_score = archive_score_matches_pe_default_mmap();
      break;
    case paired_matches_class_unique:
      archive_search_pe_compute_predictors(archive_search_end1,archive_search_end2,paired_matches,&predictors);
      paired_map[0].mapq_score = archive_score_matches_pe_default_unique(&predictors);
      break;
    case paired_matches_class_high_quality_ends: {
      match_trace_t* const match_end1 = paired_map_get_match_end1(paired_matches,paired_map);
      match_trace_t* const match_end2 = paired_map_get_match_end2(paired_matches,paired_map);
      paired_map[0].mapq_score = archive_score_matches_pe_default_high_quality_ends(match_end1,match_end2);
      break;
    }
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
  if (num_paired_map>0 && paired_map[0].pair_relation==pair_relation_discordant) return;
  /*
   * Classify
   */
  matches_predictors_t predictors;
  archive_search_pe_compute_predictors(archive_search_end1,archive_search_end2,paired_matches,&predictors);
  if (paired_map[0].pair_relation==pair_relation_discordant) {
    paired_map[0].mapq_score = 1;
    return;
  }
  switch (paired_matches->paired_matches_class) {
    case paired_matches_class_unmapped:
      break;
    case paired_matches_class_subdominant_end:
      paired_map[0].mapq_score = 2;
      break;
    case paired_matches_class_tie_d0: {
      if (paired_map[0].template_length == paired_map[1].template_length) {
        paired_map[0].mapq_score = 3;
      } else {
        paired_map[0].mapq_score = 4;
      }
      break;
    }
    case paired_matches_class_tie_d1: {
//      double pr = matches_classify_logit_ties(&predictors,&logit_model_paired_end_default);
//      paired_map[0].mapq_score = (pr >= 0.95) ? 80 + (uint8_t)((pr-0.95)*1000.0) : 79;
      paired_map[0].mapq_score = 6;
      break;
    }
    case paired_matches_class_mmap: {
//      double pr = matches_classify_logit_mmaps(&predictors,&logit_model_paired_end_default);
//      paired_map[0].mapq_score = (pr >= 0.80) ? 140 + (uint8_t)((pr-0.80)*250.0) : 139;
      paired_map[0].mapq_score = 8;
      break;
    }
    case paired_matches_class_unique: {
//      double pr = matches_classify_logit_unique(&predictors,&logit_model_paired_end_default);
//      paired_map[0].mapq_score = (pr >= 0.95) ? 200 + (uint8_t)((pr-0.95)*1000.0) : 199;
      paired_map[0].mapq_score = 10;
      break;
    }
    case paired_matches_class_high_quality_ends:
      paired_map[0].mapq_score = 11;
      break;
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
  // Sort by distance (whichever it's selected)
  paired_matches_sort_by_distance(paired_matches);
  // Classify
  paired_matches->paired_matches_class = paired_matches_classify(paired_matches);
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
      matches_predictors_t predictors;
      archive_search_pe_compute_predictors(archive_search_end1,archive_search_end2,paired_matches,&predictors);
      // Dump predictors
      const char* const tag = sequence_get_tag(&archive_search_end1->sequence);
      matches_predictors_pe_print(stdout,tag,paired_matches->paired_matches_class,&predictors);
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROFILE_STOP(GP_ARCHIVE_SCORE_PE_MATCHES,PROFILE_LEVEL);
}
