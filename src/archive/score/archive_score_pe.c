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
 *   Archive MAPQ-scoring module produces MAPQ scores for PE-searches
 */

#include "archive/score/archive_score_pe.h"
#include "archive/score/archive_score_se.h"
#include "archive/search/archive_search_pe.h"
#include "matches/classify/matches_classify.h"
#include "matches/classify/matches_classify_logit.h"
#include "matches/classify/matches_classify_logit_models.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * PE Score Categories
 */
uint8_t archive_score_matches_pe_default_mmaps_d1(const matches_predictors_t* const matches_predictors) {
  const double pr = matches_classify_logit_mmaps_d1(
      &logit_model_paired_end_default,matches_predictors);
  if (pr < PAIRED_MATCHES_MIN_CI) return 5;
  if (pr < PAIRED_MATCHES_TIES_CI) return 7;
  return archive_score_probability_scale(pr-PAIRED_MATCHES_TIES_CI,1.-PAIRED_MATCHES_TIES_CI,10,28);
}
uint8_t archive_score_matches_pe_default_mmap(const matches_predictors_t* const matches_predictors) {
  const double pr = matches_classify_logit_mmaps(
      &logit_model_paired_end_default,matches_predictors);
  if (pr < PAIRED_MATCHES_MIN_CI) return 5;
  if (pr < PAIRED_MATCHES_MMAPS_CI) return 8;
  return archive_score_probability_scale(pr-PAIRED_MATCHES_MMAPS_CI,1.-PAIRED_MATCHES_MMAPS_CI,30,48);
}
uint8_t archive_score_matches_pe_default_unique(const matches_predictors_t* const matches_predictors) {
  const double pr = matches_classify_logit_unique(
      &logit_model_paired_end_default,matches_predictors);
  if (pr < PAIRED_MATCHES_MIN_CI) return 5;
  if (pr < PAIRED_MATCHES_UNIQUE_CI) return 9;
  return archive_score_probability_scale(pr-PAIRED_MATCHES_UNIQUE_CI,1.-PAIRED_MATCHES_UNIQUE_CI,50,59);
}
/*
 * PE Scoring Models
 */
void archive_score_matches_pe_default(paired_matches_t* const paired_matches) {
  // Unmapped
  if (paired_matches->paired_matches_class==paired_matches_class_unmapped) return;
  // Score subdominant matches (MAPQ=0)
  const uint64_t num_paired_map = paired_matches_get_num_maps(paired_matches);
  paired_map_t* const paired_map = paired_matches_get_maps(paired_matches);
  uint64_t i;
  for (i=0;i<num_paired_map;++i) paired_map[i].mapq_score = 0;
  // Classify
  match_trace_t* const primary_end1 = paired_map->match_trace_end1;
  match_trace_t* const primary_end2 = paired_map->match_trace_end2;
  // Discordant
  if (paired_map[0].pair_relation==pair_relation_discordant ||
      primary_end1->type == match_type_local || primary_end2->type == match_type_local) {
    paired_map[0].mapq_score = 2;
    return;
  }
  const bool high_quality_ends =
      primary_end1->mapq_score >= MAPQ_CONFIDENCE_SCORE_MIN &&
      primary_end2->mapq_score >= MAPQ_CONFIDENCE_SCORE_MIN;
  matches_predictors_t matches_predictors;
  switch (paired_matches->paired_matches_class) {
    case paired_matches_class_tie_perfect:
      if (num_paired_map==2) {
        paired_map[0].mapq_score = 1;
      }
      break;
    case paired_matches_class_tie:
      paired_map[0].mapq_score = (high_quality_ends) ? 4 : 3;
      break;
    case paired_matches_class_mmap_d1:
      matches_predictors_compute_pe(&matches_predictors,paired_matches);
      paired_map[0].mapq_score = (high_quality_ends) ? 29 : archive_score_matches_pe_default_mmaps_d1(&matches_predictors);
      break;
    case paired_matches_class_mmap:
      matches_predictors_compute_pe(&matches_predictors,paired_matches);
      paired_map[0].mapq_score = (high_quality_ends) ? 49 : archive_score_matches_pe_default_mmap(&matches_predictors);
      break;
    case paired_matches_class_unique:
      matches_predictors_compute_pe(&matches_predictors,paired_matches);
      paired_map[0].mapq_score = (high_quality_ends) ? 60 : archive_score_matches_pe_default_unique(&matches_predictors);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
void archive_score_matches_pe_stratify(paired_matches_t* const paired_matches) {
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
  match_trace_t* const primary_end1 = paired_map->match_trace_end1;
  match_trace_t* const primary_end2 = paired_map->match_trace_end2;
  const bool high_quality_ends =
      primary_end1->mapq_score >= MAPQ_CONFIDENCE_SCORE_MIN &&
      primary_end2->mapq_score >= MAPQ_CONFIDENCE_SCORE_MIN;
//  const bool quality_one = (primary_end1->mapq_score > 0 || primary_end2->mapq_score > 0);
//  const bool subdominant_both = (primary_end1->mapq_score==0 && primary_end2->mapq_score==0);
//  const bool subdominant_one = (primary_end1->mapq_score==0 || primary_end2->mapq_score==0);
//  const bool d0_limited_end1 = paired_matches->matches_end1->matches_class == matches_class_tie_perfect &&
//                               paired_matches->matches_end1->metrics.limited_candidates;
//  const bool d0_limited_end2 = paired_matches->matches_end2->matches_class == matches_class_tie_perfect &&
//                               paired_matches->matches_end2->metrics.limited_candidates;
  matches_predictors_t matches_predictors;
  matches_predictors_compute_pe(&matches_predictors,paired_matches);
  switch (paired_matches->paired_matches_class) {
    case paired_matches_class_tie_perfect: {
      paired_map[0].mapq_score = 1;
      break;
    }
    case paired_matches_class_tie: {
      paired_map[0].mapq_score = 2;
      break;
    }
    case paired_matches_class_mmap_d1: {
      if (high_quality_ends) {
        paired_map[0].mapq_score = 38;
      } else {
        const double pr = matches_classify_logit_mmaps_d1(&logit_model_paired_end_default,&matches_predictors);
        if (pr < 0.90) {
          paired_map[0].mapq_score = 39;
        } else {
          paired_map[0].mapq_score = archive_score_probability_scale(pr-0.90,1.0-0.90,40,90);
        }
      }
      break;
    }
    case paired_matches_class_mmap: {
      if (high_quality_ends) {
        paired_map[0].mapq_score = 98;
      } else {
        const double pr = matches_classify_logit_mmaps(&logit_model_paired_end_default,&matches_predictors);
        if (pr < 0.90) {
          paired_map[0].mapq_score = 99;
        } else {
          paired_map[0].mapq_score = archive_score_probability_scale(pr-0.90,1.0-0.90,100,150);
        }
      }
      break;
    }
    case paired_matches_class_unique: {
      if (high_quality_ends) {
        paired_map[0].mapq_score = 198;
      } else {
        const double pr = matches_classify_logit_unique(&logit_model_paired_end_default,&matches_predictors);
        if (pr < 0.90) {
          paired_map[0].mapq_score = 199;
        } else {
          paired_map[0].mapq_score = archive_score_probability_scale(pr-0.90,1.0-0.90,200,250);
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
  if (search_parameters->match_alignment_model==match_alignment_model_none) return;
  // Classify
  paired_matches_classify(paired_matches);
  /*
   * Select scoring model
   */
  PROFILE_START(GP_ARCHIVE_SCORE_PE_MATCHES,PROFILE_LEVEL);
  switch (search_parameters->mapq_model_pe) {
    case mapq_model_none: break;
    case mapq_model_classify:
      archive_score_matches_pe_stratify(paired_matches);
      break;
    case mapq_model_gem:
      archive_score_matches_pe_default(paired_matches);
      break;
    case mapq_model_dump_predictors: {
      // Score default
      archive_score_matches_pe_stratify(paired_matches);
      // Compute predictors
      matches_predictors_t matches_predictors;
      matches_predictors_compute_pe(&matches_predictors,paired_matches);
      // Dump predictors
      matches_predictors_pe_print(
          stdout,sequence_get_tag(&archive_search_end1->sequence),
          paired_matches->paired_matches_class,&matches_predictors);
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROFILE_STOP(GP_ARCHIVE_SCORE_PE_MATCHES,PROFILE_LEVEL);
}
