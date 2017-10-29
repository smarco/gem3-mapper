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
#include "archive/score/archive_score_logit.h"
#include "archive/score/archive_score_logit_models.h"
#include "matches/classify/matches_classify.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

///*
// * Archive Scoring Logit Classes
// */
//const archive_score_logit_model_t* archive_score_matches_pe_get_model(
//    search_parameters_t* const search_parameters) {
//  if (search_parameters->mapping_mode==mapping_adaptive_filtering_fast) {
//    return &logit_model_paired_end_fast;
//  } else {
//    return &logit_model_paired_end_sensitive;
//  }
//}
//uint8_t archive_score_matches_pe_logit_unique(
//    search_parameters_t* const search_parameters,
//    matches_predictors_t* const matches_predictors,
//    matches_classification_t* const matches_classification) {
//  // Fetch logit model
//  const archive_score_logit_model_t* const classify_model =
//      archive_score_matches_pe_get_model(search_parameters);
//  // Compute probability
//  const double pr = archive_score_logit_unique(
//      classify_model,matches_predictors,matches_classification);
//  if (pr >= 1.0) return 59;
//  const uint64_t mapq_score = archive_score_probability_to_mapq(pr);
//  return (uint8_t)MIN(mapq_score,59);
//}
//uint8_t archive_score_matches_pe_logit_mmap(
//    search_parameters_t* const search_parameters,
//    matches_predictors_t* const matches_predictors,
//    matches_classification_t* const matches_classification) {
//  // Fetch logit model
//  const archive_score_logit_model_t* const classify_model =
//      archive_score_matches_pe_get_model(search_parameters);
//  // Compute probability
//  const double pr = archive_score_logit_mmap(
//      classify_model,matches_predictors,matches_classification);
//  if (pr >= 1.0) return 59;
//  const uint64_t mapq_score = archive_score_probability_to_mapq(pr);
//  return (uint8_t)MIN(mapq_score,59);
//}
/*
 * PE Scoring Models
 */
void archive_score_matches_pe_default(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Parameters
  const uint64_t num_maps = paired_matches_get_num_maps(paired_matches);
  if (num_maps==0) return; // Unmmapped
  // Score subdominant matches (MAPQ=0)
  paired_map_t** const paired_maps = paired_matches_get_maps(paired_matches);
  uint64_t i;
  for (i=0;i<num_maps;++i) paired_maps[i]->mapq_score = 0;
  // Classify
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  paired_map_t* const primary_map = paired_maps[0];
  //  match_trace_t* const primary_end1 = primary_map->match_trace_end1;
  //  match_trace_t* const primary_end2 = primary_map->match_trace_end2;
  // Discordant
  if (primary_map->pair_relation==pair_relation_discordant) {
    primary_map->mapq_score = 0;
    return;
  }
  // Classify
  matches_predictors_t matches_predictors;
  matches_predictors_compute_pe(&matches_predictors,paired_matches);
  // Score
  const int64_t delta = paired_matches->classification.delta_group;
  const int64_t wdelta_group = paired_matches->classification.wdelta_group;
  if (matches_end1->classification.delta_group==0 ||
      matches_end2->classification.delta_group==0) {
    if (delta==0) {
      paired_map_t* const subdominant_map = paired_maps[1];
      primary_map->mapq_score = MIN(subdominant_map->template_length_sigma,4);
    } else {
      primary_map->mapq_score = 5 + ((wdelta_group <= 0) ? 0 : MIN(wdelta_group,10));
    }
  } else {
    if (delta == 0) { // Tie
      primary_map->mapq_score = 0;
    } else {
      if (wdelta_group > 1) {
        primary_map->mapq_score = 60;
      } else {
        primary_map->mapq_score =
            (primary_map->match_trace_end1->mapq_score +
             primary_map->match_trace_end2->mapq_score)/2;
      }
    }
  }
}
void archive_score_matches_pe_classify(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  // Parameters
  const uint64_t num_maps = paired_matches_get_num_maps(paired_matches);
  if (num_maps==0) return; // Unmmapped
  // Discordant
  paired_map_t* const primary_map = paired_matches_get_primary_map(paired_matches);
  if (primary_map->pair_relation==pair_relation_discordant) {
    primary_map->mapq_score = 0;
    return;
  }
  if (primary_map->match_trace_end1->type == match_type_local ||
      primary_map->match_trace_end2->type == match_type_local) {
    primary_map->mapq_score = 1;
    return;
  }
  // Classify
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  matches_predictors_t matches_predictors;
  matches_predictors_compute_pe(&matches_predictors,paired_matches);
  // Score
  paired_map_t** const paired_maps = paired_matches_get_maps(paired_matches);
  const int64_t delta = paired_matches->classification.delta_group;
  const int64_t wdelta_group = paired_matches->classification.wdelta_group;
  if (matches_end1->classification.delta_group==0 ||
      matches_end2->classification.delta_group==0) {
    if (delta==0) {
      paired_map_t* const subdominant_map = paired_maps[1];
      primary_map->mapq_score = 10 + MIN(subdominant_map->template_length_sigma,30);
    } else {
      primary_map->mapq_score = 50 + ((wdelta_group <= 0) ? 0 : MIN(wdelta_group,30));
    }
  } else {
    if (delta == 0) { // Tie
      primary_map->mapq_score = 2;
    } else {
      if (wdelta_group > 1) {
        primary_map->mapq_score = 161;
      } else {
        primary_map->mapq_score = 100 +
            (primary_map->match_trace_end1->mapq_score +
             primary_map->match_trace_end2->mapq_score)/2;
      }
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
  search_parameters_t* const search_parameters = &archive_search_end1->search_parameters;
  if (search_parameters->match_alignment_model==match_alignment_model_none) return;
  /*
   * Select scoring model
   */
  PROFILE_START(GP_ARCHIVE_SCORE_PE_MATCHES,PROFILE_LEVEL);
  switch (search_parameters->mapq_model_pe) {
    case mapq_model_none: break;
    case mapq_model_classify:
      archive_score_matches_pe_classify(archive_search_end1,archive_search_end2,paired_matches);
      break;
    case mapq_model_gem:
      archive_score_matches_pe_default(archive_search_end1,archive_search_end2,paired_matches);
      break;
    case mapq_model_dump_predictors: {
      // Score default
      archive_score_matches_pe_classify(archive_search_end1,archive_search_end2,paired_matches);
      // Compute predictors
      matches_predictors_t matches_predictors;
      matches_predictors_compute_pe(&matches_predictors,paired_matches);
      // Dump predictors
      char* const sequence_tag = sequence_get_tag(archive_search_end1->sequence);
      matches_predictors_pe_print(stdout,sequence_tag,&matches_predictors,&paired_matches->classification);
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  PROFILE_STOP(GP_ARCHIVE_SCORE_PE_MATCHES,PROFILE_LEVEL);
}
