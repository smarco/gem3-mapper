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

#include "archive/score/archive_score_logit_models.h"

/*
 * SE-Model DEFAULT
 */
const archive_score_logit_model_t logit_model_single_end_fast = {
//    Coefficients:
//    (Intercept)         edit        event          swg        equal        acand
//     -222.88906     30.22311    169.06895     16.16693     -0.02592      1.09730
//           mcs1        arlen        mrlen         kmer
//        6.32597      1.18261     -1.64119      1.06775
  .unique_logit_coeff = {
      .coeff_intercept                            = -222.88906,
      .coeff_primary_edit_distance_norm           = 30.22311,
      .coeff_primary_event_distance_norm          = 169.06895,
      .coeff_primary_swg_score_norm               = 16.16693,
      .coeff_primary_error_quality                = -0.02592,
      .coeff_primary_template_size_sigma_norm     = 0.0,
      .coeff_mapq_end1                            = 0.0,
      .coeff_mapq_end2                            = 0.0,
      .coeff_subdominant_edit_distance_norm       = 0.0,
      .coeff_subdominant_event_distance_norm      = 0.0,
      .coeff_subdominant_swg_score_norm           = 0.0,
      .coeff_subdominant_template_size_sigma_norm = 0.0,
      .coeff_candidates_subdominant               = 0.0,
      .coeff_candidates_accepted                  = 1.09730,
      .coeff_mcs_end1                             = 6.32597,
      .coeff_mcs_end2                             = 0.0,
      .coeff_delta_group                          = 0.0,
      .coeff_wdelta_group                         = 0.0,
      .coeff_sequence_length_norm                 = 0.0,
      .coeff_avg_region_length_norm               = 1.18261,
      .coeff_max_region_length_norm               = -1.64119,
      .coeff_kmer_frequency                       = 1.06775,
  },
//  Coefficients:
//  (Intercept)       wdelta         edit        event          swg        equal
//    -83.10499      1.25320    -69.58770    183.12527     13.63890     -0.02174
//     sub_edit      sub_swg        scand        acand         mcs1        arlen
//    -37.64175     -9.02421     -0.18200      0.06244      1.77801      0.54612
//         kmer
//      0.12584
  .mmap_logit_coeff = {
      .coeff_intercept                            = -83.10499,
      .coeff_primary_edit_distance_norm           = -69.58770,
      .coeff_primary_event_distance_norm          = 183.12527,
      .coeff_primary_swg_score_norm               = 13.63890,
      .coeff_primary_error_quality                = -0.02174,
      .coeff_primary_template_size_sigma_norm     = 0.0,
      .coeff_mapq_end1                            = 0.0,
      .coeff_mapq_end2                            = 0.0,
      .coeff_subdominant_edit_distance_norm       = -37.64175,
      .coeff_subdominant_event_distance_norm      = 0.0,
      .coeff_subdominant_swg_score_norm           = -9.02421,
      .coeff_subdominant_template_size_sigma_norm = 0.0,
      .coeff_candidates_subdominant               = -0.18200,
      .coeff_candidates_accepted                  = 0.06244,
      .coeff_mcs_end1                             = 1.77801,
      .coeff_mcs_end2                             = 0.0,
      .coeff_delta_group                          = 0.0,
      .coeff_wdelta_group                         = 1.25320,
      .coeff_sequence_length_norm                 = 0.0,
      .coeff_avg_region_length_norm               = 0.54612,
      .coeff_max_region_length_norm               = 0.0,
      .coeff_kmer_frequency                       = 0.12584,
  },
};
const archive_score_logit_model_t logit_model_single_end_sensitive = {
//    Coefficients:
//    (Intercept)          swg        equal        acand         mcs1        arlen
//      -35.76634     30.18393     -0.09304      1.18547      1.85906     -3.10824
//          mrlen         kmer
//        2.89171      1.21661
  .unique_logit_coeff = {
      .coeff_intercept                            = -35.76634,
      .coeff_primary_edit_distance_norm           = 0.0,
      .coeff_primary_event_distance_norm          = 0.0,
      .coeff_primary_swg_score_norm               = 30.18393,
      .coeff_primary_error_quality                = -0.09304,
      .coeff_primary_template_size_sigma_norm     = 0.0,
      .coeff_mapq_end1                            = 0.0,
      .coeff_mapq_end2                            = 0.0,
      .coeff_subdominant_edit_distance_norm       = 0.0,
      .coeff_subdominant_event_distance_norm      = 0.0,
      .coeff_subdominant_swg_score_norm           = 0.0,
      .coeff_subdominant_template_size_sigma_norm = 0.0,
      .coeff_candidates_subdominant               = 0.0,
      .coeff_candidates_accepted                  = 1.18547,
      .coeff_mcs_end1                             = 1.85906,
      .coeff_mcs_end2                             = 0.0,
      .coeff_delta_group                          = 0.0,
      .coeff_wdelta_group                         = 0.0,
      .coeff_sequence_length_norm                 = 0.0,
      .coeff_avg_region_length_norm               = -3.10824,
      .coeff_max_region_length_norm               = 2.89171,
      .coeff_kmer_frequency                       = 1.21661,
  },
//  Coefficients:
//  (Intercept)         edit        event          swg        equal      sub_swg
//   -1.364e+02    8.979e+01    4.058e+01    1.528e+01   -6.194e-03   -1.032e+01
//        scand        acand         mcs1        arlen        mrlen         kmer
//   -1.313e-01   -9.171e-04    1.750e+00   -3.890e-01    6.126e-01    2.264e-01
  .mmap_logit_coeff = {
      .coeff_intercept                            = -1.364e+02,
      .coeff_primary_edit_distance_norm           = 8.979e+01,
      .coeff_primary_event_distance_norm          = 4.058e+01,
      .coeff_primary_swg_score_norm               = 1.528e+01,
      .coeff_primary_error_quality                = -6.194e-03,
      .coeff_primary_template_size_sigma_norm     = 0.0,
      .coeff_mapq_end1                            = 0.0,
      .coeff_mapq_end2                            = 0.0,
      .coeff_subdominant_edit_distance_norm       = 0.0,
      .coeff_subdominant_event_distance_norm      = 0.0,
      .coeff_subdominant_swg_score_norm           = -1.032e+01,
      .coeff_subdominant_template_size_sigma_norm = 0.0,
      .coeff_candidates_subdominant               = -1.313e-01,
      .coeff_candidates_accepted                  = -9.171e-04,
      .coeff_mcs_end1                             = 1.750e+00,
      .coeff_mcs_end2                             = 0.0,
      .coeff_delta_group                          = 0.0,
      .coeff_wdelta_group                         = 0.0,
      .coeff_sequence_length_norm                 = 0.0,
      .coeff_avg_region_length_norm               = -3.890e-01,
      .coeff_max_region_length_norm               = 6.126e-01,
      .coeff_kmer_frequency                       = 2.264e-01,
  },
};
/*
 * PE-Model DEFAULT
 */
const archive_score_logit_model_t logit_model_paired_end_fast = {
  .unique_logit_coeff = {
      .coeff_intercept                            = 0.0,
      .coeff_primary_edit_distance_norm           = 0.0,
      .coeff_primary_event_distance_norm          = 0.0,
      .coeff_primary_swg_score_norm               = 0.0,
      .coeff_primary_error_quality                = 0.0,
      .coeff_primary_template_size_sigma_norm     = 0.0,
      .coeff_mapq_end1                            = 0.0,
      .coeff_mapq_end2                            = 0.0,
      .coeff_subdominant_edit_distance_norm       = 0.0,
      .coeff_subdominant_event_distance_norm      = 0.0,
      .coeff_subdominant_swg_score_norm           = 0.0,
      .coeff_subdominant_template_size_sigma_norm = 0.0,
      .coeff_candidates_subdominant               = 0.0,
      .coeff_candidates_accepted                  = 0.0,
      .coeff_mcs_end1                             = 0.0,
      .coeff_mcs_end2                             = 0.0,
      .coeff_delta_group                          = 0.0,
      .coeff_wdelta_group                         = 0.0,
      .coeff_sequence_length_norm                 = 0.0,
      .coeff_avg_region_length_norm               = 0.0,
      .coeff_max_region_length_norm               = 0.0,
      .coeff_kmer_frequency                       = 0.0,

  },
//  Coefficients:
//  (Intercept)       wdelta         edit        event          swg       mapq_1
//      3.97843      0.17960    138.46925    345.56584     27.90942      0.98446
//       mapq_2     sub_edit    sub_event      sub_swg        scand        arlen
//      0.90512   -149.63586   -339.78928    -26.80168     -0.01971     -0.20189
  .mmap_logit_coeff = {
      .coeff_intercept                            = 3.97843,
      .coeff_primary_edit_distance_norm           = 138.46925,
      .coeff_primary_event_distance_norm          = 345.56584,
      .coeff_primary_swg_score_norm               = 27.90942,
      .coeff_primary_error_quality                = 0.0,
      .coeff_primary_template_size_sigma_norm     = 0.0,
      .coeff_mapq_end1                            = 0.98446,
      .coeff_mapq_end2                            = 0.90512,
      .coeff_subdominant_edit_distance_norm       = -149.63586,
      .coeff_subdominant_event_distance_norm      = -339.78928,
      .coeff_subdominant_swg_score_norm           = -26.80168,
      .coeff_subdominant_template_size_sigma_norm = 0.0,
      .coeff_candidates_subdominant               = -0.01971,
      .coeff_candidates_accepted                  = 0.0,
      .coeff_mcs_end1                             = 0.0,
      .coeff_mcs_end2                             = 0.0,
      .coeff_delta_group                          = 0.0,
      .coeff_wdelta_group                         = 0.17960,
      .coeff_sequence_length_norm                 = 0.0,
      .coeff_avg_region_length_norm               = -0.20189,
      .coeff_max_region_length_norm               = 0.0,
      .coeff_kmer_frequency                       = 0.0,
  },
};
const archive_score_logit_model_t logit_model_paired_end_sensitive = {
  .unique_logit_coeff = {
      .coeff_intercept                            = 0.0,
      .coeff_primary_edit_distance_norm           = 0.0,
      .coeff_primary_event_distance_norm          = 0.0,
      .coeff_primary_swg_score_norm               = 0.0,
      .coeff_primary_error_quality                = 0.0,
      .coeff_primary_template_size_sigma_norm     = 0.0,
      .coeff_mapq_end1                            = 0.0,
      .coeff_mapq_end2                            = 0.0,
      .coeff_subdominant_edit_distance_norm       = 0.0,
      .coeff_subdominant_event_distance_norm      = 0.0,
      .coeff_subdominant_swg_score_norm           = 0.0,
      .coeff_subdominant_template_size_sigma_norm = 0.0,
      .coeff_candidates_subdominant               = 0.0,
      .coeff_candidates_accepted                  = 0.0,
      .coeff_mcs_end1                             = 0.0,
      .coeff_mcs_end2                             = 0.0,
      .coeff_delta_group                          = 0.0,
      .coeff_wdelta_group                         = 0.0,
      .coeff_sequence_length_norm                 = 0.0,
      .coeff_avg_region_length_norm               = 0.0,
      .coeff_max_region_length_norm               = 0.0,
      .coeff_kmer_frequency                       = 0.0,
  },
  .mmap_logit_coeff = {
      .coeff_intercept                            = 3.97843,
      .coeff_primary_edit_distance_norm           = 138.46925,
      .coeff_primary_event_distance_norm          = 345.56584,
      .coeff_primary_swg_score_norm               = 27.90942,
      .coeff_primary_error_quality                = 0.0,
      .coeff_primary_template_size_sigma_norm     = 0.0,
      .coeff_mapq_end1                            = 0.98446,
      .coeff_mapq_end2                            = 0.90512,
      .coeff_subdominant_edit_distance_norm       = -149.63586,
      .coeff_subdominant_event_distance_norm      = -339.78928,
      .coeff_subdominant_swg_score_norm           = -26.80168,
      .coeff_subdominant_template_size_sigma_norm = 0.0,
      .coeff_candidates_subdominant               = -0.01971,
      .coeff_candidates_accepted                  = 0.0,
      .coeff_mcs_end1                             = 0.0,
      .coeff_mcs_end2                             = 0.0,
      .coeff_delta_group                          = 0.0,
      .coeff_wdelta_group                         = 0.17960,
      .coeff_sequence_length_norm                 = 0.0,
      .coeff_avg_region_length_norm               = -0.20189,
      .coeff_max_region_length_norm               = 0.0,
      .coeff_kmer_frequency                       = 0.0,
  },
};

/*

.coeff_intercept                            = 0.0,
.coeff_primary_edit_distance_norm           = 0.0,
.coeff_primary_event_distance_norm          = 0.0,
.coeff_primary_swg_score_norm               = 0.0,
.coeff_primary_error_quality                = 0.0,
.coeff_primary_template_size_sigma_norm     = 0.0,
.coeff_mapq_end1                            = 0.0,
.coeff_mapq_end2                            = 0.0,
.coeff_subdominant_edit_distance_norm       = 0.0,
.coeff_subdominant_event_distance_norm      = 0.0,
.coeff_subdominant_swg_score_norm           = 0.0,
.coeff_subdominant_template_size_sigma_norm = 0.0,
.coeff_candidates_subdominant               = 0.0,
.coeff_candidates_accepted                  = 0.0,
.coeff_mcs_end1                             = 0.0,
.coeff_mcs_end2                             = 0.0,
.coeff_delta_group                          = 0.0,
.coeff_wdelta_group                         = 0.0,
.coeff_sequence_length_norm                 = 0.0,
.coeff_avg_region_length_norm               = 0.0,
.coeff_max_region_length_norm               = 0.0,
.coeff_kmer_frequency                       = 0.0,

 */
