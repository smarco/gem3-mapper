/*
 * PROJECT: GEMMapper
 * FILE: approximate_search_metrics.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "approximate_search/approximate_search_metrics.h"

/*
 * Setup
 */
void approximate_search_metrics_init(
    approximate_search_metrics_t* const search_metrics,
    const double proper_length,
    const uint64_t read_length,
    const int32_t swg_match_score) {
  // Search Magnitudes
  search_metrics->proper_length = proper_length;
  search_metrics->read_length = read_length;
  search_metrics->swg_match_score = swg_match_score;
  // Mappability
  search_metrics->num_zero_regions = 0;
  search_metrics->max_region_length = UINT32_MAX;
  search_metrics->mappability_p = 0.0;
  search_metrics->mappability_2p = 0.0;

}
/*
 * Accessors
 */
void approximate_search_metrics_set_max_region_length(
    approximate_search_metrics_t* const search_metrics,
    const uint64_t max_region_length) {
  search_metrics->max_region_length = max_region_length;
}
void approximate_search_metrics_set_num_zero_regions(
    approximate_search_metrics_t* const search_metrics,
    const uint64_t num_zero_regions) {
  search_metrics->num_zero_regions = num_zero_regions;
}
void approximate_search_metrics_set_mappability(
    approximate_search_metrics_t* const search_metrics,
    const double mappability_p,
    const double mappability_2p) {
  search_metrics->mappability_p = mappability_p;
  search_metrics->mappability_2p = mappability_2p;
}
/*
 * Display
 */
void approximate_search_metrics_print(
    FILE* const stream,
    approximate_search_metrics_t* const search_metrics) {
  tab_fprintf(stream,"[GEM]>Approximate.Search.Metrics\n");
  tab_fprintf(stream,"  => Search.Magnitudes\n");
  tab_fprintf(stream,"    => Read.length     %lu\n",search_metrics->read_length);
  tab_fprintf(stream,"    => Proper.Length   %2.3f\n",search_metrics->proper_length);
  tab_fprintf(stream,"    => SWG.Match.Score %lu\n",search_metrics->swg_match_score);
  tab_fprintf(stream,"  => Mappability\n");
  tab_fprintf(stream,"    => Max.Region.length  %lu\n",search_metrics->max_region_length);
  tab_fprintf(stream,"    => Mappability.p      %2.3f\n",search_metrics->mappability_p);
  tab_fprintf(stream,"    => Mappability.2p     %2.3f\n",search_metrics->mappability_2p);
}
