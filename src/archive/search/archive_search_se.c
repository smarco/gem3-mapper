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
 *   Archive-Search Single-End main module
 */

#include "approximate_search/approximate_search_stages.h"
#include "archive/search/archive_search_se.h"
#include "archive/search/archive_select.h"
#include "archive/search/archive_check.h"
#include "archive/score/archive_score_se.h"
#include "approximate_search/approximate_search_filtering_adaptive.h"
#include "filtering/region_profile/region_profile_fixed.h"
#include "matches/classify/matches_classify.h"

/*
 * Debug
 */
#define DEBUG_ARCHIVE_SEARCH_SE GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * Archive Search SE Continue
 */
void archive_search_se_continue(
    archive_search_t* const archive_search,
    matches_t* const matches) {
  // Run the search
  approximate_search(&archive_search->approximate_search,matches);
}
/*
 * Single-End Indexed Search (SE Online Approximate String Search)
 */
void archive_search_se(
    archive_search_t* const archive_search,
    matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE,PROFILE_LEVEL);
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.SE\n");
    tab_fprintf(gem_log_get_stream(),"  => Tag %s\n",archive_search->sequence->tag.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence %s\n",archive_search->sequence->read.buffer);
    tab_global_inc();
  }
  // Search the pattern(s)
  search_parameters_t* const search_parameters = &archive_search->search_parameters;
  // Compute the full search
  approximate_search(&archive_search->approximate_search,matches);
  // Select Matches
  select_parameters_t* const select_parameters = &search_parameters->select_parameters;
  if (!search_parameters->search_paired_parameters.paired_end_search) {
    archive_select_se_matches(select_parameters,matches);
  }
  // Select alignment-Model and process accordingly
  archive_score_matches_se(archive_search,matches);
  // Check matches
  if (search_parameters->check_type!=archive_check_nothing) {
    archive_check_se_matches(
        archive_search->archive,search_parameters->match_alignment_model,
        &search_parameters->swg_penalties,archive_search->sequence,
        matches,search_parameters->check_type,archive_search->mm_allocator);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE) {
    tab_global_inc();
    archive_search_se_print(gem_log_get_stream(),archive_search,matches);
    tab_global_dec();
    tab_global_dec();
  }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE,PROFILE_LEVEL);
}
/*
 * Display
 */
void archive_search_se_print(
    FILE* const stream,
    archive_search_t* const archive_search,
    matches_t* const matches) {
  tab_fprintf(stream,"[GEM]>ArchiveSearch.SE\n");
  tab_global_inc();
  tab_fprintf(stream,"=> Approximate.Search\n");
  tab_global_inc();
  approximate_search_print(stream,&archive_search->approximate_search);
  tab_global_dec();
  if (matches!=NULL) {
    tab_fprintf(stream,"=> Matches.end\n");
    tab_global_inc();
    matches_print(stream,matches);
    tab_global_dec();
  }
  tab_global_dec();
}
