/*
 * PROJECT: GEMMapper
 * FILE: archive_search_pe_stepwise.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_search_pe.h"
#include "archive_search_se.h"
#include "approximate_search_control.h"
#include "archive_search_se_stepwise.h"
#include "archive_select.h"
#include "archive_score_pe.h"
#include "filtering_candidates_extend.h"
#include "paired_matches_classify.h"

/*
 * Debug
 */
#define DEBUG_ARCHIVE_SEARCH_PE GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * PE Extension
 */
void archive_search_pe_generate_extension_candidates(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,const sequence_end_t candidate_end) {
//  // Check
//  archive_t* const archive = archive_search_end1->archive;
//  gem_check(!archive->indexed_complement || archive_search_end1->emulate_rc_search,ARCHIVE_SEARCH_INDEX_COMPLEMENT_REQUIRED);
//  // Generate extension candidates
//  search_parameters_t* const search_parameters = archive_search_end1->as_parameters.search_parameters;
//  mapper_stats_t* const mapper_stats = archive_search_end1->mapper_stats;
//  approximate_search_t* const forward_asearch_end1 = &archive_search_end1->forward_search_state;
//  approximate_search_t* const forward_asearch_end2 = &archive_search_end2->forward_search_state;
//  if (candidate_end==paired_end2) {
//    filtering_candidates_extend_generate_candidates(
//        forward_asearch_end1->filtering_candidates,forward_asearch_end2->filtering_candidates,
//        archive,archive_search_end1->text_collection,&forward_asearch_end1->pattern,&forward_asearch_end2->pattern,
//        search_parameters,mapper_stats,paired_matches,archive_search_end1->mm_stack);
//  } else {
//    filtering_candidates_extend_generate_candidates(
//        forward_asearch_end2->filtering_candidates,forward_asearch_end1->filtering_candidates,
//        archive,archive_search_end1->text_collection,&forward_asearch_end2->pattern,&forward_asearch_end1->pattern,
//        search_parameters,mapper_stats,paired_matches,archive_search_end1->mm_stack);
//  }
}
/*
 * Archive Search PE Stepwise
 */
void archive_search_pe_stepwise_generate_candidates(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
//  PROFILE_START(GP_ARCHIVE_SEARCH_PE,PROFILE_LEVEL);
//  PROFILE_START(GP_ARCHIVE_SEARCH_PE_GENERATE_CANDIDATES,PROFILE_LEVEL);
//  // Beginning of the search (Init)
//  archive_search_end1->pe_search_state = archive_search_pe_begin;
//  archive_search_end1->pair_searched = false;
//  archive_search_end1->pair_extended = false;
//  archive_search_end1->pair_extended_shortcut = false;
//  archive_search_end2->pe_search_state = archive_search_pe_begin;
//  archive_search_end2->pair_searched = false;
//  archive_search_end2->pair_extended = false;
//  // Generate Candidates (End/1)
//  archive_search_se_stepwise_generate_candidates(archive_search_end1);
//  archive_search_end1->pair_searched = true;
//  archive_search_end1->pair_extended = archive_search_pe_is_extension_feasible(archive_search_end1);
//  if (archive_search_end1->pair_extended) {
//    // Extension of end/1
//    PROFILE_START(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT,PROFILE_LEVEL);
//    archive_search_reset(archive_search_end2); // Init end/2
//    archive_search_pe_generate_extension_candidates(
//        archive_search_end1,archive_search_end2,paired_matches,paired_end2); // Generate candidates (from extension)
//    archive_search_end1->pe_search_state = archive_search_pe_find_pairs;
//    PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT,PROFILE_LEVEL);
//  } else {
//    // Generate Candidates (End/2)
//    archive_search_se_stepwise_generate_candidates(archive_search_end2);
//    archive_search_end2->pair_searched = true;
//    archive_search_end2->pair_extended = archive_search_pe_is_extension_feasible(archive_search_end2);
//    if (archive_search_end2->pair_extended) {
//      // Extension of end/1
//      PROFILE_START(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT,PROFILE_LEVEL);
//      archive_search_pe_generate_extension_candidates(
//          archive_search_end1,archive_search_end2,paired_matches,paired_end1); // Generate candidates (from extension)
//      archive_search_end1->pe_search_state = archive_search_pe_find_pairs;
//      PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT,PROFILE_LEVEL);
//    }
//  }
//  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_GENERATE_CANDIDATES,PROFILE_LEVEL);
//  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE,PROFILE_LEVEL);
}
void archive_search_pe_stepwise_finish_search(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
//  // PE search (continue search until the end)
//  PROFILE_START(GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH,PROFILE_LEVEL);
//  archive_search_pe_continue(archive_search_end1,archive_search_end2,paired_matches);
//  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH,PROFILE_LEVEL);
}

