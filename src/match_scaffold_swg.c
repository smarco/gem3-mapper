/*
 * PROJECT: GEMMapper
 * FILE: match_scaffold_swg.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "match_scaffold_swg.h"
#include "match_scaffold_compose.h"
#include "match_scaffold_levenshtein.h"
#include "align.h"

/*
 * Debug
 */
#define DEBUG_MATCH_SCAFFOLD_SWG   GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PLOW

///*
// * Compose the scaffolding (determine the matching regions based on the SWG score)
// */
//void match_scaffold_smith_waterman_gotoh_compose_matching_regions(
//    matches_t* const matches,match_align_input_t* const align_input,
//    uint64_t key_pos,uint64_t text_pos,const uint64_t cigar_offset,
//    const uint64_t cigar_length,match_align_parameters_t* const align_parameters,
//    match_scaffold_t* const match_scaffold,mm_stack_t* const mm_stack) {
//  // Parameters
//  const swg_penalties_t* const swg_penalties = align_parameters->swg_penalties;
//  cigar_element_t* const cigar_buffer = vector_get_elm(matches->cigar_vector,cigar_offset,cigar_element_t);
//  // Init scaffolding
//  match_scaffold->scaffold_regions = mm_stack_calloc(mm_stack,cigar_length,region_matching_t,false);
//  match_scaffold->num_scaffold_regions = 0;
//  match_scaffold->scaffolding_coverage = 0;
//  match_scaffold->scaffold_type = scaffold_swg;
//  // Traverse CIGAR elements and compose scaffolding
//  int64_t i, local_score = 0, max_local_score = 0;
//  uint64_t max_local_begin_key_pos = key_pos, max_local_end_key_pos = key_pos;
//  uint64_t max_local_begin_text_pos = text_pos, max_local_end_text_pos = text_pos;
//  uint64_t max_local_begin_cigar = 0, max_local_end_cigar = 0;
//  for (i=0;i<cigar_length;) {
//    // Add CIGAR operation
//    switch (cigar_buffer[i].type) {
//      case cigar_match: {
//        const int32_t match_length = cigar_buffer[i].length;
//        local_score += align_swg_score_match(swg_penalties,match_length);
//        text_pos += match_length;
//        key_pos += match_length;
//        break;
//      }
//      case cigar_mismatch:
//        local_score += align_swg_score_mismatch(swg_penalties);
//        ++text_pos;
//        ++key_pos;
//        break;
//      case cigar_ins: {
//        const int32_t indel_length = cigar_buffer[i].length;
//        local_score += align_swg_score_insertion(swg_penalties,indel_length);
//        text_pos += indel_length;
//        break;
//      }
//      case cigar_del: {
//        const int32_t indel_length = cigar_buffer[i].length;
//        local_score += align_swg_score_deletion(swg_penalties,indel_length);
//        key_pos += indel_length;
//        break;
//      }
//      default:
//        GEM_INVALID_CASE();
//        break;
//    }
//    ++i;
//    // Check local score
//    if (local_score > max_local_score) {
//      // Max local score
//      max_local_score = local_score;
//      max_local_end_key_pos = key_pos;
//      max_local_end_text_pos = text_pos;
//      max_local_end_cigar = i;
//    } else if (local_score < 0) {
//      // Store local-max alignment chunk
//      if (max_local_score >= align_parameters->swg_threshold) {
//        match_scaffold_compose_add_matching_approximate(
//            max_local_begin_key_pos,max_local_end_key_pos,
//            max_local_begin_text_pos,max_local_end_text_pos,
//            cigar_offset+max_local_begin_cigar,max_local_end_cigar-max_local_begin_cigar,
//            max_local_score,match_scaffold); // Add approximate matching region
//      }
//      // Reset
//      local_score = 0;
//      max_local_score = 0;
//      max_local_begin_key_pos = key_pos;
//      max_local_begin_text_pos = text_pos;
//      max_local_begin_cigar = i;
//    }
//  }
//  // Check local score
//  if (max_local_score >= align_parameters->swg_threshold) {
//    match_scaffold_compose_add_matching_approximate(
//        max_local_begin_key_pos,max_local_end_key_pos,
//        max_local_begin_text_pos,max_local_end_text_pos,
//        cigar_offset+max_local_begin_cigar,max_local_end_cigar-max_local_begin_cigar,
//        max_local_score,match_scaffold); // Add approximate matching region
//  }
//}
