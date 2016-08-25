/*
 * PROJECT: GEMMapper
 * FILE: archive_search_se_stepwise.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include <approximate_search/approximate_search_stages.h>
#include "archive/archive_search_se_stepwise.h"
#include "archive/archive_search_se.h"
#include "archive/archive_select.h"
#include "archive/archive_score_se.h"
#include "archive/archive_check.h"
#include "approximate_search/approximate_search_stepwise.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * Debug
 */
#define DEBUG_ARCHIVE_SEARCH_SE_STEPWISE GEM_DEEP_DEBUG

void archive_search_se_stepwise_debug_preface(
    archive_search_t* const archive_search,
    const char* const label) {
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: %s (stage=%s,%s) (Tag=%s)\n",
        label,asearch_stage_label[archive_search->approximate_search.search_stage],
        asearch_processing_state_label[archive_search->approximate_search.processing_state],
        archive_search->sequence.tag.buffer);
    tab_global_inc();
  }
}
void archive_search_se_stepwise_debug_epilogue() {
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) { tab_global_dec(); }
}
/*
 * Stepwise: Init Search
 *   Reset initial values (Prepare pattern(s), instantiate parameters values, ...)
 */
void archive_search_se_stepwise_init_search(archive_search_t* const archive_search) {
  archive_search_se_stepwise_debug_preface(archive_search,"Init");
  archive_search_reset(archive_search); // Init
  archive_search_se_stepwise_debug_epilogue();
}
/*
 * Stepwise: Region-Profile
 */
void archive_search_se_stepwise_region_profile_static_generate(archive_search_t* const archive_search) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_preface(archive_search,"Region.Profile.Static.Generate.Static");
  // Region-Profile Generate
  approximate_search_stepwise_region_profile_static_generate(&archive_search->approximate_search);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
}
void archive_search_se_stepwise_region_profile_static_copy(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_COPY,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_preface(archive_search,"Region.Profile.Static.Copy");
  // Region-Profile Copy
  approximate_search_stepwise_region_profile_static_copy(
      &archive_search->approximate_search,gpu_buffer_fmi_ssearch);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_COPY,PROFILE_LEVEL);
}
void archive_search_se_stepwise_region_profile_static_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_ssearch_t* const gpu_buffer_fmi_ssearch) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_preface(archive_search,"Region.Profile.Static.Retrieve");
  // Region-Profile Retrieve
  approximate_search_stepwise_region_profile_static_retrieve(
      &archive_search->approximate_search,gpu_buffer_fmi_ssearch);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
}
void archive_search_se_stepwise_region_profile_adaptive_generate(archive_search_t* const archive_search) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_preface(archive_search,"Region.Profile.Adaptive.Generate.Adaptive");
  // Region-Profile Generate
  approximate_search_stepwise_region_profile_adaptive_generate(&archive_search->approximate_search);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
}
void archive_search_se_stepwise_region_profile_adaptive_copy(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_COPY,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_preface(archive_search,"Region.Profile.Adaptive.Copy");
  // Region-Profile Copy
  approximate_search_stepwise_region_profile_adaptive_copy(
      &archive_search->approximate_search,gpu_buffer_fmi_asearch);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_COPY,PROFILE_LEVEL);
}
void archive_search_se_stepwise_region_profile_adaptive_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_asearch_t* const gpu_buffer_fmi_asearch) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_preface(archive_search,"Region.Profile.Adaptive.Retrieve");
  // Region-Profile Retrieve
  approximate_search_stepwise_region_profile_adaptive_retrieve(
      &archive_search->approximate_search,gpu_buffer_fmi_asearch);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: Decode-Candidates
 */
void archive_search_se_stepwise_decode_candidates_copy(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_DECODE_CANDIDATES_COPY,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_preface(archive_search,"Decode.Candidates.Copy");
  // Decode-Candidates Copy
  approximate_search_stepwise_decode_candidates_copy(
      &archive_search->approximate_search,gpu_buffer_fmi_decode);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_DECODE_CANDIDATES_COPY,PROFILE_LEVEL);
}
void archive_search_se_stepwise_decode_candidates_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_DECODE_CANDIDATES_RETRIEVE,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_preface(archive_search,"Decode.Candidates.Retrieve");
  // Decode-Candidates Retrieve
  approximate_search_stepwise_decode_candidates_retrieve(
      &archive_search->approximate_search,gpu_buffer_fmi_decode);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_DECODE_CANDIDATES_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: Verify-Candidates
 */
void archive_search_se_stepwise_verify_candidates_copy(
    archive_search_t* const archive_search,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES_COPY,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_preface(archive_search,"Verify.Candidates.Copy");
  // Verify-Candidates Copy
  approximate_search_stepwise_verify_candidates_copy(
      &archive_search->approximate_search,gpu_buffer_align_bpm);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES_COPY,PROFILE_LEVEL);
}
void archive_search_se_stepwise_verify_candidates_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES_RETRIEVE,PROFILE_LEVEL);
  archive_search_se_stepwise_debug_preface(archive_search,"Verify.Candidates.Retrieve");
  // Verify-Candidates Retrieve
  approximate_search_stepwise_verify_candidates_retrieve(
      &archive_search->approximate_search,gpu_buffer_align_bpm,matches);
  archive_search_se_stepwise_debug_epilogue();
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: Finish Search
 */
void archive_search_se_stepwise_finish_search(
    archive_search_t* const archive_search,
    matches_t* const matches,
    const bool paired_end_search) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH,PROFILE_LEVEL);
  // DEBUG
  archive_search_se_stepwise_debug_preface(archive_search,"Finish");
  // Finish Search
  approximate_search(&archive_search->approximate_search,matches);
  // Select Matches
  search_parameters_t* const search_parameters = &archive_search->search_parameters;
  select_parameters_t* const select_parameters = (paired_end_search) ?
      &search_parameters->select_parameters_align:
      &search_parameters->select_parameters_report;
  archive_select_se_matches(archive_search,select_parameters,matches,paired_end_search);
  // Score Matches (Select alignment-Model and process accordingly)
  archive_score_matches_se(archive_search,matches);
  // Check matches
  if (search_parameters->check_type!=archive_check_nothing) {
    archive_check_se_matches(
        archive_search->archive,search_parameters->alignment_model,
        &search_parameters->swg_penalties,&archive_search->sequence,
        matches,search_parameters->check_type,archive_search->mm_stack);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_global_inc();
    archive_search_se_print(gem_log_get_stream(),archive_search,matches);
    tab_global_dec();
    tab_global_dec();
  }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH,PROFILE_LEVEL);
}
