/*
 * PROJECT: GEMMapper
 * FILE: archive_search_se_stepwise.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive/archive_search_se_stepwise.h"
#include "archive/archive_search_se.h"
#include "archive/archive_select.h"
#include "archive/archive_score_se.h"
#include "archive/archive_check.h"
#include "approximate_search/approximate_search_filtering_stages.h"
#include "approximate_search/approximate_search_stepwise.h"

/*
 * Debug
 */
#define DEBUG_ARCHIVE_SEARCH_SE_STEPWISE GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH

/*
 * Stepwise: Init Search
 *   Reset initial values (Prepare pattern(s), instantiate parameters values, ...)
 */
void archive_search_se_stepwise_init_search(archive_search_t* const archive_search) {
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: Init\n");
    tab_fprintf(gem_log_get_stream(),"  => Tag %s\n",archive_search->sequence.tag.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence %s\n",archive_search->sequence.read.buffer);
    tab_global_inc();
  }
  archive_search_reset(archive_search); // Init
}
/*
 * Stepwise: Region-Profile
 */
void archive_search_se_stepwise_region_profile_generate_static(archive_search_t* const archive_search) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: Region.Profile.Generate.Static (stage=%s,%s) (Tag=%s)\n",
        asearch_stage_label[archive_search->forward_search_state.search_stage],
        asearch_processing_state_label[archive_search->forward_search_state.processing_state],
        archive_search->sequence.tag.buffer);
    tab_global_inc();
  }
  // Region-Profile Generate
  approximate_search_stepwise_region_profile_generate_static(&archive_search->forward_search_state);
  if (archive_search->emulate_rc_search) {
    approximate_search_stepwise_region_profile_generate_static(&archive_search->reverse_search_state);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
}
void archive_search_se_stepwise_region_profile_generate_adaptive(archive_search_t* const archive_search) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: Region.Profile.Generate.Adaptive (stage=%s,%s) (Tag=%s)\n",
        asearch_stage_label[archive_search->forward_search_state.search_stage],
        asearch_processing_state_label[archive_search->forward_search_state.processing_state],
        archive_search->sequence.tag.buffer);
    tab_global_inc();
  }
  // Region-Profile Generate
  approximate_search_stepwise_region_profile_generate_adaptive(&archive_search->forward_search_state);
  if (archive_search->emulate_rc_search) {
    approximate_search_stepwise_region_profile_generate_adaptive(&archive_search->reverse_search_state);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
}
void archive_search_se_stepwise_region_profile_copy(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_COPY,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: Region.Profile.Copy (stage=%s,%s) (Tag=%s)\n",
        asearch_stage_label[archive_search->forward_search_state.search_stage],
        asearch_processing_state_label[archive_search->forward_search_state.processing_state],
        archive_search->sequence.tag.buffer);
    tab_global_inc();
  }
  // Region-Profile Copy
  approximate_search_stepwise_region_profile_copy(
      &archive_search->forward_search_state,gpu_buffer_fmi_search);
  if (archive_search->emulate_rc_search) {
    approximate_search_stepwise_region_profile_copy(
        &archive_search->reverse_search_state,gpu_buffer_fmi_search);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_COPY,PROFILE_LEVEL);
}
void archive_search_se_stepwise_region_profile_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: Region.Profile.Retrieve (stage=%s,%s) (Tag=%s)\n",
        asearch_stage_label[archive_search->forward_search_state.search_stage],
        asearch_processing_state_label[archive_search->forward_search_state.processing_state],
        archive_search->sequence.tag.buffer);
    tab_global_inc();
  }
  // Region-Profile Retrieve
  approximate_search_stepwise_region_profile_retrieve(
      &archive_search->forward_search_state,gpu_buffer_fmi_search);
  if (archive_search->emulate_rc_search) {
    approximate_search_stepwise_region_profile_retrieve(
        &archive_search->reverse_search_state,gpu_buffer_fmi_search);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: Decode-Candidates
 */
void archive_search_se_stepwise_decode_candidates_copy(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_DECODE_CANDIDATES_COPY,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: Decode.Candidates.Copy (stage=%s,%s) (Tag=%s)\n",
        asearch_stage_label[archive_search->forward_search_state.search_stage],
        asearch_processing_state_label[archive_search->forward_search_state.processing_state],
        archive_search->sequence.tag.buffer);
    tab_global_inc();
  }
  // Decode-Candidates Copy
  approximate_search_stepwise_decode_candidates_copy(
      &archive_search->forward_search_state,gpu_buffer_fmi_decode);
  if (archive_search->emulate_rc_search) {
    approximate_search_stepwise_decode_candidates_copy(
        &archive_search->reverse_search_state,gpu_buffer_fmi_decode);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_DECODE_CANDIDATES_COPY,PROFILE_LEVEL);
}
void archive_search_se_stepwise_decode_candidates_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_DECODE_CANDIDATES_RETRIEVE,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: Decode.Candidates.Retrieve (stage=%s,%s) (Tag=%s)\n",
        asearch_stage_label[archive_search->forward_search_state.search_stage],
        asearch_processing_state_label[archive_search->forward_search_state.processing_state],
        archive_search->sequence.tag.buffer);
    tab_global_inc();
  }
  // Decode-Candidates Retrieve
  approximate_search_stepwise_decode_candidates_retrieve(
      &archive_search->forward_search_state,gpu_buffer_fmi_decode);
  if (archive_search->emulate_rc_search) {
    approximate_search_stepwise_decode_candidates_retrieve(
        &archive_search->reverse_search_state,gpu_buffer_fmi_decode);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_DECODE_CANDIDATES_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: Verify-Candidates
 */
void archive_search_se_stepwise_verify_candidates_copy(
    archive_search_t* const archive_search,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES_COPY,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: Verify.Candidates.Copy (stage=%s,%s) (Tag=%s)\n",
        asearch_stage_label[archive_search->forward_search_state.search_stage],
        asearch_processing_state_label[archive_search->forward_search_state.processing_state],
        archive_search->sequence.tag.buffer);
    tab_global_inc();
  }
  // Verify-Candidates Copy
  approximate_search_stepwise_verify_candidates_copy(
      &archive_search->forward_search_state,gpu_buffer_align_bpm);
  if (archive_search->emulate_rc_search) {
    approximate_search_stepwise_verify_candidates_copy(
        &archive_search->reverse_search_state,gpu_buffer_align_bpm);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES_COPY,PROFILE_LEVEL);
}
void archive_search_se_stepwise_verify_candidates_retrieve(
    archive_search_t* const archive_search,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES_RETRIEVE,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: Verify.Candidates.Retrieve (stage=%s,%s) (Tag=%s)\n",
        asearch_stage_label[archive_search->forward_search_state.search_stage],
        asearch_processing_state_label[archive_search->forward_search_state.processing_state],
        archive_search->sequence.tag.buffer);
    tab_global_inc();
  }
  // Verify-Candidates Retrieve
  approximate_search_stepwise_verify_candidates_retrieve(
      &archive_search->forward_search_state,gpu_buffer_align_bpm,matches);
  if (archive_search->emulate_rc_search) {
    approximate_search_stepwise_verify_candidates_retrieve(
        &archive_search->reverse_search_state,gpu_buffer_align_bpm,matches);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES_RETRIEVE,PROFILE_LEVEL);
}
/*
 * Stepwise: Finish Search
 */
void archive_search_se_stepwise_finish_search(
    archive_search_t* const archive_search,
    matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_SE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE.Finish (stage=%s,%s) (Tag=%s)\n",
        asearch_stage_label[archive_search->forward_search_state.search_stage],
        asearch_processing_state_label[archive_search->forward_search_state.processing_state],
        archive_search->sequence.tag.buffer);
    tab_global_inc();
  }
  // Finish Search
  approximate_search_stepwise_finish(&archive_search->forward_search_state,matches);
  if (archive_search->emulate_rc_search) {
    approximate_search_stepwise_finish(&archive_search->reverse_search_state,matches);
  }
  // Select Matches
  search_parameters_t* const search_parameters = &archive_search->search_parameters;
  archive_select_se_matches(archive_search,&search_parameters->select_parameters_report,matches);
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
