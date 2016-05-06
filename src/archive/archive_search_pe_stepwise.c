/*
 * PROJECT: GEMMapper
 * FILE: archive_search_pe_stepwise.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive/archive_search_pe_stepwise.h"
#include "archive/archive_search.h"
#include "archive/archive_search_parameters.h"
#include "archive/archive_search_pe.h"
#include "archive/archive_select.h"
#include "archive/archive_score_se.h"
#include "archive/archive_score_pe.h"
#include "archive/archive_check.h"
#include "approximate_search/approximate_search_stepwise.h"

/*
 * Debug
 */
#define DEBUG_ARCHIVE_SEARCH_PE_STEPWISE GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PHIGH


/*
 * Stepwise: Init Search
 */
void archive_search_pe_stepwise_init_search(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_INIT,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.PE :: Region.Profile.Generate\n");
    tab_fprintf(gem_log_get_stream(),"  => Tag %s\n",archive_search_end1->sequence.tag.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence.End1 %s\n",archive_search_end1->sequence.read.buffer);
    tab_fprintf(gem_log_get_stream(),"  => Sequence.End2 %s\n",archive_search_end2->sequence.read.buffer);
  }
  // Init
  archive_search_end1->pair_searched = false;
  archive_search_end1->pair_extended = false;
  archive_search_end1->pair_extended_shortcut = false;
  archive_search_end2->pair_searched = false;
  archive_search_end2->pair_extended = false;
  archive_search_reset(archive_search_end1);
  archive_search_reset(archive_search_end2);
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_INIT,PROFILE_LEVEL);
}
/*
 * Stepwise: Region-Profile
 */
void archive_search_pe_stepwise_region_profile_generate_static(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: Region.Profile.Generate (Tag=%s)\n",
        archive_search_end1->sequence.tag.buffer);
    tab_global_inc();
  }
  // Region-Profile Generate
  approximate_search_stepwise_region_profile_generate_static(&archive_search_end1->forward_search_state);
  approximate_search_stepwise_region_profile_generate_static(&archive_search_end2->forward_search_state);
  if (archive_search_end1->emulate_rc_search) {
    approximate_search_stepwise_region_profile_generate_static(&archive_search_end1->reverse_search_state);
    approximate_search_stepwise_region_profile_generate_static(&archive_search_end2->reverse_search_state);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
}
void archive_search_pe_stepwise_region_profile_generate_adaptive(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: Region.Profile.Generate (Tag=%s)\n",
        archive_search_end1->sequence.tag.buffer);
    tab_global_inc();
  }
  // Region-Profile Generate
  approximate_search_stepwise_region_profile_generate_adaptive(&archive_search_end1->forward_search_state);
  approximate_search_stepwise_region_profile_generate_adaptive(&archive_search_end2->forward_search_state);
  if (archive_search_end1->emulate_rc_search) {
    approximate_search_stepwise_region_profile_generate_adaptive(&archive_search_end1->reverse_search_state);
    approximate_search_stepwise_region_profile_generate_adaptive(&archive_search_end2->reverse_search_state);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_GENERATE,PROFILE_LEVEL);
}
void archive_search_pe_stepwise_region_profile_copy(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_COPY,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.PE :: Region.Profile.Copy (Tag=%s)\n",
        archive_search_end1->sequence.tag.buffer);
    tab_global_inc();
  }
  // Region-Profile Copy
  approximate_search_stepwise_region_profile_copy(
      &archive_search_end1->forward_search_state,gpu_buffer_fmi_search);
  approximate_search_stepwise_region_profile_copy(
      &archive_search_end2->forward_search_state,gpu_buffer_fmi_search);
  if (archive_search_end1->emulate_rc_search) {
    approximate_search_stepwise_region_profile_copy(
        &archive_search_end1->reverse_search_state,gpu_buffer_fmi_search);
    approximate_search_stepwise_region_profile_copy(
        &archive_search_end2->reverse_search_state,gpu_buffer_fmi_search);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_COPY,PROFILE_LEVEL);
}
void archive_search_pe_stepwise_region_profile_retrieve(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_search_t* const gpu_buffer_fmi_search) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.PE :: Region.Profile.Retrieve (Tag=%s)\n",
        archive_search_end1->sequence.tag.buffer);
    tab_global_inc();
  }
  // Region-Profile Retrieve
  approximate_search_stepwise_region_profile_retrieve(
      &archive_search_end1->forward_search_state,gpu_buffer_fmi_search);
  approximate_search_stepwise_region_profile_retrieve(
      &archive_search_end2->forward_search_state,gpu_buffer_fmi_search);
  if (archive_search_end1->emulate_rc_search) {
    approximate_search_stepwise_region_profile_retrieve(
        &archive_search_end1->reverse_search_state,gpu_buffer_fmi_search);
    approximate_search_stepwise_region_profile_retrieve(
        &archive_search_end2->reverse_search_state,gpu_buffer_fmi_search);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_RETRIEVE,PROFILE_LEVEL);
}

/*
 * Stepwise: Decode-Candidates
 */
void archive_search_pe_stepwise_decode_candidates_copy(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_DECODE_CANDIDATES_COPY,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: Decode.Candidates.Copy (Tag=%s)\n",
        archive_search_end1->sequence.tag.buffer);
    tab_global_inc();
  }
  // Decode-Candidates Copy
  approximate_search_stepwise_decode_candidates_copy(
      &archive_search_end1->forward_search_state,gpu_buffer_fmi_decode);
  approximate_search_stepwise_decode_candidates_copy(
      &archive_search_end2->forward_search_state,gpu_buffer_fmi_decode);
  if (archive_search_end1->emulate_rc_search) {
    approximate_search_stepwise_decode_candidates_copy(
        &archive_search_end1->reverse_search_state,gpu_buffer_fmi_decode);
    approximate_search_stepwise_decode_candidates_copy(
        &archive_search_end2->reverse_search_state,gpu_buffer_fmi_decode);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_DECODE_CANDIDATES_COPY,PROFILE_LEVEL);
}
void archive_search_pe_stepwise_decode_candidates_retrieve(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_fmi_decode_t* const gpu_buffer_fmi_decode) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_DECODE_CANDIDATES_RETRIEVE,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: Decode.Candidates.Retrieve (Tag=%s)\n",
        archive_search_end1->sequence.tag.buffer);
    tab_global_inc();
  }
  // Decode-Candidates Retrieve
  approximate_search_stepwise_decode_candidates_retrieve(
      &archive_search_end1->forward_search_state,gpu_buffer_fmi_decode);
  approximate_search_stepwise_decode_candidates_retrieve(
      &archive_search_end2->forward_search_state,gpu_buffer_fmi_decode);
  if (archive_search_end1->emulate_rc_search) {
    approximate_search_stepwise_decode_candidates_retrieve(
        &archive_search_end1->reverse_search_state,gpu_buffer_fmi_decode);
    approximate_search_stepwise_decode_candidates_retrieve(
        &archive_search_end2->reverse_search_state,gpu_buffer_fmi_decode);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_DECODE_CANDIDATES_RETRIEVE,PROFILE_LEVEL);
}

/*
 * Stepwise: Verify-Candidates
 */
void archive_search_pe_stepwise_verify_candidates_copy(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_VERIFY_CANDIDATES_COPY,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: Verify.Candidates.Copy (Tag=%s)\n",
        archive_search_end1->sequence.tag.buffer);
    tab_global_inc();
  }
  // Verify-Candidates Copy
  approximate_search_stepwise_verify_candidates_copy(
      &archive_search_end1->forward_search_state,gpu_buffer_align_bpm);
  approximate_search_stepwise_verify_candidates_copy(
      &archive_search_end2->forward_search_state,gpu_buffer_align_bpm);
  if (archive_search_end1->emulate_rc_search) {
    approximate_search_stepwise_verify_candidates_copy(
        &archive_search_end1->reverse_search_state,gpu_buffer_align_bpm);
    approximate_search_stepwise_verify_candidates_copy(
        &archive_search_end2->reverse_search_state,gpu_buffer_align_bpm);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_VERIFY_CANDIDATES_COPY,PROFILE_LEVEL);
}
void archive_search_pe_stepwise_verify_candidates_retrieve(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    gpu_buffer_align_bpm_t* const gpu_buffer_align_bpm,
    matches_t* const matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_VERIFY_CANDIDATES_RETRIEVE,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE :: Verify.Candidates.Retrieve (Tag=%s)\n",
        archive_search_end1->sequence.tag.buffer);
    tab_global_inc();
  }
  // Verify-Candidates Retrieve
  approximate_search_stepwise_verify_candidates_retrieve(
      &archive_search_end1->forward_search_state,gpu_buffer_align_bpm,matches);
  approximate_search_stepwise_verify_candidates_retrieve(
      &archive_search_end2->forward_search_state,gpu_buffer_align_bpm,matches);
  if (archive_search_end1->emulate_rc_search) {
    approximate_search_stepwise_verify_candidates_retrieve(
        &archive_search_end1->reverse_search_state,gpu_buffer_align_bpm,matches);
    approximate_search_stepwise_verify_candidates_retrieve(
        &archive_search_end2->reverse_search_state,gpu_buffer_align_bpm,matches);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) { tab_global_dec(); }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_VERIFY_CANDIDATES_RETRIEVE,PROFILE_LEVEL);
}

/*
 * Stepwise: Finish Search
 */
void archive_search_pe_stepwise_finish_search(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches) {
  PROFILE_START(GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH,PROFILE_LEVEL);
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) {
    tab_fprintf(stderr,"[GEM]>ArchiveSearch.STEPWISE.SE.Finish (Tag=%s)\n",
        archive_search_end1->sequence.tag.buffer);
    tab_global_inc();
  }
  // Finish SE-Search
  approximate_search_stepwise_finish(&archive_search_end1->forward_search_state,paired_matches->matches_end1);
  approximate_search_stepwise_finish(&archive_search_end2->forward_search_state,paired_matches->matches_end2);
  if (archive_search_end1->emulate_rc_search) {
    approximate_search_stepwise_finish(&archive_search_end1->reverse_search_state,paired_matches->matches_end1);
    approximate_search_stepwise_finish(&archive_search_end2->reverse_search_state,paired_matches->matches_end2);
  }
  // Select Matches
  search_parameters_t* const search_parameters = &archive_search_end1->search_parameters;
  archive_select_se_matches(archive_search_end1,&search_parameters->select_parameters_report,paired_matches->matches_end1);
  archive_select_se_matches(archive_search_end2,&search_parameters->select_parameters_report,paired_matches->matches_end2);
  // Score Matches (Select alignment-Model and process accordingly)
  archive_score_matches_se(archive_search_end1,paired_matches->matches_end1);
  archive_score_matches_se(archive_search_end2,paired_matches->matches_end2);
  // Finish PE-Search
  archive_search_end1->pe_search_state = archive_search_pe_recovery;
  archive_search_end1->pair_searched = true;
  archive_search_end2->pair_searched = true;
  archive_search_pe_continue(archive_search_end1,archive_search_end2,paired_matches);
  // PE Select Matches
  archive_select_pe_matches(archive_search_end1,archive_search_end2,
      &search_parameters->select_parameters_report,paired_matches);
  // PE Score (Select alignment-Model and process accordingly)
  archive_score_matches_pe(archive_search_end1,archive_search_end2,paired_matches);
  // PE Check matches
  if (search_parameters->check_type!=archive_check_nothing) {
    archive_check_pe_matches(
        archive_search_end1->archive,search_parameters->alignment_model,
        &search_parameters->swg_penalties,&archive_search_end1->sequence,
        &archive_search_end2->sequence,paired_matches,
        search_parameters->check_type,archive_search_end1->mm_stack);
  }
  // DEBUG
  gem_cond_debug_block(DEBUG_ARCHIVE_SEARCH_PE_STEPWISE) {
    tab_global_inc();
    archive_search_pe_print(gem_log_get_stream(),archive_search_end1,archive_search_end2,paired_matches);
    tab_global_dec();
    tab_global_dec();
  }
  PROFILE_STOP(GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH,PROFILE_LEVEL);
}


