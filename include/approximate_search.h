/*
 * PROJECT: GEMMapper
 * FILE: approximate_search.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef APPROXIMATE_SEARCH_H_
#define APPROXIMATE_SEARCH_H_

#include "essentials.h"
#include "dna_string.h"
#include "quality_model.h"

#include "locator.h"
#include "graph_text.h"
#include "fm_index.h"
#include "matches.h"
#include "sequence.h"
#include "pattern.h"

#include "region_profile.h"
#include "interval_set.h"
#include "filtering_candidates.h"


// Approximate Search Internals
typedef enum { check_none, check_pmatches_correctness, check_imatches_correctness, check_correctness, check_completness } check_matches_t;
typedef enum { mapping_incremental_mapping, mapping_adaptive_filtering, mapping_fixed_filtering, mapping_fast, mapping_neighborhood_search } mapping_mode_t;
typedef struct {
  /*
   * Search parameters
   */
  /* Mapping strategy (Mapping mode + properties) */
  mapping_mode_t mapping_mode;
  float fast_mapping_degree;
  uint64_t fast_mapping_degree_nominal;
  /* Qualities */
  quality_model_t quality_model;
  quality_format_t quality_format;
  uint64_t quality_threshold;
  /* Error Model (Regulates the number of Mismatch/Indels) */
  float max_search_error;
  uint64_t max_search_error_nominal;           // Maximum number of error/differences while searching (edit distance)
  float max_filtering_error;
  uint64_t max_filtering_error_nominal;        // Maximum tolerated error at filtering (verifying candidates)
  float complete_strata_after_best;
  uint64_t complete_strata_after_best_nominal; // Maximum complete strata from first matching stratum
  float min_matching_length;
  uint64_t min_matching_length_nominal;        // Minimum mapping segment size (verifying candidates)
  /* Matches search (Regulates the number of matches) */
  uint64_t max_matches;
  /* Replacements (Regulates the bases that can be replaced/mismatched) */
  char mismatch_alphabet[DNA_EXT_RANGE];
  uint64_t mismatch_alphabet_length;
  bool allowed_chars[256];
  bool allowed_enc[DNA_EXT_RANGE];
  /*
   * Internals
   */
  /* Soft Region Profile Parameters */
  uint64_t srp_region_th; // Max. number of candidates allowed per region
  uint64_t srp_max_steps; // Max. number of characters to explore to improve the region
  uint64_t srp_dec_factor; // Decreasing factor per step in region exploration
  uint64_t srp_region_type_th; // Threshold to classify regions {ZERO,NON_ZERO}
  /* Hard Region Profile Parameters */
  uint64_t hrp_region_th;
  uint64_t hrp_max_steps;
  uint64_t hrp_dec_factor;
  uint64_t hrp_region_type_th;
  /* Read recovery parameters */
  uint64_t rrp_region_th;
  uint64_t rrp_max_steps;
  uint64_t rrp_dec_factor;
  uint64_t rrp_region_type_th;
  /* Filtering parameters */
  uint64_t filtering_threshold;
  float filtering_region_factor;
  uint64_t pa_filtering_threshold;
  /* Checkers */
  check_matches_t check_matches;
} approximate_search_parameters_t;
// Approximate Search
typedef enum { asearch_init, asearch_filtering, asearch_neighborhood, asearch_end } approximate_search_stage_t;
typedef struct {
  /* Index Structures & Parameters */
  locator_t* locator;                                 // Sequence Locator
  graph_text_t* graph;                                // Graph (text + graph = hypertext)
  dna_text_t* enc_text;                               // Index-Text (Encoded)
  fm_index_t* fm_index;                               // FM-Index
  pattern_t pattern;                                  // Search Pattern
  approximate_search_parameters_t* search_parameters; // Search Parameters
  /* Search State */
  strand_t search_strand;                             // Current search strand
  bool do_quality_search;                             // Quality search
  approximate_search_stage_t current_search_stage;    // Current Stage of the search
  approximate_search_stage_t stop_search_stage;       // Search stage to stop at
  uint64_t max_differences;
  uint64_t max_complete_stratum;
  bool max_matches_reached;                           // Quick abandon due to maximum matches reached
  /* Search Auxiliary Structures */
  region_profile_t region_profile;             // Region Profile
  filtering_candidates_t filtering_candidates; // Filtering Candidates
  interval_set_t intervals_result;             // Interval Set (to hold intervals from neighborhood searches)
  // TODO
  /* MM */
  mm_stack_t* mm_stack;  // Memory stack
} approximate_search_t;

/*
 * Approximate Search Parameters
 */
GEM_INLINE void approximate_search_parameters_init(approximate_search_parameters_t* const search_parameters);

GEM_INLINE void approximate_search_configure_mapping_strategy(
    approximate_search_parameters_t* const search_parameters,
    const mapping_mode_t mapping_mode,const float mapping_degree);
GEM_INLINE void approximate_search_configure_quality_model(
    approximate_search_parameters_t* const search_parameters,
    const quality_model_t quality_model,const quality_format_t quality_format,const uint64_t quality_threshold);
GEM_INLINE void approximate_search_configure_error_model(
    approximate_search_parameters_t* const search_parameters,
    float max_search_error,float max_filtering_error,
    float complete_strata_after_best,float min_matching_length);
GEM_INLINE void approximate_search_configure_replacements(
    approximate_search_parameters_t* const search_parameters,
    char* const mismatch_alphabet,const uint64_t mismatch_alphabet_length);
GEM_INLINE void approximate_search_configure_matches(
    approximate_search_parameters_t* const search_parameters,const uint64_t max_matches);

GEM_INLINE void approximate_search_instantiate_values(
    approximate_search_parameters_t* const search_parameters,const uint64_t pattern_length);

/*
 * Approximate Search
 */
GEM_INLINE approximate_search_t* approximate_search_new(
    locator_t* const locator,graph_text_t* const graph,dna_text_t* const enc_text,fm_index_t* const fm_index,
    approximate_search_parameters_t* const search_parameters,mm_stack_t* const mm_stack);
GEM_INLINE void approximate_search_clear(approximate_search_t* const approximate_search);
GEM_INLINE void approximate_search_delete(approximate_search_t* const approximate_search);

/*
 * Approximate Search Pattern
 */
GEM_INLINE void approximate_search_prepare_pattern(
    approximate_search_t* const approximate_search,
    const approximate_search_parameters_t* const search_parameters,sequence_t* const sequence);

// ASM-Search!!
GEM_INLINE void approximate_search(approximate_search_t* const approximate_search,matches_t* const matches);



//GEM_INLINE uint64_t fmi_mismatched_search_extend(
//    const fm_index_t* const fm_index,
//    fmi_extend_parameters* const extend_parameters,
//    const idx_t init_position,const idx_t end_position,
//    matches_t* const matches,vector_pool* const mpool);

//void fmi_check_pos_correctness(
//    const _FMI_* const fmi,fmi_search_parameters* const search_params,
//    matches* const matches,bool dump_results,vector_pool* const mpool);
//void fmi_check_int_correctness(
//    const _FMI_* const fmi,fmi_search_parameters* const search_params,
//    matches* const matches,bool dump_results,vector_pool* const mpool);
//void fmi_check_completeness(
//    const _FMI_* const fmi,fmi_search_parameters* const search_params,
//    matches* const matches,bool dump_results,vector_pool* const mpool);

/*
 * Error Msg
 */
#define GEM_ERROR_ASM_REPLACEMENT_EMPTY "Approximate Search. Valid replacements set cannot be empty"


#endif /* APPROXIMATE_SEARCH_H_ */




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////






//typedef struct {
//  ch_t *forward_key;
//  ch_t *reverse_key;
//  uint64_t key_len;
//  ch_t *mismatch_mask;
//  direction_t search_direction;
//  /* Mismatch constraints */
//  uint64_t max_distance;
//  uint64_t min_ins_size;
//  uint64_t max_ins_size;
//  uint64_t base_key_len;
//  uint64_t max_extended_matches;
//  uint64_t min_anchor_size;
//  bool check_duplicates;
//  uint64_t unique_pairing;
//  /* Replacements */
//  slch_t repls[CH_T_RANGE];
//  uint64_t repls_len;
//  bool allowed_chars[CH_T_RANGE];
//  bool allowed_repl[CH_T_RANGE];
//  uint64_t num_wildcards;
//  /* Bit vector for Myers-DP */
//  uint64_t *forward_peq;
//  uint64_t *reverse_peq;
//  /* Key tracking id */
//  uint64_t key_id;
//  /* Internals */
//  uint64_t check_alignments;
//} fm_asm_pe_extend_parameters;
