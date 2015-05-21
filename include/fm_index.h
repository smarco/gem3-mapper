/*
 * PROJECT: GEMMapper
 * FILE: fm_index.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides FM-Index functionalities
 */

#ifndef FM_INDEX_H_
#define FM_INDEX_H_

#include "essentials.h"

#include "bwt.h"
#include "sampled_sa.h"
#include "rank_mtable.h"

/*
 * Checks
 */
#define FM_INDEX_CHECK(fm_index) GEM_CHECK_NULL(fm_index)

/*
 * FM-Index
 */
typedef struct {
  /* Meta-info */
  uint64_t text_length;                   // Text length
  uint64_t proper_length;                 // Pl=log(text_length,4)
  /* Sampled SA */
  sampled_sa_t* sampled_sa;               // Sampled SuffixArray positions
  /* BWT */
  rank_mtable_t* rank_table;              // Memoizated intervals
  bwt_t* bwt;                             // BWT forward text
  bwt_reverse_t* bwt_reverse;             // BWT reverse text
} fm_index_t;
// Bi-Directional queries
typedef struct {
  /* Forward Interval (searching forward) */
  uint64_t forward_lo;
  uint64_t forward_hi;
  /* Backward Interval (searching backwards) */
  uint64_t backward_lo;
  uint64_t backward_hi;
} fm_2interval_t;
typedef struct {
  uint64_t pranks[DNA_EXT_RANGE];         // Precomputed eranks
} fm_2erank_elms_t;

/*
 * Builder
 */
GEM_INLINE void fm_index_write(
    fm_t* const file_manager,dna_text_t* const bwt_text,uint64_t* const character_occurrences,
    sampled_sa_builder_t* const sampled_sa,const bool check,const bool verbose);
GEM_INLINE void fm_index_reverse_write(
    fm_t* const file_manager,dna_text_t* const bwt_reverse_text,
    uint64_t* const character_occurrences,const bool check,const bool verbose);

/*
 * Loader
 */
GEM_INLINE fm_index_t* fm_index_read_mem(mm_t* const memory_manager,const bool check);
GEM_INLINE bool fm_index_check(const fm_index_t* const fm_index,const bool verbose);
GEM_INLINE void fm_index_delete(fm_index_t* const fm_index);

/*
 * Accessors
 */
GEM_INLINE uint64_t fm_index_get_length(const fm_index_t* const fm_index);
GEM_INLINE double fm_index_get_proper_length(const fm_index_t* const fm_index);
GEM_INLINE uint64_t fm_index_get_size(const fm_index_t* const fm_index);

/*
 * FM-Index Bidirectional Operators
 */
GEM_INLINE void fm_index_precompute_2erank_forward(
    const fm_index_t* const fm_index,fm_2erank_elms_t* const fm_2erank_elms,const uint64_t position);
GEM_INLINE void fm_index_precompute_2erank_backward(
    const fm_index_t* const fm_index,fm_2erank_elms_t* const fm_2erank_elms,const uint64_t position);
GEM_INLINE void fm_index_precomputed_2query_forward(
    fm_2interval_t* const fm_2interval,fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms,const uint8_t char_enc);
GEM_INLINE void fm_index_precomputed_2query_backward(
    fm_2interval_t* const fm_2interval,fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms,const uint8_t char_enc);
GEM_INLINE void fm_index_2query_forward(
    const fm_index_t* const fm_index,fm_2interval_t* const fm_2interval,const uint8_t char_enc);
GEM_INLINE void fm_index_2query_backward(
    const fm_index_t* const fm_index,fm_2interval_t* const fm_2interval,const uint8_t char_enc);

/*
 * FM-Index Operators
 */
// Compute SA[i]
GEM_INLINE uint64_t fm_index_lookup(const fm_index_t* const fm_index,const uint64_t bwt_position);
// Compute SA^(-1)[i]
GEM_INLINE uint64_t fm_index_inverse_lookup(const fm_index_t* const fm_index,const uint64_t text_position);
// Compute Psi[i]
GEM_INLINE uint64_t fm_index_psi(const fm_index_t* const fm_index,const uint64_t bwt_position);
// Decode fm_index->text[bwt_position..bwt_position+length-1] into @buffer.
GEM_INLINE uint64_t fm_index_decode(
    const fm_index_t* const fm_index,const uint64_t bwt_position,const uint64_t length,char* const buffer);

/*
 * Basic FM-Index search (backwards/forward)
 */
GEM_INLINE void fm_index_bsearch_pure(
    const fm_index_t* const fm_index,const uint8_t* const key,
    uint64_t key_length,uint64_t* const hi_out,uint64_t* const lo_out);
GEM_INLINE void fm_index_reverse_bsearch_pure(
    const fm_index_t* const fm_index,const uint8_t* const key,
    const uint64_t key_length,uint64_t* const hi_out,uint64_t* const lo_out);
GEM_INLINE void fm_index_bsearch(
    const fm_index_t* const fm_index,const uint8_t* const key,
    uint64_t key_length,uint64_t* const hi_out,uint64_t* const lo_out);
GEM_INLINE uint64_t fm_index_bsearch_continue(
    const fm_index_t* const fm_index,const char* const key,const uint64_t key_length,
    const bool* const allowed_repl,uint64_t last_hi,uint64_t last_lo,uint64_t begin_pos,
    const uint64_t end_pos,uint64_t* const res_hi,uint64_t* const res_lo);

/*
 * Display
 */
GEM_INLINE void fm_index_print(FILE* const stream,const fm_index_t* const fm_index);

/*
 * Errors
 */
#define GEM_ERROR_FM_INDEX_WRONG_MODEL_NO "FM-index error. Wrong FM-Index Model %lu (Expected model %lu)"
#define GEM_ERROR_FM_INDEX_INDEX_OOR "FM-index error. Index position %lu out-of-range BWT[0,%lu)"

#endif /* FM_INDEX_H_ */
