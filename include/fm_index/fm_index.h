/*
 * PROJECT: GEMMapper
 * FILE: fm_index.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides FM-Index functionalities
 */

#ifndef FM_INDEX_H_
#define FM_INDEX_H_

#include "utils/essentials.h"
#include "fm_index/bwt.h"
#include "fm_index/sampled_sa.h"
#include "fm_index/rank_mtable.h"

/*
 * FM-Index
 */
typedef struct {
  /* Meta-info */
  uint64_t text_length;                   // Text length
  uint64_t proper_length;                 // Pl=log(text_length,4)
  bool fm_index_reverse;                  // FM-Index of the reverse text included
  /* Sampled SA */
  sampled_sa_t* sampled_sa;               // Sampled SuffixArray positions
  /* BWT */
  rank_mtable_t* rank_table;              // Memoizated intervals
  bwt_t* bwt;                             // BWT forward text
  rank_mtable_t* rank_table_reverse;      // Memoizated reverse intervals
  bwt_reverse_t* bwt_reverse;             // BWT reverse text
} fm_index_t;

/*
 * Builder
 */
bwt_builder_t* fm_index_write(
    fm_t* const file_manager,
    const bool fm_index_reverse,
    dna_text_t* const bwt_text,
    uint64_t* const character_occurrences,
    sampled_sa_builder_t* const sampled_sa,
    const bool check,
    const bool verbose);
bwt_reverse_builder_t* fm_index_reverse_write(
    fm_t* const file_manager,
    dna_text_t* const bwt_reverse_text,
    uint64_t* const character_occurrences,
    const bool check,
    const bool verbose);

/*
 * Loader
 */
fm_index_t* fm_index_read_mem(mm_t* const memory_manager,const bool check);
void fm_index_delete(fm_index_t* const fm_index);

/*
 * Accessors
 */
uint64_t fm_index_get_length(const fm_index_t* const fm_index);
double fm_index_get_proper_length(const fm_index_t* const fm_index);
uint64_t fm_index_get_size(const fm_index_t* const fm_index);

/*
 * FM-Index Operators
 */
uint64_t fm_index_decode(const fm_index_t* const fm_index,const uint64_t bwt_position);
uint64_t fm_index_encode(const fm_index_t* const fm_index,const uint64_t text_position);
uint64_t fm_index_psi(const fm_index_t* const fm_index,const uint64_t bwt_position);

void fm_index_retrieve_bwt_sampled(
    const fm_index_t* const fm_index,
    uint64_t bwt_position,
    uint64_t* const sampled_bwt_position,
    uint64_t* const lf_dist);
void fm_index_retrieve_sa_sample(
    const fm_index_t* const fm_index,
    const uint64_t sampled_bwt_position,
    const uint64_t lf_dist,
    uint64_t* const text_position);

/*
 * Display
 */
void fm_index_print(
    FILE* const stream,
    const fm_index_t* const fm_index);
void fm_index_decode_print_benchmark(
    FILE* const stream,
    const fm_index_t* const fm_index,
    uint64_t bwt_position);

#endif /* FM_INDEX_H_ */
