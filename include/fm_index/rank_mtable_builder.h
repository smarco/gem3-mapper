/*
 * PROJECT: GEMMapper
 * FILE: rank_mtable_builder.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Rank memoizated queries
 */

#ifndef RANK_MTABLE_BUILDER_H_
#define RANK_MTABLE_BUILDER_H_

#include "utils/essentials.h"
#include "fm_index/bwt.h"
#include "fm_index/rank_mtable.h"

/*
 * Write mtable
 */
void rank_mtable_builder_write(
    fm_t* const file_manager,
    rank_mtable_t* const rank_mtable);

/*
 * Generation
 */
rank_mtable_t* rank_mtable_builder_new(
    const bwt_builder_t* const bwt_builder,
    const bool verbose);
rank_mtable_t* rank_mtable_reverse_builder_new(
    const bwt_reverse_builder_t* const bwt_reverse_builder,
    const bool verbose);

/*
 * Delete
 */
void rank_mtable_builder_delete(rank_mtable_t* const rank_mtable);

#endif /* RANK_MTABLE_BUILDER_H_ */
