/*
 * PROJECT: GEMMapper
 * FILE: archive_builder_index.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Build BWT (SA)
 */

#ifndef ARCHIVE_BUILDER_INDEX_H_
#define ARCHIVE_BUILDER_INDEX_H_

#include "utils/essentials.h"
#include "archive/archive_builder.h"

/*
 * Build BWT (SA)
 */
void archive_builder_index_build_bwt(
    archive_builder_t* const archive_builder,
    const bool gpu_index,
    const bool dump_bwt,
    const bool dump_explicit_sa,
    const bool verbose);
void archive_builder_index_build_bwt_reverse(
    archive_builder_t* const archive_builder,
    const bool dump_reverse_indexed_text,
    const bool dump_bwt,
    const bool dump_explicit_sa,
    const bool verbose);

/*
 * Display
 */
void archive_builder_index_print_explicit_sa(
    archive_builder_t* const archive_builder,
    const char* const extension);
void archive_builder_index_print_bwt(
    archive_builder_t* const archive_builder,
    const char* const extension,
    const bool verbose);

#endif /* ARCHIVE_BUILDER_INDEX_H_ */
