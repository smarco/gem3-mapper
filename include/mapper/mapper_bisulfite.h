/*
 * PROJECT: GEMMapper
 * FILE: mapper_bisulfite.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MAPPER_BISULFITE_H_
#define MAPPER_BISULFITE_H_

#include "utils/essentials.h"
#include "archive/search/archive_search.h"

void mapper_bisulfite_process_sequence_se(
    archive_search_t* const archive_search,
    search_parameters_t* const search_parameters);
void mapper_bisulfite_process_sequence_pe(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);

void mapper_bisulfite_restore_sequence_se(
    archive_search_t* const archive_search);
void mapper_bisulfite_restore_sequence_pe(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2);

#endif /* MAPPER_BISULFITE_H_ */
