/*
 * PROJECT: GEMMapper
 * FILE: archive_check.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef ARCHIVE_CHECK_H_
#define ARCHIVE_CHECK_H_

#include "archive.h"
#include "sequence.h"
#include "match_align.h"
#include "paired_matches.h"

/*
 * Check type
 */
typedef enum {
  archive_check_nothing,
  archive_check_correct,
  archive_check_correct__first_optimum,
  archive_check_correct__all_optimum,
  archive_check_correct__complete
} archive_check_type;

/*
 * Check Matches
 */
void archive_check_se_matches(
    archive_t* const archive,const alignment_model_t alignment_model,
    swg_penalties_t* swg_penalties,sequence_t* const sequence,
    matches_t* const matches,const archive_check_type check_type,
    mm_stack_t* const mm_stack);
void archive_check_pe_matches(
    archive_t* const archive,const alignment_model_t alignment_model,
    swg_penalties_t* swg_penalties,sequence_t* const sequence_end1,
    sequence_t* const sequence_end2,paired_matches_t* const paired_matches,
    const archive_check_type check_type,mm_stack_t* const mm_stack);

#endif /* ARCHIVE_CHECK_H_ */
