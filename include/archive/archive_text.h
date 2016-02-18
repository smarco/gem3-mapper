/*
 * PROJECT: GEMMapper
 * FILE: archive_text.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef ARCHIVE_TEXT_H_
#define ARCHIVE_TEXT_H_

#include "utils/essentials.h"
#include "archive/sampled_rl.h"
#include "data_structures/text_collection.h"
#include "data_structures/dna_text.h"

/*
 * Archive Text
 */
typedef enum {
  archive_text_regular=0,
  archive_text_run_length=1,
  archive_hypertext=2
} archive_text_type;
typedef struct {
  // Meta-information
  bool run_length;                  // Archive-Text type
  bool explicit_complement;         // Stores explicit RC-text
  uint64_t forward_text_length;     // Total length of the forward text
  // Indexed text
  dna_text_t* enc_text;             // Index-Text
  sampled_rl_t* sampled_rl;         // Sampled RL-Index Positions
} archive_text_t;

/*
 * Builder
 */
void archive_text_write(
    fm_t* const file_manager,
    dna_text_t* const enc_text,
    const bool explicit_complement,
    const uint64_t forward_text_length,
    sampled_rl_t* const sampled_rl,
    const bool verbose);

/*
 * Setup/Loader
 */
archive_text_t* archive_text_read_mem(mm_t* const memory_manager);
void archive_text_delete(archive_text_t* const archive_text);

/*
 * Accessors
 */
uint64_t archive_text_get_size(archive_text_t* const archive_text);
strand_t archive_text_get_position_strand(
    archive_text_t* const archive_text,
    const uint64_t index_position);
uint64_t archive_text_get_unitary_projection(
    archive_text_t* const archive_text,
    const uint64_t index_position);
uint64_t archive_text_get_projection(
    archive_text_t* const archive_text,
    const uint64_t index_position,
    const uint64_t length);

/*
 * Text Retriever
 */
void archive_text_retrieve(
    archive_text_t* const archive_text,
    const uint64_t text_position,
    const uint64_t text_length,
    const bool reverse_complement_text,
    const bool run_length_text,
    text_trace_t* const text_trace,
    mm_stack_t* const mm_stack);
uint64_t archive_text_retrieve_collection(
    archive_text_t* const archive_text,
    const text_collection_t* const text_collection,
    const uint64_t text_position,
    const uint64_t text_length,
    const bool reverse_complement_text,
    const bool run_length_text,
    mm_stack_t* const mm_stack);


/*
 * Display
 */
void archive_text_print(
    FILE* const stream,
    const archive_text_t* const archive_text);

/*
 * Errors
 */
#define GEM_ERROR_ARCHIVE_TEXT_WRONG_MODEL_NO "Archive-Text error. Wrong Archive-Text Model %"PRIu64" (Expected model %"PRIu64")"

#endif /* ARCHIVE_TEXT_H_ */
