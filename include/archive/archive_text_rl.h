/*
 * PROJECT: GEMMapper
 * FILE: archive_text_rl.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef ARCHIVE_TEXT_RL_H_
#define ARCHIVE_TEXT_RL_H_

#include "utils/essentials.h"
#include "archive/archive_text.h"

/*
 * Run-Length MAX-Length
 */
#define TEXT_RL_MAX_RUN_LENGTH  10

/*
 * Encode RL-Text
 */
void archive_text_rl_encode(
    const uint8_t* const text,
    const uint64_t text_length,
    uint8_t* const rl_text,
    uint64_t* const rl_text_length,
    uint32_t* const rl_runs_acc);

/*
 * Translate position
 */
uint64_t archive_text_rl_position_translate(
    archive_text_t* const archive_text,
    const uint64_t position_rl,
    mm_stack_t* const mm_stack);

/*
 * Utils
 */
uint64_t archive_text_rl_get_run_length(
    const uint32_t* const rl_runs_acc,
    const uint64_t rl_position);

uint64_t archive_text_rl_get_decoded_offset_inc(
    const uint32_t* const rl_runs_acc,
    const uint64_t rl_position);
uint64_t archive_text_rl_get_decoded_offset_exl(
    const uint32_t* const rl_runs_acc,
    const uint64_t rl_position);
uint64_t archive_text_rl_get_decoded_length(
    const uint32_t* const rl_runs_acc,
    const uint64_t rl_position,
    const uint64_t length);

uint64_t archive_text_rl_get_encoded_offset(
    const uint32_t* const rl_runs_acc,
    const uint64_t rl_text_length,
    const uint64_t text_position);


#endif /* ARCHIVE_TEXT_RL_H_ */
