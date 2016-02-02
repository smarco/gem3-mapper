/*
 * PROJECT: GEMMapper
 * FILE: archive_text_rl.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef ARCHIVE_TEXT_RL_H_
#define ARCHIVE_TEXT_RL_H_

#include "essentials.h"
#include "archive_text.h"

/*
 * Run-Length MAX-Length
 */
#define TEXT_RL_MAX_RUN_LENGTH  10

/*
 * Encode RL-Text
 */
void archive_text_rl_encode(
    const uint8_t* const text,const uint64_t text_length,
    uint8_t* const rl_text,uint8_t* const rl_runs,
    uint64_t* const rl_text_length);

/*
 * Translate position
 */
uint64_t archive_text_rl_translate(
    archive_text_t* const archive_text,const uint64_t position_rl,
    mm_stack_t* const mm_stack);

/*
 * Utils
 */
uint64_t archive_text_rl_decode_run_length(uint8_t* const rl_runs,const uint64_t rl_position);
uint64_t archive_text_rl_position_decode_inc(uint8_t* const rl_runs,const uint64_t rl_position);
uint64_t archive_text_rl_position_decode_exl(uint8_t* const rl_runs,const uint64_t rl_position);
uint64_t archive_text_rl_position_encode(
    uint8_t* const rl_runs,const uint64_t rl_text_length,
    const uint64_t text_position);


#endif /* ARCHIVE_TEXT_RL_H_ */
