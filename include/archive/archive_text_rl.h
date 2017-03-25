/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Archive text module provides data structures and functions to
 *   manipulate and query a genomic-text (ACGTN) encoded run-length (RL)
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
    mm_allocator_t* const mm_allocator);

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
