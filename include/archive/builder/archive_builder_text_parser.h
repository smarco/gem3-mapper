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
 *   Archive builder module providing basic functions to process a
 *   MFASTA (multi-FASTA) input file and compose the text-index
 *   (raw text index) and determine the individual locator-intervals
 *   (meta-info determining the sequences attributes)
 */

#ifndef ARCHIVE_BUILDER_TEXT_PARSER_H_
#define ARCHIVE_BUILDER_TEXT_PARSER_H_

#include "utils/essentials.h"
#include "archive/builder/archive_builder.h"

/*
 * Archive Builder. Text Generation Building-Blocks
 */
void archive_builder_generate_text_add_character(
    archive_builder_t* const archive_builder,
    const uint8_t char_enc);
void archive_builder_generate_text_add_separator(archive_builder_t* const archive_builder);
void archive_builder_generate_text_add_Ns(archive_builder_t* const archive_builder);

void archive_builder_generate_text_close_sequence(archive_builder_t* const archive_builder);
void archive_builder_generate_text_process_unknowns(archive_builder_t* const archive_builder);

/*
 * Inspect Text
 */
void archive_builder_inspect_text(
    archive_builder_t* const archive_builder,
    input_file_t* const input_multifasta,
    const bool verbose);

/*
 * Generate Text
 */
void archive_builder_generate_forward_text(
    archive_builder_t* const archive_builder,
    input_file_t* const input_multifasta,
    const bool verbose);

/*
 * Generate RC-Text
 */
void archive_builder_generate_rc_text(
    archive_builder_t* const archive_builder,
    const bool verbose);

/*
 * Generate C2T & G2A Texts
 */
void archive_builder_generate_bisulfite_text(
    archive_builder_t* const archive_builder,
    const bool verbose);

#endif /* ARCHIVE_BUILDER_TEXT_PARSER_H_ */
