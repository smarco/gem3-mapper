/*
 * PROJECT: GEMMapper
 * FILE: archive_builder_text_parser.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef ARCHIVE_BUILDER_TEXT_PARSER_H_
#define ARCHIVE_BUILDER_TEXT_PARSER_H_

#include "utils/essentials.h"
#include "archive/archive_builder.h"

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
