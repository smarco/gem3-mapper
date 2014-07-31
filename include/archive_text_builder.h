/*
 * PROJECT: GEMMapper
 * FILE: archive_text_builder.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef ARCHIVE_TEXT_BUILDER_H_
#define ARCHIVE_TEXT_BUILDER_H_

#include "essentials.h"
#include "archive_builder.h"
#include "input_multifasta_parser.h"

/*
 * DEBUG
 */
GEM_INLINE void archive_builder_dump_index_text(archive_builder_t* const archive_builder);

/*
 * Archive Builder. Text Generation Building-Blocks
 */
GEM_INLINE void archive_builder_generate_text_add_character(archive_builder_t* const archive_builder,const uint8_t char_enc);
GEM_INLINE void archive_builder_generate_text_add_separator(archive_builder_t* const archive_builder);
GEM_INLINE void archive_builder_generate_text_add_Ns(archive_builder_t* const archive_builder);

GEM_INLINE void archive_builder_generate_text_close_sequence(archive_builder_t* const archive_builder);
GEM_INLINE void archive_builder_generate_text_process_unknowns(archive_builder_t* const archive_builder);

#endif /* ARCHIVE_TEXT_BUILDER_H_ */
