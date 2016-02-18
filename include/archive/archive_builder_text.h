/*
 * PROJECT: GEMMapper
 * FILE: archive_builder_text.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Process MultiFASTA file
 *                1. MultiFASTA Read cycle
 *                  1.1 Filter UIPAC-DNA bases
 *                  1.2 Strip Ns
 *                2. Generate Locator
 *                3. Generate Index-Text
 *                4. Write Header & Locator
 */
#ifndef ARCHIVE_BUILDER_TEXT_H_
#define ARCHIVE_BUILDER_TEXT_H_

#include "utils/essentials.h"
#include "archive/archive_builder.h"

/*
 * Process DNA-Text
 */
// Generate standard DNA-Text
void archive_builder_text_process(
    archive_builder_t* const archive_builder,
    input_file_t* const input_multifasta,
    const bool verbose);
// Apply RL to the text
void archive_builder_text_apply_run_length(
    archive_builder_t* const archive_builder,
    const bool verbose);

/*
 * Display
 */
void archive_builder_text_dump(
    archive_builder_t* const archive_builder,
    const char* const extension);

#endif /* ARCHIVE_BUILDER_TEXT_H_ */
