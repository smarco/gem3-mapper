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
 *   Archive builder module to parse a MFASTA (multi-FASTA) input file
 *   and generate the text-index (raw text index) and meta-info
 *     1. MultiFASTA Read cycle
 *       1.1 Filter UIPAC-DNA bases
 *       1.2 Strip Ns
 *     2. Generate Locator
 *     3. Generate Index-Text
 *     4. Write Header & Locator
 */

#ifndef ARCHIVE_BUILDER_TEXT_H_
#define ARCHIVE_BUILDER_TEXT_H_

#include "utils/essentials.h"
#include "archive/builder/archive_builder.h"

/*
 * Process DNA-Text
 */
// Generate standard DNA-Text
void archive_builder_text_process(
    archive_builder_t* const archive_builder,
    const bool verbose);
// Apply RL to the text
void archive_builder_text_generate_run_length(
    archive_builder_t* const archive_builder,
    const bool verbose);

/*
 * Display
 */
void archive_builder_text_dump(
    archive_builder_t* const archive_builder,
    const char* const extension);

#endif /* ARCHIVE_BUILDER_TEXT_H_ */
