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
 *   Archive module provides data structures and functions to
 *   handle a GEM3 Index
 */

#ifndef ARCHIVE_H_
#define ARCHIVE_H_

#include "utils/essentials.h"
#include "archive/locator.h"
#include "archive/archive_text.h"
#include "fm_index/fm_index.h"

/*
 * Archive Model & Version
 */
#define ARCHIVE_MODEL_NO 5012ull

/*
 * Archive
 */
typedef enum {
  archive_dna_full         = 0,          // Standard DNA Text (forward & RC)
  archive_dna_forward      = 1,          // DNA Forward Text (Only forward DNA Strand)
  archive_dna_bisulfite    = UINT64_MAX, // Bisulfite Text
} archive_type;
typedef struct {
  // Meta-information
  archive_type type;                // Archive type
  bool gpu_index;                   // GPU Index
  uint64_t ns_threshold;            // Stretches of Ns (|Ns| >= ns_threshold) are not indexed
  bool indexed_reverse_text;        // Indexed reverse text (backwards text)
  // Locator
  locator_t* locator;               // Sequence Locator
  // Text
  archive_text_t* text;             // Archive Text
  // Index
  fm_index_t* fm_index;             // FM-Index
  // MM
  mm_t* mm;                         // MM
} archive_t;

/*
 * Archive Loader
 */
archive_t* archive_read(
    char* const file_name,
    const bool read_text_only);
void archive_delete(archive_t* const archive);

/*
 * Archive Accessors
 */
uint64_t archive_get_index_length(const archive_t* const archive);

/*
 * Display
 */
void archive_print(
    FILE* const stream,
    const archive_t* const archive);

#endif /* ARCHIVE_H_ */
