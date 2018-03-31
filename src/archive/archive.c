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

#include "archive/archive.h"

/*
 * DEBUG
 */
#define ARCHIVE_DEBUG_CHECK_INDEX false

/*
 * Error Messages
 */
#define GEM_ERROR_ARCHIVE_WRONG_MODEL_NO "Archive error. Wrong GEM-Index Model %"PRIu64" (Expected model %"PRIu64"). Please rebuild index (gem-indexer)"
#define GEM_ERROR_ARCHIVE_INDEX_OOB "Archive error. Index position (%"PRIu64") out-of-bounds [0,%"PRIu64")]"

/*
 * Setup/Loader
 */
archive_t* archive_read_mem(mm_t* const memory_manager,const bool read_text_only) {
  // Allocate handler
  archive_t* const archive = mm_alloc(archive_t);
  // Set the memory source
  archive->mm = memory_manager;
  // Load archive meta-data
  const uint64_t archive_model_no = mm_read_uint64(memory_manager);
  gem_cond_error(archive_model_no!=ARCHIVE_MODEL_NO,
      ARCHIVE_WRONG_MODEL_NO,archive_model_no,(uint64_t)ARCHIVE_MODEL_NO);
  archive->type = mm_read_uint64(archive->mm);
  archive->gpu_index = mm_read_uint64(archive->mm);
  archive->ns_threshold = mm_read_uint64(archive->mm);
  // Load archive::locator
  archive->locator = locator_read_mem(archive->mm);
  // Load archive::text
  archive->text = archive_text_read_mem(archive->mm);
  if (read_text_only) return archive;
  // Load archive::fm-index
  archive->fm_index = fm_index_read_mem(archive->mm,ARCHIVE_DEBUG_CHECK_INDEX);
  // Return
  return archive;
}
archive_t* archive_read(char* const file_name,const bool read_text_only) {
  // Load the whole archive in memory at once
  mm_t* const mm = mm_bulk_mload_file(file_name);
  // Return the loaded archive
  return archive_read_mem(mm,read_text_only);
}
void archive_delete(archive_t* const archive) {
  // Delete Locator
  locator_delete(archive->locator);
  // Delete Text
  archive_text_delete(archive->text);
  // Delete FM-index
  fm_index_delete(archive->fm_index);
  // Free MM
  if (archive->mm) mm_bulk_free(archive->mm);
  // Free handler
  mm_free(archive);
}
/*
 * Archive Accessors
 */
uint64_t archive_get_index_length(const archive_t* const archive) {
  return fm_index_get_length(archive->fm_index);
}
/*
 * Display
 */
void archive_print(FILE* const stream,const archive_t* const archive) {
  GEM_CHECK_NULL(stream);
  tab_fprintf(stream,"[GEM]>Archive\n");
  // Compute Sizes
  const uint64_t locator_size = locator_get_size(archive->locator);
  const uint64_t fm_index_size = fm_index_get_size(archive->fm_index);
  const uint64_t archive_text_size = archive_text_get_size(archive->text);
  const uint64_t archive_size = locator_size + fm_index_size + archive_text_size;
  // Archive type
  tab_fprintf(stream,"  => Index.Type GEM.FM-Index.DNA");
  if (archive->text->run_length) {
    fprintf(stream,".RunLength");
  } else {
    fprintf(stream,".Classic");
  }
  switch (archive->type) {
    case archive_dna_full:
      fprintf(stream,".full\n");
      break;
    case archive_dna_forward:
      fprintf(stream,".forwardOnly\n");
      break;
    case archive_dna_bisulfite:
      fprintf(stream,".bisulfite\n");
      break;
  }
  // Ns-Threshold
  tab_fprintf(stream,"  => Ns.Threshold %"PRIu64"\n",archive->ns_threshold);
  // Index-Reverse-Text
  tab_fprintf(stream,"  => Indexed.Reverse.Text %s\n",(archive->indexed_reverse_text) ? "YES" : "NO");
  /*
   * Size Display
   */
  tab_fprintf(stream,"  => Archive.Size %"PRIu64" MB\n",CONVERT_B_TO_MB(archive_size));
  tab_fprintf(stream,"    => Locator.Size %"PRIu64" MB (%2.3f%%)\n",
      CONVERT_B_TO_MB(locator_size),PERCENTAGE(locator_size,archive_size));
  tab_fprintf(stream,"    => Text.Size %"PRIu64" MB (%2.3f%%)\n",
      CONVERT_B_TO_MB(archive_text_size),PERCENTAGE(archive_text_size,archive_size));
  tab_fprintf(stream,"    => FM.Index.Size %"PRIu64" MB (%2.3f%%)\n",
      CONVERT_B_TO_MB(fm_index_size),PERCENTAGE(fm_index_size,archive_size));
  /*
   * Components Display
   */
  // Locator
  locator_print(stream,archive->locator,true);
  // Archive Text
  archive_text_print(stream,archive->text);
  // FM-Index
  fm_index_print(stream,archive->fm_index);
  // Flush
  fflush(stream);
}

