/*
 * PROJECT: GEMMapper
 * FILE: archive.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive.h"

/*
 * Setup/Loader
 */
GEM_INLINE archive_t* archive_read_mem(mm_t* const memory_manager,const bool do_tests,const bool verbose) {
  // Allocate handler
  archive_t* const archive = mm_alloc(archive_t);
  // Set the memory source
  archive->mm = memory_manager;
  // Load archive meta-data
  const uint64_t archive_model_no = mm_read_uint64(memory_manager);
  gem_cond_fatal_error(archive_model_no!=ARCHIVE_MODEL_NO,ARCHIVE_WRONG_MODEL_NO,archive_model_no,ARCHIVE_MODEL_NO);
  archive->filter_type = mm_read_uint64(archive->mm);
  archive->indexed_complement = mm_read_uint64(archive->mm);
  archive->ns_threshold = mm_read_uint64(archive->mm);
  // Load archive::locator
  archive->locator = locator_read_mem(archive->mm);
  // Load archive::text
  archive->text = archive_text_read_mem(archive->mm);
  // Load archive::fm-index
  archive->fm_index = fm_index_read_mem(archive->mm,do_tests);
  // Verbose
  if (verbose) archive_print(gem_info_get_stream(),archive);
  // Return
  return archive;
}
GEM_INLINE archive_t* archive_read(char* const file_name,const bool do_tests,const bool verbose) {
  // Load the whole archive in memory at once
  mm_t* const mm = mm_bulk_mload_file(file_name,1);
  // Return the loaded archive
  return archive_read_mem(mm,do_tests,verbose);
}
GEM_INLINE void archive_delete(archive_t* const archive) {
  ARCHIVE_CHECK(archive);
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
GEM_INLINE uint64_t archive_get_size(const archive_t* const archive) {
  GEM_NOT_IMPLEMENTED(); // TODO
  return 0;
}
GEM_INLINE uint64_t archive_get_index_length(const archive_t* const archive) {
  ARCHIVE_CHECK(archive);
  return fm_index_get_length(archive->fm_index);
}
/*
 * Display
 */
GEM_INLINE void archive_print(FILE* const stream,const archive_t* const archive) {
  GEM_CHECK_NULL(stream);
  ARCHIVE_CHECK(archive);
  tab_fprintf(stream,"[GEM]>Archive\n");
  // Compute Sizes
  const uint64_t locator_size = locator_get_size(archive->locator);
  const uint64_t fm_index_size = fm_index_get_size(archive->fm_index);
  const uint64_t archive_text_size = archive_text_get_size(archive->text);
  const uint64_t archive_size = locator_size + fm_index_size + archive_text_size;
  // Index type
  if (archive->text->hypertext) {
    tab_fprintf(stream,"  => Index.Type GEM.FM-Index.DNA.Graph\n");
  } else if (archive->text->run_length) {
    tab_fprintf(stream,"  => Index.Type GEM.FM-Index.DNA.RunLength\n");
  } else {
    tab_fprintf(stream,"  => Index.Type GEM.FM-Index.DNA.Classic\n");
  }
  // Filter type
  switch (archive->filter_type) {
    case Iupac_dna:
      tab_fprintf(stream,"  => Index.Filter 'Iupac_dna'\n"); break;
    case Iupac_colorspace_dna:
      tab_fprintf(stream,"  => Index.Filter 'Iupac_colorspace_dna'\n"); break;
    default: GEM_INVALID_CASE(); break;
  }
  // Index-Complement
  tab_fprintf(stream,"  => Indexed.Complement %s\n",(archive->indexed_complement) ? "YES" : "NO");
  // Ns-Threshold
  tab_fprintf(stream,"  => Ns.Threshold %"PRIu64"\n",archive->ns_threshold);
  /*
   * Size Display
   */
  tab_fprintf(stream,"  => Archive.Size %"PRIu64" MB\n",CONVERT_B_TO_MB(archive_size));
  tab_fprintf(stream,"    => Locator.Size %"PRIu64" MB (%2.3f%%)\n",CONVERT_B_TO_MB(locator_size),PERCENTAGE(locator_size,archive_size));
  tab_fprintf(stream,"    => Text.Size %"PRIu64" MB (%2.3f%%)\n",CONVERT_B_TO_MB(archive_text_size),PERCENTAGE(archive_text_size,archive_size));
  tab_fprintf(stream,"    => FM.Index.Size %"PRIu64" MB (%2.3f%%)\n",CONVERT_B_TO_MB(fm_index_size),PERCENTAGE(fm_index_size,archive_size));
  /*
   * Components Display
   */
  // Locator
  locator_print(stream,archive->locator,false);
  // Archive Text
  archive_text_print(stream,archive->text);
  // FM-Index
  fm_index_print(stream,archive->fm_index);
  // Flush
  fflush(stream);
}

