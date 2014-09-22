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
  archive->index_type = mm_read_uint64(archive->mm);
  archive->filter_type = mm_read_uint64(archive->mm);
  archive->indexed_complement = mm_read_uint64(archive->mm);
  archive->ns_threshold = mm_read_uint64(archive->mm);
  // Load archive::locator
  archive->locator = locator_read_mem(archive->mm);
  // Load archive::RL-Samples
  if (archive->index_type == fm_dna_run_length) {
    // TODO Load RL-Samples
  }
  // Load archive::graph & archive
  archive->graph = (archive->index_type == fm_dna_graph) ? graph_text_read_mem(archive->mm) : NULL;
  // Load archive::text
  archive->enc_text = dna_text_read_mem(archive->mm);
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
  // Delete RL-Samples
  // TODO
  // Delete Graph
  if (archive->index_type == fm_dna_graph) graph_text_delete(archive->graph);
  // Delete Index-Text
  dna_text_delete(archive->enc_text);
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
GEM_INLINE bool archive_is_indexed_complement(const archive_t* const archive) {
  return archive->indexed_complement;
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
  const uint64_t graph_size = (archive->graph!=NULL) ? graph_text_get_size(archive->graph) : 0;
  const uint64_t text_size = dna_text_get_length(archive->enc_text);
  const uint64_t archive_size = locator_size + fm_index_size + graph_size + text_size;
  // Index type
  switch (archive->index_type) {
    case fm_dna_classic: tab_fprintf(stream,"  => Index.Type GEM.FM-Index.DNA.Classic\n"); break;
    case fm_dna_graph:   tab_fprintf(stream,"  => Index.Type GEM.FM-Index.DNA.Graph\n"); break;
    default: GEM_INVALID_CASE(); break;
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
  if (archive->index_type == fm_dna_graph) {
    tab_fprintf(stream,"  => Indexed.Complement YES (Mandatory)\n");
  } else {
    tab_fprintf(stream,"  => Indexed.Complement %s\n",(archive->indexed_complement) ? "YES" : "NO");
  }
  // Ns-Threshold
  tab_fprintf(stream,"  => Ns.Threshold %lu\n",archive->ns_threshold);
  /*
   * Size Display
   */
  tab_fprintf(stream,"  => Archive.Size %lu MB\n",CONVERT_B_TO_MB(archive_size));
  tab_fprintf(stream,"    => Locator.Size %lu MB (%2.3f%%)\n",
        CONVERT_B_TO_MB(locator_size),PERCENTAGE(locator_size,archive_size));
  tab_fprintf(stream,"    => Text.Size %lu MB (%2.3f%%)\n",
      CONVERT_B_TO_MB(text_size),PERCENTAGE(text_size,archive_size));
  if (archive->index_type == fm_dna_run_length) {
    // TODO STH RL-Samples
  }
  if (archive->index_type == fm_dna_graph) {
    tab_fprintf(stream,"    => Graph.Size %lu MB (%2.3f%%)\n",
        CONVERT_B_TO_MB(graph_size),PERCENTAGE(graph_size,archive_size));
  }
  tab_fprintf(stream,"    => FM.Index.Size %lu GB (%2.3f%%)\n",
      CONVERT_B_TO_GB(fm_index_size),PERCENTAGE(fm_index_size,archive_size));
  /*
   * Components Display
   */
  // Locator
  locator_print(stream,archive->locator,false);
  // Archive Text
  dna_text_print(stream,archive->enc_text);
  // RL-Samples
  if (archive->index_type == fm_dna_run_length) {
    // TODO STH RL-Samples
  }
  // Graph
  if (archive->index_type == fm_dna_graph) {
    graph_text_print(stream,archive->graph,true); // FIXME -> false
  }
  // FM-Index
  fm_index_print(stream,archive->fm_index);
  // Flush
  fflush(stream);
}

