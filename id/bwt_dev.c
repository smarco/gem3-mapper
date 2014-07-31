///*
// * bwt_dev.c
// *
// *  Created on: 20/05/2014
// *      Author: smarco
// */
///*
// * Develop
// */
//GEM_INLINE void indexer_develop_build_bwt() {
//  char* buffer = NULL;
//  size_t buffer_size = 0;
//  ssize_t line_length;
//
//  // Create SA-Builder (SA-sorting) & Text
//  sa_builder_t* const sa_builder = sa_builder_new(parameters.max_memory);
//  cdna_text_t* const text = cdna_text_new(mm_pool_get_slab(mm_pool_32MB));
//
//  // Open input file
//  FILE* const input_file = fopen(parameters.input_multifasta_file_name, "r");
//  gem_cond_fatal_error_msg(input_file==NULL,"Error opening the input file '%s'",input_file);
//
//  // Read input & Count Suffixes
//  sa_builder_count_begin(sa_builder);
//  while ((line_length = getline(&buffer,&buffer_size,input_file)) != -1) {
//    // Check MULTIFASTA tag
//    if (buffer[0]==FASTA_TAG_BEGIN) {
//      sa_builder_count_suffix(sa_builder,ENC_DNA_CHAR_SEP);
//      cdna_text_add_char(text,ENC_DNA_CHAR_SEP);
//      continue;
//    }
//    // Read line
//    uint64_t i;
//    for (i=0;i<line_length-1;++i) {
//      const uint8_t char_enc = dna_encode(buffer[i]);
//      // Add character
//      sa_builder_count_suffix(sa_builder,char_enc);
//      cdna_text_add_char(text,char_enc);
//    }
//  }
//  sa_builder_count_end(sa_builder);
//
//  // Layout suffix positions
//  sa_builder_layout_kmer_blocks(sa_builder);
//  gem_debug_block() {
//    sa_builder_display_stats(stderr,sa_builder);
//  }
//
//  // Write suffix positions
//  char* const input_prefix = gem_strrmext(gem_strdup(parameters.input_multifasta_file_name));
//  char* const sa_file_name = gem_strcat(input_prefix,".sa");
//  sa_builder_store_begin(sa_builder,sa_file_name); // Prepare storage
//  cdna_text_iterator_t iterator;
//  cdna_text_iterator_init(&iterator,text,0);
//  while (!cdna_text_iterator_eoi(&iterator)) {
//    const uint8_t enc = cdna_text_iterator_get_char_encoded(&iterator);
//    sa_builder_store_suffix(sa_builder,enc); // Store SA-pos
//    cdna_text_iterator_next_char(&iterator); // Next
//  }
//  sa_builder_store_end(sa_builder,false); // Flush storage
//
//
//
//  // Free
//  if (buffer) free(buffer);
//  free(input_prefix);
//  free(sa_file_name);
//}
//
