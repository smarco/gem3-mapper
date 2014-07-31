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
  // Free handler
  mm_free(archive);
}
/*
 * Archive Accessors
 */
GEM_INLINE uint64_t archive_get_index_length(const archive_t* const archive) {
  ARCHIVE_CHECK(archive);
  return fm_index_get_length(archive->fm_index);
}
/*
 * Archive Search Setup
 */
GEM_INLINE archive_search_t* archive_search_new(archive_t* const archive) {
  ARCHIVE_CHECK(archive);
  // Allocate handler
  archive_search_t* const archive_search = mm_alloc(archive_search_t);
  // Archive
  archive_search->archive = archive;
  // Sequence
  archive_search->sequence = NULL; // Input
  archive_search->rc_sequence = sequence_new(); // Generated
  // MM
  archive_search->mm_stack = mm_stack_new(mm_pool_get_slab(mm_pool_2MB));
  // Approximate Search
  approximate_search_parameters_t* const search_parameters = &(archive_search->search_parameters);
  approximate_search_parameters_init(search_parameters); // Init parameters (defaults)
  archive_search->forward_search_state = approximate_search_new(
      archive->locator,archive->graph,archive->enc_text,archive->fm_index,
      search_parameters,archive_search->mm_stack);
  archive_search->reverse_search_state = approximate_search_new(
      archive->locator,archive->graph,archive->enc_text,archive->fm_index,
      search_parameters,archive_search->mm_stack);
  // Archive search control (Flow control) [DEFAULTS]
  archive_search->probe_strand = true;
  archive_search->search_reverse = !archive->indexed_complement;
  // Matches
  archive_search->matches = matches_new(mm_pool_get_slab(mm_pool_8MB));
  // Return
  return archive_search;
}
GEM_INLINE void archive_search_clear(archive_search_t* const archive_search) {
  // Clear F/R search state
  approximate_search_clear(archive_search->forward_search_state);
  approximate_search_clear(archive_search->reverse_search_state);
  // Clear matches
  matches_clear(archive_search->matches);
}
GEM_INLINE void archive_search_delete(archive_search_t* const archive_search) {
  // Delete Sequence
  sequence_delete(archive_search->rc_sequence);
  // Delete Approximate Search
  approximate_search_delete(archive_search->forward_search_state);
  approximate_search_delete(archive_search->reverse_search_state);
  // Delete mm_stack
  mm_stack_delete(archive_search->mm_stack);
  // Free handler
  mm_free(archive_search);
}
// [Initialize]
GEM_INLINE void archive_search_prepare_sequence(archive_search_t* const archive_search,sequence_t* const sequence) {
  // Set the sequence of the search
  archive_search->sequence = sequence;
  // Check the index characteristics & generate reverse-complement (if needed)
  if (archive_search->archive->indexed_complement) {
    archive_search->search_reverse = false;
  } else {
    if (archive_search->archive->filter_type == Iupac_colorspace_dna) {
      sequence_generate_reverse_complement(archive_search->sequence,archive_search->rc_sequence);
    } else {
      sequence_generate_reverse(archive_search->sequence,archive_search->rc_sequence);
    }
    archive_search->search_reverse =
        !(string_equals(archive_search->sequence->read,archive_search->rc_sequence->read));
  }
  // Instantiate Search Parameters Values
  approximate_search_instantiate_values(&(archive_search->search_parameters),sequence_get_length(sequence));
  // Generate the pattern(s)
  approximate_search_prepare_pattern(
      archive_search->forward_search_state,&(archive_search->search_parameters),archive_search->sequence);
  if (archive_search->search_reverse) {
    approximate_search_prepare_pattern(
        archive_search->reverse_search_state,&(archive_search->search_parameters),archive_search->rc_sequence);
  }
}
// [Accessors]
GEM_INLINE approximate_search_parameters_t* archive_search_get_search_parameters(archive_search_t* const archive_search) {
  return &(archive_search->search_parameters);
}
GEM_INLINE matches_t* archive_search_get_matches(archive_search_t* const archive_search) {
  return archive_search->matches;
}
/*
 * SingleEnd Indexed Search (SE Online Approximate String Search)
 */
GEM_INLINE void archive_search_single_end(archive_search_t* const archive_search,sequence_t* const sequence) {
  ARCHIVE_SEARCH_CHECK(archive_search);
  SEQUENCE_CHECK(sequence);
  const archive_t* const archive = archive_search->archive;
  approximate_search_parameters_t* const search_parameters = &(archive_search->search_parameters);
  // Prepare pattern(s)
  archive_search_prepare_sequence(archive_search,sequence);
  // Clean Matches
  archive_search_clear(archive_search);
  // Search the pattern(s)
  approximate_search_t* const forward_search_state = archive_search->forward_search_state;
  if (archive->indexed_complement) {
    // Compute the full search
    forward_search_state->stop_search_stage = asearch_end; // Don't stop until search is done
    forward_search_state->search_strand = Forward;
    approximate_search(forward_search_state,archive_search->matches);
  } else {
    // Configure search stage to stop at
    forward_search_state->stop_search_stage =
        (search_parameters->complete_strata_after_best_nominal <
            forward_search_state->max_differences && archive_search->probe_strand) ? asearch_neighborhood : asearch_end;
    // Run the search (FORWARD)
    forward_search_state->search_strand = Forward; // Configure forward search
    approximate_search(forward_search_state,archive_search->matches);
    // Check the number of matches & keep searching
    if (matches_get_num_matches(archive_search->matches) > search_parameters->max_matches) {
      // Give up searching (More matches than requested)
      forward_search_state->max_complete_stratum = 0;
    } else {
      // Keep on searching
      approximate_search_t* const reverse_asearch_state = archive_search->reverse_search_state;
      reverse_asearch_state->stop_search_stage = asearch_end; // Force a full search
      // Run the search (REVERSE)
      if (archive_search->search_reverse) {
        reverse_asearch_state->search_strand = Reverse; // Configure reverse search
        approximate_search(reverse_asearch_state,archive_search->matches);
      }
      // Resume forward search (if not completed before)
      if (forward_search_state->current_search_stage != asearch_end) {
        if (matches_get_num_matches(archive_search->matches) > search_parameters->max_matches) {
          // Give up searching (More matches than requested)
          forward_search_state->max_complete_stratum = 0;
        } else {
          forward_search_state->stop_search_stage = asearch_end;
          forward_search_state->search_strand = Forward; // Configure forward search
          approximate_search(forward_search_state,archive_search->matches);
        }
      }
    }
  }
//  #ifdef GEM_MAPPER_DEBUG
//   const uint64_t check_level = search_params->internal_parameters.check_alignments;
//  extern bool keepon;
//  if (keepon && (check_level&CHECK_ALG_P_CORRECTNESS)) { // Check p-correctness
//     const uint64_t old_misms = search_params->max_mismatches;
//    search_params->max_mismatches = caller_max_misms; // Full check
//    fmi_check_pos_correctness(archive->index,search_params,matches,false,mpool);
//    search_params->max_mismatches = old_misms;
//  }
//  if (keepon && (check_level&CHECK_ALG_I_CORRECTNESS)) { // Check i-correctness
//     const uint64_t old_misms = search_params->max_mismatches;
//    search_params->max_mismatches = caller_max_misms; // Full check
//    fmi_check_int_correctness(archive->index,search_params,matches,false,mpool);
//    search_params->max_mismatches = old_misms;
//  }
//  if (keepon && (check_level&CHECK_ALG_COMPLETENESS) && search_params->num_wildcards==0) { // Check completeness
//     const uint64_t old_misms = search_params->max_mismatches;
//    search_params->max_mismatches = (search_params->fast_mapping_degree>0) ?
//        matches->max_complete_stratum-1 : old_misms; /*caller_max_misms*/ // Full check
//    fmi_check_completeness(archive->index,search_params,matches,false,mpool);
//    search_params->max_mismatches = old_misms;
//  }
//  #endif
}
/*
 * Realigning Matches
 */
GEM_INLINE void archive_realign_match(
    archive_t* const archive,matches_t* const matches,
    match_trace_t* const match_trace,
    const char* const text,const uint64_t length) {
  // TODO
//  match_trace->distance = ?;
//  match_trace->cigar_length = ?;
//  match_trace->cigar_buffer_offset = ?;
}
/*
 * Decoding Matches (Retrieving & Processing matches)
 *   - 1. Expand interval-matches (compacted)
 *   - 2. Sort matches wrt distance
 *   - 3. Transform CIGAR of reverse matches
 */
GEM_INLINE void archive_search_calculate_matches_to_decode(
    matches_t* const matches,
    const uint64_t max_decoded_matches,const uint64_t min_decoded_strata,
    const uint64_t min_reported_matches,const uint64_t max_reported_matches,
    uint64_t* const total_strata_to_decode,uint64_t* const matches_to_decode_last_stratum_out) {
  // Compact counters (Shrink the counters to the last non-zero stratum)
  const uint64_t max_nz_stratum = matches_counters_compact(matches); // Maximum reachable stratum w.r.t counters
  if (max_nz_stratum==0) return;
  // Maximum stratum to decode (increased by @max_decoded_matches)
  const uint64_t* const counters = vector_get_mem(matches->counters,uint64_t);
  uint64_t strata_to_decode, total_matches;
  for (strata_to_decode=0,total_matches=0;strata_to_decode<max_nz_stratum;++strata_to_decode) {
    total_matches += counters[strata_to_decode];
    if (total_matches > max_decoded_matches) {
      total_matches -= counters[strata_to_decode];
      break;
    }
  }
  // Maximum stratum to decode (increased by @min_decoded_strata)
  const uint64_t min_nz_stratum = matches_counters_get_min_matching_stratum(matches);
  const uint64_t mandatory_strata = min_nz_stratum + min_decoded_strata;
  for (;strata_to_decode<max_nz_stratum && strata_to_decode<mandatory_strata;++strata_to_decode) {
    total_matches += counters[strata_to_decode];
  }
  // Maximum stratum to decode (increased by @min_reported_matches)
  for (;strata_to_decode<max_nz_stratum && total_matches<min_reported_matches;++strata_to_decode) {
    total_matches += counters[strata_to_decode];
  }
  // Maximum stratum to decode (lowered by @max_reported_matches & @min_reported_matches)
  for (;strata_to_decode>0;--strata_to_decode) {
    const uint64_t prev_acc = total_matches - counters[strata_to_decode-1];
    if (total_matches <= max_reported_matches || prev_acc < min_reported_matches) break;
    total_matches = prev_acc;
  }
  // Decode matches
  if (total_matches!=0) { // => (strata_to_decode > 0)
    *total_strata_to_decode = strata_to_decode;
    *matches_to_decode_last_stratum_out = UINT64_MAX;
    if (total_matches > max_reported_matches) {
      const uint64_t prev_acc = total_matches - counters[strata_to_decode-1];
      if (prev_acc >= min_reported_matches) {
        *matches_to_decode_last_stratum_out = max_reported_matches - prev_acc;
      } else {
        *matches_to_decode_last_stratum_out = max_reported_matches - prev_acc + (min_reported_matches - prev_acc);
      }
    }
  } else {
    *total_strata_to_decode = 0;
    *matches_to_decode_last_stratum_out = 0;
  }
}
GEM_INLINE void archive_decode_matches(
    archive_t* const archive,matches_t* const matches,
    const uint64_t maximum_match_distance,const uint64_t matches_to_decode_last_stratum) {
  // Count already decoded matches & discard unwanted matches
  uint64_t num_matches_last_stratum = 0;
  match_trace_t* global_matches_it = vector_get_mem(matches->global_matches,match_trace_t);
  VECTOR_ITERATE(matches->global_matches,match_trace,match_trace_num,match_trace_t) {
    if (match_trace->distance <= maximum_match_distance) {
      // Count matches last stratum
      if (match_trace->distance == maximum_match_distance) {
        if (num_matches_last_stratum < matches_to_decode_last_stratum) {
          ++num_matches_last_stratum;
        } else {
          continue; // Too many matches decoded in the last stratum
        }
      }
      // Keep match
      if (global_matches_it != match_trace) {
        *global_matches_it = *match_trace;
      }
      ++global_matches_it;
    }
  }
  vector_update_used(matches->global_matches,global_matches_it); // Update used
  // Decode interval matches until we reach @max_decoded_matches
  VECTOR_ITERATE(matches->interval_matches,match_interval,match_interval_num,match_interval_t) {
    if (num_matches_last_stratum >= matches_to_decode_last_stratum) break;
    // Get next match-interval
    if ((match_interval->lo >= match_interval->hi) || (match_interval->distance > maximum_match_distance)) continue;

    // 0.- Decode position (Use the first position of the interval, lo)
    match_trace_t match_trace;
    match_trace.strand = match_interval->strand;
    match_trace.position = fm_index_lookup(archive->fm_index,match_interval->lo);
    match_trace.trace_offset = text_collection_get_num_traces(matches->text_collection); // TODO: Undestand SANTI??

    // 1.- (Re)Align the matching-text retrieved with the seq-read(pattern) [Retrieve CIGAR string]
    if (match_interval->text == NULL) {
      archive_text_retrieve(archive,match_trace.position,
          match_interval->length /* + delta TODO*/,matches->text_collection);
      // match_interval->text = *** // TODO
    }
    // Realign (if needed)
    if (match_trace.distance > 0) {
      archive_realign_match(archive,matches,&match_trace,match_interval->text,match_interval->length);
    } else {
      match_trace.cigar_buffer_offset = UINT64_MAX;
      match_trace.cigar_length = 0;
    }
    // Store CIGAR reference
    const uint64_t cigar_buffer_offset = match_trace.cigar_buffer_offset;
    const uint64_t cigar_length = match_trace.cigar_length;
    uint64_t match_effective_length = UINT64_MAX;

    // 2.- Reverse CIGAR (If the search was performed in the reverse strand, emulated)
    if (match_trace.strand == Reverse) {
      if (archive->filter_type==Iupac_dna) {
        matches_reverse_CIGAR(matches,cigar_buffer_offset,cigar_length);
      } else { // Iupac_colorspace_dna
        matches_reverse_CIGAR_colorspace(matches,cigar_buffer_offset,cigar_length);
      }
    }
//    /* In this case, we also have to re-encode the mismatches */
//    if (in_colorspace) { // TODO
//      for (i=0;i<mism_num;++i) {
//        misms[i].mismatch=misms_colorspace[misms[i].mismatch];
//          // const ch_t misms_colorspace[CH_T_RANGE]= { ['A']='a', ['C']='b', ['G']='c', ['T']='d', ['N']='e' };
//      }
//    }

    // 3.- Locate-map the match
    location_t location;
    locator_map(archive->locator,match_trace.position,&location);
    match_trace.position = location.position;
    if (location.direction == Reverse) { // Adjust position by the effective length
      match_effective_length = matches_get_effective_length(matches,??,cigar_buffer_offset,cigar_length);
      match_trace.position -= (match_effective_length-1);
    }

    // 4.- Add the match
    matches_add_match_trace_t(matches,&match_trace);
    const bool last_stratum_match = (match_interval->distance == maximum_match_distance);
    if (last_stratum_match) ++num_matches_last_stratum;

    // 5.- Build the rest of the interval
    uint64_t sa_position;
    for (sa_position=match_interval->lo+1;sa_position<match_interval->hi;++sa_position) {
      // Check number of matches (only last stratum)
      if (last_stratum_match && ++num_matches_last_stratum >= matches_to_decode_last_stratum) break;
      // 0.- Decode position
      match_trace.position = fm_index_lookup(archive->fm_index,sa_position);
      // 3.- Locate-map the match
      locator_map(archive->locator,match_trace.position,&location);
      match_trace.position = location.position;
      if (location.direction == Reverse) { // Adjust position by the effective length
        if (match_effective_length==UINT64_MAX) {
          match_effective_length = matches_get_effective_length(matches,??,cigar_buffer_offset,cigar_length);
        }
        match_trace.position -= (match_effective_length-1);
      }
      // 4.- Add the match
      matches_add_match_trace_t(matches,&match_trace);
    }
  }
  // Sort all matches
  matches_sort_by_distance(matches);
}
GEM_INLINE void archive_search_decode_matches(
    archive_search_t* const archive_search,
    const uint64_t max_decoded_matches,const uint64_t min_decoded_strata,
    const uint64_t min_reported_matches,const uint64_t max_reported_matches) {
  // Check if we need to decode sth
  if (max_decoded_matches==0 && min_decoded_strata==0 && min_reported_matches==0) return;
  if (min_reported_matches==0 && max_reported_matches==0) return;
  // Calculate the number of matches to decode wrt input parameters
  matches_t* const matches = archive_search->matches;
  uint64_t strata_to_decode = 0, matches_to_decode_last_stratum = 0;
  archive_search_calculate_matches_to_decode(
      matches,max_decoded_matches,min_decoded_strata,
      min_reported_matches,max_reported_matches,
      &strata_to_decode,&matches_to_decode_last_stratum);
  // Decode matches
  if (strata_to_decode > 0) {
    archive_decode_matches(archive_search->archive,matches,strata_to_decode,matches_to_decode_last_stratum);
  } else {
    vector_set_used(matches->global_matches,0); // Remove all matches
  }
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

