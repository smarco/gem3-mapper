/*
 * PROJECT: GEMMapper
 * FILE: archive_builder.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "archive_builder.h"
#include "archive_text_builder.h"

#include "input_graph_parser.h"

#include "graph_text_builder.h"
#include "graph_text.h"

/*
 * Constants
 */
#define ARCHIVE_GRAPH_TAG_INITIAL_LENGTH          100
#define ARCHIVE_GRAPH_SEQUENCE_TAG_INITIAL_LENGTH 100
#define ARCHIVE_GRAPH_VARIANT_TAG_INITIAL_LENGTH  100
#define ARCHIVE_GRAPH_VARIANT_INITIAL_LENGTH      BUFFER_SIZE_1M


/*
 * Errors
 */
#define GEM_ERROR_ARCHIVE_GRAPH_BUILDER_LINK_NULL "Processing GRAPH error(%"PRI_input_file"). Link is meaningless (ignored)"


/*
 * Archive Builder. HyperText Generation Building-Blocks
 */
GEM_INLINE void archive_builder_generate_hypertext_add_jump(archive_builder_t* const archive_builder) {
  archive_builder_generate_text_add_character(archive_builder,ENC_DNA_CHAR_JUMP);
}
GEM_INLINE void archive_builder_generate_hypertext_check_jump_at_the_end(
    archive_builder_t* const archive_builder,graph_sorted_link_locator_t* const graph_link_locator) {
  if (graph_text_builder_link_locator_cmp(graph_link_locator,archive_builder->parsing_state.text_position)) {
    // Handle pending Ns
    archive_builder_generate_text_process_unknowns(archive_builder);
    // Handle non-overlapping jumps
    if (graph_text_builder_link_locator_has_jump_non_overlapping(graph_link_locator)) {
      // Solve Jump index-position
      graph_text_builder_link_locator_solve_jump_non_overlapping(graph_link_locator,
          archive_builder->parsing_state.index_position);
      // Issue jump
      archive_builder_generate_hypertext_add_jump(archive_builder);
    }
  }
}
GEM_INLINE int64_t archive_builder_generate_hypertext_add_sequence(
    archive_builder_t* const archive_builder,graph_sorted_link_locator_t* const graph_link_locator,
    input_file_t* const input_multifasta,string_t* const tag) {
  // Close sequence
  if (archive_builder->parsing_state.multifasta_read_state==Reading_sequence) {
    // Check jumps at the end of the sequence
    archive_builder_generate_hypertext_check_jump_at_the_end(archive_builder,graph_link_locator);
    // Check empty sequence
    gem_cond_fatal_error(input_multifasta_get_text_sequence_length(&archive_builder->parsing_state)==0,
        MULTIFASTA_SEQ_EMPTY,PRI_input_file_content(input_multifasta));
    // Close sequence
    archive_builder_generate_text_process_unknowns(archive_builder); // Check Ns
    archive_builder_generate_text_close_sequence(archive_builder); // Close last sequence
  }
  // Parse TAG
  input_multifasta_parse_tag(input_multifasta,tag);
  // Add to locator
  const int64_t tag_id = locator_builder_add_sequence(
      archive_builder->locator,string_get_buffer(tag),string_get_length(tag));
  // Open interval
  locator_builder_open_interval(archive_builder->locator,tag_id);
  // Begin new text-sequence (Expect sequence after TAG)
  input_multifasta_state_begin_sequence(&archive_builder->parsing_state);
  // Return
  return tag_id;
}
GEM_INLINE void archive_builder_generate_hypertext_process_character(
    archive_builder_t* const archive_builder,
    graph_sorted_link_locator_t* const graph_link_locator,const char current_char) {
  /*
   * Check Graph Pending Links for Text-Position
   * Types of Jump:
   *   - TYPE I (SNV Overlapping)
   *       Overlaps on character of the reference {ACGTE} over {ACGTN}
   *   - TYPE II (SNV Non-Overlapping)
   *       Creates a one-lasso at the position. {ACGTE} over {E}
   *   - TYPE III (General Jump)
   *       General jump to another position and/or sequence
   *  Rules of merging jumps
   *    Jump of the same type can be merged
   *    Jumps at a same position that modify the reference length can be merged (TYPE II and III can be merged)
   *    Jumps overlapping a character in the reference covers it in the index
   */
  if (graph_text_builder_link_locator_cmp(graph_link_locator,archive_builder->parsing_state.text_position)) {
    // Handle pending Ns
    archive_builder_generate_text_process_unknowns(archive_builder);
    // Handle non-overlapping jumps
    if (graph_text_builder_link_locator_has_jump_non_overlapping(graph_link_locator)) {
      // Solve Jump index-position
      graph_text_builder_link_locator_solve_jump_non_overlapping(graph_link_locator,
          archive_builder->parsing_state.index_position);
      // Issue jump
      archive_builder_generate_hypertext_add_jump(archive_builder);
    }
    // Handle overlapping jumps
    if (graph_text_builder_link_locator_has_jump_overlapping(graph_link_locator)) {
      // Solve Jump index-position
      graph_text_builder_link_locator_solve_jump_overlapping(graph_link_locator,
          archive_builder->parsing_state.index_position,current_char);
      // Issue jump
      archive_builder_generate_hypertext_add_jump(archive_builder);
    } else {
      // Add character
      archive_builder_generate_text_add_character(archive_builder,dna_encode(current_char));
    }
    ++(archive_builder->parsing_state.text_interval_length); // Inc Interval-Text Length
  } else {
    // Check Character
    if (current_char==DNA_CHAR_N || !is_extended_dna(current_char)) { // Handle Ns
      ++(archive_builder->parsing_state.ns_pending);
    } else { // Other character
      // Handle pending Ns
      archive_builder_generate_text_process_unknowns(archive_builder);
      // Add character (k-mer counting)
      archive_builder_generate_text_add_character(archive_builder,dna_encode(current_char));
      ++(archive_builder->parsing_state.text_interval_length); // Inc Interval-Text Length
    }
  }
}
/*
 * Archive Builder: Generate Graph
 */
GEM_INLINE void archive_builder_generate_graph_locate_link(
    archive_builder_t* const archive_builder,graph_link_t* const graph_link,string_t* const sequence_name) {
  graph_link->position_index = GRAPH_TEXT_BUILDER_POSITION_UNKNOWN; // FIXME Not necessarily unknown
  graph_link->tag_id *= locator_builder_add_sequence(
      archive_builder->locator,string_get_buffer(sequence_name),string_get_length(sequence_name));
}
GEM_INLINE void archive_builder_generate_graph_add_snv(
    archive_builder_t* const archive_builder,
    graph_link_t* const graph_snv,string_t* const sequence_name,
    const char character,const bool overlaps_reference) {
  // Locate SNV
  archive_builder_generate_graph_locate_link(archive_builder,graph_snv,sequence_name);
  // Add link to the graph builder
  graph_text_builder_add_snv_link(archive_builder->graph,graph_snv,character,overlaps_reference);
}
GEM_INLINE void archive_builder_generate_graph_add_simple_link(
    archive_builder_t* const archive_builder,
    graph_link_t* const graph_link_from,string_t* const sequence_name_from,
    graph_link_t* const graph_link_to,string_t* const sequence_name_to) {
  // Locate links
  archive_builder_generate_graph_locate_link(archive_builder,graph_link_from,sequence_name_from);
  archive_builder_generate_graph_locate_link(archive_builder,graph_link_to,sequence_name_to);
  // Add link to the graph builder
  graph_text_builder_add_general_link(archive_builder->graph,graph_link_from,graph_link_to);
}
GEM_INLINE int64_t archive_builder_generate_graph_add_variant_text(
    archive_builder_t* const archive_builder,string_t* const variant_tag,string_t* const variant) {
  // Add the sequence to the index_text & locator
  const int64_t tag_id = locator_builder_add_sequence(
      archive_builder->locator,string_get_buffer(variant_tag),string_get_length(variant_tag));
  locator_builder_open_interval(archive_builder->locator,tag_id);
  archive_builder_generate_hypertext_add_jump(archive_builder); // Jump BEGIN
  STRING_ITERATE(variant,seq_buffer,pos) {
    archive_builder_generate_text_add_character(archive_builder,dna_encode(seq_buffer[pos]));
  }
  archive_builder_generate_hypertext_add_jump(archive_builder); // Jump END
  locator_builder_close_interval(archive_builder->locator,
      string_get_length(variant),string_get_length(variant)+2,locator_interval_variant); // (+2) Account jumps
  archive_builder_generate_text_add_separator(archive_builder); // Add Separator
  locator_builder_skip_index(archive_builder->locator,1); // Skip Separator
  // Return Tag-ID
  return tag_id;
}
GEM_INLINE void archive_builder_generate_graph_add_variant_link(
    archive_builder_t* const archive_builder,string_t* const variant_tag,string_t* const variant,
    graph_link_t* const graph_link_from,string_t* const sequence_name_from,
    graph_link_t* const graph_link_to,string_t* const sequence_name_to) {
  const uint64_t variant_length = string_get_length(variant);
  graph_link_t graph_link_sequence_in, graph_link_sequence_out;
  // Configure the pos-links to the variant
  graph_link_sequence_in.position_text = 0;
  graph_link_sequence_out.position_text = variant_length;
  graph_link_sequence_in.position_index = archive_builder->parsing_state.index_position;
  graph_link_sequence_out.position_index = archive_builder->parsing_state.index_position+variant_length+1; // (+1) Account jumps
  // Index the variant
  const int64_t tag_id = archive_builder_generate_graph_add_variant_text(archive_builder,variant_tag,variant);
  // Configure the tag-links to the variant
  graph_link_sequence_in.tag_id = tag_id;  // Variant is forward by definition
  graph_link_sequence_out.tag_id = tag_id; // Variant is forward by definition
  // Add links to the graph builder
  archive_builder_generate_graph_locate_link(archive_builder,graph_link_from,sequence_name_from);
  archive_builder_generate_graph_locate_link(archive_builder,graph_link_to,sequence_name_to);
  graph_text_builder_add_general_link(archive_builder->graph,graph_link_from,&graph_link_sequence_in);
  graph_text_builder_add_general_link(archive_builder->graph,&graph_link_sequence_out,graph_link_to);
}
/*
 * STEP0 Archive Build :: Process GRAPH file
 *   1. GRAPH Read cycle
 *   2. Add graph-edge sequences to Locator
 *   3. Add graph-edge sequences to Index-Text
 *   4. Count k-mers (Classify SA positions)
 */
GEM_INLINE void archive_builder_process_graph(
    archive_builder_t* const archive_builder,input_file_t* const input_graph,
    const bool dump_graph_links,const bool verbose) {
  // Allocate archive graph and
  archive_builder->graph = graph_text_builder_new(mm_pool_get_slab(mm_pool_8MB));
  // Allocate temporal variables
  string_t sequence_name_from, sequence_name_to, variant_tag, variant;
  string_init(&sequence_name_from,ARCHIVE_GRAPH_SEQUENCE_TAG_INITIAL_LENGTH);
  string_init(&sequence_name_to,ARCHIVE_GRAPH_SEQUENCE_TAG_INITIAL_LENGTH);
  string_init(&variant_tag,ARCHIVE_GRAPH_VARIANT_TAG_INITIAL_LENGTH);
  string_init(&variant,ARCHIVE_GRAPH_VARIANT_INITIAL_LENGTH);
  graph_link_t graph_link_from, graph_link_to;
  // Read graph links loop
  ticker_t ticker;
  ticker_count_reset(&ticker,verbose,"Reading GRAPH",0,1000000,true);
  ticker_add_process_label(&ticker,"","graph-links parsed");
  ticker_add_finish_label(&ticker,"Total","graph-links parsed");
  while (!input_file_eof(input_graph)) {
    // Parse graph link
    error_code_t error_code;
    if ((error_code=input_graph_link_parse(input_graph,&variant_tag,
        &sequence_name_from,&graph_link_from,&sequence_name_to,&graph_link_to,&variant))) {
      input_graph_link_parse_prompt_error(input_graph,error_code);
      input_file_skip_eol(input_graph);
      continue;
    }
    input_file_skip_eol(input_graph);
    // Process link information
    ticker_update(&ticker,1);
    const uint64_t variant_length = string_get_length(&variant);
    /*
     * SNV
     */
    if (string_cmp(&sequence_name_from,&sequence_name_to)==0 && variant_length<=1) {
      int64_t abs_difference = ((int64_t)graph_link_from.position_text-(int64_t)graph_link_to.position_text);
      abs_difference = ABS(abs_difference);
      if (abs_difference <= 1) { // Is a SNV
        // Check Null link
        if (abs_difference==0 && variant_length==0) {
          gem_error(ARCHIVE_GRAPH_BUILDER_LINK_NULL,PRI_input_file_content(input_graph));
          continue; // Next
        }
        // Add SNV
        graph_link_from.position_text = MIN(graph_link_from.position_text,graph_link_to.position_text);
        archive_builder_generate_graph_add_snv(
            archive_builder,&graph_link_from,&sequence_name_from,
            (variant_length>0)?(*string_get_buffer(&variant)):0,abs_difference==1);
        // Next
        continue;
      }
    }
    /*
     * General graph link
     */
    if (variant_length==0) {
      // Add simple link
      archive_builder_generate_graph_add_simple_link(archive_builder,
          &graph_link_from,&sequence_name_from,&graph_link_to,&sequence_name_to);
      // Next
      continue;
    }
    /*
     * General graph link (through a variant)
     */
    archive_builder_generate_graph_add_variant_link(
        archive_builder,&variant_tag,&variant,
        &graph_link_from,&sequence_name_from,&graph_link_to,&sequence_name_to);
  }
  ticker_finish(&ticker);
  // Check if graph is needed
  if (graph_text_builder_get_num_links(archive_builder->graph) > 0) {
    // Set graph type
//    archive_builder->index_type = fm_dna_graph; // This index has a graph // FIXME: Graph
    // DEBUG: Print parsed links
    graph_text_builder_print(gem_info_get_stream(),archive_builder->graph,archive_builder->locator,dump_graph_links);
    // Sort table
    if (verbose) tfprintf(gem_log_get_stream(),"[Sorting Graph Links]\n");
    graph_text_builder_link_table_sort(archive_builder->graph,locator_builder_get_num_tags(archive_builder->locator));
  }
  // Free
  string_destroy(&variant_tag);
  string_destroy(&sequence_name_from);
  string_destroy(&sequence_name_to);
  string_destroy(&variant);
}
/*
 * Archive Builder. Generate HyperText & RC-HyperText
 */
GEM_INLINE void archive_builder_generate_hypertext(
    archive_builder_t* const archive_builder,input_file_t* const input_multifasta,const bool verbose) {
  string_t tag;
  string_init(&tag,100);
  // Check MultiFASTA (Empty/wrong first character)
  input_file_check_buffer(input_multifasta);
  gem_cond_fatal_error(input_file_eof(input_multifasta),MULTIFASTA_EMPTY,input_file_get_file_name(input_multifasta));
  gem_cond_fatal_error((char)input_file_get_current_char(input_multifasta)!=FASTA_TAG_BEGIN,
      MULTIFASTA_BEGINNING_TAG,PRI_input_file_content(input_multifasta));
  // Prepare ticket
  ticker_t ticker;
  ticker_count_reset(&ticker,verbose,"Reading MultiFASTA & Annotating using GRAPH",0,10000000,true);
  ticker_add_process_label(&ticker,"","bases parsed");
  ticker_add_finish_label(&ticker,"Total","bases parsed");
  // MultiFASTA Read cycle
  graph_sorted_link_locator_t graph_link_locator;
  input_multifasta_state_t* const parsing_state = &(archive_builder->parsing_state);
  while (!input_file_eof(input_multifasta)) {
    char current_char = input_file_get_current_char(input_multifasta);
    if (IS_ANY_EOL(current_char)) { // Empty line (silent skip)
      input_file_skip_eol(input_multifasta);
    } else if (current_char==FASTA_TAG_BEGIN) { // Tag
      // Archive builder parse Tag
      const int64_t tag_id = archive_builder_generate_hypertext_add_sequence(
          archive_builder,&graph_link_locator,input_multifasta,&tag);
      // Search Tag in the graph-builder (looking for unsolved jumps)
      graph_text_builder_link_locator_iterate_forward(archive_builder->graph,&graph_link_locator,tag_id);
    } else if (is_iupac_code(current_char)) { // Nucleotide -> A,C,G,T,N
      do {
        // Process Character
        archive_builder_generate_hypertext_process_character(archive_builder,&graph_link_locator,current_char);
        ticker_update(&ticker,1); // Tick
        // Next Character
        ++(parsing_state->text_position); // Inc Text Position
        if (gem_expect_false(!input_file_next_char(input_multifasta))) break;
        current_char = input_file_get_current_char(input_multifasta);
      } while (!IS_ANY_EOL(current_char));
      input_file_skip_eol(input_multifasta);
      parsing_state->multifasta_read_state = Reading_sequence;
    } else { // Wrong character
      gem_fatal_error(MULTIFASTA_INVALID_CHAR,PRI_input_file_content(input_multifasta),current_char);
    }
  }
  // Check jumps at the end of the sequence
  archive_builder_generate_hypertext_check_jump_at_the_end(archive_builder,&graph_link_locator);
  // Close sequence
  archive_builder_generate_text_close_sequence(archive_builder);
  // Ticker banner
  ticker_finish(&ticker);
  // Check sequence not null
  gem_cond_fatal_error(input_multifasta_get_text_sequence_length(parsing_state)==0,
      MULTIFASTA_SEQ_EMPTY,PRI_input_file_content(input_multifasta));
  // Free
  string_destroy(&tag);
}
GEM_INLINE void archive_builder_generate_rc_hypertext_inteval(
    archive_builder_t* const archive_builder,locator_interval_t* const locator_interval) {
  const int64_t tag_id = locator_interval->tag_id;
  const uint64_t index_interval_length = locator_interval_get_index_length(locator_interval);
  const uint64_t text_interval_length = locator_interval_get_text_length(locator_interval);
  // Locate Reverse Jumps
  graph_sorted_link_locator_t forward_graph_link_locator, reverse_graph_link_locator;
  graph_text_builder_link_locator_iterate_forward(archive_builder->graph,&reverse_graph_link_locator,-tag_id);
  graph_text_builder_link_locator_iterate_backward(archive_builder->graph,&forward_graph_link_locator,tag_id);
  // Append reverse/complement of the interval-text
  cdna_text_iterator_t rev_iterator;
//  cdna_text_reverse_iterator_init(&rev_iterator,archive_builder->text,locator_interval->end_position-1); // FIXME FIXME FIXME FIXME
  int64_t forward_text_position = text_interval_length-1;
  int64_t reverse_text_position = 0;
  // Reset Index Interval Length
  archive_builder->parsing_state.index_interval_length = 0;
  archive_builder->parsing_state.text_interval_length = 0;
  uint64_t interval_index_position = 0;
  // Skip Jump at the end (over the text length)
  if (cdna_text_reverse_iterator_get_char_encoded(&rev_iterator) == ENC_DNA_CHAR_JUMP) {
    GEM_INTERNAL_CHECK(
        graph_text_builder_link_locator_cmp(&forward_graph_link_locator,forward_text_position+1),
        "Inconsistent RC-graph generation (No jump corresponding to graph_link_locator)");
    cdna_text_reverse_iterator_next_char(&rev_iterator);
    interval_index_position++;
  }
  // Reverse traverse all interval
  while ((interval_index_position++) < index_interval_length) {
    // Check Jumps in the forward interval
    uint8_t enc_char = cdna_text_reverse_iterator_get_char_encoded(&rev_iterator);
    if (enc_char == ENC_DNA_CHAR_JUMP) { // Jumps emitted in the forward strand
      if (forward_text_position < 0) {
        // Skip Jump at the beginning
        cdna_text_reverse_iterator_next_char(&rev_iterator); continue;
      } else if (graph_text_builder_link_locator_cmp(&forward_graph_link_locator,forward_text_position)) {
        if (graph_text_builder_link_locator_has_jump_overlapping(&forward_graph_link_locator)) {
          // Extract reference character from overlapping jump
          enc_char = dna_encode(graph_text_builder_link_locator_jump_overlapping_get_reference_char(&forward_graph_link_locator));
          cdna_text_reverse_iterator_next_char(&rev_iterator);
          if (cdna_text_reverse_iterator_get_char_encoded(&rev_iterator) == ENC_DNA_CHAR_JUMP) {
            cdna_text_reverse_iterator_next_char(&rev_iterator); // Skip non-overlapping jump
            ++interval_index_position;
          }
        } else {
          // Skip non-overlapping jump
          cdna_text_reverse_iterator_next_char(&rev_iterator); continue;
        }
      } else if (graph_text_builder_link_locator_cmp(&forward_graph_link_locator,forward_text_position+1)) {
        GEM_INTERNAL_CHECK(!graph_text_builder_link_locator_has_jump_overlapping(&forward_graph_link_locator),
            "Inconsistent RC-graph generation (No jump corresponding to graph_link_locator)");
        // Skip non-overlapping jumps
        cdna_text_reverse_iterator_next_char(&rev_iterator); continue;
      } else {
        GEM_INTERNAL_CHECK(false,"Inconsistent RC-graph generation (No jump corresponding to graph_link_locator)");
      }
    } else {
      cdna_text_reverse_iterator_next_char(&rev_iterator);
    }
    ++(archive_builder->parsing_state.text_interval_length);
    // Check Graph Pending Links for Text-Position
    if (graph_text_builder_link_locator_cmp(&reverse_graph_link_locator,reverse_text_position)) {
      // Handle non-overlapping jumps
      if (graph_text_builder_link_locator_has_jump_non_overlapping(&reverse_graph_link_locator)) {
        // Solve Jump index-position
        graph_text_builder_link_locator_solve_jump_non_overlapping(&reverse_graph_link_locator,
            archive_builder->parsing_state.index_position);
        // Issue jump
        archive_builder_generate_hypertext_add_jump(archive_builder);
      }
      // Handle overlapping jumps
      if (graph_text_builder_link_locator_has_jump_overlapping(&reverse_graph_link_locator)) {
        // Solve Jump index-position
        graph_text_builder_link_locator_solve_jump_overlapping(&reverse_graph_link_locator,
            archive_builder->parsing_state.index_position,dna_complement(dna_decode(enc_char)));
        // Issue jump
        archive_builder_generate_hypertext_add_jump(archive_builder);
      } else {
        // Add character
        archive_builder_generate_text_add_character(archive_builder,enc_char);
      }
    } else {
      archive_builder_generate_text_add_character(archive_builder,dna_encoded_complement(enc_char));
    }
    // Next
    --forward_text_position;
    ++reverse_text_position;
  }
}
GEM_INLINE void archive_builder_generate_rc_hypertext(archive_builder_t* const archive_builder,const bool verbose) {
  const uint64_t num_intervals = locator_builder_get_num_intervals(archive_builder->locator);
  // Prepare ticker
  ticker_t ticker_rc;
  ticker_percentage_reset(&ticker_rc,verbose,"Generating Reverse-Complement & RC graph-links",num_intervals,20,true);
  // Traverse all reference intervals (num_base_intervals are those from the graph file; not RC)
  uint64_t i;
  for (i=0;i<num_intervals;++i) {
    // Get interval
    locator_interval_t* const locator_interval = locator_builder_get_interval(archive_builder->locator,i);
    switch (locator_interval->type) {
      case locator_interval_variant:
        continue;
        break;
      case locator_interval_unknown: {
        // Add RC interval to locator
        const int64_t rc_tag_id = -(locator_interval->tag_id);
        locator_builder_add_interval(archive_builder->locator,
            rc_tag_id,locator_interval->sequence_offset,
            locator_interval->sequence_length,0,locator_interval_unknown);
        }
        break;
      case locator_interval_regular: {
        // Generate RC of the interval & graph links
        archive_builder_generate_rc_hypertext_inteval(archive_builder,locator_interval);
        // Add RC interval to locator
        const int64_t rc_tag_id = -(locator_interval->tag_id);
        locator_builder_add_interval(archive_builder->locator,
            rc_tag_id,locator_interval->sequence_offset,
            archive_builder->parsing_state.text_interval_length,
            archive_builder->parsing_state.index_interval_length,locator_interval_regular);
        }
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
    // Add Separator
    archive_builder_generate_text_add_separator(archive_builder);
    // Tick
    ticker_update(&ticker_rc,1);
  }
  ticker_finish(&ticker_rc);
  // Add extra separator (Close text)
  archive_builder_generate_text_add_separator(archive_builder);
}
/*
 * STEP1 Archive Build :: Process MultiFASTA+GRAPH file
 *   1. MultiFASTA Read cycle
 *     1.1 Filter UIPAC-DNA bases
 *     1.2 Strip Ns
 *   2. Generate Locator
 *   3. Generate Index-Text
 *     3.1 Issue graph jumps
 *     3.2 Locate graph pending conections
 *   4. Count k-mers (Classify SA positions)
 *   5. Write (1/4) :: Header & Locator
 */
GEM_INLINE void archive_builder_process_multifasta__graph(
    archive_builder_t* const archive_builder,input_file_t* const input_multifasta,
    const bool dump_locator_intervals,const bool dump_indexed_text,const bool dump_graph_links,const bool verbose) {
  // Check number of links
  if (graph_text_builder_get_num_links(archive_builder->graph) == 0) {
    archive_builder_process_multifasta(archive_builder,
        input_multifasta,dump_locator_intervals,dump_indexed_text,verbose);
    return;
  }
  // Generate Text & Graph (Forward)
  archive_builder_generate_hypertext(archive_builder,input_multifasta,verbose);
  // Generate Text & Graph (Reverse/Complement)
  archive_builder->indexed_complement = true;
  archive_builder_generate_rc_hypertext(archive_builder,verbose);
  /*
   * DEBUG
   */
  // DEBUG locator
  locator_builder_print(gem_info_get_stream(),archive_builder->locator,dump_locator_intervals); // Locator
  // DEBUG index_text
  if (dump_indexed_text) archive_builder_dump_index_text(archive_builder,".text");
  /*
   * Write (1/3) :: Header & Locator
   */
  archive_builder_write_header(archive_builder);
  archive_builder_write_locator(archive_builder);
  // FIXME archive_builder_write_graph__jump_table(archive_builder,dump_graph_links);
}
