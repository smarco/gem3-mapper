/*
 * PROJECT: GEMMapper
 * FILE: input_map_graph_parser.c
 * DATE: 01/02/2014
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#include "input_graph_parser.h"
#include "input_parser.h"

/*
 * Error handler
 */
GEM_INLINE void input_graph_link_parse_prompt_error(input_file_t* const input_graph,const error_code_t error_code) {
  // Display textual error msg
  switch (error_code) {
    case 0: /* No error */ break;
    case GRAPH_ERROR_PREMATURE_EOF: gem_error(PARSE_GRAPH_PREMATURE_EOF,PRI_input_file_content(input_graph)); break;
    case GRAPH_ERROR_LINK_TAG_EMPTY: gem_error(PARSE_GRAPH_LINK_TAG_EMPTY,PRI_input_file_content(input_graph)); break;
    case GRAPH_ERROR_TAG_EMPTY: gem_error(PARSE_GRAPH_TAG_EMPTY,PRI_input_file_content(input_graph)); break;
    case GRAPH_ERROR_STRAND_INVALID: gem_error(PARSE_GRAPH_STRAND_INVALID,PRI_input_file_content(input_graph)); break;
    case GRAPH_ERROR_BAD_CHAR: gem_error(PARSE_GRAPH_BAD_CHAR,PRI_input_file_content(input_graph)); break;
    default:
      gem_error(PARSE_GRAPH_WRONG_FILE_FORMAT,PRI_input_file_content(input_graph));
      break;
  }
}
/*
 * Components Parsers
 */
GEM_INLINE error_code_t input_graph_link_parse_tag(input_file_t* const input_graph,string_t* const tag) {
  // Clear string
  string_clear(tag);
  // Read Tag until TAB
  input_file_parse_field(input_graph,TAB,tag);
  return 0;
}
GEM_INLINE error_code_t input_graph_link_parse_link(
    input_file_t* const input_graph,
    string_t* const sequence_name,graph_link_t* const graph_link) {
  /*
   * Graph Link Parser
   *   Chr1  +  10090
   */
  error_code_t error_code;
  // Read sequence name
  string_clear(sequence_name);
  input_file_parse_field(input_graph,TAB,sequence_name);
  if (gem_expect_false(string_get_length(sequence_name)==0)) return GRAPH_ERROR_TAG_EMPTY;
  // Skip SEPARATORS
  if ((error_code=input_file_parse_skip_separators(input_graph))) return error_code;
  if (input_file_parse_is_eol(input_graph)) return GRAPH_ERROR_PREMATURE_EOF;
  // Read Strand
  const char current_char = input_file_get_current_char(input_graph);
  switch (current_char) { // Encode strand in the @seq_id
    case '+': graph_link->tag_id =  1; break;
    case '-': graph_link->tag_id = -1; break;
    default:
      return GRAPH_ERROR_STRAND_INVALID;
      break;
  }
  input_file_next_char(input_graph); // Next
  // Skip SEPARATORS
  if ((error_code=input_file_parse_skip_separators(input_graph))) return error_code;
  if (input_file_parse_is_eol(input_graph)) return GRAPH_ERROR_PREMATURE_EOF;
  // Read Position
  int64_t number = 0;
  if ((error_code=input_file_parse_integer(input_graph,&number))) return error_code;
  graph_link->position_text = number;
  return 0;
}
GEM_INLINE error_code_t input_graph_link_parse_sequence(input_file_t* const input_graph,string_t* const sequence) {
  // Clear
  string_clear(sequence);
  // Read loop
  do {
    // Read character
    const char current_char = input_file_get_current_char(input_graph);
    if (IS_ANY_EOL(current_char) || current_char==TAB || current_char==SPACE) break;
    // Store Base
    if (!is_iupac_code(current_char)) return GRAPH_ERROR_BAD_CHAR;
    string_append_char(sequence,current_char);
    // Next
  } while (input_file_next_char(input_graph));
  // Append EOS
  string_append_eos(sequence);
  // Return
  return 0;
}
/*
 * Graph Parser
 *   rs373063181     1       +       28591   1       +       28591   GGGGG
 *   rs373063181     1       -       28591   1       -       28591   CCCCC
 *   rs369923773     1       +       52184   1       +       52186
 *   rs369923773     1       -       52186   1       -       52184
 *   rs12123356      1       +       724701  1       +       724702  A
 *   rs12123356      1       -       724702  1       -       724701  T
 *   rs12123356      1       +       724701  1       +       724702  C
 *   rs12123356      1       -       724702  1       -       724701  G
 *   rs12123356      1       +       724701  1       +       724702  T
 *   rs12123356      1       -       724702  1       -       724701  A
 */
GEM_INLINE error_code_t input_graph_link_parse_(
    input_file_t* const input_graph,string_t* const tag,
    string_t* const sequence_name_from,graph_link_t* const graph_link_from,
    string_t* const sequence_name_to,graph_link_t* const graph_link_to,string_t* const sequence) {
  error_code_t error_code;
  // Check EOF
  input_file_check_buffer(input_graph);
  if (gem_expect_false(input_file_eof(input_graph))) return INPUT_STATUS_EOF;
  // Read tag
  if ((error_code=input_graph_link_parse_tag(input_graph,tag))) return error_code;
  if (gem_expect_false(string_get_length(tag)==0)) return GRAPH_ERROR_LINK_TAG_EMPTY;
  // Skip SEPARATORS
  if ((error_code=input_file_parse_skip_separators(input_graph))) return error_code;
  if (input_file_parse_is_eol(input_graph)) return GRAPH_ERROR_PREMATURE_EOF;
  // Read Source Graph Link
  if ((error_code=input_graph_link_parse_link(input_graph,sequence_name_from,graph_link_from))) return error_code;
  if (input_file_parse_is_eol(input_graph)) return GRAPH_ERROR_PREMATURE_EOF;
  // Skip SEPARATORS
  if ((error_code=input_file_parse_skip_separators(input_graph))) return error_code;
  if (input_file_parse_is_eol(input_graph)) return GRAPH_ERROR_PREMATURE_EOF;
  // Read Destiny Graph Link
  if ((error_code=input_graph_link_parse_link(input_graph,sequence_name_to,graph_link_to))) return error_code;
  if (!input_file_parse_is_eol(input_graph)) {
    // Skip SEPARATORS
    if ((error_code=input_file_parse_skip_separators(input_graph))) return error_code;
    if (!input_file_parse_is_eol(input_graph)) { // Read Sequence (text)
      if ((error_code=input_graph_link_parse_sequence(input_graph,sequence))) return error_code;
    }
  }
  return 0;
}
GEM_INLINE error_code_t input_graph_link_parse(
    input_file_t* const input_graph,string_t* const tag,
    string_t* const sequence_name_from,graph_link_t* const graph_link_from,
    string_t* const sequence_name_to,graph_link_t* const graph_link_to,string_t* const sequence) {
  error_code_t error_code = 0;
  // Clean
  string_clear(sequence);
  /*
   * Parse the graph link
   */
  error_code = input_graph_link_parse_(input_graph,tag,
      sequence_name_from,graph_link_from,
      sequence_name_to,graph_link_to,sequence);
  // Skip the rest of the line
  input_file_parse_skip_line(input_graph);
  // Return
  return error_code;
}
