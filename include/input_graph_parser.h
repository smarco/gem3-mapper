/*
 * PROJECT: GEMMapper
 * FILE: input_map_graph_parser.h
 * DATE: 01/02/2014
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: TODO
 */

#ifndef INPUT_GRAPH_PARSER_H_
#define INPUT_GRAPH_PARSER_H_

#include "essentials.h"

#include "input_file.h"
#include "graph_text.h"
#include "cdna_text.h"

#include "graph_text_builder.h"

/*
 * Error handler
 */
void input_graph_link_parse_prompt_error(input_file_t* const input_graph,const error_code_t error_code);
/*
 * Graph Parser
 */
error_code_t input_graph_link_parse(
    input_file_t* const input_graph,string_t* const tag,
    string_t* const sequence_name_from,graph_link_t* const graph_link_from,
    string_t* const sequence_name_to,graph_link_t* const graph_link_to,string_t* const sequence);

/*
 * Error codes
 */
#define GRAPH_ERROR_WRONG_FILE_FORMAT 10
#define GRAPH_ERROR_PREMATURE_EOF     11

#define GRAPH_ERROR_LINK_TAG_EMPTY    20
#define GRAPH_ERROR_TAG_EMPTY         30
#define GRAPH_ERROR_STRAND_INVALID    40
#define GRAPH_ERROR_BAD_CHAR          50

/*
 * Error Messages
 */
#define GEM_ERROR_PARSE_GRAPH_WRONG_FILE_FORMAT "Parsing GRAPH error(%"PRI_input_file"). Wrong file format"
#define GEM_ERROR_PARSE_GRAPH_PREMATURE_EOF "Parsing GRAPH error(%"PRI_input_file"). Premature EOF"
#define GEM_ERROR_PARSE_GRAPH_LINK_TAG_EMPTY "Parsing GRAPH error(%"PRI_input_file"). Empty link TAG"
#define GEM_ERROR_PARSE_GRAPH_TAG_EMPTY "Parsing GRAPH error(%"PRI_input_file"). Empty sequence TAG"
#define GEM_ERROR_PARSE_GRAPH_STRAND_INVALID "Parsing GRAPH error(%"PRI_input_file"). Invalid strand"
#define GEM_ERROR_PARSE_GRAPH_BAD_CHAR "Parsing GRAPH error(%"PRI_input_file"). Invalid base/character in sequence found"

#endif /* INPUT_GRAPH_PARSER_H_ */
