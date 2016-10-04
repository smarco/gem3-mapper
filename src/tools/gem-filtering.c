/*
 * PROJECT: GEMMapper
 * FILE: gem-constructor.c
 * DATE:5/12/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Loads/test things
 */

#include "text/cdna_text.h"
#include "utils/essentials.h"
#include "utils/packed_integer_array.h"
#include "utils/priority_queue.h"
#include "utils/segmented_vector.h"
#include "utils/sparse_array_locator.h"
#include "utils/sparse_bitmap.h"
#include "utils/vector.h"
#include "utils/options_menu.h"
#include "stats/stats_vector.h"
#include "align/align_ond.h"
#include "archive/locator_builder.h"
#include "archive/archive.h"
#include "archive/archive_text.h"
#include "filtering/filtering_candidates.h"
#include "filtering/region_profile.h"
#include "io/input_parser.h"
#include "matches/match_align_dto.h"
#include "matches/match_alignment.h"
#include "neighborhood_search/nsearch_hamming.h"
#include "neighborhood_search/nsearch_levenshtein.h"

/*
 * Version
 */
#define GEM_VERSION_STRING(version) QUOTE(version)
char* const gem_version = GEM_VERSION_STRING(GEM_VERSION);

/*
 * Generic parameters
 */
typedef struct {
  char *name_input_file;
  char *name_output_file;
  char *option;
  uint64_t number;
  uint64_t param1;
  uint64_t param2;
  uint64_t param3;
  uint64_t param4;
  uint64_t param5;
  uint64_t param6;
  uint64_t param7;
  uint64_t num_threads;
  uint64_t max_memory;
  char* tmp_folder;
} gem_map_filter_args;

gem_map_filter_args parameters = {
    .name_input_file=NULL,
    .name_output_file=NULL,
    .option="",
    .number=0,
    .num_threads=1,
    .max_memory=0,
    .tmp_folder=NULL
};

/*
 * Filter testing
 */
char alpha = "ACGT";
void constructor_debruijn_build_pattern(
    char* const pattern,
    const uint64_t position,
    const uint64_t length) {
  // Check posiiton
  if (position==length) {

  } else {
    uint64_t a;
    for (a=0;a<4;++a) {
      // Set char
      pattern[position] = alpha[a];
      // Continue generation of pattern
      constructor_debruijn_build_pattern(pattern,position+1,length);
    }
  }
}
void constructor_debruijn_filter(const uint64_t length) {
  // Parameters
  char* const pattern = malloc(length*sizeof(char));
  // Generate all words in alpha
  constructor_debruijn_build_pattern(pattern,0,length);
  // Free
  free(pattern);
}
/*
 * Generic Menu
 */
void usage() {
  fprintf(stderr, "USE: ./gem-tools-examples -i input -o output \n"
                  "      Options::\n"
                  "        --input|i      <File>\n"
                  "        --output|o     <File>\n"
                  "        --select|s     <Number>\n"
                  "        --number|n     <Number>\n"
                  "        --tmp-folder   <Path>\n"
                  "        --help|h\n");
}
void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "select", required_argument, 0, 's' },
    { "number", required_argument, 0, 'n' },
    // Parameters
    { "param1", required_argument, 0, '1' },
    { "param2", required_argument, 0, '2' },
    { "param3", required_argument, 0, '3' },
    { "param4", required_argument, 0, '4' },
    { "param5", required_argument, 0, '5' },
    { "param6", required_argument, 0, '6' },
    { "param7", required_argument, 0, '7' },
    // Temporary folder
    { "tmp-folder", required_argument, 0, 100},
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:o:s:n:1:2:3:4:5:6:7:T:h",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    case 's':
      parameters.option = optarg;
      break;
    case 'n':
     parameters.number = atol(optarg);
     break;
    /* Numbered Params */
    case '1':
     parameters.param1 = atol(optarg);
     break;
    case '2':
     parameters.param2 = atol(optarg);
     break;
    case '3':
     parameters.param3 = atol(optarg);
     break;
    case '4':
     parameters.param4 = atol(optarg);
     break;
    case '5':
     parameters.param5 = atol(optarg);
     break;
    case '6':
     parameters.param6 = atol(optarg);
     break;
    case '7':
     parameters.param7 = atol(optarg);
     break;
    case '8':
      mm_set_tmp_folder(optarg);
      break;
    case 'T':
      parameters.num_threads = atol(optarg);
      break;
    case 'M': // --max-memory
      gem_cond_fatal_error_msg(input_text_parse_size(optarg,&(parameters.max_memory)),"Wrong Size '%s'",optarg);
      break;
    case 100: // --tmp-folder
      parameters.tmp_folder = optarg;
      break;
    case 'h':
      usage();
      exit(1);
    case '?': default:
      fprintf(stderr, "Option not recognized \n"); exit(1);
    }
  }
  /* System */
  if (parameters.tmp_folder!=NULL) mm_set_tmp_folder(parameters.tmp_folder);
}
int main(int argc,char** argv) {
  // GT error handler
  gem_handle_error_signals();

  // Parsing command-line options
  parse_arguments(argc,argv);

  constructor_debruijn_filter();

  return 0;
}

