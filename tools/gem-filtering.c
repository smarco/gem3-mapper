/*
 * PROJECT: GEMMapper
 * FILE: gem-filtering.c
 * DATE: 1/10/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include <pthread.h>
#include "gem_core.h"

/*
 * GEM-filtering Operators
 */
typedef enum {
  filtering_generate_dataset,
  filtering_counting,
  filtering_qsamples,
  filtering_qgrams
} filtering_operators_t;
/*
 * GEM-filtering Parameters
 */
typedef struct {
  /* I/O */
  char* input_file_name;
  /* Filtering */
  filtering_operators_t filtering_operators;
  /* Dataset Generation */
  uint64_t error;
  /* System */
  uint64_t num_threads;
  uint64_t max_memory;
  char* tmp_folder;
  /* Miscellaneous */
  bool verbose;
} gem_filtering_parameters;
// Defaults
gem_filtering_parameters parameters = {
  /* I/O */
  .input_file_name=NULL,
  /* Filtering */
  /* Dataset Generation */
  /* System */
  .num_threads=1,
  .max_memory=0,
  .tmp_folder=NULL,
  /* Miscellaneous */
  .verbose=true,
};
/*
 * Filtering functions
 */
GEM_INLINE uint64_t filtering_solve_dp(
    const char* const pattern,const uint64_t pattern_length,
    const uint64_t* const dp,uint64_t* const dp_next,
    const uint64_t pos,const char character) {
  uint64_t i, min = UINT64_MAX;
  const uint64_t sub = (pos>0 ? pos-1 : 0) + (character == pattern[pos] ? 0 : 1);
  const uint64_t ins = dp[0] + 1;
  const uint64_t del = pos + 1;
  dp_next[0] = MIN(sub,MIN(ins,del));
  for (i=1;i<pattern_length;++i) {
    const uint64_t sub = dp[i-1] + (character == pattern[pos] ? 0 : 1);
    const uint64_t ins = dp[i] + 1;
    const uint64_t del = dp_next[i-1] + 1;
    dp_next[i] = MIN(sub,MIN(ins,del));
    min = MIN(min,dp_next[i]);
  }
  return min;
}
GEM_INLINE void filtering_generate_neighborhood(
    const char* const pattern,const uint64_t pattern_length,
    char* const text,const uint64_t max_text_length,
    const uint64_t pos,const uint64_t* const dp,
    const uint64_t error,const bool expand) {
  // Allocate DP
  uint64_t min_val;
  uint64_t* dp_next = calloc(pattern_length,sizeof(uint64_t));
  char* alphabet = "ACGT";
  uint64_t i;
  for (i=0;i<4;++i) {
    // Solve for 'c'
    text[pos] = alphabet[i];
    min_val = filtering_solve_dp(pattern,pattern_length,dp,dp_next,pos,alphabet[i]);
    if (dp_next[pattern_length-1]==error) {
      text[pos+1] = '\0';
      fprintf(stdout,">%s\n",text);
//    } else if (min_val==error && pos+1<max_text_length) {
//      filtering_generate_neighborhood_exact(pattern,pattern_length,text,max_text_length,pos+1,dp_next,error,true);
    } else if (min_val<=error && pos+1<max_text_length) {
      filtering_generate_neighborhood(pattern,pattern_length,text,max_text_length,pos+1,dp_next,error,true);
    }
  }
  free(dp_next);
}
GEM_INLINE void filtering_test_generate_dataset() {
  char* base_read = NULL;
  uint64_t len = 0, read_length, i;
  // Read all
  while ((read_length = getline(&base_read,&len,stdin)) != -1) {
    // Allocate DP
    const uint64_t max_text_length = read_length+parameters.error;
    char* const text = calloc(max_text_length,sizeof(char));
    memset(text,0,max_text_length*sizeof(char));
    uint64_t* dp = calloc(read_length,sizeof(uint64_t));
    // Print
    fprintf(stdout,"@%.*s",(int)read_length,base_read);
    // Init
    for (i=0;i<read_length;++i) dp[i] = i;
    // Generate neighborhood
    filtering_generate_neighborhood(
        base_read,read_length,text,max_text_length,0,dp,parameters.error,parameters.error>0);
    // Free
    free(dp);
    free(text);
  }
  free(base_read);
}
GEM_INLINE void filtering_test_counting() {

}
//GEM_INLINE void filtering_test_counting() {
//
//}
//GEM_INLINE void filtering_test_counting() {
//
//}
/*
 * GEM-filtering options Menu
 */
option_t gem_filtering_options[] = {
  /* I/O */
  { 'i', "input", REQUIRED, TYPE_STRING, 2 , true, "<input_file>" , "(Multi-FASTA, default=stdin)" },
  { 'o', "output", REQUIRED, TYPE_STRING, 2 , true, "<output_prefix>" , "" },
  /* Filtering */
  { 300, "filter", REQUIRED, TYPE_STRING, 3, true, "'generate'|'count'|'q-gram'|'q-sample'" , "" },
  /* Dataset Generation */
  { 'e', "error", REQUIRED, TYPE_INT, 4, true, "<number>" , "(default=0)" },
  /* System */
  { 't', "threads", REQUIRED, TYPE_INT, 5 , true, "<number>" , "(default=1)" },
  { 500, "max-memory", REQUIRED, TYPE_STRING, 5 , true, "<maximum-memory>" , "(Eg 2GB)" },
  { 501, "tmp-folder", REQUIRED, TYPE_STRING, 5 , true, "<temporal_dir_path>" , "(/tmp/)" },
  /* Miscellaneous */
  { 'v', "verbose", NO_ARGUMENT, TYPE_NONE, 7 ,true, "", "" },
  { 'q', "quiet", NO_ARGUMENT, TYPE_NONE, 7 ,true, "", "" },
  { 'h', "help", NO_ARGUMENT, TYPE_NONE, 7 , true, "" , "(print usage)" },
  { 'H', "help", NO_ARGUMENT, TYPE_NONE, 7 , false, "" , "(print usage + extras)" },
  /* Extras/Develop */
  {  0, "", 0, 0, 0, false, "", ""}
};
char* gem_filtering_groups[] = {
  /* 0 */ "Null",
  /* 1 */ "Unclassified",
  /* 2 */ "I/O",
  /* 3 */ "Filtering",
  /* 4 */ "Dataset Generation",
  /* 5 */ "System",
  /* 6 */ "Miscellaneous",
  /* 7 */ "Extras"
};
void usage(const bool print_inactive) {
  fprintf(stderr, "USAGE: ./gem-filtering [ARGS]...\n");
  options_fprint_menu(stderr,gem_filtering_options,gem_filtering_groups,true,print_inactive);
}
void parse_arguments(int argc,char** argv) {
  struct option* getopt_options = options_adaptor_getopt(gem_filtering_options);
  string_t* const getopt_short_string = options_adaptor_getopt_short(gem_filtering_options);
  char* const getopt_short = string_get_buffer(getopt_short_string);
  int option, option_index;
  while (true) {
    // Get option &  Select case
    if ((option=getopt_long(argc,argv,getopt_short,getopt_options,&option_index))==-1) break;
    switch (option) {
    /* I/O */
    case 'i': // --input
      parameters.input_file_name = optarg;
      break;
    case 'o': // --output
      break;
    /* Filtering */
    case 300: // --filter in {'generate'|'count'|'q-gram'|'q-sample'}
      if (gem_streq(optarg,"generate")) {
        parameters.filtering_operators = filtering_generate_dataset;
      } else if (gem_streq(optarg,"q-gram")) {
        parameters.filtering_operators = filtering_counting;
      } else if (gem_streq(optarg,"q-sample")) {
        parameters.filtering_operators = filtering_qsamples;
      } else if (gem_streq(optarg,"q-gram")) {
        parameters.filtering_operators = filtering_qgrams;
      } else {
        GEM_INVALID_CASE();
      }
      break;
    /* Dataset Generation */
    case 'e':
      parameters.error = atol(optarg);
      break;
    /* System */
    case 't': // --threads
      parameters.num_threads = atol(optarg);
      break;
    case 500: // --max-memory
      gem_cond_fatal_error(input_text_parse_size(optarg,&(parameters.max_memory)),PARSING_SIZE,"--max-memory",optarg);
      break;
    case 501: // --tmp-folder
      parameters.tmp_folder = optarg;
      break;
    /* Miscellaneous */
    case 'h':
      usage(false);
      exit(1);
    case 'H':
      usage(true);
      exit(1);
    case 'v':
      parameters.verbose = true;
      break;
    case 'q':
      parameters.verbose = false;
      break;
    /* Extras/Develop */
    case '?':
    default:
      gem_fatal_error_msg("Option not recognized");
    }
  }
  /*
   * Parameters Check
   */
  // System
  if (parameters.max_memory==0) {
    parameters.max_memory = mm_get_available_mem();
  }
  /*
   * Free
   */
  string_destroy(getopt_short_string);
  mm_free(getopt_short_string);
  mm_free(getopt_options);
}
/*
 * Main()
 */
int main(int argc,char** argv) {
  // Global Filtering Timer
  gem_timer_t gem_filtering_timer;

  // Parsing command-line options
  parse_arguments(argc,argv);

  // GEM Runtime setup
  gem_runtime_init(parameters.num_threads,parameters.max_memory,parameters.tmp_folder);
  TIMER_RESTART(&gem_filtering_timer); // Start global timer

  /*
   * Select proper operator
   */
  switch (parameters.filtering_operators) {
    case filtering_generate_dataset:
      filtering_test_generate_dataset();
      break;
    case filtering_counting:
      filtering_test_counting();
      break;
    case filtering_qsamples:
      filtering_test_counting();
      break;
    case filtering_qgrams:
      filtering_test_counting();
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }

  // Display end banner
  TIMER_STOP(&gem_filtering_timer);
  if (parameters.verbose) {
    tfprintf(stderr,"[... was successfully done in %2.3f min.]\n",TIMER_GET_TOTAL_M(&gem_filtering_timer));
  }

  // Free
  gem_runtime_destroy();

  return 0;
}

