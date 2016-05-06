/*
 * PROJECT: GEMMapper
 * FILE: gem-indexer.c
 * DATE: 1/10/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Encodes a Multi-FASTA file (.fa) into a GEM index
 */

#include "utils/essentials.h"
#include "utils/options_menu.h"
#include "io/input_parser.h"
#include "archive/archive_builder.h"
#include "archive/archive_builder_text.h"
#include "archive/archive_builder_index.h"
#include "system/profiler_timer.h"

/*
 * GEM-indexer Debug
 */
#define GEM_INDEXER_DEBUG_DUMP_EXPLICIT_CHECKED_SA  false
#define GEM_INDEXER_DEBUG_DUMP_CHECKED_BWT          false

/*
 * GEM-mapper Error Handling
 */
#define gem_indexer_error_msg(error_msg,args...) \
  fprintf(stderr,"GEM-Mapper error:\n> "error_msg"\n",##args); \
  exit(1)
#define gem_indexer_cond_error_msg(condition,error_msg,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gem_indexer_error_msg(error_msg,##args); \
    } \
  } while (0)

/*
 * GEM-indexer Parameters
 */
typedef struct {
  /* I/O */
  char* input_multifasta_file_name;
  char* input_graph_file_name;
  char* output_index_file_name;
  char* output_index_file_name_prefix;
  uint64_t ns_threshold;
  indexed_complement_t index_complement;
  uint64_t complement_size_threshold;
  bool indexed_reverse_text;
  /* Index */
  sampling_rate_t sa_sampling_rate;
  sampling_rate_t text_sampling_rate;
  bool run_length_index;
  bool bisulfite_index;
  bool gpu_index;
  /* Debug */
  bool dump_locator_intervals;
  bool dump_indexed_text;
  bool dump_explicit_sa;
  bool dump_bwt;
  bool dump_run_length_text;
  bool dump_graph_links;
  bool check_index;
  bool dev_generate_only_bwt;
  /* System */
  uint64_t num_threads;
  uint64_t max_memory;
  char* tmp_folder;
  /* Miscellaneous */
  bool info_file_name_provided;
  char* info_file_name;
  bool verbose;
} indexer_parameters_t;
// Defaults
void indexer_parameters_set_defaults(indexer_parameters_t* const parameters) {
  /* I/O */
  parameters->input_multifasta_file_name=NULL;
  parameters->input_graph_file_name=NULL;
  parameters->output_index_file_name=NULL;
  parameters->output_index_file_name_prefix=NULL;
  /* MultiFasta Processing */
  parameters->ns_threshold=50;
  parameters->index_complement=index_complement_yes;
  parameters->complement_size_threshold=BUFFER_SIZE_8G;
  parameters->indexed_reverse_text = false;
  /* Index */
  parameters->sa_sampling_rate=SAMPLING_RATE_8;
  parameters->text_sampling_rate=SAMPLING_RATE_4;
  parameters->run_length_index = false;
  parameters->bisulfite_index = false;
#ifdef HAVE_CUDA
  parameters->gpu_index = true;
#else
  parameters->gpu_index = false;
#endif
  /* Debug */
  parameters->dump_locator_intervals=true;
  parameters->dump_indexed_text=false;
  parameters->dump_explicit_sa=false;
  parameters->dump_bwt=false;
  parameters->dump_run_length_text=false;
  parameters->dump_graph_links=false;
  parameters->check_index=false;
  parameters->dev_generate_only_bwt=false;
  /* System */
  parameters->num_threads=system_get_num_processors();
  parameters->max_memory=0;
  parameters->tmp_folder=NULL;
  /* Miscellaneous */
  parameters->info_file_name_provided=false;
  parameters->info_file_name=NULL;
  parameters->verbose=true;
}
/*
 * Debug
 *   Simple function to build the SA of a text using qsort
 *   and dump it to a file (to check against the indexer)
 */
char* debug_char_text;
uint64_t debug_text_length;
uint64_t* debug_SA = NULL;
#define INDEXER_GET_BWT_CHAR_FROM_SA(text,text_length,sa_pos) ((sa_pos==0) ? text[text_length-1] : text[sa_pos-1])
int indexer_debug_suffix_cmp(const uint64_t* const a,const uint64_t* const b) {
  // Same SA-pos
  if (gem_expect_false(*a==*b)) return 0;
  // Return cmp result
  const bool a_lt_b = (*a < *b);
  const uint64_t cmp_length = (a_lt_b) ? debug_text_length-*b : debug_text_length-*a;
  const int cmp = memcmp(debug_char_text+*a,debug_char_text+*b,cmp_length);
  if (cmp) return cmp;
  if (a_lt_b) {
    uint64_t offset_a = *a+cmp_length;
    uint64_t offset_b = 0;
    return indexer_debug_suffix_cmp(&offset_a,&offset_b);
  } else {
    uint64_t offset_a = 0;
    uint64_t offset_b = *b+cmp_length;
    return indexer_debug_suffix_cmp(&offset_a,&offset_b);
  }
}
void indexer_debug_generate_sa(dna_text_t* const index_text) {
  // Allocate memory for SA
  debug_char_text = (char*) dna_text_get_text(index_text);
  debug_text_length = dna_text_get_length(index_text);
  debug_SA = mm_malloc(sizeof(uint64_t)*(debug_text_length));
  uint64_t i;
  for (i=0;i<debug_text_length;++i) debug_SA[i] = i;
  // Build SA
  qsort(debug_SA,debug_text_length,sizeof(uint64_t),(int (*)(const void *,const void *))indexer_debug_suffix_cmp);
}
void indexer_debug_check_sa(char* const file_name,dna_text_t* const index_text) {
  // Open file
  FILE* const bwt_file = fopen(file_name,"w");
  // Generate SA
  if (debug_SA==NULL) indexer_debug_generate_sa(index_text);
  // Dump full SA
  uint64_t i, t, num_printed_chars;
  for (i=0;i<debug_text_length;++i) {
    num_printed_chars=0;
    fprintf(bwt_file,"Suffix=%011"PRIu64"\t\t",debug_SA[i]);
    for (t=debug_SA[i];t<debug_text_length;++t) {
      fprintf(bwt_file,"%c",dna_decode(debug_char_text[t]));
      if (++num_printed_chars >= ARCHIVE_BUILDER_DEBUG_SUFFIX_LENGTH) break;
    }
    for (t=0;t<debug_SA[i] && num_printed_chars<ARCHIVE_BUILDER_DEBUG_SUFFIX_LENGTH;++t) {
      fprintf(bwt_file,"%c",dna_decode(debug_char_text[t]));
      if (++num_printed_chars >= ARCHIVE_BUILDER_DEBUG_SUFFIX_LENGTH) break;
    }
    fprintf(bwt_file,"\n");
  }
  // Close
  fclose(bwt_file);
}
void indexer_debug_check_bwt(char* const file_name,dna_text_t* const index_text) {
  // Open file
  FILE* const bwt_file = fopen(file_name,"w");
  // Generate SA
  if (debug_SA==NULL) indexer_debug_generate_sa(index_text);
  // Dump BWT to check file
  uint64_t i;
  for (i=0;i<debug_text_length;) {
    fprintf(bwt_file,"%c",dna_decode(INDEXER_GET_BWT_CHAR_FROM_SA(debug_char_text,debug_text_length,debug_SA[i])));
    if (++i%80==0) fprintf(bwt_file,"\n"); // Print column-wise
  }
  if (i%80!=0) fprintf(bwt_file,"\n");
  // Free
  fclose(bwt_file);
}
/*
 * GEM-Indexer Build Archive
 */
void indexer_process_multifasta(archive_builder_t* const archive_builder,indexer_parameters_t* const parameters) {
  // Process input MultiFASTA
  input_file_t* const input_multifasta = input_file_open(parameters->input_multifasta_file_name,BUFFER_SIZE_32M,false);
  // Process MFASTA
  archive_builder_text_process(archive_builder,input_multifasta,parameters->verbose);
  // RL-Index
  if (parameters->run_length_index) {
    archive_builder_text_apply_run_length(archive_builder,parameters->verbose);
    // TODO if (parameters->dump_run_length_text) archive_builder_text_dump(archive_builder,".text.rl");
  }
  // DEBUG
  locator_builder_print(gem_info_get_stream(),archive_builder->locator, parameters->dump_locator_intervals);
  if (parameters->dump_indexed_text) archive_builder_text_dump(archive_builder,".text");
  // Write Metadata
  archive_builder_write_header(archive_builder);
  archive_builder_write_locator(archive_builder);
  locator_builder_delete(archive_builder->locator); // Free Locator
  // Close MultiFASTA
  input_file_close(input_multifasta);
}
void indexer_generate_bwt(archive_builder_t* const archive_builder,indexer_parameters_t* const parameters) {
  // VERBOSE
  tfprintf(gem_log_get_stream(),"[Generating BWT Forward-Text]\n");
  // Build BWT
  archive_builder_index_build_bwt(
      archive_builder,parameters->gpu_index,parameters->dump_bwt,
      parameters->dump_explicit_sa,parameters->verbose);
  // DEBUG. Print Explicit Checked-SA => (SApos,SA[SApos...SApos+SAFixLength])
  gem_cond_debug_block(GEM_INDEXER_DEBUG_DUMP_EXPLICIT_CHECKED_SA) {
    indexer_debug_check_sa(gem_strcat(parameters->output_index_file_name_prefix,".check.sa"),archive_builder->enc_text);
  }
  // DEBUG. Print Checked-BWT (Generate SA-BWT (to compare with))
  gem_cond_debug_block(GEM_INDEXER_DEBUG_DUMP_CHECKED_BWT) {
    indexer_debug_check_bwt(gem_strcat(parameters->output_index_file_name_prefix,".check.bwt"),archive_builder->enc_text);
    if (debug_SA) mm_free(debug_SA); // Free
  }
  // DEBUG. Skip the FM-index generation
  if (parameters->dev_generate_only_bwt) exit(0);
}
void indexer_generate_bwt_reverse(archive_builder_t* const archive_builder,indexer_parameters_t* const parameters) {
  // Build BWT of the reverse text
  tfprintf(gem_log_get_stream(),"[Generating BWT Reverse-Text]\n");
  archive_builder_index_build_bwt_reverse(archive_builder,parameters->dump_indexed_text,
      parameters->dump_bwt,parameters->dump_explicit_sa,parameters->verbose);
}
void indexer_write_index(archive_builder_t* const archive_builder,indexer_parameters_t* const parameters) {
  // Write Text & FM-Index
  archive_builder_write_index(archive_builder,parameters->gpu_index,parameters->check_index,parameters->verbose);
}
void indexer_write_index_reverse(archive_builder_t* const archive_builder,indexer_parameters_t* const parameters) {
  // Write FM-Index Reverse
  archive_builder_write_index_reverse(archive_builder,parameters->check_index,parameters->verbose);
}
void indexer_cleanup(archive_builder_t* const archive_builder,indexer_parameters_t* const parameters) {
  // Archive Builder
  archive_builder_delete(archive_builder);
  // Output file name
  mm_free(parameters->output_index_file_name);
  // Info File Name
  if (!parameters->info_file_name_provided)  {
    mm_free(parameters->info_file_name);
  }
}
/*
 * GEM-Indexer options Menu
 */
// TODO { 'g', "graph", REQUIRED, TYPE_STRING, 2 , VISIBILITY_USER, "<graph_file>" , "(GRAPH-File)" },
// TODO { 'c', "index-colorspace", NO_ARGUMENT, TYPE_NONE, 2 , true, "" , "(default=false)" },
// TODO { 200, "index-complement", OPTIONAL, TYPE_STRING, 2 , true, "" , "(default=false, simulated at mapping)" },
// TODO { 201, "complement-size-threshold", REQUIRED, TYPE_INT, 2 , true, "<integer>" , "(default=8GB)" },
// TODO { 202, "generate-reverse", NO_ARGUMENT, TYPE_NONE, 3 , true, "", "(default=true)"},
// TODO { 405, "dump-graph-links", OPTIONAL, TYPE_NONE, 4 , true, "" , "" },
option_t gem_indexer_options[] = {
  /* I/O */
  { 'i', "input", REQUIRED, TYPE_STRING, 2 , VISIBILITY_USER, "<input_file>" , "(Multi-FASTA)" },
  { 'o', "output", REQUIRED, TYPE_STRING, 2 , VISIBILITY_USER, "<output_prefix>" , "" },
  { 'N', "strip-unknown-bases-threshold", REQUIRED, TYPE_INT, 2 , VISIBILITY_ADVANCED, "'disable'|<integer>" , "(default=50)" },
  /* Index */
  { 's', "text-sampling-rate", REQUIRED, TYPE_INT, 3 , VISIBILITY_ADVANCED, "<sampling_rate>" , "(default=4)" },
  { 'S', "SA-sampling-rate", REQUIRED, TYPE_INT, 3 , VISIBILITY_ADVANCED, "<sampling_rate>" , "(default=8)" },
  { 'r', "run-length-index", OPTIONAL, TYPE_NONE, 3 , VISIBILITY_USER, "" , "(default=false)" },
  { 'b', "bisulfite-index", OPTIONAL, TYPE_NONE, 3 , VISIBILITY_USER, "" , "(default=false)" },
#ifdef HAVE_CUDA
  { 300, "gpu-index", OPTIONAL, TYPE_NONE, 3 , VISIBILITY_ADVANCED, "" , "(default=true)" },
#endif
  { 301, "index-backward-text", OPTIONAL, TYPE_NONE, 3 , VISIBILITY_ADVANCED, "" , "(default=false)" },
  // TODO { 300, "autotune-index-size", REQUIRED, TYPE_STRING, 5 , VISIBILITY_USER, "<index-size>" , "(Eg 2GB)" },
  /* Debug */
  { 400, "dump-locator-intervals", OPTIONAL, TYPE_NONE, 4 , VISIBILITY_ADVANCED, "" , "" },
  { 401, "dump-indexed-text", OPTIONAL, TYPE_NONE, 4 , VISIBILITY_ADVANCED, "" , "" },
  { 402, "dump-explicit-sa", OPTIONAL, TYPE_NONE, 4 , VISIBILITY_ADVANCED, "" , "" },
  { 403, "dump-bwt", OPTIONAL, TYPE_NONE, 4 , VISIBILITY_ADVANCED, "" , "" },
  { 404, "dump-run-length-text", OPTIONAL, TYPE_NONE, 4 , VISIBILITY_ADVANCED, "" , "" },
  { 405, "debug", NO_ARGUMENT, TYPE_NONE, 4 , VISIBILITY_ADVANCED, "" , "" },
  //{ 406, "check-index", OPTIONAL, TYPE_NONE, 4 , VISIBILITY_ADVANCED, "", "(default=false)"},
  { 407, "bwt", NO_ARGUMENT, TYPE_NONE, 4 ,VISIBILITY_ADVANCED, "", "(only generate BWT for benchmarking)" },
  /* System */
  { 't', "threads", REQUIRED, TYPE_INT, 5 , VISIBILITY_USER, "<number>" , "(default=#cores)" },
  //{ 500, "max-memory", REQUIRED, TYPE_STRING, 5 , VISIBILITY_ADVANCED, "<maximum-memory>" , "(Eg 2GB)" },
  { 501, "tmp-folder", REQUIRED, TYPE_STRING, 5 , VISIBILITY_ADVANCED, "<temporal_dir_path>" , "(/tmp/)" },
  /* Miscellaneous */
  { 'v', "verbose", NO_ARGUMENT, TYPE_NONE, 6 ,VISIBILITY_USER, "", "" },
  { 'q', "quiet", NO_ARGUMENT, TYPE_NONE, 6 ,VISIBILITY_USER, "", "" },
  { 600, "info-file", REQUIRED, TYPE_STRING, 6 ,VISIBILITY_ADVANCED, "<info_file_path>", "" },
  { 'h', "help", OPTIONAL, TYPE_NONE, 6 , VISIBILITY_USER, "" , "(print usage)" },
  {  0, "", 0, 0, 0, false, "", ""}
};
char* gem_indexer_groups[] = {
  /* 0 */ "Null",
  /* 1 */ "Unclassified",
  /* 2 */ "I/O",
  /* 3 */ "Index",
  /* 4 */ "Debug",
  /* 5 */ "System",
  /* 6 */ "Miscellaneous"
};
void usage(const option_visibility_t visibility_level) {
  fprintf(stderr, "USAGE: ./gem-indexer [ARGS]...\n");
  options_fprint_menu(stderr,gem_indexer_options,gem_indexer_groups,true,visibility_level);
}
void parse_arguments(int argc,char** argv,indexer_parameters_t* const parameters) {
  struct option* getopt_options = options_adaptor_getopt(gem_indexer_options);
  string_t* const getopt_short_string = options_adaptor_getopt_short(gem_indexer_options);
  char* const getopt_short = string_get_buffer(getopt_short_string);
  int option, option_index;
  while (true) {
    // Get option &  Select case
    if ((option=getopt_long(argc,argv,getopt_short,getopt_options,&option_index))==-1) break;
    switch (option) {
    /* I/O */
    case 'i': // --input
      parameters->input_multifasta_file_name = optarg;
      break;
    case 'o': // --output
      parameters->output_index_file_name = optarg;
      break;
    case 'N': // --strip-unknown-bases-threshold
      if (gem_strcaseeq(optarg,"disable")) {
        parameters->ns_threshold = UINT64_MAX;
      } else {
        parameters->ns_threshold = atol(optarg);
      }
      break;
    /* Index */
    case 's': { // --text_sampling-rate
      const uint64_t sampling = atol(optarg);
      switch (sampling) {
        case 1:   parameters->text_sampling_rate=SAMPLING_RATE_1; break;
        case 2:   parameters->text_sampling_rate=SAMPLING_RATE_2; break;
        case 4:   parameters->text_sampling_rate=SAMPLING_RATE_4; break;
        case 8:   parameters->text_sampling_rate=SAMPLING_RATE_8; break;
        case 16:  parameters->text_sampling_rate=SAMPLING_RATE_16; break;
        case 32:  parameters->text_sampling_rate=SAMPLING_RATE_32; break;
        case 64:  parameters->text_sampling_rate=SAMPLING_RATE_64; break;
        case 128: parameters->text_sampling_rate=SAMPLING_RATE_128; break;
        case 256: parameters->text_sampling_rate=SAMPLING_RATE_256; break;
        default:
          gem_indexer_error_msg("Sampling rate argument not valid (--text_sampling-rate)");
          break;
      }
    }
    break;
    case 'S': { // --SA-sampling-rate
      const uint64_t sampling = atol(optarg);
      switch (sampling) {
        case 1:   parameters->sa_sampling_rate=SAMPLING_RATE_1; break;
        case 2:   parameters->sa_sampling_rate=SAMPLING_RATE_2; break;
        case 4:   parameters->sa_sampling_rate=SAMPLING_RATE_4; break;
        case 8:   parameters->sa_sampling_rate=SAMPLING_RATE_8; break;
        case 16:  parameters->sa_sampling_rate=SAMPLING_RATE_16; break;
        case 32:  parameters->sa_sampling_rate=SAMPLING_RATE_32; break;
        case 64:  parameters->sa_sampling_rate=SAMPLING_RATE_64; break;
        case 128: parameters->sa_sampling_rate=SAMPLING_RATE_128; break;
        case 256: parameters->sa_sampling_rate=SAMPLING_RATE_256; break;
        default:
          gem_indexer_error_msg("Sampling rate argument not valid (--SA-sampling-rate)");
          break;
      }
    }
    break;
    case 'r': // --run-length-index
      parameters->run_length_index = (optarg) ? options_parse_bool(optarg) : true;
      break;
    case 'b': // --bisulfite-index
      parameters->bisulfite_index = (optarg) ? options_parse_bool(optarg) : true;
      break;
#ifdef HAVE_CUDA
    case 300: // --gpu-index
      parameters->gpu_index = (optarg) ? options_parse_bool(optarg) : true;
      break;
#endif
    case 301: // --index-backward-text
      parameters->indexed_reverse_text = (optarg) ? options_parse_bool(optarg) : true;
      break;
    /* Debug/Temporal */
    case 400: // --dump-locator-intervals
      parameters->dump_locator_intervals = (optarg) ? options_parse_bool(optarg) : true;
      break;
    case 401: // --dump-indexed-text
      parameters->dump_indexed_text = (optarg) ? options_parse_bool(optarg) : true;
      break;
    case 402: // --dump-explicit-sa
      parameters->dump_explicit_sa = (optarg) ? options_parse_bool(optarg) : true;
      break;
    case 403: // --dump-bwt
      parameters->dump_bwt = (optarg) ? options_parse_bool(optarg) : true;
      break;
    case 404:
      parameters->dump_run_length_text = (optarg) ? options_parse_bool(optarg) : true;
      break;
//    case 405: // --dump-graph-links
//      parameters->dump_graph_links = (optarg) ? options_parse_bool(optarg) : true;
//      break;
    case 405: // --debug
      parameters->dump_locator_intervals = true;
      parameters->dump_indexed_text = true;
      parameters->dump_explicit_sa = true;
      parameters->dump_bwt = true;
      parameters->dump_graph_links = true;
      parameters->dump_run_length_text = true;
      parameters->verbose = true;
      break;
    case 406: // --check-index
      parameters->check_index = (optarg) ? options_parse_bool(optarg) : true;
      break;
    case 407: // --bwt
      parameters->dev_generate_only_bwt = true;
      parameters->dump_bwt = true;
      parameters->verbose = true;
      break;
    /* System */
    case 't': // --threads
      parameters->num_threads = atol(optarg);
      break;
    case 500: // --max-memory
      gem_indexer_cond_error_msg(input_text_parse_size(optarg,&(parameters->max_memory)),
          "Error parsing --max-memory. '%s' not a valid size (Eg. 2GB)",optarg);
      break;
    case 501: // --tmp-folder
      parameters->tmp_folder = optarg;
      break;
    /* Miscellaneous */
    case 'v':
      parameters->verbose = true;
      break;
    case 'q':
      parameters->verbose = false;
      break;
    case 600: // --info-file
      parameters->info_file_name_provided = true;
      parameters->info_file_name = optarg;
      break;
    case 'h':
      if (optarg==NULL || gem_strcaseeq(optarg,"user")) {
        usage(VISIBILITY_USER);
      } else if (gem_strcaseeq(optarg,"advanced")) {
        usage(VISIBILITY_ADVANCED);
      } else if (gem_strcaseeq(optarg,"developer")) {
        usage(VISIBILITY_DEVELOPER);
      } else {
        gem_indexer_error_msg("Help argument not valid {'user','advanced'}");
      }
      exit(0);
    /* */
    case '?':
    default:
      gem_indexer_error_msg("Option not recognized");
    }
  }
  /*
   * Parameters Check
   */
//  // Index type incompatibility list
//  if (parameters->input_graph_file_name!=NULL) {
//    gem_indexer_cond_error_msg_msg(parameters->index_colorspace,
//            "Index-Type. Graph generation is not compatible with colorspace");
//    gem_indexer_cond_error_msg_msg(parameters->index_run_length,
//            "Index-Type. Graph generation is not compatible with RL-index");
//  }
  // Input file name
  gem_indexer_cond_error_msg(parameters->input_multifasta_file_name==NULL,
      "Parsing arguments. Please specify an input filename (--input)");
  // Output file name
  if (parameters->output_index_file_name==NULL) {
    parameters->output_index_file_name_prefix = gem_strrmext(gem_strbasename(parameters->input_multifasta_file_name));
    parameters->output_index_file_name = gem_strcat(parameters->output_index_file_name_prefix,".gem");
  } else {
    parameters->output_index_file_name_prefix = parameters->output_index_file_name;
    parameters->output_index_file_name = gem_strcat(parameters->output_index_file_name_prefix,".gem");
  }
  // Info file
  if (!parameters->info_file_name_provided)  {
    parameters->info_file_name = gem_strcat(parameters->output_index_file_name_prefix,".info");
  }
  // System
  if (parameters->max_memory==0) {
    parameters->max_memory = mm_get_available_mem();
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
  // Parsing command-line options
  gem_timer_t gem_indexer_timer; // Global Indexer Timer
  indexer_parameters_t parameters;
  indexer_parameters_set_defaults(&parameters);
  parse_arguments(argc,argv,&parameters);

  // GEM Runtime setup
  gruntime_init(parameters.num_threads,parameters.max_memory,parameters.tmp_folder);
  gem_info_set_stream(fopen(parameters.info_file_name,"wb")); // Set INFO file
  TIMER_RESTART(&gem_indexer_timer); // Start global timer

  // GEM Archive Builder
  fm_t* const index_file = fm_open_file(parameters.output_index_file_name,FM_WRITE);
  const archive_type type = (parameters.bisulfite_index) ? archive_dna_bisulfite : archive_dna;
  archive_builder_t* const archive_builder = archive_builder_new(index_file,
      parameters.output_index_file_name_prefix,type,
      parameters.index_complement,parameters.complement_size_threshold,
      parameters.ns_threshold,parameters.sa_sampling_rate,
      parameters.text_sampling_rate,parameters.indexed_reverse_text,
      parameters.num_threads,parameters.max_memory);

  // Process MultiFASTA
  indexer_process_multifasta(archive_builder,&parameters);

  // Generate BWT
  indexer_generate_bwt(archive_builder,&parameters);
  indexer_write_index(archive_builder,&parameters); // Write Index

  // Generate BWT reverse
  if (parameters.indexed_reverse_text) {
    indexer_generate_bwt_reverse(archive_builder,&parameters);
    indexer_write_index_reverse(archive_builder,&parameters); // Write Index Reverse
  }

  /*
   * Display end banner
   */
  TIMER_STOP(&gem_indexer_timer);
  if (parameters.verbose) {
    tfprintf(stderr,"[GEM Index '%s' was successfully built in %2.3f min.] (see '%s.info' for further info)\n",
        parameters.output_index_file_name,TIMER_GET_TOTAL_M(&gem_indexer_timer),
        parameters.output_index_file_name_prefix);
  }

  // Free
  indexer_cleanup(archive_builder,&parameters);
  gruntime_destroy();

  return 0;
}

