/*
 * PROJECT: GEMMapper
 * FILE: gem-indexer.c
 * DATE: 1/10/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Encodes a Multi-FASTA file (.fa) into a GEM index
 *   // TODO
 *   - Pending 2-step by no sorting subdominant kmers
 *   --force-fmi-general-index       (default: deduced from content)\n"
 *   --filter-function 'iupac-dna'|'iupac-colorspace-dna'|'none'\n"
 */

#include <pthread.h>
#include "gem_core.h"

/*
 * GEM-indexer Debug
 */
#define GEM_INDEXER_DEBUG_DUMP_EXPLICIT_CHECKED_SA  false
#define GEM_INDEXER_DEBUG_DUMP_CHECKED_BWT          false

/*
 * GEM-indexer Parameters
 */
//  /* 0 */ "Null",
//  /* 1 */ "Unclassified",
//  /* 2 */ "I/O",
//  /* 3 */ "MultiFasta Processing",
//  /* 4 */ "Index Generation",
//  /* 5 */ "System"
//  /* 6 */ "Miscellaneous",
//  /* 7 */ "Extras"
typedef struct {
  /* I/O */
  char* input_multifasta_file_name;
  char* input_graph_file_name;
  char* output_index_file_name;
  char* output_index_file_name_prefix;
  /* MultiFasta Processing */
  bool index_colorspace;
  bool index_run_length;
  indexed_complement_t index_complement;
  uint64_t ns_threshold;
  uint64_t complement_size_threshold;
  /* FM-Index */
  sampling_rate_t sampling_rate;
  bool check_index;
  /* System */
  uint64_t num_threads;
  uint64_t max_memory;
  char* tmp_folder;
  /* Debug */
  bool dump_locator_intervals;
  bool dump_indexed_text;
  bool dump_explicit_sa;
  bool dump_bwt;
  bool dump_run_length_text;
  bool dump_graph_links;
  /* Miscellaneous */
  bool info_file_name_provided;
  char* info_file_name;
  bool verbose;
  /* Extras */
  bool dev_generate_only_bwt;
} gem_indexer_parameters;
// Defaults
gem_indexer_parameters parameters = {
  /* I/O */
  .input_multifasta_file_name=NULL,
  .input_graph_file_name=NULL,
  .output_index_file_name=NULL,
  .output_index_file_name_prefix=NULL,
  /* MultiFasta Processing */
  .index_colorspace=false,
  .index_run_length=false,
  .index_complement=index_complement_auto,
  .ns_threshold=50,
  .complement_size_threshold=BUFFER_SIZE_1G,
  /* FM-Index */
  .sampling_rate=SAMPLING_RATE_16,
  .check_index = false,
  /* System */
  .num_threads=1,
  .max_memory=0,
  .tmp_folder=NULL,
  /* Debug */
  .dump_locator_intervals=true,
  .dump_indexed_text=false,
  .dump_explicit_sa=false,
  .dump_bwt=false,
  .dump_run_length_text=false,
  .dump_graph_links=false,
  /* Miscellaneous */
  .info_file_name_provided=false,
  .info_file_name=NULL,
  .verbose=false,
  /* Extras */
  .dev_generate_only_bwt=false
};
gem_timer_t gem_indexer_timer; // Global Indexer Timer
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
GEM_INLINE void indexer_debug_generate_sa(dna_text_t* const index_text) {
  // Allocate memory for SA
  debug_char_text = (char*) dna_text_get_buffer(index_text);
  debug_text_length = dna_text_get_length(index_text);
  debug_SA = mm_malloc(sizeof(uint64_t)*(debug_text_length));
  uint64_t i;
  for (i=0;i<debug_text_length;++i) debug_SA[i] = i;
  // Build SA
  qsort(debug_SA,debug_text_length,sizeof(uint64_t),(int (*)(const void *,const void *))indexer_debug_suffix_cmp);
}
GEM_INLINE void indexer_debug_check_sa(char* const file_name,dna_text_t* const index_text) {
  // Open file
  FILE* const bwt_file = fopen(file_name,"w");
  // Generate SA
  if (debug_SA==NULL) indexer_debug_generate_sa(index_text);
  // Dump full SA
  uint64_t i, t, num_printed_chars;
  for (i=0;i<debug_text_length;++i) {
    num_printed_chars=0;
    fprintf(bwt_file,"Suffix=%011lu\t\t",debug_SA[i]);
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
GEM_INLINE void indexer_debug_check_bwt(char* const file_name,dna_text_t* const index_text) {
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
GEM_INLINE void indexer_process_multifasta(archive_builder_t* const archive_builder) {
  /*
   * Process input MultiFASTA
   */
  input_file_t* const input_multifasta = (parameters.input_multifasta_file_name==NULL) ?
      input_stream_open(stdin) : input_file_open(parameters.input_multifasta_file_name,false);
  switch (archive_builder->index_type) {
    case fm_dna_classic:
      archive_builder_process_multifasta(archive_builder,input_multifasta,
          parameters.dump_locator_intervals,parameters.dump_indexed_text,parameters.verbose);
      break;
    case fm_dna_run_length:
      archive_builder_process_multifasta(archive_builder,input_multifasta,
          parameters.dump_locator_intervals,parameters.dump_indexed_text,parameters.verbose);
      archive_builder_process_run_length_text(
          archive_builder,parameters.dump_run_length_text,parameters.verbose);
      break;
    case fm_dna_graph: {
      input_file_t* const input_graph = input_file_open(parameters.input_graph_file_name,false);
      archive_builder_process_graph(archive_builder,input_graph,parameters.dump_graph_links,parameters.verbose);
      archive_builder_process_multifasta__graph(archive_builder,input_multifasta,
          parameters.dump_locator_intervals,parameters.dump_indexed_text,parameters.dump_graph_links,parameters.verbose);
      input_file_close(input_graph); // Close MultiFASTA
      break;
    }
    default:
      GEM_INVALID_CASE();
      break;
  }
  input_file_close(input_multifasta); // Close MultiFASTA
}
GEM_INLINE void indexer_generate_bwt(archive_builder_t* const archive_builder) {
  // Build BWT
  archive_builder_build_bwt(archive_builder,parameters.dump_explicit_sa,parameters.dump_bwt,parameters.verbose);
  // DEBUG: Print Explicit Checked-SA => (SApos,SA[SApos...SApos+SAFixLength])
  gem_cond_debug_block(GEM_INDEXER_DEBUG_DUMP_EXPLICIT_CHECKED_SA) {
    indexer_debug_check_sa(gem_strcat(parameters.output_index_file_name_prefix,".check.sa"),archive_builder->enc_text);
  }
  // DEBUG: Print Checked-BWT
  gem_cond_debug_block(GEM_INDEXER_DEBUG_DUMP_CHECKED_BWT) {
    // Generate SA-BWT (to compare with)
    indexer_debug_check_bwt(gem_strcat(parameters.output_index_file_name_prefix,".check.bwt"),archive_builder->enc_text);
    // Free
    if (debug_SA) mm_free(debug_SA);
  }
  // DEVEL: Skip the FM-index generation
  if (parameters.dev_generate_only_bwt) exit(0);
}
GEM_INLINE void indexer_build_index(archive_builder_t* const archive_builder) {
  /*
   * Build Index
   */
  archive_builder_build_index(archive_builder,parameters.check_index,parameters.verbose);
}
/*
 * GEM-Indexer options Menu
 */
option_t gem_indexer_options[] = {
  /* I/O */
  { 'i', "input", REQUIRED, TYPE_STRING, 2 , true, "<input_file>" , "(Multi-FASTA, default=stdin)" },
  { 'g', "graph", REQUIRED, TYPE_STRING, 2 , true, "<graph_file>" , "(GRAPH-File)" },
  { 'o', "output", REQUIRED, TYPE_STRING, 2 , true, "<output_prefix>" , "" },
  /* MultiFasta Processing */
  { 'c', "index-colorspace", NO_ARGUMENT, TYPE_NONE, 3 , true, "" , "(default=false)" },
  { 'r', "index-run-length", NO_ARGUMENT, TYPE_NONE, 3 , true, "" , "(default=false)" },
  { 300, "index-complement", OPTIONAL, TYPE_STRING, 3 , true, "" , "(default=false, simulated at mapping)" },
  { 'N', "strip-unknown-bases-threshold", REQUIRED, TYPE_INT, 3 , true, "'disable'|<integer>" , "(default=50)" },
  { 301, "complement-size-threshold", REQUIRED, TYPE_INT, 3 , true, "<integer>" , "(default=2GB)" },
  /* FM-Index */
  { 's', "sampling-rate", REQUIRED, TYPE_INT, 4 , true, "<sampling_rate>" , "{} (default=16)" },
  { 400, "check-index", NO_ARGUMENT, TYPE_NONE, 4 , true, "", "(default=false)"},
  /* System */
  { 't', "threads", REQUIRED, TYPE_INT, 5 , true, "<number>" , "(default=1)" },
  { 500, "max-memory", REQUIRED, TYPE_STRING, 5 , true, "<maximum-memory>" , "(Eg 2GB)" },
  { 501, "tmp-folder", REQUIRED, TYPE_STRING, 5 , true, "<temporal_dir_path>" , "(/tmp/)" },
  /* Debug/Temporal */
  { 600, "dump-locator-intervals", OPTIONAL, TYPE_NONE, 6 , true, "" , "" },
  { 601, "dump-indexed-text", OPTIONAL, TYPE_NONE, 6 , true, "" , "" },
  { 602, "dump-explicit-sa", OPTIONAL, TYPE_NONE, 6 , true, "" , "" },
  { 603, "dump-bwt", OPTIONAL, TYPE_NONE, 6 , true, "" , "" },
  { 604, "dump-run-length-text", OPTIONAL, TYPE_NONE, 6 , true, "" , "" },
  { 605, "dump-graph-links", OPTIONAL, TYPE_NONE, 6 , true, "" , "" },
  { 606, "debug", NO_ARGUMENT, TYPE_NONE, 6 , true, "" , "" },
  /* Miscellaneous */
  { 'v', "verbose", NO_ARGUMENT, TYPE_NONE, 7 ,true, "", "" },
  { 700, "info-file", REQUIRED, TYPE_STRING, 7 ,false, "<info_file_path>", "" },
  { 'h', "help", NO_ARGUMENT, TYPE_NONE, 7 , true, "" , "(print usage)" },
  { 'H', "help", NO_ARGUMENT, TYPE_NONE, 7 , false, "" , "(print usage + extras)" },
  { 701, "show-license", NO_ARGUMENT, TYPE_NONE, 7 ,true, "", "(print license and exit)"},
  /* Extras/Develop */
  { 800, "bwt", NO_ARGUMENT, TYPE_NONE, 8 ,false, "", "(only generate BWT for benchmarking)" },
  {  0, "", 0, 0, 0, false, "", ""}
};
char* gem_indexer_groups[] = {
  /* 0 */ "Null",
  /* 1 */ "Unclassified",
  /* 2 */ "I/O",
  /* 3 */ "MultiFasta Processing",
  /* 4 */ "Index Generation",
  /* 5 */ "System",
  /* 6 */ "Debug",
  /* 7 */ "Miscellaneous",
  /* 8 */ "Extras"
};
void usage(const bool print_inactive) {
  fprintf(stderr, "USAGE: ./gem-indexer [ARGS]...\n");
  options_fprint_menu(stderr,gem_indexer_options,gem_indexer_groups,false,print_inactive);
}
void parse_arguments(int argc,char** argv) {
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
      parameters.input_multifasta_file_name = optarg;
      break;
    case 'g': // --graph
      parameters.input_graph_file_name = optarg;
      break;
    case 'o': // --output
      parameters.output_index_file_name = optarg;
      break;
    /* MultiFasta Processing */
    case 'c': // --index-colorspace
      parameters.index_colorspace = true;
      break;
    case 'r': // --index-run-length
      parameters.index_run_length = true;
      break;
    case 300: // --index-complement
      parameters.index_complement = (optarg) ? (options_parse_bool(optarg) ? index_complement_yes : index_complement_no ) : index_complement_yes;
      break;
    case 'N': // --strip-unknown-bases-threshold
      parameters.ns_threshold = atol(optarg); // FIXME Disable
      break;
    case 301: // --complement-size-threshold
      gem_cond_fatal_error(input_text_parse_size(optarg,&(parameters.complement_size_threshold)),PARSING_SIZE,"-complement-size-threshold",optarg);
      break;
    /* FM-Index */
    case 's': { // --sampling-rate
      const uint64_t sampling = atol(optarg);
      switch (sampling) {
        case 1:   parameters.sampling_rate=SAMPLING_RATE_1; break;
        case 2:   parameters.sampling_rate=SAMPLING_RATE_2; break;
        case 4:   parameters.sampling_rate=SAMPLING_RATE_4; break;
        case 8:   parameters.sampling_rate=SAMPLING_RATE_8; break;
        case 16:  parameters.sampling_rate=SAMPLING_RATE_16; break;
        case 32:  parameters.sampling_rate=SAMPLING_RATE_32; break;
        case 64:  parameters.sampling_rate=SAMPLING_RATE_64; break;
        case 128: parameters.sampling_rate=SAMPLING_RATE_128; break;
        case 256: parameters.sampling_rate=SAMPLING_RATE_256; break;
        default:
          gem_error_msg("Sampling rate argument not valid. Reset to default (--sampling-rate 16)");
          break;
      }
    }
    break;
    case 400: // --check-index
      parameters.check_index = true;
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
    /* Debug/Temporal */
    case 600: // --dump-locator-intervals
      parameters.dump_locator_intervals = (optarg) ? options_parse_bool(optarg) : true;
      break;
    case 601: // --dump-indexed-text
      parameters.dump_indexed_text = (optarg) ? options_parse_bool(optarg) : true;
      break;
    case 602: // --dump-explicit-sa
      parameters.dump_explicit_sa = (optarg) ? options_parse_bool(optarg) : true;
      break;
    case 603: // --dump-bwt
      parameters.dump_bwt = (optarg) ? options_parse_bool(optarg) : true;
      break;
    case 604:
      parameters.dump_run_length_text = (optarg) ? options_parse_bool(optarg) : true;
      break;
    case 605: // --dump-graph-links
      parameters.dump_graph_links = (optarg) ? options_parse_bool(optarg) : true;
      break;
    case 606: // --debug
      parameters.dump_locator_intervals = true;
      parameters.dump_indexed_text = true;
      parameters.dump_explicit_sa = true;
      parameters.dump_bwt = true;
      parameters.dump_graph_links = true;
      parameters.dump_run_length_text = true;
      parameters.verbose = true;
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
    case 700: // --info-file
      parameters.info_file_name_provided = true;
      parameters.info_file_name = optarg;
      break;
    /* Extras/Develop */
    case 800: // --bwt
      parameters.dev_generate_only_bwt = true;
      parameters.dump_bwt = true;
      parameters.verbose = true;
      break;
    case '?':
    default:
      gem_fatal_error_msg("Option not recognized");
    }
  }
  /*
   * Parameters Check
   */
  // Index type incompatibility list
  if (parameters.input_graph_file_name!=NULL) {
    gem_cond_fatal_error_msg(parameters.index_colorspace,
            "Index-Type. Graph generation is not compatible with colorspace");
    gem_cond_fatal_error_msg(parameters.index_run_length,
            "Index-Type. Graph generation is not compatible with RL-index");
  } else if (parameters.index_colorspace!=parameters.index_run_length) {
    gem_cond_fatal_error_msg(parameters.index_run_length,
            "Index-Type. Colorspace is not compatible with RL-index");
  }
  // Output file name
  if (parameters.output_index_file_name==NULL) {
    gem_cond_fatal_error_msg(parameters.input_multifasta_file_name==NULL,
        "Parsing arguments. Please specify an output file name (--output)");
    parameters.output_index_file_name_prefix = gem_strrmext(gem_strbasename(parameters.input_multifasta_file_name));
    parameters.output_index_file_name = gem_strcat(parameters.output_index_file_name_prefix,".gem");
  } else {
    parameters.output_index_file_name_prefix = parameters.output_index_file_name;
    parameters.output_index_file_name = gem_strcat(parameters.output_index_file_name_prefix,".gem");
  }
  // Info file
  if (!parameters.info_file_name_provided)  {
    parameters.info_file_name = gem_strcat(parameters.output_index_file_name_prefix,".info");
  }
  // System
  if (parameters.max_memory==0) {
    parameters.max_memory = mm_get_available_mem();
  }
  /*
   * Free
   */
  string_delete(getopt_short_string);
  mm_free(getopt_options);
}
void indexer_cleanup(archive_builder_t* const archive_builder) {
  // Archive Builder
  archive_builder_delete(archive_builder);
  // Output file name
//  if (parameters.output_index_file_name_prefix!=parameters.output_index_file_name) {
//    mm_free(parameters.output_index_file_name_prefix);
//  }
  mm_free(parameters.output_index_file_name);
  // Info File Name
  if (!parameters.info_file_name_provided)  {
    mm_free(parameters.info_file_name);
  }
}
/*
 * Main()
 */
int main(int argc,char** argv) {
  // Parsing command-line options
  parse_arguments(argc,argv);

  // GEM Runtime setup
  gem_runtime_init(parameters.max_memory,parameters.tmp_folder,NULL);
  gem_info_set_stream(fopen(parameters.info_file_name,"wb")); // Set INFO file
  PROF_NEW(parameters.num_threads);
  TIMER_RESTART(&gem_indexer_timer); // Start global timer

  // GEM Archive Builder
  fm_t* const index_file = fm_open_file(parameters.output_index_file_name,FM_WRITE);
  const filter_t filter_type = (parameters.index_colorspace) ? Iupac_colorspace_dna : Iupac_dna;
  const index_t index_type =
      (parameters.input_graph_file_name!=NULL) ? fm_dna_graph :
          (parameters.index_run_length ? fm_dna_run_length : fm_dna_classic);
  archive_builder_t* const archive_builder = archive_builder_new(
      index_file,parameters.output_index_file_name_prefix,index_type,filter_type,
      parameters.index_complement,parameters.complement_size_threshold,
      parameters.ns_threshold,parameters.sampling_rate,parameters.num_threads,parameters.max_memory);

  /*
   * Process MultiFASTA
   */
  indexer_process_multifasta(archive_builder);

  /*
   * Generate BWT
   */
  indexer_generate_bwt(archive_builder);

  /*
   * Create FM-Index
   */
  indexer_build_index(archive_builder);

  /*
   * Display end banner
   */
  TIMER_STOP(&gem_indexer_timer);
  if (parameters.verbose) {
    tfprintf(stderr,"[GEM Index '%s' was successfully built in %2.3f min.] (see '%s.info' for further info)\n",
        parameters.output_index_file_name,TIMER_GET_S(&gem_indexer_timer)/60.0,
        parameters.output_index_file_name_prefix);
  }

  // Free
  indexer_cleanup(archive_builder);
  PROF_DELETE();

  return 0;
}

