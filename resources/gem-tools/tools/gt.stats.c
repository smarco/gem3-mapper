/*
 * PROJECT: GEM-Tools library
 * FILE: gt.stats.c
 * DATE: 02/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Utility to retrieve very naive stats from {MAP,SAM,FASTQ} files
 */

#include <getopt.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "gem_tools.h"

#define GT_STATS_OUT_FILE stdout

typedef struct {
  /* [Input] */
  char *name_input_file;
  char *name_reference_file;
  char *name_output_file;
  FILE* output_file;
  bool mmap_input;
  bool paired_end;
  uint64_t num_reads;
  /* [Tests] */
  bool first_map;
  bool maps_profile;
  bool mismatch_transitions;
  bool mismatch_quality;
  bool splitmaps_profile;
  bool indel_profile;
  bool population_profile;
  /* [MAP Specific] */
  bool use_only_decoded_maps;
  /* [Output] */
  bool verbose;
  bool compact; // FIXME Deleteme
  /* [Misc] */
  uint64_t num_threads;
} gt_stats_args;

gt_stats_args parameters = {
    /* [Input] */
    .name_input_file=NULL,
    .name_reference_file=NULL,
    .name_output_file=NULL,
    .mmap_input=false,
    .paired_end=false,
    .num_reads=0,
    .output_file=NULL,
    /* [Tests] */
    .first_map=false,
    .maps_profile = false,
    .mismatch_transitions = false,
    .mismatch_quality = false,
    .splitmaps_profile = false,
    .indel_profile = false,
    .population_profile = false,
    /* [MAP Specific] */
    .use_only_decoded_maps = false,
    /* [Output] */
    .verbose=false,
    .compact = false,
    /* [Misc] */
    .num_threads=1,
};
/*
 * STATS Print results
 */
void gt_stats_print_stats(gt_stats* const stats,uint64_t num_reads,const bool paired_end) {
  /*
   * General.Stats (Reads,Alignments,...)
   */
  fprintf(parameters.output_file,"[GENERAL.STATS]\n");
  gt_stats_print_general_stats(parameters.output_file,stats,num_reads,paired_end);
  /*
   * Maps
   */
  if (parameters.maps_profile) {
    fprintf(parameters.output_file,"[MAPS.PROFILE]\n");
    gt_stats_print_maps_stats(parameters.output_file,stats,num_reads,paired_end);
  }
  /*
   * Print Quality Scores vs Errors/Misms
   */
  if (parameters.mismatch_quality && num_reads>0) {
    const gt_maps_profile* const maps_profile = stats->maps_profile;
    if (maps_profile->total_mismatches > 0) {
      fprintf(parameters.output_file,"[MISMATCH.QUALITY]\n");
      gt_stats_print_qualities_error_distribution(
          parameters.output_file,maps_profile->qual_score_misms,maps_profile->total_mismatches);
    }
    if (maps_profile->total_errors_events > 0) {
      fprintf(parameters.output_file,"[ERRORS.QUALITY]\n");
      gt_stats_print_qualities_error_distribution(
          parameters.output_file,maps_profile->qual_score_errors,maps_profile->total_errors_events);
    }
  }
  /*
   * Print Mismatch transition table
   */
  if (parameters.mismatch_transitions && num_reads>0) {
    const gt_maps_profile* const maps_profile = stats->maps_profile;
    if (maps_profile->total_mismatches > 0) {
      fprintf(parameters.output_file,"[MISMATCH.TRANSITIONS]\n");
      fprintf(parameters.output_file,"MismsTransitions\n");
      gt_stats_print_misms_transition_table(
          parameters.output_file,maps_profile->misms_transition,maps_profile->total_mismatches);
      fprintf(parameters.output_file,"MismsTransitions.1-Nucleotide.Context\n");
      gt_stats_print_misms_transition_table_1context(
          parameters.output_file,maps_profile->misms_1context,maps_profile->total_mismatches);
    }
  }
  /*
   * Print Splitmaps profile
   */
  if (parameters.splitmaps_profile) {
    fprintf(parameters.output_file,"[SPLITMAPS.PROFILE]\n");
    gt_stats_print_split_maps_stats(parameters.output_file,stats,parameters.paired_end);
  }
  /*
   * Print Population profile
   */
  if (parameters.population_profile) {
    fprintf(parameters.output_file,"[POPULATION.PROFILE]\n");
    gt_stats_print_population_stats(parameters.output_file,stats,num_reads,parameters.paired_end);
  }
}
void gt_stats_print_stats_compact(gt_stats* const stats,uint64_t num_reads,const bool paired_end) {
  // #mapped, %mapped
  const uint64_t num_templates = paired_end ? num_reads>>1 : num_reads; // SE => 1 template. PE => 1 template
  fprintf(parameters.output_file,"%" PRIu64 ",",stats->num_mapped);
  fprintf(parameters.output_file,"%2.3f,",num_templates?100.0*(float)stats->num_mapped/(float)num_templates:0.0);
  // #unmapped, %unmapped
  const uint64_t unmapped = num_templates-stats->num_mapped;
  fprintf(parameters.output_file,"%" PRIu64 ",",unmapped);
  fprintf(parameters.output_file,"%2.3f,",num_templates?100.0*(float)unmapped/(float)num_templates:0.0);
  // MMap(maps/alg)
  fprintf(parameters.output_file,"%2.3f,",stats->num_mapped?(float)stats->num_maps/(float)stats->num_mapped:0.0);
  // Bases.aligned(%)
  fprintf(parameters.output_file,"%2.3f,",GT_GET_PERCENTAGE(stats->maps_profile->total_bases_matching,stats->maps_profile->total_bases));
  // Bases.trimmed(%)
  fprintf(parameters.output_file,"%2.3f,",GT_GET_PERCENTAGE(stats->maps_profile->total_bases_trimmed,stats->maps_profile->total_bases));
  // #Uniq-0, %Uniq-0
  const uint64_t all_uniq = stats->uniq[GT_STATS_UNIQ_RANGE_0]+
      stats->uniq[GT_STATS_UNIQ_RANGE_1]+stats->uniq[GT_STATS_UNIQ_RANGE_2]+
      stats->uniq[GT_STATS_UNIQ_RANGE_3]+stats->uniq[GT_STATS_UNIQ_RANGE_10]+
      stats->uniq[GT_STATS_UNIQ_RANGE_50]+stats->uniq[GT_STATS_UNIQ_RANGE_100]+
      stats->uniq[GT_STATS_UNIQ_RANGE_500]+stats->uniq[GT_STATS_UNIQ_RANGE_BEHOND];
  fprintf(parameters.output_file,"%" PRIu64 ",",all_uniq);
  fprintf(parameters.output_file,"%2.3f\n",num_templates?100.0*(float)all_uniq/(float)num_templates:0.0);
}

/*
 * CORE functions
 */
void gt_stats_parallel_generate_stats() {
  // Stats info
  gt_stats_analysis stats_analysis = GT_STATS_ANALYSIS_DEFAULT();
  gt_stats** stats = gt_calloc(parameters.num_threads,gt_stats*,false);

  // Select analysis
  stats_analysis.first_map = parameters.first_map;
  stats_analysis.maps_profile = parameters.maps_profile|parameters.mismatch_quality|parameters.mismatch_transitions;
  stats_analysis.nucleotide_stats = true;
  stats_analysis.splitmap_profile = parameters.splitmaps_profile;
  stats_analysis.indel_profile = parameters.indel_profile;
  stats_analysis.population_profile = parameters.population_profile;
  stats_analysis.use_map_counters = !parameters.use_only_decoded_maps;

  // Open file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,parameters.mmap_input);
  parameters.output_file = (parameters.name_output_file==NULL) ? stdout : fopen(parameters.name_output_file,"w+");

  gt_sequence_archive* sequence_archive = NULL;
  if (stats_analysis.indel_profile) {
    sequence_archive = gt_sequence_archive_new(GT_CDNA_ARCHIVE);
    gt_input_file* const reference_file = gt_input_file_open(parameters.name_reference_file,false);
    if (gt_input_multifasta_parser_get_archive(reference_file,sequence_archive)!=GT_IFP_OK) {
      gt_fatal_error_msg("Error parsing reference file '%s'\n",parameters.name_reference_file);
    }
    gt_input_file_close(reference_file);
  }

  // Parallel reading+process
#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(parameters.num_threads)
#endif
  {
#ifdef HAVE_OPENMP
    uint64_t tid = omp_get_thread_num();
#else
    uint64_t tid = 0;
#endif

    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);

    gt_status error_code;
    gt_template *template = gt_template_new();
    stats[tid] = gt_stats_new();
    gt_generic_parser_attributes* generic_parser_attribute = gt_input_generic_parser_attributes_new(parameters.paired_end);
    while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,generic_parser_attribute))) {
      if (error_code!=GT_INPUT_STATUS_OK) {
        gt_error_msg("Fatal error parsing file '%s'\n",parameters.name_input_file);
      }

      // Extract stats
      gt_stats_calculate_template_stats(stats[tid],template,sequence_archive,&stats_analysis);
    }

    // Clean
    gt_template_delete(template);
    gt_buffered_input_file_close(buffered_input);
  }

  // Merge stats
  gt_stats_merge(stats,parameters.num_threads);

  /*
   * Print Statistics
   *   Use stats[0]->num_blocks as the number of blocks in a MAP/SAM/FASTA/FASTQ file
   *   is the number of reads in a FASTA/FASTQ
   */
  if (!parameters.compact) {
    gt_stats_print_stats(stats[0],(parameters.num_reads>0)?
        parameters.num_reads:stats[0]->num_blocks,parameters.paired_end);
  } else {
    gt_stats_print_stats_compact(stats[0],(parameters.num_reads>0)?
        parameters.num_reads:stats[0]->num_blocks,parameters.paired_end);
  }

  // Clean
  gt_stats_delete(stats[0]); gt_free(stats);
  gt_input_file_close(input_file);
}

void parse_arguments(int argc,char** argv) {
  struct option* gt_stats_getopt = gt_options_adaptor_getopt(gt_stats_options);
  gt_string* const gt_stats_short_getopt = gt_options_adaptor_getopt_short(gt_stats_options);
  int option, option_index;
  while (true) {
    // Get option &  Select case
    if ((option=getopt_long(argc,argv,
        gt_string_get_string(gt_stats_short_getopt),gt_stats_getopt,&option_index))==-1) break;
    switch (option) {
    /* I/O */
    case 'i': // input
      parameters.name_input_file = optarg;
      break;
    case 200: // mmap-input
      parameters.mmap_input = true;
      gt_fatal_error(NOT_IMPLEMENTED);
      break;
    case 'r': // reference
      parameters.name_reference_file = optarg;
      gt_fatal_error(NOT_IMPLEMENTED);
      break;
    case 'I': // gem-index
      gt_fatal_error(NOT_IMPLEMENTED);
      break;
    case 'p': // paired-end
      parameters.paired_end = true;
      break;
    case 'n': // num-reads
      parameters.num_reads = atol(optarg);
      break;
    case 'o': // output
      parameters.name_output_file = optarg;
      break;
    /* Analysis */
    case 300: // first-map
      parameters.first_map = true;
      break;
    case 'a': // all-tests
      parameters.maps_profile = true;
      parameters.mismatch_transitions = true;
      parameters.mismatch_quality = true;
      parameters.splitmaps_profile = true;
      parameters.population_profile = true;
      break;
    case 'M': // maps-profile
      parameters.maps_profile = true;
      break;
    case 'T': // mismatch-transitions
      parameters.mismatch_transitions = true;
      break;
    case 'Q': // mismatch-quality
      parameters.mismatch_quality = true;
      break;
    case 'R': // rna-profile // FIXME name
      parameters.splitmaps_profile = true;
      break;
    case 'P': // population-profile
      parameters.population_profile = true;
      break;
    case 'D': // indel-profile
      gt_fatal_error(NOT_IMPLEMENTED);
      parameters.indel_profile = true;
      break;
    /* MAP Specific */
    case 400:
      parameters.use_only_decoded_maps = true;
      break;
    /* Misc */
#ifdef HAVE_OPENMP
    case 't':
      parameters.num_threads = atol(optarg);
      break;
#endif
    case 'v':
      parameters.verbose = true;
      break;
    case 'h':
      fprintf(stderr, "USE: ./gt.stats [ARGS]...\n");
      gt_options_fprint_menu(stderr,gt_stats_options,gt_stats_groups,false,false);
      exit(1);
    case 'H':
      fprintf(stderr, "USE: ./gt.stats [ARGS]...\n");
      gt_options_fprint_menu(stderr,gt_stats_options,gt_stats_groups,false,true);
      exit(1);
    case 'J':
      gt_options_fprint_json_menu(stderr,gt_stats_options,gt_stats_groups,false,true);
      exit(1);
      break;
    case '?':
    default:
      gt_fatal_error_msg("Option not recognized");
    }
  }
  /*
   * Checks
   */
  if (parameters.indel_profile && parameters.name_reference_file==NULL) {
    gt_error_msg("To generate the indel-profile, a reference file(.fa/.fasta) or GEMindex(.gem) is required");
  }
  // Free
  gt_string_delete(gt_stats_short_getopt);
}

int main(int argc,char** argv) {
  // GT error handler
  gt_handle_error_signals();

  // Parsing command-line options
  parse_arguments(argc,argv);

  // Extract stats
  gt_stats_parallel_generate_stats();

  return 0;
}


