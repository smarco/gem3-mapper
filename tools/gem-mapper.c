/*
 * PROJECT: GEMMapper
 * FILE: gem-mapper.c
 * DATE:5/12/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Genomic Read Mapper
 */

#include "gem_core.h"

/*
 * GEM-mapper I/O related functions
 */
GEM_INLINE void gem_mapper_load_index(mapper_parameters_t* const parameters) {
  PROF_START_TIMER(GP_MAPPER_LOAD_INDEX);
  // Load archive
  gem_cond_log(parameters->misc.verbose_user,"[Loading GEM index '%s']",parameters->io.index_file_name);
  parameters->archive = archive_read(parameters->io.index_file_name,
      parameters->io.check_index,parameters->misc.verbose_dev);
  PROF_STOP_TIMER(GP_MAPPER_LOAD_INDEX);
}
GEM_INLINE void gem_mapper_open_input(mapper_parameters_t* const parameters) {
  // Open input file
  const mapper_parameters_cuda_t* const cuda = &parameters->cuda;
  const uint64_t input_block_size = (cuda->cuda_enabled) ? cuda->input_block_size : parameters->io.input_block_size;
  if (parameters->io.input_file_name==NULL) {
    gem_cond_log(parameters->misc.verbose_user,"[Reading input file from stdin]");
    switch (parameters->io.input_compression) {
      case FM_GZIPPED_FILE:
        parameters->input_file = input_gzip_stream_open(stdin,input_block_size);
        break;
      case FM_BZIPPED_FILE:
        parameters->input_file = input_bzip_stream_open(stdin,input_block_size);
        break;
      default:
        parameters->input_file = input_stream_open(stdin,input_block_size);
        break;
    }
  } else {
    gem_cond_log(parameters->misc.verbose_user,"[Opening input file '%s']",parameters->io.input_file_name);
    parameters->input_file = input_file_open(parameters->io.input_file_name,input_block_size,false);
  }
  // TODO: Checks
}
GEM_INLINE void gem_mapper_close_input(mapper_parameters_t* const parameters) {
  input_file_close(parameters->input_file);
}
GEM_INLINE void gem_mapper_open_output(mapper_parameters_t* const parameters) {
  // Open output stream
  if (parameters->io.output_file_name==NULL) {
    gem_cond_log(parameters->misc.verbose_user,"[Outputting to stdout]");
    parameters->output_stream = stdout;
  } else {
    gem_cond_log(parameters->misc.verbose_user,"[Outputting to '%s']",parameters->io.output_file_name);
    parameters->output_stream = fopen(parameters->io.output_file_name,"w");
  }
  // Open output file
  const mapper_parameters_cuda_t* const cuda = &parameters->cuda;
  const uint64_t max_output_buffers = (cuda->cuda_enabled) ? cuda->output_num_buffers : parameters->io.output_num_buffers;
  const uint64_t output_buffer_size = (cuda->cuda_enabled) ? cuda->output_buffer_size : parameters->io.output_buffer_size;
  switch (parameters->io.output_compression) {
    case FM_GZIPPED_FILE:
      parameters->output_file = output_gzip_stream_new(parameters->output_stream,max_output_buffers,output_buffer_size);
      break;
    case FM_BZIPPED_FILE:
      parameters->output_file = output_bzip_stream_new(parameters->output_stream,max_output_buffers,output_buffer_size);
      break;
    default:
      parameters->output_file = output_stream_new(parameters->output_stream,max_output_buffers,output_buffer_size);
      break;
  }
}
GEM_INLINE void gem_mapper_close_output(mapper_parameters_t* const parameters) {
  output_file_close(parameters->output_file);
}
GEM_INLINE void gem_mapper_print_profile(mapper_parameters_t* const parameters) {
  // Reduce Stats
  switch (parameters->misc.profile_reduce_type) {
    case reduce_sum: PROF_REDUCE_SUM(); break;
    case reduce_max: PROF_REDUCE_MAX(); break;
    case reduce_min: PROF_REDUCE_MIN(); break;
    case reduce_mean: PROF_REDUCE_MEAN(); break;
    case reduce_sample: PROF_REDUCE_SAMPLE(); break;
    default: GEM_INVALID_CASE(); break;
  }
  // Sys
  system_print_info(gem_info_get_stream());
  // Parameters
  mapper_parameters_print(gem_info_get_stream(),parameters);
  // I/O
  mapper_profile_print_io(gem_info_get_stream());
  // Mapper
  if (!parameters->cuda.cuda_enabled) {
    /*
     * CPU Mapper
     */
    switch (parameters->mapper_type) {
      case mapper_se:
        switch (parameters->search_parameters.mapping_mode) {
          case mapping_adaptive_filtering:
            mapper_profile_print_mapper_adaptive(gem_info_get_stream());
            if (parameters->system.num_threads==1) {
              mapper_profile_print_mapper_adaptive_ranks(gem_info_get_stream());
            }
            mapper_profile_print_region_profile_soft(gem_info_get_stream());
            mapper_profile_print_filtering_verifying(gem_info_get_stream(),true);
            break;
          case mapping_incremental_mapping:
          case mapping_fixed_filtering:
          case mapping_fast:
          case mapping_neighborhood_search:
            break;
          default:
            GEM_INVALID_CASE();
            break;
        }
        break;
      case mapper_pe:
        break;
      case mapper_graph:
        break;
      default:
        break;
    }
  } else {
    /*
     * CUDA Mapper
     */
    switch (parameters->mapper_type) {
      case mapper_se:
        switch (parameters->search_parameters.mapping_mode) {
          case mapping_adaptive_filtering:
            mapper_profile_print_mapper_cuda_adaptive(gem_info_get_stream());
            mapper_profile_print_region_profile_soft(gem_info_get_stream());
            mapper_profile_print_filtering_verifying(gem_info_get_stream(),false);
            mapper_profile_print_archive_search_group(gem_info_get_stream());
            break;
          case mapping_incremental_mapping:
          case mapping_fixed_filtering:
          case mapping_fast:
          case mapping_neighborhood_search:
          default:
            GEM_INVALID_CASE();
            break;
        }
        break;
      case mapper_pe:
        break;
      case mapper_graph:
        break;
      default:
        break;
    }
  }
}
/*
 * GEM-Mapper options Menu
 */
option_t gem_mapper_options[] = {
  /* I/O */
  { 'I', "index", REQUIRED, TYPE_STRING, 2, true, "<index_file.gem>" , "" },
  { 200, "check-index", NO_ARGUMENT, TYPE_NONE, 2, false, "" , "" },
  { 'i', "input", REQUIRED, TYPE_STRING, 2, true, "<file>" , "(FASTA/FASTQ, default=stdin)" },
  { '1', "i1", REQUIRED, TYPE_STRING, 2, true, "<file>" , "(paired-end, end-1)" },
  { '2', "i2", REQUIRED, TYPE_STRING, 2, true, "<file>" , "(paired-end, end-2)" },
  { 201, "gzip-input", NO_ARGUMENT, TYPE_NONE, 2, true, "", "(gzip input)" },
  { 202, "bzip-input", NO_ARGUMENT, TYPE_NONE, 2, true, "", "(bzip input)" },
  { 203, "input-model", REQUIRED, TYPE_STRING, 2, false, "<input_block_size,num_buffers,num_records>", "(default=64M,2c,5K)" },
  { 'o', "output", REQUIRED, TYPE_STRING, 2, true, "<output_prefix>" , "(default=stdout)" },
  { 205, "gzip-output", NO_ARGUMENT, TYPE_NONE, 2, true, "", "(gzip output)" },
  { 206, "bzip-output", NO_ARGUMENT, TYPE_NONE, 2, true, "", "(bzip output)" },
  { 207, "output-model", REQUIRED, TYPE_STRING, 2, false, "<buffer_size,num_buffers>", "(default=4M,3c)" },
  /* Qualities */
  { 'q', "quality-format", REQUIRED, TYPE_STRING, 3, true, "'ignore'|'offset-33'|'offset-64'", "(default=offset-33)" },
  { 'Q', "quality-model", REQUIRED, TYPE_STRING, 3, false, "'gem'|'flat'", "(default=gem)" },
  { 300, "gem-quality-threshold", REQUIRED, TYPE_INT, 3, true, "<number>", "(default=26, that is e<=2e-3)" },
  /* Single-end Alignment */
  { 400, "mapping-mode", REQUIRED, TYPE_STRING, 4, false, "'incremental'|'adaptive'|'fixed'|'fast'|'brute-force'" , "(default=fast)" },
  { 401, "mapping-degree", REQUIRED, TYPE_FLOAT, 4, false, "<number|percentage>" , "(default=0)" },
#ifdef HAVE_CUDA
  { 402, "cuda", NO_ARGUMENT, TYPE_NONE, 4, true, "", ""},
#endif
  { 'e', "max-search-error", REQUIRED, TYPE_FLOAT, 4, true, "<number|percentage>" , "(default=0.04, 4%)" },
  { 'E', "max-filtering-error", REQUIRED, TYPE_FLOAT, 4, false, "<number|percentage>" , "(default=0.2, 20%)" },
  { 's', "complete-strata-after-best", REQUIRED, TYPE_FLOAT, 4, true, "<number|percentage>" , "(default=0)" },
  { 403, "min-matching-length", REQUIRED, TYPE_FLOAT, 4, false, "<number|percentage>" , "(default=0.20, 20%)" },
  { 404, "max-search-matches", REQUIRED, TYPE_INT, 4, true, "<number>" , "(unlimited by default)" },
  { 405, "mismatch-alphabet", REQUIRED, TYPE_STRING, 4, false, "<symbols>" , "(default='ACGT')" },
  /* Paired-end Alignment */
  { 'p', "paired-end-alignment", NO_ARGUMENT, TYPE_NONE, 5, true, "" , "" },
  // { 500, "mate-pair-alignment", NO_ARGUMENT, TYPE_NONE, 5, true, "" , "" }, // TODO
  { 501, "min-template-length", REQUIRED, TYPE_INT, 5, true, "<number>" , "(default=0)" },
  { 502, "max-template-length", REQUIRED, TYPE_INT, 5, true, "<number>" , "(default=inf)" },
  { 503, "pair-orientation", REQUIRED, TYPE_STRING, 5, true, "'FR'|'RF'|'FF'|'RR'" , "(default=FR)" },
  { 504, "discordant-pair-orientation", REQUIRED, TYPE_STRING, 5, true, "'FR'|'RF'|'FF'|'RR'" , "(default=RF,FF,RR)" },
  { 505, "pair-layout", REQUIRED, TYPE_STRING, 5, true, "'separate'|'overlap'|'contain'|'dovetail'" , "(default=separated,overlap,contain,dovetail)" },
  { 506, "map-both-ends", NO_ARGUMENT, TYPE_NONE, 5, false, "" , "(default=false)" },
  { 507, "max-extendable-candidates", REQUIRED, TYPE_INT, 5, false, "<number>" , "(default=20)" },
  { 508, "max-matches-per-extension", REQUIRED, TYPE_INT, 5, false, "<number>" , "(default=2)" },
  /* Alignment Score */
  { 600, "alignment-model", REQUIRED, TYPE_STRING, 6, false, "'none'|'hamming'|'edit'|'gap-affine'" , "(default=edit)" },
  { 601, "gap-affine-penalties", REQUIRED, TYPE_STRING, 6, false, "A,B,O,X" , "(default=1,4,6,1)" },
  { 'A', "matching-score", REQUIRED, TYPE_INT, 6, false, "" , "(default=1)" },
  { 'B', "mismatch-penalty", REQUIRED, TYPE_INT, 6, false, "" , "(default=4)" },
  { 'O', "gap-open-penalty", REQUIRED, TYPE_INT, 6, false, "" , "(default=6)" },
  { 'X', "gap-extension-penalty", REQUIRED, TYPE_INT, 6, false, "" , "(default=1)" },
  /* Mapping Quality */
  /* Reporting */
  { 'F', "output-format", REQUIRED, TYPE_STRING, 8, true, "'MAP'|'SAM'" , "(default=MAP)" },
  { 'D', "min-decoded-strata", REQUIRED, TYPE_FLOAT, 8, false, "<number|percentage>" , "(stratum-wise, default=0)" },
  { 'd', "max-decoded-matches", REQUIRED, TYPE_INT, 8, false, "<number>|'all'" , "(stratum-wise, default=20)" },
  { 'm', "min-reported-matches", REQUIRED, TYPE_INT, 8, false, "<number>|'all'" , "(default=1)" },
  { 'M', "max-reported-matches", REQUIRED, TYPE_INT, 8, true, "<number>|'all'" , "(default=all)" },
  /* System */
  { 't', "threads", REQUIRED, TYPE_STRING, 9, true, "<number>" , "(default=c)" },
  { 900, "max-memory", REQUIRED, TYPE_STRING, 9, true, "<maximum-memory>" , "(Eg 2GB)" },
  { 901, "tmp-folder", REQUIRED, TYPE_STRING, 9, true, "<temporal_dir_path>" , "(default=/tmp/)" },
  /* CUDA Settings */
  { 1000, "threads-cuda", REQUIRED, TYPE_STRING, 10, false, "<generating>,<selecting>" , "(default=1c,1c)" },
  { 1001, "cuda-search-groups", REQUIRED, TYPE_STRING, 10, false, "<num_groups,buffer_size>" , "(default=3c,16M)" },
  /* Presets/Hints */
  { 1100, "technology", REQUIRED, TYPE_STRING, 11, false, "'hiseq'|'miseq'|'454'|'ion-torrent'|'pacbio'|'nanopore'|'moleculo'" , "(default=hiseq)" },
  { 1101, "reads-model", REQUIRED, TYPE_STRING, 11, false, "<average_length>[,<std_length>]" , "(default=150,50)" },
  /* Debug */
  { 1200, "check-alignments", REQUIRED, TYPE_STRING, 12, false, "'check-correct'|'check-best'|'check-complete'" , "(default=check-correct,check-best)" },
  /* Miscellaneous */
  { 1300, "profile", OPTIONAL, TYPE_NONE, 13, false, "'sum'|'min'|'max'|'mean'|'sample'" , "(disabled)" },
  { 'v', "verbose", OPTIONAL, TYPE_NONE, 13, true, "'quiet'|'user'|'dev'" , "(default=user)" },
  { 'h', "help", NO_ARGUMENT, TYPE_NONE, 13, true, "" , "(print usage)" },
  { 'H', "help", NO_ARGUMENT, TYPE_NONE, 13, false, "" , "(print usage + extras)" },
  { 0, NULL, 0, 0, 0, 0, NULL, NULL}
};
char* gem_mapper_groups[] = {
  /*  0 */ "Null",
  /*  1 */ "Unclassified",
  /*  2 */ "I/O",
  /*  3 */ "Qualities",
  /*  4 */ "Single-end Alignment",
  /*  5 */ "Paired-end Alignment",
  /*  6 */ "Alignment Score",
  /*  7 */ "Mapping Quality",
  /*  8 */ "Reporting",
  /*  9 */ "System",
  /* 10 */ "CUDA Settings",
  /* 11 */ "Presets/Hints",
  /* 12 */ "Debug",
  /* 13 */ "Miscellaneous",
};
void usage(const bool print_inactive) {
  fprintf(stderr, "USAGE: ./gem-mapper [ARGS]...\n");
  options_fprint_menu(stderr,gem_mapper_options,gem_mapper_groups,true,print_inactive);
}
GEM_INLINE int parse_arguments_system_integer(char* const argument,uint64_t* const value) {
  char* argument_ptr = argument;
  // Textual
  if (gem_strcaseeq(argument_ptr,"c")) {
    *value = system_get_num_processors();
    return 0;
  }
  // Number
  double number;
  const int error = input_text_parse_double((const char** const)&argument_ptr,&number);
  if (error) return error;
  // Textual
  if (gem_strcaseeq(argument_ptr,"c")) {
    *value = number*(double)system_get_num_processors();
  } else {
    *value = number;
  }
  return 0;
}
void parse_arguments(int argc,char** argv,mapper_parameters_t* const parameters) {
  // Check number of parameters (quick usage & exit)
  if (argc<=1) {
    usage(false);
    exit(0);
  }
  // Set CMD line
  parameters->argc = argc;
  parameters->argv = argv;
  // Parse parameters
  struct option* getopt_options = options_adaptor_getopt(gem_mapper_options);
  string_t* const getopt_short_string = options_adaptor_getopt_short(gem_mapper_options);
  char* const getopt_short = string_get_buffer(getopt_short_string);
  int option, option_index;
  while (true) {
    // Get option &  Select case
    if ((option=getopt_long(argc,argv,getopt_short,getopt_options,&option_index))==-1) break;
    switch (option) {
    /* I/O */
    case 'I': // --index
      parameters->io.index_file_name = optarg;
      break;
    case 200: // --check-index
      parameters->io.check_index = true;
      break;
    case 'i': // --input
      parameters->io.input_file_name = optarg; // TODO Multiple input files
      break;
    case '1': // --i1
      parameters->io.input_file_name_end1 = optarg;
      break;
    case '2': // --i2
      parameters->io.input_file_name_end2 = optarg;
      break;
    case 201: // --gzip-input
      parameters->io.input_compression = FM_GZIPPED_FILE;
      break;
    case 202: // --bzip-input
      parameters->io.input_compression = FM_BZIPPED_FILE;
      break;
    case 203: { // --input-model=64M,2c,5K
      char *input_block_size=NULL, *num_buffers=NULL, *num_records=NULL;
      const int num_arguments = input_text_parse_csv_arguments(optarg,3,&input_block_size,&num_buffers,&num_records);
      gem_cond_fatal_error_msg(num_arguments!=3,"Option '--input-model' wrong number of arguments");
      // Parse input-buffer size
      gem_cond_fatal_error_msg(input_text_parse_size(input_block_size,&parameters->io.input_block_size),
          "Option '--input-model'. Error parsing 'input_block_size'");
      // Parse number of buffers
      gem_cond_fatal_error_msg(parse_arguments_system_integer(num_buffers,&parameters->io.input_num_buffers),
          "Option '--input-model'. Error parsing 'num_buffers'");
      // Parse number of records per buffer
      gem_cond_fatal_error_msg(input_text_parse_size(num_records,&parameters->io.input_buffer_lines),
          "Option '--input-model'. Error parsing 'num_records'");
      parameters->io.input_buffer_lines *= (2*4); // 2l-Paired x 4l-FASTQRecord
      // Propagate settings to CUDA
      parameters->cuda.input_block_size = parameters->io.input_block_size;
      parameters->cuda.input_num_buffers = parameters->io.input_num_buffers;
      parameters->cuda.input_buffer_lines = parameters->io.input_buffer_lines;
      break;
    }
    case 'o': // --output
      parameters->io.output_file_name = optarg;
      break;
    case 205: // --gzip-output
      parameters->io.output_compression = FM_GZIPPED_FILE;
      break;
    case 206: // --bzip-output
      parameters->io.output_compression = FM_BZIPPED_FILE;
      break;
    case 207: { // --output-model=4M,4c
      char *buffer_size=NULL, *num_buffers=NULL;
      const int num_arguments = input_text_parse_csv_arguments(optarg,2,&buffer_size,&num_buffers);
      gem_cond_fatal_error_msg(num_arguments!=2,"Option '--output-model' wrong number of arguments");
      // Parse output-buffer size
      gem_cond_fatal_error_msg(input_text_parse_size(buffer_size,&parameters->io.output_buffer_size),
          "Option '--output-model'. Error parsing 'buffer_size'");
      // Parse number of output-buffers
      gem_cond_fatal_error_msg(parse_arguments_system_integer(num_buffers,&parameters->io.output_num_buffers),
          "Option '--output-model'. Error parsing 'num_buffers'");
      // Propagate settings to CUDA
      parameters->cuda.output_buffer_size = parameters->io.output_buffer_size;
      parameters->cuda.output_num_buffers = parameters->io.output_num_buffers;
      break;
    }
    /* Qualities */
    case 'q': // --quality-format
      if (gem_strcaseeq(optarg,"ignore")) {
        parameters->search_parameters.quality_format=qualities_ignore;
      } else if (gem_strcaseeq(optarg,"offset-33")) {
        parameters->search_parameters.quality_format=qualities_offset_33;
      } else if (gem_strcaseeq(optarg,"offset-64")) {
        parameters->search_parameters.quality_format=qualities_offset_64;
      } else {
        gem_fatal_error_msg("Option '-q|--quality-format' must be 'ignore', 'offset-33' or 'offset-64'");
      }
      break;
    case 'Q': // --quality-model
      if (gem_strcaseeq(optarg,"gem")) {
        parameters->search_parameters.quality_model=quality_model_type_gem;
        // TODO parameters->quality_score=gem_quality_score;
      } else if (gem_strcaseeq(optarg,"flat")) {
        parameters->search_parameters.quality_model=quality_model_type_flat;
        // TODO parameters->quality_score=flat_quality_score;
      } else
      gem_fatal_error_msg("Option '-Q|--quality-model' must be 'gem' or 'flat'");
      break;
    case 300: // --gem-quality-threshold
      input_text_parse_extended_uint64(optarg,&parameters->search_parameters.quality_threshold);
      break;
    /* Single-end Alignment */
    case 400: // --mapping-mode
      if (gem_strcaseeq(optarg,"incremental")) {
        parameters->search_parameters.mapping_mode = mapping_incremental_mapping;
        break;
      }
      if (gem_strcaseeq(optarg,"adaptive")) {
        parameters->search_parameters.mapping_mode = mapping_adaptive_filtering;
        break;
      }
      if (gem_strcaseeq(optarg,"fixed")) {
        parameters->search_parameters.mapping_mode = mapping_fixed_filtering;
        break;
      }
      if (gem_strcaseeq(optarg,"fast")) {
        parameters->search_parameters.mapping_mode = mapping_fast;
        break;
      }
      if (gem_strcaseeq(optarg,"brute-force")) {
        parameters->search_parameters.mapping_mode = mapping_neighborhood_search;
        break;
      }
      gem_fatal_error_msg("Option '--mapping-mode' must be 'incremental'|'adaptive'|'fixed'|'fast'|'brute-force'");
      break;
    case 401: // --mapping-degree
      parameters->search_parameters.mapping_degree = atof(optarg);
      break;
    case 402: // --cuda
      if (!bpm_gpu_support()) GEM_CUDA_NOT_SUPPORTED();
      parameters->cuda.cuda_enabled = true;
      break;
    case 'e': // --max-search-error
      parameters->search_parameters.max_search_error = atof(optarg);
      break;
    case 'E': // --max-filtering-error
      parameters->search_parameters.max_filtering_error = atof(optarg);
      break;
    case 's': // --complete-strata-after-best
      input_text_parse_extended_double(optarg,(double*)&parameters->search_parameters.complete_strata_after_best);
      break;
    case 403: // --min-matching-length
      parameters->search_parameters.min_matching_length = atof(optarg);
      break;
    case 404: // --max-search-matches
      parameters->search_parameters.max_search_matches = atol(optarg);
      break;
    case 405: { // --mismatch-alphabet
      const char* const mismatch_alphabet = optarg;
      approximate_search_configure_replacements(
          &parameters->search_parameters,mismatch_alphabet,gem_strlen(mismatch_alphabet));
      break;
    }
    /* Paired-end Alignment */
    case 'p': // --paired-end-alignment
      parameters->paired_end.paired_end = true;
      break;
//    case 500: // --mate-pair-alignment
//      input_text_parse_extended_uint64(optarg,&parameters->min_template_size);
//      break;
    case 501: // --min-template-length
      input_text_parse_extended_uint64(optarg,&parameters->paired_end.min_template_length);
      break;
    case 502: // --max-template-length
      input_text_parse_extended_uint64(optarg,&parameters->paired_end.max_template_length);
      break;
    case 503: { // --pair-orientation in {'FR'|'RF'|'FF'|'RR'}
      // Start parsing
      char *pair_orientation = strtok(optarg,",");
      while (pair_orientation!=NULL) {
        if (gem_strcaseeq(pair_orientation,"FR")) { parameters->paired_end.pair_orientation_FR = true; continue; }
        if (gem_strcaseeq(pair_orientation,"RF")) { parameters->paired_end.pair_orientation_RF = true; continue; }
        if (gem_strcaseeq(pair_orientation,"FF")) { parameters->paired_end.pair_orientation_FF = true; continue; }
        if (gem_strcaseeq(pair_orientation,"RR")) { parameters->paired_end.pair_orientation_RR = true; continue; }
        gem_fatal_error_msg("Option '--pair-orientation' must be 'FR'|'RF'|'FF'|'RR'");
        pair_orientation = strtok(NULL,",");
      }
      break;
    }
    case 504: { // --discordant-pair-orientation in {'FR'|'RF'|'FF'|'RR'}
      // Start parsing
      char *discordant_pair_orientation = strtok(optarg,",");
      while (discordant_pair_orientation!=NULL) {
        if (gem_strcaseeq(discordant_pair_orientation,"FR")) { parameters->paired_end.discordant_pair_orientation_FR = true; continue; }
        if (gem_strcaseeq(discordant_pair_orientation,"RF")) { parameters->paired_end.discordant_pair_orientation_RF = true; continue; }
        if (gem_strcaseeq(discordant_pair_orientation,"FF")) { parameters->paired_end.discordant_pair_orientation_FF = true; continue; }
        if (gem_strcaseeq(discordant_pair_orientation,"RR")) { parameters->paired_end.discordant_pair_orientation_RR = true; continue; }
        gem_fatal_error_msg("Option '--discordant-pair-orientation' must be 'FR'|'RF'|'FF'|'RR'");
        discordant_pair_orientation = strtok(NULL,",");
      }
      break;
    }
    case 505: { // --pair-layout in {'separate'|'overlap'|'contain'|'dovetail'}
      // Start parsing
      char *pair_layout = strtok(optarg,",");
      while (pair_layout!=NULL) {
        if (gem_strcaseeq(pair_layout,"separate")) { parameters->paired_end.pair_layout_separate = true; continue; }
        if (gem_strcaseeq(pair_layout,"overlap"))  { parameters->paired_end.pair_layout_overlap = true; continue; }
        if (gem_strcaseeq(pair_layout,"contain"))  { parameters->paired_end.pair_layout_contain = true; continue; }
        if (gem_strcaseeq(pair_layout,"dovetail")) { parameters->paired_end.pair_layout_dovetail = true; continue; }
        gem_fatal_error_msg("Option '--pair-layout' must be 'separate'|'overlap'|'contain'|'dovetail'");
        pair_layout = strtok(NULL,",");
      }
      break;
    }
    case 506: // --map-both-ends
      parameters->paired_end.map_both_ends = true;
      break;
    case 507: // --max-extendable-candidates
      input_text_parse_extended_uint64(optarg,&parameters->paired_end.max_extendable_candidates);
      break;
    case 508: // --max-matches-per-extension
      input_text_parse_extended_uint64(optarg,&parameters->paired_end.max_matches_per_extension);
      break;
    /* Alignment Score */
    case 600: // --alignment-model
      if (gem_strcaseeq(optarg,"none")) {
        parameters->search_parameters.alignment_model = alignment_model_none;
      } else if (gem_strcaseeq(optarg,"hamming")) {
        parameters->search_parameters.alignment_model = alignment_model_hamming;
      } else if (gem_strcaseeq(optarg,"edit") || gem_strcaseeq(optarg,"levenshtein") ) {
        parameters->search_parameters.alignment_model = alignment_model_levenshtein;
      } else if (gem_strcaseeq(optarg,"gap-affine'")) {
        parameters->search_parameters.alignment_model = alignment_model_gap_affine;
      } else {
        gem_fatal_error_msg("Option '--alignment-model' must be 'none'|'hamming'|'edit'|'gap-affine'");
      }
      break;
    /* Mapping Quality */
    /* Reporting */
    case 'F': // --output-format
      if (gem_strcaseeq(optarg,"MAP")) {
        parameters->io.output_format = MAP;
      } else if (gem_strcaseeq(optarg,"SAM")) {
        parameters->io.output_format = SAM;
      } else {
        gem_fatal_error_msg("Option '-F|--output-format' must be 'MAP' or 'SAM'");
      }
      break;
    case 'D': // --min-decoded-strata
      input_text_parse_extended_double(optarg,&parameters->select_parameters.min_decoded_strata);
      break;
    case 'd': // --max-decoded-matches
      input_text_parse_extended_uint64(optarg,&parameters->select_parameters.max_decoded_matches);
      break;
    case 'm': // --min-reported-matches
      input_text_parse_extended_uint64(optarg,&parameters->select_parameters.min_reported_matches);
      break;
    case 'M': // --max-reported-matches
      input_text_parse_extended_uint64(optarg,&parameters->select_parameters.max_reported_matches);
      break;
    /* System */
    case 't': // --threads
      gem_cond_fatal_error_msg(parse_arguments_system_integer(optarg,&parameters->system.num_threads),
          "Option '--threads-cuda'. Error parsing 'num_selecting_threads'");
      break;
    case 900: // --max-memory
      gem_cond_fatal_error(input_text_parse_size(optarg,&(parameters->system.max_memory)),PARSING_SIZE,"--max-memory",optarg);
      break;
    case 901: // --tmp-folder
      parameters->system.tmp_folder = optarg;
      break;
    /* CUDA Settings */
    case 1000: { // --threads-cuda=1c,1c
      char *generating_threads=NULL, *selecting_threads=NULL;
      const int num_arguments = input_text_parse_csv_arguments(optarg,2,&generating_threads,&selecting_threads);
      gem_cond_fatal_error_msg(num_arguments!=2,"Option '--threads-cuda' wrong number of arguments");
      // Generating threads
      gem_cond_fatal_error_msg(parse_arguments_system_integer(generating_threads,&parameters->cuda.num_generating_threads),
          "Option '--threads-cuda'. Error parsing 'num_generating_threads'");
      // Selecting threads
      gem_cond_fatal_error_msg(parse_arguments_system_integer(selecting_threads,&parameters->cuda.num_selecting_threads),
          "Option '--threads-cuda'. Error parsing 'num_selecting_threads'");
      break;
    }
    case 1001: { // --cuda-search-groups=3c,16M
      char *num_groups=NULL, *buffer_size=NULL;
      const int num_arguments = input_text_parse_csv_arguments(optarg,2,&num_groups,&buffer_size);
      gem_cond_fatal_error_msg(num_arguments!=2,"Option '--cuda-search-groups' wrong number of arguments");
      // Number of groups
      gem_cond_fatal_error_msg(parse_arguments_system_integer(num_groups,&parameters->cuda.num_search_groups),
          "Option '--cuda-search-groups'. Error parsing 'num_groups'");
      // Buffer size
      gem_cond_fatal_error_msg(input_text_parse_size(buffer_size,&parameters->cuda.bpm_buffer_size),
          "Option '--cuda-search-groups'. Error parsing 'buffer_size'");
      break;
    }
    /* Presets/Hints */
    case 1100: // --technology in {'hiseq'|'miseq'|'454'|'ion-torrent'|'pacbio'|'nanopore'|'moleculo'}
      if (gem_strcaseeq(optarg,"hiseq")) {
        // TODO Presets for alignment mode, read-length, candidates, ...
      } else if (gem_strcaseeq(optarg,"miseq")) {
        // TODO Presets for alignment mode, read-length, candidates, ...
      } else if (gem_strcaseeq(optarg,"454")) {
        // TODO Presets for alignment mode, read-length, candidates, ...
      } else if (gem_strcaseeq(optarg,"ion-torrent")) {
        // TODO Presets for alignment mode, read-length, candidates, ...
      } else if (gem_strcaseeq(optarg,"pacbio")) {
        // TODO Presets for alignment mode, read-length, candidates, ...
      } else if (gem_strcaseeq(optarg,"nanopore")) {
        // TODO Presets for alignment mode, read-length, candidates, ...
      } else if (gem_strcaseeq(optarg,"moleculo")) {
        // TODO Presets for alignment mode, read-length, candidates, ...
      } else {
        gem_fatal_error_msg("Option '--technology' must be in {'hiseq'|'miseq'|'454'|'ion-torrent'|'pacbio'|'nanopore'|'moleculo'}");
      }
      break;
    case 1101: { // --reads-model <average_length>[,<std_length>]
      char *average_length, *std_length;
      const int num_arguments = input_text_parse_csv_arguments(optarg,2,&average_length,&std_length);
      gem_cond_fatal_error_msg(num_arguments==0,"Option '--reads-model' wrong number of arguments");
      // Average read length
      gem_cond_fatal_error_msg(input_text_parse_size(average_length,&parameters->hints.avg_read_length),
          "Option '--reads-model'. Error parsing 'average_length'");
      if (num_arguments > 1) {
        // Standard deviation read length
        gem_cond_fatal_error_msg(input_text_parse_size(std_length,&parameters->hints.std_read_length),
            "Option '--reads-model'. Error parsing 'std_length'");
      }
      break;
    }
    /* Debug */
    case 1200: { // --check-alignments in {'check-correct'|'check-best'|'check-complete'}
      select_parameters_t* const select_parameters = &parameters->select_parameters;
      select_parameters->check_matches_mask = check_none; // Init
      // Start parsing
      char *check = strtok(optarg,",");
      while (check!=NULL) {
        if (gem_strcaseeq(check,"check-none"))     { select_parameters->check_matches_mask |= check_none; continue; }
        if (gem_strcaseeq(check,"check-correct"))  { select_parameters->check_matches_mask |= check_correctness; continue; }
        if (gem_strcaseeq(check,"check-best"))     { select_parameters->check_matches_mask |= check_optimum; continue; }
        if (gem_strcaseeq(check,"check-complete")) { select_parameters->check_matches_mask |= check_completness; continue; }
        gem_fatal_error_msg("Option '--check-alignments' must be 'check-correct'|'check-best'|'check-complete'");
        check = strtok(NULL,",");
      }
      break;
    }
    /* Misc */
    case 1300: // --profile in {'sum'|'min'|'max'|'mean'|'sample'}
      parameters->misc.profile = true;
      if (optarg) {
        if (gem_strcaseeq(optarg,"SUM")) {
          parameters->misc.profile_reduce_type = reduce_sum;
        } else if (gem_strcaseeq(optarg,"MIN")) {
          parameters->misc.profile_reduce_type = reduce_min;
        } else if (gem_strcaseeq(optarg,"MAX")) {
          parameters->misc.profile_reduce_type = reduce_max;
        } else if (gem_strcaseeq(optarg,"MEAN")) {
          parameters->misc.profile_reduce_type = reduce_mean;
        } else if (gem_strcaseeq(optarg,"SAMPLE")) {
          parameters->misc.profile_reduce_type = reduce_sample;
        } else {
          gem_fatal_error_msg("Option '--profile' must be 'sum'|'min'|'max'|'mean'|'sample'");
        }
      }
      break;
    case 'v': // -v|--verbose in {'quiet'|'user'|'dev'}
      if (optarg) {
        if (gem_strcaseeq(optarg,"quiet")) {
          parameters->misc.verbose_user = false;
          parameters->misc.verbose_dev = false;
        } else if (gem_strcaseeq(optarg,"user")) {
          parameters->misc.verbose_user = true;
          parameters->misc.verbose_dev = false;
        } else if (gem_strcaseeq(optarg,"dev")) {
          parameters->misc.verbose_user = true;
          parameters->misc.verbose_dev = true;
        } else {
          gem_fatal_error_msg("Option '-v|--verbose' must be 'quiet'|'user'|'dev'");
        }
      } else {
        parameters->misc.verbose_user = true;
        parameters->misc.verbose_dev = false;
      }
      break;
    case 'h':
      usage(false);
      exit(1);
    case 'H':
      usage(true);
      exit(1);
    case '?':
    default:
      gem_fatal_error_msg("Option not recognized");
    }
  }
  /*
   * Parameters Check
   */
  // I/O Parameters
  gem_cond_fatal_error_msg(parameters->io.index_file_name==NULL,"Index file required");
  gem_cond_fatal_error_msg(parameters->io.output_file_name!=NULL && parameters->io.input_file_name!=NULL &&
      gem_streq(parameters->io.input_file_name,parameters->io.output_file_name),
      "Input-file and output-file must be different");
  gem_cond_fatal_error_msg(
      parameters->io.input_file_name!=NULL &&
      gem_streq(parameters->io.index_file_name,parameters->io.input_file_name),
      "Index-file and input-file must be different");
  gem_cond_fatal_error_msg(parameters->io.output_file_name!=NULL &&
      gem_streq(parameters->io.index_file_name,parameters->io.output_file_name),
      "Index-file and output-file must be different");
  /*
   * Search Parameters
   */
  /* Mapping strategy (Mapping mode + properties) */
  /* Qualities */
//  uint64_t quality_threshold;
  gem_cond_fatal_error_msg(parameters->search_parameters.quality_threshold > 94,
      "Quality threshold is too high (please consider lowering it or run ignoring qualities)");
  gem_cond_fatal_error_msg(parameters->search_parameters.quality_threshold == 0,
      "Quality threshold is zero (all base-calls will be considered wrong))");
  /* Error Model */
//  gem_cond_fatal_error_msg(parameters->search_parameters.max_search_error >
//      parameters->search_parameters.max_filtering_error,
//      "Option '--max-search-error' is greater than 'max-filtering-error' [Inefficiency]"); // FIXME
  /* Matches search */
  /* Replacements (Regulates the bases that can be replaced/mismatched) */
  /* Alignment Score */
  /*
   * Select Parameters
   */
  /* Single-end Alignment */
  /* Paired-end Alignment */
  /* Alignment Score */
  /* Reporting */
  gem_cond_fatal_error_msg(parameters->select_parameters.min_reported_matches >
      parameters->select_parameters.max_reported_matches,
      "Option '--max-reported-matches' must be greater or equal than 'min-reported-matches'");
  // System
  /*
   * TODO
   *   check path (gem_access(char* const path,const fm_mode mode))
   *   check number of cpus (force to increase it)
   */
  // Free
  string_destroy(getopt_short_string);
  mm_free(getopt_short_string);
  mm_free(getopt_options);
}
/*
 * Main
 */
int main(int argc,char** argv) {
  // Parsing command-line options
  mapper_parameters_t parameters;
  mapper_parameters_set_defaults(&parameters); // Set defaults
  parse_arguments(argc,argv,&parameters); // Parse cmd-line

  // Runtime setup
  gem_timer_t mapper_time;
  const mapper_parameters_cuda_t* const cuda = &parameters.cuda;
  const uint64_t total_threads = (cuda->cuda_enabled) ?
      (cuda->num_generating_threads + cuda->num_selecting_threads + 1) : parameters.system.num_threads + 1;
  gem_runtime_init(total_threads,parameters.system.max_memory,parameters.system.tmp_folder);
  PROF_START(GP_MAPPER_ALL); TIMER_RESTART(&mapper_time);

  // Open Input/Output File(s)
  gem_mapper_open_input(&parameters);
  gem_mapper_open_output(&parameters);

  // Load GEM-Index
  gem_mapper_load_index(&parameters);

  // Launch mapper
  if (!cuda->cuda_enabled) {
    switch (parameters.mapper_type) {
      case mapper_se:
        mapper_SE_run(&parameters);
        break;
      case mapper_pe:
        mapper_PE_run(&parameters);
        break;
      case mapper_graph:
        GEM_NOT_IMPLEMENTED(); // TODO
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  } else {
    switch (parameters.mapper_type) {
      case mapper_se:
        mapper_SE_CUDA_run(&parameters); // SE-CUDA mapping threads (Producer-Consumer)
        break;
      case mapper_pe:
        GEM_NOT_IMPLEMENTED(); // TODO
        break;
      case mapper_graph:
        GEM_NOT_IMPLEMENTED(); // TODO
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  PROF_STOP(GP_MAPPER_ALL); TIMER_STOP(&mapper_time);

  // Profile
  if (parameters.misc.profile) gem_mapper_print_profile(&parameters);

  // CleanUP
  archive_delete(parameters.archive); // Delete archive
  gem_mapper_close_input(&parameters); // Close I/O files
  gem_mapper_close_output(&parameters); // Close I/O files
  gem_runtime_destroy();

  // Display end banner
  gem_cond_log(parameters.misc.verbose_user,
      "[GEMMapper terminated successfully in %lu s.]\n",(uint64_t)TIMER_GET_TOTAL_S(&mapper_time));

  // Done!
  return 0;
}
