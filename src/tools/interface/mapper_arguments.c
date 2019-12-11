/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "tools/interface/mapper_arguments.h"
#include "stats/report_stats.h"

/*
 * Mapper options Menu
 */
option_t gem_mapper_options[] = {
  /* I/O */
  { 'I', "index", REQUIRED, TYPE_STRING, 2, VISIBILITY_USER, "<index_file.gem>", "" },
  { 'i', "input", REQUIRED, TYPE_STRING, 2, VISIBILITY_USER, "<file>", "(FASTA/FASTQ, default=stdin)" },
  { '1', "i1", REQUIRED, TYPE_STRING, 2, VISIBILITY_USER, "<file>", "(paired-end, end-1)" },
  { '2', "i2", REQUIRED, TYPE_STRING, 2, VISIBILITY_USER, "<file>", "(paired-end, end-2)" },
  { 'z', "gzip-input", NO_ARGUMENT, TYPE_NONE, 2, VISIBILITY_USER, "", "(gzip input)" },
  { 'j', "bzip-input", NO_ARGUMENT, TYPE_NONE, 2, VISIBILITY_USER, "", "(bzip input)" },
  { 201, "input-model", REQUIRED, TYPE_STRING, 2, VISIBILITY_DEVELOPER, "<block_size,num_blocks,buffer_size>", "(default=32M,c,4M)" },
  { 'o', "output", REQUIRED, TYPE_STRING, 2, VISIBILITY_USER, "<output_prefix>" , "(default=stdout)" },
  { 202, "gzip-output", NO_ARGUMENT, TYPE_NONE, 2, VISIBILITY_USER, "", "(gzip output)" },
  { 203, "bzip-output", NO_ARGUMENT, TYPE_NONE, 2, VISIBILITY_USER, "", "(bzip output)" },
  { 204, "output-model", REQUIRED, TYPE_STRING, 2, VISIBILITY_DEVELOPER, "<buffer_size,num_buffers>", "(default=4M,5c)" },
  { 205, "report-file", REQUIRED, TYPE_STRING, 2, VISIBILITY_USER, "<file_name>", "(default=disabled)" },
  { 206, "clipping", OPTIONAL, TYPE_STRING, 2, VISIBILITY_ADVANCED, "'none'|'uncalled'|'masked'|'fixed,<left_clip>,<right_clip>'", "(default=uncalled)" },
  { 207, "benchmark-mode", NO_ARGUMENT, TYPE_NONE, 2, VISIBILITY_ADVANCED, "", "" },
  { 'q', "quality-format", REQUIRED, TYPE_STRING, 2, VISIBILITY_ADVANCED, "'ignore'|'offset-33'|'offset-64'", "(default=offset-33)" },
  /* Single-end Alignment */
  { 300, "mapping-mode", REQUIRED, TYPE_STRING, 3, VISIBILITY_USER, "'fast'|'sensitive'|'customed'" , "(default=fast)" },
  { 'E', "complete-search-error", REQUIRED, TYPE_FLOAT, 3, VISIBILITY_ADVANCED, "<number|percentage>" , "(default=0.04, 4%)" },
  { 's', "complete-strata-after-best", REQUIRED, TYPE_FLOAT, 3, VISIBILITY_ADVANCED, "<number|percentage>" , "(default=1)" },
  { 'e', "alignment-max-error", REQUIRED, TYPE_FLOAT, 3, VISIBILITY_USER, "<number|percentage>" , "(default=0.12, 12%)" },
  { 301, "alignment-max-bandwidth", REQUIRED, TYPE_FLOAT, 3, VISIBILITY_ADVANCED, "<number|percentage>" , "(default=0.20, 20%)" },
  { 303, "alignment-global-min-identity", REQUIRED, TYPE_FLOAT, 3, VISIBILITY_USER, "<number|percentage>" , "(default=80%)" },
  { 304, "alignment-global-min-score", REQUIRED, TYPE_FLOAT, 3, VISIBILITY_USER, "<number|percentage>" , "(default=0.20)" },
  { 305, "alignment-local", REQUIRED, TYPE_FLOAT, 3, VISIBILITY_USER, "'never'|'if-unmapped'" , "(default=if-unmapped)" },
  { 306, "alignment-local-min-identity", REQUIRED, TYPE_FLOAT, 3, VISIBILITY_USER, "<number|percentage>" , "(default=40)" },
  { 307, "alignment-local-min-score", REQUIRED, TYPE_FLOAT, 3, VISIBILITY_USER, "<number|percentage>" , "(default=20)" },
  { 308, "alignment-force-full-swg", OPTIONAL, TYPE_STRING, 3, VISIBILITY_ADVANCED, "" , "(default=false)" },
  { 309, "alignment-scaffolding-min-matching_length", REQUIRED, TYPE_FLOAT, 3, VISIBILITY_DEVELOPER, "<number|percentage>" , "(default=10)" },
  { 310, "alignment-curation", OPTIONAL, TYPE_STRING, 3, VISIBILITY_ADVANCED, "" , "(default=true)" },
  { 311, "alignment-curation-min-end-context", REQUIRED, TYPE_FLOAT, 3, VISIBILITY_DEVELOPER, "<number|percentage>" , "(default=2)" },
  { 312, "candidate-generation", REQUIRED, TYPE_STRING, 3, VISIBILITY_DEVELOPER, "<strategy>[,<arguments>]" , "" },
  { 313, "candidate-generation-adaptive", REQUIRED, TYPE_STRING, 3, VISIBILITY_DEVELOPER, "<app_threshold>,<app_steps>,<app_dec>" , "" },
  { 314, "candidate-verification", REQUIRED, TYPE_STRING, 3, VISIBILITY_DEVELOPER, "'BPM'|'chained'" , "" },
  { 315, "qgram-filter", REQUIRED, TYPE_STRING, 3, VISIBILITY_DEVELOPER, "<kmer_tiles>,<qgram_length>" , "" },
  /* Paired-end Alignment */
  { 'p', "paired-end-alignment", NO_ARGUMENT, TYPE_NONE, 4, VISIBILITY_USER, "" , "" },
  { 'l', "min-template-length", REQUIRED, TYPE_INT, 4, VISIBILITY_USER, "<number>" , "(default=disabled)" },
  { 'L', "max-template-length", REQUIRED, TYPE_INT, 4, VISIBILITY_USER, "<number>" , "(default=10000)" },
  { 400, "discordant-pair-search", REQUIRED, TYPE_STRING, 4, VISIBILITY_USER, "'always'|'if-no-concordant'|'never'" , "(default=if-no-concordant)" },
  { 401, "pair-orientation", REQUIRED, TYPE_STRING, 4, VISIBILITY_ADVANCED, "'FR'|'RF'|'FF'|'RR'" , "(default=FR)" },
  { 402, "discordant-pair-orientation", REQUIRED, TYPE_STRING, 4, VISIBILITY_ADVANCED, "'FR'|'RF'|'FF'|'RR'" , "(default=RF)" },
  { 403, "pair-layout", REQUIRED, TYPE_STRING, 4, VISIBILITY_ADVANCED, "'separate'|'overlap'|'contain'" , "(default=separated,overlap)" },
  { 404, "discordant-pair-layout", REQUIRED, TYPE_STRING, 4, VISIBILITY_ADVANCED, "'separate'|'overlap'|'contain'" , "(default=contain)" },
  { 406, "pe-template-length", REQUIRED, TYPE_STRING, 4, VISIBILITY_ADVANCED, "<min>,<max>,<samples>" , "(default=0,800,100)" },
  /* Bisulfite and Hi-C Alignment */
  { 500, "bisulfite-conversion", REQUIRED, TYPE_STRING, 5, VISIBILITY_USER, "'inferred-C2T-G2A','inferred-G2A-C2T','C2T','G2A','non-stranded'",  "(default=inferred-C2T-G2A)" },
  { 501, "underconversion-sequence", REQUIRED, TYPE_STRING, 5, VISIBILITY_USER, "<sequence name>",  "(default=" UNDERCONVERSION_CONTROL ")" },
  { 502, "overconversion-sequence", REQUIRED, TYPE_STRING, 5, VISIBILITY_USER, "<sequence name>",  "(default=" OVERCONVERSION_CONTROL ")" },
  { 503, "control-sequence", REQUIRED, TYPE_STRING, 5, VISIBILITY_USER, "<sequence name>",  "(default=" SEQUENCING_CONTROL ")" },
  { 504, "restriction-site", REQUIRED, TYPE_STRING, 5, VISIBILITY_ADVANCED, "<restriction site> (i.e., 'C-CGG')", "(default = NULL)" },
  { 505, "rrbs", NO_ARGUMENT, TYPE_NONE, 5, VISIBILITY_ADVANCED, "", "" },
  /* Alignment Score */
  { 600, "alignment-model", REQUIRED, TYPE_STRING, 6, VISIBILITY_ADVANCED, "'pseudoalignment'|'hamming'|'edit'|'gap-affine'" , "(default=gap-affine)" },
  { 601, "gap-affine-penalties", REQUIRED, TYPE_STRING, 6, VISIBILITY_USER, "A,B,O,X" , "(default=1,4,6,1)" },
  { 'A', "matching-score", REQUIRED, TYPE_INT, 6, VISIBILITY_ADVANCED, "" , "(default=1)" },
  { 'B', "mismatch-penalty", REQUIRED, TYPE_INT, 6, VISIBILITY_ADVANCED, "" , "(default=4)" },
  { 'O', "gap-open-penalty", REQUIRED, TYPE_INT, 6, VISIBILITY_ADVANCED, "" , "(default=6)" },
  { 'X', "gap-extension-penalty", REQUIRED, TYPE_INT, 6, VISIBILITY_ADVANCED, "" , "(default=1)" },
  /* MAQ Score */
  { 700, "mapq-model", REQUIRED, TYPE_STRING, 7, VISIBILITY_ADVANCED, "'none'|'gem'" , "(default=gem)" },
  /* Reporting */
  { 'm', "min-reported-strata", REQUIRED, TYPE_FLOAT, 8, VISIBILITY_ADVANCED, "<number|percentage>|'all'" , "(stratum-wise, default=0)" },
  { 'M', "max-reported-matches", REQUIRED, TYPE_INT,  8, VISIBILITY_USER, "<number>|'all'" , "(default=5)" },
  /* Output Format */
  { 'F',  "output-format", REQUIRED, TYPE_STRING, 9, VISIBILITY_USER, "'MAP'|'SAM'" , "(default=SAM)" },
  { 900, "sam-compact", OPTIONAL, TYPE_STRING, 9, VISIBILITY_USER, "'true'|'false'" , "(default=true)" },
  { 'r',  "sam-read-group-header", REQUIRED, TYPE_STRING, 9, VISIBILITY_USER, "<read_group_header> (i.e. '@RG\\tID:xx\\tSM:yy')" , "(default=NULL)" },
  { 901, "sam-gem-compatible", OPTIONAL, TYPE_STRING, 9, VISIBILITY_DEVELOPER, "'true'|'false'" , "(default=false)" },
  { 902, "map-format", REQUIRED, TYPE_INT, 9, VISIBILITY_DEVELOPER, "'1'|'2'|'3'" , "(default=2)" },
  /* System */
  { 't',  "threads", REQUIRED, TYPE_STRING, 10, VISIBILITY_USER, "<number>" , "(default=#cores)" },
  { 1000, "max-memory", REQUIRED, TYPE_STRING, 10, VISIBILITY_DEVELOPER, "<maximum-memory>" , "(Eg 2GB)" },
  { 1001, "tmp-folder", REQUIRED, TYPE_STRING, 10, VISIBILITY_DEVELOPER, "<temporal_dir_path>" , "(default=/tmp/)" },
#ifdef HAVE_CUDA
  /* CUDA Settings */
  { 1100, "gpu", OPTIONAL, TYPE_STRING, 11, VISIBILITY_USER, "", "(default=disabled)"},
  { 1101, "gpu-devices", REQUIRED, TYPE_STRING, 11, VISIBILITY_USER, "", "(default=all)"},
  { 1102, "gpu-buffers-model", REQUIRED, TYPE_STRING, 11, VISIBILITY_DEVELOPER, "<#BSearch,#BDecode,#BKmer,#BBPMDistance,#BBPMAlign,BufferSize>" , "(default=2,3,3,3,3,1M)" },
#endif /* HAVE_CUDA */
  /* Miscellaneous */
  { 'c',  "check-alignments", REQUIRED, TYPE_STRING, 12, VISIBILITY_DEVELOPER, "'correct'|'best'|'complete'" , "" },
  { 1200, "profile", OPTIONAL, TYPE_STRING, 12, VISIBILITY_DEVELOPER, "'sum'|'min'|'max'|'mean'|'sample'" , "(disabled)" },
  { 'v',  "verbose", OPTIONAL, TYPE_STRING, 12, VISIBILITY_USER, "'quiet'|'user'|'dev'" , "(default=user)" },
  { 1201, "progress-step", REQUIRED, TYPE_INT, 12, VISIBILITY_DEVELOPER, "" , "(default=100000)" },
  { 1202, "version", NO_ARGUMENT, TYPE_STRING, 12, VISIBILITY_USER, "" , "" },
  { 'h',  "help", OPTIONAL, TYPE_NONE, 12, VISIBILITY_USER, "" , "(print usage)" },
  {  0, "", 0, 0, 0, false, "", ""}
};
char* gem_mapper_groups[] = {
  /*  0 */ "Null",
  /*  1 */ "Unclassified",
  /*  2 */ "I/O",
  /*  3 */ "Single-end Alignment",
  /*  4 */ "Paired-end Alignment",
  /*  5 */ "Bisulfite and Hi-C Alignment",
  /*  6 */ "Alignment Score",
  /*  7 */ "MAPQ Score",
  /*  8 */ "Reporting",
  /*  9 */ "Output-format",
  /* 10 */ "System",
  /* 11 */ "CUDA Settings",
  /* 12 */ "Miscellaneous",
};

/*
 * Mapper Parsing
 */
int parse_integer_system(char* const argument,uint64_t* const value) {
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
/*
 * Mapper Usage
 */
void gem_mapper_print_usage(const option_visibility_t visibility_level) {
  fprintf(stderr, "USAGE: ./gem-mapper [ARGS]...\n");
  options_fprint_menu(stderr,gem_mapper_options,gem_mapper_groups,true,visibility_level);
}
/*
 * Mapper Parameters Checks
 */
void gem_mapper_parameters_check(mapper_parameters_t* const parameters) {
  // Parameters
  search_parameters_t* const search = &parameters->search_parameters;
  search_paired_parameters_t* const paired_search = &search->search_paired_parameters;
  // I/O Parameters
  const char* const pindex = parameters->io.index_file_name;
  const char* const poutput = parameters->io.output_file_name;
  mapper_cond_error_msg(pindex==NULL,"Index file required");
  if (!parameters->io.separated_input_files) {
    const char* const pinput = parameters->io.input_file_name;
    if (pinput!=NULL) {
      mapper_cond_error_msg(gem_streq(pindex,pinput), "Index and Input-file must be different");
      if (poutput!=NULL) {
        mapper_cond_error_msg(gem_streq(pinput,poutput),"Input-file and Output-file must be different");
        mapper_cond_error_msg(gem_streq(pindex,poutput),"Index and Output-file must be different");
      }
    }
  } else {
    const char* const pinput_1 = parameters->io.input_file_name_end1;
    const char* const pinput_2 = parameters->io.input_file_name_end2;
    mapper_cond_error_msg(pinput_1==NULL, "Missing Input-End1 (--i1)");
    mapper_cond_error_msg(pinput_2==NULL, "Missing Input-End2 (--i2)");
    mapper_cond_error_msg(gem_streq(pinput_1,pinput_2), "Input-End1 and Input-End2 must be different");
    mapper_cond_error_msg(gem_streq(pindex,pinput_1), "Index and Input-End1 must be different");
    mapper_cond_error_msg(gem_streq(pindex,pinput_2), "Index and Input-End2 must be different");
    if (poutput!=NULL) {
      mapper_cond_error_msg(gem_streq(pinput_1,poutput),"Input-End1 and Output-file must be different");
      mapper_cond_error_msg(gem_streq(pinput_2,poutput),"Input-End2 and Output-file must be different");
      mapper_cond_error_msg(gem_streq(pindex,poutput),"Index and Output-file must be different");
    }
  }
  /*
   * Search Parameters
   */
  /* Mapping strategy (Mapping mode + properties) */
  if (paired_search->paired_end_search) {
    parameters->mapper_type = mapper_pe;
    gem_cond_warn_msg((search->bisulfite_read != bisulfite_inferred_C2T_G2A) && (search->bisulfite_read != bisulfite_inferred_G2A_C2T) && (search->bisulfite_read != bisulfite_non_stranded),
        "Option '--bisulfite_read' can only be set in paired end mode to {'inferred-C2T-G2A'|'inferred-G2A-C2T'|'non-stranded'}");
    search->mapq_model_pe = search->mapq_model_se;
    search->mapq_model_se = mapq_model_gem;
  } else {
    parameters->mapper_type = mapper_se;
  }
  if (parameters->cuda.gpu_enabled) {
    if (search->mapping_mode != mapping_adaptive_filtering_fast &&
        search->mapping_mode != mapping_hybrid_sensitive) {
      mapper_error_msg("Option '--gpu' can only operate with '--mapping-mode' in {'fast'|'sensitive'}");
    }
  }
  /* Reporting */
  mapper_cond_error_msg(
      parameters->min_reported_strata_set && parameters->max_reported_matches_set,
      "Options '--max-reported-matches' and 'min-reported-strata' cannot be used at the same time");
  mapper_cond_error_msg(
      search->select_parameters.max_reported_matches == 0,
      "Option '--max-reported-matches' must be greater than zero'");
  /* RRBS */
  if(search->rrbs) {
    if(search->restriction_sites == NULL) {
      restriction_t *rest = restriction_new("C-CGG"); // If no restriction site has been set, default to MspI site
      if(rest != NULL) {
        search->restriction_sites = vector_new(1, restriction_t *);
        vector_insert(search->restriction_sites, rest, restriction_t *);
      }
    }
  }
}
/*
 * Mapper Arguments Parsing
 */
bool gem_mapper_parse_arguments_io(
    mapper_parameters_t* const parameters,
    const int option,
    char* optarg) {
  // Parameters
  mapper_parameters_io_t* const io = &parameters->io;
  search_parameters_t* const search = &parameters->search_parameters;
  search_paired_parameters_t* const paired_search = &search->search_paired_parameters;
  mapper_parameters_cuda_t* const cuda = &parameters->cuda;
  // I/O
  switch (option) {
    case 'I': // --index
      io->index_file_name = optarg;
      return true;
    case 'i': // --input
      io->separated_input_files = false;
      io->input_file_name = optarg; // TODO Multiple input files
      return true;
    case '1': // --i1
      paired_search->paired_end_search = true;
      io->separated_input_files = true;
      io->input_file_name_end1 = optarg;
      return true;
    case '2': // --i2
      paired_search->paired_end_search = true;
      io->separated_input_files = true;
      io->input_file_name_end2 = optarg;
      return true;
    case 'z': // --gzip-input
      io->input_compression = FM_GZIPPED_FILE;
      return true;
    case 'j': // --bzip-input
      io->input_compression = FM_BZIPPED_FILE;
      return true;
    case 201: { // --input-model=64M,2c,5K
      char *block_size = NULL, *num_blocks = NULL, *buffer_size = NULL;
      const int num_arguments = input_text_parse_csv_arguments(optarg, 3,
          &block_size, &num_blocks, &buffer_size);
      mapper_cond_error_msg(num_arguments != 3,
          "Option '--input-model' wrong number of arguments");
      // Parse input-buffer size
      int error = input_text_parse_size(block_size, &io->input_block_size);
      mapper_cond_error_msg(error,
          "Option '--input-model'. Error parsing 'block_size'");
      // Parse number of buffers
      error = parse_integer_system(num_blocks, &io->input_num_blocks);
      mapper_cond_error_msg(error,
          "Option '--input-model'. Error parsing 'num_blocks'");
      // Parse number of records per buffer
      error = input_text_parse_size(buffer_size, &io->input_buffer_size);
      mapper_cond_error_msg(error,
          "Option '--input-model'. Error parsing 'buffer_size'");
      // Propagate settings to CUDA
      cuda->input_block_size = io->input_block_size;
      cuda->input_buffer_size = io->input_buffer_size;
      return true;
    }
    case 'o': // --output
      io->output_file_name = optarg;
      return true;
    case 202: // --gzip-output
      io->output_compression = FM_GZIPPED_FILE;
      return true;
    case 203: // --bzip-output
      io->output_compression = FM_BZIPPED_FILE;
      return true;
    case 204: { // --output-model=4M,4c
      char *buffer_size = NULL, *num_buffers = NULL;
      const int num_arguments = input_text_parse_csv_arguments(optarg, 2,
          &buffer_size, &num_buffers);
      mapper_cond_error_msg(num_arguments != 2,
          "Option '--output-model' wrong number of arguments");
      // Parse output-buffer size
      int error = input_text_parse_size(buffer_size, &io->output_buffer_size);
      mapper_cond_error_msg(error,
          "Option '--output-model'. Error parsing 'buffer_size'");
      // Parse number of output-buffers
      error = parse_integer_system(num_buffers, &io->output_num_buffers);
      mapper_cond_error_msg(error,
          "Option '--output-model'. Error parsing 'num_buffers'");
      // Propagate settings to CUDA
      cuda->output_buffer_size = io->output_buffer_size;
      cuda->output_num_buffers = io->output_num_buffers;
      return true;
    }
    case 205: // --report-file
      io->report_file_name = optarg;
      return true;
    case 206: // --clipping 'none'|'uncalled'|'masked'|'fixed,<left_clip>,<right_clip>' (default=uncalled)
      if (optarg == NULL) {
        search->clipping = clipping_uncalled;
      } else {
        // Parse optargs
        char *clipping_mode = NULL, *left_clip = NULL, *right_clip = NULL;
        const int num_arguments = input_text_parse_csv_arguments(optarg, 3,
            &clipping_mode, &left_clip, &right_clip);
        if (num_arguments == 1 && gem_strcaseeq(clipping_mode, "none")) {
          search->clipping = clipping_disabled;
        } else if (num_arguments == 1
            && gem_strcaseeq(clipping_mode, "uncalled")) {
          search->clipping = clipping_uncalled;
        } else if (num_arguments == 1 && gem_strcaseeq(clipping_mode, "masked")) {
          search->clipping = clipping_masked;
        } else if (num_arguments == 3 && gem_strcaseeq(clipping_mode, "fixed")) {
          search->clipping = clipping_fixed;
          input_text_parse_extended_uint64(left_clip, &search->clip_left);
          input_text_parse_extended_uint64(right_clip, &search->clip_right);
        } else {
          mapper_error_msg(
              "Option '--clipping' must be 'none', 'uncalled', 'masked' or 'fixed,<left_clip>,<right_clip>'");
        }
      }
      return true;
    case 207: // benchmark-mode
    	parameters->io.sam_parameters.benchmark_mode = true;
    	return true;
    case 'q': // --quality-format
      if (gem_strcaseeq(optarg, "ignore")) {
        search->qualities_format = sequence_qualities_ignore;
      } else if (gem_strcaseeq(optarg, "offset-33")) {
        search->qualities_format = sequence_qualities_offset_33;
      } else if (gem_strcaseeq(optarg, "offset-64")) {
        search->qualities_format = sequence_qualities_offset_64;
      } else {
        mapper_error_msg(
            "Option '-q|--quality-format' must be 'ignore', 'offset-33' or 'offset-64'");
      }
      return true;
    default:
      return false; // Not found
  }
  // Not found
  return false;
}
bool gem_mapper_parse_arguments_single_end(
    mapper_parameters_t* const parameters,
    const int option,
    char* optarg) {
  // Parameters
  search_parameters_t* const search = &parameters->search_parameters;
  // Single-End
  switch (option) {
    case 300: // --mapping-mode in {'fast'|'sensitive'|'customed'} (default=fast)
      // Filtering Modes
      if (gem_strcaseeq(optarg,"fast")) {
        search->mapping_mode = mapping_adaptive_filtering_fast;
      // NS Modes
      } else if (gem_strcaseeq(optarg,"complete-brute-force")) {
        search->mapping_mode = mapping_neighborhood_search_brute_force;
      } else if (gem_strcaseeq(optarg,"complete-partition")) {
        search->nsearch_parameters.filtering_cutoff = false;
        search->nsearch_parameters.matches_accuracy_cutoff = false;
        search->nsearch_parameters.matches_max_searched_cutoff = false;
        search->mapping_mode = mapping_neighborhood_search_partition;
      } else if (gem_strcaseeq(optarg,"complete-partition-filter")) {
        search->nsearch_parameters.filtering_cutoff = true;
        search->nsearch_parameters.matches_accuracy_cutoff = false;
        search->nsearch_parameters.matches_max_searched_cutoff = false;
        search->mapping_mode = mapping_neighborhood_search_partition;
      // Hybrid Modes
      } else if (gem_strcaseeq(optarg,"sensitive") ||
                 gem_strcaseeq(optarg,"hybrid")) {
        search->mapping_mode = mapping_hybrid_sensitive;
      // Complete Modes
      } else if (gem_strcaseeq(optarg,"complete-pure")) {
        search->nsearch_parameters.filtering_cutoff = false;
        search->nsearch_parameters.matches_accuracy_cutoff = false;
        search->nsearch_parameters.matches_max_searched_cutoff = false;
        search->mapping_mode = mapping_hybrid_complete;
      } else if (gem_strcaseeq(optarg,"customed") ||
                 gem_strcaseeq(optarg,"complete")) {
        search->nsearch_parameters.matches_accuracy_cutoff = false;
        search->nsearch_parameters.matches_max_searched_cutoff = false;
        search->mapping_mode = mapping_hybrid_complete;
      } else {
        mapper_error_msg("Option '--mapping-mode' must be 'fast'|'sensitive'|'customed'");
      }
      return true;
    case 'E': // --complete-search-error (default=0.04, 4%)
      search->complete_search_error = atof(optarg);
      return true;
    case 's': // --complete-strata-after-best (default=1)
      input_text_parse_extended_double(optarg,(double*)&search->complete_strata_after_best);
      return true;
    case 'e': // --alignment-max-error (default=0.12, 12%)
      search->alignment_max_error = atof(optarg);
      return true;
    case 301: // --alignment-max-bandwidth (default=0.20, 20%)
      input_text_parse_extended_double(optarg,(double*)&search->alignment_max_bandwidth);
      return true;
    case 303: // --alignment-global-min-identity (default=80%)
      input_text_parse_extended_double(optarg,(double*)&search->alignment_global_min_identity);
      return true;
    case 304: // --alignment-global-min-score (default=40%)
      input_text_parse_extended_double(optarg,(double*)&search->alignment_global_min_swg_threshold);
      return true;
    case 305: // --alignment-local in {'never'|'if-unmapped'} (default=if-unmapped)
      if (gem_strcaseeq(optarg,"never")) {
        search->alignment_local = local_alignment_never;
      } else if (gem_strcaseeq(optarg,"if-unmapped")) {
        search->alignment_local = local_alignment_if_unmapped;
      } else {
        mapper_error_msg("Option '--alignment-local' must be 'never'|'if-unmapped'");
      }
      return true;
    case 306: // --alignment-local-min-identity (default=40)
      input_text_parse_extended_double(optarg,(double*)&search->alignment_local_min_identity);
      return true;
    case 307: // --alignment-local-min-score (default=20)
      input_text_parse_extended_double(optarg,(double*)&search->alignment_local_min_swg_threshold);
      return true;
    case 308: // --alignment-force-full-swg (default=true)
      search->alignment_force_full_swg = input_text_parse_extended_bool(optarg);
      return true;
    case 309: // --alignment-scaffolding-min-matching_length (default=10)
      input_text_parse_extended_double(optarg,(double*)&search->alignment_scaffolding_min_matching_length);
      return true;
    case 310: // --alignment-curation (default=true)
      search->cigar_curation = input_text_parse_extended_bool(optarg);
      return true;
    case 311: // --alignment-curation-min-end-context (default=2)
      input_text_parse_extended_double(optarg,(double*)&search->cigar_curation_min_end_context);
      return true;
    case 312: { // --candidate-generation <strategy>,[...]
      char *strategy=NULL, *num_regions=NULL, *region_length=NULL;
      char *region_step=NULL, *region_error=NULL, *max_candidates=NULL;
      if (strncasecmp(optarg,"fixed,",6)==0) {
        // --candidate-generation fixed,<region_length>,<region_step>,<region_error>,<max_candidates>
        const int num_arguments = input_text_parse_csv_arguments(
            optarg,5,&strategy,&region_length,&region_step,&region_error,&max_candidates);
        mapper_cond_error_msg(num_arguments!=5,
            "Invalid usage '--candidate-generation fixed,<region_length>,<region_step>,<region_error>,<max_candidates>'");
        search->region_profile_model.strategy = region_profile_fixed;
        input_text_parse_extended_uint64(region_length,&search->region_profile_model.region_length);
        input_text_parse_extended_uint64(region_step,&search->region_profile_model.region_step);
        input_text_parse_extended_uint64(region_error,&search->region_profile_model.region_error);
        input_text_parse_extended_uint64(max_candidates,&search->region_profile_model.max_candidates);
      } else if (strncasecmp(optarg,"CKS,",4)==0) {
        // --candidate-generation CKS,<num_regions>,<region_length>,<max_candidates>
        const int num_arguments = input_text_parse_csv_arguments(
            optarg,4,&strategy,&num_regions,&region_length,&max_candidates);
        mapper_cond_error_msg(num_arguments!=4,
            "Invalid usage '--candidate-generation CKS,<num_regions>,<region_length>,<max_candidates>'");
        search->region_profile_model.strategy = region_profile_CKS;
        input_text_parse_extended_uint64(num_regions,&search->region_profile_model.num_regions);
        input_text_parse_extended_uint64(region_length,&search->region_profile_model.region_length);
        input_text_parse_extended_uint64(max_candidates,&search->region_profile_model.max_candidates);
      } else if (strncasecmp(optarg,"OPS,",4)==0) {
        // --candidate-generation OPS,<num_regions>,<region_length>,<max_candidates>
        const int num_arguments = input_text_parse_csv_arguments(
            optarg,4,&strategy,&num_regions,&region_length,&max_candidates);
        mapper_cond_error_msg(num_arguments!=4,
            "Invalid usage '--candidate-generation OPS,<num_regions>,<region_length>,<max_candidates>'");
        search->region_profile_model.strategy = region_profile_OPS;
        input_text_parse_extended_uint64(num_regions,&search->region_profile_model.num_regions);
        input_text_parse_extended_uint64(region_length,&search->region_profile_model.region_length);
        input_text_parse_extended_uint64(max_candidates,&search->region_profile_model.max_candidates);
      } else if (strncasecmp(optarg,"factors,",8)==0) {
        // --candidate-generation factors,<num_regions>,<region_error>,<max_candidates>
        const int num_arguments = input_text_parse_csv_arguments(
            optarg,4,&strategy,&num_regions,&region_error,&max_candidates);
        mapper_cond_error_msg(num_arguments!=4,
            "Invalid usage '--candidate-generation factors,<num_regions>,<region_error>,<max_candidates>'");
        search->region_profile_model.strategy = region_profile_factor;
        input_text_parse_extended_uint64(num_regions,&search->region_profile_model.num_regions);
        input_text_parse_extended_uint64(region_error,&search->region_profile_model.region_error);
        input_text_parse_extended_uint64(max_candidates,&search->region_profile_model.max_candidates);
      } else if (strncasecmp(optarg,"adaptive,",9)==0) {
        // --candidate-generation adaptive
        mapper_cond_error_msg(strcasecmp(optarg,"adaptive")!=0,"Invalid usage '--candidate-generation adaptive'");
        search->region_profile_model.strategy = region_profile_adaptive;
      } else if (strncasecmp(optarg,"adaptive-limited,",17)==0) {
        // --candidate-generation adaptive-limited
        mapper_cond_error_msg(strcasecmp(optarg,"adaptive-limited")!=0,
            "Invalid usage '--candidate-generation adaptive-limited'");
        search->region_profile_model.strategy = region_profile_adaptive_limited;
      } else if (strncasecmp(optarg,"MEM",3)==0) {
        // --candidate-generation MEM
        const int num_arguments = input_text_parse_csv_arguments(optarg,2,&strategy,&max_candidates);
        mapper_cond_error_msg(num_arguments!=2,"Invalid usage '--candidate-generation MEM,<max_candidates>'");
        search->region_profile_model.strategy = region_profile_MEM;
        input_text_parse_extended_uint64(max_candidates,&search->region_profile_model.max_candidates);
      } else if (strncasecmp(optarg,"SMEM",4)==0) {
        // --candidate-generation SMEM
        const int num_arguments = input_text_parse_csv_arguments(optarg,2,&strategy,&max_candidates);
        mapper_cond_error_msg(num_arguments!=2,"Invalid usage '--candidate-generation SMEM,<max_candidates>'");
        search->region_profile_model.strategy = region_profile_SMEM;
        input_text_parse_extended_uint64(max_candidates,&search->region_profile_model.max_candidates);
      } else if (strncasecmp(optarg,"OPP",3)==0) {
        // --candidate-generation OPP,<num_regions>
        const int num_arguments = input_text_parse_csv_arguments(optarg,2,&strategy,&num_regions);
        mapper_cond_error_msg(num_arguments!=2,"Invalid usage '--candidate-generation OPP,<num_regions>'");
        search->region_profile_model.strategy = region_profile_OPP;
        input_text_parse_extended_uint64(num_regions,&search->region_profile_model.num_regions);
      } else if (strncasecmp(optarg,"test,",4)==0) {
        // --candidate-generation <strategy>,<num_regions>,<region_length>,<region_step>,<region_error>,<max_candidates>
        const int num_arguments = input_text_parse_csv_arguments(optarg,6,
            &strategy,&num_regions,&region_length,&region_step,&region_error,&max_candidates);
        mapper_cond_error_msg(num_arguments!=6,
            "Invalid usage '--candidate-generation test,<num_regions>,<region_length>,<region_step>,<region_error>,<max_candidates>'");
        search->region_profile_model.strategy = region_profile_test;
        input_text_parse_extended_uint64(num_regions,&search->region_profile_model.num_regions);
        input_text_parse_extended_uint64(region_length,&search->region_profile_model.region_length);
        input_text_parse_extended_uint64(region_step,&search->region_profile_model.region_step);
        input_text_parse_extended_uint64(region_error,&search->region_profile_model.region_error);
        input_text_parse_extended_uint64(max_candidates,&search->region_profile_model.max_candidates);
      } else {
        mapper_error_msg(
            "Invalid <strategy> from option '--candidate-generation'. "
            "Select from 'fixed','CKS','OPS','factors',"
            "'adaptive','adaptive-limited','MEM','SMEM','OPP','test'");
      }
      return true;
    }
    case 313: { // --candidate-generation-adaptive <app_threshold>,<app_steps>,<app_dec>
      char *region_th=NULL, *max_steps=NULL, *dec_factor=NULL;
      const int num_arguments = input_text_parse_csv_arguments(
          optarg,3,&region_th,&max_steps,&dec_factor);
      mapper_cond_error_msg(num_arguments!=3,
          "Option '--candidate-generation-adaptive' wrong arguments (<app_threshold>,<app_steps>,<app_dec>)");
      input_text_parse_extended_uint64(region_th,&search->region_profile_model.region_th);
      input_text_parse_extended_uint64(max_steps,&search->region_profile_model.max_steps);
      input_text_parse_extended_uint64(dec_factor,&search->region_profile_model.dec_factor);
      search->region_profile_model.strategy = region_profile_adaptive;
      return true;
    }
    case 314: // --candidate-verification in 'BPM'|'chained'
      if (gem_strcaseeq(optarg,"BPM")) {
        search->candidate_verification.verification_strategy = verification_BPM;
      } else if (gem_strcaseeq(optarg,"chained")) {
        search->candidate_verification.verification_strategy = verification_chained;
      } else {
        mapper_error_msg("Option '--candidate-verification' must be 'BPM'|'chained'");
      }
      return true;
    case 315: { // --qgram-filter <num_slices>,<qgram_length>
      char *kmer_tiles=NULL, *qgram_length=NULL;
      const int num_arguments = input_text_parse_csv_arguments(optarg,2,&kmer_tiles,&qgram_length);
      mapper_cond_error_msg(num_arguments!=2,"Option '--qgram-filter' wrong arguments (<kmer_tiles>,<qgram_length>)");
      input_text_parse_extended_uint64(kmer_tiles,&search->candidate_verification.kmer_tiles);
      input_text_parse_extended_uint64(qgram_length,&search->candidate_verification.kmer_length);
      return true;
    }
    case 316: // --candidate-drop-off in {'true'|'false'} (default=true)
      search->candidate_verification.candidate_local_drop_off = (optarg==NULL) ? true : input_text_parse_extended_bool(optarg);
      return true;
    default:
      return false; // Not found
  }
  // Not found
  return false;
}
bool gem_mapper_parse_arguments_paired_end(
    mapper_parameters_t* const parameters,
    const int option,
    char* optarg) {
  // Parameters
  search_parameters_t* const search = &parameters->search_parameters;
  search_paired_parameters_t* const paired_search = &search->search_paired_parameters;
  // Paired-End
  switch (option) {
    case 'p': // --paired-end-alignment
      paired_search->paired_end_search = true;
      return true;
    case 'l': // --min-template-length (default=disabled)
      input_text_parse_extended_uint64(optarg,&paired_search->min_template_length);
      return true;
    case 'L': // --max-template-length (default=10000)
      input_text_parse_extended_uint64(optarg,&paired_search->max_template_length);
      return true;
    case 400: // --discordant-pair-search in 'always'|'if-no-concordant'|'never' (default=if-no-concordant)
      if (gem_strcaseeq(optarg,"always")) {
        paired_search->pair_discordant_search = pair_discordant_search_always;
      } else if (gem_strcaseeq(optarg,"if-no-concordant")) {
        paired_search->pair_discordant_search = pair_discordant_search_only_if_no_concordant;
      } else if (gem_strcaseeq(optarg,"never")) {
        paired_search->pair_discordant_search = pair_discordant_search_never;
      } else {
        mapper_error_msg("Option '--search-discordant' must be 'always'|'if-no-concordant'|'never'");
      }
      return true;
    case 401: { // --pair-orientation in {'FR'|'RF'|'FF'|'RR'} (default=FR)
      // Init null
      paired_search->pair_orientation[pair_orientation_FR] = pair_relation_invalid;
      paired_search->pair_orientation[pair_orientation_RF] = pair_relation_invalid;
      paired_search->pair_orientation[pair_orientation_FF] = pair_relation_invalid;
      paired_search->pair_orientation[pair_orientation_RR] = pair_relation_invalid;
      // Start parsing
      char *pair_orientation = strtok(optarg,",");
      while (pair_orientation!=NULL) {
        if (gem_strcaseeq(pair_orientation,"FR")) {
          paired_search->pair_orientation[pair_orientation_FR] = pair_relation_concordant; continue;
        }
        if (gem_strcaseeq(pair_orientation,"RF")) {
          paired_search->pair_orientation[pair_orientation_RF] = pair_relation_concordant; continue;
        }
        if (gem_strcaseeq(pair_orientation,"FF")) {
          paired_search->pair_orientation[pair_orientation_FF] = pair_relation_concordant; continue;
        }
        if (gem_strcaseeq(pair_orientation,"RR")) {
          paired_search->pair_orientation[pair_orientation_RR] = pair_relation_concordant; continue;
        }
        mapper_error_msg("Option '--pair-orientation' must be 'FR'|'RF'|'FF'|'RR'");
        pair_orientation = strtok(NULL,",");
      }
      return true;
    }
    case 402: { // --discordant-pair-orientation in {'FR'|'RF'|'FF'|'RR'} (default=RF)
      // Start parsing
      char *discordant_pair_orientation = strtok(optarg,",");
      while (discordant_pair_orientation!=NULL) {
        if (gem_strcaseeq(discordant_pair_orientation,"FR")) {
          paired_search->pair_orientation[pair_orientation_FR] = pair_relation_discordant; continue;
        }
        if (gem_strcaseeq(discordant_pair_orientation,"RF")) {
          paired_search->pair_orientation[pair_orientation_RF] = pair_relation_discordant; continue;
        }
        if (gem_strcaseeq(discordant_pair_orientation,"FF")) {
          paired_search->pair_orientation[pair_orientation_FF] = pair_relation_discordant; continue;
        }
        if (gem_strcaseeq(discordant_pair_orientation,"RR")) {
          paired_search->pair_orientation[pair_orientation_RR] = pair_relation_discordant; continue;
        }
        mapper_error_msg("Option '--discordant-pair-orientation' must be 'FR'|'RF'|'FF'|'RR'");
        discordant_pair_orientation = strtok(NULL,",");
      }
      return true;
    }
    case 403: { // --pair-layout in {'separate'|'overlap'|'contain'} (default=separated,overlap,contain)
      paired_search->pair_layout[pair_layout_separate] = pair_relation_invalid;
      paired_search->pair_layout[pair_layout_overlap] = pair_relation_invalid;
      paired_search->pair_layout[pair_layout_contain] = pair_relation_invalid;
      char *pair_layout = strtok(optarg,","); // Start parsing
      while (pair_layout!=NULL) {
        if (gem_strcaseeq(pair_layout,"separate")) {
          paired_search->pair_layout[pair_layout_separate] = pair_relation_concordant; continue;
        }
        if (gem_strcaseeq(pair_layout,"overlap"))  {
          paired_search->pair_layout[pair_layout_overlap] = pair_relation_concordant; continue;
        }
        if (gem_strcaseeq(pair_layout,"contain"))  {
          paired_search->pair_layout[pair_layout_contain] = pair_relation_concordant; continue;
        }
        mapper_error_msg("Option '--pair-layout' must be 'separate'|'overlap'|'contain'");
        pair_layout = strtok(NULL,",");
      }
      return true;
    }
    case 404: { // --discordant-pair-layout in {'separate'|'overlap'|'contain'} (default=none)
      char *pair_layout = strtok(optarg,","); // Start parsing
      while (pair_layout!=NULL) {
        if (gem_strcaseeq(pair_layout,"separate")) {
          paired_search->pair_layout[pair_layout_separate] = pair_relation_discordant; continue;
        }
        if (gem_strcaseeq(pair_layout,"overlap"))  {
          paired_search->pair_layout[pair_layout_overlap] = pair_relation_discordant; continue;
        }
        if (gem_strcaseeq(pair_layout,"contain"))  {
          paired_search->pair_layout[pair_layout_contain] = pair_relation_discordant; continue;
        }
        mapper_error_msg("Option '--pair-layout' must be 'separate'|'overlap'|'contain'");
        pair_layout = strtok(NULL,",");
      }
      return true;
    }
    case 405: // --pe-extension in {'none'|'recovery'|'shortcut'|'all'} (default=none)
      if (gem_strcaseeq(optarg,"none")) {
        paired_search->paired_end_extension_shortcut = false;
        paired_search->paired_end_extension_recovery = false;
      } else if (gem_strcaseeq(optarg,"shortcut")) {
        paired_search->paired_end_extension_shortcut = true;
        paired_search->paired_end_extension_recovery = false;
      } else if (gem_strcaseeq(optarg,"recovery")) {
        paired_search->paired_end_extension_shortcut = false;
        paired_search->paired_end_extension_recovery = true;
      } else if (gem_strcaseeq(optarg,"all")) {
        paired_search->paired_end_extension_shortcut = true;
        paired_search->paired_end_extension_recovery = true;
      } else {
        mapper_error_msg("Option '--pe-extension' must be 'none'|'recovery'|'shortcut'|'all'");
      }
      return true;
    case 406: { // --pe-template-length=<min>,<max>,<samples> (default=0,800,100)
      char *min=NULL, *max=NULL, *num_samples=NULL;
      const int num_arguments = input_text_parse_csv_arguments(optarg,3,&min,&max,&num_samples);
      mapper_cond_error_msg(num_arguments!=2,"Option '--pe-template-length' wrong number of arguments");
      // Parse minimum
      input_text_parse_extended_uint64(min,&paired_search->template_length_estimation_min);
      // Parse maximum
      input_text_parse_extended_uint64(max,&paired_search->template_length_estimation_max);
      // Parse num-samples
      input_text_parse_extended_uint64(num_samples,&paired_search->template_length_estimation_samples);
      return true;
    }
    default:
      return false; // Not found
  }
  // Not found
  return false;
}
bool gem_mapper_parse_arguments_bisulfite(
    mapper_parameters_t* const parameters,
    const int option,
    char* optarg) {
  // Parameters
  search_parameters_t* const search = &parameters->search_parameters;
  restriction_t *rest = NULL;
  // Bisulfite
  switch (option) {
    case 500: // --bisulfite_read
      if (gem_strcaseeq(optarg,"inferred-C2T-G2A")) {
        search->bisulfite_read = bisulfite_inferred_C2T_G2A;
        return true;
      }
      if (gem_strcaseeq(optarg,"inferred-G2A-C2T")) {
        search->bisulfite_read = bisulfite_inferred_G2A_C2T;
        return true;
      }
      if (gem_strcaseeq(optarg,"C2T")) {
        search->bisulfite_read = bisulfite_C2T;
        return true;
      }
      if (gem_strcaseeq(optarg,"G2A")) {
        search->bisulfite_read = bisulfite_G2A;
        return true;
      }
      if (gem_strcaseeq(optarg,"non-stranded")) {
        search->bisulfite_read = bisulfite_non_stranded;
        return true;
      }
      gem_fatal_error_msg("Option '--bisulfite-conversion' must be 'inferred_C2T_G2A'|'inferred_G2A_C2T'|'C2T'|'G2A'|'non-stranded'");
      return true;
    case 501: // --underconversion_sequence
      search->control_sequences[1] = strdup(optarg);
      return true;
    case 502: // --overconversion_sequence
      search->control_sequences[2] = strdup(optarg);
      return true;
    case 503: // --control_sequence
      search->control_sequences[0] = strdup(optarg);
      return true;
    case 504: // -- restriction-site
      rest = restriction_new(optarg);
      if(rest != NULL) {
        if(search->restriction_sites == NULL) {
          search->restriction_sites = vector_new(1, restriction_t *);
        }
        vector_insert(search->restriction_sites, rest, restriction_t *);
      } else {
        gem_fatal_error_msg("Error setting --restriction-site option");
      }
    return true;
    case 505: // --rrbs
      search->rrbs = true;
      return true;
    default:
      return false; // Not found
  }
  // Not found
  return false;
}
bool gem_mapper_parse_arguments_alignment_score(
    mapper_parameters_t* const parameters,
    const int option,
    char* optarg) {
  // Parameters
  search_parameters_t* const search = &parameters->search_parameters;
  // Alignment score
  switch (option) {
    case 600: // --alignment-model
      if (gem_strcaseeq(optarg,"none") || gem_strcaseeq(optarg,"pseudoalignment")) {
        search->match_alignment_model = match_alignment_model_none;
      } else if (gem_strcaseeq(optarg,"hamming")) {
        search->match_alignment_model = match_alignment_model_hamming;
      } else if (gem_strcaseeq(optarg,"edit") || gem_strcaseeq(optarg,"levenshtein") ) {
        search->match_alignment_model = match_alignment_model_levenshtein;
      } else if (gem_strcaseeq(optarg,"gap-affine")) {
        search->match_alignment_model = match_alignment_model_gap_affine;
      } else {
        mapper_error_msg("Option '--alignment-model' must be 'none'|'hamming'|'edit'|'gap-affine'");
      }
      return true;
    case 601: { // --gap-affine-penalties (A,B,O,X) (default=1,4,6,1)
      char *matching=NULL, *mismatch=NULL, *gap_open=NULL, *gap_extension=NULL;
      const int num_arguments = input_text_parse_csv_arguments(optarg,4,&matching,&mismatch,&gap_open,&gap_extension);
      mapper_cond_error_msg(num_arguments!=4,"Option '--gap-affine-penalties' wrong number of arguments");
      uint64_t matching_score, mismatch_penalty, gap_open_penalty, gap_extension_penalty;
      // Parse matching-score
      input_text_parse_extended_uint64(matching,&matching_score);
      // Parse mismatch-penalty
      input_text_parse_extended_uint64(mismatch,&mismatch_penalty);
      // Parse gap-open-penalty
      input_text_parse_extended_uint64(gap_open,&gap_open_penalty);
      // Parse gap-extension-penalty
      input_text_parse_extended_uint64(gap_extension,&gap_extension_penalty);
      // Configure scores
      search_configure_alignment_match_scores(search,matching_score);
      search_configure_alignment_mismatch_scores(search,mismatch_penalty);
      search->swg_penalties.gap_open_score = -((int32_t)gap_open_penalty);
      search->swg_penalties.gap_extension_score = -((int32_t)gap_extension_penalty);
      return true;
    }
    case 'A': { // --matching-score (default=1)
      uint64_t matching_score;
      input_text_parse_extended_uint64(optarg,&matching_score);
      search_configure_alignment_match_scores(search,matching_score);
      return true;
    }
    case 'B': { // --mismatch-penalty (default=4)
      uint64_t mismatch_penalty;
      input_text_parse_extended_uint64(optarg,&mismatch_penalty);
      search_configure_alignment_mismatch_scores(search,mismatch_penalty);
      return true;
    }
    case 'O': { // --gap-open-penalty (default=6)
      uint64_t gap_open_penalty;
      input_text_parse_extended_uint64(optarg,&gap_open_penalty);
      search->swg_penalties.gap_open_score = -((int32_t)gap_open_penalty);
      return true;
    }
    case 'X': { // --gap-extension-penalty (default=1)
      uint64_t gap_extension_penalty;
      input_text_parse_extended_uint64(optarg,&gap_extension_penalty);
      search->swg_penalties.gap_extension_score = -((int32_t)gap_extension_penalty);
      return true;
    }
    default:
      return false; // Not found
  }
  // Not found
  return false;
}
bool gem_mapper_parse_arguments_mapq_score(
    mapper_parameters_t* const parameters,
    const int option,
    char* optarg) {
  // Parameters
  search_parameters_t* const search = &parameters->search_parameters;
  // MAPQ score
  switch (option) {
    case 700: // --mapq-model in {'none'|'gem'|'classify'|'dump-predictors'} (default=gem)
      if (gem_strcaseeq(optarg,"none")) {
        search->mapq_model_se = mapq_model_none;
      } else if (gem_strcaseeq(optarg,"gem")) {
        search->mapq_model_se = mapq_model_gem;
      } else if (gem_strcaseeq(optarg,"classify")) {
        search->mapq_model_se = mapq_model_classify;
      } else if (gem_strcaseeq(optarg,"dump-predictors")) {
        search->mapq_model_se = mapq_model_dump_predictors;
      } else {
        mapper_error_msg("Option '--mapq-model' must be in {'none'|'gem'}");
      }
      return true;
    default:
      return false; // Not found
  }
  // Not found
  return false;
}
bool gem_mapper_parse_arguments_reporting(
    mapper_parameters_t* const parameters,
    const int option,
    char* optarg) {
  // Parameters
  search_parameters_t* const search = &parameters->search_parameters;
  // Reporting
  switch (option) {
    case 'm': // --min-reported-strata
      parameters->min_reported_strata_set = true;
      input_text_parse_extended_double(optarg,&search->select_parameters.min_reported_strata);
      return true;
    case 'M': // --max-reported-matches
      parameters->max_reported_matches_set = true;
      input_text_parse_extended_uint64(optarg,&search->select_parameters.max_reported_matches);
      return true;
    default:
      return false; // Not found
  }
  // Not found
  return false;
}
bool gem_mapper_parse_arguments_output_format(
    mapper_parameters_t* const parameters,
    const int option,
    char* optarg) {
  // Output format
  switch (option) {
    case 'F': // --output-format
      if (gem_strcaseeq(optarg,"MAP")) {
        parameters->io.output_format = MAP;
      } else if (gem_strcaseeq(optarg,"SAM")) {
        parameters->io.output_format = SAM;
      } else {
        mapper_error_msg("Option '-F|--output-format' must be 'MAP' or 'SAM'");
      }
      return true;
    case 900: // --sam-compact in {'true'|'false'} (default=true)
      parameters->io.sam_parameters.compact_xa = (optarg==NULL) ? true : input_text_parse_extended_bool(optarg);
      return true;
    case 'r': // --sam-read-group-header
      output_sam_parse_read_group_header(optarg,&parameters->io.sam_parameters);
      return true;
    case 901: // --sam-gem-compatible in {'true'|'false'} (default=true)
      parameters->io.sam_parameters.print_gem_fields = (optarg==NULL) ? true : input_text_parse_extended_bool(optarg);
      return true;
    case 902: { // --map-format in {'1'|'2'|'3'} (default=1)
      const uint64_t format_version = atol(optarg);
      if (format_version==1) {
        parameters->io.map_parameters.format_version = map_format_v1;
      } else if (format_version==2) {
        parameters->io.map_parameters.format_version = map_format_v2;
      } else if (format_version==3) {
        parameters->io.map_parameters.format_version = map_format_v3;
      } else {
        mapper_error_msg("Option '--map-format' must be in {'1'|'2'|'3'}");
      }
      return true;
    }
    default:
      return false; // Not found
  }
  // Not found
  return false;
}
bool gem_mapper_parse_arguments_system(
    mapper_parameters_t* const parameters,
    const int option,
    char* optarg) {
  // System
  switch (option) {
    case 't': // --threads
      mapper_cond_error_msg(parse_integer_system(optarg,&parameters->system.num_threads),
          "Option '--threads'. Error parsing 'num_threads'");
      return true;
    case 1000: // --max-memory
      mapper_cond_error_msg(
          input_text_parse_size(optarg,&(parameters->system.max_memory)),
          "Error parsing --max-memory. '%s' not a valid size (Eg. 2GB)",optarg);
      return true;
    case 1001: // --tmp-folder
      parameters->system.tmp_folder = optarg;
      return true;
    default:
      return false; // Not found
  }
  // Not found
  return false;
}
bool gem_mapper_parse_arguments_gpu(
    mapper_parameters_t* const parameters,
    const int option,
    char* optarg) {
  // GPU
  switch (option) {
    case 1100: // --gpu
      if (!gpu_supported()) GEM_CUDA_NOT_SUPPORTED();
      parameters->cuda.gpu_enabled = true;
      if (optarg && gem_strcaseeq(optarg,"emulated")) {
        parameters->cuda.cpu_emulation = true;
      }
      return true;
    case 1101: // --gpu-devices (default=all)
      if (gem_strcaseeq(optarg,"all")) {
        parameters->cuda.gpu_devices = UINT64_MAX;
      } else {
        char* devices_active[64];
        const int num_devices = input_text_parse_csv_array(optarg,devices_active,64);
        // Parse devices active
        uint64_t i;
        parameters->cuda.gpu_devices = 0;
        for (i=0;i<num_devices;++i) {
          int64_t device_no;
          if (!input_text_parse_integer((const char ** const)(devices_active + i),&device_no)) {
            parameters->cuda.gpu_devices |= ((uint64_t)1 << device_no);
          }
        }
      }
      break;
    case 1102: { // --gpu-buffers-model=2,3,3,3,3,1M
      if (!gpu_supported()) GEM_CUDA_NOT_SUPPORTED();
      char* num_fmi_bsearch_buffers=NULL;
      char* num_fmi_decode_buffers=NULL;
      char* num_bpm_kmer_filter_buffers=NULL;
      char* num_bpm_distance_buffers=NULL;
      char* num_bpm_align_buffers=NULL;
      char* buffer_size=NULL;
      const int num_arguments = input_text_parse_csv_arguments(optarg,6,
          &num_fmi_bsearch_buffers,&num_fmi_decode_buffers,
          &num_bpm_kmer_filter_buffers,&num_bpm_distance_buffers,
          &num_bpm_align_buffers,&buffer_size);
      mapper_cond_error_msg(num_arguments!=6,"Option '--gpu-buffers-model' wrong number of arguments");
      // Number of region-profile buffers per thread
      mapper_cond_error_msg(input_text_parse_integer(
          (const char** const)&num_fmi_bsearch_buffers,
          (int64_t*)&parameters->cuda.num_fmi_bsearch_buffers),
          "Option '--gpu-buffers-model'. Error parsing 'region-profile buffers'");
      // Number of decode-candidates buffers per thread
      mapper_cond_error_msg(input_text_parse_integer(
          (const char** const)&num_fmi_decode_buffers,
          (int64_t*)&parameters->cuda.num_fmi_decode_buffers),
          "Option '--gpu-buffers-model'. Error parsing 'decode-candidates buffers'");
      // Number of Kmer-Filter candidates buffers per thread
      mapper_cond_error_msg(input_text_parse_integer(
          (const char** const)&num_bpm_kmer_filter_buffers,
          (int64_t*)&parameters->cuda.num_kmer_filter_buffers),
          "Option '--gpu-buffers-model'. Error parsing 'Kmer-Filter candidates buffers'");
      // Number of BPM-Distance candidates buffers per thread
      mapper_cond_error_msg(input_text_parse_integer(
          (const char** const)&num_bpm_distance_buffers,
          (int64_t*)&parameters->cuda.num_bpm_distance_buffers),
          "Option '--gpu-buffers-model'. Error parsing 'BPM-Distance candidates buffers'");
      // Number of BPM-Align candidates buffers per thread
      mapper_cond_error_msg(input_text_parse_integer(
          (const char** const)&num_bpm_align_buffers,
          (int64_t*)&parameters->cuda.num_bpm_align_buffers),
          "Option '--gpu-buffers-model'. Error parsing 'BPM-Align candidates buffers'");
      // Buffer size
      mapper_cond_error_msg(input_text_parse_size(buffer_size,&parameters->cuda.gpu_buffer_size),
          "Option '--gpu-buffers-model'. Error parsing 'buffer_size'");
      return true;
    }
    default:
      return false; // Not found
  }
  // Not found
  return false;
}
bool gem_mapper_parse_arguments_misc(
    mapper_parameters_t* const parameters,
    const int option,
    char* optarg) {
  // Parameters
  search_parameters_t* const search = &parameters->search_parameters;
  // Debug/Misc
  switch (option) {
    case 'c': { // --check-alignments in {'correct'|'best'|'complete'}
      if (!optarg) {
        search->check_type = archive_check_correct;
      } else if (gem_strcaseeq(optarg,"none")) {
        search->check_type = archive_check_nothing;
      } else if (gem_strcaseeq(optarg,"correct")) {
        search->check_type = archive_check_correct;
      } else if (gem_strcaseeq(optarg,"first-best") || gem_strcaseeq(optarg,"best")) {
        search->check_type = archive_check_correct__first_optimum;
      } else if (gem_strcaseeq(optarg,"all-best")) {
        search->check_type = archive_check_correct__all_optimum;
      } else if (gem_strcaseeq(optarg,"complete")) {
        search->check_type = archive_check_correct__complete;
      } else {
        mapper_error_msg("Option '--check-alignments' must be 'correct'|'first-best'|'all-best'|'complete'");
      }
      return true;
    }
    case 1200: // --profile in {'sum'|'min'|'max'|'mean'|'sample'}
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
          mapper_error_msg("Option '--profile' must be 'sum'|'min'|'max'|'mean'|'sample'");
        }
      }
      return true;
    case 'v': // -v|--verbose in {'quiet'|'user'|'dev'}
      if (optarg) {
        if (gem_strcaseeq(optarg,"quiet")) {
          parameters->misc.verbose_user = false;
          parameters->misc.verbose_dev = false;
        } else if (gem_strcaseeq(optarg,"user")) {
          parameters->misc.verbose_user = true;
          parameters->misc.verbose_dev = false;
        } else if (gem_strcaseeq(optarg,"dev") || gem_strcaseeq(optarg,"debug")) {
          parameters->misc.verbose_user = true;
          parameters->misc.verbose_dev = true;
        } else {
          mapper_error_msg("Option '-v|--verbose' must be 'quiet'|'user'|'dev'");
        }
      } else {
        parameters->misc.verbose_user = true;
        parameters->misc.verbose_dev = false;
      }
      return true;
    case 1201: // --progress-step
      input_text_parse_extended_uint64(optarg,&parameters->io.mapper_ticker_step);
      return true;
    case 1202: // --version
      fprintf(stderr,"%s\n",parameters->gem_version);
      exit(0);
    case 'h':
      if (optarg==NULL || gem_strcaseeq(optarg,"user")) {
        gem_mapper_print_usage(VISIBILITY_USER);
      } else if (gem_strcaseeq(optarg,"advanced")) {
        gem_mapper_print_usage(VISIBILITY_ADVANCED);
      } else if (gem_strcaseeq(optarg,"developer")) {
        gem_mapper_print_usage(VISIBILITY_DEVELOPER);
      } else {
        mapper_error_msg("Help argument not valid {'user','advanced'}");
      }
      exit(0);
    default:
      return false; // Not found
  }
  // Not found
  return false;
}
void gem_mapper_parse_arguments(
    int argc,
    char** argv,
    mapper_parameters_t* const parameters,
    char* const gem_version) {
  // Check number of parameters (quick usage & exit)
  if (argc <= 1) {
    gem_mapper_print_usage(VISIBILITY_USER);
    exit(0);
  }
  // Set CMD line
  parameters->argc = argc;
  parameters->argv = argv;
  parameters->gem_version = gem_version;
  // Parse parameters
  struct option* getopt_options = options_adaptor_getopt(gem_mapper_options);
  string_t* const getopt_short_string = options_adaptor_getopt_short(gem_mapper_options);
  char* const getopt_short = string_get_buffer(getopt_short_string);
  int option, option_index;
  while (true) {
    // Get option &  Select case
    if ((option=getopt_long(argc,argv,getopt_short,getopt_options,&option_index))==-1) break;
    /* I/O */
    if (gem_mapper_parse_arguments_io(parameters,option,optarg)) continue;
    /* Single-end Alignment */
    if (gem_mapper_parse_arguments_single_end(parameters,option,optarg)) continue;
    /* Paired-end Alignment */
    if (gem_mapper_parse_arguments_paired_end(parameters,option,optarg)) continue;
    /* Bisulfite Alignment */
    if (gem_mapper_parse_arguments_bisulfite(parameters,option,optarg)) continue;
    /* Alignment Score */
    if (gem_mapper_parse_arguments_alignment_score(parameters,option,optarg)) continue;
    /* MAQ Score */
    if (gem_mapper_parse_arguments_mapq_score(parameters,option,optarg)) continue;
    /* Reporting */
    if (gem_mapper_parse_arguments_reporting(parameters,option,optarg)) continue;
    /* Output-format */
    if (gem_mapper_parse_arguments_output_format(parameters,option,optarg)) continue;
    /* System */
    if (gem_mapper_parse_arguments_system(parameters,option,optarg)) continue;
    /* GPU Settings */
    if (gem_mapper_parse_arguments_gpu(parameters,option,optarg)) continue;
    /* Debug/Misc */
    if(gem_mapper_parse_arguments_misc(parameters,option,optarg)) continue;
    /* Unrecognized option */
    mapper_error_msg("Option not recognized");
  }
  /*
   * Parameters Check
   */
  gem_mapper_parameters_check(parameters);
  // Free
  string_destroy(getopt_short_string);
  mm_free(getopt_short_string);
  mm_free(getopt_options);
}
