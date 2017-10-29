/*
 * PROJECT: GEM-Tools library
 * FILE: gt.mapset.c
 * DATE: 08/11/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Utility to perform set operations {UNION,INTERSECTION,DIFFERENCE} over alignment files {MAP,SAM}
 */

#include <getopt.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "gem_tools.h"

typedef enum { GT_MAP_SET_UNKNOWN,
               GT_MAP_SET_INTERSECTION, GT_MAP_SET_UNION, GT_MAP_SET_DIFFERENCE,
               GT_MAP_SET_JOIN, GT_MAP_SET_COMPARE,
               GT_MERGE_MAP, GT_DISPLAY_COMPACT_MAP} gt_operation;

typedef struct {
  gt_operation operation;
  char* name_input_file_1;
  char* name_input_file_2;
  char* name_output_file;
  bool mmap_input;
  bool paired_end;
  bool files_contain_same_reads;
  double eq_threshold;
  bool strict;
  bool verbose;
  uint64_t num_threads;
} gt_stats_args;

gt_stats_args parameters = {
    .operation=GT_MAP_SET_UNKNOWN,
    .name_input_file_1=NULL,
    .name_input_file_2=NULL,
    .name_output_file=NULL,
    .mmap_input=false,
    .paired_end=false,
    .files_contain_same_reads=false,
    .eq_threshold=0.5,
    .strict=false,
    .verbose=false,
    .num_threads=1
};
/*
 * CMP
 */
uint64_t current_read_length;
int64_t gt_mapset_map_cmp(gt_map* const map_1,gt_map* const map_2) {
  const uint64_t eq_threshold = (parameters.eq_threshold <= 1.0) ?
      parameters.eq_threshold*current_read_length: parameters.eq_threshold;
  return parameters.strict ? gt_map_cmp(map_1,map_2) : gt_map_range_cmp(map_1,map_2,eq_threshold);
}
int64_t gt_mapset_mmap_cmp(gt_map** const map_1,gt_map** const map_2,const uint64_t num_maps) {
  const uint64_t eq_threshold = (parameters.eq_threshold <= 1.0) ?
      parameters.eq_threshold*current_read_length: parameters.eq_threshold;
  return parameters.strict ? gt_mmap_cmp(map_1,map_2,num_maps) : gt_mmap_range_cmp(map_1,map_2,num_maps,eq_threshold);
}
/*
 * Synch
 */
GT_INLINE gt_status gt_mapset_read_template_get_commom_map(
    gt_buffered_input_file* const buffered_input_master,gt_buffered_input_file* const buffered_input_slave,
    gt_template* const template_master,gt_template* const template_slave) {
  gt_status error_code_master, error_code_slave;
  gt_generic_parser_attributes* generic_parser_attr = gt_input_generic_parser_attributes_new(parameters.paired_end);
  // Read master
  if ((error_code_master=gt_input_generic_parser_get_template(
      buffered_input_master,template_master,generic_parser_attr))==GT_INPUT_STATUS_FAIL) {
    gt_fatal_error_msg("Fatal error parsing file <<Master>>");
  }
  if (error_code_master==GT_INPUT_STATUS_EOF) return GT_INPUT_STATUS_EOF;
  // Read slave
  if ((error_code_slave=gt_input_generic_parser_get_template(
      buffered_input_slave,template_slave,generic_parser_attr))==GT_INPUT_STATUS_FAIL) {
    gt_fatal_error_msg("Fatal error parsing file <<Slave>>");
  }
  if (error_code_slave==GT_INPUT_STATUS_EOF) { // Check EOF conditions
    gt_fatal_error_msg("<<Slave>> is not contained in master <<Master>> (looking for '"PRIgts"')",
        PRIgts_content(gt_template_get_string_tag(template_master)));
  }
  // Synch loop
  while (gt_string_cmp(gt_template_get_string_tag(template_master),gt_template_get_string_tag(template_slave))) {
    // Fetch next slave's template
    if ((error_code_master=gt_input_generic_parser_get_template(
        buffered_input_slave,template_slave,generic_parser_attr))!=GT_INPUT_STATUS_OK) {
      gt_fatal_error_msg("<<Slave>> is not contained in master <<Master>> (looking for '"PRIgts"')",
          PRIgts_content(gt_template_get_string_tag(template_master)));
    }
  }
  return GT_INPUT_STATUS_OK;
}
/*
 * Core
 */
void gt_get_subdominant() {
  // File IN/OUT
  gt_input_file* input_file_1 = gt_input_file_open(parameters.name_input_file_1,parameters.mmap_input);
  gt_input_file* input_file_2 = (parameters.name_input_file_2==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file_2,parameters.mmap_input);
  if (parameters.name_input_file_2==NULL) GT_SWAP(input_file_1,input_file_2);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);

  // Buffered I/O
  gt_buffered_input_file* buffered_input_1 = gt_buffered_input_file_new(input_file_1);
  gt_buffered_input_file* buffered_input_2 = gt_buffered_input_file_new(input_file_2);
  gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
  gt_buffered_input_file_attach_buffered_output(buffered_input_1,buffered_output);

  // Template I/O (synch)
  gt_template *template_1 = gt_template_new();
  gt_template *template_2 = gt_template_new();
  gt_output_map_attributes* output_map_attributes = gt_output_map_attributes_new();
  while (gt_mapset_read_template_get_commom_map(buffered_input_1,buffered_input_2,template_1,template_2)) {
    // Record current read length
    current_read_length = gt_template_get_total_length(template_1);

    GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template_1,alignment_1) {
      GT_TEMPLATE_REDUCTION(template_2,alignment_2);
      // TYPE A
      if (gt_alignment_get_num_maps(alignment_1) > 1) { // Multimap
        // Do the intersection
        gt_alignment* const alignment_intersection = gt_alignment_intersect_alignment_maps_fx(gt_mapset_map_cmp,alignment_1,alignment_2);
        if (gt_alignment_get_num_maps(alignment_intersection) >= 1) { // Intersection non-empty
          gt_map* const mmap_first = gt_alignment_get_map(alignment_1,0);
          gt_map* const mmap_intersection = gt_alignment_get_map(alignment_intersection,0);

          const uint64_t distance_mmap_first = gt_map_get_distance(mmap_first);
          const uint64_t distance_mmap_intersection = gt_map_get_distance(mmap_intersection);
          if (distance_mmap_first < distance_mmap_intersection) { // Sub-dominant
            //gt_output_map_bofprint_alignment(buffered_output,alignment_1,output_map_attributes);
            gt_output_map_bofprint_tag(buffered_output,template_1->tag,template_1->attributes,output_map_attributes);
            gt_bofprintf(buffered_output,"\t");
            gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_1),
                template_1->attributes,output_map_attributes); // Master's Counters
            gt_bofprintf(buffered_output,"\t");
            gt_output_map_bofprint_alignment_maps(buffered_output,alignment_intersection,output_map_attributes); // Slave's Maps
            gt_bofprintf(buffered_output,"\t");
            gt_output_map_bofprint_template_maps(buffered_output,template_1,output_map_attributes); // Master's Maps
            gt_bofprintf(buffered_output,"\n");
          }
        }
        gt_alignment_delete(alignment_intersection);
      }
    } GT_TEMPLATE_END_REDUCTION;


  }
  // Clean
  gt_template_delete(template_1);
  gt_template_delete(template_2);
  gt_buffered_input_file_close(buffered_input_1);
  gt_buffered_input_file_close(buffered_input_2);
  gt_buffered_output_file_close(buffered_output);
  gt_input_file_close(input_file_1);
  gt_input_file_close(input_file_2);
  gt_output_file_close(output_file);
//    // Apply operation
//    switch (parameters.operation) {
//      case GT_MAP_SET_JOIN:
//        // Print Master's TAG+Counters+Maps
//        gt_output_map_bofprint_tag(buffered_output,template_1->tag,template_1->attributes,output_map_attributes);
//        gt_bofprintf(buffered_output,"\t");
//        gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_1),
//            template_1->attributes,output_map_attributes); // Master's Counters
//        gt_bofprintf(buffered_output,"\t");
//        gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_2),
//            template_1->attributes,output_map_attributes); // Slave's Counters
//        gt_bofprintf(buffered_output,"\t");
//        gt_output_map_bofprint_template_maps(buffered_output,template_1,output_map_attributes); // Master's Maps
//        gt_bofprintf(buffered_output,"\t");
//        gt_output_map_bofprint_template_maps(buffered_output,template_2,output_map_attributes); // Slave's Maps
//        gt_bofprintf(buffered_output,"\n");
//        break;
//      case GT_MAP_SET_COMPARE: {
//        // Perform simple cmp operations
//        gt_template *template_master_minus_slave=gt_template_subtract_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_1,template_2);
//        gt_template *template_slave_minus_master=gt_template_subtract_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_2,template_1);
//        gt_template *template_intersection=gt_template_intersect_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_1,template_2);
//        /*
//         * Print results :: (TAG (Master-Slave){COUNTER MAPS} (Slave-Master){COUNTER MAPS} (Intersection){COUNTER MAPS})
//         */
//        gt_output_map_bofprint_tag(buffered_output,template_1->tag,template_1->attributes,output_map_attributes);
//        // Counters
//        gt_bofprintf(buffered_output,"\t");
//        gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_master_minus_slave),
//            template_master_minus_slave->attributes,output_map_attributes); // (Master-Slave){COUNTER}
//        gt_bofprintf(buffered_output,"\t");
//        gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_slave_minus_master),
//            template_slave_minus_master->attributes,output_map_attributes); // (Slave-Master){COUNTER}
//        gt_bofprintf(buffered_output,"\t");
//        gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_intersection),
//            template_intersection->attributes,output_map_attributes); // (Intersection){COUNTER}
//        // Maps
//        gt_bofprintf(buffered_output,"\t");
//        gt_output_map_bofprint_template_maps(buffered_output,template_master_minus_slave,output_map_attributes); // (Master-Slave){COUNTER}
//        gt_bofprintf(buffered_output,"\t");
//        gt_output_map_bofprint_template_maps(buffered_output,template_slave_minus_master,output_map_attributes); // (Slave-Master){COUNTER}
//        gt_bofprintf(buffered_output,"\t");
//        gt_output_map_bofprint_template_maps(buffered_output,template_intersection,output_map_attributes); // (Intersection){COUNTER}
//        gt_bofprintf(buffered_output,"\n");
//        // Delete templates
//        gt_template_delete(template_master_minus_slave);
//        gt_template_delete(template_slave_minus_master);
//        gt_template_delete(template_intersection);
//        }
//        break;
//      default:
//        gt_fatal_error(SELECTION_NOT_VALID);
//        break;
//    }
//  }
}
/*
 * Menu
 */
gt_option gt_constructor_options[] = {
  /* I/O */
  { 300, "i1", GT_OPT_REQUIRED, GT_OPT_STRING, 3 , true, "<file>" , "" },
  { 301, "i2", GT_OPT_REQUIRED, GT_OPT_STRING, 3 , true, "<file>" , "" },
  { 'p', "paired-end", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , true, "" , "" },
  { 302, "mmap-input", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3 , false, "" , "" },
  { 'o', "output", GT_OPT_REQUIRED, GT_OPT_STRING, 3 , true, "<file>" , "" },
  /* Compare Function */
  { 400, "eq-th", GT_OPT_REQUIRED, GT_OPT_FLOAT, 4 , true, "<integer>|<float> (Difference tolerated between positions)" , "" },
  { 401, "strict", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4 , true, "(Strict comparison of mappings)" , "" },
  /* Misc */
  { 'v', "verbose", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 5, true, "", ""},
#ifdef HAVE_OPENMP
  { 't', "threads", GT_OPT_REQUIRED, GT_OPT_INT, 5, true, "", ""},
#endif
  {  0, "", 0, 0, 0, false, "", ""}
};
char* gt_constructor_groups[] = {
  /*  0 */ "Null",
  /*  1 */ "Unclassified",
  /*  2 */ "Operations",
  /*  3 */ "I/O",
  /*  4 */ "Compare Function",
  /*  5 */ "Misc"
};
void parse_arguments(int argc,char** argv) {
  struct option* gt_constructor_getopt = gt_options_adaptor_getopt(gt_constructor_options);
  gt_string* const gt_constructor_short_getopt = gt_options_adaptor_getopt_short(gt_constructor_options);
  int option, option_index;
  while (true) {
    // Get option & Select case
    if ((option=getopt_long(argc,argv,
        gt_string_get_string(gt_constructor_short_getopt),gt_constructor_getopt,&option_index))==-1) break;
    switch (option) {
    /* I/O */
    case 300:
      parameters.name_input_file_1 = optarg;
      break;
    case 301:
      parameters.name_input_file_2 = optarg;
      break;
    case 'p':
      parameters.paired_end = true;
      break;
    case 302:
      parameters.mmap_input = true;
      gt_fatal_error(NOT_IMPLEMENTED);
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    /* Compare Function */
    case 400: // eq-th
      parameters.eq_threshold = atof(optarg);
      break;
    case 401: // strict
      parameters.strict = true;
      break;
    /* Misc */
    case 'v':
      parameters.verbose = true;
      break;
    case 't':
#ifdef HAVE_OPENMP
      parameters.num_threads = atol(optarg);
#endif
      break;
    case 'h':
      fprintf(stderr, "USE: ./gt.constructor [OPERATION] [ARGS]...\n");
      gt_options_fprint_menu(stderr,gt_constructor_options,gt_constructor_groups,false,false);
      exit(1);
    case '?':
    default:
      gt_fatal_error_msg("Option not recognized");
    }
  }
  // Free
  gt_string_delete(gt_constructor_short_getopt);
}

int main(int argc,char** argv) {
  // GT error handler
  gt_handle_error_signals();

  // Parsing command-line options
  parse_arguments(argc,argv);

  // Do it !
  gt_get_subdominant();

  return 0;
}


