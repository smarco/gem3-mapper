/*
 * PROJECT: GEM-Tools library
 * FILE: gt.construct.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Constructor. Used to load examples/debug/testing/deploy/develop/...
 */

#include <getopt.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "gem_tools.h"

#define GT_EXAMPLE_MMAP_FILE false

typedef struct {
  char *name_input_file;
  char *name_output_file;
  uint64_t num_threads;
} gem_map_filter_args;

gem_map_filter_args parameters = {
    .name_input_file=NULL,
    .name_output_file=NULL,
    .num_threads=1
};

void gt_map_PE_2_SE() {
  // Open input file
  gt_input_file* input_file = (parameters.name_input_file==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file,GT_EXAMPLE_MMAP_FILE);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);

  gt_map_parser_attributes map_in_attributes = GT_MAP_PARSER_ATTR_DEFAULT(true);
  gt_output_map_attributes map_out_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();

  // Buffered reading of the file
  gt_status error_code;

  gt_buffered_input_file* buffered_input_file = gt_buffered_input_file_new(input_file);
  gt_buffered_output_file* buffered_output_file = gt_buffered_output_file_new(output_file);
  gt_buffered_input_file_attach_buffered_output(buffered_input_file,buffered_output_file);

#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(parameters.num_threads)
#endif
  {
    gt_template* template = gt_template_new();
    while ((error_code=gt_input_map_parser_get_template(buffered_input_file,template,&map_in_attributes))) {
      if (error_code==GT_INPUT_STATUS_FAIL) continue;

      const uint64_t num_blocks = gt_template_get_num_blocks(template);
      gt_cond_fatal_error_msg(num_blocks!=2,"Error processing. Not PE mapping found. Number of ends is %lu",num_blocks);
      if (gt_template_get_num_mmaps(template)>0) { // Select paired-end mappings
        gt_alignment* end1;
        gt_alignment_recalculate_counters(gt_template_get_block(template,0));
        end1 = gt_alignment_intersect_alignment_maps(gt_template_get_block(template,0),gt_template_get_block(template,0));
        gt_output_map_bofprint_alignment(buffered_output_file,end1,&map_out_attributes);
        gt_alignment_delete(end1);

        gt_alignment* end2;
        gt_alignment_recalculate_counters(gt_template_get_block(template,1));
        end2 = gt_alignment_intersect_alignment_maps(gt_template_get_block(template,1),gt_template_get_block(template,1));
        gt_output_map_bofprint_alignment(buffered_output_file,end2,&map_out_attributes);
        gt_alignment_delete(end2);
      }
    }
    gt_template_delete(template);
  }

  // Cleanup
  gt_buffered_input_file_close(buffered_input_file);
  gt_buffered_output_file_close(buffered_output_file);
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);

}

void usage() {
  fprintf(stderr, "USE: ./gt.pe2se [-i <input>] [-o <output>] [-t <threads>] [-h]\n"
                  "      Options::\n"
                  "        --input|i     <File>\n"
                  "        --output|o    <File>\n"
                  "        --threads|t    <Number>\n"
                  "        --help|h\n");
}
void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "threads", required_argument, 0, 't' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:o:t:h",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    case 't':
#ifdef HAVE_OPENMP
      parameters.num_threads = atol(optarg);
#endif
      break;
    case 'h':
      usage();
      exit(1);
    case '?': default:
      fprintf(stderr, "Option not recognized \n"); exit(1);
    }
  }
}

int main(int argc,char** argv) {
  // GT error handler
  gt_handle_error_signals();

  // Parsing command-line options
  parse_arguments(argc,argv);

  gt_map_PE_2_SE();

  return 0;
}


