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

#include <text/text_trace.h>
#include "text/dna_text.h"
#include "utils/essentials.h"
#include "utils/options_menu.h"
#include "utils/string_buffer.h"
#include "archive/archive.h"
#include "archive/archive_text.h"
#include "archive/locator.h"
#include "io/input_text.h"

/*
 * Version
 */
#define GEM_VERSION_STRING(version) QUOTE(version)
char* const gem_version = GEM_VERSION_STRING(GEM_VERSION);

/*
 * Debug
 */
//#define GEM_RETRIEVER_DEBUG_
/*
 * Error Handling
 */
#define gem_retriever_error_msg(error_msg,args...) \
  fprintf(stderr,"GEM-Retriever error:\n> "error_msg"\n",##args); \
  exit(1)
#define gem_retriever_cond_error_msg(condition,error_msg,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gem_retriever_error_msg(error_msg,##args); \
    } \
  } while (0)

/*
 * GEM-retriever parameters
 */
typedef struct {
  /* I/O */
  char* index_file_name;
  char* input_file_name;
  char* output_file_name;
  /* Verbose */
  bool verbose;
} retriever_parameters_t;
typedef struct {
  /* Retriever parameters */
  retriever_parameters_t parameters;
  /* I/O */
  FILE* input_stream;
  FILE* output_stream;
  /* Archive */
  archive_t* archive;
  /* Misc */
  mm_allocator_t* mm_allocator;
} retriever_data_t;
// Defaults
void retriever_parameters_set_defaults(retriever_parameters_t* const parameters) {
  /* I/O */
  parameters->index_file_name=NULL;
  parameters->input_file_name=NULL;
  parameters->output_file_name=NULL;
  /* Verbose */
  parameters->verbose = true;
}
/*
 * GEM-Retriever options Menu
 */
option_t gem_retriever_options[] = {
  /* I/O */
  { 'I', "index", REQUIRED, TYPE_STRING, 2, VISIBILITY_USER, "<index_file.gem>", "(GEM archive)" },
  { 'i', "input", REQUIRED, TYPE_STRING, 2 , VISIBILITY_USER, "<input_file>" , "(default=stdin)" },
  { 'o', "output", REQUIRED, TYPE_STRING, 2 , VISIBILITY_USER, "<output_prefix>" , "(default=stdout)" },
  /* Miscellaneous */
  { 'h',  "help", OPTIONAL, TYPE_NONE, 3, VISIBILITY_USER, "" , "(print usage)" },
  { 300, "version", NO_ARGUMENT, TYPE_STRING, 3, VISIBILITY_USER, "" , "" },
  {  0, "", 0, 0, 0, false, "", ""}
};
char* gem_retriever_groups[] = {
  /* 0 */ "Null",
  /* 1 */ "Unclassified",
  /* 2 */ "I/O",
  /* 3 */ "Miscellaneous",
};
void usage(const option_visibility_t visibility_level) {
  fprintf(stderr, "USAGE: ./gem-retriever [ARGS]...\n");
  options_fprint_menu(stderr,gem_retriever_options,gem_retriever_groups,true,visibility_level);
}
void parse_arguments(int argc,char** argv,retriever_parameters_t* const parameters) {
  struct option* getopt_options = options_adaptor_getopt(gem_retriever_options);
  string_t* const getopt_short_string = options_adaptor_getopt_short(gem_retriever_options);
  char* const getopt_short = string_get_buffer(getopt_short_string);
  int option, option_index;
  while (true) {
    // Get option &  Select case
    if ((option=getopt_long(argc,argv,getopt_short,getopt_options,&option_index))==-1) break;
    switch (option) {
    /* I/O */
    case 'I': // --index
      parameters->index_file_name = optarg;
      break;
    case 'i': // --input
      parameters->input_file_name = optarg;
      break;
    case 'o': // --output
      parameters->output_file_name = optarg;
      break;
    case 'h':
      if (optarg==NULL || gem_strcaseeq(optarg,"user")) {
        usage(VISIBILITY_USER);
      } else if (gem_strcaseeq(optarg,"advanced")) {
        usage(VISIBILITY_ADVANCED);
      } else if (gem_strcaseeq(optarg,"developer")) {
        usage(VISIBILITY_DEVELOPER);
      } else {
        gem_retriever_error_msg("Help argument not valid {'user','advanced'}");
      }
      exit(0);
    case 300: // --version
      fprintf(stderr,"%s\n",gem_version);
      exit(0);
    /* */
    case '?':
    default:
      gem_retriever_error_msg("Option not recognized");
    }
  }
  /*
   * Parameters Checks
   */
  // I/O Parameters
  gem_retriever_cond_error_msg(parameters->index_file_name==NULL,"Index file required");
  parameters->verbose = (parameters->input_file_name==NULL) && isatty(0);
  /*
   * Free
   */
  string_destroy(getopt_short_string);
  mm_free(getopt_short_string);
  mm_free(getopt_options);
}
/*
 * Query
 */
typedef struct {
  uint64_t id;
  uint8_t* tag;
  strand_t strand;
  bs_strand_t bs_strand;
  uint64_t text_position;
  uint64_t text_length;
} retriever_query_t;
/*
 * Parsing
 */
int retriever_query_parse(
    retriever_query_t* const retriever_query,
    char* text_line,
    const uint64_t text_line_length) {
  // Parse sequence name
  if (text_line_length==0) return -1;
  retriever_query->tag = (uint8_t*)text_line;
  while (*text_line != '\0') {
    if (*text_line==':' || *text_line=='\t') break;
    ++text_line;
  }
  if (*text_line!=':' && *text_line!='\t') return -1;
  *text_line = '\0';
  ++text_line;
  // Parse strand
  switch (*text_line) {
    case '+': retriever_query->strand = Forward; break;
    case 'F': retriever_query->strand = Forward; break;
    case '-': retriever_query->strand = Reverse; break;
    case 'R': retriever_query->strand = Reverse; break;
    default: return -1;
  }
  ++text_line;
  // Parse bs-strand (optional)
  if (strcmp(text_line,"C2T")==0) {
    retriever_query->bs_strand = bs_strand_C2T;
    text_line += 3;
  } else if (strcmp(text_line,"G2A")==0) {
    retriever_query->bs_strand = bs_strand_G2A;
    text_line += 3;
  } else {
    retriever_query->bs_strand = bs_strand_none;
  }
  if (*text_line!=':' && *text_line!='\t') return -1;
  ++text_line;
  // Parse position
  int64_t value;
  input_text_parse_integer((const char** const)&text_line,&value);
  retriever_query->text_position = value;
  if (*text_line!=':' && *text_line!='\t') return -1;
  ++text_line;
  // Parse length
  input_text_parse_integer((const char** const)&text_line,&value);
  retriever_query->text_length = value;
  return 0;
}
void retriever_query_location(
    retriever_data_t* const retriever_data,retriever_query_t* const retriever_query) {
  // Locate the sequence
  locator_interval_t* const locator_interval =
      locator_inverse_map(retriever_data->archive->locator,retriever_query->tag,
          Forward,retriever_query->bs_strand,retriever_query->text_position);
  if (locator_interval==NULL) {
    gem_retriever_error_msg("[#%"PRIu64"]Sequence location (%s:%c%s:%"PRIu64") not found in the archive\n",
        retriever_query->id,retriever_query->tag,
        (retriever_query->strand==Forward) ? '+' : '-',
        (retriever_query->bs_strand==bs_strand_C2T) ? "C2T" :
            ((retriever_query->bs_strand==bs_strand_G2A) ? "G2A" :
                ((retriever_query->bs_strand==bs_strand_none) ? "" : "ERROR")),
        retriever_query->text_position);
    fprintf(retriever_data->output_stream,"-\n");
    return;
  }
  if (locator_interval->type == locator_interval_uncalled) {
    fprintf(retriever_data->output_stream,"N\n");
    return;
  }
  // Adjust the sequence boundaries
  const uint64_t index_begin_position = locator_interval->begin_position +
      retriever_query->text_position-locator_interval->sequence_offset;
  uint64_t index_end_position = index_begin_position + retriever_query->text_length;
  if (index_end_position > locator_interval->end_position) {
    index_end_position = locator_interval->end_position;
  }
  const uint64_t text_length = index_end_position-index_begin_position;
  // Retrieve the sequence
  text_trace_t text_trace;
  archive_text_retrieve(
      retriever_data->archive->text,index_begin_position,
      text_length,retriever_query->strand==Reverse,false,
      &text_trace,retriever_data->mm_allocator);
  const uint8_t* const text = text_trace.text; // Candidate
  // Output the sequence
  uint64_t i;
  for (i=0;i<text_length;++i) {
    fprintf(retriever_data->output_stream,"%c",dna_decode(text[i]));
  }
  fprintf(retriever_data->output_stream,"\n");
}
/*
 * Main()
 */
int main(int argc,char** argv) {
  // Parsing command-line options
  retriever_data_t retriever_data;
  retriever_parameters_t* const parameters = &retriever_data.parameters;
  retriever_parameters_set_defaults(parameters);
  parse_arguments(argc,argv,parameters);

  // Open I/O files
  retriever_data.input_stream = (parameters->input_file_name==NULL) ? stdin : fopen(parameters->input_file_name,"r");
  retriever_data.output_stream = (parameters->output_file_name==NULL) ? stdout : fopen(parameters->output_file_name,"w");

  // Load archive
  gem_cond_log(parameters->verbose,"[Loading GEM index '%s']",parameters->index_file_name);
  retriever_data.archive = archive_read(parameters->index_file_name,false);
  gem_cond_log(parameters->verbose,"... done");

  // Allocate
  mm_slab_t* const mm_slab = mm_slab_new_(BUFFER_SIZE_64M,BUFFER_SIZE_512M,MM_UNLIMITED_MEM);
  retriever_data.mm_allocator = mm_allocator_new(mm_slab);

  // Read all retriever queries
  retriever_query_t retriever_query;
  char* input_buffer = NULL;
  size_t input_buffer_size = 0;
  retriever_query.id = 0;
  while (getline(&input_buffer,&input_buffer_size,retriever_data.input_stream)>=0) {
    // Next query
    ++(retriever_query.id);
    // Parse query
    if (retriever_query_parse(&retriever_query,input_buffer,input_buffer_size)!=-1) {
      // Query & output
      retriever_query_location(&retriever_data,&retriever_query);
    }
    mm_allocator_clear(retriever_data.mm_allocator);
  }

  // Clean-up
  if (input_buffer!=NULL) free(input_buffer);
  mm_allocator_delete(retriever_data.mm_allocator);
  mm_slab_delete(mm_slab);
  archive_delete(retriever_data.archive); // Delete archive
  return 0;
}

