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
 * PROJECT: GEM-Mappability v3 (GEM3)
 * AUTHOR(S): Paolo Ribeca <paolo.ribeca@gmail.com>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "utils/essentials.h"
#include "archive/archive.h"
#include "archive/archive_text.h"
#include "archive/locator.h"
#include "archive/search/archive_search.h"
#include "archive/search/archive_search_handlers.h"
#include "archive/search/archive_search_se.h"
#include "archive/search/archive_search_se_parameters.h"
#include "io/output_map.h"
#include "matches/match_trace.h"
#include "matches/matches.h"
#include "matches/matches_counters.h"
#include <getopt.h>

/*
 * Version
 */
#define GEM_VERSION_STRING(version) QUOTE(version)
char* const gem_version = GEM_VERSION_STRING(GEM_VERSION);

/*
 * Errors
 */
#define GEM_ERROR_MAPP_LOG "Integer %lu is out of logarithm range"
#define GEM_ERROR_MAPP_READ_LEN "The specified read length exceeds the index size"
#define GEM_ERROR_FILE_WRITE "Couldn't write output"

/*
 * Debug
 */
#define DEBUG_MAPPABILITY_C false

/*
 * Constants
 */
#define VLOGP32_0 ' ' // vlogp32(0)
#define VLOGP32_1 '!' // vlogp32(1)

/*
 * Mappability Parameters
 */
typedef struct {
  /* CMD line */
  int argc;
  char** argv;
  char* gem_version;
  /* I/O */
  char* index_file_name;
  char* output_file_name;
  uint64_t output_line_width;
  /* Mappability structures */
  archive_t* archive;
  uint64_t remaining;               // Current position
  pthread_mutex_t remaining_mutex;
  vector_t* bitmap_buf;             // The logarithm bitmap
  uint64_t done_new;
  uint64_t done_unique;
  pthread_mutex_t bitmap_mutex;
  uint64_t log_table[95];           // Logarithm table (128-33[\000-\032]-2[\126\127])[intervals]+2
  uint64_t eff_text_len;            // Optimization
  /* Search parameters */
  search_parameters_t search_parameters;
  uint64_t read_length;
  uint64_t approx_threshold;
  uint64_t max_errors;
  /* Misc */
  bool verbose;
  bool debug;
  uint64_t granularity;
  uint64_t num_threads;
} mappability_parameters_t;
mappability_parameters_t parameters = {
    /* I/O */
    .output_line_width = 70,
    /* Search */
    .approx_threshold = 0,
    .max_errors = 0,
    /* Misc */
    .verbose = true,
    .debug = false,
    .granularity = 10000,
    .num_threads = 1,
};
/*
 * Bitmap Data Structures
 */
typedef struct {
  uint64_t position;
  uint8_t vlog;
} bitmap_t;
typedef struct {
  uint64_t position;
  UT_hash_handle hh;
} hashed_t;
/*
 * Index loader
 */
archive_t* mappability_load_index(
    char* const index_file_name,
    const bool verbose) {
  gem_cond_log(verbose,"[Loading GEM index '%s']",index_file_name);
  return archive_read(index_file_name,false);
}
/*
 * Output
 */
FILE* mappability_open_output(
    char* const output_file_name,
    const uint64_t read_length,
    const uint64_t effective_max_errors,
    const uint64_t approx_threshold,
    uint64_t* const log_table) {
  // Compose output file name
  char* const output_file_name__ext = calloc(strlen(output_file_name)+14,sizeof(char));
  gem_cond_fatal_error_msg(!output_file_name__ext,"Couldn't allocate memory");
  sprintf(output_file_name__ext,"%s.mappability",output_file_name);
  FILE* output = fopen(output_file_name__ext,"w");
  gem_cond_fatal_error_msg(!output,"Couldn't open output file");
  // Write header
  fprintf(output,"~~K-MER LENGTH\n%lu\n~~APPROXIMATION THRESHOLD\n",(uint64_t)read_length);
  if (approx_threshold == UINT64_MAX) {
    fprintf(output,"disabled\n");
  } else {
    fprintf(output,"%lu\n",(uint64_t)approx_threshold);
  }
  fprintf(output,
      "~~MAX ERRORS\n%lu\n"
      "~~ENCODING\n",
      effective_max_errors);
  // Dump logarithmic table
  uint64_t i;
  for (i=0;i<94;++i) {
    fprintf(output,"'%c'~[%lu-%lu]\n",(uint8_t)(32+i),log_table[i],log_table[i+1]-1);
  }
  fflush(output);
  // Return output file handler
  return output;
}
/*
 * Logaritmic table computation
 */
void mappability_compute_log_table(
    const uint64_t read_length,
    uint64_t* const log_table) {
  // Compute lengths
  const uint64_t text_len = archive_get_index_length(parameters.archive);
  gem_cond_fatal_error(text_len+1<read_length,MAPP_READ_LEN);
  parameters.eff_text_len = text_len-read_length+1;
  /*
   * Setup of the logarithm table.
   * The maximum has to be doubled, since both strands may contribute.
   */
  uint64_t table_max = 2*parameters.eff_text_len; // Account for complement
  double log_table_max = log(table_max);
  log_table[0] = 0;
  log_table[1] = 1;
  double multiplier = exp(ceil(log_table_max*100)/100./93.);
  double prod = 1.;
  int64_t i;
  for (i=2;i<=94;++i) {
    prod *= multiplier;
    if ((uint64_t)prod == log_table[i-1]) {
      prod = ceil(prod);
      multiplier = exp(ceil((log_table_max-log(prod))*100)/100./(93.-i));
    }
    log_table[i] = prod;
    if (!parameters.approx_threshold && log_table[i]>log_table[i-1]+1) {
      parameters.approx_threshold = log_table[i-1];
    }
  }
  assert(log_table[94] > table_max);
}
void mappability_retrieve_sequence(
    archive_text_t* const archive_text,
    const uint64_t idx_position,
    const uint64_t read_length,
    sequence_t* const sequence,
    uint64_t* const num_invalid_chars,
    mm_allocator_t* const mm_allocator) {
  // Decode the needed reference text
  text_trace_t text_trace;
  archive_text_retrieve(archive_text,idx_position,read_length,false,false,&text_trace,mm_allocator);
  // Prepare sequence
  sequence->tag.buffer = "Mappability";
  sequence->tag.length = strlen(sequence->tag.buffer);
  sequence->read.buffer = mm_allocator_calloc(mm_allocator,read_length+1,char,false);
  uint64_t i, invalid_chars = 0;
  for (i=0;i<read_length;++i) {
    const char base = dna_decode(text_trace.text[i]);
    sequence->read.buffer[i] = base;
    if (!is_dna_canonical(base)) ++invalid_chars; // No matches containing Ns can ever be found
  }
  sequence->read.buffer[read_length] = '\0';
  sequence->read.length = read_length;
  sequence->qualities.buffer = "";
  sequence->qualities.length = 0;
  sequence->has_qualities = false;
  sequence->end_info = single_end;
  *num_invalid_chars = invalid_chars;
}
/*
 * Mappability mapping sequence
 */
void mappability_print_countersaccount_mcs(
    FILE* const stream,
    const uint64_t current_counter_pos,
    const uint64_t num_zeros) {
  uint64_t i = 0;
  if (current_counter_pos==0) {
    fprintf(stream,"0");
    i=1;
  }
  for (;i<num_zeros;++i) {
    fprintf(stream,":0");
  }
  fprintf(stream,"+0");
}
void mappability_print_counters(
    FILE* const stream,
    matches_counters_t* const matches_counter,
    const uint64_t mcs) {
  const uint64_t num_counters = matches_counters_get_num_counters(matches_counter);
  // Zero counters
  if (gem_expect_false(num_counters==0)) {
    if (mcs==0 || mcs==ALL) {
      fprintf(stream,"0");
    } else {
      mappability_print_countersaccount_mcs(stream,0,mcs);
    }
    return;
  }
  // Print counters
  uint64_t i = 0;
  const uint64_t* counters = matches_counters_get_counts(matches_counter);
  while (i < num_counters) {
    // Print Counter
    if (i>0) fprintf(stream,"%c",(mcs==i?'+':':'));
    fprintf(stream,"%lu",*counters);
    i++; // Next (+1)
    ++counters;
  }
  // Account for MCS
  if (i<=mcs) {
    mappability_print_countersaccount_mcs(stream,i,mcs-i);
  }
}
void mappability_print_cigar(
    FILE* const stream,
    const cigar_element_t* cigar_array,
    const uint64_t cigar_length) {
  // Traverse all CIGAR elements
  uint64_t i;
  for (i=0;i<cigar_length;++i,++cigar_array) {
    // Print CIGAR element
    switch (cigar_array->type) {
      case cigar_match:
        fprintf(stream,"%d",(uint32_t)cigar_array->length);
        break;
      case cigar_mismatch:
        fprintf(stream,"%c",dna_decode(cigar_array->mismatch));
        break;
      case cigar_ins:
        fprintf(stream,"%c",'>');
        fprintf(stream,"%d",cigar_array->length);
        fprintf(stream,"%c",'+');
        break;
      case cigar_del:
        if (cigar_array->attributes == cigar_attr_trim) {
          fprintf(stream,"%c",'(');
          fprintf(stream,"%d",cigar_array->length);
          fprintf(stream,"%c",')');
        } else {
          fprintf(stream,"%c",'>');
          fprintf(stream,"%d",cigar_array->length);
          fprintf(stream,"%c",'-');
        }
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
void mappability_mapping_sequence(
    archive_search_handlers_t* const archive_search_handlers,
    archive_search_t* const archive_search,
    const uint64_t idx_position,
    sequence_t* const sequence,
    matches_t* const matches) {
  // Clear matches & handlers
  archive_search_handlers_clear(archive_search_handlers);
  matches_clear(matches);
  // Prepare Search
  archive_search_handlers_prepare_se(archive_search,sequence,archive_search_handlers);
  // Archive Search
  archive_search_se(archive_search,matches);
  // DEBUG
  if (parameters.debug) {
    uint64_t i;
    // Searched sequence
    fprintf(stderr,"SEQ=");
    for (i=0;i<sequence->read.length;++i) {
      fprintf(stderr,"%c",sequence->read.buffer[i]);
    }
    fprintf(stderr,"\tpos=%lu\t",idx_position);
    mappability_print_counters(stderr,matches->counters,matches->max_complete_stratum);
    fprintf(stderr,"\n");
    // Matches
    // Traverse found matches
    const uint64_t num_matches = matches_get_num_match_traces(matches);
    match_trace_t** const match_traces = matches_get_match_traces(matches);
    for (i=0;i<num_matches;++i) {
      match_trace_t* const match_trace = match_traces[i];
      fprintf(stderr,"[#%lu](e=%lu,%lu:+:",i,match_trace->edit_distance,
          match_trace->match_alignment.match_position);
      const cigar_element_t* const cigar_buffer = match_trace_get_cigar_buffer(match_trace,matches->cigar_vector);
      const uint64_t cigar_length = match_trace_get_cigar_length(match_trace);
      mappability_print_cigar(stderr,cigar_buffer,cigar_length);
      fprintf(stderr,")\n");
    }
  }
}
/*
 * Mappability process found matches
 */
uint8_t vlogp32(
    const uint64_t* const log_table,
    const uint64_t x) {
  uint64_t index=0,new_index,increment=64;
  while (increment>0) {
    new_index = index+increment;
    if (new_index<94 && log_table[new_index]<=x) {
      index = new_index;
    }
    increment >>= 1;
  }
  assert(log_table[index]<=x);
  gem_cond_fatal_error(x>=log_table[index+1],MAPP_LOG,(uint64_t)x);
  return 32 + index;
}
void mappability_process_matches(
    vector_t* const queue_buf,
    const uint64_t current_idx_position,
    hashed_t** const hits,
    matches_t* const matches) {
  // Parameters
  const uint64_t max_errors = parameters.max_errors;
  const uint64_t approx_threshold = parameters.approx_threshold;
  const uint64_t eff_text_len = parameters.eff_text_len;
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  assert(num_matches!=0);
  // Count matches
  const uint64_t found = matches_get_num_match_traces(matches);
  const uint64_t found_exact = matches_counters_get_count(matches->counters,0);
  const uint8_t vlogp32_found = vlogp32(parameters.log_table,found);
  // Traverse found matches
  match_trace_t** const match_traces = matches_get_match_traces(matches);
  uint64_t i;
  for (i=0;i<num_matches;++i) {
    match_trace_t* const match_trace = match_traces[i];
    const uint64_t idx_position = match_trace->match_alignment.match_position;
    if ((idx_position < eff_text_len) && (idx_position >= current_idx_position)) {
      // Check position & exact mappability (then matches can be propagated)
      const bool exact_multiple_matches =
          (match_trace->edit_distance==0 && found_exact>=approx_threshold);
      const bool exact_forward_propagated =
          ((idx_position > current_idx_position) && (max_errors==0 || exact_multiple_matches));
      // Check position
      if (idx_position == current_idx_position || exact_forward_propagated) {
        // Add position to the queue
        bitmap_t* queue;
        vector_alloc_new(queue_buf,bitmap_t,queue);
        queue->position = idx_position;
        queue->vlog = (found>=approx_threshold) ? VLOGP32_0 : vlogp32_found;
        gem_cond_debug_block(DEBUG_MAPPABILITY_C) {
          fprintf(stderr,"Assigning value '%c' (%lu%s) to position %lu (mismatches=%ld,max_mismatches=%ld)\n",
              queue->vlog,found,found>=approx_threshold?"=inf":"",idx_position,
              match_trace->edit_distance,max_errors);
        }
        if (exact_forward_propagated) {
          // As for exact matches an equivalence relation applies, we can
          // safely assume that all the exact matches are also above the threshold
          hashed_t* const hit = calloc(1,sizeof(hashed_t));
          hit->position = idx_position;
          HASH_ADD_INT((*hits),position,hit); // Add the position to the hash table
        }
      }
    }
  }
  gem_cond_debug_block(DEBUG_MAPPABILITY_C) {
    fprintf(stderr,"There are %d hashed hits\n",HASH_COUNT((*hits)));
  }
}
/*
 * Mappability update global bitmap
 */
void mappability_update_bitmap(
    vector_t* const queue_buf,
    const uint64_t start_pos,
    const uint64_t end_pos) {
  // Parameters
  uint8_t* const bitmap = vector_get_mem(parameters.bitmap_buf,uint8_t);
  const uint64_t eff_text_len = parameters.eff_text_len;
  // Bitmap update
  pthread_mutex_lock(&parameters.bitmap_mutex);
  const uint64_t num_queue_elements = vector_get_used(queue_buf);
  bitmap_t* pos__vlog = vector_get_mem(queue_buf,bitmap_t);
  uint64_t i;
  for (i=0;i<num_queue_elements;++i,++pos__vlog) {
    const uint64_t pos = pos__vlog->position;
    const uint8_t vlog = pos__vlog->vlog;
    const uint8_t bm = bitmap[pos];
    if (!bm) {
      assert(pos < eff_text_len && vlog);
      gem_cond_debug_block(DEBUG_MAPPABILITY_C) {
        fprintf(stderr,"Assigning value '%c' to bitmap position %lu\n",vlog,pos);
      }
      // We just found out a new position
      bitmap[pos] = vlog;
      ++(parameters.done_new);
      if (vlog == VLOGP32_1) {
        ++(parameters.done_unique);
      }
    } else {
      assert(vlog == bm);
    }
  }
  /* We print updated statistics */
  gem_cond_debug_block(DEBUG_MAPPABILITY_C) {
    fprintf(stderr,"eff_text_len=%lu\tdone=%lu\tdone_new=%lu\tdone_unique=%lu\n",
        eff_text_len,end_pos,parameters.done_new,parameters.done_unique);
  }
  if (parameters.num_threads==1) {
    assert(parameters.done_new >= end_pos);
  }
  tprintf("Pos=%.3g%% Done=%lu(%.3g%%) Uniq=%.3g%%\n",
    (100.*end_pos)/eff_text_len,
    (uint64_t)parameters.done_new,
    (100.*parameters.done_new)/eff_text_len,
    (100.*parameters.done_unique)/parameters.done_new);
  pthread_mutex_unlock(&parameters.bitmap_mutex);
}
/*
 * Mappability thread
 */
uint64_t go_on(uint64_t* const start_pos,uint64_t* const end_pos) {
  pthread_mutex_lock(&parameters.remaining_mutex);
  const uint64_t res =
      (parameters.remaining>=parameters.granularity) ?
          parameters.granularity : parameters.remaining;
  *start_pos = parameters.eff_text_len - parameters.remaining;
  *end_pos = *start_pos + res;
  parameters.remaining -= res;
  pthread_mutex_unlock(&parameters.remaining_mutex);
  return res;
}
bool already_done(
    hashed_t* const hits,
    uint64_t pos) {
  hashed_t* found;
  HASH_FIND_INT(hits,&pos,found);
  return (found!=NULL);
}
void* mappability_thread(void* thread_id) {
  gem_thread_register_id(*((int*)thread_id)+1);
  // Parameters
  locator_t* const locator = parameters.archive->locator;
  archive_text_t* const text = parameters.archive->text;
  const uint64_t read_length = parameters.read_length;
  const uint64_t granularity = parameters.granularity;
  uint8_t* const bitmap = vector_get_mem(parameters.bitmap_buf,uint8_t);
  // Search structures
  archive_search_handlers_t* archive_search_handlers;
  archive_search_t* archive_search;
  matches_t* matches;
  mm_slab_t* mm_slab;
  mm_allocator_t* mm_allocator;
  // Init search structures
  archive_search_handlers = archive_search_handlers_new(parameters.archive);
  archive_search_se_new(&parameters.search_parameters,false,&archive_search);
  matches = matches_new();
  mm_slab = mm_slab_new(BUFFER_SIZE_8M);
  mm_allocator = mm_allocator_new(mm_slab);
  // Init position queue
  vector_t* const queue_buf = vector_new(granularity,bitmap_t);
  uint64_t start_pos, end_pos;
  // Mapping reference sequences
  while (go_on(&start_pos,&end_pos)) {
    // Parameters
    hashed_t* hits = NULL; // Hash table for hits
    vector_clear(queue_buf); // No new positions in the queue so far
    uint64_t current_idx_position = start_pos;
    while (current_idx_position < end_pos) {
      // Check character (Separators)
      if ((*dna_text_retrieve_sequence(text->enc_text,current_idx_position))==ENC_DNA_CHAR_SEP) {
        ++current_idx_position;
        continue;
      }
      // Fetch location
      locator_interval_t* const locator_interval = locator_lookup_interval(locator,current_idx_position);
      if (locator_interval->type == locator_interval_uncalled) {
        current_idx_position += (locator_interval->end_position - current_idx_position);
        continue;
      }
      // Process within current locator interval
      while (current_idx_position < end_pos &&
             current_idx_position < locator_interval->end_position) {
        // Retrieve reference sequence
        sequence_t sequence;
        uint64_t num_invalid_chars;
        mm_allocator_clear(mm_allocator);
        mappability_retrieve_sequence(
            text,current_idx_position,read_length,
            &sequence,&num_invalid_chars,mm_allocator);
        if (num_invalid_chars > 0) {
          // Make room for the position itself
          bitmap_t* queue;
          vector_alloc_new(queue_buf,bitmap_t,queue);
          queue->position = current_idx_position;
          queue->vlog = VLOGP32_0;
          gem_cond_debug_block(DEBUG_MAPPABILITY_C) {
            fprintf(stderr,"Assigning value '%c' (0) to position %lu (invalid=%lu)\n",
                    VLOGP32_0,current_idx_position,num_invalid_chars);
          }
        } else {
          // The position might have been already computed if approximation mode is on
          if (!bitmap[current_idx_position] && !already_done(hits,current_idx_position)) {
            // Mapping sequence
            mappability_mapping_sequence(
                archive_search_handlers,archive_search,
                current_idx_position,&sequence,matches);
            // Process matches
            mappability_process_matches(queue_buf,current_idx_position,&hits,matches);
          }
        }
        // Next
        ++current_idx_position;
      }
    }
    // Update global bitmap
    mappability_update_bitmap(queue_buf,start_pos,end_pos);
    // Delete hashtable
    hashed_t* hit, *tmp;
    HASH_ITER(hh,hits,hit,tmp) {
      HASH_DEL(hits,hit);
      free(hit);
    }
  }
  // Clean up
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(mm_slab);
  matches_delete(matches);
  archive_search_delete(archive_search);
  archive_search_handlers_delete(archive_search_handlers);
  vector_delete(queue_buf);
  // Exit
  pthread_exit(0);
}
/*
 * Mappability output frequencies
 */
void mappability_output_blanks(
    FILE* const output_file,
    const uint64_t total_bytes,
    const uint64_t output_line_width,
    uint64_t* const open_line_len) {
  size_t remaining_bytes = total_bytes;
  while (remaining_bytes) {
    const size_t eff_line_width = output_line_width - *open_line_len;
    size_t to_be_written;
    if (remaining_bytes >= eff_line_width) {
      to_be_written = eff_line_width;
      *open_line_len = 0;
    } else {
      to_be_written = remaining_bytes;
      *open_line_len += remaining_bytes;
    }
    int64_t j;
    for (j=0;j<to_be_written;++j) {
      gem_cond_fatal_error(fprintf(output_file,"%c",VLOGP32_0)!=1,FILE_WRITE);
    }
    gem_cond_fatal_error(!(*open_line_len) && fprintf(output_file,"\n")!=1,FILE_WRITE);
    remaining_bytes -= to_be_written;
  }
}
void mappability_output_write_by_columns(
    FILE* const output_file,
    uint8_t* const bitmap,
    const uint64_t total_bytes,
    const uint64_t output_line_width,
    uint64_t* const open_line_len) {
  size_t remaining_bytes = total_bytes;
  while (remaining_bytes) {
    const size_t eff_line_width = output_line_width - *open_line_len;
    size_t to_be_written;
    if (remaining_bytes >= eff_line_width) {
      to_be_written = eff_line_width;
      *open_line_len = 0;
    } else {
      to_be_written = remaining_bytes;
      *open_line_len += remaining_bytes;
    }
    gem_cond_fatal_error(
        fwrite(bitmap,sizeof(uint8_t),to_be_written,output_file)!=to_be_written,FILE_WRITE);
    gem_cond_fatal_error(!(*open_line_len) && fprintf(output_file,"\n")!=1,FILE_WRITE);
    remaining_bytes -= to_be_written;
  }
}
void mappability_output_frequencies(
  FILE* const output_file,
  vector_t* const bitmap_buf,
  const uint64_t read_length) {
  // Parameters
  uint8_t* bitmap = vector_get_mem(bitmap_buf,uint8_t);
  locator_t* const locator = parameters.archive->locator;
  char* last_tag = NULL;
  uint64_t last_end_position = 0;
  const uint64_t separators = 1; // Separator-related stuff
  uint64_t open_line_len = 0; // Size of the beginning of the line
  // Interval cycle
  tprintf("Writing frequencies to disk...\n");
  uint64_t i;
  for (i=0;i<locator->num_intervals;++i) {
    // Fetch interval
    locator_interval_t* const interval = locator->intervals + i;
    char* const interval_tag = locator_interval_get_tag(locator,interval);
    const uint64_t end_position = interval->end_position;
    uint64_t position = interval->begin_position;
    // If the interval is not the first one, we have to correct for the separator(s)
    if (i > 0) {
      position += separators;
      bitmap += separators;
    }
    const uint64_t interval_length = end_position - position;
    switch (interval->strand) {
    case Forward:
      if (interval_tag == last_tag) {
        // We insert blanks to connect segments
        mappability_output_blanks(
            output_file,position-last_end_position,
            parameters.output_line_width,&open_line_len);
      } else {
        // Nothing to insert here, in case we just truncate
        last_tag = interval_tag;
        if (open_line_len) {
          open_line_len = 0;
          gem_cond_fatal_error(fprintf(output_file,"\n")!=1,FILE_WRITE);
        }
        gem_cond_fatal_error(fprintf(output_file,"~%s\n",interval_tag)!=strlen(interval_tag)+2,FILE_WRITE);
        // There might be trailing 'N's here
        assert(position>=0);
        if (position > 0) {
          mappability_output_blanks(
              output_file,position,
              parameters.output_line_width,&open_line_len);
        }
      }
      last_end_position = end_position;
      // If the segment is too short w.r.t. the read length, we just pad it
      if (interval_length < read_length-1) {
        mappability_output_blanks(
            output_file,interval_length,
            parameters.output_line_width,&open_line_len);
        bitmap += interval_length;
      } else {
        // First we write part of the bitmap, then we pad with blanks to preserve the length
        const uint64_t eff_bitmap_len = interval_length+1-read_length;
        mappability_output_write_by_columns(
            output_file,bitmap,eff_bitmap_len,
            parameters.output_line_width,&open_line_len);
        bitmap += eff_bitmap_len;
        // We output the last part of the sequence as blanks
        mappability_output_blanks(
            output_file,read_length-1,
            parameters.output_line_width,&open_line_len);
        bitmap += (read_length-1);
      }
      break;
    case Reverse:
      // No need to write anything: the segment is same as the one before
      bitmap += interval_length;
      break;
    default:
      assert(0);
      break;
    }
  }
  assert(bitmap==((uint8_t*)bitmap_buf->memory)+archive_get_index_length(parameters.archive));
  gem_cond_fatal_error(open_line_len && fprintf(output_file,"\n")!=1,FILE_WRITE);
  tprintf("...done.\n");
}
/*
 * GEM-Mappability options Menu
 */
void gem_mappability_usage() {
  printf(
    "Usage:\n"
    " gem-mappability\n"
    "  -I <index_prefix>                            (mandatory)\n"
    "  --granularity <number>                       (default=%lu)\n"
    "  -o <output_prefix>                           (mandatory)\n"
    "  --output-line-width <number>                 (default=%lu)\n"
    "  -l <read_length>                             (mandatory)\n"
    "  --approximation-threshold <number>|'disable' (default: first multiple bin)\n"
    "  -e <max_edit_distance>                       (default=%lu)\n"
    "  -t <thread_number>                           (default=%lu)\n"
    "  -h|--help                                    (print usage)\n",
    parameters.granularity,
    parameters.output_line_width,
    parameters.max_errors,
    parameters.num_threads);
}
void gem_mappability_parse_arguments(
    const int argc,
    char** const argv,
    char* const gem_version) {
  bool done_I=false, done_o=false, done_l=false;
  struct option long_options[]={
    { "granularity", required_argument, 0, 0 },
    { "output-line-width", required_argument, 0, 1 },
    { "approximation-threshold", required_argument, 0, 2 },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 }
  };
  int c,option_index;
  while (1) {
    c = getopt_long(argc,argv,"I:o:l:m:e:t:gvh",long_options,&option_index);
    if (c == -1) break;
    switch (c) {
    case 0:
      parameters.granularity = atoll(optarg);
      break;
    case 1:
      parameters.output_line_width = atoll(optarg);
      break;
    case 'I':
      parameters.index_file_name = optarg;
      done_I = true;
      break;
    case 'o':
      parameters.output_file_name = optarg;
      done_o = true;
      break;
    case 'l':
      parameters.read_length = atoll(optarg);
      done_l = true;
      break;
    case 2:
      if (!strcmp(optarg,"disable")) {
        parameters.approx_threshold = UINT64_MAX;
      } else {
        parameters.approx_threshold = atoll(optarg);
      }
      break;
    case 'e':
      parameters.max_errors = atoll(optarg);
      gem_cond_fatal_error_msg(parameters.max_errors < 0,"-e must be positive");
      break;
    case 't':
      parameters.num_threads = atoll(optarg);
      break;
    case 'g':
      parameters.debug = true;
      break;
    case 'v':
      parameters.verbose = true;
      break;
    case 'h':
      gem_mappability_usage();
      exit(0);
      break;
    case '?':
    default:
      gem_error_msg("Option '%s' not recognized",optarg);
      gem_mappability_usage();
      exit(1);
      break;
    }
  }
  gem_cond_fatal_error_msg(!done_I||!done_o||!done_l,"Options -I, -o, -l are mandatory");
}
void mappability_parameters_init() {
  // Set defaults
  search_parameters_t* const search_parameters = &parameters.search_parameters;
  search_parameters_init(search_parameters);
  // Search Mode
  search_parameters->mapping_mode = mapping_hybrid_complete;
  search_parameters->match_alignment_model = match_alignment_model_levenshtein;
  search_parameters->alignment_local = local_alignment_never;
  // Search Error
  search_parameters->alignment_max_error = parameters.max_errors;
  search_parameters->alignment_max_bandwidth = parameters.max_errors;
  search_parameters->complete_search_error = parameters.max_errors;
  search_parameters->complete_strata_after_best = parameters.max_errors;
  // Matches
  search_parameters->select_parameters.max_searched_matches = 1000000;
  search_parameters->select_parameters.max_reported_matches = 1000000;
}
/*
 * GEM-Mappability
 */
int main(int argc,char** argv) {
  // Parsing command-line options
  gem_mappability_parse_arguments(argc,argv,gem_version); // Parse cmd-line
  gruntime_init(parameters.num_threads+1,NULL);
  // Init parameters
  mappability_parameters_init(&parameters);
  // Load index
  parameters.archive = mappability_load_index(parameters.index_file_name,parameters.verbose);
  // Compute logarithm table
  mappability_compute_log_table(parameters.read_length,parameters.log_table);
  // Output file
  FILE* const output = mappability_open_output(
      parameters.output_file_name,
      parameters.read_length,
      parameters.max_errors,
      parameters.approx_threshold,
      parameters.log_table);
  // Creating and zeroing the logarithm bitmap(s)
  parameters.bitmap_buf = vector_new(parameters.eff_text_len,uint8_t);
  memset(parameters.bitmap_buf->memory,0,parameters.eff_text_len*sizeof(uint8_t));
  // Setting global counters
  parameters.done_new = 0;
  parameters.done_unique = 0;
  parameters.remaining = parameters.eff_text_len;
  // Generate mappability threads
  tprintf("Starting (%lu positions to go)...\n",parameters.eff_text_len);
  const uint64_t num_threads = parameters.num_threads;
  pthread_t* threads = calloc(num_threads,sizeof(pthread_t));
  int* thread_ids = calloc(num_threads,sizeof(int));
  uint64_t i;
  for (i=0;i<num_threads;++i) {
    thread_ids[i] = i;
    gem_cond_fatal_error_msg(
        pthread_create(threads+i,0,mappability_thread,(void*)(thread_ids+i)),
        "Couldn't create mappability threads");
  }
  for (i=0;i<num_threads;++i) {
    pthread_join(threads[i],0);
  }
  // Threads work check & clean-up
  free(threads);
  free(thread_ids);
  assert(parameters.done_new==parameters.eff_text_len);
  tprintf("...done.\n");
  // Output mappability results
  mappability_output_frequencies(output,parameters.bitmap_buf,parameters.read_length);
  // Close files
  gem_cond_fatal_error_msg(fclose(output),"Couldn't close output file");
  archive_delete(parameters.archive);
  // Clean-up bitmap
  vector_delete(parameters.bitmap_buf);
  return 0;
}
