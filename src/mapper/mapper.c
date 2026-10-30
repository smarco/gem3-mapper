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
 * DESCRIPTION:
 *   Mapper module provides high-level functions to run
 *   the standard mapper workflow (SE/PE)
 */

#include "mapper/mapper.h"
#include "mapper/mapper_io.h"
#include "archive/search/archive_search_se.h"
#include "archive/search/archive_search_pe.h"
#include "text/sequence_bisulfite.h"
#include "io/output_sam.h"
#include "io/input_fasta.h"
#include "filtering/region_profile/region_profile_adaptive.h"
#include "profiler/profiler.h"

/*
 * Profile Level
 */
#define PROFILE_LEVEL PHIGH

/*
 * Debug
 */
//#define DEBUG_MAPPER_DISPLAY_EACH_READ_TIME

/*
 * Error Report
 */
mapper_search_t* g_mapper_searches; // Global searches on going
pthread_mutex_t mapper_error_report_mutex = PTHREAD_MUTEX_INITIALIZER;
void mapper_error_report_cmd(
    FILE* stream,
    mapper_parameters_t* const mapper_parameters) {
  // Display header
  fprintf(stream,"GEM::Version %s\n",mapper_parameters->gem_version);
  fprintf(stream,"GEM::CMD");
  // Print CMD line used
  uint64_t i;
  for (i=0;i<mapper_parameters->argc;++i) {
    fprintf(stream," %s",mapper_parameters->argv[i]);
  }
  fprintf(stream,"\n");
}
void mapper_error_report_input_state(
    FILE* stream,
    buffered_input_file_t* const buffered_fasta_input,
    const sequence_t* const sequence) {
  // Display header
  fprintf(stream,"GEM::Input.State\n");
  // Check NULL
  if (sequence==NULL) { fprintf(stream,"Sequence is NULL\n"); return; }
  if (buffered_fasta_input==NULL) { fprintf(stream,"Buffered_fasta_input is NULL\n"); return; }
  // Dump FASTA/FASTQ read
  if (!string_is_null(&sequence->tag) && !string_is_null(&sequence->read)) {
    const bool has_qualities = sequence_has_qualities(sequence);
    char* const end_tag =
        (sequence->end_info == paired_end1) ? "/1" :
      ( (sequence->end_info == paired_end2) ? "/2" : " " );
    fprintf(stream,"Sequence (File '%s' Line '%"PRIu64"')\n",
        buffered_input_file_get_file_name(buffered_fasta_input),
        buffered_fasta_input->current_buffer_line_no - (has_qualities ? 4 : 2));
    if (has_qualities) {
      if (!string_is_null(&sequence->qualities)) {
        fprintf(stream,"@%"PRIs"%s\n%"PRIs"\n+\n%"PRIs"\n",
            PRIs_content(&sequence->tag),end_tag,
            PRIs_content(&sequence->read),
            PRIs_content(&sequence->qualities));
      } else {
        fprintf(stream,"@%"PRIs"%s\n%"PRIs"\n+\n<<Null Qualities>>\n",
            PRIs_content(&sequence->tag),end_tag,
            PRIs_content(&sequence->read));
      }
    } else {
      fprintf(stream,">%"PRIs"%s\n%"PRIs"\n",
          PRIs_content(&sequence->tag),end_tag,
          PRIs_content(&sequence->read));
    }
  } else {
    fprintf(stream,"Current sequence is <<Empty>>\n");
  }
}
void mapper_error_report(FILE* stream) {
  // Select thread
  const uint64_t thread_id = gem_thread_get_thread_id();
  if (thread_id==0) {
    mapper_parameters_t* const mapper_parameters = g_mapper_searches->mapper_parameters;
    MUTEX_BEGIN_SECTION(mapper_parameters->error_report_mutex) {
      fprintf(stream,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
      fprintf(stream,"GEM::Unexpected error occurred. Sorry for the inconvenience\n"
                     "     Feedback and bug reporting it's highly appreciated,\n"
                     "     => Please report or email (gem.mapper.dev@gmail.com)\n");
      fprintf(stream,"GEM::Running-Thread (threadID = MASTER)\n");
      fprintf(stream,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
      mapper_error_report_cmd(stream,mapper_parameters); // Display CMD used
    } MUTEX_END_SECTION(mapper_parameters->error_report_mutex);
  } else {
    mapper_search_t* const mapper_search = g_mapper_searches + (thread_id-1); // Thread
    mapper_parameters_t* const mapper_parameters = mapper_search->mapper_parameters;
    MUTEX_BEGIN_SECTION(mapper_parameters->error_report_mutex) {
      fprintf(stream,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
      fprintf(stream,"GEM::Unexpected error occurred. Sorry for the inconvenience\n"
                     "     Feedback and bug reporting it's highly appreciated,\n"
                     "     => Please report or email (gem.mapper.dev@gmail.com)\n");
      fprintf(stream,"GEM::Running-Thread (threadID = %"PRIu64")\n",thread_id);
      fprintf(stream,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
      mapper_error_report_cmd(stream,mapper_parameters); // Display CMD used
      fprintf(stream,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
      // Display Input State
      mapper_io_handler_t* const mapper_io_handler = mapper_search->mapper_io_handler;
      if (!mapper_search->mapper_io_handler->paired_end) {
        const sequence_t* const sequence = *mapper_search->sequence_end1;
        mapper_error_report_input_state(stream,mapper_io_handler->buffered_fasta_input_end1,sequence);
      } else {
        const sequence_t* const sequence_end1 = *mapper_search->sequence_end1;
        const sequence_t* const sequence_end2 = *mapper_search->sequence_end2;
        mapper_error_report_input_state(stream,mapper_io_handler->buffered_fasta_input_end1,sequence_end1);
        mapper_error_report_input_state(stream,mapper_io_handler->buffered_fasta_input_end2,sequence_end2);
      }
    } MUTEX_END_SECTION(mapper_parameters->error_report_mutex);
  }
}
/*
 * Region profiles computed ahead of the reads' searches
 *   The keys of the reads that follow in the current input block are read from the block
 *   without consuming it, and profiled in lockstep with the current read's key. Each search
 *   replays its profile only when its own key turns out to be the one profiled.
 */
#define MAPPER_AHEAD_KEYS 16 // Keys profiled in lockstep (SE: reads; PE: ends of pairs)
typedef struct {
  bool enabled;
  uint64_t num_aheads;             // Profiles of the last batch, in input order (PE: end1,end2 per pair)
  uint64_t next_ahead;             // Next profile to hand out
  region_profile_ahead_t* aheads;  // MAPPER_AHEAD_KEYS profiles
  string_t line;                   // Peeked line
  string_t read_end1;              // Peeked reads
  string_t read_end2;
} mapper_ahead_t;
static void mapper_ahead_init(
    mapper_ahead_t* const ahead,
    mapper_parameters_t* const parameters,
    const bool paired_end) {
  const search_parameters_t* const search_parameters = &parameters->search_parameters;
  const archive_t* const archive = parameters->archive;
  ahead->enabled =
      MAPPER_AHEAD_KEYS >= 2 && MAPPER_AHEAD_KEYS <= REGION_PROFILE_AHEAD_MAX_KEYS &&
      search_parameters->mapping_mode == mapping_adaptive_filtering_fast &&
      search_parameters->region_profile_model.strategy == region_profile_adaptive &&
      !archive->text->run_length &&
      (paired_end || archive->type != archive_dna_bisulfite);
  ahead->num_aheads = 0;
  ahead->next_ahead = 0;
  if (!ahead->enabled) return;
  ahead->aheads = mm_calloc(MAPPER_AHEAD_KEYS,region_profile_ahead_t,true);
  string_init(&ahead->line,256,NULL);
  string_init(&ahead->read_end1,256,NULL);
  string_init(&ahead->read_end2,256,NULL);
}
static void mapper_ahead_destroy(mapper_ahead_t* const ahead) {
  if (!ahead->enabled) return;
  mm_free(ahead->aheads);
  string_destroy(&ahead->line);
  string_destroy(&ahead->read_end1);
  string_destroy(&ahead->read_end2);
}
static bool mapper_ahead_peek_read(
    buffered_input_file_t* const cursor,
    string_t* const line,
    string_t* const read) {
  // One FASTA/FASTQ record (false at the end of the block or on anything unexpected)
  if (buffered_input_file_eob(cursor)) return false;
  buffered_input_file_get_line(cursor,line);
  const char tag_begin = string_get_buffer(line)[0];
  if (tag_begin != FASTA_TAG_BEGIN && tag_begin != FASTQ_TAG_BEGIN) return false;
  // Read (normalized as input_fasta_parse_sequence_read does)
  if (buffered_input_file_eob(cursor)) return false;
  buffered_input_file_get_line(cursor,read);
  const uint64_t read_length = string_get_length(read)-1; // No EOL
  if (read_length == 0) return false;
  char* const buffer = string_get_buffer(read);
  uint64_t i;
  for (i=0;i<read_length-1;++i) {
    if (!is_dna(buffer[i])) buffer[i] = dna_normalized(buffer[i]);
  }
  string_set_length(read,read_length);
  // Separator & qualities
  if (tag_begin == FASTQ_TAG_BEGIN) {
    if (buffered_input_file_eob(cursor)) return false;
    buffered_input_file_get_line(cursor,line);
    if (string_get_buffer(line)[0] != FASTQ_SEP) return false;
    if (buffered_input_file_eob(cursor)) return false;
    buffered_input_file_get_line(cursor,line);
  }
  return true;
}
static void mapper_ahead_bisulfite_convert(
    string_t* const read,
    const char* const bs_table) {
  const uint64_t length = string_get_length(read);
  char* const buffer = string_get_buffer(read);
  uint64_t pos;
  for (pos=0;pos<length;pos++) buffer[pos] = bs_table[(int)buffer[pos]];
}
static bool mapper_ahead_bisulfite_flip(
    string_t* const read_end1,
    string_t* const read_end2) {
  // As sequence_bisulfite_check_cg_depletion_pe
  uint64_t counts1[7] = {0, 0, 0, 0, 0, 0, 0};
  uint64_t counts2[7] = {0, 0, 0, 0, 0, 0, 0};
  uint64_t length = string_get_length(read_end1);
  const char* buffer = string_get_buffer(read_end1);
  uint64_t pos;
  for (pos=0;pos<length;pos++) counts1[(int)dna_encode_table[(int)buffer[pos]]]++;
  length = string_get_length(read_end2);
  buffer = string_get_buffer(read_end2);
  for (pos=0;pos<length;pos++) counts2[(int)dna_encode_table[(int)buffer[pos]]]++;
  return (counts1[1] + counts2[2]) > (counts1[2] + counts2[1]);
}
static void mapper_ahead_prepare_key(
    region_profile_ahead_t* const ahead,
    const search_parameters_t* const parameters,
    string_t* const read) {
  // Key as pattern_init_encode computes it
  const uint64_t sequence_length = string_get_length(read);
  const char* const buffer = string_get_buffer(read);
  int64_t clip_left, clip_right, i;
  switch (parameters->clipping) {
    case clipping_disabled:
      clip_left = 0;
      break;
    case clipping_uncalled:
      i = 0;
      while (i<sequence_length && !is_dna_canonical(buffer[i])) ++i;
      clip_left = i;
      break;
    case clipping_masked:
      i = 0;
      while (i<sequence_length && !is_unmasked_dna(buffer[i])) ++i;
      clip_left = i;
      break;
    case clipping_fixed:
      clip_left = MIN(parameters->clip_left,sequence_length);
      break;
    default:
      ahead->ready = false;
      return;
  }
  switch (parameters->clipping) {
    case clipping_disabled:
      clip_right = 0;
      break;
    case clipping_uncalled:
      i = sequence_length-1;
      while (i>=clip_left && !is_dna_canonical(buffer[i])) --i;
      clip_right = sequence_length - (i+1);
      break;
    case clipping_masked:
      i = sequence_length-1;
      while (i>=clip_left && !is_unmasked_dna(buffer[i])) --i;
      clip_right = sequence_length - (i+1);
      break;
    case clipping_fixed:
      clip_right = MIN(parameters->clip_right,sequence_length - clip_left);
      break;
    default:
      ahead->ready = false;
      return;
  }
  const uint64_t key_length = sequence_length - clip_right - clip_left;
  if (!region_profile_ahead_prepare(ahead,key_length)) return;
  const char* const clipped = buffer + clip_left;
  uint64_t j;
  for (j=0;j<key_length;++j) ahead->key[j] = dna_encode(clipped[j]);
}
static void mapper_ahead_prepare_pattern(
    region_profile_ahead_t* const ahead,
    const pattern_t* const pattern) {
  if (!region_profile_ahead_prepare(ahead,pattern->key_length)) return;
  memcpy(ahead->key,pattern->key,pattern->key_length);
}
static void mapper_ahead_generate(
    mapper_ahead_t* const ahead,
    const uint64_t num_aheads,
    archive_search_t* const archive_search) {
  // Profile every key prepared, in lockstep
  region_profile_ahead_t* batch[MAPPER_AHEAD_KEYS];
  uint64_t i, num_batch = 0;
  for (i=0;i<num_aheads;++i) {
    if (ahead->aheads[i].key_length > 0 && ahead->aheads[i].region_profile.filtering_region != NULL) {
      batch[num_batch++] = ahead->aheads + i;
    }
  }
  region_profile_ahead_generate(batch,num_batch,
      archive_search->archive->fm_index,&archive_search->search_parameters.region_profile_model);
  ahead->num_aheads = num_aheads;
  ahead->next_ahead = 0;
}
static void mapper_ahead_batch_se(
    mapper_ahead_t* const ahead,
    archive_search_t* const archive_search,
    mapper_io_handler_t* const mapper_io_handler) {
  const search_parameters_t* const parameters = &archive_search->search_parameters;
  uint64_t i;
  for (i=0;i<MAPPER_AHEAD_KEYS;++i) {
    ahead->aheads[i].key_length = 0;
    ahead->aheads[i].ready = false;
  }
  // The current read
  mapper_ahead_prepare_pattern(ahead->aheads,&archive_search->approximate_search.pattern);
  // The reads that follow in the block
  buffered_input_file_t cursor = *(mapper_io_handler->buffered_fasta_input_end1);
  uint64_t num_aheads = 1;
  while (num_aheads < MAPPER_AHEAD_KEYS) {
    if (!mapper_ahead_peek_read(&cursor,&ahead->line,&ahead->read_end1)) break;
    mapper_ahead_prepare_key(ahead->aheads+num_aheads,parameters,&ahead->read_end1);
    ++num_aheads;
  }
  mapper_ahead_generate(ahead,num_aheads,archive_search);
}
static void mapper_ahead_batch_pe(
    mapper_ahead_t* const ahead,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    mapper_io_handler_t* const mapper_io_handler) {
  uint64_t i;
  for (i=0;i<MAPPER_AHEAD_KEYS;++i) {
    ahead->aheads[i].key_length = 0;
    ahead->aheads[i].ready = false;
  }
  // The current pair
  mapper_ahead_prepare_pattern(ahead->aheads,&archive_search_end1->approximate_search.pattern);
  mapper_ahead_prepare_pattern(ahead->aheads+1,&archive_search_end2->approximate_search.pattern);
  // The pairs that follow in the block(s)
  buffered_input_file_t cursor_end1 = *(mapper_io_handler->buffered_fasta_input_end1);
  buffered_input_file_t cursor_end2 = *(mapper_io_handler->buffered_fasta_input_end2);
  buffered_input_file_t* const cursor2 =
      mapper_io_handler->separated_input_files ? &cursor_end2 : &cursor_end1;
  uint64_t num_aheads = 2;
  while (num_aheads+2 <= MAPPER_AHEAD_KEYS) {
    if (!mapper_ahead_peek_read(&cursor_end1,&ahead->line,&ahead->read_end1)) break;
    if (!mapper_ahead_peek_read(cursor2,&ahead->line,&ahead->read_end2)) break;
    if (mapper_io_handler->archive_bisulfite) { // As sequence_bisulfite_process_pe
      bool flip = false;
      if (mapper_io_handler->bisulfite_read==bisulfite_non_stranded) {
        flip = mapper_ahead_bisulfite_flip(&ahead->read_end1,&ahead->read_end2);
      }
      mapper_ahead_bisulfite_convert(&ahead->read_end1,flip ? dna_bisulfite_G2A_table : dna_bisulfite_C2T_table);
      mapper_ahead_bisulfite_convert(&ahead->read_end2,flip ? dna_bisulfite_C2T_table : dna_bisulfite_G2A_table);
    }
    mapper_ahead_prepare_key(ahead->aheads+num_aheads,
        &archive_search_end1->search_parameters,&ahead->read_end1);
    mapper_ahead_prepare_key(ahead->aheads+num_aheads+1,
        &archive_search_end2->search_parameters,&ahead->read_end2);
    num_aheads += 2;
  }
  mapper_ahead_generate(ahead,num_aheads,archive_search_end1);
}
/*
 * SE Mapper
 */
void* mapper_se_thread(mapper_search_t* const mapper_search) {
  // GEM-thread error handler
  gem_thread_register_id(mapper_search->thread_id+1);

  // Parameters
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  search_parameters_t* const search_parameters = &mapper_search->mapper_parameters->search_parameters;
  archive_t* const archive = parameters->archive;
  mapper_io_handler_t* mapper_io_handler;
  archive_search_handlers_t* archive_search_handlers;
  archive_search_t* archive_search;
  matches_t* matches;
  sequence_t* sequence;

  // Init search structures
  archive_search_handlers = archive_search_handlers_new(archive);
  archive_search_se_new(search_parameters,false,&archive_search);
  matches = matches_new();

  // Init I/O handler
  mapper_io_handler = mapper_io_handler_new_se(
      parameters,parameters->io.input_buffer_size,
      archive_search_handlers->mm_allocator);

  // Set structures pointers (DEBUG)
  mapper_search->mapper_io_handler = mapper_io_handler;
  mapper_search->sequence_end1 = &sequence;

  // Region profiles computed ahead
  mapper_ahead_t ahead;
  mapper_ahead_init(&ahead,parameters,false);

  // FASTA/FASTQ reading loop
  uint64_t reads_processed = 0;
  while (mapper_read_sequence(mapper_io_handler,true,&sequence)) {
//    // DEBUG
//    if (gem_streq(sequence->tag.buffer,"Sim.Illumina.l100.0000009233")) {
//      printf("HERE\n");
//    }

    // Prepare Search
    archive_search_handlers_prepare_se(archive_search,sequence,archive_search_handlers);
    if (ahead.enabled) {
      if (ahead.next_ahead >= ahead.num_aheads) {
        mapper_ahead_batch_se(&ahead,archive_search,mapper_io_handler);
      }
      archive_search->approximate_search.region_profile_ahead = ahead.aheads + ahead.next_ahead;
      ++(ahead.next_ahead);
    }

    // Search into the archive
#ifdef DEBUG_MAPPER_DISPLAY_EACH_READ_TIME
    gem_timer_t timer;
    TIMER_RESTART(&timer); archive_search_se(archive_search,matches); TIMER_STOP(&timer);
    fprintf(stderr,"Done %s in %2.4f ms.\n",sequence->tag.buffer,TIMER_GET_TOTAL_MS(&timer));
#else
    archive_search_se(archive_search,matches);
#endif

    // Output matches
    mapper_io_handler_output_matches(mapper_io_handler,
        archive_search,matches,mapper_search->mapping_stats);

    // Update processed
    if (++reads_processed == parameters->io.mapper_ticker_step) {
      ticker_update_mutex(mapper_search->ticker,reads_processed);
      reads_processed=0;
    }

    // Clear (these are not needed => just wipe clean all mm_allocator)
    // archive_search_destroy(archive_search);
    // filtering_candidates_clear(&archive_search_handlers->filtering_candidates_end1,true);
    // Clear
    archive_search_handlers_clear(archive_search_handlers);
    matches_clear(matches);
  }
  // Update processed
  ticker_update_mutex(mapper_search->ticker,reads_processed);
  // Clean up
  mapper_ahead_destroy(&ahead);
  matches_delete(matches);
  archive_search_delete(archive_search);
  archive_search_handlers_delete(archive_search_handlers);
  mapper_io_handler_delete(mapper_io_handler);
  pthread_exit(0);
}
/*
 * PE Mapper
 */
void* mapper_pe_thread(mapper_search_t* const mapper_search) {
  // GEM-thread error handler
  gem_thread_register_id(mapper_search->thread_id+1);

  // Parameters
  mapper_parameters_t* const parameters = mapper_search->mapper_parameters;
  search_parameters_t* const search_parameters = &mapper_search->mapper_parameters->search_parameters;
  archive_t* const archive = parameters->archive;
  mapper_io_handler_t* mapper_io_handler;
  archive_search_handlers_t* archive_search_handlers;
  archive_search_t* archive_search_end1;
  archive_search_t* archive_search_end2;
  paired_matches_t* paired_matches;
  sequence_t* sequence_end1;
  sequence_t* sequence_end2;

  // Init search structures
  archive_search_handlers = archive_search_handlers_new(archive);
  archive_search_pe_new(search_parameters,false,&archive_search_end1,&archive_search_end2);
  paired_matches = paired_matches_new();

  // Init I/O handler
  mapper_io_handler = mapper_io_handler_new_pe(
      parameters,parameters->io.input_buffer_size,
      archive_search_handlers->mapper_stats,
      archive_search_handlers->mm_allocator);

  // Set structures pointers (DEBUG)
  mapper_search->mapper_io_handler = mapper_io_handler;
  mapper_search->sequence_end1 = &sequence_end1;
  mapper_search->sequence_end2 = &sequence_end2;

  // Region profiles computed ahead
  mapper_ahead_t ahead;
  mapper_ahead_init(&ahead,parameters,true);

  // FASTA/FASTQ reading loop
  uint64_t reads_processed = 0;
  while (mapper_read_paired_sequence(mapper_io_handler,true,&sequence_end1,&sequence_end2)) {
//    // DEBUG
//    if (gem_streq(sequence_end1->tag.buffer,"Sim.Illumina.l100.0000491385/1")) {
//      printf("HERE\n");
//    }

    // Prepare Search
    archive_search_handlers_prepare_pe(
        archive_search_end1,archive_search_end2,
        sequence_end1,sequence_end2,
        archive_search_handlers);
    if (ahead.enabled) {
      if (ahead.next_ahead >= ahead.num_aheads) {
        mapper_ahead_batch_pe(&ahead,archive_search_end1,archive_search_end2,mapper_io_handler);
      }
      archive_search_end1->approximate_search.region_profile_ahead = ahead.aheads + ahead.next_ahead;
      archive_search_end2->approximate_search.region_profile_ahead = ahead.aheads + ahead.next_ahead + 1;
      ahead.next_ahead += 2;
    }

    // Search into the archive
#ifdef DEBUG_MAPPER_DISPLAY_EACH_READ_TIME
    gem_timer_t timer;
    TIMER_RESTART(&timer); archive_search_pe(archive_search_end1,archive_search_end2,paired_matches); TIMER_STOP(&timer);
    fprintf(stderr,"Done %s in %2.4f ms.\n",sequence_end1->tag.buffer,TIMER_GET_TOTAL_MS(&timer));
#else
    archive_search_pe(archive_search_end1,archive_search_end2,paired_matches);
#endif

    // Output matches
    mapper_io_handler_output_paired_matches(
        mapper_io_handler,archive_search_end1,archive_search_end2,
        paired_matches,mapper_search->mapping_stats);
    //output_fastq(mapper_search->buffered_output_file,&archive_search_end1->sequence);
    //output_fastq(mapper_search->buffered_output_file,&archive_search_end2->sequence);

    // Update processed
    if (++reads_processed == parameters->io.mapper_ticker_step) {
      ticker_update_mutex(mapper_search->ticker,reads_processed);
      reads_processed=0;
    }

    // Clear (these are not needed => just wipe clean all mm_allocator)
//    archive_search_destroy(archive_search_end1);
//    archive_search_destroy(archive_search_end2);
//    filtering_candidates_clear(&archive_search_handlers->filtering_candidates_end1,true);
//    filtering_candidates_clear(&archive_search_handlers->filtering_candidates_end2,true);
    // Clear
    archive_search_handlers_clear(archive_search_handlers);
    paired_matches_clear(paired_matches,true);
  }
  // Update processed
  ticker_update_mutex(mapper_search->ticker,reads_processed);

  // Clean up
  mapper_ahead_destroy(&ahead);
  archive_search_delete(archive_search_end1);
  archive_search_delete(archive_search_end2);
  paired_matches_delete(paired_matches);
  archive_search_handlers_delete(archive_search_handlers);
  mapper_io_handler_delete(mapper_io_handler);
  pthread_exit(0);
}
/*
 * SE/PE runnable
 */
void mapper_run(mapper_parameters_t* const mapper_parameters,const bool paired_end) {
  // Load GEM-Index
  mapper_load_index(mapper_parameters);
  gem_cond_fatal_error_msg(
      paired_end && mapper_parameters->archive->text->run_length,
      "Archive RL-text not supported for paired-end mode (use standard index)");
  gem_cond_fatal_error_msg(
      paired_end && (mapper_parameters->archive->type == archive_dna_forward),
      "Archive no-complement not supported for paired-end mode (use standard index)");
  if (mapper_parameters->archive->gpu_index) {
    gem_warn_msg(
        "Running CPU-mode using GPU-index; "
        "CPU-index can achieve better performance "
        "(gem-indexer --gpu-index=false)");
  }
  // Setup threads
  const uint64_t num_threads = mapper_parameters->system.num_threads;
  mapper_search_t* const mapper_search = mm_calloc(num_threads,mapper_search_t,false);
  // Set error-report function
  g_mapper_searches = mapper_search;
  MUTEX_INIT(mapper_parameters->error_report_mutex);
  gem_error_set_report_function(mapper_error_report);
  // Prepare output file/parameters (SAM headers)
  archive_t* const archive = mapper_parameters->archive;
  const bool bisulfite_index = (archive->type == archive_dna_bisulfite);
  if (mapper_parameters->io.output_format==SAM) {
    output_sam_print_header(
        mapper_parameters->output_file,archive,&mapper_parameters->io.sam_parameters,
        mapper_parameters->argc,mapper_parameters->argv,mapper_parameters->gem_version);
    mapper_parameters->io.sam_parameters.bisulfite_output = bisulfite_index;
  }
  // Setup Ticker
  ticker_t ticker;
  ticker_count_reset(&ticker,mapper_parameters->misc.verbose_user,
      paired_end ? "PE::Mapping Sequences" : "SE::Mapping Sequences",0,
      mapper_parameters->io.mapper_ticker_step,false);
  ticker_add_process_label(&ticker,"#","sequences processed");
  ticker_add_finish_label(&ticker,"Total","sequences processed");
  ticker_mutex_enable(&ticker);
	// Allocate per thread mapping stats
  mapping_stats_t* const mstats = mapper_parameters->global_mapping_stats ?
      mm_calloc(num_threads,mapping_stats_t,false) : NULL;
  // Launch threads
  PROF_START(GP_MAPPER_MAPPING);
  pthread_handler_t mapper_thread;
  if (paired_end) {
    mapper_thread = (pthread_handler_t) mapper_pe_thread;
  } else {
    mapper_thread = (pthread_handler_t) mapper_se_thread;
  }
  uint64_t i;
  PROFILE_VTUNE_START(); // Vtune
  for (i=0;i<num_threads;++i) {
    // Setup thread
    mapper_search[i].thread_id = i;
    mapper_search[i].thread_data = mm_alloc(pthread_t);
    mapper_search[i].mapper_parameters = mapper_parameters;
    mapper_search[i].ticker = &ticker;
		if (mstats)	{
			 mapper_search[i].mapping_stats = mstats + i;
			 init_mapping_stats(mstats + i);
		} else {
		  mapper_search[i].mapping_stats = NULL;
		}
    // Launch thread
    gem_cond_fatal_error(
        pthread_create(mapper_search[i].thread_data,0,
            mapper_thread,(void*)(mapper_search+i)),SYS_THREAD_CREATE);
  }
  // Join all threads
  for (i=0;i<num_threads;++i) {
    gem_cond_fatal_error(pthread_join(*(mapper_search[i].thread_data),0),SYS_THREAD_JOIN);
    mm_free(mapper_search[i].thread_data);
  }
  PROFILE_VTUNE_STOP(); // Vtune
  ticker_finish(&ticker);
  ticker_mutex_cleanup(&ticker);
  PROF_STOP(GP_MAPPER_MAPPING);
	// Merge report stats
	if (mstats) {
		 merge_mapping_stats(mapper_parameters->global_mapping_stats,mstats,num_threads);
		 mm_free(mstats);
	}
  // Clean up
  mm_free(mapper_search);
  //MUTEX_DESTROY(mapper_parameters->error_report_mutex);
}
void mapper_se_run(mapper_parameters_t* const mapper_parameters) {
  mapper_run(mapper_parameters,false);
}
void mapper_pe_run(mapper_parameters_t* const mapper_parameters) {
  mapper_run(mapper_parameters,true);
}
