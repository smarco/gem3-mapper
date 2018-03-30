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
 *   Mapper module provides high-level functions for I/O
 */

#include "io/input_fasta.h"
#include "mapper/mapper_io.h"
#include "text/sequence_bisulfite.h"
#include "stats/report_stats.h"

/*
 * Profile Level
 */
#define PROFILE_LEVEL PHIGH

/*
 * Debug/Profile
 */
#define MAPPER_OUTPUT

/*
 * Error Messages
 */
#define GEM_ERROR_MAPPER_IO_RELOAD_BUFFER_SINGLE "Parsing Input File. Error reloading buffer from file '%s'"
#define GEM_ERROR_MAPPER_IO_RELOAD_BUFFER_PAIR "Parsing Input Files. Error reloading buffer from files '%s,%s'"

#define GEM_ERROR_MAPPER_IO_PARSE_SEQUENCE_SINGLE "Parsing Input File. Error parsing sequence from file '%s'"
#define GEM_ERROR_MAPPER_IO_PARSE_SEQUENCE_PAIR "Parsing Input Files. Error parsing sequence from files '%s,%s'"

#define GEM_ERROR_MAPPER_IO_UNSYNCH_PAIR "Parsing Input Files. Files '%s,%s' doesn't contain the same number of reads (cannot pair)"
#define GEM_ERROR_MAPPER_IO_UNSYNCH_SINGLE "Parsing Input File. File '%s' doesn't contain a pair number of reads (cannot pair)"
#define GEM_ERROR_MAPPER_IO_UNSYNCH_EOF "Parsing Input Files. File '%s' could not read second end (unexpected end-of-file)"
#define GEM_ERROR_MAPPER_IO_UNSYNCH_NOT_EOF "Parsing Input Files. File '%s' has too many reads (expected end-of-file)"

/*
 * Error Macros
 */
#define MAPPER_ERROR_IO_RELOAD_BUFFER(mapper_io_handler) \
  if (mapper_io_handler->separated_input_files) { \
    gem_fatal_error(MAPPER_IO_RELOAD_BUFFER_PAIR, \
        input_file_sliced_get_file_name(mapper_io_handler->buffered_fasta_input_end1->input_file_sliced), \
        input_file_sliced_get_file_name(mapper_io_handler->buffered_fasta_input_end2->input_file_sliced)); \
  } else { \
    gem_fatal_error(MAPPER_IO_RELOAD_BUFFER_SINGLE, \
        input_file_sliced_get_file_name(mapper_io_handler->buffered_fasta_input_end1->input_file_sliced)); \
  }
#define MAPPER_ERROR_IO_UNSYNCH(mapper_io_handler) \
  if (mapper_io_handler->separated_input_files) { \
    gem_fatal_error(MAPPER_IO_UNSYNCH_PAIR, \
        input_file_sliced_get_file_name(mapper_io_handler->buffered_fasta_input_end1->input_file_sliced), \
        input_file_sliced_get_file_name(mapper_io_handler->buffered_fasta_input_end2->input_file_sliced)); \
  } else { \
    gem_fatal_error(MAPPER_IO_UNSYNCH_SINGLE, \
        input_file_sliced_get_file_name(mapper_io_handler->buffered_fasta_input_end1->input_file_sliced)); \
  }
#define MAPPER_ERROR_IO_PARSE_SEQUENCE(mapper_io_handler) \
  if (mapper_io_handler->separated_input_files) { \
    gem_fatal_error(MAPPER_IO_PARSE_SEQUENCE_PAIR, \
        input_file_sliced_get_file_name(mapper_io_handler->buffered_fasta_input_end1->input_file_sliced), \
        input_file_sliced_get_file_name(mapper_io_handler->buffered_fasta_input_end2->input_file_sliced)); \
  } else { \
    gem_fatal_error(MAPPER_IO_PARSE_SEQUENCE_SINGLE, \
        input_file_sliced_get_file_name(mapper_io_handler->buffered_fasta_input_end1->input_file_sliced)); \
  }

/*
 * Index loader
 */
void mapper_load_index(mapper_parameters_t* const parameters) {
  TIMER_START(&parameters->loading_time);
  PROFILE_START(GP_MAPPER_LOAD_INDEX,PROFILE_LEVEL);
  // Load archive
  gem_cond_log(parameters->misc.verbose_user,"[Loading GEM index '%s']",parameters->io.index_file_name);
  parameters->archive = archive_read(parameters->io.index_file_name,false);
  if (parameters->misc.verbose_dev) archive_print(gem_error_get_stream(),parameters->archive);
  PROFILE_STOP(GP_MAPPER_LOAD_INDEX,PROFILE_LEVEL);
  TIMER_STOP(&parameters->loading_time);
}
/*
 * Setup
 */
mapper_io_handler_t* mapper_io_handler_new_se(
    mapper_parameters_t* const parameters,
    const uint64_t input_buffer_size,
    mm_allocator_t* const mm_allocator) {
  // Allocate
  mapper_io_handler_t* const mapper_io_handler = mm_alloc(mapper_io_handler_t);
  // Attributes
  mapper_io_handler->paired_end = false;
  mapper_io_handler->mapper_parameters_io = &parameters->io;
  mapper_io_handler->separated_input_files = false;
  // Bisulfite
  mapper_io_handler->archive_bisulfite = (parameters->archive->type==archive_dna_bisulfite);
  mapper_io_handler->bisulfite_read = parameters->search_parameters.bisulfite_read;
  // Create buffered I/O files
  mapper_io_handler->buffered_fasta_input_end1 =
      buffered_input_file_new(parameters->input_file,input_buffer_size);
  mapper_io_handler->buffered_output_file = buffered_output_file_new(parameters->output_file);
  buffered_input_file_attach_buffered_output(
      mapper_io_handler->buffered_fasta_input_end1,mapper_io_handler->buffered_output_file);
  // MM
  mapper_io_handler->mm_allocator = mm_allocator;
  // Return
  return mapper_io_handler;
}
mapper_io_handler_t* mapper_io_handler_new_pe(
    mapper_parameters_t* const parameters,
    const uint64_t input_buffer_size,
    mapper_stats_t* const mapper_stats,
    mm_allocator_t* const mm_allocator) {
  // Allocate
  mapper_io_handler_t* const mapper_io_handler = mm_alloc(mapper_io_handler_t);
  // Attributes
  mapper_io_handler->paired_end = true;
  mapper_io_handler->archive_bisulfite = (parameters->archive->type==archive_dna_bisulfite);
  mapper_io_handler->mapper_parameters_io = &parameters->io;
  search_paired_parameters_t* const paired_parameters =
      &parameters->search_parameters.search_paired_parameters;
  mapper_io_handler->separated_input_files = parameters->io.separated_input_files;
  mapper_io_handler->template_length_estimation_samples = paired_parameters->template_length_estimation_samples;
  mapper_io_handler->template_length_estimation_min = paired_parameters->template_length_estimation_min;
  mapper_io_handler->template_length_estimation_max = paired_parameters->template_length_estimation_max;
  // Stats
  mapper_io_handler->mapper_stats = mapper_stats;
  // Create buffered I/O files
  if (mapper_io_handler->separated_input_files) {
    mapper_io_handler->buffered_fasta_input_end1 =
        buffered_input_file_new(parameters->input_file_end1,input_buffer_size);
    mapper_io_handler->buffered_fasta_input_end2 =
        buffered_input_file_new(parameters->input_file_end2,input_buffer_size);
  } else {
    mapper_io_handler->buffered_fasta_input_end1 =
        buffered_input_file_new(parameters->input_file,input_buffer_size);
    mapper_io_handler->buffered_fasta_input_end2 = mapper_io_handler->buffered_fasta_input_end1;
  }
  mapper_io_handler->buffered_output_file = buffered_output_file_new(parameters->output_file);
  buffered_input_file_attach_buffered_output(
      mapper_io_handler->buffered_fasta_input_end1,mapper_io_handler->buffered_output_file);
  // Mutex
  mapper_io_handler->input_mutex = &parameters->input_file_mutex;
  // MM
  mapper_io_handler->mm_allocator = mm_allocator;
  // Return
  return mapper_io_handler;
}
void mapper_io_handler_delete(
    mapper_io_handler_t* const mapper_io_handler) {
  // Close  buffered I/O files
  buffered_input_file_close(mapper_io_handler->buffered_fasta_input_end1);
  if (mapper_io_handler->separated_input_files) {
    buffered_input_file_close(mapper_io_handler->buffered_fasta_input_end2);
  }
  buffered_output_file_close(mapper_io_handler->buffered_output_file);
  // Free handler
  mm_free(mapper_io_handler);
}
/*
 * Pair-End Reload/Parse Helpers
 */
uint64_t mapper_se_reload_buffer(
    mapper_io_handler_t* const mapper_io_handler) {
  return buffered_input_file_reload(mapper_io_handler->buffered_fasta_input_end1,0);
}
error_code_t mapper_se_parse_sequence(
    mapper_io_handler_t* const mapper_io_handler,
    sequence_t** const sequence) {
  return input_fasta_parse_sequence(mapper_io_handler->buffered_fasta_input_end1,*sequence,false);
}
uint64_t mapper_pe_reload_buffer(
    mapper_io_handler_t* const mapper_io_handler) {
  // Check end-of-block
  error_code_t error_code;
  // Reset template-length estimation
  mapper_stats_template_init(
      mapper_io_handler->mapper_stats,
      mapper_io_handler->template_length_estimation_samples,
      mapper_io_handler->template_length_estimation_min,
      mapper_io_handler->template_length_estimation_max);
  // Read new input-block
  if (!mapper_io_handler->separated_input_files) {
    error_code = buffered_input_file_reload(mapper_io_handler->buffered_fasta_input_end1,0);
    if (error_code==INPUT_STATUS_EOF) return INPUT_STATUS_EOF;
    const uint64_t buffered_lines = mapper_io_handler->buffered_fasta_input_end1->num_lines;
    if ((buffered_lines % 8) != 0) {
      MAPPER_ERROR_IO_UNSYNCH(mapper_io_handler);
    }
  } else {
    // Dump buffer (explicitly before synch-reload)
    buffered_output_file_dump_buffer(
        mapper_io_handler->buffered_fasta_input_end1->attached_buffered_output_file);
    // Reload buffers (in synch)
    gem_cond_fatal_error(pthread_mutex_lock(mapper_io_handler->input_mutex),SYS_MUTEX);
    // Check synch
    if (!buffered_input_file_eob(mapper_io_handler->buffered_fasta_input_end2)) {
      MAPPER_ERROR_IO_UNSYNCH(mapper_io_handler);
    }
    // Reload end1
    error_code = buffered_input_file_reload(mapper_io_handler->buffered_fasta_input_end1,0);
    // Reload end2
    const uint64_t num_lines_read = mapper_io_handler->buffered_fasta_input_end1->num_lines;
    error_code = buffered_input_file_reload(
        mapper_io_handler->buffered_fasta_input_end2,num_lines_read);
    if (mapper_io_handler->buffered_fasta_input_end2->num_lines != num_lines_read) {
      MAPPER_ERROR_IO_UNSYNCH(mapper_io_handler);
    }
    if (num_lines_read==0) {
      gem_cond_fatal_error(pthread_mutex_unlock(mapper_io_handler->input_mutex),SYS_MUTEX);
      return INPUT_STATUS_EOF;
    }
    gem_cond_fatal_error(pthread_mutex_unlock(mapper_io_handler->input_mutex),SYS_MUTEX);
  }
  // OK
  return INPUT_STATUS_OK;
}
error_code_t mapper_pe_parse_sequence(
    mapper_io_handler_t* const mapper_io_handler,
    sequence_t** const sequence_end1,
    sequence_t** const sequence_end2) {
  // Parameters
  mm_allocator_t* const mm_allocator = mapper_io_handler->mm_allocator;
  // Allocate sequences
  *sequence_end1 = mm_allocator_alloc(mm_allocator,sequence_t);
  *sequence_end2 = mm_allocator_alloc(mm_allocator,sequence_t);
  sequence_init(*sequence_end1,mapper_io_handler->archive_bisulfite,mm_allocator);
  sequence_init(*sequence_end2,mapper_io_handler->archive_bisulfite,mm_allocator);
  // Read end1
  error_code_t error_code;
  error_code = input_fasta_parse_sequence(
      mapper_io_handler->buffered_fasta_input_end1,*sequence_end1,false);
  if (gem_expect_false(error_code!=INPUT_STATUS_OK)) return error_code;
  // Read end2
  error_code = input_fasta_parse_sequence(
      mapper_io_handler->buffered_fasta_input_end2,*sequence_end2,false);
  if (gem_expect_false(error_code!=INPUT_STATUS_OK)) return error_code;
  // OK
  return INPUT_STATUS_OK;
}
/*
 * Sequence readers
 */
error_code_t mapper_read_sequence(
    mapper_io_handler_t* const mapper_io_handler,
    const bool reload_input_buffer,
    sequence_t** const sequence) {
  // Parameters
  mm_allocator_t* const mm_allocator = mapper_io_handler->mm_allocator;
  // Check Buffer
  error_code_t error_code;
  if (buffered_input_file_eob(mapper_io_handler->buffered_fasta_input_end1)) {
    // Check reload-buffer flag
    if (!reload_input_buffer) return 0;
    // Reload buffer
    error_code = mapper_se_reload_buffer(mapper_io_handler);
    if (error_code==INPUT_STATUS_EOF) return INPUT_STATUS_EOF;
    if (error_code==INPUT_STATUS_FAIL) {
      MAPPER_ERROR_IO_RELOAD_BUFFER(mapper_io_handler);
    }
  }
  // Allocate sequence
  *sequence = mm_allocator_alloc(mm_allocator,sequence_t);
  sequence_init(*sequence,mapper_io_handler->archive_bisulfite,mm_allocator);
  // Parse sequence
  error_code = mapper_se_parse_sequence(mapper_io_handler,sequence);
  if (error_code==INPUT_STATUS_FAIL) {
    MAPPER_ERROR_IO_PARSE_SEQUENCE(mapper_io_handler);
  }
  // Bisulfite: Fully convert reads before searching into archive, making a copy of the original
  if (mapper_io_handler->archive_bisulfite) {
    sequence_bisulfite_process_se(*sequence,mapper_io_handler->bisulfite_read);
  }
  // Return
  PROF_INC_COUNTER(GP_MAPPER_NUM_READS);
  return error_code;
}
error_code_t mapper_read_paired_sequence(
    mapper_io_handler_t* const mapper_io_handler,
    const bool reload_input_buffer,
    sequence_t** const sequence_end1,
    sequence_t** const sequence_end2) {
  error_code_t error_code;
  // Check buffers
  if (buffered_input_file_eob(mapper_io_handler->buffered_fasta_input_end1)) {
    // Check reload-buffer flag
    if (!reload_input_buffer) return 0;
    // Reload buffer
    error_code = mapper_pe_reload_buffer(mapper_io_handler);
    if (error_code==INPUT_STATUS_EOF) return INPUT_STATUS_EOF;
    if (error_code==INPUT_STATUS_FAIL) {
      MAPPER_ERROR_IO_RELOAD_BUFFER(mapper_io_handler);
    }
  }
  // Read end1 & end2
  error_code = mapper_pe_parse_sequence(mapper_io_handler,sequence_end1,sequence_end2);
  if (error_code==INPUT_STATUS_FAIL) {
    MAPPER_ERROR_IO_PARSE_SEQUENCE(mapper_io_handler);
  }
  // Bisulfite: Fully convert reads before searching into archive, making a copy of the original
  if (mapper_io_handler->archive_bisulfite) {
    sequence_bisulfite_process_pe(*sequence_end1,*sequence_end2,mapper_io_handler->bisulfite_read);
  }
  // Return
  PROF_ADD_COUNTER(GP_MAPPER_NUM_READS,2);
  return INPUT_STATUS_OK;
}
/*
 * Output
 */
void mapper_io_handler_output_matches(
    mapper_io_handler_t* const mapper_io_handler,
    archive_search_t* const archive_search,
    matches_t* const matches,
    mapping_stats_t* const mstats) {
  // Bisulfite: Copy back original read
  if (mapper_io_handler->archive_bisulfite) {
    sequence_bisulfite_restore_se(archive_search->sequence);
  }
  // Output
#ifdef MAPPER_OUTPUT
  mapper_parameters_io_t* const mapper_parameters_io = mapper_io_handler->mapper_parameters_io;
  switch (mapper_parameters_io->output_format) {
    case MAP:
      output_map_single_end_matches(
          mapper_io_handler->buffered_output_file,
          archive_search,matches,&mapper_parameters_io->map_parameters);
      break;
    case SAM:
      output_sam_single_end_matches(
          mapper_io_handler->buffered_output_file,
          archive_search,matches,&mapper_parameters_io->sam_parameters);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
#endif
  // Stats
  if (mstats) {
    collect_se_mapping_stats(archive_search,matches,mstats);
  }
}
void mapper_io_handler_output_paired_matches(
    mapper_io_handler_t* const mapper_io_handler,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,
    mapping_stats_t* const mstats) {
  // Bisulfite: Copy back original read
  if (mapper_io_handler->archive_bisulfite) {
    sequence_bisulfite_restore_pe(
        archive_search_end1->sequence,archive_search_end2->sequence);
  }
  // Output
#ifdef MAPPER_OUTPUT
  mapper_parameters_io_t* const mapper_parameters_io = mapper_io_handler->mapper_parameters_io;
  switch (mapper_parameters_io->output_format) {
    case MAP:
      output_map_paired_end_matches(
          mapper_io_handler->buffered_output_file,archive_search_end1,
          archive_search_end2,paired_matches,&mapper_parameters_io->map_parameters);
      break;
    case SAM:
      output_sam_paired_end_matches(
          mapper_io_handler->buffered_output_file,archive_search_end1,
          archive_search_end2,paired_matches,&mapper_parameters_io->sam_parameters);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
#endif
  // Stats
  if (mstats) {
    collect_pe_mapping_stats(archive_search_end1,archive_search_end2,paired_matches,mstats);
  }
}
