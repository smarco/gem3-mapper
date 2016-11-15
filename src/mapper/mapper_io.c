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

#include "mapper/mapper_io.h"
#include "stats/report_stats.h"
#include "io/input_fasta_parser.h"

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
#define GEM_ERROR_MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_PAIR "Parsing Input Files. Files '%s,%s' doesn't contain the same number of reads (cannot pair)"
#define GEM_ERROR_MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_SINGLE "Parsing Input File. File '%s' doesn't contain a pair number of reads (cannot pair)"
#define GEM_ERROR_MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_EOF "Parsing Input Files. File '%s' could not read second end (unexpected end-of-file)"
#define GEM_ERROR_MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_NOT_EOF "Parsing Input Files. File '%s' has too many reads (expected end-of-file)"

/*
 * Error Macros
 */
#define MAPPER_ERROR_PE_PARSE_UNSYNCH_INPUT_FILES(mapper_parameters) \
  if (mapper_parameters->io.separated_input_files) { \
    gem_fatal_error(MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_PAIR, \
        input_file_sliced_get_file_name(mapper_parameters->input_file_end1), \
        input_file_sliced_get_file_name(mapper_parameters->input_file_end2)); \
  } else { \
    gem_fatal_error(MAPPER_PE_PARSE_UNSYNCH_INPUT_FILES_SINGLE, \
        input_file_sliced_get_file_name(mapper_parameters->input_file)); \
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
 * Input (Low-level)
 */
void mapper_se_prepare_io_buffers(
    const mapper_parameters_t* const parameters,
    const uint64_t input_buffer_size,
    buffered_input_file_t** const buffered_fasta_input,
    buffered_output_file_t** const buffered_output_file) {
  *buffered_fasta_input = buffered_input_file_new(parameters->input_file,input_buffer_size);
  *buffered_output_file = buffered_output_file_new(parameters->output_file);
  buffered_input_file_attach_buffered_output(*buffered_fasta_input,*buffered_output_file);
}
uint64_t mapper_pe_reload_buffers(
    mapper_parameters_t* const parameters,
    buffered_input_file_t* const buffered_fasta_input_end1,
    buffered_input_file_t* const buffered_fasta_input_end2,
    mapper_stats_t* const mapper_stats) {
  // Check end-of-block
  error_code_t error_code;
  if (buffered_input_file_eob(buffered_fasta_input_end1)) {
    // Reset template-length estimation
    search_paired_parameters_t* const search_paired_parameters =
        &parameters->search_parameters.search_paired_parameters;
    mapper_stats_template_init(mapper_stats,
        search_paired_parameters->template_length_estimation_samples,
        search_paired_parameters->template_length_estimation_min,
        search_paired_parameters->template_length_estimation_max);
    // Read new input-block
    if (!parameters->io.separated_input_files) {
      error_code = buffered_input_file_reload(buffered_fasta_input_end1,0);
      if (error_code==INPUT_STATUS_EOF) return INPUT_STATUS_EOF;
      const uint64_t buffered_lines = buffered_fasta_input_end1->num_lines;
      if ((buffered_lines % 8) != 0) {
        MAPPER_ERROR_PE_PARSE_UNSYNCH_INPUT_FILES(parameters);
      }
    } else {
      // Dump buffer (explicitly before synch-reload)
      buffered_output_file_dump_buffer(buffered_fasta_input_end1->attached_buffered_output_file);
      // Reload buffers (in synch)
      MUTEX_BEGIN_SECTION(parameters->input_file_mutex) {
        // Check synch
        if (!buffered_input_file_eob(buffered_fasta_input_end2)) {
          MAPPER_ERROR_PE_PARSE_UNSYNCH_INPUT_FILES(parameters);
        }
        // Reload end1
        error_code = buffered_input_file_reload(buffered_fasta_input_end1,0);
        // Reload end2
        const uint64_t num_lines_read = buffered_fasta_input_end1->num_lines;
        error_code = buffered_input_file_reload(buffered_fasta_input_end2,num_lines_read);
        if (buffered_fasta_input_end2->num_lines != num_lines_read) {
          MAPPER_ERROR_PE_PARSE_UNSYNCH_INPUT_FILES(parameters);
        }
        if (num_lines_read==0) {
          MUTEX_END_SECTION(parameters->input_file_mutex);
          return INPUT_STATUS_EOF;
        }
      } MUTEX_END_SECTION(parameters->input_file_mutex);
    }
  }
  // OK
  return INPUT_STATUS_OK;
}
void mapper_pe_prepare_io_buffers(
    const mapper_parameters_t* const parameters,
    const uint64_t input_buffer_size,
    buffered_input_file_t** const buffered_fasta_input_end1,
    buffered_input_file_t** const buffered_fasta_input_end2,
    buffered_output_file_t** const buffered_output_file) {
  if (parameters->io.separated_input_files) {
    *buffered_fasta_input_end1 = buffered_input_file_new(parameters->input_file_end1,input_buffer_size);
    *buffered_fasta_input_end2 = buffered_input_file_new(parameters->input_file_end2,input_buffer_size);
  } else {
    *buffered_fasta_input_end1 = buffered_input_file_new(parameters->input_file,input_buffer_size);
    *buffered_fasta_input_end2 = *buffered_fasta_input_end1;
  }
  *buffered_output_file = buffered_output_file_new(parameters->output_file);
  buffered_input_file_attach_buffered_output(*buffered_fasta_input_end1,*buffered_output_file);
}
error_code_t mapper_pe_parse_paired_sequences(
    const mapper_parameters_t* const parameters,
    buffered_input_file_t* const buffered_fasta_input_end1,
    buffered_input_file_t* const buffered_fasta_input_end2,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  error_code_t error_code;
  // Read end1
  error_code = input_fasta_parse_sequence(buffered_fasta_input_end1,
      archive_search_get_sequence(archive_search_end1),false);
  if (gem_expect_false(error_code!=INPUT_STATUS_OK)) return error_code;
  // Read end2
  error_code = input_fasta_parse_sequence(buffered_fasta_input_end2,
      archive_search_get_sequence(archive_search_end2),false);
  if (gem_expect_false(error_code!=INPUT_STATUS_OK)) return error_code;
  // OK
  PROF_ADD_COUNTER(GP_MAPPER_NUM_READS,2);
  return INPUT_STATUS_OK;
}
/*
 * Input (High-level)
 */
error_code_t mapper_se_read_single_sequence(
    archive_search_t* const archive_search,
    buffered_input_file_t* const buffered_fasta_input) {
  sequence_t* const sequence = archive_search_get_sequence(archive_search);
  const error_code_t error_code = input_fasta_parse_sequence(buffered_fasta_input,sequence,true);
  if (gem_expect_false(error_code==INPUT_STATUS_FAIL)) pthread_exit(0); // Abort
  PROF_INC_COUNTER(GP_MAPPER_NUM_READS);
  return error_code; // OK
}
error_code_t mapper_pe_read_paired_sequences(
    mapper_parameters_t* const mapper_parameters,
    buffered_input_file_t* const buffered_fasta_input_end1,
    buffered_input_file_t* const buffered_fasta_input_end2,
    archive_search_handlers_t* const archive_search_handlers,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  error_code_t error_code;
  // Check/Reload buffers manually
  error_code = mapper_pe_reload_buffers(
      mapper_parameters,buffered_fasta_input_end1,
      buffered_fasta_input_end2,archive_search_handlers->mapper_stats);
  if (error_code==INPUT_STATUS_EOF) return INPUT_STATUS_EOF;
  if (gem_expect_false(error_code==INPUT_STATUS_FAIL)) pthread_exit(0); // Abort
  // Read end1 & end2
  error_code = mapper_pe_parse_paired_sequences(mapper_parameters,
      buffered_fasta_input_end1,buffered_fasta_input_end2,
      archive_search_end1,archive_search_end2);
  if (gem_expect_false(error_code==INPUT_STATUS_FAIL)) pthread_exit(0); // Abort
  // OK
  return INPUT_STATUS_OK;
}
/*
 * Output
 */
void mapper_se_output_matches(
    mapper_parameters_t* const parameters,
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search,
    matches_t* const matches,
    mapping_stats_t *mstats) {

  if (archive_search->approximate_search.ns) {
    output_fastq(buffered_output_file,&archive_search->sequence);
  }

//#ifdef MAPPER_OUTPUT
//  switch (parameters->io.output_format) {
//    case MAP:
//      output_map_single_end_matches(buffered_output_file,archive_search,matches,&parameters->io.map_parameters);
//      break;
//    case SAM:
//      output_sam_single_end_matches(buffered_output_file,archive_search,matches,&parameters->io.sam_parameters);
//      break;
//    default:
//      GEM_INVALID_CASE();
//      break;
//  }
//  if (mstats) {
//    collect_se_mapping_stats(archive_search,matches,mstats);
//  }
//#endif
}
void mapper_pe_output_matches(
    mapper_parameters_t* const parameters,
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,
    mapping_stats_t* const mstats) {
#ifdef MAPPER_OUTPUT
  switch (parameters->io.output_format) {
    case MAP:
      output_map_paired_end_matches(buffered_output_file,archive_search_end1,
          archive_search_end2,paired_matches,&parameters->io.map_parameters);
      break;
    case SAM:
      output_sam_paired_end_matches(buffered_output_file,
          archive_search_end1,archive_search_end2,
          paired_matches,&parameters->io.sam_parameters);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  if (mstats) collect_pe_mapping_stats(archive_search_end1,archive_search_end2,paired_matches,mstats);
#endif
}
