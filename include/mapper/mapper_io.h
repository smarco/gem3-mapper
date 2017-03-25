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

#ifndef MAPPER_IO_H_
#define MAPPER_IO_H_

#include "utils/essentials.h"
#include "mapper/mapper_parameters.h"

/*
 * Mapper I/O handler
 */
typedef struct {
  // Parameters
  bool paired_end;
  bool separated_input_files;
  uint64_t template_length_estimation_samples;
  uint64_t template_length_estimation_min;
  uint64_t template_length_estimation_max;
  mapper_parameters_io_t* mapper_parameters_io;
  // Bisulfite
  bool archive_bisulfite;
  bisulfite_read_t bisulfite_read;
  // Buffered Input
  buffered_input_file_t* buffered_fasta_input_end1;
  buffered_input_file_t* buffered_fasta_input_end2;
  // Buffered Output
  buffered_output_file_t* buffered_output_file;
  // Stats
  mapper_stats_t* mapper_stats;
  // Mutex
  pthread_mutex_t* input_mutex;
  // MM
  mm_allocator_t* mm_allocator;
} mapper_io_handler_t;

/*
 * Index loader
 */
void mapper_load_index(mapper_parameters_t* const parameters);

/*
 * Setup
 */
mapper_io_handler_t* mapper_io_handler_new_se(
    mapper_parameters_t* const parameters,
    const uint64_t input_buffer_size,
    mm_allocator_t* const mm_allocator);
mapper_io_handler_t* mapper_io_handler_new_pe(
    mapper_parameters_t* const parameters,
    const uint64_t input_buffer_size,
    mapper_stats_t* const mapper_stats,
    mm_allocator_t* const mm_allocator);
void mapper_io_handler_delete(
    mapper_io_handler_t* const mapper_io_handler);

/*
 * Reload/Parse Helpers
 */
uint64_t mapper_se_reload_buffer(
    mapper_io_handler_t* const mapper_io_handler);
error_code_t mapper_se_parse_sequence(
    mapper_io_handler_t* const mapper_io_handler,
    sequence_t** const sequence_end);

uint64_t mapper_pe_reload_buffer(
    mapper_io_handler_t* const mapper_io_handler);
error_code_t mapper_pe_parse_sequence(
    mapper_io_handler_t* const mapper_io_handler,
    sequence_t** const sequence_end1,
    sequence_t** const sequence_end2);

/*
 * Sequence readers
 */
error_code_t mapper_read_sequence(
    mapper_io_handler_t* const mapper_io_handler,
    const bool reload_input_buffer,
    sequence_t** const sequence);
error_code_t mapper_read_paired_sequence(
    mapper_io_handler_t* const mapper_io_handler,
    const bool reload_input_buffer,
    sequence_t** const sequence_end1,
    sequence_t** const sequence_end2);

/*
 * Output
 */
void mapper_io_handler_output_matches(
    mapper_io_handler_t* const mapper_io_handler,
    archive_search_t* const archive_search,
    matches_t* const matches,
    mapping_stats_t* const mstats);
void mapper_io_handler_output_paired_matches(
    mapper_io_handler_t* const mapper_io_handler,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,
    mapping_stats_t* const mstats);

#endif /* MAPPER_IO_H_ */
