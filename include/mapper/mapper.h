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

#ifndef MAPPER_H_
#define MAPPER_H_

#include "utils/essentials.h"
#include "archive/search/archive_search.h"
#include "archive/search/archive_search_handlers.h"
#include "mapper/mapper_parameters.h"
#include "mapper/mapper_io.h"
#include "matches/matches.h"
#include "matches/paired_matches.h"
#include "stats/report_stats.h"
#include "stats/report_stats_mstats.h"

/*
 * Mapper Search
 */
typedef struct {
  /* Thread Info */
  uint64_t thread_id;
  pthread_t* thread_data;
  /* Mapper parameters */
  mapper_parameters_t* mapper_parameters;
	/* Stats */
	mapping_stats_t* mapping_stats; // Per thread stats report structures
  /* Debug/Error structures */
  mapper_io_handler_t* mapper_io_handler;
  sequence_t** sequence_end1;
  sequence_t** sequence_end2;
  /* Ticker */
  ticker_t* ticker;
} mapper_search_t;

/*
 * Report
 */
void mapper_display_input_state(
    FILE* stream,
    buffered_input_file_t* const buffered_fasta_input,
    const sequence_t* const sequence);

/*
 * SE Mapper
 */
void mapper_se_run(mapper_parameters_t* const mapper_parameters);

/*
 * PE Mapper
 */
void mapper_pe_run(mapper_parameters_t* const mapper_parameters);

#endif /* MAPPER_H_ */
