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
 *   Input module allows fast parsing of FASTA/FASTQ buffered-input files
 */

#ifndef INPUT_FASTA_H_
#define INPUT_FASTA_H_

#include "system/commons.h"
#include "io/buffered_input_file.h"
#include "text/sequence.h"

/*
 * Input FASTA/FASTQ Constants
 */
#define FASTA_TAG_BEGIN '>'
#define FASTQ_TAG_BEGIN '@'
#define FASTQ_SEP '+'

/*
 * Read Parser
 */
int input_fasta_parse_sequence(
    buffered_input_file_t* const buffered_fasta_input,
    sequence_t* const sequence,
    const bool check_input_buffer);

/*
 * Display
 */
void input_fasta_parser_prompt_error(
    buffered_input_file_t* const buffered_input,
    const error_code_t error_code);

#endif /* INPUT_FASTA_H_ */
