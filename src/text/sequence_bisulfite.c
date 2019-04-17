/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *  Copyright (c) 2011-2017 by Simon Heath  <simon.heath@gmail.com>
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
 *            Simon Heath <simon.heath@gmail.com>
 * DESCRIPTION:
 */

#include "text/sequence_bisulfite.h"

bool sequence_bisulfite_check_cg_depletion_se(
    sequence_t* const read_sequence) {
  uint64_t counts[7] = {0, 0, 0, 0, 0, 0, 0};
  const uint64_t length = sequence_get_length(read_sequence);
  char* const buffer = sequence_get_read(read_sequence);
  uint64_t pos;
  for (pos=0;pos<length;pos++) {
    counts[(int)dna_encode_table[(int)buffer[pos]]]++;
  }
  return counts[1] > counts[2];
}
bool sequence_bisulfite_check_cg_depletion_pe(
    sequence_t* const read_sequence1,
    sequence_t* const read_sequence2) {
  uint64_t counts1[7] = {0, 0, 0, 0, 0, 0, 0};
  uint64_t counts2[7] = {0, 0, 0, 0, 0, 0, 0};
  uint64_t length = sequence_get_length(read_sequence1);
  char* buffer = sequence_get_read(read_sequence1);
  uint64_t pos;
  for (pos=0;pos<length;pos++) {
    counts1[(int)dna_encode_table[(int)buffer[pos]]]++;
  }
  length = sequence_get_length(read_sequence2);
  buffer = sequence_get_read(read_sequence2);
  for (pos=0;pos<length;pos++) {
    counts2[(int)dna_encode_table[(int)buffer[pos]]]++;
  }
  return (counts1[1] + counts2[2]) > (counts1[2] + counts2[1]);
  }
