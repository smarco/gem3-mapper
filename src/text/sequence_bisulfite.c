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

/*
 * Process Sequence
 */
void sequence_bisulfite_convert(
    string_t* const input_sequence,
    string_t* const backup_sequence,
    const char* const BS_table) {
  // Backup original sequence
  string_copy(backup_sequence,input_sequence);
  // Convert sequence
  const uint64_t length = string_get_length(input_sequence);
  char* const buffer = string_get_buffer(input_sequence);
  uint64_t pos;
  for (pos=0;pos<length;pos++) {
    buffer[pos] = BS_table[(int)buffer[pos]];
  }
}
bool sequence_bisulfite_check_cg_depletion_se(
    string_t* const read_sequence) {
  uint64_t counts[7] = {0, 0, 0, 0, 0, 0, 0};
  const uint64_t length = string_get_length(read_sequence);
  char* const buffer = string_get_buffer(read_sequence);
  uint64_t pos;
  for (pos=0;pos<length;pos++) {
    counts[(int)dna_encode_table[(int)buffer[pos]]]++;
  }
  return counts[1] > counts[2];
}
bool sequence_bisulfite_check_cg_depletion_pe(
    string_t* const read_sequence1,
    string_t* const read_sequence2) {
  uint64_t counts1[7] = {0, 0, 0, 0, 0, 0, 0};
  uint64_t counts2[7] = {0, 0, 0, 0, 0, 0, 0};
  uint64_t length = string_get_length(read_sequence1);
  char* buffer = string_get_buffer(read_sequence1);
  uint64_t pos;
  for (pos=0;pos<length;pos++) {
    counts1[(int)dna_encode_table[(int)buffer[pos]]]++;
  }
  length = string_get_length(read_sequence2);
  buffer = string_get_buffer(read_sequence2);
  for (pos=0;pos<length;pos++) {
    counts2[(int)dna_encode_table[(int)buffer[pos]]]++;
  }
  return (counts1[1] + counts2[2]) > (counts1[2] + counts2[1]);
}
void sequence_bisulfite_process_se(
    sequence_t* const sequence,
    const bisulfite_read_t bisulfite_read_mode) {
  // Fully convert reads before searching into archive, making a copy of the original
  if (bisulfite_read_mode==bisulfite_read_inferred) {
    sequence->bs_sequence_end = sequence->end_info;
  } else if (bisulfite_read_mode==bisulfite_non_stranded) {
    sequence->bs_sequence_end = sequence_bisulfite_check_cg_depletion_se(&sequence->read) ? paired_end2 : paired_end1;
  }
  string_clear(&sequence->bs_original_sequence);
  switch(sequence->bs_sequence_end) {
  case paired_end1:
    sequence_bisulfite_convert(&sequence->read,
        &sequence->bs_original_sequence,dna_bisulfite_C2T_table);
    break;
  case paired_end2:
    sequence_bisulfite_convert(&sequence->read,
        &sequence->bs_original_sequence,dna_bisulfite_G2A_table);
    break;
  default:
    break;
  }
  if (bisulfite_read_mode==bisulfite_read_interleaved) {
    sequence->bs_sequence_end = (sequence->bs_sequence_end==paired_end1) ? paired_end2 : paired_end1;
  }
}
void sequence_bisulfite_process_pe(
    sequence_t* const sequence_end1,
    sequence_t* const sequence_end2,
    const bisulfite_read_t bisulfite_read_mode) {
  // Fully convert reads before searching into archive, making a copy of the original
  bool flip = bisulfite_read_mode==bisulfite_non_stranded ?
    sequence_bisulfite_check_cg_depletion_pe(&sequence_end1->read, &sequence_end2->read) : false;
  if(flip) {
    sequence_bisulfite_convert(&sequence_end1->read,
      &sequence_end1->bs_original_sequence,dna_bisulfite_G2A_table);
    sequence_bisulfite_convert(&sequence_end2->read,
      &sequence_end2->bs_original_sequence,dna_bisulfite_C2T_table);
  } else {
    sequence_bisulfite_convert(&sequence_end1->read,
      &sequence_end1->bs_original_sequence,dna_bisulfite_C2T_table);
    sequence_bisulfite_convert(&sequence_end2->read,
      &sequence_end2->bs_original_sequence,dna_bisulfite_G2A_table);
  }
}
/*
 * Restore Sequence
 */
void sequence_bisulfite_restore_se(sequence_t* const sequence) {
  if (!string_is_null(&sequence->bs_original_sequence)) {
    string_copy(&sequence->read,&sequence->bs_original_sequence);
  }
}
void sequence_bisulfite_restore_pe(
    sequence_t* const sequence_end1,
    sequence_t* const sequence_end2) {
  string_copy(&sequence_end1->read,&sequence_end1->bs_original_sequence);
  string_copy(&sequence_end2->read,&sequence_end2->bs_original_sequence);
}
