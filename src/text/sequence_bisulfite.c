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
void sequence_bisulfite_process_se(
    sequence_t* const sequence,
    const bisulfite_read_t bisulfite_read_mode) {
  // Fully convert reads before searching into archive, making a copy of the original
  if (bisulfite_read_mode==bisulfite_read_inferred) {
    sequence->bs_sequence_end = sequence->end_info;
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
    sequence_t* const sequence_end2) {
  // Fully convert reads before searching into archive, making a copy of the original
  sequence_bisulfite_convert(&sequence_end1->read,
      &sequence_end1->bs_original_sequence,dna_bisulfite_C2T_table);
  sequence_bisulfite_convert(&sequence_end2->read,
      &sequence_end2->bs_original_sequence,dna_bisulfite_G2A_table);
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
