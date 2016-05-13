/*
 * PROJECT: GEMMapper
 * FILE: mapper_bisulfite.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper/mapper_bisulfite.h"
#include "archive/archive_search_se.h"
#include "archive/archive_search_pe.h"

#define BISULFITE_SEQUENCE_INITIAL_LENGTH 200

void mapper_bisulfite_convert_sequence(
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
void mapper_bisulfite_process_sequence_se(
    archive_search_t* const archive_search,
    search_parameters_t* const search_parameters) {
  // Fully convert reads before searching into archive, making a copy of the original
  const bisulfite_read_t bs_read_mode = search_parameters->bisulfite_read;
  if (bs_read_mode==bisulfite_read_inferred) {
    archive_search->bs_sequence_end = sequence_get_end_info(&archive_search->sequence);
  }
  string_clear(archive_search->bs_original_sequence);
  switch(archive_search->bs_sequence_end) {
  case paired_end1:
    mapper_bisulfite_convert_sequence(&archive_search->sequence.read,
        archive_search->bs_original_sequence,dna_bisulfite_C2T_table);
    break;
  case paired_end2:
    mapper_bisulfite_convert_sequence(&archive_search->sequence.read,
        archive_search->bs_original_sequence,dna_bisulfite_G2A_table);
    break;
  default:
    break;
  }
  if (bs_read_mode==bisulfite_read_interleaved) {
    archive_search->bs_sequence_end = (archive_search->bs_sequence_end==paired_end1) ? paired_end2 : paired_end1;
  }
}
void mapper_bisulfite_process_sequence_pe(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  // Fully convert reads before searching into archive, making a copy of the original
  mapper_bisulfite_convert_sequence(&archive_search_end1->sequence.read,
      archive_search_end1->bs_original_sequence,dna_bisulfite_C2T_table);
  mapper_bisulfite_convert_sequence(&archive_search_end2->sequence.read,
      archive_search_end2->bs_original_sequence,dna_bisulfite_G2A_table);
}
void mapper_bisulfite_restore_sequence_se(
    archive_search_t* const archive_search) {
  if (!string_is_null(archive_search->bs_original_sequence)) {
    string_copy(&archive_search->sequence.read,archive_search->bs_original_sequence);
  }
}
void mapper_bisulfite_restore_sequence_pe(
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2) {
  string_copy(&archive_search_end1->sequence.read,archive_search_end1->bs_original_sequence);
  string_copy(&archive_search_end2->sequence.read,archive_search_end2->bs_original_sequence);
}
