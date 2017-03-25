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
 */

#include "text/dna_text.h"
#include "matches/align/match_alignment.h"
#include "matches/matches_cigar.h"

/*
 * Display
 */
void match_alignment_print_pretty(
    FILE* const stream,
    match_alignment_t* const match_alignment,
    vector_t* const cigar_vector,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length,
    mm_allocator_t* const mm_allocator) {
  mm_allocator_push_state(mm_allocator);
  tab_fprintf(stream,"[GEM]>Match.alignment\n");
  tab_global_inc();
  tab_fprintf(stream,"=> CIGAR  ");
  const uint64_t max_buffer_length = text_length+key_length+1;
  char* const key_alg = mm_allocator_calloc(mm_allocator,max_buffer_length,char,true);
  char* const ops_alg = mm_allocator_calloc(mm_allocator,max_buffer_length,char,true);
  char* const text_alg = mm_allocator_calloc(mm_allocator,max_buffer_length,char,true);
  cigar_element_t* cigar_element = vector_get_elm(cigar_vector,match_alignment->cigar_offset,cigar_element_t);
  uint64_t i, j, alg_pos = 0, read_pos = 0, text_pos = 0;
  for (i=0;i<match_alignment->cigar_length;++i,++cigar_element) {
    switch (cigar_element->type) {
      case cigar_match:
        fprintf(stream,"%d",(uint32_t)cigar_element->length);
        for (j=0;j<cigar_element->length;++j) {
          if (key[read_pos] != text[text_pos]) {
            key_alg[alg_pos] = dna_decode(key[read_pos]);
            ops_alg[alg_pos] = '*';
            text_alg[alg_pos++] = dna_decode(text[text_pos]);
          } else {
            key_alg[alg_pos] = dna_decode(key[read_pos]);
            ops_alg[alg_pos] = '|';
            text_alg[alg_pos++] = dna_decode(text[text_pos]);
          }
          read_pos++; text_pos++;
        }
        break;
      case cigar_mismatch:
        fprintf(stream,"%c",dna_decode(cigar_element->mismatch));
        if (key[read_pos] != text[text_pos]) {
          key_alg[alg_pos] = dna_decode(key[read_pos++]);
          ops_alg[alg_pos] = 'M';
          text_alg[alg_pos++] = dna_decode(text[text_pos++]);
        } else {
          key_alg[alg_pos] = dna_decode(key[read_pos++]);
          ops_alg[alg_pos] = '*';
          text_alg[alg_pos++] = dna_decode(text[text_pos++]);
        }
        break;
      case cigar_ins:
        fprintf(stream,">%u+",cigar_element->length);
        for (j=0;j<cigar_element->length;++j) {
          key_alg[alg_pos] = '-';
          ops_alg[alg_pos] = ' ';
          text_alg[alg_pos++] = dna_decode(text[text_pos++]);
        }
        break;
      case cigar_del:
        for (j=0;j<cigar_element->length;++j) {
          key_alg[alg_pos] = dna_decode(key[read_pos++]);
          ops_alg[alg_pos] = ' ';
          text_alg[alg_pos++] = '-';
        }
        if (cigar_element->attributes==cigar_attr_trim) {
          fprintf(stream,"(%u)",cigar_element->length);
        } else {
          fprintf(stream,">%u-",cigar_element->length);
        }
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
  key_alg[alg_pos] = '\0';
  ops_alg[alg_pos] = '\0';
  text_alg[alg_pos] = '\0';
  fprintf(stream,"\n");
  tab_fprintf(stream,"=> Pretty.Alignment\n");
  tab_fprintf(stream,"   KEY--%s--\n",key_alg);
  tab_fprintf(stream,"        %s  \n",ops_alg);
  tab_fprintf(stream,"   TXT--%s--\n",text_alg);
  tab_global_dec();
  mm_allocator_pop_state(mm_allocator);
}
