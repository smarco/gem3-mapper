/*
 * PROJECT: GEMMapper
 * FILE: match_alignment_region.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "matches/align/match_alignment_region.h"
#include "matches/matches_cigar.h"

/*
 * Key/Text Operators
 */
uint64_t match_alignment_region_text_coverage(
    match_alignment_region_t* const match_alignment_region) {
  return match_alignment_region->text_end - match_alignment_region->text_begin;
}
uint64_t match_alignment_region_text_distance(
    match_alignment_region_t* const match_alignment_region_a,
    match_alignment_region_t* const match_alignment_region_b) {
  return match_alignment_region_b->text_begin - match_alignment_region_a->text_end;
}
bool match_alignment_region_text_overlap(
    match_alignment_region_t* const match_alignment_region_a,
    match_alignment_region_t* const match_alignment_region_b) {
  return !(match_alignment_region_a->text_end <= match_alignment_region_b->text_begin ||
           match_alignment_region_b->text_end <= match_alignment_region_a->text_begin);
}
/*
 * Compare
 */
int match_alignment_region_key_cmp(
    match_alignment_region_t* const match_alignment_region_a,
    match_alignment_region_t* const match_alignment_region_b) {
  return (int)match_alignment_region_a->key_begin - (int)match_alignment_region_b->key_begin;
}
int match_alignment_region_cmp_text_position(
    const match_alignment_region_t* const a,
    const match_alignment_region_t* const b) {
  return a->text_begin - b->text_begin;
}
/*
 * Display
 */
void match_alignment_region_print(
    FILE* const stream,
    match_alignment_region_t* const match_alignment_region,
    const uint64_t match_alignment_region_id,
    matches_t* const matches) {
  // Print alignment-region
  switch (match_alignment_region->region_type) {
    case match_alignment_region_exact:
      tab_fprintf(stream,"    %"PRIu64"[exact]\t",match_alignment_region_id);
      break;
    case match_alignment_region_approximate:
      tab_fprintf(stream,"    %"PRIu64"[approximate]\t",match_alignment_region_id);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  tab_fprintf(stream,"-> [%"PRIu64",%"PRIu64") ~> [+%"PRIu64",+%"PRIu64")",
      match_alignment_region->key_begin,match_alignment_region->key_end,
      match_alignment_region->text_begin,match_alignment_region->text_end);
  // Print CIGAR
  if (matches!=NULL &&
      match_alignment_region->region_type==match_alignment_region_approximate &&
      match_alignment_region->cigar_length>0) {
    tab_fprintf(stream,"\tCIGAR=");
    match_cigar_print(stream,matches->cigar_vector,
        match_alignment_region->cigar_buffer_offset,match_alignment_region->cigar_length);
  }
  tab_fprintf(stream,"\n");
}
void match_alignment_region_print_pretty(
    FILE* const stream,
    match_alignment_region_t* const match_alignment_region,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text) {
  // Compute offset
  uint64_t i, offset_text = 0, offset_key = 0;
  if (match_alignment_region->key_begin > match_alignment_region->text_begin) {
    offset_key = match_alignment_region->key_begin-match_alignment_region->text_begin;
  } else {
    offset_text = match_alignment_region->text_begin-match_alignment_region->key_begin;
    for (i=0;i<offset_text;++i) fprintf(stream," ");
  }
  // Print Key
  for (i=offset_key;i<key_length;++i) {
    if (match_alignment_region->key_begin <= i && i < match_alignment_region->key_end) {
      if (text[offset_text+i]==key[i]) {
        fprintf(stream,"%c",dna_decode(key[i]));
      } else {
        fprintf(stream,"%c",tolower(dna_decode(key[i])));
      }
    } else {
      fprintf(stream,"-");
    }
  }
  fprintf(stream,"\n");
}
