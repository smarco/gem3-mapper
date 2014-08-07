/*
 * PROJECT: GEMMapper
 * FILE: output_map.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "output_map.h"

/*
 * Constants
 */
#define OUTPUT_MAP_COUNTERS_MIN_COMPACT_ZEROS 5

GEM_INLINE void output_map_print_tag(
    buffered_output_file_t* const buffered_output_file,const sequence_t* const seq_read) {
  // Print Tag + End-Info
  const uint64_t tag_length = string_get_length(seq_read->tag);
  switch (sequence_get_end_info(seq_read)) {
    case SINGLE_END:
      bofprintf_fixed(buffered_output_file,tag_length+1,"%"PRIs,PRIs_content(seq_read->tag));
      break;
    case PAIRED_END1:
      bofprintf_fixed(buffered_output_file,tag_length+3,"%"PRIs"/1",PRIs_content(seq_read->tag));
      break;
    case PAIRED_END2:
      bofprintf_fixed(buffered_output_file,tag_length+3,"%"PRIs"/2",PRIs_content(seq_read->tag));
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Print CASAVA Tag (if present)
  if (sequence_has_casava_tag(seq_read)) {
    const uint64_t casava_tag_length = string_get_length(seq_read->attributes.casava_tag);
    bofprintf_fixed(buffered_output_file,casava_tag_length+1,
        " %"PRIs,casava_tag_length,string_get_buffer(seq_read->attributes.casava_tag));
  }
  // Print Extra Tag (if present)
  if (sequence_has_extra_tag(seq_read)) {
    const uint64_t extra_tag_length = string_get_length(seq_read->attributes.extra_tag);
    bofprintf_fixed(buffered_output_file,extra_tag_length+1,
        " %"PRIs,extra_tag_length,string_get_buffer(seq_read->attributes.extra_tag));
  }
}
GEM_INLINE void output_map_print_read__qualities(
    buffered_output_file_t* const buffered_output_file,const sequence_t* const seq_read) {
  // Select proper case
  const uint64_t read_length = string_get_length(seq_read->read);
  if (sequence_has_qualities(seq_read)) {
    const uint64_t qualities_length = string_get_length(seq_read->qualities);
    bofprintf_fixed(buffered_output_file,read_length+qualities_length+2,"\t%"PRIs"\t%"PRIs,
        read_length,string_get_buffer(seq_read->read),
        qualities_length,string_get_buffer(seq_read->qualities));
  } else {
    bofprintf_fixed(buffered_output_file,read_length+1,"\t%"PRIs,
        read_length,string_get_buffer(seq_read->read));
  }
}
GEM_INLINE void output_map_print_counters(
    buffered_output_file_t* const buffered_output_file,matches_t* const matches) {
  const uint64_t* counters = vector_get_mem(matches->counters,uint64_t);
  const uint64_t num_counters = vector_get_used(matches->counters);
  uint64_t i = 0;
  while (i < num_counters) {
    // Check zero counter
    if (*counters==0) {
      // Check ahead number of zeros
      uint64_t j = 1;
      while (i+j < num_counters && counters[j]==0) ++j;
      // Check if zero-counters must be compacted
      if (j >= OUTPUT_MAP_COUNTERS_MIN_COMPACT_ZEROS) {
        bofprintf_fixed(buffered_output_file,INT_MAX_LENGTH+3,"%c0x%lu",(i > 0)?':':'\t',j);
      } else {
        // Print a series of zeros
        uint64_t k;
        for (k=0;k<j;++k) {
          bofprintf_fixed(buffered_output_file,INT_MAX_LENGTH+1,"%c0",(k==0 && i==0)?'\t':':');
        }
      }
      i+=j; // Next (+j)
    } else {
      // Print Counter
      bofprintf_fixed(buffered_output_file,INT_MAX_LENGTH+1,"%c%lu",(i > 0)?':':'\t',*counters);
      i++; // Next (+1)
    }
  }
}
GEM_INLINE void output_map_print_cigar(
    buffered_output_file_t* const buffered_output_file,
    cigar_element_t* cigar_array,const uint64_t cigar_length) {
  // Traverse all CIGAR elements
  uint64_t i;
  for (i=0;i<cigar_length;++i,++cigar_array) {
    // Print CIGAR element
    switch (cigar_array->type) {
      case cigar_match:
        bofprintf_fixed(buffered_output_file,INT_MAX_LENGTH,"%lu",(uint32_t)cigar_array->length);
        break;
      case cigar_mismatch:
        bofprintf_fixed(buffered_output_file,1,"%c",(char)cigar_array->length);
        break;
      case cigar_ins:
        bofprintf_fixed(buffered_output_file,INT_MAX_LENGTH+2,">%lu+",cigar_array->length);
        break;
      case cigar_del:
        bofprintf_fixed(buffered_output_file,INT_MAX_LENGTH+2,">%lu+",cigar_array->length);
        break;
      case cigar_soft_trim:
        bofprintf_fixed(buffered_output_file,INT_MAX_LENGTH+2,"(%lu)",cigar_array->length);
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
GEM_INLINE void output_map_print_match(
    buffered_output_file_t* const buffered_output_file,const matches_t* const matches,
    const uint64_t match_position,const match_trace_t* const match_trace) {
  // Print Sequence Name + Strand + Position
  bofprintf_fixed(buffered_output_file,gem_strlen(match_trace->sequence_name)+INT_MAX_LENGTH+5,
      "%c%s:%c:%lu:",
      (match_position > 0) ? ',' : '\t',
      match_trace->sequence_name,
      (match_trace->strand==Forward) ? '+' : '-',
      match_trace->position);
  // Print CIGAR
  output_map_print_cigar(buffered_output_file,
      match_trace_get_cigar_array(matches,match_trace),
      match_trace_get_cigar_length(match_trace));
}
GEM_INLINE void output_map_single_end_matches(
    buffered_output_file_t* const buffered_output_file,
    const sequence_t* const seq_read,matches_t* const matches) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output_file);
  SEQUENCE_CHECK(seq_read);
  MATCHES_CHECK(matches);
  // Print TAG
  output_map_print_tag(buffered_output_file,seq_read);
  // Print READ & QUALITIES
  output_map_print_read__qualities(buffered_output_file,seq_read);
  // Print COUNTERS
  output_map_print_counters(buffered_output_file,matches);
  // Print MATCHES (Traverse all matches (Position-matches))
  VECTOR_ITERATE(matches->global_matches,match_trace,match_position,match_trace_t) {
    // Print match
    output_map_print_match(buffered_output_file,matches,match_position,match_trace);
  }
  // Next
  bofprintf_fixed(buffered_output_file,1,"\n");
}
