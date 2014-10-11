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
    buffered_output_file_t* const buffered_output_file,sequence_t* const seq_read) {
  // Print Tag + End-Info
  const uint64_t tag_length = string_get_length(&seq_read->tag);
  buffered_output_file_reserve(buffered_output_file,tag_length+3);
  switch (sequence_get_end_info(seq_read)) {
    case SINGLE_END:
      bofprintf_string(buffered_output_file,PRIs_content(&seq_read->tag));
      break;
    case PAIRED_END1:
      bofprintf_string(buffered_output_file,PRIs_content(&seq_read->tag));
      bofprintf_string_literal(buffered_output_file,"/1");
      break;
    case PAIRED_END2:
      bofprintf_string(buffered_output_file,PRIs_content(&seq_read->tag));
      bofprintf_string_literal(buffered_output_file,"/2");
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Print CASAVA Tag (if present)
  if (sequence_has_casava_tag(seq_read)) {
    const uint64_t casava_tag_length = string_get_length(&seq_read->attributes.casava_tag);
    buffered_output_file_reserve(buffered_output_file,casava_tag_length);
    bofprintf_string(buffered_output_file,casava_tag_length,string_get_buffer(&seq_read->attributes.casava_tag));
  }
  // Print Extra Tag (if present)
  if (sequence_has_extra_tag(seq_read)) {
    const uint64_t extra_tag_length = string_get_length(&seq_read->attributes.extra_tag);
    buffered_output_file_reserve(buffered_output_file,extra_tag_length);
    bofprintf_string(buffered_output_file,extra_tag_length,string_get_buffer(&seq_read->attributes.extra_tag));
  }
}
GEM_INLINE void output_map_print_read__qualities(
    buffered_output_file_t* const buffered_output_file,sequence_t* const seq_read) {
  // Select proper case
  const uint64_t read_length = string_get_length(&seq_read->read);
  if (sequence_has_qualities(seq_read)) {
    const uint64_t qualities_length = string_get_length(&seq_read->qualities);
    buffered_output_file_reserve(buffered_output_file,read_length+qualities_length+2);
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_string(buffered_output_file,read_length,string_get_buffer(&seq_read->read));
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_string(buffered_output_file,qualities_length,string_get_buffer(&seq_read->qualities));
  } else {
    buffered_output_file_reserve(buffered_output_file,read_length+1);
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_string(buffered_output_file,read_length,string_get_buffer(&seq_read->read));
  }
}
GEM_INLINE void output_map_print_counters(
    buffered_output_file_t* const buffered_output_file,matches_t* const matches,const bool compact) {
  // TODO Simple -> Use MCS to add zeros -> Or add them at search expanding counters !!!
  const uint64_t num_counters = vector_get_used(matches->counters);
  buffered_output_file_reserve(buffered_output_file,num_counters*(INT_MAX_LENGTH+1));
  // Zero counters
  if (gem_expect_false(num_counters==0)) {
    bofprintf_string_literal(buffered_output_file,"\t0");
    return;
  }
  // Print counters
  uint64_t i = 0;
  const uint64_t* counters = vector_get_mem(matches->counters,uint64_t);
  while (i < num_counters) {
    // Check zero counter
    if (compact && *counters==0) {
      // Check ahead number of zeros
      uint64_t j = 1;
      while (i+j < num_counters && counters[i+j]==0) ++j;
      // Check if zero-counters must be compacted
      if (j >= OUTPUT_MAP_COUNTERS_MIN_COMPACT_ZEROS) {
        bofprintf_char(buffered_output_file,(i > 0)?':':'\t');
        bofprintf_string_literal(buffered_output_file,"0x");
        bofprintf_uint64(buffered_output_file,j);
      } else {
        // Print a series of zeros
        uint64_t k;
        for (k=0;k<j;++k) {
          bofprintf_char(buffered_output_file,(k==0 && i==0)?'\t':':');
          bofprintf_char(buffered_output_file,'0');
        }
      }
      i+=j; // Next (+j)
      counters+=j;
    } else {
      // Print Counter
      bofprintf_char(buffered_output_file,(i > 0)?':':'\t');
      bofprintf_uint64(buffered_output_file,*counters);
      i++; // Next (+1)
      ++counters;
    }
  }
}
GEM_INLINE void output_map_print_cigar(
    buffered_output_file_t* const buffered_output_file,
    cigar_element_t* cigar_array,const uint64_t cigar_length) {
  // Reserve (upper-bound)
  buffered_output_file_reserve(buffered_output_file,cigar_length*(INT_MAX_LENGTH+2));
  // Traverse all CIGAR elements
  uint64_t i;
  for (i=0;i<cigar_length;++i,++cigar_array) {
    // Print CIGAR element
    switch (cigar_array->type) {
      case cigar_match:
        bofprintf_uint64(buffered_output_file,(uint32_t)cigar_array->length);
        break;
      case cigar_mismatch:
        bofprintf_char(buffered_output_file,dna_decode(cigar_array->mismatch));
        break;
      case cigar_ins:
        bofprintf_char(buffered_output_file,'>');
        bofprintf_int64(buffered_output_file,cigar_array->length);
        bofprintf_char(buffered_output_file,'+');
        break;
      case cigar_del:
        bofprintf_char(buffered_output_file,'>');
        bofprintf_int64(buffered_output_file,cigar_array->length);
        bofprintf_char(buffered_output_file,'-');
        break;
      case cigar_soft_trim:
        bofprintf_char(buffered_output_file,'(');
        bofprintf_int64(buffered_output_file,cigar_array->length);
        bofprintf_char(buffered_output_file,')');
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
GEM_INLINE void output_map_print_match(
    buffered_output_file_t* const buffered_output_file,const matches_t* const matches,
    const uint64_t match_number,const match_trace_t* const match_trace) {
  // Reserve
  const uint64_t sequence_length = gem_strlen(match_trace->sequence_name);
  buffered_output_file_reserve(buffered_output_file,sequence_length+INT_MAX_LENGTH+5);
  // Print Sequence Name
  bofprintf_char(buffered_output_file,(match_number > 0) ? ',' : '\t');
  bofprintf_string(buffered_output_file,sequence_length,match_trace->sequence_name);
  // Print Strand
  bofprintf_char(buffered_output_file,':');
  bofprintf_char(buffered_output_file,(match_trace->strand==Forward) ? '+' : '-');
  // Print Position
  bofprintf_char(buffered_output_file,':');
  bofprintf_uint64(buffered_output_file,match_trace->position+1); /* Base-1 */
  // Print CIGAR
  bofprintf_char(buffered_output_file,':');
  output_map_print_cigar(buffered_output_file,
      match_trace_get_cigar_array(matches,match_trace),
      match_trace_get_cigar_length(match_trace));
}
GEM_INLINE void output_map_single_end_matches(
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read,matches_t* const matches) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output_file);
  SEQUENCE_CHECK(seq_read);
  MATCHES_CHECK(matches);
  PROF_START_TIMER(GP_OUTPUT_MAP_SE);
  // Print TAG
  output_map_print_tag(buffered_output_file,seq_read);
  // Print READ & QUALITIES
  output_map_print_read__qualities(buffered_output_file,seq_read);
  // Print COUNTERS
  output_map_print_counters(buffered_output_file,matches,false);
  // Print MATCHES (Traverse all matches (Position-matches))
  if (gem_expect_false(vector_get_used(matches->global_matches)==0)) {
    buffered_output_file_reserve(buffered_output_file,3);
    bofprintf_string_literal(buffered_output_file,"\t-\n");
  } else {
    VECTOR_ITERATE(matches->global_matches,match_trace,match_number,match_trace_t) {
      // Print match
      output_map_print_match(buffered_output_file,matches,match_number,match_trace);
    }
    // Next
    buffered_output_file_reserve(buffered_output_file,1);
    bofprintf_string_literal(buffered_output_file,"\n");
  }
  PROF_STOP_TIMER(GP_OUTPUT_MAP_SE);
}
