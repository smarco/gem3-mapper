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
 *   Output module provides MAP output generation
 */

#include "text/sequence.h"
#include "io/output_map.h"

/*
 * Constants
 */
#define OUTPUT_MAP_COUNTERS_MIN_COMPACT_ZEROS 5
//#define OUTPUT_MAP_SWG_SCORE_SUCCINT

/*
 * Utils
 */
#define OUTPUT_MAP_COUNTER_SEPARATOR(counter_position,mcs) \
  if (counter_position>0) { bofprintf_char(buffered_output_file,(mcs==counter_position?'+':':')); }

/*
 * Setup
 */
void output_map_parameters_set_defaults(output_map_parameters_t* const output_map_parameters) {
  output_map_parameters->format_version = map_format_v2;
}
/*
 * Output MAP
 */
void output_map_print_cigar_mapv2(
    buffered_output_file_t* const buffered_output_file,
    const cigar_element_t* cigar_array,
    const uint64_t cigar_length) {
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
        if (cigar_array->attributes == cigar_attr_trim) {
          bofprintf_char(buffered_output_file,'(');
          bofprintf_int64(buffered_output_file,cigar_array->length);
          bofprintf_char(buffered_output_file,')');
        } else {
          bofprintf_char(buffered_output_file,'>');
          bofprintf_int64(buffered_output_file,cigar_array->length);
          bofprintf_char(buffered_output_file,'-');
        }
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
void output_map_print_cigar_mapv3(
    buffered_output_file_t* const buffered_output_file,
    const cigar_element_t* cigar_array,
    const uint64_t cigar_length) {
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
        if (cigar_array->attributes == cigar_attr_homopolymer) {
          bofprintf_char(buffered_output_file,'&');
          bofprintf_int64(buffered_output_file,cigar_array->length);
          bofprintf_char(buffered_output_file,'+');
        } else {
          bofprintf_char(buffered_output_file,'>');
          bofprintf_int64(buffered_output_file,cigar_array->length);
          bofprintf_char(buffered_output_file,'+');
        }
        break;
      case cigar_del:
        if (cigar_array->attributes == cigar_attr_trim) {
          bofprintf_char(buffered_output_file,'(');
          bofprintf_int64(buffered_output_file,cigar_array->length);
          bofprintf_char(buffered_output_file,')');
        } else if (cigar_array->attributes == cigar_attr_homopolymer) {
          bofprintf_char(buffered_output_file,'&');
          bofprintf_int64(buffered_output_file,cigar_array->length);
          bofprintf_char(buffered_output_file,'-');
        } else {
          bofprintf_char(buffered_output_file,'>');
          bofprintf_int64(buffered_output_file,cigar_array->length);
          bofprintf_char(buffered_output_file,'-');
        }
        break;
      default:
        GEM_INVALID_CASE();
        break;
    }
  }
}
void output_map_print_match(
    buffered_output_file_t* const buffered_output_file,
    const matches_t* const matches,
    const match_trace_t* const match_trace,
    const bool print_mapq,
    const output_map_format_t output_map_format) {
  // Reserve
  const uint64_t sequence_length = gem_strlen(match_trace->sequence_name);
  buffered_output_file_reserve(buffered_output_file,sequence_length+INT_MAX_LENGTH+7);
  // Print Sequence Name
  bofprintf_string(buffered_output_file,sequence_length,match_trace->sequence_name);
  // Print Strand
  bofprintf_char(buffered_output_file,':');
  bofprintf_char(buffered_output_file,(match_trace->strand==Forward) ? '+' : '-');
  // Print Position
  bofprintf_char(buffered_output_file,':');
  bofprintf_uint64(buffered_output_file,match_trace->text_position+1); /* Base-1 */
#ifdef OUTPUT_MAP_SWG_SCORE_SUCCINT
  // SWG Score
  bofprintf_char(buffered_output_file,':');
  bofprintf_int64(buffered_output_file,match_trace->swg_score);
  bofprintf_char(buffered_output_file,':');
  bofprintf_uint64(buffered_output_file,match_trace->match_alignment.effective_length);
#else
  // Print CIGAR
  bofprintf_char(buffered_output_file,':');
  const cigar_element_t* const cigar_buffer = match_trace_get_cigar_buffer(match_trace,matches->cigar_vector);
  const uint64_t cigar_length = match_trace_get_cigar_length(match_trace);
  switch (output_map_format) {
    case map_format_v2:
      output_map_print_cigar_mapv2(buffered_output_file,cigar_buffer,cigar_length);
      break;
    case map_format_v3:
      output_map_print_cigar_mapv3(buffered_output_file,cigar_buffer,cigar_length);
      break;
    default:
      GEM_NOT_IMPLEMENTED();
      break;
  }
#endif
  if (print_mapq) {
    bofprintf_char(buffered_output_file,':');
    bofprintf_char(buffered_output_file,':');
    bofprintf_char(buffered_output_file,':');
    bofprintf_uint64(buffered_output_file,match_trace->mapq_score);
  }
}
void output_map_print_paired_match(
    buffered_output_file_t* const buffered_output_file,
    const matches_t* const matches_end1,
    const matches_t* const matches_end2,
    const paired_map_t* const paired_map,
    const output_map_format_t output_map_format) {
  // Map end/1
  output_map_print_match(buffered_output_file,matches_end1,paired_map->match_trace_end1,false,output_map_format);
  buffered_output_file_reserve(buffered_output_file,2);
  // Paired-end Separator
  bofprintf_char(buffered_output_file,':');
  bofprintf_char(buffered_output_file,':');
  // Map end/2
  output_map_print_match(buffered_output_file,matches_end2,paired_map->match_trace_end2,false,output_map_format);
  // MAPQ Score Separator
  bofprintf_char(buffered_output_file,':');
  bofprintf_char(buffered_output_file,':');
  bofprintf_char(buffered_output_file,':');
  // Print Paired-end MAPQ Score
  bofprintf_uint64(buffered_output_file,paired_map->mapq_score);
}
/*
 * MAP Alignment pretty
 */
void output_map_alignment_pretty(
    FILE* const stream,
    match_trace_t* const match_trace,
    matches_t* const matches,
    uint8_t* const key,
    const uint64_t key_length,
    uint8_t* const text,
    const uint64_t text_length,
    mm_allocator_t* const mm_allocator) {
  tab_fprintf(stream,"%s:%"PRIu64":%c\n",match_trace->sequence_name,
      match_trace->text_position,(match_trace->strand==Forward)?'+':'-');
  match_alignment_print_pretty(
      stream,&match_trace->match_alignment,
      key,key_length,text,text_length,
      matches->cigar_vector,mm_allocator);
}
/*
 * Separator
 */
void output_map_print_separator(
    buffered_output_file_t* const buffered_output_file,
    const char separator) {
  buffered_output_file_reserve(buffered_output_file,1);
  bofprintf_char(buffered_output_file,separator);
}
/*
 * MAP Tag
 */
void output_map_single_end_print_tag(
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read) {
  // Print Tag + End-Info
  const uint64_t tag_length = string_get_length(&seq_read->tag);
  buffered_output_file_reserve(buffered_output_file,tag_length+3);
  switch (sequence_get_end_info(seq_read)) {
    case single_end:
      bofprintf_string(buffered_output_file,PRIs_content(&seq_read->tag));
      break;
    case paired_end1:
      bofprintf_string(buffered_output_file,PRIs_content(&seq_read->tag));
      bofprintf_string_literal(buffered_output_file,"/1");
      break;
    case paired_end2:
      bofprintf_string(buffered_output_file,PRIs_content(&seq_read->tag));
      bofprintf_string_literal(buffered_output_file,"/2");
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Print CASAVA Tag (if present)
  if (sequence_has_casava_tag(seq_read)) {
    const uint64_t casava_tag_length = string_get_length(&seq_read->casava_tag);
    buffered_output_file_reserve(buffered_output_file,casava_tag_length);
    bofprintf_string(buffered_output_file,casava_tag_length,string_get_buffer(&seq_read->casava_tag));
  }
  // Print Extra Tag (if present)
  if (sequence_has_extra_tag(seq_read)) {
    const uint64_t extra_tag_length = string_get_length(&seq_read->extra_tag);
    buffered_output_file_reserve(buffered_output_file,extra_tag_length);
    bofprintf_string(buffered_output_file,extra_tag_length,string_get_buffer(&seq_read->extra_tag));
  }
}
void output_map_paired_end_print_tag(
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read_end1) {
  // TODO Redefine behavior in case of mismatching tags
  // Print Tag + End-Info
  const uint64_t tag_length = string_get_length(&seq_read_end1->tag);
  buffered_output_file_reserve(buffered_output_file,tag_length+3);
  bofprintf_string(buffered_output_file,PRIs_content(&seq_read_end1->tag));
  // Print CASAVA Tag (if present)
  if (sequence_has_casava_tag(seq_read_end1)) {
    const uint64_t casava_tag_length = string_get_length(&seq_read_end1->casava_tag);
    buffered_output_file_reserve(buffered_output_file,casava_tag_length);
    bofprintf_string(buffered_output_file,casava_tag_length,string_get_buffer(&seq_read_end1->casava_tag));
  }
  // Print Extra Tag (if present)
  if (sequence_has_extra_tag(seq_read_end1)) {
    const uint64_t extra_tag_length = string_get_length(&seq_read_end1->extra_tag);
    buffered_output_file_reserve(buffered_output_file,extra_tag_length);
    bofprintf_string(buffered_output_file,extra_tag_length,string_get_buffer(&seq_read_end1->extra_tag));
  }
}
/*
 * MAP Sequence Read/Qualities
 */
void output_map_single_end_print_read__qualities(
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read) {
  // Select proper case
  const uint64_t read_length = string_get_length(&seq_read->read);
  if (sequence_has_qualities(seq_read)) {
    buffered_output_file_reserve(buffered_output_file,read_length+read_length+1);
    bofprintf_string(buffered_output_file,read_length,string_get_buffer(&seq_read->read));
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_string(buffered_output_file,read_length,string_get_buffer(&seq_read->qualities));
  } else {
    buffered_output_file_reserve(buffered_output_file,read_length);
    bofprintf_string(buffered_output_file,read_length,string_get_buffer(&seq_read->read));
  }
}
void output_map_paired_end_print_read__qualities(
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read_end1,
    sequence_t* const seq_read_end2) {
  // Select proper case
  const uint64_t read_length_end1 = string_get_length(&seq_read_end1->read);
  const uint64_t read_length_end2 = string_get_length(&seq_read_end2->read);
  if (sequence_has_qualities(seq_read_end1)) {
    buffered_output_file_reserve(buffered_output_file,2*read_length_end1+2*read_length_end2+3);
    bofprintf_string(buffered_output_file,read_length_end1,string_get_buffer(&seq_read_end1->read));
    bofprintf_char(buffered_output_file,' ');
    bofprintf_string(buffered_output_file,read_length_end2,string_get_buffer(&seq_read_end2->read));
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_string(buffered_output_file,read_length_end1,string_get_buffer(&seq_read_end1->qualities));
    bofprintf_char(buffered_output_file,' ');
    bofprintf_string(buffered_output_file,read_length_end2,string_get_buffer(&seq_read_end2->qualities));
  } else {
    buffered_output_file_reserve(buffered_output_file,read_length_end1+read_length_end2+1);
    bofprintf_string(buffered_output_file,read_length_end1,string_get_buffer(&seq_read_end1->read));
    bofprintf_char(buffered_output_file,' ');
    bofprintf_string(buffered_output_file,read_length_end2,string_get_buffer(&seq_read_end2->read));
  }
}
/*
 * MAP Counters
 */
void output_map_print_counters_account_mcs(
    buffered_output_file_t* const buffered_output_file,
    const uint64_t current_counter_pos,
    const uint64_t num_zeros,
    const bool compact) {
  if (compact && num_zeros>0) {
    if (current_counter_pos>0) bofprintf_char(buffered_output_file,':');
    bofprintf_string_literal(buffered_output_file,"0x");
    bofprintf_uint64(buffered_output_file,num_zeros);
    bofprintf_string_literal(buffered_output_file,"+0");
  } else {
    uint64_t i = 0;
    if (current_counter_pos==0) {
      bofprintf_string_literal(buffered_output_file,"0");
      i=1;
    }
    for (;i<num_zeros;++i) {
      bofprintf_string_literal(buffered_output_file,":0");
    }
    bofprintf_string_literal(buffered_output_file,"+0");
  }
}
void output_map_print_counters(
    buffered_output_file_t* const buffered_output_file,
    matches_counters_t* const matches_counter,
    const uint64_t mcs,
    const bool compact) {
  const uint64_t num_counters = matches_counters_get_num_counters(matches_counter);
  const uint64_t max_length = MAX(num_counters,mcs);
  buffered_output_file_reserve(buffered_output_file,(max_length+1)*(INT_MAX_LENGTH+1));
  // Zero counters
  if (gem_expect_false(num_counters==0)) {
    if (mcs==0 || mcs==ALL) {
      bofprintf_string_literal(buffered_output_file,"0");
    } else {
      output_map_print_counters_account_mcs(buffered_output_file,0,mcs,compact);
    }
    return;
  }
  // Print counters
  uint64_t i = 0;
  const uint64_t* counters = matches_counters_get_counts(matches_counter);
  while (i < num_counters) {
    // Check zero counter
    if (compact && *counters==0) {
      // Check ahead number of zeros
      uint64_t j = 1;
      while (i+j<num_counters && i+j!=mcs && counters[i+j]==0) ++j;
      // Check if zero-counters must be compacted
      if (j >= OUTPUT_MAP_COUNTERS_MIN_COMPACT_ZEROS) {
        OUTPUT_MAP_COUNTER_SEPARATOR(i,mcs);
        bofprintf_string_literal(buffered_output_file,"0x");
        bofprintf_uint64(buffered_output_file,j);
      } else {
        // Print a series of zeros
        uint64_t k;
        for (k=0;k<j;++k) {
          OUTPUT_MAP_COUNTER_SEPARATOR(i+k,mcs);
          bofprintf_char(buffered_output_file,'0');
        }
      }
      i+=j; // Next (+j)
      counters+=j;
    } else {
      // Print Counter
      OUTPUT_MAP_COUNTER_SEPARATOR(i,mcs);
      bofprintf_uint64(buffered_output_file,*counters);
      i++; // Next (+1)
      ++counters;
    }
  }
  // Account for MCS
  if (i<=mcs) {
    output_map_print_counters_account_mcs(buffered_output_file,i,mcs-i,compact);
  }
}
/*
 * MAP SingleEnd
 */
void output_map_single_end_matches(
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search,
    matches_t* const matches,
    output_map_parameters_t* const output_map_parameters) {
  PROF_START_TIMER(GP_OUTPUT_MAP_SE);
  // Print TAG
  output_map_single_end_print_tag(buffered_output_file,archive_search->sequence);
  // Print READ & QUALITIES
  output_map_print_separator(buffered_output_file,'\t'); // Separator
  output_map_single_end_print_read__qualities(buffered_output_file,archive_search->sequence);
  // Print COUNTERS
  output_map_print_separator(buffered_output_file,'\t'); // Separator
  output_map_print_counters(buffered_output_file,matches->counters,matches->max_complete_stratum,false);
  // Print MATCHES (Traverse all matches (Position-matches))
  if (gem_expect_false(matches_get_num_match_traces(matches)==0)) {
    buffered_output_file_reserve(buffered_output_file,3);
    bofprintf_string_literal(buffered_output_file,"\t-\n");
  } else {
    const uint64_t num_match_traces = matches_get_num_match_traces(matches);
    match_trace_t** const match_traces = matches_get_match_traces(matches);
    uint64_t i;
    for (i=0;i<num_match_traces;++i) {
      // Separator
      if (i==0) output_map_print_separator(buffered_output_file,'\t');
      else output_map_print_separator(buffered_output_file,',');
      // Print Match
      output_map_print_match(buffered_output_file,matches,match_traces[i],true,output_map_parameters->format_version);
    }
    // Next
    output_map_print_separator(buffered_output_file,'\n'); // Separator
  }
  PROF_STOP_TIMER(GP_OUTPUT_MAP_SE);
}
/*
 * MAP PairedEnd
 */
void output_map_paired_end_matches(
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,
    output_map_parameters_t* const output_map_parameters) {
  PROF_START_TIMER(GP_OUTPUT_MAP_PE);
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  if (!paired_matches_is_mapped(paired_matches)) {
    output_map_single_end_matches(buffered_output_file,archive_search_end1,matches_end1,output_map_parameters);
    output_map_single_end_matches(buffered_output_file,archive_search_end2,matches_end2,output_map_parameters);
  } else {
    // Print TAG
    output_map_paired_end_print_tag(buffered_output_file,archive_search_end1->sequence);
    // Print READ & QUALITIES
    output_map_print_separator(buffered_output_file,'\t'); // Separator
    output_map_paired_end_print_read__qualities(
        buffered_output_file,
        archive_search_end1->sequence,
        archive_search_end2->sequence);
    // Print COUNTERS
    output_map_print_separator(buffered_output_file,'\t'); // Separator
    output_map_print_counters(
        buffered_output_file,paired_matches->counters,
        paired_matches_get_max_complete_stratum(paired_matches),false);
    if (paired_matches_get_num_maps(paired_matches)==0) {
      buffered_output_file_reserve(buffered_output_file,3);
      bofprintf_string_literal(buffered_output_file,"\t-\n");
    } else {
      // Print PAIRED MATCHES (Traverse all matches (Position-matches))
      const uint64_t num_paired_map = paired_matches_get_num_maps(paired_matches);
      paired_map_t** const paired_map = paired_matches_get_maps(paired_matches);
      uint64_t i;
      for (i=0;i<num_paired_map;++i) {
        // Separator
        output_map_print_separator(buffered_output_file,(i==0) ? '\t' : ',');
        // Paired-Map
        output_map_print_paired_match(buffered_output_file,matches_end1,
            matches_end2,paired_map[i],output_map_parameters->format_version);
      }
    }
    output_map_print_separator(buffered_output_file,'\n');
  }
  PROF_STOP_TIMER(GP_OUTPUT_MAP_PE);
}
/*
 * FASTA/FASTQ
 */
void output_fastq(
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const sequence) {
  // Print TAG
  buffered_output_file_reserve(buffered_output_file,1);
  output_map_print_separator(buffered_output_file,sequence->has_qualities?'@':'>');
  output_map_single_end_print_tag(buffered_output_file,sequence);
  buffered_output_file_reserve(buffered_output_file,1);
  output_map_print_separator(buffered_output_file,'\n');
  // Print Sequence
  const uint64_t seq_length = sequence_get_length(sequence);
  if (sequence->has_qualities) {
    buffered_output_file_reserve(buffered_output_file,2*seq_length+10);
    bofprintf_string(buffered_output_file,seq_length,string_get_buffer(&sequence->read));
    bofprintf_string(buffered_output_file,3,"\n+\n");
    bofprintf_string(buffered_output_file,seq_length,string_get_buffer(&sequence->qualities));
  } else {
    buffered_output_file_reserve(buffered_output_file,seq_length+10);
    bofprintf_string(buffered_output_file,seq_length,string_get_buffer(&sequence->read));
  }
  output_map_print_separator(buffered_output_file,'\n');
}
