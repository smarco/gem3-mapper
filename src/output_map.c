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

/*
 * Utils
 */
#define OUTPUT_MAP_COUNTER_SEPARATOR(counter_position,mcs) (counter_position==0)?'\t':(mcs==counter_position?'+':':')

/*
 * Output MAP
 */
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
        bofprintf_uint64(buffered_output_file,(uint32_t)cigar_array->match_length);
        break;
      case cigar_mismatch:
        bofprintf_char(buffered_output_file,dna_decode(cigar_array->mismatch));
        break;
      case cigar_ins:
        bofprintf_char(buffered_output_file,'>');
        bofprintf_int64(buffered_output_file,cigar_array->indel.indel_length);
        bofprintf_char(buffered_output_file,'+');
        break;
      case cigar_del:
        bofprintf_char(buffered_output_file,'>');
        bofprintf_int64(buffered_output_file,cigar_array->indel.indel_length);
        bofprintf_char(buffered_output_file,'-');
        break;
      case cigar_soft_trim:
        bofprintf_char(buffered_output_file,'(');
        bofprintf_int64(buffered_output_file,cigar_array->indel.indel_length);
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
    const uint64_t match_number,const match_trace_t* const match_trace,const bool print_mapq) {
  // Reserve
  const uint64_t sequence_length = gem_strlen(match_trace->sequence_name);
  buffered_output_file_reserve(buffered_output_file,sequence_length+INT_MAX_LENGTH+7);
  // Print Sequence Name
  if (gem_expect_true(match_number!=UINT64_MAX)) {
    bofprintf_char(buffered_output_file,(match_number>0) ? ',' : '\t');
  }
  bofprintf_string(buffered_output_file,sequence_length,match_trace->sequence_name);
  // Print Strand
  bofprintf_char(buffered_output_file,':');
  bofprintf_char(buffered_output_file,(match_trace->strand==Forward) ? '+' : '-');
  // Print Position
  bofprintf_char(buffered_output_file,':');
  bofprintf_uint64(buffered_output_file,match_trace->text_position+1); /* Base-1 */
  // Print CIGAR
  bofprintf_char(buffered_output_file,':');
  output_map_print_cigar(buffered_output_file,
      match_trace_get_cigar_array(matches,match_trace),
      match_trace_get_cigar_length(match_trace));
  if (print_mapq) {
    bofprintf_char(buffered_output_file,':');
    bofprintf_char(buffered_output_file,':');
    bofprintf_char(buffered_output_file,':');
    bofprintf_uint64(buffered_output_file,match_trace->mapq_score);
  }
}
GEM_INLINE void output_map_print_paired_match(
    buffered_output_file_t* const buffered_output_file,
    const matches_t* const matches_end1,const matches_t* const matches_end2,
    const uint64_t match_number,const paired_match_t* const paired_match) {
  output_map_print_match(buffered_output_file,matches_end1,match_number,paired_match->match_end1,false);
  buffered_output_file_reserve(buffered_output_file,2);
  bofprintf_char(buffered_output_file,':');
  bofprintf_char(buffered_output_file,':');
  output_map_print_match(buffered_output_file,matches_end2,UINT64_MAX,paired_match->match_end2,false);
  bofprintf_char(buffered_output_file,':');
  bofprintf_char(buffered_output_file,':');
  bofprintf_char(buffered_output_file,':');
  bofprintf_uint64(buffered_output_file,paired_match->mapq_score);
}
/*
 * MAP Alignment pretty
 */
GEM_INLINE void output_map_alignment_pretty(
    FILE* const stream,match_trace_t* const match_trace,matches_t* const matches,
    uint8_t* const key,const uint64_t key_length,uint8_t* const text,
    const uint64_t text_length,mm_stack_t* const mm_stack) {
  mm_stack_push_state(mm_stack);
  fprintf(stream,"%s:%lu:%c:",match_trace->sequence_name,
      match_trace->text_position,(match_trace->strand==Forward)?'+':'-');
  char* const key_alg = mm_stack_calloc(mm_stack,2*key_length,char,true);
  char* const ops_alg = mm_stack_calloc(mm_stack,2*key_length,char,true);
  char* const text_alg = mm_stack_calloc(mm_stack,2*key_length,char,true);
  cigar_element_t* cigar_element = vector_get_elm(matches->cigar_buffer,match_trace->cigar_buffer_offset,cigar_element_t);
  uint64_t i, j, alg_pos = 0, read_pos = 0, text_pos = 0;
  for (i=0;i<match_trace->cigar_length;++i,++cigar_element) {
    switch (cigar_element->type) {
      case cigar_match:
        fprintf(stream,"%d",(uint32_t)cigar_element->match_length);
        for (j=0;j<cigar_element->match_length;++j) {
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
        fprintf(stream,">%u+",cigar_element->indel.indel_length);
        for (j=0;j<cigar_element->indel.indel_length;++j) {
          key_alg[alg_pos] = '-';
          ops_alg[alg_pos] = ' ';
          text_alg[alg_pos++] = dna_decode(text[text_pos++]);
        }
        break;
      case cigar_del:
      case cigar_soft_trim:
        for (j=0;j<cigar_element->indel.indel_length;++j) {
          key_alg[alg_pos] = dna_decode(key[read_pos++]);
          ops_alg[alg_pos] = ' ';
          text_alg[alg_pos++] = '-';
        }
        if (cigar_element->type==cigar_del) {
          fprintf(stream,">%u-",cigar_element->indel.indel_length);
        } else {
          fprintf(stream,"(%u)",cigar_element->indel.indel_length);
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
  fprintf(stream,"\nKEY--%s--\n",key_alg);
  fprintf(stream,"     %s  \n",ops_alg);
  fprintf(stream,"TXT--%s--\n",text_alg);
  mm_stack_pop_state(mm_stack,false);
}
/*
 * MAP Tag
 */
GEM_INLINE void output_map_single_end_print_tag(
    buffered_output_file_t* const buffered_output_file,sequence_t* const seq_read) {
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
GEM_INLINE void output_map_paired_end_print_tag(
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read_end1,sequence_t* const seq_read_end2) {
  // TODO Redefine behavior in case of mismatching tags
  // Print Tag + End-Info
  const uint64_t tag_length = string_get_length(&seq_read_end1->tag);
  buffered_output_file_reserve(buffered_output_file,tag_length+3);
  bofprintf_string(buffered_output_file,PRIs_content(&seq_read_end1->tag));
  // Print CASAVA Tag (if present)
  if (sequence_has_casava_tag(seq_read_end1)) {
    const uint64_t casava_tag_length = string_get_length(&seq_read_end1->attributes.casava_tag);
    buffered_output_file_reserve(buffered_output_file,casava_tag_length);
    bofprintf_string(buffered_output_file,casava_tag_length,string_get_buffer(&seq_read_end1->attributes.casava_tag));
  }
  // Print Extra Tag (if present)
  if (sequence_has_extra_tag(seq_read_end1)) {
    const uint64_t extra_tag_length = string_get_length(&seq_read_end1->attributes.extra_tag);
    buffered_output_file_reserve(buffered_output_file,extra_tag_length);
    bofprintf_string(buffered_output_file,extra_tag_length,string_get_buffer(&seq_read_end1->attributes.extra_tag));
  }
}
/*
 * MAP Sequence Read/Qualities
 */
GEM_INLINE void output_map_single_end_print_read__qualities(
    buffered_output_file_t* const buffered_output_file,sequence_t* const seq_read) {
  // Select proper case
  const uint64_t read_length = string_get_length(&seq_read->read);
  if (sequence_has_qualities(seq_read)) {
    buffered_output_file_reserve(buffered_output_file,read_length+read_length+2);
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_string(buffered_output_file,read_length,string_get_buffer(&seq_read->read));
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_string(buffered_output_file,read_length,string_get_buffer(&seq_read->qualities));
  } else {
    buffered_output_file_reserve(buffered_output_file,read_length+1);
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_string(buffered_output_file,read_length,string_get_buffer(&seq_read->read));
  }
}
GEM_INLINE void output_map_paired_end_print_read__qualities(
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read_end1,sequence_t* const seq_read_end2) {
  // Select proper case
  const uint64_t read_length_end1 = string_get_length(&seq_read_end1->read);
  const uint64_t read_length_end2 = string_get_length(&seq_read_end2->read);
  if (sequence_has_qualities(seq_read_end1)) {
    buffered_output_file_reserve(buffered_output_file,2*read_length_end1+2*read_length_end2+4);
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_string(buffered_output_file,read_length_end1,string_get_buffer(&seq_read_end1->read));
    bofprintf_char(buffered_output_file,' ');
    bofprintf_string(buffered_output_file,read_length_end2,string_get_buffer(&seq_read_end2->read));
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_string(buffered_output_file,read_length_end1,string_get_buffer(&seq_read_end1->qualities));
    bofprintf_char(buffered_output_file,' ');
    bofprintf_string(buffered_output_file,read_length_end2,string_get_buffer(&seq_read_end2->qualities));
  } else {
    buffered_output_file_reserve(buffered_output_file,read_length_end1+read_length_end2+1);
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_string(buffered_output_file,read_length_end1,string_get_buffer(&seq_read_end1->read));
    bofprintf_char(buffered_output_file,' ');
    bofprintf_string(buffered_output_file,read_length_end2,string_get_buffer(&seq_read_end2->read));
  }
}
/*
 * MAP Counters
 */
GEM_INLINE void output_map_print_counters_account_mcs(
    buffered_output_file_t* const buffered_output_file,
    const uint64_t current_counter_pos,const uint64_t num_zeros,const bool compact) {
  if (compact && num_zeros>0) {
    bofprintf_char(buffered_output_file,(current_counter_pos==0)?'\t':':');
    bofprintf_string_literal(buffered_output_file,"0x");
    bofprintf_uint64(buffered_output_file,num_zeros);
    bofprintf_string_literal(buffered_output_file,"+0");
  } else {
    uint64_t i = 0;
    if (current_counter_pos==0) {
      bofprintf_string_literal(buffered_output_file,"\t0");
      i=1;
    }
    for (;i<num_zeros;++i) {
      bofprintf_string_literal(buffered_output_file,":0");
    }
    bofprintf_string_literal(buffered_output_file,"+0");
  }
}
GEM_INLINE void output_map_print_counters(
    buffered_output_file_t* const buffered_output_file,
    const vector_t* const counters_vector,const uint64_t mcs,const bool compact) {
  const uint64_t num_counters = vector_get_used(counters_vector);
  buffered_output_file_reserve(buffered_output_file,(num_counters+1)*(INT_MAX_LENGTH+1));
  // Zero counters
  if (gem_expect_false(num_counters==0)) {
    if (mcs==0 || mcs==ALL) {
      bofprintf_string_literal(buffered_output_file,"\t0");
    } else {
      output_map_print_counters_account_mcs(buffered_output_file,0,mcs,compact);
    }
    return;
  }
  // Print counters
  uint64_t i = 0;
  const uint64_t* counters = vector_get_mem(counters_vector,uint64_t);
  while (i < num_counters) {
    // Check zero counter
    if (compact && *counters==0) {
      // Check ahead number of zeros
      uint64_t j = 1;
      while (i+j<num_counters && i+j!=mcs && counters[i+j]==0) ++j;
      // Check if zero-counters must be compacted
      if (j >= OUTPUT_MAP_COUNTERS_MIN_COMPACT_ZEROS) {
        bofprintf_char(buffered_output_file,OUTPUT_MAP_COUNTER_SEPARATOR(i,mcs));
        bofprintf_string_literal(buffered_output_file,"0x");
        bofprintf_uint64(buffered_output_file,j);
      } else {
        // Print a series of zeros
        uint64_t k;
        for (k=0;k<j;++k) {
          bofprintf_char(buffered_output_file,OUTPUT_MAP_COUNTER_SEPARATOR(i+k,mcs));
          bofprintf_char(buffered_output_file,'0');
        }
      }
      i+=j; // Next (+j)
      counters+=j;
    } else {
      // Print Counter
      bofprintf_char(buffered_output_file,OUTPUT_MAP_COUNTER_SEPARATOR(i,mcs));
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
GEM_INLINE void output_map_single_end_matches(
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search,matches_t* const matches) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output_file);
  MATCHES_CHECK(matches);
  PROF_START_TIMER(GP_OUTPUT_MAP_SE);
  // Sort matches
  matches_sort_by_mapq_score(matches); // FIXME matches_sort_by_distance(matches);
  // Print TAG
  output_map_single_end_print_tag(buffered_output_file,&archive_search->sequence);
  // Print READ & QUALITIES
  output_map_single_end_print_read__qualities(buffered_output_file,&archive_search->sequence);
  // Print COUNTERS
  output_map_print_counters(buffered_output_file,matches->counters,matches->max_complete_stratum,false);
  // Print MATCHES (Traverse all matches (Position-matches))
  if (gem_expect_false(vector_get_used(matches->global_matches)==0)) {
    buffered_output_file_reserve(buffered_output_file,3);
    bofprintf_string_literal(buffered_output_file,"\t-\n");
  } else {
    VECTOR_ITERATE(matches->global_matches,match_trace,match_number,match_trace_t) {
      output_map_print_match(buffered_output_file,matches,match_number,match_trace,true); // Print match
    }
    // Next
    buffered_output_file_reserve(buffered_output_file,1);
    bofprintf_string_literal(buffered_output_file,"\n");
  }
  PROF_STOP_TIMER(GP_OUTPUT_MAP_SE);
}
/*
 * MAP PairedEnd
 */
GEM_INLINE void output_map_paired_end_matches(
    buffered_output_file_t* const buffered_output_file,archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,paired_matches_t* const paired_matches) {
  BUFFERED_OUTPUT_FILE_CHECK(buffered_output_file);
  PROF_START_TIMER(GP_OUTPUT_MAP_PE);
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  if (!paired_matches_is_mapped(paired_matches)) {
    output_map_single_end_matches(buffered_output_file,archive_search_end1,matches_end1);
    output_map_single_end_matches(buffered_output_file,archive_search_end2,matches_end2);
  } else {
    // Sort matches
    paired_matches_sort_by_distance(paired_matches);
    // Print TAG
    output_map_paired_end_print_tag(buffered_output_file,
        &archive_search_end1->sequence,&archive_search_end2->sequence);
    // Print READ & QUALITIES
    output_map_paired_end_print_read__qualities(buffered_output_file,
        &archive_search_end1->sequence,&archive_search_end2->sequence);
    // Print COUNTERS
    output_map_print_counters(buffered_output_file,
        paired_matches->concordant_counters,paired_matches->max_complete_stratum,false);
    if (gem_expect_false(vector_get_used(paired_matches->concordant_matches)==0)) {
      buffered_output_file_reserve(buffered_output_file,3);
      bofprintf_string_literal(buffered_output_file,"\t-\n");
    } else {
      // Print PAIRED MATCHES (Traverse all matches (Position-matches))
      VECTOR_ITERATE_CONST(paired_matches->concordant_matches,paired_match,n,paired_match_t) {
        output_map_print_paired_match(buffered_output_file,matches_end1,matches_end2,n,paired_match);
      }
    }
    buffered_output_file_reserve(buffered_output_file,1);
    bofprintf_string_literal(buffered_output_file,"\n");
  }
  PROF_STOP_TIMER(GP_OUTPUT_MAP_PE);
}
