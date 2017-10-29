/*
 * PROJECT: GEM-Tools library
 * FILE: gt_dna_read.c
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple data structure to store genomic reads
 */

#include "gt_dna_read.h"
#include "gt_dna_string.h"

#define GT_DNA_READ_TAG_INITIAL_LENGTH 40
#define GT_DNA_READ_INITIAL_LENGTH 100

/*
 * Constructor
 */
GT_INLINE gt_dna_read* gt_dna_read_new(void) {
  gt_dna_read* read = gt_alloc(gt_dna_read);
  read->tag = gt_string_new(GT_DNA_READ_TAG_INITIAL_LENGTH);
  read->read = gt_string_new(GT_DNA_READ_INITIAL_LENGTH);
  read->qualities = gt_string_new(GT_DNA_READ_INITIAL_LENGTH);
  read->attributes = gt_attributes_new();
  return read;
}
GT_INLINE void gt_dna_read_clear(gt_dna_read* read) {
  GT_DNA_READ_CHECK(read);
  gt_string_clear(read->tag);
  gt_string_clear(read->read);
  gt_string_clear(read->qualities);
  gt_attributes_clear(read->attributes);
}
GT_INLINE void gt_dna_read_delete(gt_dna_read* read) {
  GT_DNA_READ_CHECK(read);
  gt_string_delete(read->tag);
  gt_string_delete(read->read);
  gt_string_delete(read->qualities);
  gt_attributes_delete(read->attributes);
}

/*
 * Accessors
 */
GT_INLINE void gt_dna_read_set_ntag(gt_dna_read* read,char* const text,const uint64_t length) {
  GT_DNA_READ_CHECK(read);
  gt_string_set_nstring(read->tag,text,length);
}
GT_INLINE void gt_dna_read_set_tag(gt_dna_read* read,char* const text) {
  GT_DNA_READ_CHECK(read);
  GT_NULL_CHECK(text);
  gt_string_set_string(read->tag,text);
}
GT_INLINE char* gt_dna_read_get_tag(gt_dna_read* const read) {
  GT_DNA_READ_CHECK(read);
  return gt_string_get_string(read->tag);
}

GT_INLINE void gt_dna_read_set_nread(gt_dna_read* read,char* const text,const uint64_t length) {
  GT_DNA_READ_CHECK(read);
  gt_string_set_nstring(read->read,text,length);
}
GT_INLINE void gt_dna_read_set_read(gt_dna_read* read,char* const text) {
  GT_DNA_READ_CHECK(read);
  GT_NULL_CHECK(text);
  gt_string_set_string(read->read,text);
}
GT_INLINE char* gt_dna_read_get_read(gt_dna_read* const read) {
  GT_DNA_READ_CHECK(read);
  return gt_string_get_string(read->read);
}

GT_INLINE void gt_dna_read_set_nqualities(gt_dna_read* read,char* const text,const uint64_t length) {
  GT_DNA_READ_CHECK(read);
  gt_string_set_nstring(read->qualities,text,length);
}
GT_INLINE void gt_dna_read_set_qualities(gt_dna_read* read,char* const text) {
  GT_DNA_READ_CHECK(read);
  GT_NULL_CHECK(text);
  gt_string_set_string(read->qualities,text);
}
GT_INLINE char* gt_dna_read_get_qualities(gt_dna_read* const read) {
  GT_DNA_READ_CHECK(read);
  return gt_string_get_string(read->qualities);
}

/*
 * Handlers/Utils
 */
GT_INLINE void gt_qualities_get_min__max(gt_string* const qualities,uint8_t* const min,uint8_t* const max) {
  GT_STRING_CHECK(qualities);
  *min = UINT8_MAX;
  *max = 0;
  GT_STRING_ITERATE(qualities,string,pos) {
    if (string[pos]<*min) *min = string[pos];
    if (string[pos]>*max) *max = string[pos];
  }
}
GT_INLINE gt_status gt_qualities_deduce_offset(gt_string* const qualities,gt_qualities_offset_t* qualities_offset_type) {
  GT_STRING_CHECK(qualities);
  GT_NULL_CHECK(qualities_offset_type);
  uint8_t min, max;
  gt_qualities_get_min__max(qualities,&min,&max);
  if (min >= 64) {*qualities_offset_type=GT_QUALS_OFFSET_64; return GT_STATUS_OK;}
  if (min >= 33) {*qualities_offset_type=GT_QUALS_OFFSET_33; return GT_STATUS_OK;}
  return GT_STATUS_FAIL;
}
GT_INLINE gt_status gt_qualities_adapt_from_offset33_to_offset64(gt_string* const qualities) {
  GT_STRING_CHECK(qualities);
  uint8_t min, max;
  gt_qualities_get_min__max(qualities,&min,&max);
  if (max > 255-64) return -1;
  GT_STRING_ITERATE(qualities,string,pos) {
    string[pos] += 33;
  }
  return 0;
}
GT_INLINE gt_status gt_qualities_adapt_from_offset64_to_offset33(gt_string* const qualities) {
  GT_STRING_CHECK(qualities);
  uint8_t min, max;
  gt_qualities_get_min__max(qualities,&min,&max);
  if (min < 64) return -1;
  GT_STRING_ITERATE(qualities,string,pos) {
    string[pos] -= 33;
  }
  return 0;
}
GT_INLINE gt_string* gt_qualities_dup__adapt_offset64_to_offset33(gt_string* const qualities) {
  gt_string* qualities_dst = gt_string_new(gt_string_get_length(qualities)+1);
  gt_string_set_length(qualities_dst,gt_string_get_length(qualities)+1);
  char* const qualities_dst_buf = gt_string_get_string(qualities_dst);
  GT_STRING_ITERATE(qualities,qualities_buf,pos) {
    qualities_dst_buf[pos] = qualities_buf[pos]-64+33;
  }
  return qualities_dst;
}
GT_INLINE void gt_dna_read_uniform_content_fix_length(
    gt_string* const read,gt_string* const qualities,
    const uint64_t length,const uint64_t uncall_begin,const uint64_t uncall_end) {
  // Set final length
  const uint64_t final_length = length-uncall_begin-uncall_end;
  gt_string_set_length(read,final_length);
  // Uniform qualities (if needed)
  if (qualities!=NULL) {
    char* const qualities_buffer = gt_string_get_string(qualities);
    int64_t i=0;
    for (i=0;i<final_length;++i) {
      qualities_buffer[i] = qualities_buffer[i+uncall_begin];
    }
    gt_string_set_length(qualities,final_length);
  }
}
GT_INLINE void gt_dna_read_uniform_content(gt_string* const read,gt_string* const qualities) {
  GT_STRING_CHECK(read);
  const uint64_t length = gt_string_get_length(read);
  char* const read_buffer = gt_string_get_string(read);
  int64_t i=0, uncall_begin=0, uncall_end=0;
  // Find all uncall bases at the beginning
  while (i<length && gt_get_dna_normalized(read_buffer[i])=='N') {
    ++i; ++uncall_begin;
  }
  // Convert the remaining bases
  while (i<length) {
    read_buffer[i-uncall_begin] = gt_get_dna_normalized(read_buffer[i]);
    ++i;
  }
  // Find all uncall bases at the end
  while (i>=0 && read_buffer[i]=='N') {
    --i; ++uncall_end;
  }
  // Fix length
  gt_dna_read_uniform_content_fix_length(read,qualities,length,uncall_begin,uncall_end);
}
GT_INLINE void gt_dna_read_uniform_strict_content(gt_string* const read,gt_string* const qualities) {
  GT_STRING_CHECK(read);
  const uint64_t length = gt_string_get_length(read);
  char* const read_buffer = gt_string_get_string(read);
  int64_t i=0, uncall_begin=0, uncall_end=0;
  // Find all uncall bases at the beginning
  while (i<length && gt_get_dna_strictly_normalized(read_buffer[i])=='N') {
    ++i; ++uncall_begin;
  }
  // Convert the remaining bases
  while (i<length) {
    read_buffer[i-uncall_begin] = gt_get_dna_strictly_normalized(read_buffer[i]);
    ++i;
  }
  // Find all uncall bases at the end
  while (i>=0 && read_buffer[i]=='N') {
    --i; ++uncall_end;
  }
  // Fix length
  gt_dna_read_uniform_content_fix_length(read,qualities,length,uncall_begin,uncall_end);
}
