/*
 * PROJECT: GEMMapper
 * FILE: string_buffer.c
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple string-buffer implementation.
 *   Static sting string_new(0), which share memory across instances (stores memory pointer)
 *   Dynamic string string_new(n>0), which handle their own memory and hold copy of the string
 */

#include "utils/string_buffer.h"
#include "system/errors.h"
#include "system/mm.h"

/*
 * Setup
 */
void string_init_(string_t* const string,char* const buffer,const uint64_t length,mm_stack_t* const mm_stack) {
  // Initialize
  if (buffer!=NULL) {
    // Static String
    string->buffer = buffer;
    string->allocated = 0;
    string->length = gem_strlen(buffer);
    string->mm_stack = NULL;
  } else if (mm_stack!=NULL) {
    // Dynamic-mmStack String
    string->buffer = mm_stack_malloc(mm_stack,length+1);
    string->buffer[0] = EOS;
    string->allocated = length+1;
    string->length = 0;
    string->mm_stack = mm_stack;
  } else {
    // Dynamic-heap String
    string->buffer = mm_malloc(length+1);
    string->buffer[0] = EOS;
    string->allocated = length+1;
    string->length = 0;
    string->mm_stack = NULL;
  }
}
void string_init(string_t* const string,const uint64_t length) {
  // Initialize Dynamic-heap String
  string_init_(string,NULL,length,NULL);
}
void string_init_static(string_t* const string,char* const buffer) {
  // Initialize Static String
  string_init_(string,buffer,0,NULL);
}
void string_init_mm(string_t* const string,const uint64_t length,mm_stack_t* const mm_stack) {
  // Initialize Dynamic-mmStack String
  string_init_(string,NULL,length,mm_stack);
}
void string_resize(string_t* const string,const uint64_t length) {
  const uint64_t new_buffer_size = length+1;
  if (string->allocated < new_buffer_size) {
    if (string->mm_stack!=NULL) {
      const char* const old_buffer = string->buffer;
      string->buffer = mm_stack_malloc(string->mm_stack,new_buffer_size);
      strncpy(string->buffer,old_buffer,string->length);
    } else {
      string->buffer = mm_realloc(string->buffer,new_buffer_size);
    }
    string->allocated = new_buffer_size;
  }
}
void string_clear(string_t* const string) {
  if (string->allocated) {
    // Dynamic String
    string->buffer[0] = EOS;
  } else {
    // Static String
    string->buffer = NULL;
  }
  string->length = 0;
}
void string_destroy(string_t* const string) {
  if (string->allocated && string->mm_stack==NULL) {
    mm_free(string->buffer);
  }
}
/*
 * Accessors
 */
void string_set_buffer(string_t* const string,char* const buffer,const uint64_t length) {
  if (gem_expect_true(string->allocated)) {
    // Dynamic String
    string_resize(string,length);
    gem_strncpy(string->buffer,buffer,length);
    string->length = length;
  } else {
    // Static String
    string->buffer = buffer;
    string->length = length;
  }
}
char* string_get_buffer(string_t* const string) {
  return string->buffer;
}
void string_set_length(string_t* const string,const uint64_t length) {
  string->length = length;
}
uint64_t string_get_length(string_t* const string) {
  return string->length;
}
char* string_char_at(string_t* const string,const uint64_t pos) {
  gem_fatal_check(pos>string->length,POSITION_OUT_OF_RANGE,pos,(uint64_t)0,string->length);
  return string->buffer+pos;
}
/*
 * Basic editing
 */
void string_append_char(string_t* const string,const char character) {
  string_resize(string,string->length);
  string->buffer[string->length] = character; // NOTE: No EOS appended
  ++string->length;
}
void string_append_eos(string_t* const string) {
  string_resize(string,string->length);
  string->buffer[string->length] = EOS;
}
/*
 * Append & trimming
 */
void string_left_append_buffer(string_t* const string,const char* const buffer,const uint64_t length) {
  const uint64_t base_length = string->length;
  const uint64_t final_length = base_length+length;
  string_resize(string,final_length);
  // Shift dst characters to the left
  char* const string_buffer = string->buffer;
  int64_t i;
  for (i=base_length;i>=0;--i) {
    string_buffer[i+length] = string_buffer[i];
  }
  // Left append the src string
  for (i=0;i<length;++i) {
    string_buffer[i] = buffer[i];
  }
  string->length = final_length;
}
void string_left_append_string(string_t* const string_dst,string_t* const string_src) {
  const uint64_t base_src_length = string_src->length;
  const uint64_t base_dst_length = string_dst->length;
  const uint64_t final_length = base_dst_length+base_src_length;
  string_resize(string_dst,final_length);
  // Shift dst characters to the left
  char* const buffer_src = string_src->buffer;
  char* const buffer_dst = string_dst->buffer;
  int64_t i;
  for (i=base_dst_length;i>=0;--i) {
    buffer_dst[i+base_src_length] = buffer_dst[i];
  }
  // Left append the src string
  for (i=0;i<base_src_length;++i) {
    buffer_dst[i] = buffer_src[i];
  }
  string_dst->length = final_length;
}
void string_right_append_buffer(string_t* const string,const char* const buffer,const uint64_t length) {
  const uint64_t final_length = string->length+length;
  string_resize(string,final_length);
  strncpy(string->buffer+string->length,buffer,length);
  string->buffer[final_length] = EOS;
  string->length = final_length;
}
void string_right_append_string(string_t* const string_dst,string_t* const string_src) {
  const uint64_t final_length = string_dst->length+string_src->length;
  string_resize(string_dst,final_length);
  strncpy(string_dst->buffer+string_dst->length,string_src->buffer,string_src->length);
  string_dst->buffer[final_length] = EOS;
  string_dst->length = final_length;
}
void string_trim_left(string_t* const string,const uint64_t length) {
  if (length > 0) {
    if (gem_expect_false(length >= string->length)) {
      string_clear(string);
    } else {
      const uint64_t new_length = string->length-length;
      uint64_t i;
      for (i=0;i<new_length;++i) string->buffer[i] = string->buffer[i+length];
      string->buffer[new_length] = EOS;
      string->length = new_length;
    }
  }
}
void string_trim_right(string_t* const string,const uint64_t length) {
  if (length > 0) {
    if (gem_expect_false(length >= string->length)) {
      string_clear(string);
    } else {
      string->length -= length;
      string->buffer[string->length] = EOS;
    }
  }
}
void string_copy_reverse(string_t* const string_dst,string_t* const string_src) {
  // Prepare Buffer
  const uint64_t length = string_get_length(string_src);
  string_resize(string_dst,length);
  string_clear(string_dst);
  // Reverse Read
  int64_t pos;
  const char* const buffer = string_get_buffer(string_src);
  for (pos=length-1;pos>=0;--pos) {
    string_append_char(string_dst,buffer[pos]);
  }
  string_append_eos(string_dst);
}
/*
 * Compare functions
 */
bool string_is_null(const string_t* const string) {
  if (gem_expect_false(string==NULL || string->length==0)) return true;
  return gem_expect_true(string->allocated > 0) ? string->buffer[0]==EOS : false;
}
int64_t string_cmp(string_t* const string_a,string_t* const string_b) {
  char* const buffer_a = string_a->buffer;
  char* const buffer_b = string_b->buffer;
  const uint64_t min_length = MIN(string_a->length,string_b->length);
  const uint64_t cmp = gem_strncmp(buffer_a,buffer_b,min_length);
  if (cmp) {
    return cmp;
  } else if (string_a->length == string_b->length) {
    return 0;
  } else {
    if (string_a->length < string_b->length) {
      return min_length+1;
    } else {
      return -min_length-1;
    }
  }
}
int64_t string_ncmp(string_t* const string_a,string_t* const string_b,const uint64_t length) {
  char* const buffer_a = string_a->buffer;
  char* const buffer_b = string_b->buffer;
  const uint64_t min_length = MIN(MIN(string_a->length,string_b->length),length);
  return gem_strncmp(buffer_a,buffer_b,min_length);
}
bool string_equals(string_t* const string_a,string_t* const string_b) {
  return string_cmp(string_a,string_b)==0;
}
bool string_nequals(string_t* const string_a,string_t* const string_b,const uint64_t length) {
  return string_ncmp(string_a,string_b,length)==0;
}
/*
 * Handlers
 */
string_t* string_dup(string_t* const string) {
  string_t* const string_cpy = mm_alloc(string_t);
  // Duplicate
  if (string->allocated==0) {
    // Static String
    string_init_(string_cpy,string->buffer,0,NULL);
  } else {
    // Dynamic String
    string_init_(string_cpy,NULL,string->length,string->mm_stack);
    strncpy(string_cpy->buffer,string->buffer,string->length);
  }
  string_cpy->length = string->length;
  return string_cpy;
}
void string_copy(string_t* const string_dst,string_t* const string_src) {
  string_resize(string_dst,string_src->length);
  strncpy(string_dst->buffer,string_src->buffer,string_src->length);
  string_dst->length = string_src->length;
}
/*
 * String Printers
 */
int sbprintf_v(string_t* const string,const char *template,va_list v_args) {
  int chars_printed;
  if (string->allocated>0) { // Allocate memory
    const uint64_t mem_required = calculate_memory_required_v(template,v_args);
    string_resize(string,mem_required);
  }
  chars_printed=vsprintf(string_get_buffer(string),template,v_args);
  string_set_length(string,chars_printed);
  string_get_buffer(string)[chars_printed] = EOS;
  return chars_printed;
}
int sbprintf(string_t* const string,const char *template,...) {
  va_list v_args;
  va_start(v_args,template);
  const int chars_printed = sbprintf_v(string,template,v_args);
  va_end(v_args);
  return chars_printed;
}
int sbprintf_append_v(string_t* const string,const char *template,va_list v_args) {
  int chars_printed = string_get_length(string);
  if (string->allocated>0) { // Allocate memory
    const uint64_t mem_required = calculate_memory_required_v(template,v_args);
    string_resize(string,mem_required+chars_printed);
  }
  chars_printed+=vsprintf(string_get_buffer(string)+chars_printed,template,v_args);
  string_set_length(string,chars_printed);
  string_get_buffer(string)[chars_printed] = EOS;
  return chars_printed;
}
int sbprintf_append(string_t* const string,const char *template,...) {
  va_list v_args;
  va_start(v_args,template);
  const int chars_printed = sbprintf_append_v(string,template,v_args);
  va_end(v_args);
  return chars_printed;
}

/*
 * String-Buffer functions
 */
void gem_strncpy(char* const buffer_dst,const char* const buffer_src,const uint64_t length) {
  memcpy(buffer_dst,buffer_src,length);
  buffer_dst[length] = EOS;
}
char* gem_strndup(const char* const buffer,const uint64_t length) {
  char* const buffer_cpy = mm_malloc(length+1);
  strncpy(buffer_cpy,buffer,length);
  return buffer_cpy;
}
char* gem_strdup(const char* const buffer) {
  return gem_strndup(buffer,strlen(buffer));
}
int gem_strcmp(const char* const buffer_a,const char* const buffer_b) {
  return strcmp(buffer_a,buffer_b);
}
int gem_strcasecmp(const char* const buffer_a,const char* const buffer_b) {
  return strcasecmp(buffer_a,buffer_b);
}
bool gem_streq(const char* const buffer_a,const char* const buffer_b) {
  return strcmp(buffer_a,buffer_b)==0;
}
bool gem_strcaseeq(const char* const buffer_a,const char* const buffer_b) {
  return strcasecmp(buffer_a,buffer_b)==0;
}
int gem_strncmp(const char* const buffer_a,const char* const buffer_b,const uint64_t length) {
  return strncmp(buffer_a,buffer_b,length);
}
int gem_strncasecmp(const char* const buffer_a,const char* const buffer_b,const uint64_t length) {
  return strncasecmp(buffer_a,buffer_b,length);
}
bool gem_strneq(const char* const buffer_a,const char* const buffer_b,const uint64_t length) {
  return strncmp(buffer_a,buffer_b,length)==0;
}
char* gem_strcat(const char* const buffer_a,const char* const buffer_b) {
  const uint64_t total_length = strlen(buffer_a) + strlen(buffer_b);
  char* const buffer_dst = mm_malloc(total_length+1);
  buffer_dst[0] = EOS;
  strcat(buffer_dst,buffer_a);
  strcat(buffer_dst,buffer_b);
  return buffer_dst;
}
uint64_t gem_strlen(const char* const buffer) {
  return strlen(buffer);
}
void gem_strrev(char* const buffer,const uint64_t length) {
  const uint64_t last_idx = length-1;
  const uint64_t middle = length/2;
  uint64_t i;
  for (i=0;i<middle;++i) buffer[i] = buffer[last_idx-i];
}
void gem_encrev(uint8_t* const buffer,const uint64_t length) {
  const uint64_t last_idx = length-1;
  const uint64_t middle = length/2;
  uint64_t i;
  for (i=0;i<middle;++i) buffer[i] = buffer[last_idx-i];
}
char* gem_strrmext(char* const buffer) {
  const int64_t total_length = strlen(buffer);
  int64_t i = total_length-1;
  while (i>=0) {
    if (buffer[i]==DOT) {buffer[i]=EOS; break;}
    --i;
  }
  return buffer;
}
char* gem_strbasename(char* const buffer) {
  const int64_t total_length = strlen(buffer);
  int64_t i = total_length-1;
  while (i>=0) {
    if (buffer[i]==SLASH) {
      return gem_strdup(buffer+i+1);
    }
    --i;
  }
  return gem_strdup(buffer);
}

