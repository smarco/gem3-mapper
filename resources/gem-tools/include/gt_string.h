/*
 * PROJECT: GEM-Tools library
 * FILE: gt_string.h
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple string implementation.
 *   Static stings gt_string_new(0), which share memory across instances (stores mem ptr)
 *   Dynamic strings gt_string_new(n>0), which handle their own memory and hold copy of the string
 */

#ifndef GT_STRING_H_
#define GT_STRING_H_

#include "gt_commons.h"
#include "gt_error.h"

typedef struct {
  char* buffer;
  uint64_t allocated;
  uint64_t length;
} gt_string;

/*
 * Checkers
 */
#define GT_STRING_CHECK(string) GT_NULL_CHECK(string)
#define GT_STRING_CHECK_BUFFER(string) \
  GT_STRING_CHECK(string); \
  GT_NULL_CHECK(string->buffer)
#define GT_STRING_CHECK_NO_STATIC(string) \
  GT_STRING_CHECK(string); \
  gt_fatal_check(string->allocated==0,STRING_STATIC)

/*
 * Printers
 */
#define PRIgts "%.*s"
#define PRIgts_content(string) (int)gt_string_get_length(string),gt_string_get_string(string)
#define PRIgts_range_content(string,start_pos,length) (length),(gt_string_get_string(string)+start_pos)
#define PRIgts_trimmed_content(string,left_trim,right_trim) ((int)gt_string_get_length(string)-((left_trim)+(right_trim))),(gt_string_get_string(string)+(left_trim))

/*
 * Constructor & Accessors
 */
GT_INLINE gt_string* gt_string_new(const uint64_t initial_buffer_size);
GT_INLINE gt_string* gt_string_set_new(const char* const string_src);
GT_INLINE void gt_string_resize(gt_string* const string,const uint64_t new_buffer_size);
GT_INLINE void gt_string_clear(gt_string* const string);
GT_INLINE void gt_string_delete(gt_string* const string);

GT_INLINE bool gt_string_is_static(gt_string* const string);
GT_INLINE void gt_string_cast_static(gt_string* const string);
GT_INLINE void gt_string_cast_dynamic(gt_string* const string,const uint64_t initial_buffer_size);

GT_INLINE void gt_string_set_string(gt_string* const string,char* const string_src);
GT_INLINE void gt_string_set_nstring(gt_string* const string,char* const string_src,const uint64_t length);
GT_INLINE void gt_string_set_string_static(gt_string* const string,const char* const string_src);
GT_INLINE void gt_string_set_nstring_static(gt_string* const string,const char* const string_src,const uint64_t length);
GT_INLINE char* gt_string_get_string(gt_string* const string);

GT_INLINE uint64_t gt_string_get_length(gt_string* const string);
GT_INLINE void gt_string_set_length(gt_string* const string,const uint64_t length);

GT_INLINE char* gt_string_char_at(gt_string* const string,const uint64_t pos);

GT_INLINE void gt_string_append_char(gt_string* const string_dst,const char character);
GT_INLINE void gt_string_append_eos(gt_string* const string_dst);

#define gt_string_append_string    gt_string_right_append_string
#define gt_string_append_gt_string gt_string_right_append_gt_string
GT_INLINE void gt_string_left_append_string(gt_string* const string_dst,const char* const string_src,const uint64_t length);
GT_INLINE void gt_string_left_append_gt_string(gt_string* const string_dst,gt_string* const string_src);
GT_INLINE void gt_string_right_append_string(gt_string* const string_dst,const char* const string_src,const uint64_t length);
GT_INLINE void gt_string_right_append_gt_string(gt_string* const string_dst,gt_string* const string_src);

GT_INLINE void gt_string_trim_left(gt_string* const string,const uint64_t length);
GT_INLINE void gt_string_trim_right(gt_string* const string,const uint64_t length);

/*
 * Cmp functions
 */
GT_INLINE bool gt_string_is_null(gt_string* const string);
GT_INLINE int64_t gt_string_cmp(gt_string* const string_a,gt_string* const string_b);
GT_INLINE int64_t gt_string_ncmp(gt_string* const string_a,gt_string* const string_b,const uint64_t length);
GT_INLINE bool gt_string_equals(gt_string* const string_a,gt_string* const string_b);
GT_INLINE bool gt_string_nequals(gt_string* const string_a,gt_string* const string_b,const uint64_t length);

/*
 * Handlers
 */
GT_INLINE void gt_string_reverse(gt_string* const sequence);

GT_INLINE gt_string* gt_string_dup(gt_string* const sequence);
GT_INLINE void gt_string_copy(gt_string* const sequence_dst,gt_string* const sequence_src);
GT_INLINE void gt_string_reverse_copy(gt_string* const sequence_dst,gt_string* const sequence_src);
GT_INLINE gt_string* gt_string_reverse_dup(gt_string* const sequence);
GT_INLINE uint64_t gt_string_copy_substr(gt_string * const sequence_dst,gt_string * const sequence_src,uint64_t off,uint64_t len);

/*
 * String Printers
 */
GT_INLINE gt_status gt_vsprintf(gt_string* const sequence,const char *template,va_list v_args);
GT_INLINE gt_status gt_sprintf(gt_string* const sequence,const char *template,...);
GT_INLINE gt_status gt_vsprintf_append(gt_string* const sequence,const char *template,va_list v_args);
GT_INLINE gt_status gt_sprintf_append(gt_string* const sequence,const char *template,...);

/*
 * Iterator
 */
#define GT_STRING_ITERATE(string,mem,pos) \
  uint64_t pos; \
  const uint64_t __length_##mem = gt_string_get_length(string); \
  char* mem = gt_string_get_string(string); \
  for (pos=0;pos<__length_##mem;++pos) /* mem[pos] */

/*
 * String-Buffer functions
 */
GT_INLINE void gt_strncpy(char* const buffer_dst,const char* const buffer_src,const uint64_t length);
GT_INLINE char* gt_strndup(const char* const buffer,const uint64_t length);
GT_INLINE int gt_strcmp(const char* const buffer_a,const char* const buffer_b);
GT_INLINE bool gt_streq(const char* const buffer_a,const char* const buffer_b);
GT_INLINE int gt_strncmp(const char* const buffer_a,const char* const buffer_b,const uint64_t length);
GT_INLINE bool gt_strneq(const char* const buffer_a,const char* const buffer_b,const uint64_t length);
GT_INLINE uint64_t gt_strlen(const char* const buffer);

/*
 * Error Messages
 */
#define GT_ERROR_STRING_STATIC "Could not perform operation on static string"

#endif /* GT_STRING_H_ */
