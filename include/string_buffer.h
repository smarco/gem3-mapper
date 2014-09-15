/*
 * PROJECT: GEMMapper
 * FILE: string_buffer.h
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple string-buffer implementation.
 *   Static sting string_new(0), which share memory across instances (stores memory pointer)
 *   Dynamic string string_new(n>0), which handle their own memory and hold copy of the string
 */

#ifndef STRING_BUFFER_H_
#define STRING_BUFFER_H_

#include "commons.h"
#include "errors.h"
#include "mm_stack.h"

typedef struct {
  char* buffer;         // String buffer
  uint64_t allocated;   // Number of bytes allocated
  uint64_t length;      // Length of the string (not including EOS)
  /* MM */
  mm_stack_t* mm_stack; // MM-Stack
} string_t;

// Direction (traversal)
typedef enum { traversal_forward, traversal_backward } traversal_direction_t;

/*
 * Checkers
 */
#define STRING_CHECK(string) \
  GEM_CHECK_NULL((string)); \
  GEM_CHECK_NULL((string)->buffer)
#define STRING_DYNAMIC_CHECK(string) \
  STRING_CHECK((string)); \
  gem_fatal_check(!(string)->allocated,STRING_STATIC)

/*
 * Printers
 */
#define PRIs ".*s"
#define PRIs_content(string) (int)string_get_length((string_t*)string),string_get_buffer((string_t*)string)

/*
 * Constructor & Accessors
 */
GEM_INLINE void string_init(string_t* const string,const uint64_t length);
GEM_INLINE void string_init_static(string_t* const string,char* const buffer);
GEM_INLINE void string_init_mm(string_t* const string,const uint64_t length,mm_stack_t* const mm_stack);
GEM_INLINE void string_resize(string_t* const string,const uint64_t length);
GEM_INLINE void string_clear(string_t* const string);
GEM_INLINE void string_destroy(string_t* const string);

GEM_INLINE char* string_get_buffer(string_t* const string);
GEM_INLINE void string_set_buffer_const(string_t* const string,const char* const buffer,const uint64_t length);
GEM_INLINE void string_set_buffer(string_t* const string,char* const buffer_src,const uint64_t length);

GEM_INLINE char* string_char_at(string_t* const string,const uint64_t pos);
GEM_INLINE uint64_t string_get_length(string_t* const string);
GEM_INLINE void string_set_length(string_t* const string,const uint64_t length);

/*
 * Basic editing
 */
GEM_INLINE void string_append_char(string_t* const string,const char character);
GEM_INLINE void string_append_eos(string_t* const string);

/*
 * Append & trimming
 */
#define string_append_buffer string_right_append_buffer
#define string_append_string string_right_append_string
GEM_INLINE void string_left_append_buffer(string_t* const string,const char* const buffer,const uint64_t length);
GEM_INLINE void string_left_append_string(string_t* const string_dst,string_t* const string_src);
GEM_INLINE void string_right_append_buffer(string_t* const string,const char* const buffer,const uint64_t length);
GEM_INLINE void string_right_append_string(string_t* const string_dst,string_t* const string_src);
GEM_INLINE void string_trim_left(string_t* const string,const uint64_t length);
GEM_INLINE void string_trim_right(string_t* const string,const uint64_t length);

GEM_INLINE void string_copy_reverse(string_t* const string_dst,string_t* const string_src);

/*
 * Compare functions
 */
GEM_INLINE bool string_is_null(const string_t* const string);
GEM_INLINE int64_t string_cmp(string_t* const string_a,string_t* const string_b);
GEM_INLINE int64_t string_ncmp(string_t* const string_a,string_t* const string_b,const uint64_t length);
GEM_INLINE bool string_equals(string_t* const string_a,string_t* const string_b);
GEM_INLINE bool string_nequals(string_t* const string_a,string_t* const string_b,const uint64_t length);

/*
 * Handlers
 */
GEM_INLINE string_t* string_dup(string_t* const string);
GEM_INLINE void string_copy(string_t* const string_dst,string_t* const string_src);

/*
 * String Printers
 */
GEM_INLINE int sbprintf_v(string_t* const string,const char *template,va_list v_args);
GEM_INLINE int sbprintf(string_t* const string,const char *template,...);
GEM_INLINE int sbprintf_append_v(string_t* const string,const char *template,va_list v_args);
GEM_INLINE int sbprintf_append(string_t* const string,const char *template,...);

/*
 * Iterator
 */
#define STRING_ITERATE(string,mem,pos) \
  uint64_t pos; \
  const uint64_t __length_##mem = string_get_length(string); \
  char* mem = string_get_buffer(string); \
  for (pos=0;pos<__length_##mem;++pos) /* mem[pos] */

/*
 * Basic String Functions Wrappers
 */
GEM_INLINE void gem_strncpy(char* const buffer_dst,const char* const buffer_src,const uint64_t length);
GEM_INLINE char* gem_strndup(const char* const buffer,const uint64_t length);
GEM_INLINE char* gem_strdup(const char* const buffer);
GEM_INLINE int gem_strcmp(const char* const buffer_a,const char* const buffer_b);
GEM_INLINE int gem_strcasecmp(const char* const buffer_a,const char* const buffer_b);
GEM_INLINE bool gem_streq(const char* const buffer_a,const char* const buffer_b);
GEM_INLINE bool gem_strcaseeq(const char* const buffer_a,const char* const buffer_b);
GEM_INLINE int gem_strncmp(const char* const buffer_a,const char* const buffer_b,const uint64_t length);
GEM_INLINE int gem_strncasecmp(const char* const buffer_a,const char* const buffer_b,const uint64_t length);
GEM_INLINE bool gem_strneq(const char* const buffer_a,const char* const buffer_b,const uint64_t length);
GEM_INLINE char* gem_strcat(const char* const buffer_a,const char* const buffer_b);
GEM_INLINE uint64_t gem_strlen(const char* const buffer);

GEM_INLINE char* gem_strrmext(char* const buffer);
GEM_INLINE char* gem_strbasename(char* const buffer);


/*
 * Error Messages
 */
#define GEM_ERROR_STRING_STATIC "Could not perform operation on static string"

#endif /* STRING_BUFFER_H_ */
