/*
 * PROJECT: GEM-Tools library
 * FILE: gt_commons.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Base module containing general purpose functions
 */

#include "gt_commons.h"
#include "gt_error.h"

/*
 * Helper functions (OPERATIVE)
 */
GT_INLINE uint64_t gt_get_integer_proportion(const float proportion,const uint64_t length) {
  if (proportion<=0.0) return 0;
  if (proportion>=1.0) return (uint64_t)proportion;
  return (uint64_t)(proportion*(float)length);
}
/*
 * Print's template helpers
 */
GT_INLINE uint64_t gt_calculate_memory_required_v(const char *template,va_list v_args) {
  GT_NULL_CHECK(template);
  GT_NULL_CHECK(v_args);
  // Copy to avoid spoiling v_args
  va_list v_args_cpy;
  va_copy(v_args_cpy,v_args);
  // Calculate memory required to print the template{v_args}
  uint64_t mem_required = 0, precision=0;
  const char* centinel;
  for (centinel=template;*centinel!=EOS;++centinel,precision=0) {
    if (*centinel==FORMAT) {
      ++centinel;
      // Read modifiers
      while (gt_is_number(*centinel)) ++centinel;
      if (*centinel==DOT){
        ++centinel;
        if (*centinel==STAR) {
          ++centinel;
          precision = va_arg(v_args_cpy,int);
        } else {
          while (gt_is_number(*centinel)) ++centinel;
        }
      }
      gt_check(centinel==EOS,PRINT_FORMAT);
      // Check format
      switch (*centinel) {
        case 's': { // String requires fetching the argument length // FIXME: %.*s
          char* const string = va_arg(v_args_cpy,char*);
          mem_required += (precision>0) ? precision : strlen(string);
          break;
        }
        default:
          // As for the rest, we estimate the memory usage
          // Also we assume an upper bound over the possible formats (int, chars, floats, ...)
          va_arg(v_args_cpy,int);
          mem_required+=20;
          break;
      }
    } else {
      ++mem_required;
    }
  }
  return mem_required;
}
GT_INLINE uint64_t gt_calculate_memory_required_va(const char *template,...) {
  GT_NULL_CHECK(template);
  va_list v_args;
  va_start(v_args,template);
  return gt_calculate_memory_required_v(template,v_args);
}
/*
 * Random number generator
 */
GT_INLINE uint64_t gt_rand_IID(const uint64_t min,const uint64_t max) {
  int n_rand = rand(); // [0, RAND_MAX]
  const uint64_t range = max - min;
  const uint64_t rem = RAND_MAX % range;
  const uint64_t sample = RAND_MAX / range;
  // Consider the small interval within remainder of RAND_MAX
  if (n_rand < RAND_MAX - rem) {
    return min + n_rand/sample;
  } else {
    return gt_rand_IID(min,max);
  }
}
