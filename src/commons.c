/*
 * PROJECT: GEMMapper
 * FILE: commons.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Base module containing general purpose functions
 */

#include "commons.h"

/*
 * Common masks
 */
const uint64_t uint64_mask_ones[] = {
  0x0000000000000000ull,
  0x0000000000000001ull, 0x0000000000000003ull, 0x0000000000000007ull, 0x000000000000000Full,
  0x000000000000001Full, 0x000000000000003Full, 0x000000000000007Full, 0x00000000000000FFull,
  0x00000000000001FFull, 0x00000000000003FFull, 0x00000000000007FFull, 0x0000000000000FFFull,
  0x0000000000001FFFull, 0x0000000000003FFFull, 0x0000000000007FFFull, 0x000000000000FFFFull,
  0x000000000001FFFFull, 0x000000000003FFFFull, 0x000000000007FFFFull, 0x00000000000FFFFFull,
  0x00000000001FFFFFull, 0x00000000003FFFFFull, 0x00000000007FFFFFull, 0x0000000000FFFFFFull,
  0x0000000001FFFFFFull, 0x0000000003FFFFFFull, 0x0000000007FFFFFFull, 0x000000000FFFFFFFull,
  0x000000001FFFFFFFull, 0x000000003FFFFFFFull, 0x000000007FFFFFFFull, 0x00000000FFFFFFFFull,
  0x00000001FFFFFFFFull, 0x00000003FFFFFFFFull, 0x00000007FFFFFFFFull, 0x0000000FFFFFFFFFull,
  0x0000001FFFFFFFFFull, 0x0000003FFFFFFFFFull, 0x0000007FFFFFFFFFull, 0x000000FFFFFFFFFFull,
  0x000001FFFFFFFFFFull, 0x000003FFFFFFFFFFull, 0x000007FFFFFFFFFFull, 0x00000FFFFFFFFFFFull,
  0x00001FFFFFFFFFFFull, 0x00003FFFFFFFFFFFull, 0x00007FFFFFFFFFFFull, 0x0000FFFFFFFFFFFFull,
  0x0001FFFFFFFFFFFFull, 0x0003FFFFFFFFFFFFull, 0x0007FFFFFFFFFFFFull, 0x000FFFFFFFFFFFFFull,
  0x001FFFFFFFFFFFFFull, 0x003FFFFFFFFFFFFFull, 0x007FFFFFFFFFFFFFull, 0x00FFFFFFFFFFFFFFull,
  0x01FFFFFFFFFFFFFFull, 0x03FFFFFFFFFFFFFFull, 0x07FFFFFFFFFFFFFFull, 0x0FFFFFFFFFFFFFFFull,
  0x1FFFFFFFFFFFFFFFull, 0x3FFFFFFFFFFFFFFFull, 0x7FFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull
};
const uint64_t uint64_mask_inverse_ones[] =
{
  0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFEull, 0xFFFFFFFFFFFFFFFCull, 0xFFFFFFFFFFFFFFF8ull,
  0xFFFFFFFFFFFFFFF0ull, 0xFFFFFFFFFFFFFFE0ull, 0xFFFFFFFFFFFFFFC0ull, 0xFFFFFFFFFFFFFF80ull,
  0xFFFFFFFFFFFFFF00ull, 0xFFFFFFFFFFFFFE00ull, 0xFFFFFFFFFFFFFC00ull, 0xFFFFFFFFFFFFF800ull,
  0xFFFFFFFFFFFFF000ull, 0xFFFFFFFFFFFFE000ull, 0xFFFFFFFFFFFFC000ull, 0xFFFFFFFFFFFF8000ull,
  0xFFFFFFFFFFFF0000ull, 0xFFFFFFFFFFFE0000ull, 0xFFFFFFFFFFFC0000ull, 0xFFFFFFFFFFF80000ull,
  0xFFFFFFFFFFF00000ull, 0xFFFFFFFFFFE00000ull, 0xFFFFFFFFFFC00000ull, 0xFFFFFFFFFF800000ull,
  0xFFFFFFFFFF000000ull, 0xFFFFFFFFFE000000ull, 0xFFFFFFFFFC000000ull, 0xFFFFFFFFF8000000ull,
  0xFFFFFFFFF0000000ull, 0xFFFFFFFFE0000000ull, 0xFFFFFFFFC0000000ull, 0xFFFFFFFF80000000ull,
  0xFFFFFFFF00000000ull, 0xFFFFFFFE00000000ull, 0xFFFFFFFC00000000ull, 0xFFFFFFF800000000ull,
  0xFFFFFFF000000000ull, 0xFFFFFFE000000000ull, 0xFFFFFFC000000000ull, 0xFFFFFF8000000000ull,
  0xFFFFFF0000000000ull, 0xFFFFFE0000000000ull, 0xFFFFFC0000000000ull, 0xFFFFF80000000000ull,
  0xFFFFF00000000000ull, 0xFFFFE00000000000ull, 0xFFFFC00000000000ull, 0xFFFF800000000000ull,
  0xFFFF000000000000ull, 0xFFFE000000000000ull, 0xFFFC000000000000ull, 0xFFF8000000000000ull,
  0xFFF0000000000000ull, 0xFFE0000000000000ull, 0xFFC0000000000000ull, 0xFF80000000000000ull,
  0xFF00000000000000ull, 0xFE00000000000000ull, 0xFC00000000000000ull, 0xF800000000000000ull,
  0xF000000000000000ull, 0xE000000000000000ull, 0xC000000000000000ull, 0x8000000000000000ull,
  0x0000000000000000ull
};
const uint64_t uint64_mask_mod_pow2[] =
{
  0x0000000000000000ull, 0x0000000000000001ull, 0x0000000000000003ull, 0x0000000000000007ull,
  0x000000000000000Full, 0x000000000000001Full, 0x000000000000003Full, 0x000000000000007Full,
  0x00000000000000FFull, 0x00000000000001FFull, 0x00000000000003FFull, 0x00000000000007FFull,
  0x0000000000000FFFull, 0x0000000000001FFFull, 0x0000000000003FFFull, 0x0000000000007FFFull,
  0x000000000000FFFFull, 0x000000000001FFFFull, 0x000000000003FFFFull, 0x000000000007FFFFull,
  0x00000000000FFFFFull, 0x00000000001FFFFFull, 0x00000000003FFFFFull, 0x00000000007FFFFFull,
  0x0000000000FFFFFFull, 0x0000000001FFFFFFull, 0x0000000003FFFFFFull, 0x0000000007FFFFFFull,
  0x000000000FFFFFFFull, 0x000000001FFFFFFFull, 0x000000003FFFFFFFull, 0x000000007FFFFFFFull,
  0x00000000FFFFFFFFull, 0x00000001FFFFFFFFull, 0x00000003FFFFFFFFull, 0x00000007FFFFFFFFull,
  0x0000000FFFFFFFFFull, 0x0000001FFFFFFFFFull, 0x0000003FFFFFFFFFull, 0x0000007FFFFFFFFFull,
  0x000000FFFFFFFFFFull, 0x000001FFFFFFFFFFull, 0x000003FFFFFFFFFFull, 0x000007FFFFFFFFFFull,
  0x00000FFFFFFFFFFFull, 0x00001FFFFFFFFFFFull, 0x00003FFFFFFFFFFFull, 0x00007FFFFFFFFFFFull,
  0x0000FFFFFFFFFFFFull, 0x0001FFFFFFFFFFFFull, 0x0003FFFFFFFFFFFFull, 0x0007FFFFFFFFFFFFull,
  0x000FFFFFFFFFFFFFull, 0x001FFFFFFFFFFFFFull, 0x003FFFFFFFFFFFFFull, 0x007FFFFFFFFFFFFFull,
  0x00FFFFFFFFFFFFFFull, 0x01FFFFFFFFFFFFFFull, 0x03FFFFFFFFFFFFFFull, 0x07FFFFFFFFFFFFFFull,
  0x0FFFFFFFFFFFFFFFull, 0x1FFFFFFFFFFFFFFFull, 0x3FFFFFFFFFFFFFFFull, 0x7FFFFFFFFFFFFFFFull,
  0xFFFFFFFFFFFFFFFFull
};
/*
 * Helper functions (OPERATIVE)
 */
GEM_INLINE uint64_t integer_proportion(const float proportion,const uint64_t length) {
  if (proportion<=0.0) return 0;
  if (proportion>=1.0) return (uint64_t)proportion;
  return (uint64_t)(proportion*(float)length);
}
GEM_INLINE uint64_t integer_lower_power_of_two(uint64_t number) {
  uint64_t result = 0;
  while (number>>=1) {
    result++;
  }
  return result;
}
GEM_INLINE uint64_t integer_upper_power_of_two(uint64_t number) {
  const uint64_t lower_power_of_two = integer_lower_power_of_two(number);
  return (number | 1<<lower_power_of_two) ? lower_power_of_two+1 : lower_power_of_two;
}
/*
 * Random number generator
 */
GEM_INLINE uint64_t gem_rand_IID(const uint64_t min,const uint64_t max) {
  int n_rand = rand(); // [0, RAND_MAX]
  const uint64_t range = max - min;
  const uint64_t rem = RAND_MAX % range;
  const uint64_t sample = RAND_MAX / range;
  // Consider the small interval within remainder of RAND_MAX
  if (n_rand < RAND_MAX - rem) {
    return min + n_rand/sample;
  } else {
    return gem_rand_IID(min,max);
  }
}
/*
 * CheckSum & BitDisplay
 */
GEM_INLINE uint64_t checksum_uint64(uint64_t* mem,const uint64_t num_words) {
  uint64_t i, checksum = 0;
  for (i=0;i<num_words;++mem,++i) {
    checksum ^= *mem;
  }
  return checksum;
}
GEM_INLINE void checksum_incremental_uint64(uint64_t* const checksum,const uint64_t word) {
  *checksum ^= word;
}
GEM_INLINE void fprintf_uint64_binary(FILE* const stream,const uint64_t word) {
  char binary_string[65];
  uint64_t mask = UINT64_ONE_MASK;
  int64_t i;
  binary_string[64] = '\0';
  for (i=63;i>=0;--i) {
    binary_string[i] = (word & mask) ? '1' : '0';
    mask <<= 1;
  }
  // Print
  fprintf(stream,"%s",binary_string);
}
GEM_INLINE void fprintf_uint64_footprint(FILE* const stream,const uint64_t word) {
  uint8_t* components = (uint8_t*)(&word);
  // Print
  fprintf(stream,"%c%c%c%c",
      CHAR_TO_PRINTABLE(components[0]),CHAR_TO_PRINTABLE(components[1]),
      CHAR_TO_PRINTABLE(components[2]),CHAR_TO_PRINTABLE(components[3]));
}
/*
 * System
 */
GEM_INLINE uint64_t system_get_num_processors() {
  const uint64_t num_processors = sysconf(_SC_NPROCESSORS_ONLN);
  return num_processors ? num_processors : 1;
}
GEM_INLINE char* system_get_cwd() {
  char* const cwd = malloc(BUFFER_SIZE_2K);
  if (getcwd(cwd,BUFFER_SIZE_2K)==NULL) cwd[0]='\0';
  return cwd;
}
GEM_INLINE char* system_get_hostname() {
  char* const hostname = malloc(BUFFER_SIZE_2K);
  if (gethostname(hostname,1024)) hostname[0]='\0';
  return hostname;
}
GEM_INLINE char* system_get_user_name() {
  struct passwd *pw;
  uid_t uid = geteuid();
  pw = getpwuid(uid);
  return pw ? (pw->pw_name ? pw->pw_name : "") : "";
}

