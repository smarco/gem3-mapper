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
 * DESCRIPTION: Base module containing general purpose functions
 */

#include "system/commons.h"

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

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
 * Common numerical data processing/formating
 */
uint64_t integer_num_ciphers(const uint64_t number) {
  return (number >= 10000000000000000000ull) ? 20 :
         (number >= 1000000000000000000ull) ? 19 :
         (number >= 100000000000000000ull) ? 18 :
         (number >= 10000000000000000ull) ? 17 :
         (number >= 1000000000000000ull) ? 16 :
         (number >= 100000000000000ull) ? 15 :
         (number >= 10000000000000ull) ? 14 :
         (number >= 1000000000000ull) ? 13 :
         (number >= 100000000000ull) ? 12 :
         (number >= 10000000000ull) ? 11 :
         (number >= 1000000000ull) ? 10 :
         (number >= 100000000ull) ? 9 :
         (number >= 10000000ull) ? 8 :
         (number >= 1000000ull) ? 7 :
         (number >= 100000ull) ? 6 :
         (number >= 10000ull) ? 5 :
         (number >= 1000ull) ? 4 :
         (number >= 100ull) ? 3 :
         (number >= 10ull) ? 2 : 1;
}
uint64_t integer_proportion(const double proportion,const uint64_t length) {
  if (proportion<=0.0) return 0;
  if (proportion>=1.0) return (uint64_t)proportion;
  return (uint64_t)(proportion*(double)length);
}
uint64_t integer_lower_power_of_two(uint64_t number) {
  uint64_t result = 0;
  while (number>>=1) {
    result++;
  }
  return result;
}
uint64_t integer_upper_power_of_two(uint64_t number) {
  const uint64_t lower_power_of_two = integer_lower_power_of_two(number);
  return (number | 1<<lower_power_of_two) ? lower_power_of_two+1 : lower_power_of_two;
}
int integer_to_ascii(char* const buffer,uint64_t number) {
  // Calculate the number of digits of the number
  const uint64_t num_digits =
      (number >= 10000) ? 5 :
      (number >= 1000) ? 4 :
      (number >= 100) ? 3 :
      (number >= 10) ? 2 : 1;
  // Decompose number
  char* centinel = buffer;
  switch (num_digits) {
    case 5: {
      const uint64_t number_div10K = number / 10000;
      number -= number_div10K * 10000;
      centinel += integer_to_ascii(centinel,number_div10K);
    }
    // no break
    case 4: {
      const uint64_t div_1000 = (number * 8389UL) >> 23;
      *centinel = '0' + (char) div_1000;
      ++centinel;
      number -= div_1000 * 1000;
    }
    // no break
    case 3: {
      const uint64_t div_100 = (number * 5243UL) >> 19;
      *centinel = '0' + (char) div_100;
      ++centinel;
      number -= div_100 * 100;
    }
    // no break
    case 2: {
      const uint64_t div_10 = (number * 6554UL) >> 16;
      *centinel = '0' + (char) div_10;
      ++centinel;
      number -= div_10 * 10;
    }
    // no break
    case 1: {
      *centinel = '0' + (char) number;
      ++centinel;
    }
    default:
      break;
  }
  // Return number of ciphers written
  return centinel - buffer;
}
float gem_log2(float number) {
//  union { float val; int32_t x; } u = { number };
//  float log_2 = (float)(((u.x >> 23) & 255) - 128);
//  u.x &= ~(255 << 23);
//  u.x += 127 << 23;
//  log_2 += ((-0.34484843f) * u.val + 2.02466578f) * u.val - 0.67487759f;
//  return log_2;
//  /* Ankerl's inversion of Schraudolph's published algorithm, converted to explicit multiplication */
//  union { float f; int x; } u = { number };
//  return (u.x - 1064866805) * 8.262958405176314e-8f; /* 1 / 12102203.0; */
  return __builtin_log2(number);
}
float gem_loge(float number) {
  int * const exp_ptr = ((int*) &number);
  int x = *exp_ptr;
  const int log_2 = ((x >> 23) & 255) - 128;
  x &= ~(255 << 23);
  x += 127 << 23;
  *exp_ptr = x;
  number = ((-1.0f / 3) * number + 2) * number - 2.0f / 3;
  return ((number + log_2) * 0.69314718);
}
/*
 * Random number generator
 */
uint64_t gem_rand_IID(const uint64_t min,const uint64_t max) {
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
 * Statistical Utils
 */
double standard_normal_CDF(double x) {
  // CDF Constants
  const double a1 =  0.254829592;
  const double a2 = -0.284496736;
  const double a3 =  1.421413741;
  const double a4 = -1.453152027;
  const double a5 =  1.061405429;
  const double p  =  0.3275911;
  // Sign of x
  int sign = 1;
  if (x < 0) sign = -1;
  x = fabs(x)/sqrt(2.0);
  // A&S formula 7.1.26
  double t = 1.0/(1.0 + p*x);
  double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
  return 0.5*(1.0 + sign*y);
}
/*
 * CheckSum & BitDisplay
 */
uint64_t checksum_uint64(uint64_t* mem,const uint64_t num_words) {
  uint64_t i, checksum = 0;
  for (i=0;i<num_words;++mem,++i) {
    checksum ^= *mem;
  }
  return checksum;
}
void checksum_incremental_uint64(uint64_t* const checksum,const uint64_t word) {
  *checksum ^= word;
}
void fprintf_uint64_binary(FILE* const stream,const uint64_t word) {
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
void fprintf_uint64_footprint(FILE* const stream,const uint64_t word) {
  uint8_t* components = (uint8_t*)(&word);
  // Print
  fprintf(stream,"%c%c%c%c",
      CHAR_TO_PRINTABLE(components[0]),CHAR_TO_PRINTABLE(components[1]),
      CHAR_TO_PRINTABLE(components[2]),CHAR_TO_PRINTABLE(components[3]));
}
/*
 * System
 */
uint64_t system_get_num_processors(void) {
  const uint64_t num_processors = sysconf(_SC_NPROCESSORS_ONLN);
  return num_processors ? num_processors : 1;
}
char* system_get_cwd(void) {
  char* const cwd = malloc(BUFFER_SIZE_2K);
  if (getcwd(cwd,BUFFER_SIZE_2K)==NULL) cwd[0]='\0';
  return cwd;
}
char* system_get_hostname(void) {
  char* const hostname = malloc(BUFFER_SIZE_2K);
  if (gethostname(hostname,1024)) hostname[0]='\0';
  return hostname;
}
char* system_get_user_name(void) {
  struct passwd *pw;
  uid_t uid = geteuid();
  pw = getpwuid(uid);
  return pw ? (pw->pw_name ? pw->pw_name : "") : "";
}
void system_print_info(FILE* const stream) {
  fprintf(stream,"[GEM]>System.Info\n");
  // Current Date
  time_t current_time=time(0);
  struct tm local_time;
  localtime_r(&current_time,&local_time);
  fprintf(stream,"  => Date %4d/%d/%d %02d:%02d:%02d CET\n",
      1900+local_time.tm_year,local_time.tm_mon+1,local_time.tm_mday,
      local_time.tm_hour,local_time.tm_min,local_time.tm_sec);
  char* cwd = system_get_cwd();
  fprintf(stream,"  => CWD %s\n",cwd);
  free(cwd);
  // Host
  char* hostname = system_get_hostname();
  fprintf(stream,"  => Hostname %s\n",hostname);
  free(hostname);
  // User
  const char* const user_name = system_get_user_name();
  fprintf(stream,"  => Username %s\n",user_name);
  // CPU
  const uint64_t num_processors = system_get_num_processors();
  fprintf(stream,"  => CPUs %"PRIu64"\n",num_processors);
  // Memory
  // TODO
  // Disk
  // TODO
}
void system_get_time(struct timespec *ts) {
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(),CALENDAR_CLOCK,&cclock);
  clock_get_time(cclock,&mts);
  mach_port_deallocate(mach_task_self(),cclock);
  ts->tv_sec = mts.tv_sec;
  ts->tv_nsec = mts.tv_nsec;
#else
  clock_gettime(CLOCK_REALTIME,ts);
#endif
}

