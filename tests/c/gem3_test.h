/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2026 by Jacob Schreiber  <jacob.schreiber@umassmed.edu>
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
 * AUTHOR(S): Jacob Schreiber <jacob.schreiber@umassmed.edu> (unit tests, written with Claude Code)
 * DESCRIPTION:
 *   Minimal unit-test framework for the C drivers in tests/c/ (header only).
 *
 *   A driver is one file, tests/c/test_<name>.c, that includes this header,
 *   defines its cases with GT_TEST(name) and ends with
 *
 *     GT_CASES(GT_CASE(case_a), GT_CASE(case_b))
 *
 *   which defines main(). Running the driver prints TAP ("1..N", then
 *   "ok i - name" / "not ok i - name" with "#" diagnostic lines) and exits
 *   non-zero when any case fails. "-l" lists the cases; names given as
 *   arguments run only those cases.
 *
 *   GT_CHECK* record a failure and continue; GT_REQUIRE* end the current case.
 *   The runner calls gruntime_init(1,NULL) before the first case, as
 *   gem-mapper's main does (GEM_PROFILE builds need it).
 *   Only the first GT_MAX_REPORTED failures of a case are printed.
 *   gt_rng_* is a SplitMix64 generator so random cases are the same on every
 *   platform and compiler.
 */

#ifndef GEM3_TEST_H_
#define GEM3_TEST_H_

// gem3's headers set feature-test macros (_GNU_SOURCE), so they come before
// any system header.
#include "utils/essentials.h"

#include <inttypes.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GT_MAX_REPORTED 20

typedef struct {
  const char* name;
  void (*fn)(void);
} gt_case_t;

typedef struct {
  uint64_t checks;
  uint64_t failures;
  jmp_buf abort_case;
} gt_state_t;

static gt_state_t gt_state;

/*
 * Failure reporting
 */
static inline void gt_fail(const char* const file,const int line,const char* const fmt,...) {
  ++gt_state.failures;
  if (gt_state.failures > GT_MAX_REPORTED) return;
  fprintf(stdout,"#   %s:%d: ",file,line);
  va_list args;
  va_start(args,fmt);
  vfprintf(stdout,fmt,args);
  va_end(args);
  fputc('\n',stdout);
  if (gt_state.failures == GT_MAX_REPORTED) {
    fprintf(stdout,"#   (further failures in this case are counted, not printed)\n");
  }
}
static inline void gt_abort_case(void) {
  longjmp(gt_state.abort_case,1);
}

/*
 * Checks (continue after a failure)
 */
#define GT_CHECK(cond) do { \
  ++gt_state.checks; \
  if (!(cond)) gt_fail(__FILE__,__LINE__,"check failed: %s",#cond); \
} while (0)
#define GT_CHECK_MSG(cond,...) do { \
  ++gt_state.checks; \
  if (!(cond)) gt_fail(__FILE__,__LINE__,__VA_ARGS__); \
} while (0)
#define GT_CHECK_EQ_INT(a,b) do { \
  const int64_t gt_a_ = (int64_t)(a), gt_b_ = (int64_t)(b); \
  ++gt_state.checks; \
  if (gt_a_ != gt_b_) gt_fail(__FILE__,__LINE__,"%s == %s: %" PRId64 " != %" PRId64,#a,#b,gt_a_,gt_b_); \
} while (0)
#define GT_CHECK_EQ_UINT(a,b) do { \
  const uint64_t gt_a_ = (uint64_t)(a), gt_b_ = (uint64_t)(b); \
  ++gt_state.checks; \
  if (gt_a_ != gt_b_) gt_fail(__FILE__,__LINE__,"%s == %s: %" PRIu64 " != %" PRIu64,#a,#b,gt_a_,gt_b_); \
} while (0)
#define GT_CHECK_EQ_STR(a,b) do { \
  const char* const gt_a_ = (a); const char* const gt_b_ = (b); \
  ++gt_state.checks; \
  if (gt_a_ == NULL || gt_b_ == NULL || strcmp(gt_a_,gt_b_) != 0) \
    gt_fail(__FILE__,__LINE__,"%s == %s: \"%s\" != \"%s\"",#a,#b, \
        gt_a_ ? gt_a_ : "(null)",gt_b_ ? gt_b_ : "(null)"); \
} while (0)
#define GT_CHECK_EQ_MEM(a,b,n) do { \
  const size_t gt_n_ = (size_t)(n); \
  ++gt_state.checks; \
  if (memcmp((a),(b),gt_n_) != 0) gt_fail(__FILE__,__LINE__,"%s == %s (%zu bytes) differ",#a,#b,gt_n_); \
} while (0)

/*
 * Requirements (end the current case after a failure)
 */
#define GT_REQUIRE(cond) do { \
  ++gt_state.checks; \
  if (!(cond)) { gt_fail(__FILE__,__LINE__,"requirement failed: %s",#cond); gt_abort_case(); } \
} while (0)
#define GT_REQUIRE_MSG(cond,...) do { \
  ++gt_state.checks; \
  if (!(cond)) { gt_fail(__FILE__,__LINE__,__VA_ARGS__); gt_abort_case(); } \
} while (0)

/*
 * Note a diagnostic line without failing
 */
static inline void gt_note(const char* const fmt,...) {
  fputs("# ",stdout);
  va_list args;
  va_start(args,fmt);
  vfprintf(stdout,fmt,args);
  va_end(args);
  fputc('\n',stdout);
}

/*
 * Deterministic random numbers (SplitMix64)
 */
typedef struct {
  uint64_t state;
} gt_rng_t;
static inline gt_rng_t gt_rng_new(const uint64_t seed) {
  gt_rng_t rng;
  rng.state = seed;
  return rng;
}
static inline uint64_t gt_rng_u64(gt_rng_t* const rng) {
  uint64_t z = (rng->state += UINT64_C(0x9E3779B97F4A7C15));
  z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
  return z ^ (z >> 31);
}
// Uniform in [0,n) (n > 0), without modulo bias
static inline uint64_t gt_rng_below(gt_rng_t* const rng,const uint64_t n) {
  const uint64_t limit = UINT64_MAX - (UINT64_MAX % n);
  uint64_t x;
  do { x = gt_rng_u64(rng); } while (x >= limit);
  return x % n;
}
// Uniform in [lo,hi]
static inline uint64_t gt_rng_range(gt_rng_t* const rng,const uint64_t lo,const uint64_t hi) {
  return lo + gt_rng_below(rng,hi-lo+1);
}
// True with probability num/den
static inline bool gt_rng_chance(gt_rng_t* const rng,const uint64_t num,const uint64_t den) {
  return gt_rng_below(rng,den) < num;
}

/*
 * Runner
 */
#define GT_TEST(name) static void name(void)
#define GT_CASE(name) { #name, name }
static inline int gt_main(
    const int argc,
    char** const argv,
    const gt_case_t* const cases,
    const int num_cases) {
  int i, j;
  if (argc > 1 && strcmp(argv[1],"-l") == 0) {
    for (i=0;i<num_cases;++i) printf("%s\n",cases[i].name);
    return 0;
  }
  // Select cases
  int num_selected = 0;
  bool* const selected = calloc((size_t)num_cases+1,sizeof(bool));
  for (i=0;i<num_cases;++i) {
    if (argc <= 1) {
      selected[i] = true;
    } else {
      for (j=1;j<argc;++j) {
        if (strcmp(argv[j],cases[i].name) == 0) selected[i] = true;
      }
    }
    if (selected[i]) ++num_selected;
  }
  for (j=1;j<argc;++j) {
    bool found = false;
    for (i=0;i<num_cases;++i) found = found || strcmp(argv[j],cases[i].name)==0;
    if (!found) {
      fprintf(stdout,"Bail out! no case named '%s'\n",argv[j]);
      free(selected);
      return 2;
    }
  }
  // gem3 runtime, as gem-mapper's main sets it up (profiler counters in
  // GEM_PROFILE builds, thread registry, error signal handlers)
  gruntime_init(1,NULL);
  // Run
  printf("1..%d\n",num_selected);
  fflush(stdout);
  int num_failed = 0, n = 0;
  uint64_t total_checks = 0;
  for (i=0;i<num_cases;++i) {
    if (!selected[i]) continue;
    ++n;
    gt_state.checks = 0;
    gt_state.failures = 0;
    if (setjmp(gt_state.abort_case) == 0) {
      cases[i].fn();
    }
    total_checks += gt_state.checks;
    if (gt_state.failures == 0) {
      printf("ok %d - %s (%" PRIu64 " checks)\n",n,cases[i].name,gt_state.checks);
    } else {
      ++num_failed;
      printf("not ok %d - %s (%" PRIu64 " of %" PRIu64 " checks failed)\n",
          n,cases[i].name,gt_state.failures,gt_state.checks);
    }
    fflush(stdout);
  }
  printf("# %d of %d cases failed, %" PRIu64 " checks\n",num_failed,num_selected,total_checks);
  gruntime_destroy();
  free(selected);
  return num_failed ? 1 : 0;
}
#define GT_CASES(...) \
  static const gt_case_t gt_cases[] = { __VA_ARGS__ }; \
  int main(int argc,char** argv) { \
    return gt_main(argc,argv,gt_cases,(int)(sizeof(gt_cases)/sizeof(gt_cases[0]))); \
  }

#endif /* GEM3_TEST_H_ */
