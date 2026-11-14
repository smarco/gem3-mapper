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
 *   Area D: number formatting and the small utilities the writers and the
 *   command line use, each against the C library or a direct reference.
 *     - integer_to_ascii, integer_num_ciphers against snprintf("%llu")
 *     - bprintf_* (output buffer) and bofprintf_* (buffered output file)
 *       against snprintf of the same values
 *     - integer_proportion, gem_log2, gem_loge (its quadratic approximation
 *       written out independently, and its error against log()), gem_rand_IID
 *     - string_buffer (string_t), gem_str* against libc
 *     - vector_t, ihash_t against plain arrays
 *     - options_menu adaptors on a small option table
 *     - parse_integer_system ('<n>', 'c', '<x>c') for -t and the models
 *
 *   Internal signatures called (a speed iteration may change them; update
 *   this driver with them):
 *     int integer_to_ascii(char* buffer,uint64_t number)
 *     uint64_t integer_num_ciphers(uint64_t), integer_proportion(double,uint64_t)
 *     float gem_log2(float), float gem_loge(float), uint64_t gem_rand_IID(uint64_t,uint64_t)
 *     output_buffer_new/clear/delete/get_used/get_buffer/set_state/get_state,
 *       bprintf_uint64/int64/char/buffer(output_buffer_t*,...)
 *     output_file_new(char* file_name,uint64_t num_buffers,uint64_t buffer_size), output_file_close,
 *       buffered_output_file_new/request_buffer/reserve/dump_buffer/close,
 *       bofprintf_uint64/int64/char/string
 *     string_init/resize/clear/destroy/set_buffer/get_buffer/set_length/get_length/
 *       append_char/append_eos/right_append_buffer/right_append_string/is_null/copy,
 *       gem_strncpy/strcmp/streq/strcaseeq/strcat/strlen
 *     vector_new_/vector_reserve/vector_delete (+ vector_insert, vector_get_* macros)
 *     ihash_new/clear/delete/get_ihash_element/insert_element/get_element/sort_by_key
 *     options_get_num_options, options_adaptor_getopt, options_adaptor_getopt_short,
 *       options_fprint_menu
 *     int parse_integer_system(char* argument,uint64_t* value)   (mapper_arguments.c)
 *   and reads ihash_t.head and ihash_element_t (key, hh.next) to walk a
 *   sorted hash, as report_stats.c does.
 */

#include "gem3_test.h"
#include <math.h>
#include <getopt.h>
#include "system/commons.h"
#include "system/mm.h"
#include "io/output_buffer.h"
#include "io/output_file.h"
#include "io/buffered_output_file.h"
#include "utils/string_buffer.h"
#include "utils/vector.h"
#include "utils/ihash.h"
#include "utils/options_menu.h"

int parse_integer_system(char* const argument,uint64_t* const value);

/*
 * Numbers to text
 */
static void check_ascii(const uint64_t n) {
  char got[32], want[32];
  memset(got,'#',sizeof(got));
  const int len = integer_to_ascii(got,n);
  const int wlen = snprintf(want,sizeof(want),"%" PRIu64,n);
  GT_CHECK_MSG(len==wlen && memcmp(got,want,wlen)==0 && got[len]=='#',
      "integer_to_ascii(%" PRIu64 ") wrote %.*s (%d chars)",n,len>0&&len<32?len:0,got,len);
  GT_CHECK_MSG(integer_num_ciphers(n)==(uint64_t)wlen,
      "integer_num_ciphers(%" PRIu64 ") = %" PRIu64,n,integer_num_ciphers(n));
}
GT_TEST(integer_to_ascii_matches_printf) {
  uint64_t n, p;
  for (n=0;n<200000;++n) check_ascii(n);
  for (p=1;p<=UINT64_C(10000000000000000000);p*=10) {
    check_ascii(p-1); check_ascii(p); check_ascii(p+1);
    if (p > UINT64_MAX/10) break;
  }
  check_ascii(UINT64_MAX); check_ascii(UINT64_MAX-1);
  check_ascii(UINT64_C(4294967295)); check_ascii(UINT64_C(4294967296));
  gt_rng_t rng = gt_rng_new(11);
  int i;
  for (i=0;i<200000;++i) {
    const uint64_t x = gt_rng_u64(&rng);
    check_ascii(x >> gt_rng_below(&rng,64));
  }
}

GT_TEST(output_buffer_printers_match_printf) {
  output_buffer_t* const ob = output_buffer_new(1<<20);
  GT_CHECK_EQ_UINT(output_buffer_get_used(ob),0);
  GT_CHECK(output_buffer_get_state(ob)==OUTPUT_BUFFER_FREE);
  output_buffer_set_state(ob,OUTPUT_BUFFER_BUSY);
  GT_CHECK(output_buffer_get_state(ob)==OUTPUT_BUFFER_BUSY);
  char* const want = malloc(1<<20);
  size_t wlen = 0;
  gt_rng_t rng = gt_rng_new(5);
  int i;
  for (i=0;i<20000;++i) {
    switch (gt_rng_below(&rng,4)) {
      case 0: {
        const uint64_t v = gt_rng_u64(&rng) >> gt_rng_below(&rng,64);
        bprintf_uint64(ob,v);
        wlen += sprintf(want+wlen,"%" PRIu64,v);
        break;
      }
      case 1: {
        int64_t v = (int64_t)(gt_rng_u64(&rng) >> (1+gt_rng_below(&rng,63)));
        if (gt_rng_below(&rng,2)) v = -v;
        bprintf_int64(ob,v);
        wlen += sprintf(want+wlen,"%" PRId64,v);
        break;
      }
      case 2: {
        const char c = (char)('!'+gt_rng_below(&rng,90));
        bprintf_char(ob,c);
        want[wlen++] = c;
        break;
      }
      default: {
        const char* const s = "\tNM:i:";
        const int n = (int)gt_rng_below(&rng,7);
        bprintf_buffer(ob,n,s);
        memcpy(want+wlen,s,n);
        wlen += n;
        break;
      }
    }
  }
  bprintf_int64(ob,INT64_MIN+1); wlen += sprintf(want+wlen,"%" PRId64,INT64_MIN+1);
  bprintf_int64(ob,-1); wlen += sprintf(want+wlen,"%d",-1);
  bprintf_int64(ob,0); want[wlen++] = '0';
  GT_REQUIRE(output_buffer_get_used(ob)==wlen);
  GT_CHECK_EQ_MEM(output_buffer_get_buffer(ob),want,wlen);
  output_buffer_clear(ob);
  GT_CHECK_EQ_UINT(output_buffer_get_used(ob),0);
  bprintf_char(ob,'x');
  GT_CHECK(output_buffer_get_buffer(ob)[0]=='x' && output_buffer_get_used(ob)==1);
  output_buffer_delete(ob);
  free(want);
}

GT_TEST(buffered_output_file_printers_reach_the_stream) {
  char path[] = "/tmp/gem3-test-d-format-XXXXXX";
  const int fd = mkstemp(path);
  GT_REQUIRE(fd >= 0);
  close(fd);
  output_file_t* const of = output_file_new(path,2,4096);
  buffered_output_file_t* const bof = buffered_output_file_new(of);
  buffered_output_file_request_buffer(bof,0);
  static char want[65536], got[65536];
  size_t wlen = 0;
  int i;
  for (i=0;i<1200;++i) {
    buffered_output_file_reserve(bof,64);
    bofprintf_uint64(bof,(uint64_t)i*977);
    bofprintf_char(bof,'\t');
    bofprintf_int64(bof,-(int64_t)i);
    bofprintf_string(bof,3,"ab\n");
    wlen += sprintf(want+wlen,"%d\t%d", i*977, -i);
    memcpy(want+wlen,"ab\n",3); wlen += 3;
  }
  buffered_output_file_close(bof);
  output_file_close(of);
  FILE* const f = fopen(path,"rb");
  GT_REQUIRE(f!=NULL);
  const size_t n = fread(got,1,sizeof(got),f);
  fclose(f);
  unlink(path);
  GT_CHECK_EQ_UINT(n,wlen);
  GT_CHECK_EQ_MEM(got,want,wlen);
}

/*
 * Numerical helpers
 */
GT_TEST(integer_proportion_cases) {
  GT_CHECK_EQ_UINT(integer_proportion(0.0,100),0);
  GT_CHECK_EQ_UINT(integer_proportion(-0.5,100),0);
  GT_CHECK_EQ_UINT(integer_proportion(0.12,100),12);
  GT_CHECK_EQ_UINT(integer_proportion(0.04,150),6);
  GT_CHECK_EQ_UINT(integer_proportion(0.5,7),3);    // truncated
  GT_CHECK_EQ_UINT(integer_proportion(1.0,100),1);  // >=1 is an absolute value
  GT_CHECK_EQ_UINT(integer_proportion(5.0,100),5);
  GT_CHECK_EQ_UINT(integer_proportion(12.9,100),12);
}

/* gem_loge: x = m*2^e with m in [1,2); 1+log2(m) ~ (-m/3+2)*m-2/3, times ln 2 */
static double loge_reference(const float x) {
  int e;
  const double f = frexp((double)x,&e); // x = f*2^e, f in [0.5,1)
  const double m = 2.0*f;
  const double l2 = ((-1.0/3.0)*m+2.0)*m - 2.0/3.0 + (double)(e-1) - 1.0;
  return l2*0.69314718;
}
GT_TEST(gem_loge_approximation) {
  gt_rng_t rng = gt_rng_new(3);
  double max_err = 0.0;
  int i;
  for (i=0;i<100000;++i) {
    const float x = (float)(1.0 + (double)gt_rng_below(&rng,1u<<30)/16.0);
    const double got = gem_loge(x);
    const double ref = loge_reference(x);
    GT_CHECK_MSG(fabs(got-ref) <= 1e-5*fabs(ref)+1e-5,"gem_loge(%g) = %.9g, formula %.9g",x,got,ref);
    const double err = fabs(got-log((double)x));
    if (err > max_err) max_err = err;
  }
  // The quadratic's own error (0.0067 at most here on base)
  GT_CHECK_MSG(max_err < 0.007,"max |gem_loge-log| = %g",max_err);
  GT_CHECK(fabs(gem_loge(4.0f)-log(4.0)) < 1e-6);
  GT_CHECK(fabs(gem_loge(1.0f)) < 1e-6);
}

GT_TEST(gem_log2_matches_log2) {
  gt_rng_t rng = gt_rng_new(4);
  int i;
  for (i=0;i<100000;++i) {
    const float x = (float)(1.0 + (double)gt_rng_below(&rng,1u<<31));
    const double want = log2((double)x);
    GT_CHECK_MSG(fabs(gem_log2(x)-want) <= 1e-6*want+1e-6,"gem_log2(%g)",x);
  }
  GT_CHECK(gem_log2(1024.0f)==10.0f);
}

GT_TEST(gem_rand_iid_in_range) {
  srand(7);
  int i;
  for (i=0;i<20000;++i) {
    const uint64_t lo = (uint64_t)(i % 50);
    const uint64_t hi = lo + 1 + (uint64_t)(i % 977);
    const uint64_t r = gem_rand_IID(lo,hi);
    GT_CHECK_MSG(r>=lo && r<hi,"gem_rand_IID(%" PRIu64 ",%" PRIu64 ") = %" PRIu64,lo,hi,r);
  }
}

/*
 * Strings
 */
GT_TEST(string_buffer_operations) {
  int pass;
  mm_slab_t* const slab = mm_slab_new_(BUFFER_SIZE_1M,BUFFER_SIZE_1M,MM_UNLIMITED_MEM);
  mm_allocator_t* const mma = mm_allocator_new(slab);
  for (pass=0;pass<2;++pass) {
    mm_allocator_t* const a = pass ? mma : NULL;
    string_t s, t;
    string_init(&s,4,a);
    GT_CHECK(string_is_null(&s));
    GT_CHECK_EQ_UINT(string_get_length(&s),0);
    GT_CHECK_EQ_STR(string_get_buffer(&s),"");
    // Append characters past the initial allocation
    char want[2048];
    int i;
    for (i=0;i<1500;++i) {
      want[i] = (char)('a'+i%26);
      string_append_char(&s,want[i]);
    }
    string_append_eos(&s);
    want[1500] = 0;
    GT_CHECK_EQ_UINT(string_get_length(&s),1500);
    GT_CHECK_EQ_STR(string_get_buffer(&s),want);
    GT_CHECK(!string_is_null(&s));
    // Right append of a buffer and of a string
    string_right_append_buffer(&s,"XYZ123",3);
    GT_CHECK_EQ_UINT(string_get_length(&s),1503);
    GT_CHECK(memcmp(string_get_buffer(&s)+1500,"XYZ",4)==0);
    string_init(&t,1,a);
    string_set_buffer(&t,"ACGTN",5);
    GT_CHECK_EQ_STR(string_get_buffer(&t),"ACGTN");
    string_right_append_string(&s,&t);
    GT_CHECK_EQ_UINT(string_get_length(&s),1508);
    GT_CHECK_EQ_STR(string_get_buffer(&s)+1500,"XYZACGTN");
    // Copy, resize keeping content, set length, clear
    string_copy(&t,&s);
    GT_CHECK_EQ_UINT(string_get_length(&t),1508);
    GT_CHECK_EQ_MEM(string_get_buffer(&t),string_get_buffer(&s),1508);
    string_set_buffer(&t,"abc",3);
    string_resize(&t,5000,true);
    GT_CHECK_EQ_MEM(string_get_buffer(&t),"abc",3);
    string_set_length(&t,2);
    GT_CHECK_EQ_UINT(string_get_length(&t),2);
    string_clear(&t);
    GT_CHECK(string_is_null(&t) && string_get_length(&t)==0);
    string_set_buffer(&t,"",0);
    GT_CHECK(string_is_null(&t));
    string_destroy(&t);
    string_destroy(&s);
  }
  mm_allocator_delete(mma);
  mm_slab_delete(slab);
  GT_CHECK(string_is_null(NULL));
}

GT_TEST(gem_string_helpers_match_libc) {
  char buf[16];
  memset(buf,'#',sizeof(buf));
  gem_strncpy(buf,"chr10extra",5);
  GT_CHECK_EQ_STR(buf,"chr10");
  GT_CHECK(gem_strcmp("chr1","chr2") < 0 && gem_strcmp("chr2","chr1") > 0 && gem_strcmp("a","a")==0);
  GT_CHECK(gem_streq("NC_001416.1","NC_001416.1") && !gem_streq("chr1","chr1 "));
  GT_CHECK(gem_strcaseeq("Offset-33","OFFSET-33") && gem_strcaseeq("sam","SAM") && !gem_strcaseeq("SAM","MAP"));
  char* const cat = gem_strcat("gem-","mapper");
  GT_CHECK_EQ_STR(cat,"gem-mapper");
  mm_free(cat);
  GT_CHECK_EQ_UINT(gem_strlen(""),0);
  GT_CHECK_EQ_UINT(gem_strlen("chrShort1"),9);
}

/*
 * Containers
 */
GT_TEST(vector_grows_and_keeps_elements) {
  vector_t* const v = vector_new(2,uint64_t);
  GT_CHECK_EQ_UINT(vector_get_used(v),0);
  uint64_t i;
  for (i=0;i<100000;++i) vector_insert(v,i*i+7,uint64_t);
  GT_REQUIRE(vector_get_used(v)==100000);
  const uint64_t* const mem = vector_get_mem(v,uint64_t);
  for (i=0;i<100000;++i) GT_CHECK_MSG(mem[i]==i*i+7,"element %" PRIu64,i);
  GT_CHECK_EQ_UINT(*vector_get_last_elm(v,uint64_t),UINT64_C(99999)*99999+7);
  vector_clear(v);
  vector_reserve(v,300000,true);
  const uint64_t* const z = vector_get_mem(v,uint64_t);
  for (i=0;i<300000;++i) {
    if (z[i]!=0) { GT_CHECK_MSG(false,"reserve(zero_mem) left %" PRIu64 " at %" PRIu64,z[i],i); break; }
  }
  vector_delete(v);
}

GT_TEST(ihash_insert_get_sort) {
  ihash_t* const h = ihash_new(NULL);
  static int64_t vals[4000];
  int64_t keys[4000];
  gt_rng_t rng = gt_rng_new(9);
  int i, n = 0;
  for (i=0;i<4000;++i) {
    // keys as the report uses them: read lengths, distances, template lengths
    const int64_t k = (int64_t)gt_rng_below(&rng,20000) - 2000;
    vals[i] = i;
    if (ihash_get_ihash_element(h,k)==NULL) keys[n++] = k;
    ihash_insert_element(h,k,vals+i); // a repeated key replaces the element
  }
  GT_CHECK_EQ_UINT(ihash_get_num_elements(h),n);
  // Last insertion wins
  gt_rng_t rng2 = gt_rng_new(9);
  int64_t last[22000];
  for (i=0;i<22000;++i) last[i] = -1;
  for (i=0;i<4000;++i) last[(int64_t)gt_rng_below(&rng2,20000)] = i;
  for (i=0;i<20000;++i) {
    const int64_t* const e = ihash_get_element(h,(int64_t)i-2000);
    if (last[i] < 0) {
      GT_CHECK(e==NULL);
    } else {
      GT_CHECK_MSG(e!=NULL && *e==last[i],"key %d",i-2000);
    }
  }
  // Sorted iteration visits every key once, ascending
  ihash_sort_by_key(h);
  int64_t prev = INT64_MIN;
  int count = 0;
  ihash_element_t* it;
  for (it=h->head;it!=NULL;it=it->hh.next) {
    GT_CHECK(it->key > prev);
    prev = it->key;
    ++count;
  }
  GT_CHECK_EQ_INT(count,n);
  ihash_clear(h);
  GT_CHECK_EQ_UINT(ihash_get_num_elements(h),0);
  GT_CHECK(ihash_get_element(h,keys[0])==NULL);
  ihash_delete(h);
}

/*
 * Options menu
 */
static option_t test_options[] = {
  { 'I', "index", REQUIRED, TYPE_STRING, 1, VISIBILITY_USER, "<file>", "(index)" },
  { 'v', "verbose", OPTIONAL, TYPE_STRING, 1, VISIBILITY_USER, "'quiet'|'user'", "(default=user)" },
  { 300, "mapping-mode", REQUIRED, TYPE_STRING, 2, VISIBILITY_ADVANCED, "'fast'", "" },
  { 'p', "paired", NO_ARGUMENT, TYPE_NONE, 2, VISIBILITY_USER, "", "" },
  { 400, "secret", NO_ARGUMENT, TYPE_NONE, 3, VISIBILITY_DEVELOPER, "", "(dev)" },
  {  0, "", 0, 0, 0, false, "", ""}
};
static char* test_groups[] = { "Null", "I/O", "Search", "Misc" };

GT_TEST(options_menu_adaptors) {
  GT_CHECK_EQ_UINT(options_get_num_options(test_options),5);
  struct option* const lo = options_adaptor_getopt(test_options);
  GT_CHECK_EQ_STR(lo[0].name,"index");
  GT_CHECK(lo[0].has_arg==required_argument && lo[0].val=='I' && lo[0].flag==NULL);
  GT_CHECK(lo[1].has_arg==optional_argument && lo[1].val=='v');
  GT_CHECK(lo[2].has_arg==required_argument && lo[2].val==300);
  GT_CHECK(lo[3].has_arg==no_argument && lo[3].val=='p');
  GT_CHECK(lo[5].name==NULL && lo[5].val==0);
  string_t* const so = options_adaptor_getopt_short(test_options);
  GT_CHECK_EQ_STR(string_get_buffer(so),"I:v::p");
  // getopt with the adapted tables
  char* argv[] = { "prog", "-I", "x.gem", "--mapping-mode", "fast", "-vquiet", "--secret", "-p", NULL };
  int argc = 8, opt, idx;
  int seen[5] = {0,0,0,0,0};
  optind = 1;
  while ((opt=getopt_long(argc,argv,string_get_buffer(so),lo,&idx))!=-1) {
    if (opt=='I') { seen[0]=1; GT_CHECK_EQ_STR(optarg,"x.gem"); }
    else if (opt=='v') { seen[1]=1; GT_CHECK_EQ_STR(optarg,"quiet"); }
    else if (opt==300) { seen[2]=1; GT_CHECK_EQ_STR(optarg,"fast"); }
    else if (opt=='p') seen[3]=1;
    else if (opt==400) seen[4]=1;
    else GT_CHECK_MSG(false,"unexpected option %d",opt);
  }
  GT_CHECK(seen[0]&&seen[1]&&seen[2]&&seen[3]&&seen[4]);
  optind = 1;
  string_destroy(so);
  mm_free(so);
  mm_free(lo);
  // Menu text at the advanced level: no developer option
  FILE* const f = tmpfile();
  options_fprint_menu(f,test_options,test_groups,true,VISIBILITY_ADVANCED);
  rewind(f);
  char text[1024];
  const size_t n = fread(text,1,sizeof(text)-1,f);
  text[n] = 0;
  fclose(f);
  GT_CHECK_EQ_STR(text,
      "    [I/O]\n"
      "      --index|-I <file> (index)\n"
      "      --verbose|-v 'quiet'|'user' (default=user)\n"
      "    [Search]\n"
      "      --mapping-mode 'fast'\n"
      "      --paired|-p \n");
}

/*
 * -t / model counts
 */
GT_TEST(parse_integer_system_forms) {
  const uint64_t cores = system_get_num_processors();
  uint64_t v = 12345;
  char a[][16] = { "4", "c", "C", "2c", "0.5c", "1.5", "64" };
  const uint64_t want[] = { 4, cores, cores, 2*cores, (uint64_t)(0.5*(double)cores), 1, 64 };
  int i;
  for (i=0;i<7;++i) {
    v = 12345;
    GT_CHECK_MSG(parse_integer_system(a[i],&v)==0 && v==want[i],
        "'%s' -> %" PRIu64 " (want %" PRIu64 ")",a[i],v,want[i]);
  }
  GT_CHECK_EQ_UINT(cores,(uint64_t)sysconf(_SC_NPROCESSORS_ONLN));
}

GT_CASES(
  GT_CASE(integer_to_ascii_matches_printf),
  GT_CASE(output_buffer_printers_match_printf),
  GT_CASE(buffered_output_file_printers_reach_the_stream),
  GT_CASE(integer_proportion_cases),
  GT_CASE(gem_loge_approximation),
  GT_CASE(gem_log2_matches_log2),
  GT_CASE(gem_rand_iid_in_range),
  GT_CASE(string_buffer_operations),
  GT_CASE(gem_string_helpers_match_libc),
  GT_CASE(vector_grows_and_keeps_elements),
  GT_CASE(ihash_insert_get_sort),
  GT_CASE(options_menu_adaptors),
  GT_CASE(parse_integer_system_forms)
)
