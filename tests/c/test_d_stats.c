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
 *   Area D: statistics.
 *     - gem_counter_t (COUNTER_*) against a two-pass mean/variance
 *     - the template-length estimate of the paired mapper (mapper_stats_*),
 *       which starts from a prior of <samples> virtual observations with
 *       mean (min+max)/2 and standard deviation mean/6 and adds observed
 *       lengths; checked against the pooled (Chan et al.) mean and sum of
 *       squares and the formulas of the getters
 *     - the report's per-thread statistics: init_mapping_stats,
 *       update_counts, get_read_type, merge_mapping_stats/merge_ihash and the
 *       JSON printers
 *     - the timers (TIMER_*) against nanosleep
 *
 *   Internal signatures called (a speed iteration may change them; update
 *   this driver with them):
 *     COUNTER_RESET/ADD/GET_TOTAL/GET_NUM_SAMPLES/GET_MEAN/GET_VARIANCE/GET_STDDEV(gem_counter_t*,...)
 *     mapper_stats_t* mapper_stats_new(void); mapper_stats_clear/delete(mapper_stats_t*)
 *     mapper_stats_template_init(mapper_stats_t*,uint64_t samples,uint64_t min,uint64_t max)
 *     mapper_stats_template_length_sample(mapper_stats_t*,uint64_t)
 *     mapper_stats_template_length_get_num_samples/get_ci_min_samples(..,margin)/is_reliable/
 *       get_mean/get_stddev/get_expected_max/get_expected_min/get_sigma_dev(..,length)
 *     init_mapping_stats(mapping_stats_t*), merge_mapping_stats(mapping_stats_t*,mapping_stats_t*,uint64_t)
 *     void update_counts(sequence_t*,mapping_stats_t*,int end)
 *     int get_read_type(match_trace_t*,char** control_seqs)
 *     void output_json_uint_element(FILE*,char*,uint64_t,int indent,bool last)
 *     void output_json_uint_array(FILE*,char*,uint64_t*,int n,int indent,bool last)
 *     TIMER_START/STOP/PAUSE/CONTINUE/RESET/RESTART/GET_TOTAL_NS(gem_timer_t*)
 *   and reads and writes these fields: mapping_stats_t (reads, BSreads,
 *   unmapped, correct_pairs, base_counts, hist_mapq, read_length_dist,
 *   distance_dist, insert_size_dist), sequence_t.read and
 *   match_trace_t.sequence_name.
 */

#include "gem3_test.h"
#include <math.h>
#include <time.h>
#include "stats/stats_counter.h"
#include "mapper/mapper_stats.h"
#include "stats/report_stats.h"
#include "profiler/profiler_timer.h"
#include "text/sequence.h"
#include "matches/match_trace.h"

void update_counts(sequence_t* const seq_read,mapping_stats_t* const mstats,const int end);
int get_read_type(match_trace_t* const match,char** control_seqs);
void merge_ihash(ihash_t* ihash1,ihash_t* ihash2);
void output_json_uint_element(FILE* fp,char* key,uint64_t value,int indent,bool last);
void output_json_uint_array(FILE* fp,char* key,uint64_t* values,int n,int indent,bool last);

static bool close_to(const double a,const double b) {
  return fabs(a-b) <= 1e-9*fmax(1.0,fmax(fabs(a),fabs(b)));
}

/*
 * Counters
 */
GT_TEST(counter_mean_variance_min_max) {
  gt_rng_t rng = gt_rng_new(21);
  int trial;
  for (trial=0;trial<50;++trial) {
    gem_counter_t c;
    COUNTER_RESET(&c);
    GT_CHECK_EQ_UINT(COUNTER_GET_TOTAL(&c),0);
    GT_CHECK_EQ_UINT(COUNTER_GET_NUM_SAMPLES(&c),0);
    GT_CHECK(COUNTER_GET_VARIANCE(&c)==0.0);
    const int n = 1 + (int)gt_rng_below(&rng,400);
    uint64_t v[400], total = 0, mn = UINT64_MAX, mx = 0;
    int i;
    for (i=0;i<n;++i) {
      v[i] = gt_rng_below(&rng,trial%2 ? 1000 : 10000000);
      COUNTER_ADD(&c,v[i]);
      total += v[i];
      mn = MIN(mn,v[i]);
      mx = MAX(mx,v[i]);
    }
    double mean = (double)total/n, ss = 0.0;
    for (i=0;i<n;++i) ss += ((double)v[i]-mean)*((double)v[i]-mean);
    const double var = n > 1 ? ss/(n-1) : 0.0;
    GT_CHECK_EQ_UINT(COUNTER_GET_TOTAL(&c),total);
    GT_CHECK_EQ_UINT(COUNTER_GET_NUM_SAMPLES(&c),n);
    GT_CHECK_EQ_UINT(COUNTER_GET_MIN(&c),mn);
    GT_CHECK_EQ_UINT(COUNTER_GET_MAX(&c),mx);
    GT_CHECK_MSG(close_to(COUNTER_GET_MEAN(&c),mean),"mean %g vs %g",COUNTER_GET_MEAN(&c),mean);
    GT_CHECK_MSG(fabs(COUNTER_GET_VARIANCE(&c)-var) <= 1e-7*fmax(1.0,var),"var %.12g vs %.12g",COUNTER_GET_VARIANCE(&c),var);
    GT_CHECK_MSG(fabs(COUNTER_GET_STDDEV(&c)-sqrt(var)) <= 1e-7*fmax(1.0,sqrt(var)),"sd");
  }
}

/*
 * Template-length estimate
 */
typedef struct { double n, mean, ss; uint64_t total; } pooled_t;
static void pooled_add(pooled_t* const p,const uint64_t x) {
  // Chan et al.: combine (n,mean,ss) with one observation
  const double n1 = p->n + 1.0;
  const double d = (double)x - p->mean;
  p->ss += d*d*p->n/n1;
  p->mean += d/n1;
  p->n = n1;
  p->total += x;
}
GT_TEST(template_length_estimate) {
  uint64_t samples_list[] = { 100, 37 };
  uint64_t bounds[][2] = { {0,800}, {150,501}, {200,200} };
  gt_rng_t rng = gt_rng_new(31);
  int a, b;
  for (a=0;a<2;++a) for (b=0;b<3;++b) {
    const uint64_t n0 = samples_list[a], lo = bounds[b][0], hi = bounds[b][1];
    mapper_stats_t* const ms = mapper_stats_new();
    GT_CHECK_EQ_UINT(mapper_stats_template_length_get_num_samples(ms),0);
    mapper_stats_template_init(ms,n0,lo,hi);
    const double m0 = (double)(lo+hi)/2.0, d0 = m0/6.0;
    pooled_t p = { (double)n0, m0, d0*d0*(double)n0, (uint64_t)((double)n0*m0) };
    int k;
    for (k=0;k<=300;++k) {
      if (k > 0) {
        const uint64_t x = 250 + gt_rng_below(&rng,100) + (k%7==0 ? gt_rng_below(&rng,500) : 0);
        mapper_stats_template_length_sample(ms,x);
        pooled_add(&p,x);
      }
      const uint64_t n = mapper_stats_template_length_get_num_samples(ms);
      GT_REQUIRE_MSG(n==n0+(uint64_t)k,"samples %" PRIu64,n);
      // Mean is the integer total over the samples; deviation from the pooled sum of squares
      const double mean = (double)p.total/(double)n;
      const double sd = sqrt(p.ss/(double)(n-1));
      GT_CHECK_MSG(close_to(mapper_stats_template_length_get_mean(ms),mean),"mean %.12g vs %.12g",
          mapper_stats_template_length_get_mean(ms),mean);
      GT_CHECK_MSG(fabs(mapper_stats_template_length_get_stddev(ms)-sd) <= 1e-6*fmax(1.0,sd),
          "k=%d sd %.12g vs %.12g",k,mapper_stats_template_length_get_stddev(ms),sd);
      const double got_sd = mapper_stats_template_length_get_stddev(ms);
      // Truncations to integers: allow one unit for rounding differences
      // between compilers at an integer boundary
      const int64_t emax = (int64_t)(mean+6.0*got_sd);
      const int64_t emin = (int64_t)(mean > 6.0*got_sd ? mean-6.0*got_sd : 0.0);
      GT_CHECK(llabs((int64_t)mapper_stats_template_length_get_expected_max(ms)-emax) <= 1);
      GT_CHECK(llabs((int64_t)mapper_stats_template_length_get_expected_min(ms)-emin) <= 1);
      const uint64_t ci = mapper_stats_template_length_get_ci_min_samples(ms,20);
      const double f = 1.96*got_sd/20.0;
      GT_CHECK_MSG(llabs((int64_t)ci-(int64_t)((uint64_t)(f*f)+1)) <= 1,"ci %" PRIu64,ci);
      GT_CHECK(mapper_stats_template_length_is_reliable(ms)==(n > ci && n >= 100));
      const uint64_t probe = 260 + (uint64_t)k;
      GT_CHECK(close_to(mapper_stats_template_length_get_sigma_dev(ms,probe),fabs((double)probe-mean)/got_sd));
    }
    // Re-initializing resets the estimate to the prior
    mapper_stats_template_init(ms,n0,lo,hi);
    GT_CHECK_EQ_UINT(mapper_stats_template_length_get_num_samples(ms),n0);
    GT_CHECK(close_to(mapper_stats_template_length_get_mean(ms),(double)((uint64_t)((double)n0*m0))/(double)n0));
    // Invalid bounds leave the estimate as it is
    mapper_stats_template_length_sample(ms,1000);
    mapper_stats_template_init(ms,n0,hi+1,hi);
    GT_CHECK_EQ_UINT(mapper_stats_template_length_get_num_samples(ms),n0+1);
    mapper_stats_template_init(ms,n0,UINT64_MAX,hi);
    GT_CHECK_EQ_UINT(mapper_stats_template_length_get_num_samples(ms),n0+1);
    mapper_stats_clear(ms);
    GT_CHECK_EQ_UINT(mapper_stats_template_length_get_num_samples(ms),0);
    mapper_stats_delete(ms);
  }
}

/*
 * Report statistics
 */
GT_TEST(update_counts_counts_bases_and_lengths) {
  mapping_stats_t ms;
  init_mapping_stats(&ms);
  sequence_t seq;
  sequence_init(&seq,bisulfite_disabled,NULL);
  const char* reads[] = { "ACGTACGTAC", "AAAAANNNNN", "GGGCCCTTTAAA", "ACGTACGTAC" };
  uint64_t want[5] = {0,0,0,0,0};
  const char* const order = "NACGT";
  int i;
  for (i=0;i<4;++i) {
    string_set_buffer(&seq.read,(char*)reads[i],strlen(reads[i]));
    update_counts(&seq,&ms,i==2 ? 1 : 0);
    const char* p;
    if (i!=2) for (p=reads[i];*p;++p) ++want[strchr(order,*p)-order];
  }
  for (i=0;i<5;++i) GT_CHECK_MSG(ms.base_counts[0][0][i]==want[i],"end 1 base %c: %" PRIu64 " != %" PRIu64,
      order[i],ms.base_counts[0][0][i],want[i]);
  GT_CHECK_EQ_UINT(ms.base_counts[0][1][1],3);
  uint64_t* const n10 = ihash_get(ms.read_length_dist[0],10,uint64_t);
  uint64_t* const n12 = ihash_get(ms.read_length_dist[1],12,uint64_t);
  GT_CHECK(n10!=NULL && *n10==3);
  GT_CHECK(n12!=NULL && *n12==1);
  GT_CHECK(ihash_get(ms.read_length_dist[0],12,uint64_t)==NULL);
  sequence_destroy(&seq);
}

GT_TEST(read_type_by_contig_name) {
  char* controls[3] = { "NC_001422.1", "NC_001416.1", "NC_001604.1" };
  match_trace_t m;
  memset(&m,0,sizeof(m));
  char* names[] = { "chr1", "NC_001422.1", "NC_001416.1", "NC_001604.1", "x_NC_001416.1_y", "NC_00142" };
  const int want[] = { 0, 1, 2, 3, 2, 0 };
  int i;
  for (i=0;i<6;++i) {
    m.sequence_name = names[i];
    GT_CHECK_MSG(get_read_type(&m,controls)==want[i],"%s -> %d",names[i],get_read_type(&m,controls));
  }
}

static void fill_ihash(ihash_t* const h,gt_rng_t* const rng,uint64_t* const ref,const int64_t offset) {
  int i;
  for (i=0;i<30;++i) {
    const int64_t k = (int64_t)gt_rng_below(rng,60);
    const uint64_t c = 1 + gt_rng_below(rng,9);
    uint64_t* e = ihash_get(h,k+offset,uint64_t);
    if (e==NULL) {
      e = mm_alloc(uint64_t);
      *e = 0;
      ihash_insert(h,k+offset,e);
    }
    *e += c;
    ref[k] += c;
  }
}
GT_TEST(merge_mapping_stats_sums_threads) {
  enum { T = 4 };
  mapping_stats_t per[T], global;
  uint64_t reads[2][4] = {{0}}, bs[2][2] = {{0}}, un[2] = {0}, pairs = 0, mapq[256] = {0};
  uint64_t bases[7][2][5];
  uint64_t rl[2][60], dd[2][60], ins[60];
  memset(bases,0,sizeof(bases)); memset(rl,0,sizeof(rl)); memset(dd,0,sizeof(dd)); memset(ins,0,sizeof(ins));
  gt_rng_t rng = gt_rng_new(41);
  int t, e, j, k;
  for (t=0;t<T;++t) {
    init_mapping_stats(per+t);
    per[t].correct_pairs = gt_rng_below(&rng,1000); pairs += per[t].correct_pairs;
    for (e=0;e<2;++e) {
      for (j=0;j<4;++j) { per[t].reads[e][j] = gt_rng_below(&rng,1000); reads[e][j] += per[t].reads[e][j]; }
      for (j=0;j<2;++j) { per[t].BSreads[e][j] = gt_rng_below(&rng,1000); bs[e][j] += per[t].BSreads[e][j]; }
      per[t].unmapped[e] = gt_rng_below(&rng,1000); un[e] += per[t].unmapped[e];
      for (k=0;k<7;++k) for (j=0;j<5;++j) {
        per[t].base_counts[k][e][j] = gt_rng_below(&rng,100000); bases[k][e][j] += per[t].base_counts[k][e][j];
      }
      fill_ihash(per[t].read_length_dist[e],&rng,rl[e],0);
      fill_ihash(per[t].distance_dist[e],&rng,dd[e],0);
    }
    for (j=0;j<256;++j) { per[t].hist_mapq[j] = gt_rng_below(&rng,50); mapq[j] += per[t].hist_mapq[j]; }
    fill_ihash(per[t].insert_size_dist,&rng,ins,200);
  }
  merge_mapping_stats(&global,per,T);
  GT_CHECK_EQ_UINT(global.correct_pairs,pairs);
  for (e=0;e<2;++e) {
    for (j=0;j<4;++j) GT_CHECK_EQ_UINT(global.reads[e][j],reads[e][j]);
    for (j=0;j<2;++j) GT_CHECK_EQ_UINT(global.BSreads[e][j],bs[e][j]);
    GT_CHECK_EQ_UINT(global.unmapped[e],un[e]);
    for (k=0;k<7;++k) for (j=0;j<5;++j) GT_CHECK_EQ_UINT(global.base_counts[k][e][j],bases[k][e][j]);
    for (k=0;k<60;++k) {
      uint64_t* const a = ihash_get(global.read_length_dist[e],k,uint64_t);
      uint64_t* const b = ihash_get(global.distance_dist[e],k,uint64_t);
      GT_CHECK_MSG(rl[e][k] ? (a!=NULL && *a==rl[e][k]) : a==NULL,"read length %d end %d",k,e);
      GT_CHECK_MSG(dd[e][k] ? (b!=NULL && *b==dd[e][k]) : b==NULL,"distance %d end %d",k,e);
    }
  }
  for (j=0;j<256;++j) GT_CHECK_EQ_UINT(global.hist_mapq[j],mapq[j]);
  for (k=0;k<60;++k) {
    uint64_t* const a = ihash_get(global.insert_size_dist,k+200,uint64_t);
    GT_CHECK_MSG(ins[k] ? (a!=NULL && *a==ins[k]) : a==NULL,"insert %d",k+200);
  }
  // Per-thread tables are untouched
  uint64_t sum = 0;
  for (t=0;t<T;++t) sum += per[t].unmapped[0];
  GT_CHECK_EQ_UINT(sum,un[0]);
}

static char* print_to_string(void (*fn)(FILE*,void*),void* arg) {
  FILE* const f = tmpfile();
  fn(f,arg);
  const long n = ftell(f);
  rewind(f);
  char* const out = calloc((size_t)n+1,1);
  if (fread(out,1,(size_t)n,f) != (size_t)n) out[0] = 0;
  fclose(f);
  return out;
}
static uint64_t json_values[20];
static int json_n;
static void print_array(FILE* f,void* last) { output_json_uint_array(f,"HistMapq",json_values,json_n,1,*(bool*)last); }
static void print_element(FILE* f,void* last) { output_json_uint_element(f,"CorrectPairs",12345,2,*(bool*)last); }

GT_TEST(json_printers) {
  bool last = false, final = true;
  char* s = print_to_string(print_element,&last);
  GT_CHECK_EQ_STR(s,"\t\t\"CorrectPairs\": 12345,\n");
  free(s);
  s = print_to_string(print_element,&final);
  GT_CHECK_EQ_STR(s,"\t\t\"CorrectPairs\": 12345\n");
  free(s);
  int i;
  for (i=0;i<20;++i) json_values[i] = (uint64_t)i*i*1000;
  // Short arrays: one line; long arrays: eight values per line
  for (json_n=1;json_n<=20;++json_n) {
    char want[1024];
    int n = 0;
    for (i=0;i<json_n;++i) n += sprintf(want+n,"%s%" PRIu64,i ? "," : "",json_values[i]);
    s = print_to_string(print_array,&last);
    char flat[1024];
    int m = 0;
    const char* p;
    for (p=s;*p;++p) if (*p!=' ' && *p!='\t' && *p!='\n') flat[m++] = *p;
    flat[m] = 0;
    char expect[1100];
    sprintf(expect,"\"HistMapq\":[%s],",want);
    GT_CHECK_MSG(strcmp(flat,expect)==0,"n=%d: %s",json_n,s);
    int lines = 0;
    for (p=s;*p;++p) lines += *p=='\n';
    GT_CHECK_MSG(lines == (json_n <= 8 ? 1 : 2 + (json_n+7)/8),"n=%d: %d lines",json_n,lines);
    free(s);
  }
}

/*
 * Timers
 */
static void sleep_ms(const long ms) {
  struct timespec ts = { ms/1000, (ms%1000)*1000000L };
  while (nanosleep(&ts,&ts)!=0) {}
}
static uint64_t now_ns(void) {
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME,&ts);
  return (uint64_t)ts.tv_sec*1000000000ull + (uint64_t)ts.tv_nsec;
}
GT_TEST(timers_measure_elapsed_time) {
  gem_timer_t t;
  TIMER_RESET(&t);
  GT_CHECK_EQ_UINT(TIMER_GET_TOTAL_NS(&t),0);
  uint64_t w0 = now_ns();
  TIMER_START(&t);
  sleep_ms(20);
  TIMER_STOP(&t);
  uint64_t w1 = now_ns();
  const uint64_t first = TIMER_GET_TOTAL_NS(&t);
  GT_CHECK_MSG(first >= 20000000ull && first <= w1-w0,"first %" PRIu64 " ns",first);
  GT_CHECK_EQ_UINT(TIMER_GET_NUM_SAMPLES(&t),1);
  // Pause/continue: the 30 ms paused are not counted
  w0 = now_ns();
  TIMER_START(&t);
  sleep_ms(10);
  TIMER_PAUSE(&t);
  sleep_ms(30);
  TIMER_CONTINUE(&t);
  sleep_ms(10);
  TIMER_STOP(&t);
  w1 = now_ns();
  const uint64_t second = TIMER_GET_TOTAL_NS(&t) - first;
  GT_CHECK_MSG(second >= 20000000ull && second + 30000000ull <= w1-w0,"second %" PRIu64 " ns of %" PRIu64,second,w1-w0);
  GT_CHECK_EQ_UINT(TIMER_GET_NUM_SAMPLES(&t),2);
  // Restart forgets earlier samples
  TIMER_RESTART(&t);
  sleep_ms(5);
  TIMER_STOP(&t);
  GT_CHECK_EQ_UINT(TIMER_GET_NUM_SAMPLES(&t),1);
  GT_CHECK(TIMER_GET_TOTAL_NS(&t) >= 5000000ull);
  TIMER_RESET(&t);
  GT_CHECK_EQ_UINT(TIMER_GET_TOTAL_NS(&t),0);
  GT_CHECK_EQ_UINT(TIMER_GET_NUM_SAMPLES(&t),0);
}

GT_CASES(
  GT_CASE(counter_mean_variance_min_max),
  GT_CASE(template_length_estimate),
  GT_CASE(update_counts_counts_bases_and_lengths),
  GT_CASE(read_type_by_contig_name),
  GT_CASE(merge_mapping_stats_sums_threads),
  GT_CASE(json_printers),
  GT_CASE(timers_measure_elapsed_time)
)
