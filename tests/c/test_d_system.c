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
 *   Area D: runtime services.
 *     - thread ids (gem_thread_register_id / gem_thread_get_thread_id) from
 *       several threads
 *     - the error-signal handler installed by gruntime_init: a SIGSEGV in a
 *       child process runs the report function, prints the signal and exits 1
 *     - log streams, tabulated printing and the progress ticker (report.c)
 *     - system information (commons.c) against the C library
 *
 *   Internal signatures called (a speed iteration may change them; update
 *   this driver with them):
 *     void gem_thread_register_id(uint64_t); int64_t gem_thread_get_thread_id(void)
 *     void gem_handle_error_signals(void) (through gruntime_init)
 *     gem_error_get/set_report_function, gem_error_get_stream, gem_log_get/set_stream,
 *       gem_is_mute_error_stream, gem_is_mute_log_stream
 *     int tfprintf(FILE*,const char*,...), tab_fprintf(FILE*,...), fprintf_tabs(FILE*,int),
 *       tab_global_inc/dec/print
 *     ticker_count_reset(ticker_t*,bool,const char*,uint64_t top,uint64_t each,bool timed),
 *       ticker_add_process_label/ticker_add_finish_label, ticker_mutex_enable/cleanup,
 *       ticker_update_mutex, ticker_finish (ticker_count_finish, ticker_percentage_finish)
 *     system_get_num_processors/cwd/hostname/user_name, system_print_info(FILE*),
 *       system_get_time(struct timespec*)
 */

#include "gem3_test.h"
#include <pthread.h>
#include <pwd.h>
#include <signal.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
#include "system/commons.h"
#include "system/report.h"
#include "system/gthread.h"
#include "system/errors.h"

static char* file_text(FILE* const f) {
  const long n = ftell(f);
  rewind(f);
  char* const out = calloc((size_t)n+1,1);
  if (fread(out,1,(size_t)n,f) != (size_t)n) out[0] = 0;
  return out;
}

/*
 * Thread ids
 */
typedef struct { uint64_t id; int64_t seen[50]; } thread_arg_t;
// All threads stay alive until every one has registered, so no thread
// reuses the pthread_t of one that already exited.
static pthread_mutex_t live_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t live_cond = PTHREAD_COND_INITIALIZER;
static int live_count = 0, live_total = 0;
static void* register_and_read(void* p) {
  thread_arg_t* const a = p;
  a->seen[0] = gem_thread_get_thread_id(); // not registered yet
  gem_thread_register_id(a->id);
  pthread_mutex_lock(&live_mutex);
  if (++live_count == live_total) pthread_cond_broadcast(&live_cond);
  while (live_count < live_total) pthread_cond_wait(&live_cond,&live_mutex);
  pthread_mutex_unlock(&live_mutex);
  int i;
  for (i=1;i<50;++i) {
    a->seen[i] = gem_thread_get_thread_id();
    sched_yield();
  }
  return NULL;
}
GT_TEST(thread_ids_are_per_thread) {
  GT_CHECK_EQ_INT(gem_thread_get_thread_id(),0); // master, registered by gruntime_init
  enum { T = 8 };
  live_total = T;
  pthread_t th[T];
  thread_arg_t args[T];
  int i, k;
  for (i=0;i<T;++i) {
    args[i].id = (uint64_t)(i+1);
    GT_REQUIRE(pthread_create(th+i,NULL,register_and_read,args+i)==0);
  }
  for (i=0;i<T;++i) pthread_join(th[i],NULL);
  for (i=0;i<T;++i) {
    GT_CHECK_EQ_INT(args[i].seen[0],-1);
    for (k=1;k<50;++k) GT_CHECK_EQ_INT(args[i].seen[k],i+1);
  }
  GT_CHECK_EQ_INT(gem_thread_get_thread_id(),0);
  gem_thread_register_id(0);
  GT_CHECK_EQ_INT(gem_thread_get_thread_id(),0);
}

/*
 * Error signals
 */
static void custom_report(FILE* stream) {
  fprintf(stream,"CUSTOM-REPORT-MARK\n");
}
GT_TEST(error_signal_handler_reports_and_exits_1) {
  GT_CHECK(gem_error_get_report_function()!=NULL);
  report_function_t const before = gem_error_get_report_function();
  gem_error_set_report_function(custom_report);
  GT_CHECK(gem_error_get_report_function()==custom_report);
  GT_CHECK(gem_error_get_stream()==stderr);
  int fds[2];
  GT_REQUIRE(pipe(fds)==0);
  fflush(stdout);
  const pid_t pid = fork();
  GT_REQUIRE(pid >= 0);
  if (pid==0) {
    // Child: stderr into the pipe, then fault
    close(fds[0]);
    dup2(fds[1],2);
    raise(SIGSEGV);
    _exit(99); // not reached
  }
  close(fds[1]);
  char buf[4096];
  size_t n = 0;
  ssize_t r;
  while ((r=read(fds[0],buf+n,sizeof(buf)-1-n)) > 0) n += (size_t)r;
  buf[n] = 0;
  close(fds[0]);
  int status = 0;
  waitpid(pid,&status,0);
  gem_error_set_report_function(before);
  GT_CHECK_MSG(WIFEXITED(status) && WEXITSTATUS(status)==1,"child status 0x%x",status);
  GT_CHECK_MSG(strstr(buf,"CUSTOM-REPORT-MARK\n")!=NULL,"stderr: %s",buf);
  GT_CHECK_MSG(strstr(buf,">> GEM.System.Error::Signal raised (no=11)")!=NULL,"stderr: %s",buf);
}

/*
 * Log output
 */
static bool timestamp_prefix(const char* s,const char** rest) {
  // "%4d/%d/%d %02d:%02d:%02d -- " then the message
  const char* const sep = strstr(s," -- ");
  if (sep==NULL || sep-s > 19 || sep-s < 14) return false;
  int y, mo, d, h, mi, se;
  if (sscanf(s,"%4d/%d/%d %2d:%2d:%2d",&y,&mo,&d,&h,&mi,&se)!=6) return false;
  if (sep[-3]!=':' || sep[-6]!=':') return false;
  *rest = sep+4;
  return y>=2000 && mo>=1 && mo<=12 && d>=1 && d<=31 && h<24 && mi<60 && se<61;
}
GT_TEST(log_streams_and_timestamped_print) {
  GT_CHECK(gem_log_get_stream()==stderr);
  GT_CHECK(!gem_is_mute_error_stream() && !gem_is_mute_log_stream());
  FILE* const f = tmpfile();
  gem_log_set_stream(f);
  GT_CHECK(gem_log_get_stream()==f);
  const int printed = tfprintf(gem_log_get_stream(),"[%s %d]\n","value",42);
  gem_log_set_stream(NULL);
  GT_CHECK(gem_log_get_stream()==stderr);
  char* const text = file_text(f);
  fclose(f);
  const char* rest = NULL;
  GT_CHECK_MSG(timestamp_prefix(text,&rest) && strcmp(rest,"[value 42]\n")==0,"%s",text);
  GT_CHECK_EQ_INT(printed,(int)strlen(text));
  free(text);
}

GT_TEST(tabulated_printing) {
  FILE* const f = tmpfile();
  tab_fprintf(f,"a\n");
  tab_global_inc();
  tab_fprintf(f,"b%d\n",1);
  tab_global_inc();
  tab_global_print(f);
  fputs("c\n",f);
  tab_global_dec();
  tab_global_dec();
  tab_global_dec(); // never below zero
  tab_fprintf(f,"d\n");
  fprintf_tabs(f,3);
  fputs("|\n",f);
  fprintf_tabs(f,25);
  fputs("|\n",f);
  char* const text = file_text(f);
  fclose(f);
  GT_CHECK_EQ_STR(text,"a\n  b1\n    c\nd\n   |\n                         |\n");
  free(text);
}

static char* strip_timestamps(char* const text) {
  // Remove the "YYYY/M/D HH:MM:SS -- " prefix of every line
  char* const out = calloc(strlen(text)+1,1);
  char* line = text;
  size_t n = 0;
  while (*line) {
    char* const eol = strchr(line,'\n');
    const char* rest = line;
    if (eol) *eol = 0;
    if (!timestamp_prefix(line,&rest)) rest = line;
    n += (size_t)sprintf(out+n,"%s\n",rest);
    if (!eol) break;
    line = eol+1;
  }
  return out;
}
GT_TEST(count_ticker) {
  // As the mapper drives it: a step of 10, updates from several threads
  FILE* const f = tmpfile();
  gem_log_set_stream(f);
  ticker_t t;
  ticker_count_reset(&t,true,"SE::Mapping Sequences",0,10,false);
  ticker_add_process_label(&t,"#","sequences processed");
  ticker_add_finish_label(&t,"Total","sequences processed");
  ticker_mutex_enable(&t);
  ticker_update_mutex(&t,10);
  ticker_update_mutex(&t,4);
  ticker_update_mutex(&t,10);
  ticker_update_mutex(&t,3);
  ticker_finish(&t);
  ticker_finish(&t); // a second finish prints nothing
  ticker_update_mutex(&t,10); // after finishing, nothing
  ticker_mutex_cleanup(&t);
  // A disabled ticker prints nothing
  ticker_t q;
  ticker_count_reset(&q,false,"quiet",0,1,false);
  ticker_finish(&q);
  gem_log_set_stream(NULL);
  char* const text = file_text(f);
  fclose(f);
  char* const plain = strip_timestamps(text);
  GT_CHECK_EQ_STR(plain,
      "[SE::Mapping Sequences]\n"
      "# 10 sequences processed\n"
      "# 20 sequences processed\n"
      "Total 27 sequences processed... done\n");
  free(plain);
  free(text);
}

GT_TEST(percentage_ticker_finishes_at_max) {
  FILE* const f = tmpfile();
  gem_log_set_stream(f);
  ticker_t t;
  ticker_percentage_reset(&t,true,"Work",200,25,false);
  ticker_add_process_label(&t,"","done");
  ticker_add_finish_label(&t,"All","items");
  ticker_update(&t,50);
  ticker_update(&t,100);
  ticker_update(&t,50); // reaches max: finishes
  ticker_finish(&t);
  gem_log_set_stream(NULL);
  char* const text = file_text(f);
  fclose(f);
  char* const plain = strip_timestamps(text);
  GT_CHECK_EQ_STR(plain,"[Work]\n  25% done\n  75% done\nAll 100% items... done\n");
  free(plain);
  free(text);
}

/*
 * System information
 */
GT_TEST(system_information_matches_libc) {
  GT_CHECK_EQ_UINT(system_get_num_processors(),(uint64_t)sysconf(_SC_NPROCESSORS_ONLN));
  char want[4096];
  GT_REQUIRE(getcwd(want,sizeof(want))!=NULL);
  char* const cwd = system_get_cwd();
  GT_CHECK_EQ_STR(cwd,want);
  char host[1100];
  GT_REQUIRE(gethostname(host,1024)==0);
  char* const hn = system_get_hostname();
  GT_CHECK_EQ_STR(hn,host);
  struct passwd* const pw = getpwuid(geteuid());
  GT_CHECK_EQ_STR(system_get_user_name(),pw ? pw->pw_name : "");
  FILE* const f = tmpfile();
  system_print_info(f);
  char* const text = file_text(f);
  fclose(f);
  char expect[8192];
  snprintf(expect,sizeof(expect),"  => CWD %s\n  => Hostname %s\n  => Username %s\n  => CPUs %" PRIu64 "\n",
      want,host,pw ? pw->pw_name : "",system_get_num_processors());
  GT_CHECK_MSG(strncmp(text,"[GEM]>System.Info\n  => Date ",28)==0,"%s",text);
  const char* const tail = strstr(text," CET\n");
  GT_CHECK_MSG(tail!=NULL && strcmp(tail+5,expect)==0,"%s",text);
  free(text);
  free(cwd);
  free(hn);
  struct timespec ts;
  const time_t before = time(NULL);
  system_get_time(&ts);
  const time_t after = time(NULL);
  GT_CHECK(ts.tv_sec >= before && ts.tv_sec <= after && ts.tv_nsec >= 0 && ts.tv_nsec < 1000000000L);
}

GT_CASES(
  GT_CASE(thread_ids_are_per_thread),
  GT_CASE(error_signal_handler_reports_and_exits_1),
  GT_CASE(log_streams_and_timestamped_print),
  GT_CASE(tabulated_printing),
  GT_CASE(count_ticker),
  GT_CASE(percentage_ticker_finishes_at_max),
  GT_CASE(system_information_matches_libc)
)
