/*
 * PROJECT: GEM-Tools library
 * FILE: gt_error.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_error.h"
#include "gt_commons.h"

// Gem-tools error/output streams
FILE* gt_error_stream=NULL;
FILE* gt_log_stream=NULL;
FILE* gt_debug_stream=NULL;
bool mute_error_stream = false;
bool mute_report_stream = false;

/*
 * Getters/Setters ELD-streams
 */
inline FILE* gt_error_get_stream() {
  return (gt_error_stream!=NULL) ? gt_error_stream : stderr;
}
inline void gt_error_set_stream(FILE* const stream) {
  gt_error_stream = stream;
}
inline FILE* gt_log_get_stream() {
  return (gt_log_stream!=NULL) ? gt_log_stream : stderr;
}
inline void gt_log_set_stream(FILE* const stream) {
  gt_log_stream = stream;
}
inline FILE* gt_debug_get_stream() {
  return (gt_debug_stream!=NULL) ? gt_debug_stream : stderr;
}
inline void gt_debug_set_stream(FILE* const stream) {
  gt_debug_stream = stream;
}
/*
 * Mute/Articulate ELD-streams
 */
inline void gt_mute_error_stream() {mute_error_stream=true;}
inline void gt_mute_report_stream() {mute_error_stream=true;}
inline void gt_articulate_error_stream() {mute_error_stream=false;}
inline void gt_articulate_report_stream() {mute_error_stream=false;}

/*
 * StackTrace Printer
 */
#ifdef GT_DEBUG
#include <execinfo.h>
#include <stdlib.h>

#define GT_STACK_TRACE_SIZE 15

void gt_print_stack_trace() {
  void *stack[GT_STACK_TRACE_SIZE];
  size_t size = backtrace(stack,GT_STACK_TRACE_SIZE);
  fprintf(stderr,"<GT::StackTrace>\n");
  backtrace_symbols_fd(stack,size,2);
  fprintf(stderr,"<<<\n");
}

#else
inline void gt_print_stack_trace() {}
#endif

/*
 * Error Signal Handler
 */
void gt_error_signal_handler(int signal) {
  FILE* const error_stream=gt_error_get_stream();
  fprintf(error_stream,"> System.Error::Signal raised (no=%d)\n",signal);
  fflush(error_stream);
  gt_print_stack_trace();
  exit(1);
}
void gt_handle_error_signals() {
  struct sigaction error_signal_handler;
  error_signal_handler.sa_handler = gt_error_signal_handler;
  sigemptyset(&error_signal_handler.sa_mask);
  error_signal_handler.sa_flags = 0;
  sigaction(SIGBUS,&error_signal_handler,NULL);
  sigaction(SIGSEGV,&error_signal_handler,NULL);
}

/*
 * Print ErrNo
 */
void gt_perror() {
  FILE* const error_stream=gt_error_get_stream();
  fprintf(error_stream,"> System.Error::%s (errno=%d)\n",strerror(errno),errno);
  fflush(error_stream);
}

// Text weekDays & months
char* gt_label_weekdays[] = {"Sun","Mon","Tue","Wed","Thu","Fri","Sat"};
char* gt_label_months[] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};

/*
 * Time Printed Formated functions
 */
gt_status gt_vtfprintf(FILE* stream,const char* format,va_list v_args) {
  gt_status chars_printed = 0; // Base printed chars
  // Get Current Time
  time_t current_time=time(0);
  struct tm local_time;
  localtime_r(&current_time,&local_time);
  // Print Current year
  chars_printed+=fprintf(stream,"%4d/%d/%d ",1900+local_time.tm_year,local_time.tm_mon+1,local_time.tm_mday);
  // Print Current Time
  chars_printed+=fprintf(stream,"%02d:%02d:%02d -- ",local_time.tm_hour,local_time.tm_min,local_time.tm_sec);
  chars_printed+=vfprintf(stream,format,v_args);
  fflush(stream);
  //  // Print Current Day
  //  if (GT_BETWEEN(local_time.tm_wday,0,6)) {
  //    chars_printed+=fprintf(stream,"%s:",gt_label_weekdays[local_time.tm_wday]);
  //  } else {
  //    gt_fatal_error(TPRINTF);
  //  }
  //  // Print Current Month
  //  if (GT_BETWEEN(local_time.tm_mon,0,11)) {
  //    chars_printed+=fprintf(stream,"%s ",gt_label_months[local_time.tm_wday]);
  //  } else {
  //    gt_fatal_error(TPRINTF);
  //  }
  return chars_printed;
}
gt_status gt_tfprintf(FILE* stream,const char* format,...) {
  GT_NULL_CHECK(stream);
  GT_NULL_CHECK(format);
  va_list v_args;
  va_start(v_args,format);
  const gt_status chars_printed = gt_vtfprintf(stream,format,v_args);
  va_end(v_args);
  return chars_printed;
}
gt_status gt_vtprintf(const char* format,va_list v_args) {
  GT_NULL_CHECK(format);
  return gt_vtfprintf(stdout,format,v_args);
}
gt_status gt_tprintf(const char* format,...) {
  GT_NULL_CHECK(format);
  va_list v_args;
  va_start(v_args,format);
  const gt_status chars_printed = gt_vtfprintf(stdout,format,v_args);
  va_end(v_args);
  return chars_printed;
}

