/*
 * PROJECT: GEMMapper
 * FILE: errors.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Error handling module
 */

#include "system/errors.h"
#include "system/gthread.h"

/*
 * StackTrace Printer
 */
#ifdef GEM_DEBUG
#include <execinfo.h>
#include <stdlib.h>
#define STACK_TRACE_SIZE 15
pthread_mutex_t gem_print_stack_trace_mutex = PTHREAD_MUTEX_INITIALIZER;
void gem_print_stack_trace() {
  pthread_mutex_lock(&gem_print_stack_trace_mutex);
  {
    void *stack[STACK_TRACE_SIZE];
    size_t size = backtrace(stack,STACK_TRACE_SIZE);
    fprintf(stderr,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
    fprintf(stderr,"<GEM::StackTrace>\n");
    backtrace_symbols_fd(stack,size,2);
    fprintf(stderr,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
  }
  pthread_mutex_unlock(&gem_print_stack_trace_mutex);
}
#else
void gem_print_stack_trace() {}
#endif

/*
 * Error Signal Handler
 */
pthread_mutex_t gem_error_signal_handler_mutex = PTHREAD_MUTEX_INITIALIZER;
void gem_error_signal_handler(int signal) {
  pthread_mutex_lock(&gem_error_signal_handler_mutex);
    static bool reported = false;
    if (!reported) {
      FILE* const error_stream = gem_error_get_stream();
      // Print stack trace
      gem_print_stack_trace();
      // Print custom error function
      report_function_t report_function = gem_error_get_report_function();
      if (report_function!=NULL) report_function(error_stream);
      // Print signal label
      fprintf(stderr,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
      fprintf(error_stream,">> GEM.System.Error::Signal raised (no=%d) [errno=%d,%s]\n",signal,errno,strerror(errno));
      fprintf(stderr,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
      fflush(error_stream);
    }
  exit(1);
  // pthread_mutex_unlock(&gem_error_signal_handler_mutex);
}
void gem_handle_error_signals() {
  struct sigaction error_signal_handler;
  error_signal_handler.sa_handler = gem_error_signal_handler;
  sigemptyset(&error_signal_handler.sa_mask);
  error_signal_handler.sa_flags = 0;
  sigaction(SIGBUS,&error_signal_handler,NULL);
  sigaction(SIGSEGV,&error_signal_handler,NULL);
  sigaction(SIGFPE,&error_signal_handler,NULL);
  sigaction(SIGILL,&error_signal_handler,NULL);
}
/*
 * Print ErrNo
 */
void gem_perror() {
  FILE* const error_stream = gem_error_get_stream();
  fprintf(error_stream,">> GEM.System.Error::%s (errno=%d)\n",strerror(errno),errno);
  fflush(error_stream);
}

