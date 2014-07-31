/*
 * PROJECT: GEMMapper
 * FILE: errors.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Error handling module
 */

#include "errors.h"

/*
 * StackTrace Printer
 */
#ifdef GEM_DEBUG
#include <execinfo.h>
#include <stdlib.h>

#define STACK_TRACE_SIZE 15

void gem_print_stack_trace() {
  void *stack[STACK_TRACE_SIZE];
  size_t size = backtrace(stack,STACK_TRACE_SIZE);
  fprintf(stderr,"<GT::StackTrace>\n");
  backtrace_symbols_fd(stack,size,2);
  fprintf(stderr,"<<<\n");
}

#else
void gem_print_stack_trace() {}
#endif

/*
 * Error Signal Handler
 */
void gem_error_signal_handler(int signal) {
  FILE* const error_stream = gem_error_get_stream();
  fprintf(error_stream,">>GEM:System.Error::Signal raised (no=%d) [%s (errno=%d)]\n",signal,strerror(errno),errno);
  fflush(error_stream);
  gem_print_stack_trace();
  exit(1);
}
void gem_handle_error_signals() {
  struct sigaction error_signal_handler;
  error_signal_handler.sa_handler = gem_error_signal_handler;
  sigemptyset(&error_signal_handler.sa_mask);
  error_signal_handler.sa_flags = 0;
  sigaction(SIGBUS,&error_signal_handler,NULL);
  sigaction(SIGSEGV,&error_signal_handler,NULL);
}

/*
 * Print ErrNo
 */
void gem_perror() {
  FILE* const error_stream = gem_error_get_stream();
  fprintf(error_stream,"> System.Error::%s (errno=%d)\n",strerror(errno),errno);
  fflush(error_stream);
}

