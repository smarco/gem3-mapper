/*
 * PROJECT: GEMMapper
 * FILE: report.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Reporting module
 */

#ifndef REPORT_H_
#define REPORT_H_

#include "system/commons.h"

// Labels
#define GEM_LABEL_ERROR "GEM::Error"
#define GEM_LABEL_FATAL_ERROR "GEM::FatalError"
#define GEM_LABEL_DEBUG "GEM::Debug"
#define GEM_LABEL_LOG "GEM::Log"
#define GEM_LABEL_WARNING "GEM::Warning"

// Base-name of the sources
#define GEM_ERROR_BASENAME(S) \
  ({ const char* const slash=strrchr((S),'/'); \
     slash ? slash + 1 : (S); })

// Getters/Setters ELD-function
typedef void (*report_function_t)(FILE*);
report_function_t gem_error_get_report_function();
void gem_error_set_report_function(report_function_t report_function);

// Getters/Setters ELD-streams
#define GEM_DEFAULT_REPORT_STREAM stderr
FILE* gem_error_get_stream();
void gem_error_set_stream(FILE* const stream);
FILE* gem_log_get_stream();
void gem_log_set_stream(FILE* const stream);
FILE* gem_info_get_stream();
void gem_info_set_stream(FILE* const stream);
FILE* gem_debug_get_stream();
void gem_debug_set_stream(FILE* const stream);
// Mute/Articulate ELD-streams
void gem_mute_error_stream();
void gem_mute_report_stream();
void gem_articulate_error_stream();
void gem_articulate_report_stream();
bool gem_is_mute_error_stream();
bool gem_is_mute_report_stream();
/*
 * Low-level code report functions
 *
 * E.g. EXITING ERROR (FATAL)
 *  gem_report_error_begin_block(args...) {
 *    ..Code Block..
 *  } gem_report_end_block(args...);
 *
 * E.g. NOT-EXITING
 *  gem_report_error_begin_block(args...) {
 *    ..Code Block..
 *  } gem_report_end_block(args...);
 */
#define gem_report_error_begin_block(gem_label,gem_error_name,args...) \
  do { \
    FILE* const gem_stream=gem_error_get_stream(); \
    if (!gem_is_mute_report_stream()) { \
      fprintf(gem_stream,gem_label" (%s:%d,%s)\n "GEM_ERROR_##gem_error_name"\n", \
        GEM_ERROR_BASENAME(__FILE__),__LINE__,__func__, ##args); \
      fflush(gem_stream); \
    }
#define gem_report_begin_block(gem_label,stream,gem_report_msg,args...) \
  do { \
    FILE* const gem_stream=gem_##stream##_get_stream(); \
    if (!gem_is_mute_report_stream()) { \
      fprintf(gem_stream,gem_label" (%s:%d,%s)\n "gem_report_msg"\n", \
        GEM_ERROR_BASENAME(__FILE__),__LINE__,__func__, ##args); \
      fflush(gem_stream); \
    }
#define gem_report_raw_begin_block(gem_label,stream,gem_report_msg,args...) \
  do { \
    FILE* const gem_stream=gem_##stream##_get_stream(); \
    if (!gem_is_mute_report_stream()) { \
      fprintf(gem_stream,gem_label":: "gem_report_msg"\n",##args); \
      fflush(gem_stream); \
    }
#define gem_report__timestamp_begin_block(gem_label,stream,gem_report_msg,args...) \
  do { \
    FILE* const gem_stream=gem_##stream##_get_stream(); \
    if (!gem_is_mute_report_stream()) { \
      tfprintf(gem_stream,gem_label":: "gem_report_msg"\n",##args); \
      fflush(gem_stream); \
    }
// Gem error report end block (forces SEGFAULT so the IDE-debugger can capture the breakpoint -handy-)
#ifdef GEM_DEBUG_FORCE_SEGFAULT_AT_ERROR
  #define gem_report_end_block(print_gem_report,print_stack_trace,do_exit,exit_code) \
      if (print_gem_report) gem_error_get_report_function()(gem_stream); \
      if (print_stack_trace) gem_print_stack_trace(); \
      if (do_exit) { raise(SIGSEGV); exit(exit_code); } \
    } while (0)
#else
  #define gem_report_end_block(print_gem_report,print_stack_trace,do_exit,exit_code) \
      if (print_gem_report) gem_error_get_report_function()(gem_stream); \
      if (print_stack_trace) gem_print_stack_trace(); \
      if (do_exit) exit(exit_code); \
    } while (0)
#endif
/*
 * Log Utilities
 */
// Simple LOG
#define gem_slog(gem_log_msg,args...) \
  do { \
    FILE* const gem_stream=gem_log_get_stream(); \
    if (!gem_is_mute_report_stream()) { \
      fprintf(gem_stream,gem_log_msg,##args); \
      fflush(gem_stream); \
    } \
  } while (0)
// Time LOG
#define gem_xlog(gem_log_msg,args...) \
  gem_report__timestamp_begin_block(GEM_LABEL_LOG,log,gem_log_msg,##args) \
  gem_report_end_block(0,0,0,0)
#define gem_log(gem_log_msg,args...) \
  do { \
    FILE* const gem_stream=gem_log_get_stream(); \
    if (!gem_is_mute_report_stream()) { \
      tfprintf(gem_stream,gem_log_msg"\n",##args); \
      fflush(gem_stream); \
    } \
  } while (0)
// Conditional Time LOG
#define gem_cond_log(condition,gem_log_msg,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gem_log(gem_log_msg,##args); \
    } \
  } while (0)
// File LOG
#define gem_flog(stream,gem_log_msg,args...) \
  do { \
    fprintf(stream,gem_log_msg,##args); \
    fflush(stream); \
  } while (0)
// Conditional File LOG
#define gem_cond_flog(condition,stream,gem_log_msg,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gem_flog(stream,gem_log_msg,##args); \
    } \
  } while (0)
// Info LOG
#define gem_info(gem_info_msg,args...) \
  do { \
    FILE* const gem_stream=gem_info_get_stream(); \
    if (!gem_is_mute_report_stream()) { \
      fprintf(gem_stream,gem_info_msg,##args); \
      fflush(gem_stream); \
    } \
  } while (0)
/*
 * Printing functions
 */
#define PRINT_COLUMN(stream,position,column_size,separator) \
  if (position && ((position+1)%column_size==0)) fprintf(stream,"\n"); \
  else fprintf(stream,separator);
#define PRINT_COLUMN_CLOSE(stream,position,column_size) \
  if (position%column_size!=0) fprintf(stream,"\n");
// Time Printed Formated functions
int vtfprintf(FILE* stream,const char* format,va_list v_args);
int tfprintf(FILE* stream,const char* format,...);
int vtprintf(const char* format,va_list v_args);
int tprintf(const char* format,...);
// Tabulated Printed Formated functions
int tab_vfprintf(FILE* stream,const char* format,va_list v_args);
int tab_fprintf(FILE* stream,const char* format,...);
int tab_vprintf(const char* format,va_list v_args);
int tab_printf(const char* format,...);
/*
 * Tabulate Data
 */
void fprintf_tabs(FILE* const stream,const int num_spaces);
void tab_global_print(FILE* const stream);
void tab_global_inc();
void tab_global_add(const uint64_t amount);
void tab_global_dec();
void tab_global_subtract(const uint64_t amount);
void tab_global_reset();
/*
 * Ticker
 */
typedef enum { ticker_percentage, ticker_count } ticker_type_t;
typedef struct {
  // Type
  ticker_type_t ticker_type;
  // Limits
  uint64_t max;
  uint64_t global_ticks;
  uint64_t local_ticks;
  uint64_t step_ticks;
  // Timer
  bool timed;
  struct timespec begin_timer;
  struct timespec end_timer;
  // Labels
  char* process_begin;
  char* process_end;
  char* finish_begin;
  char* finish_end;
  // Status
  bool enabled;
  bool finished;
  // Mutex
  pthread_mutex_t mutex;
} ticker_t;
// Percentage based
void ticker_percentage_reset(
    ticker_t* const ticker,const bool enabled,const char* const label,
    const uint64_t max,const uint64_t step,const bool timed);
// Count based
void ticker_count_reset(
    ticker_t* const ticker,const bool enabled,const char* const label,
    const uint64_t top,const uint64_t each,const bool timed);
void ticker_update(ticker_t* const ticker,const uint64_t n);
void ticker_update_mutex(ticker_t* const ticker,const uint64_t n);
void ticker_finish(ticker_t* const ticker);
void ticker_finish_mutex(ticker_t* const ticker);
// Adding labels
void ticker_add_process_label(ticker_t* const ticker,char* const process_begin,char* const process_end);
void ticker_add_finish_label(ticker_t* const ticker,char* const finish_begin,char* const finish_end);
// Enable/Disable ticker
void ticker_set_status(ticker_t* const ticker,const bool enabled);
// Set granularity
void ticker_set_step(ticker_t* const ticker,const uint64_t step);
// Enable mutex
void ticker_mutex_enable(ticker_t* const ticker);
void ticker_mutex_cleanup(ticker_t* const ticker);
// Percentage/Count Ticker Macros
#define TICKER_COND_BEGIN(condition,maximum,msg,timed) \
  ticker_t __ticker; \
  ticker_percentage_reset(&__ticker,condition,msg,maximum,1,timed);
#define TICKER_BEGIN(maximum,msg,timed) TICKER_COND_BEGIN(true,maximum,msg,timed)
#define CTICKER_COND_BEGIN(condition,maximum,each,msg,timed) \
  ticker_t __ticker; \
  ticker_count_reset(&__ticker,condition,msg,maximum,each,timed);
#define CTICKER_BEGIN(maximum,each,msg,timed) CTICKER_COND_BEGIN(true,maximum,each,msg,timed)
// Tick Macros
#define TICK() ticker_update(&__ticker,1)
#define TICK_INC(increment) ticker_update(&__ticker,increment)
#define TICKER_END() ticker_finish(&__ticker)
/*
 * Print's template helpers
 */
uint64_t calculate_memory_required_v(const char *template,va_list v_args);
uint64_t calculate_memory_required_va(const char *template,...);

#endif /* REPORT_H_ */
