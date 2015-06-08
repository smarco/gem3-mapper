/*
 * PROJECT: GEM-Tools library
 * FILE: gt_error.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_ERROR_H_
#define GT_ERROR_H_

#include <stdio.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdarg.h>

// Internally to Gem-tools error codes are returned as gt_status
typedef int32_t gt_status;

// Codes gt_status
#define GT_STATUS_OK 1
#define GT_STATUS_EOF 0
#define GT_STATUS_FAIL -1

// Base-name of the sources
#define GT_ERROR_BASENAME(S) \
  ({ const char* const slash=strrchr((S),'/'); \
     slash ? slash + 1 : (S); })

// Gem-tools error/output ELD-streams
extern FILE* gt_error_stream;
extern FILE* gt_log_stream;
extern FILE* gt_debug_stream;
extern bool mute_error_stream;
extern bool mute_report_stream;
// Getters/Setters ELD-streams
#define GT_DEFAULT_REPORT_STREAM stderr
inline FILE* gt_error_get_stream();
inline void gt_error_set_stream(FILE* const stream);
inline FILE* gt_log_get_stream();
inline void gt_log_set_stream(FILE* const stream);
inline FILE* gt_debug_get_stream();
inline void gt_debug_set_stream(FILE* const stream);
// Mute/Articulate ELD-streams
inline void gt_mute_error_stream();
inline void gt_mute_report_stream();
inline void gt_articulate_error_stream();
inline void gt_articulate_report_stream();

// Labels
#define GT_LABEL_ERROR "Error"
#define GT_LABEL_FATAL_ERROR "Fatal error"
#define GT_LABEL_DEBUG "GTDebug"
#define GT_LABEL_LOG "GTLog"
#define GT_LABEL_WARNING "Warning"

/*
 * StackTrace Printer
 */
void gt_print_stack_trace();

/*
 * Error Signal Handler
 */
void gt_handle_error_signals();

/*
 * Print ErrNo
 */
void gt_perror();

/*
 * Low-level code report functions
 *
 * E.g. EXITING ERROR (FATAL)
 *  gt_report_error_begin_block(gt_label,gt_error_name,args...) {
 *    ..Code Block..
 *  } gt_report_end_block__exit(exit_error_code,true/false);
 *
 * E.g. NOT-EXITING
 *  gt_report_error_begin_block(gt_label,gt_report_msg,args...) {
 *    ..Code Block..
 *  } gt_report_end_block();
 */
#define gt_report_error_begin_block(gt_label,gt_error_name,args...) \
  do { \
    if (!mute_error_stream) { \
      FILE* const error_stream=gt_error_get_stream(); \
      fprintf(error_stream,gt_label" (%s:%d,%s)\n "GT_ERROR_##gt_error_name"\n", \
        GT_ERROR_BASENAME(__FILE__),__LINE__,__func__, ##args); \
      fflush(error_stream); \
    }
#define gt_report_begin_block(gt_label,stream,gt_report_msg,args...) \
  do { \
    if (!mute_report_stream) { \
      FILE* const gt_stream=gt_##stream##_get_stream(); \
      fprintf(gt_stream,gt_label" (%s:%d,%s)\n "gt_report_msg"\n", \
        GT_ERROR_BASENAME(__FILE__),__LINE__,__func__, ##args); \
      fflush(gt_stream); \
    }
#define gt_report_raw_begin_block(gt_label,stream,gt_report_msg,args...) \
  do { \
    if (!mute_report_stream) { \
      FILE* const gt_stream=gt_##stream##_get_stream(); \
      fprintf(gt_stream,gt_label":: "gt_report_msg"\n",##args); \
      fflush(gt_stream); \
    }
#define gt_report__timestamp_begin_block(gt_label,stream,gt_report_msg,args...) \
  do { \
    if (!mute_report_stream) { \
      FILE* const gt_stream=gt_##stream##_get_stream(); \
      gt_tfprintf(gt_stream,gt_label":: "gt_report_msg"\n",##args); \
      fflush(gt_stream); \
    }
#define gt_report_end_block() \
  } while (0)
#define gt_report_end_block__exit(exit_code,print_stack_trace) \
    if (print_stack_trace) gt_print_stack_trace(); \
    exit(exit_code); \
  } while (0)
/*
 * GT-Error handlers
 */
#define gt_fatal_error(gt_error_name,args...) \
  gt_report_error_begin_block(GT_LABEL_FATAL_ERROR,gt_error_name,##args) \
  gt_report_end_block__exit(1,1)
#define gt_fatal_error__perror(gt_error_name,args...) \
  gt_report_error_begin_block(GT_LABEL_FATAL_ERROR,gt_error_name,##args) \
  gt_perror(); \
  gt_report_end_block__exit(1,1)
#define gt_error(gt_error_name,args...) \
  gt_report_error_begin_block(GT_LABEL_ERROR,gt_error_name,##args) \
  gt_report_end_block()
/*
 * ErrorMsg handlers
 */
#define gt_fatal_error_msg(gt_error_msg,args...) \
  gt_report_begin_block(GT_LABEL_FATAL_ERROR,error,gt_error_msg,##args) \
  gt_report_end_block__exit(1,1)
#define gt_error_msg(gt_error_msg,args...) \
  gt_report_begin_block(GT_LABEL_ERROR,error,gt_error_msg,##args) \
  gt_report_end_block()
/*
 * GT-Warning handler
 */
#define gt_warn(gt_warning_name,args...) \
  gt_report_error_begin_block(GT_LABEL_WARNING,gt_warning_name,##args) \
  gt_report_end_block()
/*
 * GT-Exception handlers (conditional error handlers)
 */
#define gt_cond_fatal_error(condition,gt_error_name,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gt_fatal_error(gt_error_name,##args); \
    } \
  } while (0)
#define gt_cond_fatal_error__perror(condition,gt_error_name,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gt_fatal_error__perror(gt_error_name,##args); \
    } \
  } while (0)
#define gt_cond_error(condition,gt_error_name,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gt_error(gt_error_name,##args); \
    } \
  } while (0)
/*
 * Exception handlers (conditional error handlers)
 */
#define gt_cond_fatal_error_msg(condition,error_msg,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gt_fatal_error_msg(error_msg,##args); \
    } \
  } while (0)
#define gt_cond_error_msg(condition,error_msg,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gt_error_msg(error_msg,##args); \
    } \
  } while (0)
/*
 * Robust checkers
 */
#ifndef GT_NO_CONSISTENCY_CHECKS
  #define gt_fatal_check(condition,gt_error_name,args...) gt_cond_fatal_error(condition,gt_error_name,##args)
  #define gt_check(condition,gt_error_name,args...) gt_cond_error(condition,gt_error_name,##args)
  #define gt_check_block(condition) if (condition)
#else
  #define gt_fatal_check(condition,gt_error_name,args...)
  #define gt_check(condition,gt_error_name,args...)
  #define gt_check_block(condition) if (0)
#endif
/*
 * Debug Utilities
 */
#ifdef GT_DEBUG
  #define gt_debug(gt_debug_msg,args...) \
    gt_report_raw_begin_block(GT_LABEL_DEBUG,debug,gt_debug_msg,##args) \
    gt_report_end_block()
  #define gt_debug_msg(gt_debug_msg,args...) \
    gt_report_begin_block(GT_LABEL_DEBUG,debug,gt_debug_msg,##args) \
    gt_report_end_block()
  #define gt_cond_debug_msg(condition,gt_debug_msg,args...) \
    do { \
      if (__builtin_expect((condition),0)){ \
        gt_debug_msg(gt_debug_msg,##args); \
      } \
    } while (0)
  #define gt_debug_msg__stack(gt_debug_msg,args...) \
    gt_debug_msg(gt_debug_msg,args...); \
    gt_print_stack_trace();
  #define gt_debug_block(condition) if (condition)
#else
  #define gt_debug_msg(gt_debug_msg,args...)
  #define gt_debug_msg__stack(gt_debug_msg,args...)
  #define gt_cond_debug_msg(condition,gt_debug_msg,args...)
  #define gt_debug_block(condition) if (0)
#endif
/*
 * Log Utilities
 */
#define gt_log(gt_log_msg,args...) \
  gt_report__timestamp_begin_block(GT_LABEL_LOG,log,gt_log_msg,##args) \
  gt_report_end_block()
#define gt_slog(gt_log_msg,args...) \
  do { \
    if (!mute_report_stream) { \
      FILE* const gt_stream=gt_log_get_stream(); \
      fprintf(gt_stream,gt_log_msg,##args); \
      fflush(gt_stream); \
    } \
  } while (0)
#define gt_cond_log(condition,gt_log_msg,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gt_log(gt_log_msg,##args); \
    } \
  } while (0)

/*
 * Time Printed Formated functions
 */
gt_status gt_vtfprintf(FILE* stream,const char* format,va_list v_args);
gt_status gt_tfprintf(FILE* stream,const char* format,...);
gt_status gt_vtprintf(const char* format,va_list v_args);
gt_status gt_tprintf(const char* format,...);

/*
 * ERROR CODES/MSG
 *   #define GT_ERROR_<CODE> "<MSG>"
 */
// Library/Program errors
#define GT_ERROR_NOT_ZERO "Value Zero. Variable %s must be non-zero"
#define GT_ERROR_INVALID_VALUE "Invalid Value. Variable %s must be %s"
#define GT_ERROR_POSITION_OUT_OF_RANGE "Requested position (%"PRIu64") out of range [%"PRIu64",%"PRId64"]"
#define GT_ERROR_SELECTION_NOT_VALID "Library error. Selection not valid"
#define GT_ERROR_ALG_INCONSISNTENCY "Library error. Algorithmic inconsistency, check your program"
#define GT_ERROR_NOT_IMPLEMENTED "Function/Feature not implemented yet (sorry)"

#define GT_ERROR_NULL_HANDLER "Null handler %s "

// System errors
#define GT_ERROR_SYS_ERROR "System error signal raised"
#define GT_ERROR_SYS_MMAP "Mmap call error"
#define GT_ERROR_SYS_MMAP_FILE "Could not map file '%s' to memory"
#define GT_ERROR_SYS_UNMAP "Could not unmap memory"
#define GT_ERROR_SYS_THREAD "Could not create thread"
#define GT_ERROR_SYS_PIPE "Could not create pipe"
#define GT_ERROR_SYS_MUTEX "Mutex call error"
#define GT_ERROR_SYS_MUTEX_INIT "Mutex initialization error"
#define GT_ERROR_SYS_MUTEX_DESTROY "Mutex destroy call error"
#define GT_ERROR_SYS_COND_VAR "Conditional variable call error"
#define GT_ERROR_SYS_COND_VAR_INIT "Conditional variable initialization error"
#define GT_ERROR_SYS_COND_VAR_DESTROY "Conditional variable destroy call error"
#define GT_ERROR_SYS_MKSTEMP "Could not create a temporal file (mkstemp:'%s')"
#define GT_ERROR_SYS_HANDLE_TMP "Failed to handle temporal file"
#define GT_ERROR_SYS_SYSCONF "Failed called to sysconf"
#define GT_ERROR_SYS_STAT "Could not stat file '%s'"
#define GT_ERROR_SYS_OPEN "Could not open file '%s'"

#define GT_ERROR_STRDUP "strdup error"
#define GT_ERROR_PRINT_FORMAT "Incorrect print format. Expected format character"

/*
 * General purpose checkers
 */
#define GT_INVALID_CASE() gt_fatal_error(SELECTION_NOT_VALID)
#define GT_NOT_IMPLEMENTED() gt_fatal_error(NOT_IMPLEMENTED)
//#ifdef GT_ECLIPSE
  /* Eclipse debugging definitions */
//  #define GT_NULL_CHECK(object) if ((object)==NULL) { printf("%d\n",(int)(*(int*)object)); }
//  #define GT_ZERO_CHECK(object) if ((object)==0) { printf("%d\n",(*(int*)0)); }
//#else
  #define GT_NULL_CHECK(object) gt_fatal_check((object)==NULL,NULL_HANDLER,((char*)GT_QUOTE(object)))
  #define GT_ZERO_CHECK(object) gt_fatal_check((object)==0,NOT_ZERO,((char*)GT_QUOTE(object)))
//#endif


#endif /* GT_ERROR_H_ */
