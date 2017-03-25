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
 * DESCRIPTION: Error handling module
 */

#ifndef ERRORS_H_
#define ERRORS_H_

#include "system/commons.h"
#include "system/report.h"

// Internally to GEM error codes are returned as error_code
typedef int32_t error_code_t;

// Codes gem_status
#define GEM_STATUS_OK 1
#define GEM_STATUS_FAIL -1

/*
 * StackTrace Printer
 */
void gem_print_stack_trace(void);

/*
 * Error Signal Handler
 */
void gem_handle_error_signals(void);

/*
 * Print ErrNo
 */
void gem_perror(void);

/*
 * Handlers
 */
// Error Code
#define gem_error(gem_error_name,args...) \
  gem_report_error_begin_block(GEM_LABEL_ERROR,gem_error_name,##args) \
  gem_report_end_block(0,0,1)
#define gem_error_msg(gem_error_msg,args...) \
  gem_report_error_msg_begin_block(GEM_LABEL_ERROR,gem_error_msg,##args) \
  gem_report_end_block(0,0,1)
#define gem_fatal_error(gem_error_name,args...) \
  gem_report_error_begin_block(GEM_LABEL_FATAL_ERROR,gem_error_name,##args) \
  gem_report_end_block(1,1,1)
#define gem_fatal_error_msg(gem_error_msg,args...) \
  gem_report_error_msg_begin_block(GEM_LABEL_FATAL_ERROR,gem_error_msg,##args) \
  gem_report_end_block(1,1,1)
// Conditional Error Code
#define gem_cond_error(condition,gem_error_name,args...) \
  do { if (__builtin_expect((condition),0)) { gem_error(gem_error_name,##args); } } while (0)
#define gem_cond_fatal_error(condition,gem_error_name,args...) \
  do { if (__builtin_expect((condition),0)) { gem_fatal_error(gem_error_name,##args); } } while (0)
// Conditional Error Message
#define gem_cond_error_msg(condition,error_msg,args...) \
  do { if (__builtin_expect((condition),0)) { gem_error_msg(error_msg,##args); } } while (0)
#define gem_cond_fatal_error_msg(condition,error_msg,args...) \
  do { if (__builtin_expect((condition),0)) { gem_fatal_error_msg(error_msg,##args); } } while (0)
// Warning
#define gem_warn(gem_warning_name,args...) \
  gem_report_error_begin_block(GEM_LABEL_WARNING,gem_warning_name,##args) \
  gem_report_end_block(0,0,0)
#define gem_warn_msg(gem_warn_msg,args...) \
  gem_report_error_msg_begin_block(GEM_LABEL_WARNING,gem_warn_msg,##args) \
  gem_report_end_block(0,0,0)
// Conditional Warning
#define gem_cond_warn(condition,warning_name,args...) \
  do { if (__builtin_expect((condition),0)) { gem_warn(warning_name,##args); } } while (0)
#define gem_cond_warn_msg(condition,warn_msg,args...) \
  do { if (__builtin_expect((condition),0)) { gem_warn_msg(warn_msg,##args); } } while (0)
// Debug
#ifdef GEM_DEBUG
  #define gem_debug_msg(gem_debug_msg,args...) \
    gem_report_error_msg_begin_block(GEM_LABEL_DEBUG,gem_debug_msg,##args) \
    gem_report_end_block(1,0,0)
  #define gem_debug_block()
  #define gem_cond_debug_block(condition) if (condition)
#else
  #define gem_debug_msg(gem_debug_msg,args...)
  #define gem_debug_block() if (0)
  #define gem_cond_debug_block(condition) if (0)
#endif
// Conditional Debug
#define gem_cond_debug_msg(condition,debug_msg,args...) \
  do { if (__builtin_expect((condition),0)) { gem_debug_msg(debug_msg,##args); } } while (0)

/*
 * Checkers (Robustness check)
 */
#ifdef GEM_DEBUG
  #define gem_fatal_check(condition,gem_error_name,args...) gem_cond_fatal_error(condition,gem_error_name,##args)
  #define gem_fatal_check_msg(condition,error_msg,args...) gem_cond_fatal_error_msg(condition,error_msg,##args)
#else
  #define gem_fatal_check(condition,gem_error_name,args...)
  #define gem_fatal_check_msg(condition,error_msg,args...)
#endif

/*
 * ERROR CODES/MSG
 *   #define GEM_ERROR_<CODE> "<MSG>"
 */
// General checkers
#define GEM_CHECK_NULL(object) gem_fatal_check((object)==NULL,NULL_HANDLER_INFO,((char*)QUOTE(object)))
#define GEM_CHECK_ZERO(object) gem_fatal_check((object)==0,NOT_ZERO,((char*)QUOTE(object)))
#define GEM_CHECK_NOT_NEGATIVE(object) gem_fatal_check((object)<0,NOT_NEGATIVE,((char*)QUOTE(object)))
#define GEM_CHECK_POSITIVE(object) gem_fatal_check((object)<=0,POSITIVE,((char*)QUOTE(object)))
#define GEM_INVALID_CASE() gem_fatal_error(SELECTION_NOT_VALID)
#define GEM_NOT_SUPPORTED() gem_fatal_error(NOT_SUPPORTED)
#define GEM_NOT_IMPLEMENTED() gem_fatal_error(NOT_IMPLEMENTED)
#define GEM_CUDA_NOT_SUPPORTED() gem_fatal_error(CUDA_NOT_SUPPORTED)
#define GEM_UNREACHABLE_CODE() gem_fatal_error(UNREACHABLE_CODE)
#define GEM_INTERNAL_CHECK(cond,msg) gem_fatal_check(!(cond),INTERNAL_CHECK,msg)
// General errors
#define GEM_ERROR_NULL_HANDLER "Null handler or fields not properly allocated"
#define GEM_ERROR_NULL_HANDLER_INFO "Null handler %s"
#define GEM_ERROR_NOT_ZERO "Invalid Value. Variable %s must be non-zero"
#define GEM_ERROR_NOT_NEGATIVE "Invalid Value. Variable %s must be non-negative (>=0)"
#define GEM_ERROR_POSITIVE "Invalid Value. Variable %s must be strictly positive (>0)"
#define GEM_ERROR_INVALID_VALUE "Invalid value: '%s'"
#define GEM_ERROR_SELECTION_NOT_VALID "Library error. Selection not valid"
#define GEM_ERROR_NOT_IMPLEMENTED "Function/Feature not implemented yet (Sorry for the inconvenience)"
#define GEM_ERROR_NOT_SUPPORTED "Function/Feature not supported"
#define GEM_ERROR_CUDA_NOT_SUPPORTED "No CUDA support detected"
#define GEM_ERROR_UNREACHABLE_CODE "Unreachable code, please report"
#define GEM_ERROR_POSITION_OUT_OF_RANGE "Requested position (%"PRIu64") out of range [%"PRIu64",%"PRId64"]"
#define GEM_ERROR_ALG_INCONSISNTENCY "Algorithmic inconsistency, please report (Sorry for the inconvenience)"
#define GEM_ERROR_INTERNAL_CHECK "Internal check failed, please report '%s'"
// String
#define GEM_ERROR_STRDUP "strdup error"
#define GEM_ERROR_PRINT_FORMAT "Incorrect print format. Expected format character"
// System errors
#define GEM_ERROR_SYS_ERROR "System error signal raised"
#define GEM_ERROR_SYS_MMAP "Mmap call error"
#define GEM_ERROR_SYS_MMAP_FILE "Could not map file '%s' to memory"
#define GEM_ERROR_SYS_UNMAP "Could not unmap memory"
#define GEM_ERROR_SYS_THREAD_CREATE "Could not create thread"
#define GEM_ERROR_SYS_THREAD_JOIN "Could not join thread"
#define GEM_ERROR_SYS_MUTEX "Mutex call error"
#define GEM_ERROR_SYS_MUTEX_INIT "Mutex initialization error"
#define GEM_ERROR_SYS_MUTEX_DESTROY "Mutex destroy call error"
#define GEM_ERROR_SYS_COND_VAR "Conditional variable call error"
#define GEM_ERROR_SYS_COND_VAR_INIT "Conditional variable initialization error"
#define GEM_ERROR_SYS_COND_VAR_DESTROY "Conditional variable destroy call error"
#define GEM_ERROR_SYS_MKSTEMP "Could not create a temporal file (mkstemp:'%s')"
#define GEM_ERROR_SYS_HANDLE_TMP "Failed to handle temporal file"
#define GEM_ERROR_SYS_SYSCONF "Failed called to sysconf"
#define GEM_ERROR_SYS_STAT "Could not stat file '%s'"
#define GEM_ERROR_SYS_OPEN "Could not open file '%s'"

#endif /* ERRORS_H_ */
