/*
 * PROJECT: GEMMapper
 * FILE: errors.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Error handling module
 */

#ifndef ERRORS_H_
#define ERRORS_H_

#include "commons.h"
#include "report.h"

// Internally to GEM error codes are returned as error_code
typedef int32_t error_code_t;

// Codes gem_status
#define GEM_STATUS_OK 1
#define GEM_STATUS_FAIL -1

/*
 * StackTrace Printer
 */
void gem_print_stack_trace();

/*
 * Error Signal Handler
 */
void gem_handle_error_signals();

/*
 * Print ErrNo
 */
void gem_perror();

/*
 * GEM-Error handlers
 */
#define gem_fatal_error(gem_error_name,args...) \
  gem_report_error_begin_block(GEM_LABEL_FATAL_ERROR,gem_error_name,##args) \
  gem_report_end_block(1,1,1,1)
#define gem_fatal_error__perror(gem_error_name,args...) \
  gem_report_error_begin_block(GEM_LABEL_FATAL_ERROR,gem_error_name,##args) \
  gem_perror(); \
  gem_report_end_block(1,1,1,1)
#define gem_error(gem_error_name,args...) \
  gem_report_error_begin_block(GEM_LABEL_ERROR,gem_error_name,##args) \
  gem_report_end_block(1,0,0,0)
/*
 * ErrorMsg handlers
 */
#define gem_fatal_error_msg(gem_error_msg,args...) \
  gem_report_begin_block(GEM_LABEL_FATAL_ERROR,error,gem_error_msg,##args) \
  gem_report_end_block(1,1,1,1)
#define gem_error_msg(gem_error_msg,args...) \
  gem_report_begin_block(GEM_LABEL_ERROR,error,gem_error_msg,##args) \
  gem_report_end_block(1,0,0,0)
/*
 * GEM-Warning handler
 */
#define gem_warn(gem_warning_name,args...) \
  gem_report_error_begin_block(GEM_LABEL_WARNING,gem_warning_name,##args) \
  gem_report_end_block(0,0,0,0)
/*
 * GEM-Exception handlers (conditional error handlers)
 */
#define gem_cond_error(condition,gem_error_name,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gem_error(gem_error_name,##args); \
    } \
  } while (0)
#define gem_cond_fatal_error(condition,gem_error_name,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gem_fatal_error(gem_error_name,##args); \
    } \
  } while (0)
#define gem_cond_fatal_error__perror(condition,gem_error_name,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gem_fatal_error__perror(gem_error_name,##args); \
    } \
  } while (0)
/*
 * Exception handlers (conditional error handlers)
 */
#define gem_cond_fatal_error_msg(condition,error_msg,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gem_fatal_error_msg(error_msg,##args); \
    } \
  } while (0)
#define gem_cond_error_msg(condition,error_msg,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gem_error_msg(error_msg,##args); \
    } \
  } while (0)
/*
 * Robust checkers
 */
#ifdef GEM_DEBUG
  #define gem_fatal_check(condition,gem_error_name,args...) gem_cond_fatal_error(condition,gem_error_name,##args)
  #define gem_check(condition,gem_error_name,args...) gem_cond_error(condition,gem_error_name,##args)
  #define gem_check_block(condition) if (condition)
#else
  #define gem_fatal_check(condition,gem_error_name,args...)
  #define gem_check(condition,gem_error_name,args...)
  #define gem_check_block(condition) if (0)
#endif
/*
 * Debug Utilities
 */
#define gem_cond_debug_msg(condition,gem_debug_msg,args...) \
  do { \
    if (__builtin_expect((condition),0)){ \
      gem_debug_msg(gem_debug_msg,##args); \
    } \
  } while (0)
#define gem_cond_debug_block(condition) if (condition)
#ifdef GEM_DEBUG
  #define gem_debug(gem_debug_msg,args...) \
    gem_report_raw_begin_block(GEM_LABEL_DEBUG,debug,gem_debug_msg,##args) \
    gem_report_end_block(1,0,0,0)
  #define gem_debug_msg(gem_debug_msg,args...) \
    gem_report_begin_block(GEM_LABEL_DEBUG,debug,gem_debug_msg,##args) \
    gem_report_end_block(1,0,0,0)
  #define gem_debug_msg__stack(gem_debug_msg,args...) \
    gem_debug_msg(gem_debug_msg,args...); \
    gem_print_stack_trace();
  #define gem_debug_block()
#else
  #define gem_debug(gem_debug_msg,args...)
  #define gem_debug_msg(gem_debug_msg,args...)
  #define gem_debug_msg__stack(gem_debug_msg,args...)
  #define gem_debug_block() if (0)
#endif

/*
 * ERROR CODES/MSG
 *   #define GEM_ERROR_<CODE> "<MSG>"
 */
/*
 * General purpose checkers/errors
 */
#define GEM_CHECK_NULL(object) gem_fatal_check((object)==NULL,NULL_HANDLER_INFO,((char*)QUOTE(object)))
#define GEM_CHECK_ZERO(object) gem_fatal_check((object)==0,NOT_ZERO,((char*)QUOTE(object)))
#define GEM_CHECK_NOT_NEGATIVE(object) gem_fatal_check((object)<0,NOT_NEGATIVE,((char*)QUOTE(object)))
#define GEM_CHECK_POSITIVE(object) gem_fatal_check((object)<=0,POSITIVE,((char*)QUOTE(object)))
#define GEM_INVALID_CASE() gem_fatal_error(SELECTION_NOT_VALID)
#define GEM_NOT_SUPPORTED() gem_fatal_error(NOT_SUPPORTED)
#define GEM_NOT_IMPLEMENTED() gem_fatal_error(NOT_IMPLEMENTED)
#define GEM_INTERNAL_CHECK(cond,msg) gem_fatal_check(!(cond),INTERNAL_CHECK,msg)
/* Eclipse debugging definitions */
  // #define GEM_CHECK_NULL(object) if ((object)==NULL) { printf("%d\n",(*(int*)object)); }
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

//
//// String errors
//#define GEM_ERROR_STRING_STATIC "Could not perform operation on static string"
//
//// File errors
//#define GEM_ERROR_FILE_FORMAT "Could not determine file format"
//
//// Output errors
//#define GEM_ERROR_FPRINTF "Printing output. 'fprintf' call failed"
//#define GEM_ERROR_SPRINTF "Printing output. 'sprintf' call failed"
//#define GEM_ERROR_TPRINTF "Printing output. 'tprintf' call failed"
//#define GEM_ERROR_BPRINTF "Printing output. Buffer print formated 'gem_bprintf' call failed"
//#define GEM_ERROR_OFPRINTF "Printing output. Output File print formated 'gem_ofprintf' call failed"
//#define GEM_ERROR_BOFPRINTF "Printing output. Buffered Output file print formated 'gem_bofprintf' call failed"
//
//
//
//// Template/Alignment/Map/Misms errors
//#define GEM_ERROR_MISMS_TYPE "Misms incorrect type"
//#define GEM_ERROR_MISMS_TRANSITION "Incorrect mismatch transition (Same base at both sides)"
//#define GEM_ERROR_MISMS_SPLICE_POS "Splicing distance must be positive (non-zero)"
//#define GEM_ERROR_COUNTERS_POS_STRATUM "Stratum must be strictly positive (stratum>0)"
//#define GEM_ERROR_MAP_MISMS_NOT_PARSED "Map's mismatches not parsed yet"
//#define GEM_ERROR_MAP_NEG_LENGTH "Negative Map total length"
//#define GEM_ERROR_MAP_NEG_MAPPED_BASES "Negative number of bases mapped"
//#define GEM_ERROR_ALIGNMENT_READ_QUAL_LENGTH "Read and quality length differs"
//#define GEM_ERROR_ALIGNMENT_MAPS_NOT_PARSED "Alignment's maps not parsed yet"
//#define GEM_ERROR_ALIGNMENT_INCONSISTENT_COUNTERS "Alignment inconsistency. Maps inconsistent with counters values"
//#define GEM_ERROR_TEMPLATE_ZERO_BLOCKS "Zero alignment blocks (num_blocks_template==0)"
//#define GEM_ERROR_TEMPLATE_TOO_MANY_BLOCKS "Template contains already two ends"
//#define GEM_ERROR_TEMPLATE_BLOCKS_EXCESS "Template blocks exceeds two ends"
//#define GEM_ERROR_TEMPLATE_MMAP_NULL "Template mmap is null (all maps are null)"
//#define GEM_ERROR_TEMPLATE_INCONSISTENT_MMAPS_ALIGNMENT "Template inconsistency. Multimaps' members must be contained by single alignments"
//#define GEM_ERROR_TEMPLATE_INCONSISTENT_MMAPS_ATTRB_RELATION "Template inconsistency. Incorrect number of mmaps and mmaps' attributes"
//#define GEM_ERROR_TEMPLATE_INCONSISTENT_NUM_MAPS_RELATION "Template inconsistency. Incorrect number of matches' elements (check num_blocks_template)"
//#define GEM_ERROR_TEMPLATE_INCONSISTENT_NUM_BLOCKS "Template inconsistency. Number of blocks must be the same across templates"
//#define GEM_ERROR_TEMPLATE_INCONSISTENT_COUNTERS "Template inconsistency. MMaps inconsistent with counters values"
//#define GEM_ERROR_TEMPLATE_ADD_BAD_NUM_BLOCKS "Trying to add wrong number of blocks to the template"
//#define GEM_ERROR_PALIGN_BAD_NUM_BLOCKS "Invalid Paired-alignment. Wrong number of alignment blocks (%"PRIu64")"
//
//// Sequence Archive/Segmented Sequence errors
//#define GEM_ERROR_SEGMENTED_SEQ_IDX_OUT_OF_RANGE "Error accessing segmented sequence. Index %"PRIu64" out out range [0,%"PRIu64")"
//#define GEM_ERROR_CDNA_IT_OUT_OF_RANGE "Error seeking sequence. Index %"PRIu64" out out range [0,%"PRIu64")"
//#define GEM_ERROR_SEQ_ARCHIVE_WRONG_TYPE "Wrong sequence archive type"
//#define GEM_ERROR_SEQ_ARCHIVE_NOT_FOUND "Sequence '%s' not found in reference archive"
//#define GEM_ERROR_SEQ_ARCHIVE_POS_OUT_OF_RANGE "Requested position '%"PRIu64"' out of sequence boundaries"
//#define GEM_ERROR_SEQ_ARCHIVE_CHUNK_OUT_OF_RANGE "Requested sequence string [%"PRIu64",%"PRIu64") out of sequence '%s' boundaries"
//#define GEM_ERROR_GEMIDX_SEQ_ARCHIVE_NOT_FOUND "GEMIdx. Sequence '%s' not found in reference archive"
//#define GEM_ERROR_GEMIDX_INTERVAL_NOT_FOUND "GEMIdx. Interval relative to sequence '%s' not found in reference archive"
//
///*
// * Parsing FASTQ File format errors
// */
//// IFP (Input FASTA Parser). General
//#define GEM_ERROR_PARSE_FASTA "Parsing FASTA/FASTQ error(%s:%"PRIu64":%"PRIu64")"
//
///*
// * Parsing MAP File format errors
// */
//// IMP (Input MAP Parser). General
//#define GEM_ERROR_PARSE_MAP "Parsing MAP error(%s:%"PRIu64")"
//#define GEM_ERROR_PARSE_MAP_BAD_FILE_FORMAT "Parsing MAP error(%s:%"PRIu64"). Not a MAP file"
//#define GEM_ERROR_PARSE_MAP_BAD_NUMBER_FIELDS "Parsing MAP error(%s:%"PRIu64"). Wrong number of TAB separated fields (%"PRIu64")"
//#define GEM_ERROR_PARSE_MAP_BAD_READ_QUAL_LENGTH "Parsing MAP error(%s:%"PRIu64"). Mismatching Read length (%"PRIu64") and Quality length (%"PRIu64")"
//#define GEM_ERROR_PARSE_MAP_COUNTERS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Error parsing counters"
//#define GEM_ERROR_PARSE_MAP_BAD_TEMPLATE_SEP "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Read character '%c' not valid (%s)"
//#define GEM_ERROR_PARSE_MAP_DIFF_TEMPLATE_BLOCKS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Different number of template blocks {read(%"PRIu64"),qualities(%"PRIu64")}"
//#define GEM_ERROR_PARSE_MAP_NOT_AN_ALIGNMENT "Parsing MAP error(%s:%"PRIu64"). File doesn't contains simple alignments (use template)"
//#define GEM_ERROR_PARSE_MAP_MISMS_ALREADY_PARSED "Parsing MAP error(%s:%"PRIu64"). Mismatch string already parsed or null lazy-parsing handler"
//#define GEM_ERROR_PARSE_MAP_NOT_IMPLEMENTED "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Feature not implemented yet (sorry)"
//#define GEM_ERROR_PARSE_MAP_PREMATURE_EOL "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Premature End-of-line found"
//// IMP (Input MAP Parser). Parsing Read Errors
//#define GEM_ERROR_PARSE_MAP_READ_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing read, bad character found"
//// IMP (Input MAP Parser). Parsing Qualities Errors
//#define GEM_ERROR_PARSE_MAP_QUAL_BAD_SEPARATOR "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing quality string, bad block-separator found"
//#define GEM_ERROR_PARSE_MAP_QUAL_BAD_LENGTH "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing quality string, wrong length (w.r.t. read length)"
//#define GEM_ERROR_PARSE_MAP_QUAL_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing quality string, wrong quality value (bad character)"
//// IMP (Input MAP Parser). Parsing Counters Errors
//#define GEM_ERROR_PARSE_MAP_COUNTERS_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing counters, bad character found"
//// IMP (Input MAP Parser). Parsing Maps Errors
//#define GEM_ERROR_PARSE_MAP_BAD_NUMBER_OF_BLOCKS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing maps, wrong number of blocks"
//#define GEM_ERROR_PARSE_MAP_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing maps, bad character found"
//#define GEM_ERROR_PARSE_MAP_INCONSISTENT_BLOCKS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing maps, block(s) length doesn't match the read length"
//#define GEM_ERROR_PARSE_MAP_SPLIT_MAP_BAD_NUM_ACCEPTORS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing split-map, bad number of acceptors"
//#define GEM_ERROR_PARSE_MAP_SPLIT_MAP_BAD_NUM_DONORS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing split-map, bad number of donors"
//// IMP (Input MAP Parser). Parsing Mismatch String Errors
//#define GEM_ERROR_PARSE_MAP_MISMS_BAD_CHARACTER "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing mismatch string, bad character found"
//#define GEM_ERROR_PARSE_MAP_MISMS_BAD_MISMS_POS "Parsing MAP error(%s:%"PRIu64":%"PRIu64"). Parsing mismatch string, unsorted mismatches"
//
///*
// * Parsing SAM File format errors
// */
//// ISP (Input SAM Parser). General
//#define GEM_ERROR_PARSE_SAM "Parsing SAM error(%s:%"PRIu64":%"PRIu64")"
//#define GEM_ERROR_PARSE_SAM_BAD_FILE_FORMAT "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Not a SAM file"
//#define GEM_ERROR_PARSE_SAM_BAD_CHARACTER "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Bad character found"
//#define GEM_ERROR_PARSE_SAM_UNMAPPED_XA "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Unmapped read contains XA field (inconsistency)"
//#define GEM_ERROR_PARSE_SAM_PREMATURE_EOL "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Premature End-of-line found"
//#define GEM_ERROR_PARSE_SAM_EXPECTED_NUMBER "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Expected number"
//#define GEM_ERROR_PARSE_SAM_WRONG_READ_CONTENT "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Read in template doesn't match previously parse reads with same tag"
//#define GEM_ERROR_PARSE_SAM_CIGAR_PREMATURE_END "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Premature end of CIGAR string"
//#define GEM_ERROR_PARSE_SAM_WRONG_NUM_XA "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Wrong number of eXtra mAps (as to pair them)"
//#define GEM_ERROR_PARSE_SAM_UNSOLVED_PENDING_MAPS "Parsing SAM error(%s:%"PRIu64":%"PRIu64"). Failed to pair maps"
//
///*
// * Output File
// */
//#define GEM_ERROR_OUTPUT_FILE_FAIL_WRITE "Output file. Error writing to to file"
//#define GEM_ERROR_BUFFER_SAFETY_DUMP "Output buffer. Could not perform safety dump"
//
//#define GEM_ERROR_OUTPUT_SAM_NO_PRIMARY_ALG "Output SAM. No primary alignment specified"
//
///*
// * Map Alignment
// */
//#define GEM_ERROR_MAP_ALG_WRONG_ALG "(Re)Aligning Map. Wrong alignment"
//#define GEM_ERROR_MAP_RECOVER_MISMS_WRONG_BASE_ALG "Recovering mismatches from map. Wrong initial alignment"

#endif /* ERRORS_H_ */
