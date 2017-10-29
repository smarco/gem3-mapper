/*
 * PROJECT: GEM-Tools library
 * FILE: gt_generic_printer.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_GENERIC_PRINTER_H_
#define GT_GENERIC_PRINTER_H_

#include "gt_commons.h"
#include "gt_output_buffer.h"
#include "gt_output_file.h"
#include "gt_buffered_output_file.h"

typedef enum { GT_FILE_PRINTER, GT_STRING_PRINTER, GT_BUFFER_PRINTER, GT_OUTPUT_FILE_PRINTER, GT_BOF_PRINTER } gt_generic_printer_t;
typedef struct {
  gt_generic_printer_t printer_type;
  union {
    FILE* file;
    gt_string *string;
    gt_output_buffer* output_buffer;
    gt_output_file* output_file;
    gt_buffered_output_file* buffered_output_file;
  }; // Printer Object
} gt_generic_printer;

/*
 * Checkers
 */
#define GT_GENERIC_PRINTER_CHECK(generic_printer) \
  GT_NULL_CHECK(generic_printer); \
  switch ((generic_printer)->printer_type) { \
    case GT_FILE_PRINTER: GT_NULL_CHECK((generic_printer)->file); break; \
    case GT_STRING_PRINTER: GT_STRING_CHECK((generic_printer)->string); break; \
    case GT_BUFFER_PRINTER: GT_OUTPUT_BUFFER_CHECK((generic_printer)->output_buffer); break; \
    case GT_OUTPUT_FILE_PRINTER: GT_OUTPUT_FILE_CHECK(generic_printer->output_file); break;\
    case GT_BOF_PRINTER: GT_BUFFERED_OUTPUT_FILE_CHECK(generic_printer->buffered_output_file); break; \
    default: gt_fatal_check(true,SELECTION_NOT_VALID); break; \
  }

/*
 * Generic printer
 */
GT_INLINE void gt_generic_new_file_printer(gt_generic_printer* const generic_printer,FILE* const file);
GT_INLINE void gt_generic_new_string_printer(gt_generic_printer* const generic_printer,gt_string* const string);
GT_INLINE void gt_generic_new_buffer_printer(gt_generic_printer* const generic_printer,gt_output_buffer* const output_buffer);
GT_INLINE void gt_generic_new_output_file_printer(gt_generic_printer* const generic_printer,gt_output_file* const output_file);
GT_INLINE void gt_generic_new_buffered_output_file_printer(gt_generic_printer* const generic_printer,gt_buffered_output_file* const buffered_output_file);

GT_INLINE gt_status gt_vgprintf(gt_generic_printer* const generic_printer,const char *template,va_list v_args);
GT_INLINE gt_status gt_gprintf(gt_generic_printer* const generic_printer,const char *template,...);

/*
 * Automatic bindings generator
 */
#define GT_GENERIC_PRINTER_PROTOTYPE(MODULE_NAME,FUNCTION_NAME,SIGNATURE...) \
  GT_INLINE gt_status MODULE_NAME## _g   ##FUNCTION_NAME(gt_generic_printer* const gprinter,##SIGNATURE); \
  GT_INLINE gt_status MODULE_NAME## _f   ##FUNCTION_NAME(FILE* file,##SIGNATURE); \
  GT_INLINE gt_status MODULE_NAME## _s   ##FUNCTION_NAME(gt_string* const string,##SIGNATURE); \
  GT_INLINE gt_status MODULE_NAME## _b   ##FUNCTION_NAME(gt_output_buffer* const output_buffer,##SIGNATURE); \
  GT_INLINE gt_status MODULE_NAME## _of  ##FUNCTION_NAME(gt_output_file* const output_file,##SIGNATURE); \
  GT_INLINE gt_status MODULE_NAME## _bof ##FUNCTION_NAME(gt_buffered_output_file* const buffered_output_file,##SIGNATURE)
#define GT_GENERIC_PRINTER_IMPLEMENTATION(MODULE_NAME,FUNCTION_NAME,SIGNATURE...) \
  GT_INLINE gt_status MODULE_NAME##_f##FUNCTION_NAME(FILE* file,##SIGNATURE) { \
    GT_NULL_CHECK(file); \
    gt_generic_printer gprinter; \
    gt_generic_new_file_printer(&gprinter,file); \
    return MODULE_NAME##_g##FUNCTION_NAME(&gprinter,GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS); \
  } \
  GT_INLINE gt_status MODULE_NAME##_s##FUNCTION_NAME(gt_string* const string,##SIGNATURE) { \
    GT_STRING_CHECK(string); \
    gt_generic_printer gprinter; \
    gt_generic_new_string_printer(&gprinter,string); \
    return MODULE_NAME##_g##FUNCTION_NAME(&gprinter,GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS); \
  } \
  GT_INLINE gt_status MODULE_NAME##_b##FUNCTION_NAME(gt_output_buffer* const output_buffer,##SIGNATURE) { \
    GT_OUTPUT_BUFFER_CHECK(output_buffer); \
    gt_generic_printer gprinter; \
    gt_generic_new_buffer_printer(&gprinter,output_buffer); \
    return MODULE_NAME##_g##FUNCTION_NAME(&gprinter,GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS); \
  } \
  GT_INLINE gt_status MODULE_NAME##_of##FUNCTION_NAME(gt_output_file* const output_file,##SIGNATURE) { \
    GT_OUTPUT_FILE_CHECK(output_file); \
    gt_generic_printer gprinter; \
    gt_generic_new_output_file_printer(&gprinter,output_file); \
    return MODULE_NAME##_g##FUNCTION_NAME(&gprinter,GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS); \
  } \
  GT_INLINE gt_status MODULE_NAME##_bof##FUNCTION_NAME(gt_buffered_output_file* const buffered_output_file,##SIGNATURE) { \
    GT_BUFFERED_OUTPUT_FILE_CHECK(buffered_output_file); \
    gt_generic_printer gprinter; \
    gt_generic_new_buffered_output_file_printer(&gprinter,buffered_output_file); \
    return MODULE_NAME##_g##FUNCTION_NAME(&gprinter,GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS); \
  }

/*
 * Error Messages
 */
#define GT_ERROR_FPRINTF "Printing output. 'fprintf' call failed"
#define GT_ERROR_SPRINTF "Printing output. 'sprintf' call failed"
#define GT_ERROR_TPRINTF "Printing output. 'tprintf' call failed"
#define GT_ERROR_BPRINTF "Printing output. Buffer print formated 'gt_bprintf' call failed"
#define GT_ERROR_OFPRINTF "Printing output. Output File print formated 'gt_ofprintf' call failed"
#define GT_ERROR_BOFPRINTF "Printing output. Buffered Output file print formated 'gt_bofprintf' call failed"

#endif /* GT_OUTPUT_PRINTER_H_ */
