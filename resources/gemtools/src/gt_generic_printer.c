/*
 * PROJECT: GEM-Tools library
 * FILE: gt_generic_printer.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_generic_printer.h"

#define GT_GEN_PRINTER_INITIAL_STRING_SIZE GT_BUFFER_SIZE_1K

/*
 * Generic printer
 */
GT_INLINE void gt_generic_new_file_printer(gt_generic_printer* const generic_printer,FILE* const file) {
  GT_NULL_CHECK(generic_printer);
  GT_NULL_CHECK(file);
  generic_printer->printer_type = GT_FILE_PRINTER;
  generic_printer->file = file;
}
GT_INLINE void gt_generic_new_string_printer(gt_generic_printer* const generic_printer,gt_string* const string) {
  GT_NULL_CHECK(generic_printer);
  GT_STRING_CHECK(string);
  generic_printer->printer_type = GT_STRING_PRINTER;
  generic_printer->string = string;
}
GT_INLINE void gt_generic_new_buffer_printer(gt_generic_printer* const generic_printer,gt_output_buffer* const output_buffer) {
  GT_NULL_CHECK(generic_printer);
  GT_OUTPUT_BUFFER_CHECK(output_buffer);
  generic_printer->printer_type = GT_BUFFER_PRINTER;
  generic_printer->output_buffer = output_buffer;
}
GT_INLINE void gt_generic_new_output_file_printer(gt_generic_printer* const generic_printer,gt_output_file* const output_file) {
  GT_NULL_CHECK(generic_printer);
  GT_OUTPUT_FILE_CHECK(output_file);
  generic_printer->printer_type = GT_OUTPUT_FILE_PRINTER;
  generic_printer->output_file = output_file;
}
GT_INLINE void gt_generic_new_buffered_output_file_printer(gt_generic_printer* const generic_printer,gt_buffered_output_file* const buffered_output_file) {
  GT_NULL_CHECK(generic_printer);
  GT_BUFFERED_OUTPUT_FILE_CHECK(buffered_output_file);
  generic_printer->printer_type = GT_BOF_PRINTER;
  generic_printer->buffered_output_file = buffered_output_file;
}

GT_INLINE gt_status gt_vgprintf(gt_generic_printer* const generic_printer,const char *template,va_list v_args) {
  GT_GENERIC_PRINTER_CHECK(generic_printer);
  GT_NULL_CHECK(template);
  gt_status chars_printed = 0;
  switch (generic_printer->printer_type) {
    case GT_FILE_PRINTER:
      GT_NULL_CHECK(generic_printer->file);
      gt_cond_fatal_error( (chars_printed=
          vfprintf(generic_printer->file,template,v_args))<0,FPRINTF);
      break;
    case GT_STRING_PRINTER:
      GT_STRING_CHECK(generic_printer->string);
      gt_cond_fatal_error( (chars_printed=
          gt_vsprintf_append(generic_printer->string,template,v_args))<0,SPRINTF);
      break;
    case GT_BUFFER_PRINTER:
      GT_OUTPUT_BUFFER_CHECK(generic_printer->output_buffer);
      gt_cond_fatal_error( (chars_printed=
          gt_vbprintf(generic_printer->output_buffer,template,v_args))<0,BPRINTF);
      break;
    case GT_OUTPUT_FILE_PRINTER:
      GT_OUTPUT_FILE_CHECK(generic_printer->output_file);
      gt_cond_fatal_error( (chars_printed=
          gt_vofprintf(generic_printer->output_file,template,v_args))<0,OFPRINTF);
      break;
    case GT_BOF_PRINTER:
      GT_BUFFERED_OUTPUT_FILE_CHECK(generic_printer->buffered_output_file);
      gt_cond_fatal_error( (chars_printed=
          gt_vbofprintf(generic_printer->buffered_output_file,template,v_args))<0,BOFPRINTF);
      break;
    default:
      GT_INVALID_CASE();
      break;
  }
  return chars_printed;
}
GT_INLINE gt_status gt_gprintf(gt_generic_printer* const generic_printer,const char *template,...) {
  GT_GENERIC_PRINTER_CHECK(generic_printer);
  GT_NULL_CHECK(template);
  va_list v_args;
  va_start(v_args,template);
  const gt_status chars_printed = gt_vgprintf(generic_printer,template,v_args);
  va_end(v_args);
  return chars_printed;
}

