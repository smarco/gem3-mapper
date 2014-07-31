/*
 * PROJECT: GEMMapper
 * FILE: buffered_file_printer.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */


#include "buffered_file_printer.h"

#define BUFFERED_FILE_PRINTER_BUFFER_SIZE BUFFER_SIZE_16M

/*
 * Setup
 */
GEM_INLINE buffered_file_printer_t* buffered_file_printer_new(FILE* const output_file) {
  buffered_file_printer_t* const buffered_file_printer = mm_alloc(buffered_file_printer_t);
  buffered_file_printer->output_file = output_file;
  buffered_file_printer->buffer = mm_malloc(BUFFERED_FILE_PRINTER_BUFFER_SIZE);
  buffered_file_printer->cursor = 0;
  return buffered_file_printer;
}
GEM_INLINE void buffered_file_printer_close(buffered_file_printer_t* const buffered_file_printer) {
  if (buffered_file_printer->cursor) {
    fprintf(buffered_file_printer->output_file,"%.*s",(int)buffered_file_printer->cursor,buffered_file_printer->buffer);
  }
  mm_free(buffered_file_printer->buffer);
  mm_free(buffered_file_printer);
}
/*
 * Buffered File Printers
 */
GEM_INLINE void bfp_putc(buffered_file_printer_t* const buffered_file_printer,const char c) {
  buffered_file_printer->buffer[buffered_file_printer->cursor++] = c;
  if (buffered_file_printer->cursor == BUFFERED_FILE_PRINTER_BUFFER_SIZE) {
    buffered_file_printer->cursor = 0;
    fprintf(buffered_file_printer->output_file,"%.*s",(int)BUFFERED_FILE_PRINTER_BUFFER_SIZE,buffered_file_printer->buffer);
  }
}
