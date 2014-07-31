/*
 * PROJECT: GEMMapper
 * FILE: buffered_file_printer.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef BUFFERED_FILE_PRINTER_H_
#define BUFFERED_FILE_PRINTER_H_

#include "essentials.h"

typedef struct {
  FILE* output_file;
  uint8_t* buffer;
  uint64_t cursor;
} buffered_file_printer_t;

/*
 * Checkers
 */
#define BUFFERED_FILE_PRINTER_CHECK(buffered_file_printer) GEM_CHECK_NULL(buffered_file_printer)

/*
 * Setup
 */
GEM_INLINE buffered_file_printer_t* buffered_file_printer_new(FILE* const output_file);
GEM_INLINE void buffered_file_printer_close(buffered_file_printer_t* const buffered_file_printer);

/*
 * Buffered File Printers
 */
GEM_INLINE void bfp_putc(buffered_file_printer_t* const buffered_file_printer,const char c);

#endif /* BUFFERED_FILE_PRINTER_H_ */
