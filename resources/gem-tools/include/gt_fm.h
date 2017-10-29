/*
 * PROJECT: GEM-Tools library
 * FILE: gt_fm.h
 * DATE: 01/02/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   // TODO
 */

#ifndef GT_FILE_MANAGEMENT_H_
#define GT_FILE_MANAGEMENT_H_

#include "gt_commons.h"
#include "gt_error.h"
#include "gt_vector.h"
#include "gt_string.h"

#include <zlib.h>
#include <bzlib.h>

/*
 * I/O Constants/Values
 */
extern const int gt_fm_open_flags[3];
extern const mode_t gt_fm_open_mode[3];

typedef enum { GT_FM_READ=0, GT_FM_WRITE=1 } gt_fm_mode;
typedef enum { GT_FM_STREAM, GT_FM_REGULAR_FILE, GT_FM_GZIPPED_FILE, GT_FM_BZIPPED_FILE } gt_fm_type;
typedef struct {
  /* File */
  int fd;                 /* File descriptor */
  FILE* file;             /* FILE */
  gzFile gz_file;         /* GZip FILE */
  BZFILE* bz_file;        /* BZip FILE */
  /* Attributes */
  gt_fm_type file_type;   /* File type */
  gt_fm_mode mode;        /* File mode {R,W,R/W} */
  char *file_name;        /* File name */
  /* Locator */
  uint64_t byte_position; /* Current byte position */
  uint64_t file_size;     /* File size */
  bool eof;               /* End of file flag */
  /* Auxiliary Skip Buffers */
  uint8_t* skip_read_buffer;
  uint8_t* skip_write_buffer;
} gt_fm;

/*
 * Checkers
 */
#define GT_FM_CHECK(fm) GT_NULL_CHECK(fm)

/*
 * Setup
 */
GT_INLINE gt_fm* gt_fm_open_file(char* const file_name,const gt_fm_mode mode);
GT_INLINE gt_fm* gt_fm_open_FILE(FILE* const stream,const gt_fm_mode mode);
GT_INLINE gt_fm* gt_fm_open_gzFILE(FILE* const stream,const gt_fm_mode mode);
GT_INLINE gt_fm* gt_fm_open_bzFILE(FILE* const stream,const gt_fm_mode mode);
GT_INLINE void gt_fm_close(gt_fm* const fm);

/*
 * Accessors
 */
GT_INLINE uint64_t gt_fm_get_current_position(gt_fm* const fm);
GT_INLINE bool gt_fm_eof(gt_fm* const fm);
GT_INLINE char* gt_fm_get_file_name(gt_fm* const fm);

/*
 * Seek
 */
GT_INLINE void gt_fm_skip_forward(gt_fm* const fm,const uint64_t num_bytes);
GT_INLINE void gt_fm_skip_uint64(gt_fm* const fm);
GT_INLINE void gt_fm_skip_uint32(gt_fm* const fm);
GT_INLINE void gt_fm_skip_uint16(gt_fm* const fm);
GT_INLINE void gt_fm_skip_uint8(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align(gt_fm* const fm,const uint64_t num_bytes);
GT_INLINE void gt_fm_skip_align_16(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align_32(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align_64(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align_128(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align_512(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align_1024(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align_4KB(gt_fm* const fm);
GT_INLINE void gt_fm_skip_align_mempage(gt_fm* const fm);

/*
 * Read
 */
#define gt_fm_read(fm,var) gt_fm_copy_mem(fm,&var,sizeof(var))
GT_INLINE uint64_t gt_fm_read_uint64(gt_fm* const fm);
GT_INLINE uint32_t gt_fm_read_uint32(gt_fm* const fm);
GT_INLINE uint16_t gt_fm_read_uint16(gt_fm* const fm);
GT_INLINE uint8_t gt_fm_read_uint8(gt_fm* const fm);
GT_INLINE uint64_t gt_fm_read_mem(gt_fm* const fm,void* const dst,const uint64_t num_bytes);
GT_INLINE uint64_t gt_fm_read_mem_parallel(gt_fm* const fm,void* const dst,const uint64_t num_bytes,const uint64_t num_threads);

/*
 * Write
 */
#define gt_fm_write(fm,var) gt_fm_write_mem(fm,&var,sizeof(var))
GT_INLINE void gt_fm_write_uint64(gt_fm* const fm,const uint64_t data);
GT_INLINE void gt_fm_write_uint32(gt_fm* const fm,const uint32_t data);
GT_INLINE void gt_fm_write_uint16(gt_fm* const fm,const uint16_t data);
GT_INLINE void gt_fm_write_uint8(gt_fm* const fm,const uint8_t data);
GT_INLINE void gt_fm_write_mem(gt_fm* const fm,const void* const src,const uint64_t num_bytes);

/*
 * Bulk Read of file
 */
GT_INLINE void gt_fm_bulk_read_fd(const int fd,void* const dst,const uint64_t size);
GT_INLINE void gt_fm_bulk_read_file(char* const file_name,void* const dst,const uint64_t offset,const uint64_t size);
GT_INLINE void gt_fm_bulk_read_file_parallel(
    char* const file_name,void* const dst,const uint64_t offset,const uint64_t size,const uint64_t num_threads);

/*
 * FileManager Wrappers
 */
GT_INLINE void gt_stat(char* const file_name,struct stat *stat_info);
GT_INLINE int gt_open_fd(char* const file_name,const int flags,const mode_t mode);
GT_INLINE FILE* gt_open_FILE(char* const file_name,const char* opentype);

/*
 * FileManager Printers
 */
GT_INLINE gt_status gt_vfmprintf(gt_fm* const fm,const char *template,va_list v_args);
GT_INLINE gt_status gt_fmprintf(gt_fm* const fm,const char *template,...);

/*
 * Error Messages
 */
#define GT_ERROR_FM_NO_ZLIB_SUPPORT "ZLIB functionality disable for this compilation"
#define GT_ERROR_FM_NO_BZLIB_SUPPORT "BZLIB functionality disable for this compilation"

#define GT_ERROR_FM_OPEN "Could not open file '%s'"
#define GT_ERROR_FM_FDOPEN "Could not open file '%s'"
#define GT_ERROR_FM_GZOPEN "Could not open GZip file"
#define GT_ERROR_FM_BZOPEN "Could not open BZip file"
#define GT_ERROR_FM_STAT "Could not stat file '%s'"
#define GT_ERROR_FM_CLOSE "Could not close file '%s'"
#define GT_ERROR_FM_GZCLOSE "Could not close GZip file '%s'"
#define GT_ERROR_FM_BZCLOSE "Could not close BZip file '%s'"
#define GT_ERROR_FM_FILENO "Invalid stream, fileno() call failed"
#define GT_ERROR_FM_SEEK "Could not seek file '%s' to %"PRIu64

#define GT_ERROR_FM_WRITE "Could not write to file '%s'"
#define GT_ERROR_FM_GZWRITE "Could not write to GZip file '%s'"
#define GT_ERROR_FM_BZWRITE "Could not write to BZip file '%s'"

#define GT_ERROR_FM_READ "Could not read %"PRIu64" bytes from file '%s'"
#define GT_ERROR_FM_FDREAD "Could not read %"PRIu64" bytes from fd(%d)"
#define GT_ERROR_FM_GZREAD "Could not read %"PRIu64" bytes from GZip file '%s'"
#define GT_ERROR_FM_BZREAD "Could not read %"PRIu64" bytes from BZip file '%s'"

#define GT_ERROR_FM_INVALID_MODE_WRITE "Invalid file mode. File '%s' cannot be written"
#define GT_ERROR_FM_INVALID_MODE_READ "Invalid file mode. File '%s' cannot be read"
#define GT_ERROR_FM_NOT_SEEKABLE "Cannot seek file '%s' "

#endif /* GT_FILE_MANAGEMENT_H_ */

