/*
 * PROJECT: GEMMapper
 * FILE: fm.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: File management module.
 *   Provides functionality for an easy file management providing stream-based access for reading & writing
 */

#ifndef FILE_MANAGEMENT_H_
#define FILE_MANAGEMENT_H_

#include "commons.h"
#include "mm.h"

/*
 * I/O Constants/Values
 */
extern const int fm_open_flags[3];
extern const mode_t fm_open_mode[3];

typedef enum { FM_READ=0, FM_WRITE=1, FM_READ_WRITE=2 } fm_mode;
typedef enum { FM_STREAM, FM_REGULAR_FILE, FM_GZIPPED_FILE, FM_BZIPPED_FILE } fm_type;
typedef struct _fm_t fm_t;

/*
 * Checkers
 */
#define FM_CHECK(fm) GEM_CHECK_NULL(fm)

/*
 * Setup
 */
GEM_INLINE fm_t* fm_open_file(char* const file_name,const fm_mode mode);
GEM_INLINE fm_t* fm_open_FILE(FILE* const stream,const fm_mode mode);
GEM_INLINE fm_t* fm_open_gzFILE(FILE* const stream,const fm_mode mode);
GEM_INLINE fm_t* fm_open_bzFILE(FILE* const stream,const fm_mode mode);
GEM_INLINE fm_t* fm_open_temp_file();
GEM_INLINE void fm_close(fm_t* const file_manager);

/*
 * Accessors
 */
GEM_INLINE uint64_t fm_get_current_position(fm_t* const file_manager);
GEM_INLINE bool fm_eof(fm_t* const file_manager);
GEM_INLINE char* fm_get_file_name(fm_t* const file_manager);
GEM_INLINE uint64_t fm_get_file_size(fm_t* const file_manager);

/*
 * Seek
 */
GEM_INLINE void fm_seek(fm_t* const file_manager,const uint64_t position);

GEM_INLINE void fm_skip_forward(fm_t* const file_manager,const uint64_t num_bytes);
GEM_INLINE void fm_skip_uint64(fm_t* const file_manager);
GEM_INLINE void fm_skip_uint32(fm_t* const file_manager);
GEM_INLINE void fm_skip_uint16(fm_t* const file_manager);
GEM_INLINE void fm_skip_uint8(fm_t* const file_manager);
GEM_INLINE void fm_skip_align(fm_t* const file_manager,const uint64_t num_bytes);
GEM_INLINE void fm_skip_align_16(fm_t* const file_manager);
GEM_INLINE void fm_skip_align_32(fm_t* const file_manager);
GEM_INLINE void fm_skip_align_64(fm_t* const file_manager);
GEM_INLINE void fm_skip_align_128(fm_t* const file_manager);
GEM_INLINE void fm_skip_align_512(fm_t* const file_manager);
GEM_INLINE void fm_skip_align_1024(fm_t* const file_manager);
GEM_INLINE void fm_skip_align_4KB(fm_t* const file_manager);
GEM_INLINE void fm_skip_align_mempage(fm_t* const file_manager);

/*
 * Read
 */
#define fm_read(file_manager,var) fm_copy_mem(file_manager,&var,sizeof(var))
GEM_INLINE uint64_t fm_read_uint64(fm_t* const file_manager);
GEM_INLINE uint32_t fm_read_uint32(fm_t* const file_manager);
GEM_INLINE uint16_t fm_read_uint16(fm_t* const file_manager);
GEM_INLINE uint8_t fm_read_uint8(fm_t* const file_manager);
GEM_INLINE uint64_t fm_read_mem(fm_t* const file_manager,void* const dst,const uint64_t num_bytes);
GEM_INLINE uint64_t fm_read_mem_parallel(
    fm_t* const file_manager,void* const dst,const uint64_t num_bytes,const uint64_t num_threads);

GEM_INLINE mm_t* fm_load_mem(fm_t* const file_manager,const uint64_t num_bytes);

/*
 * Write
 */
#define fm_write(file_manager,var) fm_write_mem(file_manager,&var,sizeof(var))
GEM_INLINE void fm_write_uint64(fm_t* const file_manager,const uint64_t data);
GEM_INLINE void fm_write_uint32(fm_t* const file_manager,const uint32_t data);
GEM_INLINE void fm_write_uint16(fm_t* const file_manager,const uint16_t data);
GEM_INLINE void fm_write_uint8(fm_t* const file_manager,const uint8_t data);
GEM_INLINE void fm_write_mem(fm_t* const file_manager,const void* const src,const uint64_t num_bytes);

/*
 * Bulk Read of file
 */
GEM_INLINE void fm_bulk_read_fd(const int fd,void* const dst,const uint64_t size);
GEM_INLINE void fm_bulk_read_file(char* const file_name,void* const dst,const uint64_t offset,const uint64_t size);
GEM_INLINE void fm_bulk_read_file_parallel(
    char* const file_name,void* const dst,const uint64_t offset,const uint64_t size,const uint64_t num_threads);

/*
 * FileManager Wrappers
 */
GEM_INLINE void gem_stat(char* const file_name,struct stat *stat_info);
GEM_INLINE int gem_open_fd(char* const file_name,const int flags,const mode_t mode);
GEM_INLINE FILE* gem_open_FILE(char* const file_name,const char* opentype);
GEM_INLINE void gem_unlink(char* const file_name);

/*
 * Utils
 */
GEM_INLINE bool gem_access(char* const path,const fm_mode mode);

/*
 * FileManager Printers
 */
GEM_INLINE int vfmprintf(fm_t* const file_manager,const char *template,va_list v_args);
GEM_INLINE int fmprintf(fm_t* const file_manager,const char *template,...);

/*
 * Error Messages
 */
#define GEM_ERROR_FM_NO_ZLIB_SUPPORT "ZLIB functionality disable for this compilation"
#define GEM_ERROR_FM_NO_BZLIB_SUPPORT "BZLIB functionality disable for this compilation"

#define GEM_ERROR_FM_OPEN "Could not open file '%s'"
#define GEM_ERROR_FM_FDOPEN "Could not open file '%s'"
#define GEM_ERROR_FM_GZOPEN "Could not open GZip file"
#define GEM_ERROR_FM_BZOPEN "Could not open BZip file"
#define GEM_ERROR_FM_STAT "Could not stat file '%s'"
#define GEM_ERROR_FM_CLOSE "Could not close file '%s'"
#define GEM_ERROR_FM_GZCLOSE "Could not close GZip file '%s'"
#define GEM_ERROR_FM_BZCLOSE "Could not close BZip file '%s'"
#define GEM_ERROR_FM_FILENO "Invalid stream, fileno() call failed"
#define GEM_ERROR_FM_SEEK "Could not seek file '%s' to %"PRIu64
#define GEM_ERROR_FM_UNLINK "Could not unlink file '%s'"

#define GEM_ERROR_FM_WRITE "Could not write to file '%s'"
#define GEM_ERROR_FM_GZWRITE "Could not write to GZip file '%s'"
#define GEM_ERROR_FM_BZWRITE "Could not write to BZip file '%s'"

#define GEM_ERROR_FM_READ "Could not read %"PRIu64" bytes from file '%s'"
#define GEM_ERROR_FM_READ_ZERO "Zero bytes read from file '%s'"
#define GEM_ERROR_FM_FDREAD "Could not read %"PRIu64" bytes from fd(%d)"
#define GEM_ERROR_FM_GZREAD "Could not read %"PRIu64" bytes from GZip file '%s'"
#define GEM_ERROR_FM_BZREAD "Could not read %"PRIu64" bytes from BZip file '%s'"

#define GEM_ERROR_FM_INVALID_MODE_WRITE "Invalid file mode. File '%s' cannot be written"
#define GEM_ERROR_FM_INVALID_MODE_READ "Invalid file mode. File '%s' cannot be read"
#define GEM_ERROR_FM_NOT_SEEKABLE "Cannot seek file '%s' "

#define GEM_ERROR_FM_LOAD "Could not load memory chunk (size=%lu,read=%lu)"
#define GEM_ERROR_FM_DUP  "Could not duplicate file '%s' (not regular FILE or stream)"

#endif /* FILE_MANAGEMENT_H_ */

