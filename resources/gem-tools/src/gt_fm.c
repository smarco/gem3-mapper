/*
 * PROJECT: GEM-Tools library
 * FILE: gt_fm.c
 * DATE: 01/02/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   TODO
 *   1.- Incorporate the use of madvise() / readahead() / posix_fadvise()
 */

#include "gt_fm.h"
#include "gt_mm.h"
#include "gt_error.h"
#include "gt_vector.h"

#ifndef O_NOATIME
  #define O_NOATIME 0 /* FIXME: O_NOATIME is not allowed if only read rights are guaranteed */
#endif

/*
 * I/O Constants/Values
 */
// Buffer Size
#define GT_FM_BULK_COPY_BLOCK_SIZE GT_BUFFER_SIZE_2G
#define GT_FM_SKIP_BUFFER_SIZE     GT_BUFFER_SIZE_2M
// Open Flags
const int gt_fm_open_flags[3] = { O_RDONLY, O_WRONLY|O_CREAT|O_NOATIME, O_RDWR|O_CREAT|O_NOATIME };
const char* const gt_fm_file_open_flags[3] = { "rb", "wb", "wb+" };
const char* const gt_fm_gzfile_open_flags[2] = { "rb", "wb9"};

/*
 * Setup
 */
GT_INLINE void gt_fm_initialize(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  // Locator
  fm->eof = (fm->mode==GT_FM_READ) ? feof(fm->file) : false;
  fm->byte_position=0;
  // Init Special Files
  switch (fm->file_type) {
    case GT_FM_STREAM:
    case GT_FM_REGULAR_FILE:
      break;
#ifdef HAVE_ZLIB
    case GT_FM_GZIPPED_FILE: {
      fm->fd = fileno(fm->file);
      gt_cond_fatal_error__perror(fm->fd==-1,FM_FILENO);
      fm->gz_file = gzdopen(fm->fd,gt_fm_gzfile_open_flags[fm->mode]);
      gt_cond_fatal_error(!fm->gz_file,FM_GZOPEN);
      break;
    }
#endif
#ifdef HAVE_BZLIB
    case GT_FM_BZIPPED_FILE: {
      int bzerror;
      fm->bz_file = (fm->mode==GT_FM_READ) ?
          BZ2_bzReadOpen(&bzerror,fm->file,0,0,NULL,0) :
          BZ2_bzWriteOpen(&bzerror,fm->file,9,0,0);
      gt_cond_fatal_error(bzerror!=BZ_OK,FM_BZOPEN);
      break;
    }
#endif
    default:
      GT_INVALID_CASE();
      break;
  }
  // Auxiliary Skip Buffers
  fm->skip_read_buffer=NULL;
  fm->skip_write_buffer=NULL;
}
GT_INLINE void gt_fm_check_file_type(gt_fm* const fm,char* const file_name) {
  GT_NULL_CHECK(fm);
  // Check FILE stats
  struct stat stat_info;
  gt_cond_fatal_error__perror(stat(file_name,&stat_info)==-1,FM_STAT,file_name);
  // Open the FILE and check type
  FILE* file;
  gt_cond_fatal_error__perror(!(file=fopen(file_name,"r")),FM_OPEN,file_name);
  if (S_ISREG(stat_info.st_mode)) { // Regular file
    fm->file_type = GT_FM_REGULAR_FILE;
    fm->file_size = stat_info.st_size;
    // Check if GZIP/BZIP compressed
    unsigned char tbuf[4];
    if (fread(tbuf,1,4,file) == 4) {
      if(tbuf[0]==0x1f && tbuf[1]==0x8b && tbuf[2]==0x08) {
        fm->file_type = GT_FM_GZIPPED_FILE;
      } else if(tbuf[0]=='B' && tbuf[1]=='Z' && tbuf[2]=='h' && tbuf[3]>='0' && tbuf[3]<='9') {
        fm->file_type = GT_FM_BZIPPED_FILE;
      }
    }
  } else {
    fm->file_type = GT_FM_STREAM;
    fm->file_size = UINT64_MAX;
  }
  gt_cond_fatal_error__perror(fclose(file),FM_CLOSE,file_name);
}
GT_INLINE gt_fm* gt_fm_open_file(char* const file_name,const gt_fm_mode mode) {
  GT_NULL_CHECK(file_name);
  // Allocate handler
  gt_fm* fm = gt_alloc(gt_fm);
  // File
  fm->fd = open(file_name,gt_fm_open_flags[mode],S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
  gt_cond_fatal_error__perror(fm->fd==-1,FM_OPEN,file_name);
  fm->file = fdopen(fm->fd,gt_fm_file_open_flags[mode]);
  gt_cond_fatal_error__perror(fm->file==NULL,FM_FDOPEN,file_name);
  fm->gz_file = NULL;
  fm->bz_file = NULL;
  // Attributes
  fm->mode = mode;
  fm->file_name = strdup(file_name);
  gt_cond_fatal_error(fm->file_name==NULL,STRDUP);
  fm->file_size = UINT64_MAX;
  // Check BZIP/GZIP compress formats (low layer)
  if (mode==GT_FM_READ) {
    gt_fm_check_file_type(fm,file_name);
  } else {
    fm->file_type = GT_FM_REGULAR_FILE;
  }
  // Initialize file manager
  gt_fm_initialize(fm);
  // Return fm
  return fm;
}
GT_INLINE gt_fm* gt_fm_open_FILE(FILE* const stream,const gt_fm_mode mode) {
  GT_NULL_CHECK(stream);
  // Allocate handler
  gt_fm* fm = gt_alloc(gt_fm);
  // File
  fm->fd = 0;
  fm->file = stream;
  fm->gz_file = NULL;
  fm->bz_file = NULL;
  // Attributes
  fm->mode = mode;
  fm->file_name = GT_STREAM_FILE_NAME;
  fm->file_type = GT_FM_STREAM;
  fm->file_size = UINT64_MAX;
  // Initialize file manager
  gt_fm_initialize(fm);
  // Return fm
  return fm;
}
GT_INLINE gt_fm* gt_fm_open_gzFILE(FILE* const stream,const gt_fm_mode mode) {
#ifndef HAVE_BZLIB
  gt_fatal_error(FM_NO_ZLIB_SUPPORT);
  return NULL;
#else
  GT_NULL_CHECK(stream);
  // Allocate handler
  gt_fm* fm = gt_alloc(gt_fm);
  // File
  fm->fd = 0;
  fm->file = stream;
  fm->gz_file = NULL;
  fm->bz_file = NULL;
  // Attributes
  fm->mode = mode;
  fm->file_name = GT_STREAM_FILE_NAME;
  fm->file_type = GT_FM_GZIPPED_FILE;
  fm->file_size = UINT64_MAX;
  // Initialize file manager
  gt_fm_initialize(fm);
  // Return fm
  return fm;
#endif
}
GT_INLINE gt_fm* gt_fm_open_bzFILE(FILE* const stream,const gt_fm_mode mode) {
#ifndef HAVE_BZLIB
  gt_fatal_error(FM_NO_BZLIB_SUPPORT);
  return NULL;
#else
  GT_NULL_CHECK(stream);
  // Allocate handler
  gt_fm* fm = gt_alloc(gt_fm);
  // File
  fm->fd = 0;
  fm->file = stream;
  fm->gz_file = NULL;
  fm->bz_file = NULL;
  // Attributes
  fm->mode = mode;
  fm->file_name = GT_STREAM_FILE_NAME;
  fm->file_type = GT_FM_BZIPPED_FILE;
  fm->file_size = UINT64_MAX;
  // Initialize file manager
  gt_fm_initialize(fm);
  // Return fm
  return fm;
#endif
}
GT_INLINE void gt_fm_close(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  // Close fm
  switch (fm->file_type) {
    case GT_FM_STREAM:
      gt_cond_fatal_error__perror(fclose(fm->file),FM_CLOSE,fm->file_name);
      break;
    case GT_FM_REGULAR_FILE:
      gt_cond_fatal_error__perror(fclose(fm->file),FM_CLOSE,fm->file_name);
      break;
#ifdef HAVE_ZLIB
    case GT_FM_GZIPPED_FILE:
      gt_cond_fatal_error(gzclose(fm->gz_file),FM_GZCLOSE,fm->file_name);
      break;
#endif
#ifdef HAVE_BZLIB
    case GT_FM_BZIPPED_FILE: {
      int bzerr;
      if (fm->mode==GT_FM_READ) {
        BZ2_bzReadClose(&bzerr,fm->bz_file);
      } else {
        BZ2_bzWriteClose(&bzerr,fm->bz_file,0,NULL,NULL);
      }
      gt_cond_fatal_error(bzerr!=BZ_OK,FM_BZCLOSE,fm->file_name);
      break;
    }
#endif
    default:
      GT_INVALID_CASE();
      break;
  }
  // Free handler
  if (fm->skip_read_buffer!=NULL)  gt_free(fm->skip_read_buffer);
  if (fm->skip_write_buffer!=NULL) gt_free(fm->skip_write_buffer);
  gt_free(fm);
}
/*
 * Accesors
 */
GT_INLINE uint64_t gt_fm_get_current_position(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  return fm->byte_position;
}
GT_INLINE bool gt_fm_eof(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  // Shortcut eof
  if (fm->eof) return true;
  // Refresh eof
  switch (fm->file_type) {
    case GT_FM_STREAM:
    case GT_FM_REGULAR_FILE:
      fm->eof = feof(fm->file);
      return fm->eof;
      break;
#ifdef HAVE_ZLIB
    case GT_FM_GZIPPED_FILE:
      fm->eof = gzeof(fm->gz_file);
      return fm->eof;
      break;
#endif
#ifdef HAVE_BZLIB
    case GT_FM_BZIPPED_FILE:
      return fm->eof;
      break;
#endif
    default:
      GT_INVALID_CASE();
      break;
  }
}
GT_INLINE char* gt_fm_get_file_name(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  return fm->file_name;
}
/*
 * Seek
 */
GT_INLINE void gt_fm_skip_forward(gt_fm* const fm,const uint64_t num_bytes) {
  GT_FM_CHECK(fm);
  GT_ZERO_CHECK(num_bytes);
  if (fm->mode==GT_FM_READ) {
    /*
     * GT_FM_READ
     */
    if (fm->file_type==GT_FM_STREAM || fm->file_type==GT_FM_REGULAR_FILE) {
      // True Skip (if possible)
      gt_cond_fatal_error__perror(fseek(fm->file,num_bytes,SEEK_CUR),FM_SEEK,
          fm->file_name,fm->byte_position+num_bytes);
      // Update locator
      fm->byte_position += num_bytes;
    } else {
      // Skip By Reading
      if (fm->skip_read_buffer==NULL) {
        fm->skip_read_buffer = gt_malloc_(1,GT_FM_SKIP_BUFFER_SIZE,false,0);
      }
      uint64_t total_bytes_read = 0, bytes_read;
      while (total_bytes_read < num_bytes) {
        const uint64_t bytes_to_read = (total_bytes_read+GT_FM_SKIP_BUFFER_SIZE < num_bytes) ?
            GT_FM_SKIP_BUFFER_SIZE : num_bytes-total_bytes_read;
        bytes_read = gt_fm_read_mem(fm,fm->skip_read_buffer,bytes_to_read);
        total_bytes_read += bytes_to_read;
        if (bytes_read != bytes_to_read) break;
      }
      // Update locator
      fm->byte_position += total_bytes_read;
    }
  } else {
    /*
     * GT_FM_WRITE
     */
    // Skip Writing Zeros
    if (fm->skip_write_buffer==NULL) {
      fm->skip_write_buffer = gt_malloc_(1,GT_FM_SKIP_BUFFER_SIZE,true,0);
    }
    uint64_t total_bytes_written = 0;
    while (total_bytes_written < num_bytes) {
      const uint64_t bytes_to_write = (total_bytes_written+GT_FM_SKIP_BUFFER_SIZE < num_bytes) ?
          GT_FM_SKIP_BUFFER_SIZE : num_bytes-total_bytes_written;
      gt_fm_write_mem(fm,fm->skip_write_buffer,bytes_to_write);
      total_bytes_written += bytes_to_write;
    }
    // Update locator
    fm->byte_position += num_bytes;
  }
}
GT_INLINE void gt_fm_skip_uint64(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  if (fm->mode==GT_FM_READ) {
    gt_fm_read_uint64(fm);
  } else {
    gt_fm_write_uint64(fm,0);
  }
}
GT_INLINE void gt_fm_skip_uint32(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  if (fm->mode==GT_FM_READ) {
    gt_fm_read_uint32(fm);
  } else {
    gt_fm_write_uint32(fm,0);
  }
}
GT_INLINE void gt_fm_skip_uint16(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  if (fm->mode==GT_FM_READ) {
    gt_fm_read_uint16(fm);
  } else {
    gt_fm_write_uint16(fm,0);
  }
}
GT_INLINE void gt_fm_skip_uint8(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  if (fm->mode==GT_FM_READ) {
    gt_fm_read_uint8(fm);
  } else {
    gt_fm_write_uint8(fm,0);
  }
}
GT_INLINE void gt_fm_skip_align(gt_fm* const fm,const uint64_t num_bytes) {
  GT_FM_CHECK(fm);
  GT_ZERO_CHECK(num_bytes);
  const uint64_t bytes_to_skip = fm->byte_position + (num_bytes-1) % num_bytes;
  if (bytes_to_skip > 0) gt_fm_skip_forward(fm,bytes_to_skip);
}
GT_INLINE void gt_fm_skip_align_16(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  gt_fm_skip_align(fm,2);
}
GT_INLINE void gt_fm_skip_align_32(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  gt_fm_skip_align(fm,4);
}
GT_INLINE void gt_fm_skip_align_64(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  gt_fm_skip_align(fm,8);
}
GT_INLINE void gt_fm_skip_align_128(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  gt_fm_skip_align(fm,16);
}
GT_INLINE void gt_fm_skip_align_512(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  gt_fm_skip_align(fm,64);
}
GT_INLINE void gt_fm_skip_align_1024(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  gt_fm_skip_align(fm,128);
}
GT_INLINE void gt_fm_skip_align_4KB(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  gt_fm_skip_align(fm,512);
}
GT_INLINE void gt_fm_skip_align_mempage(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  int64_t sz = sysconf(_SC_PAGESIZE);
  gt_cond_fatal_error__perror(sz==-1,SYS_SYSCONF);
  gt_fm_skip_align(fm,sz);
}
/*
 * Read
 */
GT_INLINE uint64_t gt_fm_read_uint64(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  uint64_t var;
  gt_fm_read_mem(fm,&var,8);
  return var;
}
GT_INLINE uint32_t gt_fm_read_uint32(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  uint32_t var;
  gt_fm_read_mem(fm,&var,4);
  return var;
}
GT_INLINE uint16_t gt_fm_read_uint16(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  uint16_t var;
  gt_fm_read_mem(fm,&var,2);
  return var;
}
GT_INLINE uint8_t gt_fm_read_uint8(gt_fm* const fm) {
  GT_FM_CHECK(fm);
  uint8_t var;
  gt_fm_read_mem(fm,&var,1);
  return var;
}
GT_INLINE uint64_t gt_fm_read_mem(gt_fm* const fm,void* const dst,const uint64_t num_bytes) {
  GT_FM_CHECK(fm);
  GT_NULL_CHECK(dst);
  GT_ZERO_CHECK(num_bytes);
  gt_fatal_check(fm->mode!=GT_FM_READ,FM_INVALID_MODE_READ,fm->file_name);
  // Read
  int64_t num_bytes_read = 0;
  switch (fm->file_type) {
    case GT_FM_STREAM:
    case GT_FM_REGULAR_FILE:
      num_bytes_read = fread(dst,1,num_bytes,fm->file);
      fm->eof = (num_bytes_read<num_bytes);
      break;
#ifdef HAVE_ZLIB
    case GT_FM_GZIPPED_FILE:
      num_bytes_read = gzread(fm->gz_file,dst,num_bytes);
      gt_cond_fatal_error(num_bytes_read==-1,FM_GZREAD,num_bytes,fm->file_name);
      fm->eof = (num_bytes_read<num_bytes);
      break;
#endif
#ifdef HAVE_BZLIB
    case GT_FM_BZIPPED_FILE:
      if (!fm->eof) {
        int bzerror;
        num_bytes_read = BZ2_bzRead(&bzerror,fm->bz_file,dst,num_bytes);
        if (gt_expect_false(bzerror==BZ_STREAM_END)) {
          fm->eof = true;
        } else {
          gt_cond_fatal_error(bzerror!=BZ_OK,FM_BZREAD,num_bytes,fm->file_name);
          fm->eof = num_bytes_read<num_bytes;
        }
      }
      break;
#endif
    default:
      GT_INVALID_CASE();
      break;
  }
  // Update locator
  fm->byte_position += num_bytes_read;
  return num_bytes_read;
}
GT_INLINE uint64_t gt_fm_read_mem_parallel(gt_fm* const fm,void* const dst,const uint64_t num_bytes,const uint64_t num_threads) {
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
}
/*
 * Write
 */
GT_INLINE void gt_fm_write_uint64(gt_fm* const fm,const uint64_t data) {
  GT_FM_CHECK(fm);
  gt_fm_write_mem(fm,&data,8);
}
GT_INLINE void gt_fm_write_uint32(gt_fm* const fm,const uint32_t data) {
  GT_FM_CHECK(fm);
  gt_fm_write_mem(fm,&data,4);
}
GT_INLINE void gt_fm_write_uint16(gt_fm* const fm,const uint16_t data) {
  GT_FM_CHECK(fm);
  gt_fm_write_mem(fm,&data,2);
}
GT_INLINE void gt_fm_write_uint8(gt_fm* const fm,const uint8_t data) {
  GT_FM_CHECK(fm);
  gt_fm_write_mem(fm,&data,1);
}
GT_INLINE void gt_fm_write_mem(gt_fm* const fm,const void* const src,const uint64_t num_bytes) {
  GT_FM_CHECK(fm);
  gt_fatal_check(fm->mode!=GT_FM_WRITE,FM_INVALID_MODE_WRITE,fm->file_name);
  switch (fm->file_type) {
    case GT_FM_STREAM:
    case GT_FM_REGULAR_FILE:
      gt_cond_fatal_error__perror(fwrite(src,1,num_bytes,fm->file)!=num_bytes,FM_WRITE,fm->file_name);
      break;
#ifdef HAVE_ZLIB
    case GT_FM_GZIPPED_FILE: {
      gt_cond_fatal_error(gzwrite(fm->gz_file,src,num_bytes)!=num_bytes,FM_GZWRITE,fm->file_name);
      break;
    }
#endif
#ifdef HAVE_BZLIB
    case GT_FM_BZIPPED_FILE: {
      int bzerror;
      BZ2_bzWrite(&bzerror,fm->bz_file,(void*)src,num_bytes);
      gt_cond_fatal_error(bzerror!=BZ_OK,FM_BZWRITE,fm->file_name);
      break;
    }
#endif
    default:
      GT_INVALID_CASE();
      break;
  }
  // Update locator
  fm->byte_position += num_bytes;
}
/*
 * Bulk read of file
 */
GT_INLINE void gt_fm_bulk_read_fd(const int fd,void* const dst,const uint64_t size) {
  GT_NULL_CHECK(dst);
  uint64_t bytes_written = 0;
  while (bytes_written < size) {
    const uint64_t bytes_pending = size-bytes_written;
    const uint64_t chunk_size = (bytes_pending<GT_FM_BULK_COPY_BLOCK_SIZE) ? bytes_pending : GT_FM_BULK_COPY_BLOCK_SIZE;
    // Copy chunk
    gt_cond_fatal_error__perror(read(fd,dst+bytes_written,chunk_size)!=chunk_size,FM_FDREAD,chunk_size,fd);
    bytes_written += chunk_size;
  }
}
GT_INLINE void gt_fm_bulk_read_file(char* const file_name,void* const dst,const uint64_t offset,const uint64_t size) {
  GT_NULL_CHECK(file_name);
  GT_NULL_CHECK(dst);
  // Retrieve input file info
  struct stat stat_info;
  gt_cond_fatal_error__perror(stat(file_name,&stat_info)==-1,FM_STAT,file_name);
  // Open file descriptor
  const int fd = open(file_name,O_RDONLY,S_IRUSR);
  gt_cond_fatal_error__perror(fd==-1,FM_FDOPEN,file_name);
  if (offset > 0) {
    gt_cond_fatal_error__perror(lseek(fd,offset,SEEK_SET)==-1,FM_SEEK,file_name,offset); // Seek
  }
  // Read file
  gt_fm_bulk_read_fd(fd,dst,(size==0) ? stat_info.st_size-offset : size);
}
GT_INLINE void gt_fm_bulk_read_file_parallel(
    char* const file_name,void* const dst,const uint64_t offset,const uint64_t size,const uint64_t num_threads) {
//  GT_NULL_CHECK(file_name);
//  GT_NULL_CHECK(dst);
//  // Retrieve input file info
//  struct stat stat_info;
//  gt_cond_fatal_error__perror(stat(file_name,&stat_info)==-1,FM_STAT,file_name);
//  // Calculate size of each chunk
//  const uint64_t size_to_read = (size==0) ? stat_info.st_size : size;
//  const uint64_t chunk_size = size_to_read/num_threads; // FIXME: stat_info.st_size > num_threads
//  //#pragma omp parallel num_threads(num_threads) // FIXME: Pthreads
//  {
//    // Calculate offsets
//    const uint64_t tid = omp_get_thread_num();
//    const uint64_t thread_mem_offset = tid*chunk_size;
//    const uint64_t thread_file_offset = offset + thread_mem_offset;
//    const uint64_t thread_size = (tid < (num_threads-1)) ? chunk_size : stat_info.st_size-chunk_size*tid;
//    // Open file descriptor
//    const int fd = open(file_name,O_RDONLY,S_IRUSR);
//    gt_cond_fatal_error__perror(fd==-1,FM_FDOPEN,file_name);
//    gt_cond_fatal_error__perror(lseek(fd,thread_file_offset,SEEK_SET)==-1,FM_SEEK,file_name,thread_file_offset); // Seek
//    // Copy file chunk
//    fm_bulk_read_fd(fd,dst+thread_mem_offset,thread_size);
//  }
}

/*
 * FileManager Wrappers
 */
GT_INLINE void gt_stat(char* const file_name,struct stat *stat_info) {
  GT_NULL_CHECK(file_name);
  gt_cond_fatal_error__perror(stat(file_name,stat_info)==-1,FM_STAT,file_name);
}
GT_INLINE int gt_open_fd(char* const file_name,const int flags,const mode_t mode) {
  int fd = open(file_name,flags,mode);
  gt_cond_fatal_error__perror(fd==-1,FM_OPEN,file_name);
  return fd;
}
GT_INLINE FILE* gt_open_FILE(char* const file_name,const char* opentype) {
  FILE* const file = fopen(file_name,opentype);
  gt_cond_fatal_error__perror(file==NULL,FM_FDOPEN,file_name);
  return file;
}
/*
 * FileManager Printers
 */
GT_INLINE gt_status gt_vfmprintf(gt_fm* const fm,const char *template,va_list v_args) {
  GT_FM_CHECK(fm);
  GT_NULL_CHECK(template);
  GT_NULL_CHECK(v_args);
  // Depending on the type we might have to compose the data into a buffer
  gt_status num_bytes = 0;
  switch (fm->file_type) {
    case GT_FM_STREAM:
    case GT_FM_REGULAR_FILE:
      num_bytes = vfprintf(fm->file,template,v_args);
      break;
    case GT_FM_GZIPPED_FILE:
    case GT_FM_BZIPPED_FILE:
      if (fm->skip_write_buffer==NULL) {
        fm->skip_write_buffer = gt_malloc_(1,GT_FM_SKIP_BUFFER_SIZE,true,0);
      }
      num_bytes = vsprintf((char*)fm->skip_write_buffer,template,v_args);
      gt_fm_write_mem(fm,fm->skip_write_buffer,num_bytes);
      break;
    default:
      GT_INVALID_CASE();
      break;
  }
  return num_bytes;
}
GT_INLINE gt_status gt_fmprintf(gt_fm* const fm,const char *template,...) {
  GT_FM_CHECK(fm);
  GT_NULL_CHECK(template);
  va_list v_args;
  va_start(v_args,template);
  const gt_status error_code = gt_vfmprintf(fm,template,v_args);
  va_end(v_args);
  return error_code;
}
