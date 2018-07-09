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
 * DESCRIPTION:
 *   Provides functionality for an easy file management
 *   providing stream-based access for reading & writing
 */

#ifndef FILE_MANAGEMENT_H_
#define FILE_MANAGEMENT_H_

#include "system/commons.h"
#include "system/mm.h"

/*
 * I/O Constants/Values
 */
extern const int fm_open_flags[3];
extern const mode_t fm_open_mode[3];

typedef enum { FM_READ=0, FM_WRITE=1, FM_READ_WRITE=2 } fm_mode;
typedef enum { FM_STREAM, FM_REGULAR_FILE, FM_GZIPPED_FILE, FM_BZIPPED_FILE, FM_POPEN } fm_type;
typedef struct _fm_t fm_t;

/*
 * Setup
 */
fm_t* fm_open_file(char* const file_name,const fm_mode mode);
fm_t* fm_open_popen(char* const file_name,const fm_mode mode);
fm_t* fm_open_FILE(FILE* const stream,const fm_mode mode);
fm_t* fm_open_gzFILE(FILE* const stream,const fm_mode mode);
fm_t* fm_open_bzFILE(FILE* const stream,const fm_mode mode);
fm_t* fm_open_temp_file(void);
void fm_close(fm_t* const file_manager);

/*
 * Accessors
 */
uint64_t fm_get_current_position(fm_t* const file_manager);
bool fm_eof(fm_t* const file_manager);
char* fm_get_file_name(fm_t* const file_manager);
uint64_t fm_get_file_size(fm_t* const file_manager);

/*
 * Seek
 */
void fm_seek(fm_t* const file_manager,const uint64_t position);

void fm_skip_forward(fm_t* const file_manager,const uint64_t num_bytes);
void fm_skip_uint64(fm_t* const file_manager);
void fm_skip_uint32(fm_t* const file_manager);
void fm_skip_uint16(fm_t* const file_manager);
void fm_skip_uint8(fm_t* const file_manager);
void fm_skip_align(fm_t* const file_manager,const uint64_t num_bytes);
void fm_skip_align_16(fm_t* const file_manager);
void fm_skip_align_32(fm_t* const file_manager);
void fm_skip_align_64(fm_t* const file_manager);
void fm_skip_align_128(fm_t* const file_manager);
void fm_skip_align_512(fm_t* const file_manager);
void fm_skip_align_1024(fm_t* const file_manager);
void fm_skip_align_4KB(fm_t* const file_manager);
void fm_skip_align_mempage(fm_t* const file_manager);

/*
 * Read
 */
#define fm_read(file_manager,var) fm_copy_mem(file_manager,&var,sizeof(var))
uint64_t fm_read_uint64(fm_t* const file_manager);
uint32_t fm_read_uint32(fm_t* const file_manager);
uint16_t fm_read_uint16(fm_t* const file_manager);
uint8_t fm_read_uint8(fm_t* const file_manager);
ssize_t fm_getline(char **buf, size_t *bufsiz, fm_t* const file_manager);
uint64_t fm_read_mem(
    fm_t* const file_manager,
    void* const dst,
    const uint64_t num_bytes);

// TODO uint64_t fm_read_mem_parallel(
//    fm_t* const file_manager,
//    void* const dst,
//    const uint64_t num_bytes,
//    const uint64_t num_threads);

mm_t* fm_load_mem(fm_t* const file_manager,const uint64_t num_bytes);

void fm_prefetch_next(fm_t* const file_manager,const uint64_t num_bytes);

/*
 * Write
 */
#define fm_write(file_manager,var) fm_write_mem(file_manager,&var,sizeof(var))
void fm_write_uint64(fm_t* const file_manager,const uint64_t data);
void fm_write_uint32(fm_t* const file_manager,const uint32_t data);
void fm_write_uint16(fm_t* const file_manager,const uint16_t data);
void fm_write_uint8(fm_t* const file_manager,const uint8_t data);
void fm_write_mem(fm_t* const file_manager,const void* const src,const uint64_t num_bytes);

/*
 * Bulk Read of file
 */
void fm_bulk_read_fd(
    const int fd,
    void* const dst,
    const uint64_t size);
void fm_bulk_read_file(
    char* const file_name,
    void* const dst,
    const uint64_t offset,
    const uint64_t size);
// TODO void fm_bulk_read_file_parallel(
//    char* const file_name,
//    void* const dst,
//    const uint64_t offset,
//    const uint64_t size,
//    const uint64_t num_threads);

/*
 * FileManager Wrappers
 */
void gem_stat(char* const file_name,struct stat *stat_info);
int gem_open_fd(char* const file_name,const int flags,const mode_t mode);
FILE* gem_open_FILE(char* const file_name,const char* opentype);
void gem_unlink(char* const file_name);

/*
 * Utils
 */
bool gem_access(char* const path,const fm_mode mode);
uint64_t gem_file_size(const char* const file_name);

/*
 * FileManager Printers
 */
int vfmprintf(fm_t* const file_manager,const char *template,va_list v_args);
int fmprintf(fm_t* const file_manager,const char *template,...);

#endif /* FILE_MANAGEMENT_H_ */

