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
 *   Memory Manager provides memory allocation functions. Different types of memory are supported.
 *     - UnitMemory
 *         Allocate relative small chunks of memory relying on the regular memory manager,
 *         usually malloc/calloc using a BuddySystem (Helper functions)
 *     - BulkMemory
 *         Allocate big chunks of memory and resort to disk if memory is not enough
 *  CONSIDERATIONS OF PERFORMANCE:
 *    1.- Use of other underlying memory managers for unit memory
 *      TCMalloc (Thread-Caching Malloc)/nedmalloc
 *    2.- Use of hints ...
 *      madvise() / readahead() / posix_fadvise()
 */

#include "system/mm.h"
#include "system/fm.h"
#include <sys/resource.h>

/*
 * Config
 */
#ifndef MAP_HUGETLB
  #define MAP_HUGETLB 0 // In some environments MAP_HUGETLB can be undefined
#endif
#ifndef MAP_ANONYMOUS
  #define MAP_ANONYMOUS 0 // Disabled for Mac compatibility
#endif
#ifndef MAP_POPULATE
  #define MAP_POPULATE 0 // Disabled for Mac compatibility
#endif

/*
 * Error Messages
 */
#define GEM_ERROR_MEM_ALLOC "Could not allocate memory (%"PRIu64"Bytes requested)"
#define GEM_ERROR_MEM_CALLOC "Could not allocate memory (%"PRIu64"Bytes x %"PRIu64" requested)"
#define GEM_ERROR_MEM_REALLOC "Could not re-allocate memory (%"PRIu64" requested)"
#define GEM_ERROR_MEM_ALLOC_DISK "Requested %"PRIu64" Bytes. Resorting to disk"
#define GEM_ERROR_MEM_ALLOC_MMAP_DISK_FAIL "Requested %"PRIu64" Bytes. Failed to mmap memory to '%s'"
#define GEM_ERROR_MEM_ALLOC_MMAP_FAIL "Requested %"PRIu64" Bytes. Failed to mmap memory"
#define GEM_ERROR_MEM_CURSOR_OUT_OF_SEGMENT "Current memory cursor is out of boundaries (Segmentation fault)"
#define GEM_ERROR_MEM_CURSOR_SEEK "Could not seek to address %"PRIu64". Out of boundaries (Segmentation fault)"
#define GEM_ERROR_MEM_MAP_ZERO_FILE "Mapping zero-bytes File '%s' to memory"
#define GEM_ERROR_MEM_ALG_FAILED "Failed aligning the memory address to the specified boundary"
#define GEM_ERROR_MEM_STAT_MEMINFO "Could not stat '/proc/meminfo' (%s)"
#define GEM_ERROR_MEM_SYSINFO "Error calling sysinfo"
#define GEM_ERROR_MEM_PARSE_STATM "Could not parse stat-memory information"

/*
 * Memory Alignment Utils
 */
const uint64_t mm_mem_alignment_bits_mask[] = { // Check Memory Alignment Bits (Masks)
    0x0000000000000001lu, /*  [0]    16 bits aligned (   2B  @2^1)  */
    0x0000000000000003lu, /*  [1]    32 bits aligned (   4B  @2^2)  */
    0x0000000000000007lu, /*  [2]    64 bits aligned (   8B  @2^3)  */
    0x000000000000000Flu, /*  [3]   128 bits aligned (  16B  @2^4)  */
    0x000000000000001Flu, /*  [4]   256 bits aligned (  32B  @2^5)  */
    0x000000000000003Flu, /*  [5]   512 bits aligned (  64B  @2^6)  */
    0x000000000000007Flu, /*  [6]  1024 bits aligned ( 128B  @2^7)  */
    0x00000000000000FFlu, /*  [7]  2048 bits aligned ( 256B  @2^8)  */
    0x00000000000001FFlu, /*  [8]  4096 bits aligned ( 512B  @2^9)  */
    0x00000000000003FFlu, /*  [9]  8192 bits aligned (   1KB @2^10) */
    0x00000000000007FFlu, /* [10] 16384 bits aligned (   2KB @2^11) */
    0x0000000000000FFFlu, /* [11] 32768 bits aligned (   4KB @2^12) */ /* RegularPage Size */
    0x0000000000001FFFlu, /* [12] ***** bits aligned (   8KB @2^13) */
    0x0000000000003FFFlu, /* [13] ***** bits aligned (  16KB @2^14) */
    0x0000000000007FFFlu, /* [14] ***** bits aligned (  32KB @2^15) */
    0x000000000000FFFFlu, /* [15] ***** bits aligned (  64KB @2^16) */
    0x000000000001FFFFlu, /* [16] ***** bits aligned ( 128KB @2^17) */
    0x000000000003FFFFlu, /* [17] ***** bits aligned ( 256KB @2^18) */
    0x000000000007FFFFlu, /* [18] ***** bits aligned ( 512KB @2^19) */
    0x00000000000FFFFFlu, /* [19] ***** bits aligned (   1MB @2^20) */
    0x00000000001FFFFFlu, /* [20] ***** bits aligned (   2MB @2^21) */ /* RegularPageHugeTLB Size */
    0x00000000003FFFFFlu, /* [21] ***** bits aligned (   4MB @2^22) */
    0x00000000007FFFFFlu, /* [22] ***** bits aligned (   8MB @2^23) */
};

/*
 * MMap Constants/Values
 */
int mm_proc_flags[3] = { PROT_READ, PROT_READ|PROT_WRITE, PROT_READ|PROT_WRITE };
int mm_mmap_mode[3] = { MAP_PRIVATE, MAP_SHARED, MAP_SHARED };
/*
 * Temporal folder path
 */
char* mm_temp_folder_path = MM_DEFAULT_TMP_FOLDER;

char* mm_get_tmp_folder(void) {
  return mm_temp_folder_path;
}
void mm_set_tmp_folder(char* const tmp_folder_path) {
  GEM_CHECK_NULL(tmp_folder_path);
  mm_temp_folder_path = tmp_folder_path;
}
/*
 * UnitMemory
 *   Allocate relative small chunks of memory relying on the regular memory manager,
 *   usually malloc/calloc using a BuddySystem (Helper functions)
 */
void* mm_malloc_nothrow(uint64_t const num_elements,const uint64_t size_element,const bool init_mem,const int init_value) {
  // Request memory
  const uint64_t total_memory = num_elements*size_element;
  void* const allocated_mem = (gem_expect_false(init_mem && init_value==0)) ?
      calloc(num_elements,size_element) : malloc(total_memory);
  // Check memory
  if (!allocated_mem) return NULL;
  // Initialize memory
  if (gem_expect_false(init_mem && init_value!=0)) memset(allocated_mem,init_value,total_memory);
  // MM_PRINT_MEM_ALIGMENT(allocated_mem); // Debug
  return allocated_mem;
}
void* mm_malloc_(uint64_t const num_elements,const uint64_t size_element,const bool init_mem,const int init_value) {
  // Request memory
  const uint64_t total_memory = num_elements*size_element;
  void* allocated_mem;
  if (gem_expect_false(init_mem && init_value==0)) {
    allocated_mem = calloc(num_elements,size_element);
    gem_cond_fatal_error(!allocated_mem,MEM_CALLOC,num_elements,size_element);
  } else {
    allocated_mem = malloc(total_memory);
    gem_cond_fatal_error(!allocated_mem,MEM_ALLOC,total_memory);
  }
  // Initialize memory
  if (gem_expect_false(init_mem && init_value!=0)) memset(allocated_mem,init_value,total_memory);
  // MM_PRINT_MEM_ALIGMENT(allocated_mem); // Debug
  return allocated_mem;
}
void* mm_realloc_nothrow(void* mem_addr,const uint64_t num_bytes) {
  GEM_CHECK_NULL(mem_addr);
  GEM_CHECK_ZERO(num_bytes);
  return realloc(mem_addr,num_bytes);
}
void* mm_realloc(void* mem_addr,const uint64_t num_bytes) {
  GEM_CHECK_NULL(mem_addr);
  GEM_CHECK_ZERO(num_bytes);
  void* const new_mem_addr = mm_realloc_nothrow(mem_addr,num_bytes);
  gem_cond_fatal_error(!new_mem_addr,MEM_REALLOC,num_bytes);
  return new_mem_addr;
}
void mm_free(void* mem_addr) {
  GEM_CHECK_NULL(mem_addr);
  free(mem_addr);
}
/*
 * BulkMemory
 *   Allocate big chunks of memory and resort to disk if memory is not enough
 */
mm_t* mm_bulk_malloc(const uint64_t num_bytes,const bool init_mem) {
  GEM_CHECK_ZERO(num_bytes);
  void* memory = mm_malloc_nothrow(num_bytes,1,init_mem,0);
  if (gem_expect_true(memory!=NULL)) { // Fits in HEAP
    mm_t* const mem_manager = mm_alloc(mm_t);
    mem_manager->memory = memory;
    mem_manager->mem_type = MM_HEAP;
    mem_manager->mode = MM_READ_WRITE;
    mem_manager->allocated = num_bytes;
    mem_manager->cursor = mem_manager->memory;
    // MM_PRINT_MEM_ALIGMENT(mem_manager->memory); // Debug
    return mem_manager;
  } else { // Resort to MMAP in disk
    gem_warn(MEM_ALLOC_DISK,num_bytes);
    return mm_bulk_mmalloc_temp(num_bytes);
  }
}
mm_t* mm_bulk_mmalloc(const uint64_t num_bytes,const bool use_huge_pages) {
  GEM_CHECK_ZERO(num_bytes);
  // Allocate handler
  mm_t* const mem_manager = mm_alloc(mm_t);
#ifdef MM_NO_MMAP
  mem_manager->memory = mm_malloc_nothrow(num_bytes,1,0,0);
  mem_manager->cursor = mem_manager->memory;
  mem_manager->mem_type = MM_HEAP;
  mem_manager->mode = MM_READ_WRITE;
  mem_manager->allocated = num_bytes;
  mem_manager->fd = -1;
  mem_manager->file_name = NULL;
  return mem_manager;
#else
  /*
   * MMap memory (anonymous)
   *   - MAP_PRIVATE => Fits in RAM+SWAP
   *   - MAP_ANONYMOUS => The mapping is not backed by any file; its contents are initialized to zero.
   *       Map against /dev/zero (Allocate anonymous memory segment, without open)
   *   - MAP_NORESERVE to explicitly enable swap space overcommitting. (echo 1 > /proc/sys/vm/overcommit_memory)
   *     Useful when you wish to map a file larger than the amount of free memory
   *     available on your system (RAM+SWAP).
   *       In this case, the lazy swap space reservation may cause the program
   *       to consume all the free RAM and swap on the system, eventually
   *       triggering the OOM killer (Linux) or causing a SIGSEGV.
   */
  int flags = MAP_PRIVATE | MAP_ANONYMOUS; // MAP_NORESERVE (Seems that MAP_NORESERVE & MAP_HUGETLB gives a problem)
  if (use_huge_pages) flags |= MAP_HUGETLB;
  mem_manager->memory = mmap(0,num_bytes,PROT_READ|PROT_WRITE,flags,-1,0);
  gem_cond_fatal_error(mem_manager->memory==MAP_FAILED,MEM_ALLOC_MMAP_FAIL,num_bytes);
  mem_manager->cursor = mem_manager->memory;
  // Set MM
  mem_manager->mem_type = MM_MMAPPED;
  mem_manager->mode = MM_READ_WRITE;
  mem_manager->allocated = num_bytes;
  mem_manager->fd = -1;
  mem_manager->file_name = NULL;
  // MM_PRINT_MEM_ALIGMENT(mem_manager->memory); // Debug
  return mem_manager;
#endif
}
mm_t* mm_bulk_mmalloc_temp(const uint64_t num_bytes) {
  GEM_CHECK_ZERO(num_bytes);
#ifdef MM_NO_MMAP
  return mm_bulk_mmalloc(num_bytes,false);
#else
  // Allocate handler
  mm_t* const mem_manager = mm_alloc(mm_t);
  // TemporalMemory (backed by a file)
  mem_manager->file_name = mm_calloc(strlen(mm_get_tmp_folder())+22,char,true);
  sprintf(mem_manager->file_name,"%smm_temp_XXXXXX",mm_get_tmp_folder());
  // Create temporary file
  mem_manager->fd = mkstemp(mem_manager->file_name);
  gem_cond_fatal_error(mem_manager->fd==-1,SYS_MKSTEMP,mem_manager->file_name);
  gem_cond_fatal_error(unlink(mem_manager->file_name),SYS_HANDLE_TMP); // Make it temporary
  gem_log("Allocating memory mapped to disk: %s (%"PRIu64" MBytes) [PhysicalMem Available %"PRIu64" MBytes]",
      mem_manager->file_name,num_bytes/1024/1024,mm_get_mem_available_total()/1024/1024);
  // Set the size of the temporary file (disk allocation)
  gem_cond_fatal_error(lseek(mem_manager->fd,num_bytes-1,SEEK_SET)==-1,SYS_HANDLE_TMP);
  gem_cond_fatal_error(write(mem_manager->fd,"",1)<=0,SYS_HANDLE_TMP);
  gem_cond_fatal_error(lseek(mem_manager->fd,0,SEEK_SET)==-1,SYS_HANDLE_TMP);
  /*
   * Mmap file.
   *   - MAP_SHARED as we the mapping will be reflected on disk (no copy-on-write)
   *     As such, the kernel knows it can always free up memory by doing writeback.
   *   - MAP_NORESERVE to explicitly enable swap space overcommitting. (echo 1 > /proc/sys/vm/overcommit_memory)
   *     Useful when you wish to map a file larger than the amount of free memory
   *     available on your system (RAM+SWAP).
   *       In this case, the lazy swap space reservation may cause the program
   *       to consume all the free RAM and swap on the system, eventually
   *       triggering the OOM killer (Linux) or causing a SIGSEGV.
   */
  mem_manager->memory = mmap(NULL,num_bytes,PROT_READ|PROT_WRITE,MAP_SHARED|MAP_NORESERVE,mem_manager->fd,0);
  gem_cond_fatal_error(mem_manager->memory==MAP_FAILED,MEM_ALLOC_MMAP_DISK_FAIL,num_bytes,mem_manager->file_name);
  mem_manager->cursor = mem_manager->memory;
  // Set MM
  mem_manager->mem_type = MM_MMAPPED;
  mem_manager->mode = MM_READ_WRITE;
  mem_manager->allocated = num_bytes;
  // MM_PRINT_MEM_ALIGMENT(mem_manager->memory); // Debug
  return mem_manager;
#endif
}
void mm_bulk_free(mm_t* const mem_manager) {
  GEM_CHECK_NULL(mem_manager);
  GEM_CHECK_NULL(mem_manager->memory);
  GEM_CHECK_NULL(mem_manager->cursor);
  if (mem_manager->mem_type==MM_HEAP) { // Heap BulkMemory
    mm_free(mem_manager->memory);
  } else { // MMapped BulkMemory
#ifdef MM_NO_MMAP
    mm_free(mem_manager->memory);
#else
    munmap(mem_manager->memory,mem_manager->allocated);
    if (mem_manager->fd!=-1) {
      gem_cond_fatal_error(close(mem_manager->fd),SYS_HANDLE_TMP);
    }
#endif
  }
  CFREE(mem_manager->file_name);
  mm_free(mem_manager);
}
mm_t* mm_bulk_mmap_file(char* const file_name,const mm_mode mode,const bool populate_page_tables) {
  GEM_CHECK_NULL(file_name);
#ifdef MM_NO_MMAP
  return mm_bulk_mload_file(file_name);
#else
  // Allocate handler
  mm_t* const mem_manager = mm_alloc(mm_t);
  // Retrieve input file info
  struct stat stat_info;
  gem_stat(file_name,&stat_info);
  gem_cond_fatal_error(stat_info.st_size==0,MEM_MAP_ZERO_FILE,file_name);
  // Open file descriptor
  mem_manager->fd = gem_open_fd(file_name,fm_open_flags[mode],S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
  /*
   * Mmap file
   *   - @mode::
   *       MM_READ_ONLY => MAP_PRIVATE (no copy-on-write as it's not allowed)
   *       MM_WRITE_ONLY or MM_READ_WRITE => MAP_SHARED
   *   - MAP_POPULATE (since Linux 2.5.46)
   *       Populate (prefault) page tables for a mapping. For a file mapping, this causes
   *       read-ahead on the file. Later accesses to the mapping will not be blocked by page faults.
   *       MAP_POPULATE is only supported for private mappings since Linux 2.6.23.
   */
  int flags = mm_mmap_mode[mode];
  if (mode==MM_READ_ONLY && populate_page_tables) flags |= MAP_POPULATE;
  mem_manager->memory = mmap(0,stat_info.st_size,mm_proc_flags[mode],flags,mem_manager->fd,0);
  gem_cond_fatal_error(mem_manager->memory==MAP_FAILED,SYS_MMAP_FILE,file_name);
  mem_manager->cursor = mem_manager->memory;
  // Set MM
  mem_manager->mem_type = MM_MMAPPED;
  mem_manager->mode = mode;
  mem_manager->allocated = stat_info.st_size;
  mem_manager->file_name = strndup(file_name,strlen(file_name));
  // MM_PRINT_MEM_ALIGMENT(mem_manager->memory); // Debug
  return mem_manager;
#endif
}
mm_t* mm_bulk_load_file(char* const file_name) {
  GEM_CHECK_NULL(file_name);
#ifdef MM_NO_MMAP
  return mm_bulk_mload_file(file_name);
#else
  // Allocate handler
  mm_t* const mem_manager = mm_alloc(mm_t);
  // Retrieve input file info
  struct stat stat_info;
  gem_cond_fatal_error(stat(file_name,&stat_info)==-1,SYS_STAT,file_name);
  gem_cond_fatal_error(stat_info.st_size==0,MEM_MAP_ZERO_FILE,file_name);
  // Allocate memory to dump the content of the file
  mem_manager->memory = mm_malloc(stat_info.st_size);
  gem_cond_fatal_error(!mem_manager->memory,MEM_ALLOC,stat_info.st_size);
  mem_manager->mem_type = MM_HEAP;
  mem_manager->mode = MM_READ_ONLY;
  mem_manager->allocated = stat_info.st_size;
  mem_manager->cursor = mem_manager->memory;
  // MM_PRINT_MEM_ALIGMENT(mem_manager->memory); // Debug
  // Read the file and dump it into memory
  fm_bulk_read_file(file_name,mem_manager->memory,0,0);
  return mem_manager;
#endif
}
mm_t* mm_bulk_mload_file(char* const file_name) {
  GEM_CHECK_NULL(file_name);
  // Retrieve input file info
  struct stat stat_info;
  gem_cond_fatal_error(stat(file_name,&stat_info)==-1,SYS_STAT,file_name);
  gem_cond_fatal_error(stat_info.st_size==0,MEM_MAP_ZERO_FILE,file_name);
  // Allocate memory to dump the content of the file
  mm_t* const mem_manager = mm_bulk_mmalloc(stat_info.st_size,false);
  // Read the file and dump it into memory
  fm_bulk_read_file(file_name,mem_manager->memory,0,0);
  return mem_manager;
}

/*
 * Accessors
 */
void* mm_get_mem(mm_t* const mem_manager) {
  return mem_manager->cursor;
}
void* mm_get_base_mem(mm_t* const mem_manager) {
  return mem_manager->memory;
}
mm_mode mm_get_mode(mm_t* const mem_manager) {
  return mem_manager->mode;
}
uint64_t mm_get_allocated(mm_t* const mem_manager) {
  return mem_manager->allocated;
}
char* mm_get_mfile_name(mm_t* const mem_manager) {
  return mem_manager->file_name;
}
/*
 * Seek functions
 */
uint64_t mm_get_current_position(mm_t* const mem_manager) {
  return (mem_manager->cursor-mem_manager->memory);
}
bool mm_eom(mm_t* const mem_manager) {
  return mm_get_current_position(mem_manager) >= mem_manager->allocated;
}
void mm_seek(mm_t* const mem_manager,const uint64_t byte_position) {
  gem_fatal_check(byte_position>=mem_manager->allocated,MEM_CURSOR_SEEK,byte_position);
  mem_manager->cursor = mem_manager->memory + byte_position;
}
void mm_skip_forward(mm_t* const mem_manager,const uint64_t num_bytes) {
  mem_manager->cursor += num_bytes;
  MM_CHECK_SEGMENT(mem_manager);
}
void mm_skip_backward(mm_t* const mem_manager,const uint64_t num_bytes) {
  mem_manager->cursor -= num_bytes;
  MM_CHECK_SEGMENT(mem_manager);
}
void mm_skip_uint64(mm_t* const mem_manager) {
  mem_manager->cursor += 8;
  MM_CHECK_SEGMENT(mem_manager);
}
void mm_skip_uint32(mm_t* const mem_manager) {
  mem_manager->cursor += 4;
  MM_CHECK_SEGMENT(mem_manager);
}
void mm_skip_uint16(mm_t* const mem_manager) {
  mem_manager->cursor += 2;
  MM_CHECK_SEGMENT(mem_manager);
}
void mm_skip_uint8(mm_t* const mem_manager) {
  mem_manager->cursor += 1;
  MM_CHECK_SEGMENT(mem_manager);
}
#ifdef MM_NO_MMAP
void mm_skip_align(mm_t* const mem_manager,const uint64_t num_bytes) {
  GEM_CHECK_ZERO(num_bytes);
  if (gem_expect_true(num_bytes > 1)) {
    const uint64_t byte_position = mem_manager->cursor - mem_manager->memory;
    const uint64_t bytes_mod = byte_position % num_bytes;
    if (bytes_mod > 0) {
      const uint64_t bytes_to_skip = num_bytes - bytes_mod;
      mm_skip_forward(mem_manager,bytes_to_skip);
    }
  }
}
#else
void mm_skip_align(mm_t* const mem_manager,const uint64_t num_bytes) {
  GEM_CHECK_ZERO(num_bytes);
  if (gem_expect_true(num_bytes > 1)) {
    mem_manager->cursor = mem_manager->cursor+(num_bytes-1);
    mem_manager->cursor = mem_manager->cursor-(MM_CAST_ADDR(mem_manager->cursor)%num_bytes);
    MM_CHECK_SEGMENT(mem_manager);
    gem_fatal_check(MM_CAST_ADDR(mem_manager->cursor)%num_bytes!=0,MEM_ALG_FAILED);
  }
}
#endif
void mm_skip_align_16(mm_t* const mem_manager) {
#ifdef MM_NO_MMAP
  mm_skip_align(mem_manager,2);
#else
  mem_manager->cursor = MM_CAST_PTR(
      (MM_CAST_ADDR(mem_manager->cursor)+MM_MEM_ALIGNED_MASK_16b) & (~MM_MEM_ALIGNED_MASK_16b));
  MM_CHECK_SEGMENT(mem_manager);
  MM_CHECK_ALIGNMENT(mem_manager,16b);
#endif
}
void mm_skip_align_32(mm_t* const mem_manager) {
#ifdef MM_NO_MMAP
  mm_skip_align(mem_manager,4);
#else
  mem_manager->cursor = MM_CAST_PTR(
      (MM_CAST_ADDR(mem_manager->cursor)+MM_MEM_ALIGNED_MASK_32b) & (~MM_MEM_ALIGNED_MASK_32b));
  MM_CHECK_SEGMENT(mem_manager);
  MM_CHECK_ALIGNMENT(mem_manager,32b);
#endif
}
void mm_skip_align_64(mm_t* const mem_manager) {
#ifdef MM_NO_MMAP
  mm_skip_align(mem_manager,8);
#else
  mem_manager->cursor = MM_CAST_PTR(
      (MM_CAST_ADDR(mem_manager->cursor)+MM_MEM_ALIGNED_MASK_64b) & (~MM_MEM_ALIGNED_MASK_64b));
  MM_CHECK_SEGMENT(mem_manager);
  MM_CHECK_ALIGNMENT(mem_manager,64b);
#endif
}
void mm_skip_align_128(mm_t* const mem_manager) {
#ifdef MM_NO_MMAP
  mm_skip_align(mem_manager,16);
#else
  mem_manager->cursor = MM_CAST_PTR(
      (MM_CAST_ADDR(mem_manager->cursor)+MM_MEM_ALIGNED_MASK_128b) & (~MM_MEM_ALIGNED_MASK_128b));
  MM_CHECK_SEGMENT(mem_manager);
  MM_CHECK_ALIGNMENT(mem_manager,128b);
#endif
}
void mm_skip_align_512(mm_t* const mem_manager) {
#ifdef MM_NO_MMAP
  mm_skip_align(mem_manager,64);
#else
  mem_manager->cursor = MM_CAST_PTR(
      (MM_CAST_ADDR(mem_manager->cursor)+MM_MEM_ALIGNED_MASK_512b) & (~MM_MEM_ALIGNED_MASK_512b));
  MM_CHECK_SEGMENT(mem_manager);
  MM_CHECK_ALIGNMENT(mem_manager,512b);
#endif
}
void mm_skip_align_1024(mm_t* const mem_manager) {
#ifdef MM_NO_MMAP
  mm_skip_align(mem_manager,1024);
#else
  mem_manager->cursor = MM_CAST_PTR(
      (MM_CAST_ADDR(mem_manager->cursor)+MM_MEM_ALIGNED_MASK_1024b) & (~MM_MEM_ALIGNED_MASK_1024b));
  MM_CHECK_SEGMENT(mem_manager);
  MM_CHECK_ALIGNMENT(mem_manager,1KB);
#endif
}
void mm_skip_align_4KB(mm_t* const mem_manager) {
#ifdef MM_NO_MMAP
  mm_skip_align(mem_manager,4096);
#else
  mem_manager->cursor = MM_CAST_PTR(
      (MM_CAST_ADDR(mem_manager->cursor)+MM_MEM_ALIGNED_MASK_4KB) & (~MM_MEM_ALIGNED_MASK_4KB));
  MM_CHECK_SEGMENT(mem_manager);
  MM_CHECK_ALIGNMENT(mem_manager,4KB);
#endif
}
void mm_skip_align_mempage(mm_t* const mem_manager) {
  const int64_t page_size = mm_get_page_size();
  gem_cond_fatal_error(page_size==-1,SYS_SYSCONF);
  mm_skip_align(mem_manager,page_size);
}
/*
 * Read functions
 */
uint64_t mm_read_uint64(mm_t* const mem_manager) {
  MM_CHECK_SEGMENT(mem_manager);
  const uint64_t data = *((uint64_t*)mem_manager->cursor);
  mem_manager->cursor += 8;
  return data;
}
uint32_t mm_read_uint32(mm_t* const mem_manager) {
  MM_CHECK_SEGMENT(mem_manager);
  const uint32_t data = *((uint32_t*)mem_manager->cursor);
  mem_manager->cursor += 4;
  return data;
}
uint16_t mm_read_uint16(mm_t* const mem_manager) {
  MM_CHECK_SEGMENT(mem_manager);
  const uint16_t data = *((uint16_t*)mem_manager->cursor);
  mem_manager->cursor += 2;
  return data;
}
uint8_t mm_read_uint8(mm_t* const mem_manager) {
  MM_CHECK_SEGMENT(mem_manager);
  const uint8_t data = *((uint8_t*)mem_manager->cursor);
  mem_manager->cursor += 1;
  return data;
}
void* mm_read_mem(mm_t* const mem_manager,const uint64_t num_bytes) {
  MM_CHECK_SEGMENT(mem_manager);
  void* const current_cursor = mem_manager->cursor;
  mem_manager->cursor += num_bytes;
  return current_cursor;
}
void mm_copy_mem(
    mm_t* const mem_manager,
    void* const dst,
    const uint64_t num_bytes) {
  MM_CHECK_SEGMENT(mem_manager);
  memcpy(dst,mem_manager->cursor,num_bytes);
  mem_manager->cursor += num_bytes;
}
/*
 * Write functions
 */
void mm_write_uint64(mm_t* const mem_manager,const uint64_t data) {
  MM_CHECK_SEGMENT(mem_manager);
  *((uint64_t*)mem_manager->cursor) = data;
  mem_manager->cursor += 8;
}
void mm_write_uint32(mm_t* const mem_manager,const uint32_t data) {
  MM_CHECK_SEGMENT(mem_manager);
  *((uint32_t*)mem_manager->cursor) = data;
  mem_manager->cursor += 4;
}
void mm_write_uint16(mm_t* const mem_manager,const uint16_t data) {
  MM_CHECK_SEGMENT(mem_manager);
  *((uint16_t*)mem_manager->cursor) = data;
  mem_manager->cursor += 2;
}
void mm_write_uint8(mm_t* const mem_manager,const uint8_t data) {
  MM_CHECK_SEGMENT(mem_manager);
  *((uint8_t*)mem_manager->cursor) = data;
  mem_manager->cursor += 1;
}
void mm_write_mem(mm_t* const mem_manager,void* const src,const uint64_t num_bytes) {
  MM_CHECK_SEGMENT(mem_manager);
  memcpy(mem_manager->cursor,src,num_bytes);
  mem_manager->cursor += num_bytes;
}
/*
 * Status
 */
#ifdef __MACH__
#include <sys/sysctl.h>
int64_t mm_get_page_size(void) {
  return BUFFER_SIZE_4K;
}
int64_t mm_get_mem_available_virtual(void) {
  return mm_get_mem_total();
}
int64_t mm_get_mem_available_cached(void) {
  return mm_get_mem_total();
}
int64_t mm_get_mem_available_free(void) {
  return mm_get_mem_total();
}
int64_t mm_get_mem_available_total(void) {
  return mm_get_mem_total();
}
int64_t mm_get_mem_total(void) {
  int64_t x;
  size_t l = sizeof(x);
  sysctlbyname("hw.memsize", &x, &l, NULL, 0);
  return x;
}
#else
int64_t mm_get_page_size(void) {
  int64_t page_size = sysconf(_SC_PAGESIZE);
  gem_cond_fatal_error(page_size==-1,SYS_SYSCONF);
  return page_size; // Bytes
}
int64_t mm_get_stat_meminfo(const char* const label,const uint64_t label_length) {
  // Open /proc/meminfo
  FILE* const meminfo = fopen("/proc/meminfo", "r");
  gem_cond_fatal_error(meminfo == NULL,MEM_STAT_MEMINFO,"no such file");
  // Parse
  char *line = NULL;
  uint64_t chars_read=0, size=0;
  size_t line_length=0;
  while ((chars_read=getline(&line,&line_length,meminfo))!=-1) {
    if (strncmp(line,label,label_length)==0) {
      sscanf(line,"%*s %"PRIu64"",&size);
      free(line); fclose(meminfo);
      return size*1024; // Bytes
    }
  }
  // Not found
  free(line);
  fclose(meminfo);
  return -1;
}
int64_t mm_get_mem_available_virtual(void) {
  // Get Total Program Size
  uint64_t vm_size = 0;
  FILE *statm = fopen("/proc/self/statm", "r");
  gem_cond_fatal_error(fscanf(statm,"%"PRIu64,&vm_size)==-1,MEM_PARSE_STATM);
  vm_size = (vm_size + 1) * 1024;
  fclose(statm);
  // Get Virtual Memory Limit
  struct rlimit lim;
  getrlimit(RLIMIT_AS,&lim);
  // Return Virtual Memory Available
  return lim.rlim_cur - vm_size; // Bytes
}
int64_t mm_get_mem_available_cached(void) {
  const int64_t size = mm_get_stat_meminfo("Cached:",7);
  gem_cond_fatal_error(size==-1,MEM_STAT_MEMINFO,"Cached");
  return size; // Bytes
}
int64_t mm_get_mem_available_free(void) {
  const int64_t size = mm_get_stat_meminfo("MemFree:",8);
  gem_cond_fatal_error(size==-1,MEM_STAT_MEMINFO,"MemFree");
  return size; // Bytes
}
int64_t mm_get_mem_available_total(void) {
  return mm_get_mem_available_free()+mm_get_mem_available_cached(); // Bytes
}
int64_t mm_get_mem_total(void) {
  const int64_t size = mm_get_stat_meminfo("MemTotal:",9);
  gem_cond_fatal_error(size==-1,MEM_STAT_MEMINFO,"MemTotal");
  return size; // Bytes
}
#endif

