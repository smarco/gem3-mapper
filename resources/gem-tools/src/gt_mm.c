/*
 * PROJECT: GEM-Tools library
 * FILE: gt_mm.c
 * DATE: 01/02/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Memory Manager provides memory allocation functions. Different types of memory are supported.
 *     - UnitMemory
 *         Allocate relative small chunks of memory relying on the regular memory manager,
 *         usually malloc/calloc using a BuddySystem (Helper functions)
 *     - BulkMemory
 *         Allocate big chunks of memory and resort to disk if memory is not enough
 *     - SlabMemory
 *         Relative big amounts of objects allocated all at once (like the LINUX slab allocator)
 *         Objects of a certain type are ready to go inside the slab, thus reducing
 *         the overhead of malloc/setup/free cycles along the program
 *     - PoolMemory
 *         Pool of Slabs as gather all slabs needed along a program
 *         The goal is to minimize all memory malloc/setup/free overhead
 *         Offers thread safe allocation of slabs as to balance memory consumption across threads
 */

// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
// 1.- TCMalloc : Thread-Caching Malloc
// 2.- nedmalloc()
// 4.- madvise() / readahead() / posix_fadvise()
// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

#include "gt_mm.h"

// In some environments MAP_HUGETLB can be undefined
#ifndef MAP_HUGETLB
  #define MAP_HUGETLB 0
#endif
#ifndef MAP_ANONYMOUS
  #define MAP_ANONYMOUS 0
#endif
#ifndef MAP_POPULATE
  #define MAP_POPULATE 0
#endif

/*
 * Memory Alignment Utils
 */
const uint64_t gt_mm_mem_alignment_bits_mask[] = { // Check Memory Alignment Bits (Masks)
    0x0000000000000001lu, /*    16 bits aligned (  2B / 2^4)  */
    0x0000000000000003lu, /*    32 bits aligned (  4B / 2^5)  */
    0x0000000000000007lu, /*    64 bits aligned (  8B / 2^6)  */
    0x000000000000000Flu, /*   128 bits aligned ( 16B / 2^7)  */
    0x000000000000001Flu, /*   256 bits aligned ( 32B / 2^8)  */
    0x000000000000003Flu, /*   512 bits aligned ( 64B / 2^9)  */
    0x000000000000007Flu, /*  1024 bits aligned ( 1KB / 2^10) */
    0x00000000000000FFlu, /*  2048 bits aligned ( 2KB / 2^11) */
    0x00000000000001FFlu, /*  4096 bits aligned ( 4KB / 2^12) RegularPage Size*/
    0x00000000000003FFlu, /*  8192 bits aligned ( 8KB / 2^13) */
    0x00000000000007FFlu, /* 16384 bits aligned (16KB / 2^14) */
    0x0000000000000FFFlu, /* 32768 bits aligned (32KB / 2^15) */
    0x000000000003FFFFlu, /*   n/a bits aligned ( 2MB / 2^21) RegularPageHugeTLB Size */
    0x000000000007FFFFlu, /*   n/a bits aligned ( 4MB / 2^21) */
};

/*
 * MMap Constants/Values
 */
int gt_mm_proc_flags[3] = { PROT_READ, PROT_READ|PROT_WRITE, PROT_READ|PROT_WRITE };
int gt_mm_mmap_mode[3] = { MAP_PRIVATE, MAP_SHARED, MAP_SHARED };
/*
 * Temporal folder path
 */
char* gt_mm_temp_folder_path = GT_MM_DEFAULT_TMP_FOLDER;

GT_INLINE char* gt_mm_get_tmp_folder() {
  return gt_mm_temp_folder_path;
}
GT_INLINE void gt_mm_set_tmp_folder(char* const tmp_folder_path) {
  GT_NULL_CHECK(tmp_folder_path);
  gt_mm_temp_folder_path = tmp_folder_path;
}
/*
 * UnitMemory
 *   Allocate relative small chunks of memory relying on the regular memory manager,
 *   usually malloc/calloc using a BuddySystem (Helper functions)
 */
GT_INLINE void* gt_malloc_nothrow(uint64_t const num_elements,const uint64_t size_element,const bool init_mem,const int init_value) {
  // Request memory
  const uint64_t total_memory = num_elements*size_element;
  void* const allocated_mem = (gt_expect_false(init_mem && init_value==0)) ?
      calloc(num_elements,size_element) : malloc(total_memory);
  if (!allocated_mem) return NULL;
  // Initialize memory
  if (gt_expect_false(init_mem && init_value!=0)) memset(allocated_mem,init_value,total_memory);
  //GT_MM_PRINT_MEM_ALIGMENT(allocated_mem); // Debug
  return allocated_mem;
}
GT_INLINE void* gt_malloc_(uint64_t const num_elements,const uint64_t size_element,const bool init_mem,const int init_value) {
  // Request memory
  const uint64_t total_memory = num_elements*size_element;
  void* allocated_mem;
  if (gt_expect_false(init_mem && init_value==0)) {
    allocated_mem = calloc(num_elements,size_element);
    gt_cond_fatal_error(!allocated_mem,MEM_CALLOC,num_elements,size_element);
  } else {
    allocated_mem = malloc(total_memory);
    gt_cond_fatal_error(!allocated_mem,MEM_ALLOC,total_memory);
  }
  if (gt_expect_false(init_mem && init_value!=0)) memset(allocated_mem,init_value,total_memory);
  //GT_MM_PRINT_MEM_ALIGMENT(allocated_mem); // Debug
  return allocated_mem;
}
GT_INLINE void* gt_realloc_nothrow(void* mem_addr,const uint64_t num_bytes) {
  GT_NULL_CHECK(mem_addr);
  GT_ZERO_CHECK(num_bytes);
  return realloc(mem_addr,num_bytes);
}
GT_INLINE void* gt_realloc(void* mem_addr,const uint64_t num_bytes) {
  GT_NULL_CHECK(mem_addr);
  GT_ZERO_CHECK(num_bytes);
  void* const new_mem_addr = gt_realloc_nothrow(mem_addr,num_bytes);
  gt_cond_fatal_error(!new_mem_addr,MEM_REALLOC,num_bytes);
  return new_mem_addr;
}
GT_INLINE void gt_free(void* mem_addr) {
  free(mem_addr);
}

/*
 * BulkMemory
 *   Allocate big chunks of memory and resort to disk if memory is not enough
 */
GT_INLINE gt_mm* gt_mm_bulk_malloc(const uint64_t num_bytes,const bool init_mem) {
  GT_ZERO_CHECK(num_bytes);
  void* memory = gt_malloc_nothrow(num_bytes,1,init_mem,0);
  if (gt_expect_true(memory!=NULL)) { // Fits in HEAP
    gt_mm* const mm = gt_alloc(gt_mm);
    mm->memory = memory;
    mm->mem_type = GT_MM_HEAP;
    mm->mode = GT_MM_READ_WRITE;
    mm->allocated = num_bytes;
    mm->cursor = mm->memory;
    GT_MM_PRINT_MEM_ALIGMENT(mm->memory); // Debug
    return mm;
  } else { // Resort to MMAP in disk
    gt_warn(MEM_ALLOC_DISK,num_bytes);
    return gt_mm_bulk_mmalloc_temp(num_bytes);
  }
}
GT_INLINE gt_mm* gt_mm_bulk_mmalloc(const uint64_t num_bytes,const bool use_huge_pages) {
  GT_ZERO_CHECK(num_bytes);
  // Allocate handler
  gt_mm* const mm = gt_alloc(gt_mm);
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
  int flags = MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE;
  if (use_huge_pages) flags |= MAP_HUGETLB;
  mm->memory = mmap(0,num_bytes,PROT_READ|PROT_WRITE,flags,-1,0);
  gt_cond_fatal_error__perror(mm->memory==MAP_FAILED,MEM_ALLOC_MMAP_FAIL,num_bytes);
  mm->cursor = mm->memory;
  // Set MM
  mm->mem_type = GT_MM_MMAPPED;
  mm->mode = GT_MM_READ_WRITE;
  mm->allocated = num_bytes;
  mm->fd = -1;
  mm->file_name = NULL;
  // GT_MM_PRINT_MEM_ALIGMENT(mm->memory); // Debug
  return mm;
}
GT_INLINE gt_mm* gt_mm_bulk_mmalloc_temp(const uint64_t num_bytes) {
  GT_ZERO_CHECK(num_bytes);
  // Allocate handler
  gt_mm* const mm = gt_alloc(gt_mm);
  // TemporalMemory (backed by a file)
  mm->file_name = gt_calloc(strlen(gt_mm_get_tmp_folder())+22,char,true);
  sprintf(mm->file_name,"%sgt_mmalloc_temp_XXXXXX",gt_mm_get_tmp_folder());
  // Create temporary file
  mm->fd = mkstemp(mm->file_name);
  gt_cond_fatal_error__perror(mm->fd==-1,SYS_MKSTEMP,mm->file_name);
  gt_cond_fatal_error__perror(unlink(mm->file_name),SYS_HANDLE_TMP); // Make it temporary
  // Set the size of the temporary file (disk allocation)
  gt_cond_fatal_error__perror(lseek(mm->fd,num_bytes-1,SEEK_SET)==-1,SYS_HANDLE_TMP);
  gt_cond_fatal_error__perror(write(mm->fd,"",1)<=0,SYS_HANDLE_TMP);
  gt_cond_fatal_error__perror(lseek(mm->fd,0,SEEK_SET)==-1,SYS_HANDLE_TMP);
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
  mm->memory = mmap(NULL,num_bytes,PROT_READ|PROT_WRITE,MAP_SHARED|MAP_NORESERVE,mm->fd,0);
  gt_cond_fatal_error__perror(mm->memory==MAP_FAILED,MEM_ALLOC_MMAP_DISK_FAIL,num_bytes,mm->file_name);
  mm->cursor = mm->memory;
  // Set MM
  mm->mem_type = GT_MM_MMAPPED;
  mm->mode = GT_MM_READ_WRITE;
  mm->allocated = num_bytes;
  // GT_MM_PRINT_MEM_ALIGMENT(mm->memory); // Debug
  return mm;
}
GT_INLINE void gt_mm_realloc(gt_mm* const mm,const uint64_t num_bytes) {
  GT_MM_CHECK(mm);
  gt_fatal_error(NOT_IMPLEMENTED); // TODO
//  const uint64_t current_cursor_pos = gt_mm_get_current_position(mm);
//  if (mm->mem_type==GT_MM_HEAP) { // Heap BulkMemory
//    if (num_bytes > mm->allocated) {
//      mm->memory = realloc(mm->memory);
//      gt_cond_fatal_error(mm->memory==NULL,MEM_REALLOC);
//      mm->cursor = mm->memory + current_cursor_pos;
//      mm->allocated = num_bytes;
//    }
//  } else { // MMapped BulkMemory
//    if (mm->fd!=-1) {
//      if (mm->tmp_file) { // TemporalMemory
//        mremap(mm->memory,mm->allocated,num_bytes,MREMAP_MAYMOVE);
//      } else { // File mapped
//
//      }
//    } else { // Anonymous
//
//    }
//  }
}
GT_INLINE void gt_mm_free(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  if (mm->mem_type==GT_MM_HEAP) { // Heap BulkMemory
    gt_free(mm->memory);
  } else { // MMapped BulkMemory
    gt_cond_fatal_error__perror(munmap(mm->memory,mm->allocated)==-1,SYS_UNMAP);
    if (mm->fd!=-1) {
      gt_cond_fatal_error__perror(close(mm->fd),SYS_HANDLE_TMP);
    }
  }
  gt_cfree(mm->file_name);
  gt_free(mm);
}
GT_INLINE gt_mm* gt_mm_bulk_mmap_file(char* const file_name,const gt_mm_mode mode,const bool populate_page_tables) {
  GT_NULL_CHECK(file_name);
  // Allocate handler
  gt_mm* const mm = gt_alloc(gt_mm);
  // Retrieve input file info
  struct stat stat_info;
  gt_stat(file_name,&stat_info);
  // Open file descriptor
  mm->fd = gt_open_fd(file_name,gt_fm_open_flags[mode],S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
  /*
   * Mmap file
   *   - @mode::
   *       GT_MM_READ_ONLY => MAP_PRIVATE (no copy-on-write as it's not allowed)
   *       GT_MM_WRITE_ONLY or GT_MM_READ_WRITE => MAP_SHARED
   *   - MAP_POPULATE (since Linux 2.5.46)
   *       Populate (prefault) page tables for a mapping. For a file mapping, this causes
   *       read-ahead on the file. Later accesses to the mapping will not be blocked by page faults.
   *       MAP_POPULATE is only supported for private mappings since Linux 2.6.23.
   */
  int flags = gt_mm_mmap_mode[mode];
  if (populate_page_tables) flags |= MAP_POPULATE;
  mm->memory = mmap(0,stat_info.st_size,gt_mm_proc_flags[mode],flags,mm->fd,0);
  gt_cond_fatal_error__perror(mm->memory==MAP_FAILED,SYS_MMAP_FILE,file_name);
  mm->cursor = mm->memory;
  // Set MM
  mm->mem_type = GT_MM_MMAPPED;
  mm->mode = mode;
  mm->allocated = stat_info.st_size;
  mm->file_name = gt_strndup(file_name,gt_strlen(file_name));
  // GT_MM_PRINT_MEM_ALIGMENT(mm->memory); // Debug
  return mm;
}
GT_INLINE gt_mm* gt_mm_bulk_load_file(char* const file_name,const uint64_t num_threads) {
  GT_NULL_CHECK(file_name);
  // Allocate handler
  gt_mm* const mm = gt_alloc(gt_mm);
  // Retrieve input file info
  struct stat stat_info;
  gt_stat(file_name,&stat_info);
  // Allocate memory to dump the content of the file
  mm->memory = gt_malloc(stat_info.st_size);
  gt_cond_fatal_error(!mm->memory,MEM_ALLOC,stat_info.st_size);
  mm->mem_type = GT_MM_HEAP;
  mm->mode = GT_MM_READ_ONLY;
  mm->allocated = stat_info.st_size;
  mm->cursor = mm->memory;
  GT_MM_PRINT_MEM_ALIGMENT(mm->memory); // Debug
  // Read the file and dump it into memory
  if (num_threads>1 && (stat_info.st_size > num_threads*8)) {
    gt_fm_bulk_read_file_parallel(file_name,mm->memory,0,0,num_threads);
  } else {
    gt_fm_bulk_read_file(file_name,mm->memory,0,0);
  }
  return mm;
}
GT_INLINE gt_mm* gt_mm_bulk_mload_file(char* const file_name,const uint64_t num_threads) {
  GT_NULL_CHECK(file_name);
  // Retrieve input file info
  struct stat stat_info;
  gt_stat(file_name,&stat_info);
  // Allocate memory to dump the content of the file
  gt_mm* const mm = gt_mm_bulk_mmalloc(stat_info.st_size,false);
  // Read the file and dump it into memory
  if (num_threads>1 && (stat_info.st_size > num_threads*8)) {
    gt_fm_bulk_read_file_parallel(file_name,mm->memory,0,0,num_threads);
  } else {
    gt_fm_bulk_read_file(file_name,mm->memory,0,0);
  }
  return mm;
}

/*
 * Accessors
 */
GT_INLINE void* gt_mm_get_mem(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  return mm->cursor;
}
GT_INLINE void* gt_mm_get_base_mem(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  return mm->memory;
}
GT_INLINE gt_mm_mode gt_mm_get_mode(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  return mm->mode;
}
GT_INLINE uint64_t gt_mm_get_allocated(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  return mm->allocated;
}
GT_INLINE char* gt_mm_get_mfile_name(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  return mm->file_name;
}
/*
 * Seek functions
 */
GT_INLINE uint64_t gt_mm_get_current_position(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  return (mm->cursor-mm->memory);
}
GT_INLINE bool gt_mm_eom(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  return gt_mm_get_current_position(mm) >= mm->allocated;
}
GT_INLINE void gt_mm_seek(gt_mm* const mm,const uint64_t byte_position) {
  GT_MM_CHECK(mm);
  gt_fatal_check(byte_position>=mm->allocated,MEM_CURSOR_SEEK,byte_position);
  mm->cursor = mm->memory + byte_position;
}

GT_INLINE void gt_mm_skip_forward(gt_mm* const mm,const uint64_t num_bytes) {
  GT_MM_CHECK(mm);
  mm->cursor += num_bytes;
  GT_MM_CHECK_SEGMENT(mm);
}
GT_INLINE void gt_mm_skip_backward(gt_mm* const mm,const uint64_t num_bytes) {
  GT_MM_CHECK(mm);
  mm->cursor -= num_bytes;
  GT_MM_CHECK_SEGMENT(mm);
}
GT_INLINE void gt_mm_skip_uint64(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  mm->cursor += 8;
  GT_MM_CHECK_SEGMENT(mm);
}
GT_INLINE void gt_mm_skip_uint32(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  mm->cursor += 4;
  GT_MM_CHECK_SEGMENT(mm);
}
GT_INLINE void gt_mm_skip_uint16(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  mm->cursor += 2;
  GT_MM_CHECK_SEGMENT(mm);
}
GT_INLINE void gt_mm_skip_uint8(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  mm->cursor += 1;
  GT_MM_CHECK_SEGMENT(mm);
}
GT_INLINE void gt_mm_skip_align(gt_mm* const mm,const uint64_t num_bytes) {
  GT_MM_CHECK(mm);
  GT_ZERO_CHECK(num_bytes);
  if (gt_expect_true(num_bytes > 1)) {
    mm->cursor = mm->cursor+(num_bytes-1);
    mm->cursor = mm->cursor-(GT_MM_CAST_ADDR(mm->cursor)%num_bytes);
    GT_MM_CHECK_SEGMENT(mm);
    gt_fatal_check(GT_MM_CAST_ADDR(mm->cursor)%num_bytes!=0,MEM_ALG_FAILED);
  }
}
GT_INLINE void gt_mm_skip_align_16(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  mm->cursor = GT_MM_CAST_PTR(
      (GT_MM_CAST_ADDR(mm->cursor)+GT_MM_MEM_ALIGNED_MASK_16b) & (~GT_MM_MEM_ALIGNED_MASK_16b));
  GT_MM_CHECK_SEGMENT(mm);
  GT_MM_CHECK_ALIGNMENT(mm,16b);
}
GT_INLINE void gt_mm_skip_align_32(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  mm->cursor = GT_MM_CAST_PTR(
      (GT_MM_CAST_ADDR(mm->cursor)+GT_MM_MEM_ALIGNED_MASK_32b) & (~GT_MM_MEM_ALIGNED_MASK_32b));
  GT_MM_CHECK_SEGMENT(mm);
  GT_MM_CHECK_ALIGNMENT(mm,32b);
}
GT_INLINE void gt_mm_skip_align_64(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  mm->cursor = GT_MM_CAST_PTR(
      (GT_MM_CAST_ADDR(mm->cursor)+GT_MM_MEM_ALIGNED_MASK_64b) & (~GT_MM_MEM_ALIGNED_MASK_64b));
  GT_MM_CHECK_SEGMENT(mm);
  GT_MM_CHECK_ALIGNMENT(mm,64b);
}
GT_INLINE void gt_mm_skip_align_128(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  mm->cursor = GT_MM_CAST_PTR(
      (GT_MM_CAST_ADDR(mm->cursor)+GT_MM_MEM_ALIGNED_MASK_128b) & (~GT_MM_MEM_ALIGNED_MASK_128b));
  GT_MM_CHECK_SEGMENT(mm);
  GT_MM_CHECK_ALIGNMENT(mm,128b);
}
GT_INLINE void gt_mm_skip_align_512(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  mm->cursor = GT_MM_CAST_PTR(
      (GT_MM_CAST_ADDR(mm->cursor)+GT_MM_MEM_ALIGNED_MASK_512b) & (~GT_MM_MEM_ALIGNED_MASK_512b));
  GT_MM_CHECK_SEGMENT(mm);
  GT_MM_CHECK_ALIGNMENT(mm,512b);
}
GT_INLINE void gt_mm_skip_align_1024(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  mm->cursor = GT_MM_CAST_PTR(
      (GT_MM_CAST_ADDR(mm->cursor)+GT_MM_MEM_ALIGNED_MASK_1KB) & (~GT_MM_MEM_ALIGNED_MASK_1KB));
  GT_MM_CHECK_SEGMENT(mm);
  GT_MM_CHECK_ALIGNMENT(mm,1KB);
}
GT_INLINE void gt_mm_skip_align_4KB(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  mm->cursor = GT_MM_CAST_PTR(
      (GT_MM_CAST_ADDR(mm->cursor)+GT_MM_MEM_ALIGNED_MASK_4KB) & (~GT_MM_MEM_ALIGNED_MASK_4KB));
  GT_MM_CHECK_SEGMENT(mm);
  GT_MM_CHECK_ALIGNMENT(mm,4KB);
}
GT_INLINE void gt_mm_skip_align_mempage(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  uint64_t sz = sysconf(_SC_PAGESIZE);
  gt_mm_skip_align(mm,sz);
}

/*
 * Read functions
 */
GT_INLINE uint64_t gt_mm_read_uint64(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  GT_MM_CHECK_SEGMENT(mm);
  const uint64_t data = *((uint64_t*)mm->cursor);
  mm->cursor += 8;
  return data;
}
GT_INLINE uint32_t gt_mm_read_uint32(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  GT_MM_CHECK_SEGMENT(mm);
  const uint32_t data = *((uint32_t*)mm->cursor);
  mm->cursor += 4;
  return data;
}
GT_INLINE uint16_t gt_mm_read_uint16(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  GT_MM_CHECK_SEGMENT(mm);
  const uint16_t data = *((uint16_t*)mm->cursor);
  mm->cursor += 2;
  return data;
}
GT_INLINE uint8_t gt_mm_read_uint8(gt_mm* const mm) {
  GT_MM_CHECK(mm);
  GT_MM_CHECK_SEGMENT(mm);
  const uint8_t data = *((uint8_t*)mm->cursor);
  mm->cursor += 1;
  return data;
}
GT_INLINE void* gt_mm_read_mem(gt_mm* const mm,const uint64_t num_bytes) {
  GT_MM_CHECK(mm);
  GT_MM_CHECK_SEGMENT(mm);
  void* const current_cursor = mm->cursor;
  mm->cursor += num_bytes;
  return current_cursor;
}
GT_INLINE void gt_mm_copy_mem(gt_mm* const mm,void* const dst,const uint64_t num_bytes) {
  GT_MM_CHECK(mm);
  GT_MM_CHECK_SEGMENT(mm);
  memcpy(dst,mm->cursor,num_bytes);
  mm->cursor += num_bytes;
}
GT_INLINE void gt_mm_copy_mem_parallel(
    gt_mm* const mm,void* const dst,const uint64_t num_bytes,const uint64_t num_threads) {
  GT_MM_CHECK(mm);
  GT_MM_CHECK_SEGMENT(mm);
  // Calculate size of each chunk
  const uint64_t chunk_size = num_bytes/num_threads; // num_bytes > num_threads
#ifdef HAVE_OPENMP
  //#pragma omp parallel num_threads(num_threads)
#endif
  {
    // Calculate offsets
#ifdef HAVE_OPENMP
    const uint64_t tid = omp_get_thread_num();
#else
    const uint64_t tid = 0;
#endif
    const uint64_t offset = tid*chunk_size;
    const uint64_t size = (tid < (num_threads-1)) ? chunk_size : num_bytes-chunk_size*tid;
    // Copy the chunk
    memcpy(dst+offset,mm->cursor+offset,size);
  }
  mm->cursor += num_bytes;
}

/*
 * Write functions
 */
GT_INLINE void gt_mm_write_uint64(gt_mm* const mm,const uint64_t data) {
  GT_MM_CHECK(mm);
  GT_MM_CHECK_SEGMENT(mm);
  *((uint64_t*)mm->cursor) = data;
  mm->cursor += 8;
}
GT_INLINE void gt_mm_write_uint32(gt_mm* const mm,const uint32_t data) {
  GT_MM_CHECK(mm);
  GT_MM_CHECK_SEGMENT(mm);
  *((uint32_t*)mm->cursor) = data;
  mm->cursor += 4;
}
GT_INLINE void gt_mm_write_uint16(gt_mm* const mm,const uint16_t data) {
  GT_MM_CHECK(mm);
  GT_MM_CHECK_SEGMENT(mm);
  *((uint16_t*)mm->cursor) = data;
  mm->cursor += 2;
}
GT_INLINE void gt_mm_write_uint8(gt_mm* const mm,const uint8_t data) {
  GT_MM_CHECK(mm);
  GT_MM_CHECK_SEGMENT(mm);
  *((uint8_t*)mm->cursor) = data;
  mm->cursor += 1;
}
GT_INLINE void gt_mm_write_mem(gt_mm* const mm,void* const src,const uint64_t num_bytes) {
  GT_MM_CHECK(mm);
  GT_MM_CHECK_SEGMENT(mm);
  // TODO
}
