/*
 * PROJECT: GEMMapper
 * FILE: mm.h
 * DATE: 06/06/2012
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

#ifndef MEMORY_MANAGEMENT_H_
#define MEMORY_MANAGEMENT_H_

#include "commons.h"
#include "errors.h"

/******************************************************************************
 *
 * Memory Alignment Utils
 *
 ******************************************************************************/
// Check Memory Alignment
#define MM_MEM_ALIGNED_MASK_16b  mm_mem_alignment_bits_mask[0]
#define MM_MEM_ALIGNED_MASK_32b  mm_mem_alignment_bits_mask[1]
#define MM_MEM_ALIGNED_MASK_64b  mm_mem_alignment_bits_mask[2]
#define MM_MEM_ALIGNED_MASK_128b mm_mem_alignment_bits_mask[3]
#define MM_MEM_ALIGNED_MASK_256b mm_mem_alignment_bits_mask[4]
#define MM_MEM_ALIGNED_MASK_512b mm_mem_alignment_bits_mask[5]
#define MM_MEM_ALIGNED_MASK_1KB  mm_mem_alignment_bits_mask[6]
#define MM_MEM_ALIGNED_MASK_2KB  mm_mem_alignment_bits_mask[7]
#define MM_MEM_ALIGNED_MASK_4KB  mm_mem_alignment_bits_mask[8]
#define MM_MEM_ALIGNED_MASK_8KB  mm_mem_alignment_bits_mask[9]
#define MM_MEM_ALIGNED_MASK_16KB mm_mem_alignment_bits_mask[10]
#define MM_MEM_ALIGNED_MASK_32KB mm_mem_alignment_bits_mask[11]
#define MM_MEM_ALIGNED_MASK_2MB  mm_mem_alignment_bits_mask[12]
#define MM_MEM_ALIGNED_MASK_4MB  mm_mem_alignment_bits_mask[13]
// Check Memory Alignment Bits (Masks)
extern const uint64_t mm_mem_alignment_bits_mask[];
#define MM_MEM_IS_ALIGNED(memory_address,align) ((MM_CAST_ADDR(memory_address) & MM_MEM_ALIGNED_MASK_##align)==0)
// Cast to address as to get the alignment
#define MM_CAST_ADDR(memory_address) ((uintptr_t)(const void *)(memory_address))
#define MM_CAST_PTR(memory_address)  ((void *)(memory_address))
// Debug print the alignment of the memory block
#define MM_PRINT_MEM_ALIGMENT(memory_address) \
  if (MM_MEM_IS_ALIGNED(memory_address,4MB)) { \
    debug_msg("Memory aligned to 4MB (%p)",(void*)memory_address); \
  } else if (MM_MEM_IS_ALIGNED(memory_address,2MB)) { \
    debug_msg("Memory aligned to 2MB (%p)",(void*)memory_address); \
  } else if (MM_MEM_IS_ALIGNED(memory_address,32KB)) { \
    debug_msg("Memory aligned to 32KB (%p)",(void*)memory_address); \
  } else if (MM_MEM_IS_ALIGNED(memory_address,16KB)) { \
    debug_msg("Memory aligned to 16KB (%p)",(void*)memory_address); \
  } else if (MM_MEM_IS_ALIGNED(memory_address,8KB)) { \
    debug_msg("Memory aligned to 8KB (%p)",(void*)memory_address); \
  } else if (MM_MEM_IS_ALIGNED(memory_address,4KB)) { \
    debug_msg("Memory aligned to 4KB (%p)",(void*)memory_address); \
  } else if (MM_MEM_IS_ALIGNED(memory_address,2KB)) { \
    debug_msg("Memory aligned to 2KB (%p)",(void*)memory_address); \
  } else if (MM_MEM_IS_ALIGNED(memory_address,1KB)) { \
    debug_msg("Memory aligned to 1KB (%p)",(void*)memory_address); \
  } else if (MM_MEM_IS_ALIGNED(memory_address,512b)) { \
    debug_msg("Memory aligned to 512 bits (%p)",(void*)memory_address); \
  } else if (MM_MEM_IS_ALIGNED(memory_address,256b)) { \
    debug_msg("Memory aligned to 256 bits (%p)",(void*)memory_address); \
  } else if (MM_MEM_IS_ALIGNED(memory_address,128b)) { \
    debug_msg("Memory aligned to 128 bits (%p)",(void*)memory_address); \
  } else if (MM_MEM_IS_ALIGNED(memory_address,64b)) { \
    debug_msg("Memory aligned to 64 bits (%p)",(void*)memory_address); \
  } else if (MM_MEM_IS_ALIGNED(memory_address,32b)) { \
    debug_msg("Memory aligned to 32 bits (%p)",(void*)memory_address); \
  } else if (MM_MEM_IS_ALIGNED(memory_address,16b)) { \
    debug_msg("Memory aligned to 16 bits (%p)",(void*)memory_address); \
  }

/*
 * Temporal folder path
 */
#define MM_DEFAULT_TMP_FOLDER "./"

GEM_INLINE char* mm_get_tmp_folder();
GEM_INLINE void mm_set_tmp_folder(char* const tmp_folder_path);

/******************************************************************************
 * UnitMemory
 *   Allocate relative small chunks of memory relying on the regular memory manager,
 *   usually malloc/calloc using a BuddySystem (Helper functions)
 ******************************************************************************/
#define mm_alloc(type) ((type*)mm_malloc_(1,sizeof(type),false,0))
#define mm_malloc(num_bytes) (mm_malloc_(1,num_bytes,false,0))
#define mm_calloc(num_elements,type,clear_mem) ((type*)mm_malloc_(num_elements,sizeof(type),clear_mem,0))
#define mm_malloc_uint64() mm_malloc(sizeof(uint64_t))
#define mm_malloc_uint32() mm_malloc(sizeof(uint32_t))
#define mm_malloc_uint16() mm_malloc(sizeof(uint16_t))
#define mm_malloc_uint8()  mm_malloc(sizeof(uint8_t))
GEM_INLINE void* mm_malloc_(uint64_t const num_elements,const uint64_t size_element,const bool init_mem,const int init_value);
GEM_INLINE void* mm_malloc_nothrow(uint64_t const num_elements,const uint64_t size_element,const bool init_mem,const int init_value);
GEM_INLINE void* mm_realloc(void* mem_addr,const uint64_t num_bytes);
GEM_INLINE void* mm_realloc_nothrow(void* mem_addr,const uint64_t num_bytes);
GEM_INLINE void mm_free(void* mem_addr);

/******************************************************************************
 * BulkMemory
 *   Allocate big chunks of memory and resort to disk if memory is not enough
 *
 ******************************************************************************/
typedef enum { MM_READ_ONLY=0, MM_WRITE_ONLY=1, MM_READ_WRITE=2 } mm_mode;
typedef enum { MM_HEAP, MM_MMAPPED } mm_type;
typedef struct {
  mm_type mem_type;       /* Type {Heap,MMapped} */
  mm_mode mode;           /* Memory mode {R,W,R/W} */
  /* Memory Data */
  uint64_t allocated;     /* Total allocated Bytes*/
  void* memory;           /* Pointer to block of memory */
  void* cursor;           /* Pointer current position of memory */
  /* Mapped File Memory*/
  int fd;                 /* File descriptor */
  char *file_name;        /* File name */
} mm_t;

// Checkers
#define MM_CHECK(mm) \
  GEM_CHECK_NULL((mm)); \
  GEM_CHECK_NULL((mm)->memory); \
  GEM_CHECK_NULL((mm)->cursor); \
  MM_CHECK_SEGMENT(mm)

#define MM_CHECK_SEGMENT(mm) \
  gem_fatal_check( ((mm)->cursor <  (mm)->memory) || \
                   ((mm)->cursor >= (mm)->memory+(mm)->allocated),MEM_CURSOR_OUT_OF_SEGMENT)

#define MM_CHECK_ALIGNMENT(mm,align) \
  gem_fatal_check(!MM_MEM_IS_ALIGNED((mm)->cursor,align),MEM_ALG_FAILED);

// Array allocation
#define mm_malloc_array(mem_addr,num_vectors,vector_size,type) { \
  void* memory_allocated = mm_malloc(num_vectors*(sizeof(type*)+vector_size)); \
  mem_addr = memory_allocated; \
  memory_allocated += num_vectors*sizeof(type*); \
  uint64_t i; \
  for (i=0;i<num_vectors;++i) { \
    mem_addr[i] = memory_allocated; \
    memory_allocated += vector_size; \
  } \
}

// Memory/MFile Allocators/Handlers
GEM_INLINE mm_t* mm_bulk_malloc(const uint64_t num_bytes,const bool init_mem);
GEM_INLINE mm_t* mm_bulk_mmalloc(const uint64_t num_bytes,const bool use_huge_pages);
GEM_INLINE mm_t* mm_bulk_mmalloc_temp(const uint64_t num_bytes);
GEM_INLINE void mm_bulk_free(mm_t* const mem_manager);

GEM_INLINE mm_t* mm_bulk_mmap_file(char* const file_name,const mm_mode mode,const bool populate_page_tables);
GEM_INLINE mm_t* mm_bulk_load_file(char* const file_name,const uint64_t num_threads);
GEM_INLINE mm_t* mm_bulk_mload_file(char* const file_name,const uint64_t num_threads);

// Accessors
GEM_INLINE void* mm_get_mem(mm_t* const mem_manager);
GEM_INLINE void* mm_get_base_mem(mm_t* const mem_manager);
GEM_INLINE mm_mode mm_get_mode(mm_t* const mem_manager);
GEM_INLINE uint64_t mm_get_allocated(mm_t* const mem_manager);
GEM_INLINE char* mm_get_mfile_name(mm_t* const mem_manager);

// Seek
GEM_INLINE uint64_t mm_get_current_position(mm_t* const mem_manager);
GEM_INLINE bool mm_eom(mm_t* const mem_manager);
GEM_INLINE void mm_seek(mm_t* const mem_manager,const uint64_t byte_position);

GEM_INLINE void mm_skip_forward(mm_t* const mem_manager,const uint64_t num_bytes);
GEM_INLINE void mm_skip_backward(mm_t* const mem_manager,const uint64_t num_bytes);
GEM_INLINE void mm_skip_uint64(mm_t* const mem_manager);
GEM_INLINE void mm_skip_uint32(mm_t* const mem_manager);
GEM_INLINE void mm_skip_uint16(mm_t* const mem_manager);
GEM_INLINE void mm_skip_uint8(mm_t* const mem_manager);
GEM_INLINE void mm_skip_align(mm_t* const mem_manager,const uint64_t num_bytes);
GEM_INLINE void mm_skip_align_16(mm_t* const mem_manager);
GEM_INLINE void mm_skip_align_32(mm_t* const mem_manager);
GEM_INLINE void mm_skip_align_64(mm_t* const mem_manager);
GEM_INLINE void mm_skip_align_128(mm_t* const mem_manager);
GEM_INLINE void mm_skip_align_512(mm_t* const mem_manager);
GEM_INLINE void mm_skip_align_1024(mm_t* const mem_manager);
GEM_INLINE void mm_skip_align_4KB(mm_t* const mem_manager);
GEM_INLINE void mm_skip_align_mempage(mm_t* const mem_manager);

// Read
#define mm_read(mem_manager,var) mm_copy_mem(mem_manager,&var,sizeof(var))
GEM_INLINE uint64_t mm_read_uint64(mm_t* const mem_manager);
GEM_INLINE uint32_t mm_read_uint32(mm_t* const mem_manager);
GEM_INLINE uint16_t mm_read_uint16(mm_t* const mem_manager);
GEM_INLINE uint8_t mm_read_uint8(mm_t* const mem_manager);
GEM_INLINE void* mm_read_mem(mm_t* const mem_manager,const uint64_t num_bytes);
GEM_INLINE void mm_copy_mem(mm_t* const mem_manager,void* const dst,const uint64_t num_bytes);
GEM_INLINE void mm_copy_mem_parallel(mm_t* const mem_manager,void* const dst,const uint64_t num_bytes,const uint64_t num_threads);

// Write
#define mm_write(mem_manager,var) mm_write_mem(mem_manager,&var,sizeof(var))
GEM_INLINE void mm_write_uint64(mm_t* const mem_manager,const uint64_t data);
GEM_INLINE void mm_write_uint32(mm_t* const mem_manager,const uint32_t data);
GEM_INLINE void mm_write_uint16(mm_t* const mem_manager,const uint16_t data);
GEM_INLINE void mm_write_uint8(mm_t* const mem_manager,const uint8_t data);
GEM_INLINE void mm_write_mem(mm_t* const mem_manager,void* const src,const uint64_t num_bytes);

// Status
GEM_INLINE int64_t mm_get_page_size();
GEM_INLINE int64_t mm_get_available_cached_mem();
GEM_INLINE int64_t mm_get_available_free_mem();
GEM_INLINE int64_t mm_get_available_mem();
GEM_INLINE int64_t mm_get_available_virtual_mem();

/*
 * Error Messages
 */
//#define GEM_ERROR_MEM_HANDLER "Could not allocate handler"
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
#define GEM_ERROR_MEM_PARSE_STATM "Could not parse stat-memory information"

#endif /* MEMORY_MANAGEMENT_H_ */
