/*
 * PROJECT: GEM-Tools library
 * FILE: gt_mm.h
 * DATE: 01/02/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Memory Manager provides memory allocation functions.
 *   Different types of memory are supported.
 *     - UnitMemory
 *         Allocate relative small chunks of memory relying on the regular memory manager,
 *         usually malloc/calloc using a BuddySystem (Helper functions)
 *     - BulkMemory
 *         Allocate big chunks of memory and resort to disk if memory is not enough
 */

#ifndef GT_MEMORY_MANAGEMENT_H_
#define GT_MEMORY_MANAGEMENT_H_

//#include <unistd.h>
#include "gt_commons.h"
#include "gt_error.h"
#include "gt_vector.h"
#include "gt_string.h"
#include "gt_fm.h"

/*
 * Memory Alignment Utils
 */
// Check Memory Alignment
#define GT_MM_MEM_ALIGNED_MASK_16b  gt_mm_mem_alignment_bits_mask[0]
#define GT_MM_MEM_ALIGNED_MASK_32b  gt_mm_mem_alignment_bits_mask[1]
#define GT_MM_MEM_ALIGNED_MASK_64b  gt_mm_mem_alignment_bits_mask[2]
#define GT_MM_MEM_ALIGNED_MASK_128b gt_mm_mem_alignment_bits_mask[3]
#define GT_MM_MEM_ALIGNED_MASK_256b gt_mm_mem_alignment_bits_mask[4]
#define GT_MM_MEM_ALIGNED_MASK_512b gt_mm_mem_alignment_bits_mask[5]
#define GT_MM_MEM_ALIGNED_MASK_1KB  gt_mm_mem_alignment_bits_mask[6]
#define GT_MM_MEM_ALIGNED_MASK_2KB  gt_mm_mem_alignment_bits_mask[7]
#define GT_MM_MEM_ALIGNED_MASK_4KB  gt_mm_mem_alignment_bits_mask[8]
#define GT_MM_MEM_ALIGNED_MASK_8KB  gt_mm_mem_alignment_bits_mask[9]
#define GT_MM_MEM_ALIGNED_MASK_16KB gt_mm_mem_alignment_bits_mask[10]
#define GT_MM_MEM_ALIGNED_MASK_32KB gt_mm_mem_alignment_bits_mask[11]
#define GT_MM_MEM_ALIGNED_MASK_2MB  gt_mm_mem_alignment_bits_mask[12]
#define GT_MM_MEM_ALIGNED_MASK_4MB  gt_mm_mem_alignment_bits_mask[13]
// Check Memory Alignment Bits (Masks)
extern const uint64_t gt_mm_mem_alignment_bits_mask[];
#define GT_MM_MEM_IS_ALIGNED(memory_address,align) ((GT_MM_CAST_ADDR(memory_address) & GT_MM_MEM_ALIGNED_MASK_##align)==0)
// Cast to address as to get the alignment
#define GT_MM_CAST_ADDR(memory_address) ((uintptr_t)(const void *)(memory_address))
#define GT_MM_CAST_PTR(memory_address)  ((void *)(memory_address))
// Debug print the alignment of the memory block
#define GT_MM_PRINT_MEM_ALIGMENT(memory_address) \
  if (GT_MM_MEM_IS_ALIGNED(memory_address,4MB)) { \
    gt_debug_msg("Memory aligned to 4MB (%p)",(void*)memory_address); \
  } else if (GT_MM_MEM_IS_ALIGNED(memory_address,2MB)) { \
    gt_debug_msg("Memory aligned to 2MB (%p)",(void*)memory_address); \
  } else if (GT_MM_MEM_IS_ALIGNED(memory_address,32KB)) { \
    gt_debug_msg("Memory aligned to 32KB (%p)",(void*)memory_address); \
  } else if (GT_MM_MEM_IS_ALIGNED(memory_address,16KB)) { \
    gt_debug_msg("Memory aligned to 16KB (%p)",(void*)memory_address); \
  } else if (GT_MM_MEM_IS_ALIGNED(memory_address,8KB)) { \
    gt_debug_msg("Memory aligned to 8KB (%p)",(void*)memory_address); \
  } else if (GT_MM_MEM_IS_ALIGNED(memory_address,4KB)) { \
    gt_debug_msg("Memory aligned to 4KB (%p)",(void*)memory_address); \
  } else if (GT_MM_MEM_IS_ALIGNED(memory_address,2KB)) { \
    gt_debug_msg("Memory aligned to 2KB (%p)",(void*)memory_address); \
  } else if (GT_MM_MEM_IS_ALIGNED(memory_address,1KB)) { \
    gt_debug_msg("Memory aligned to 1KB (%p)",(void*)memory_address); \
  } else if (GT_MM_MEM_IS_ALIGNED(memory_address,512b)) { \
    gt_debug_msg("Memory aligned to 512 bits (%p)",(void*)memory_address); \
  } else if (GT_MM_MEM_IS_ALIGNED(memory_address,256b)) { \
    gt_debug_msg("Memory aligned to 256 bits (%p)",(void*)memory_address); \
  } else if (GT_MM_MEM_IS_ALIGNED(memory_address,128b)) { \
    gt_debug_msg("Memory aligned to 128 bits (%p)",(void*)memory_address); \
  } else if (GT_MM_MEM_IS_ALIGNED(memory_address,64b)) { \
    gt_debug_msg("Memory aligned to 64 bits (%p)",(void*)memory_address); \
  } else if (GT_MM_MEM_IS_ALIGNED(memory_address,32b)) { \
    gt_debug_msg("Memory aligned to 32 bits (%p)",(void*)memory_address); \
  } else if (GT_MM_MEM_IS_ALIGNED(memory_address,16b)) { \
    gt_debug_msg("Memory aligned to 16 bits (%p)",(void*)memory_address); \
  }

/*
 * Temporal folder path
 */
#define GT_MM_DEFAULT_TMP_FOLDER "/tmp/"

GT_INLINE char* gt_mm_get_tmp_folder();
GT_INLINE void gt_mm_set_tmp_folder(char* const tmp_folder_path);

/*
 * UnitMemory
 *   Allocate relative small chunks of memory relying on the regular memory manager,
 *   usually malloc/calloc using a BuddySystem (Helper functions)
 */
#define gt_alloc(type) ((type*)gt_malloc_(1,sizeof(type),false,0))
#define gt_malloc(num_bytes) (gt_malloc_(1,num_bytes,false,0))
#define gt_calloc(num_elements,type,clear_mem) ((type*)gt_malloc_(num_elements,sizeof(type),clear_mem,0))
#define gt_malloc_uint64() gt_malloc(sizeof(uint64_t))
#define gt_malloc_uint32() gt_malloc(sizeof(uint32_t))
#define gt_malloc_uint16() gt_malloc(sizeof(uint16_t))
#define gt_malloc_uint8()  gt_malloc(sizeof(uint8_t))
GT_INLINE void* gt_malloc_(uint64_t const num_elements,const uint64_t size_element,const bool init_mem,const int init_value);
GT_INLINE void* gt_malloc_nothrow(uint64_t const num_elements,const uint64_t size_element,const bool init_mem,const int init_value);
GT_INLINE void* gt_realloc(void* mem_addr,const uint64_t num_bytes);
GT_INLINE void* gt_realloc_nothrow(void* mem_addr,const uint64_t num_bytes);
GT_INLINE void gt_free(void* mem_addr);

/*
 * BulkMemory
 *   Allocate big chunks of memory and resort to disk if memory is not enough
 */
typedef enum { GT_MM_READ_ONLY=0, GT_MM_WRITE_ONLY=1, GT_MM_READ_WRITE=2 } gt_mm_mode;
typedef enum { GT_MM_HEAP, GT_MM_MMAPPED } gt_mm_t;
typedef struct {
  gt_mm_t mem_type;       /* Type {Heap,MMapped} */
  gt_mm_mode mode;        /* Memory mode {R,W,R/W} */
  /* Memory Data */
  uint64_t allocated;     /* Total allocated Bytes*/
  void* memory;           /* Pointer to block of memory */
  void* cursor;           /* Pointer current position of memory */
  /* Mapped File Memory*/
  int fd;                 /* File descriptor */
  char *file_name;        /* File name */
} gt_mm;

// Checkers
#define GT_MM_CHECK(mm) \
  GT_NULL_CHECK((mm)); \
  GT_NULL_CHECK((mm)->memory); \
  GT_NULL_CHECK((mm)->cursor); \
  GT_MM_CHECK_SEGMENT(mm)

#define GT_MM_CHECK_SEGMENT(mm) \
  gt_fatal_check( ((mm)->cursor < (mm)->memory) || \
                  ((mm)->cursor >= (mm)->memory+(mm)->allocated),MEM_CURSOR_OUT_OF_SEGMENT)

#define GT_MM_CHECK_ALIGNMENT(mm,align) \
  gt_fatal_check(!GT_MM_MEM_IS_ALIGNED((mm)->cursor,align),MEM_ALG_FAILED);

// Memory/MFile Allocators/Handlers
GT_INLINE gt_mm* gt_mm_bulk_malloc(const uint64_t num_bytes,const bool init_mem);
GT_INLINE gt_mm* gt_mm_bulk_mmalloc(const uint64_t num_bytes,const bool use_huge_pages);
GT_INLINE gt_mm* gt_mm_bulk_mmap_file(char* const file_name,const gt_mm_mode mode,const bool populate_page_tables);
GT_INLINE gt_mm* gt_mm_bulk_mmalloc_temp(const uint64_t num_bytes);
GT_INLINE void gt_mm_realloc(gt_mm* const mm,const uint64_t num_bytes);
GT_INLINE void gt_mm_free(gt_mm* const mm);

GT_INLINE gt_mm* gt_mm_bulk_load_file(char* const file_name,const uint64_t num_threads);
GT_INLINE gt_mm* gt_mm_bulk_mload_file(char* const file_name,const uint64_t num_threads);

// Accessors
GT_INLINE void* gt_mm_get_mem(gt_mm* const mm);
GT_INLINE void* gt_mm_get_base_mem(gt_mm* const mm);
GT_INLINE gt_mm_mode gt_mm_get_mode(gt_mm* const mm);
GT_INLINE uint64_t gt_mm_get_allocated(gt_mm* const mm);
GT_INLINE char* gt_mm_get_mfile_name(gt_mm* const mm);

// Seek
GT_INLINE uint64_t gt_mm_get_current_position(gt_mm* const mm);
GT_INLINE bool gt_mm_eom(gt_mm* const mm);
GT_INLINE void gt_mm_seek(gt_mm* const mm,const uint64_t byte_position);

GT_INLINE void gt_mm_skip_forward(gt_mm* const mm,const uint64_t num_bytes);
GT_INLINE void gt_mm_skip_backward(gt_mm* const mm,const uint64_t num_bytes);
GT_INLINE void gt_mm_skip_uint64(gt_mm* const mm);
GT_INLINE void gt_mm_skip_uint32(gt_mm* const mm);
GT_INLINE void gt_mm_skip_uint16(gt_mm* const mm);
GT_INLINE void gt_mm_skip_uint8(gt_mm* const mm);
GT_INLINE void gt_mm_skip_align(gt_mm* const mm,const uint64_t num_bytes);
GT_INLINE void gt_mm_skip_align_16(gt_mm* const mm);
GT_INLINE void gt_mm_skip_align_32(gt_mm* const mm);
GT_INLINE void gt_mm_skip_align_64(gt_mm* const mm);
GT_INLINE void gt_mm_skip_align_128(gt_mm* const mm);
GT_INLINE void gt_mm_skip_align_512(gt_mm* const mm);
GT_INLINE void gt_mm_skip_align_1024(gt_mm* const mm);
GT_INLINE void gt_mm_skip_align_4KB(gt_mm* const mm);
GT_INLINE void gt_mm_skip_align_mempage(gt_mm* const mm);

// Read
#define gt_mm_read(mm,var) gt_mm_copy_mem(mm,&var,sizeof(var))
GT_INLINE uint64_t gt_mm_read_uint64(gt_mm* const mm);
GT_INLINE uint32_t gt_mm_read_uint32(gt_mm* const mm);
GT_INLINE uint16_t gt_mm_read_uint16(gt_mm* const mm);
GT_INLINE uint8_t gt_mm_read_uint8(gt_mm* const mm);
GT_INLINE void* gt_mm_read_mem(gt_mm* const mm,const uint64_t num_bytes);
GT_INLINE void gt_mm_copy_mem(gt_mm* const mm,void* const dst,const uint64_t num_bytes);
GT_INLINE void gt_mm_copy_mem_parallel(gt_mm* const mm,void* const dst,const uint64_t num_bytes,const uint64_t num_threads);

// Write
#define gt_mm_write(mm,var) gt_mm_write_mem(mm,&var,sizeof(var))
GT_INLINE void gt_mm_write_uint64(gt_mm* const mm,const uint64_t data);
GT_INLINE void gt_mm_write_uint32(gt_mm* const mm,const uint32_t data);
GT_INLINE void gt_mm_write_uint16(gt_mm* const mm,const uint16_t data);
GT_INLINE void gt_mm_write_uint8(gt_mm* const mm,const uint8_t data);
GT_INLINE void gt_mm_write_mem(gt_mm* const mm,void* const src,const uint64_t num_bytes);

/*
 * Error Messages
 */
#define GT_ERROR_MEM_HANDLER "Could not allocate handler"
#define GT_ERROR_MEM_ALLOC "Could not allocate memory (%"PRIu64" requested)"
#define GT_ERROR_MEM_CALLOC "Could not allocate memory (%"PRIu64"x%"PRIu64" requested)"
#define GT_ERROR_MEM_REALLOC "Could not re-allocate memory (%"PRIu64" requested)"
#define GT_ERROR_MEM_ALLOC_DISK "Requested %"PRIu64" Bytes. Resorting to disk"
#define GT_ERROR_MEM_ALLOC_MMAP_DISK_FAIL "Requested %"PRIu64" Bytes. Failed to mmap memory to '%s'"
#define GT_ERROR_MEM_ALLOC_MMAP_FAIL "Requested %"PRIu64" Bytes. Failed to mmap memory"
#define GT_ERROR_MEM_CURSOR_OUT_OF_SEGMENT "Current memory cursor is out of boundaries (Segmentation fault)"
#define GT_ERROR_MEM_CURSOR_SEEK "Could not seek to address %"PRIu64". Out of boundaries (Segmentation fault)"
#define GT_ERROR_MEM_ALG_FAILED "Failed aligning the memory address to the specified boundary"

#endif /* GT_MEMORY_MANAGEMENT_H_ */
