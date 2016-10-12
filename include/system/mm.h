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

#ifndef MEMORY_MANAGEMENT_H_
#define MEMORY_MANAGEMENT_H_

#include "system/commons.h"
#include "system/errors.h"

/*
 * Config
 */
#ifdef GEM_DEBUG
#define MM_NO_MMAP
#endif

/******************************************************************************
 * Memory Alignment Utils
 ******************************************************************************/
// Check Memory Alignment
#define MM_MEM_ALIGNED_MASK_16b   mm_mem_alignment_bits_mask[0]
#define MM_MEM_ALIGNED_MASK_32b   mm_mem_alignment_bits_mask[1]
#define MM_MEM_ALIGNED_MASK_64b   mm_mem_alignment_bits_mask[2]
#define MM_MEM_ALIGNED_MASK_128b  mm_mem_alignment_bits_mask[3]
#define MM_MEM_ALIGNED_MASK_256b  mm_mem_alignment_bits_mask[4]
#define MM_MEM_ALIGNED_MASK_512b  mm_mem_alignment_bits_mask[5]
#define MM_MEM_ALIGNED_MASK_1024b mm_mem_alignment_bits_mask[6]
/* (...) */
#define MM_MEM_ALIGNED_MASK_1KB   mm_mem_alignment_bits_mask[9]
#define MM_MEM_ALIGNED_MASK_2KB   mm_mem_alignment_bits_mask[10]
#define MM_MEM_ALIGNED_MASK_4KB   mm_mem_alignment_bits_mask[11]
#define MM_MEM_ALIGNED_MASK_2MB   mm_mem_alignment_bits_mask[20]
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
// Align Memory Macros
#define MM_MEM_ALIGNED_MASK_16b   mm_mem_alignment_bits_mask[0]
#define MM_MEM_ALIGNED_MASK_32b   mm_mem_alignment_bits_mask[1]
#define MM_MEM_ALIGNED_MASK_64b   mm_mem_alignment_bits_mask[2]
#define MM_MEM_ALIGNED_MASK_128b  mm_mem_alignment_bits_mask[3]
#define MM_MEM_ALIGNED_MASK_256b  mm_mem_alignment_bits_mask[4]
#define MM_MEM_ALIGNED_MASK_512b  mm_mem_alignment_bits_mask[5]
#define MM_MEM_ALIGNED_MASK_1024b mm_mem_alignment_bits_mask[6]
/* (...) */
#define MM_MEM_ALIGNED_MASK_1KB   mm_mem_alignment_bits_mask[9]
#define MM_MEM_ALIGNED_MASK_2KB   mm_mem_alignment_bits_mask[10]
#define MM_MEM_ALIGNED_MASK_4KB   mm_mem_alignment_bits_mask[11]
#define MM_MEM_ALIGNED_MASK_2MB   mm_mem_alignment_bits_mask[20]

#define MM_MEM_ALIGN_16b(mem)   MM_CAST_PTR((MM_CAST_ADDR(mem)+MM_MEM_ALIGNED_MASK_16b)&(~MM_MEM_ALIGNED_MASK_16b))
#define MM_MEM_ALIGN_32b(mem)   MM_CAST_PTR((MM_CAST_ADDR(mem)+MM_MEM_ALIGNED_MASK_32b)&(~MM_MEM_ALIGNED_MASK_32b))
#define MM_MEM_ALIGN_64b(mem)   MM_CAST_PTR((MM_CAST_ADDR(mem)+MM_MEM_ALIGNED_MASK_64b)&(~MM_MEM_ALIGNED_MASK_64b))
#define MM_MEM_ALIGN_128b(mem)  MM_CAST_PTR((MM_CAST_ADDR(mem)+MM_MEM_ALIGNED_MASK_128b)&(~MM_MEM_ALIGNED_MASK_128b))
#define MM_MEM_ALIGN_256b(mem)  MM_CAST_PTR((MM_CAST_ADDR(mem)+MM_MEM_ALIGNED_MASK_256b)&(~MM_MEM_ALIGNED_MASK_256b))
#define MM_MEM_ALIGN_512b(mem)  MM_CAST_PTR((MM_CAST_ADDR(mem)+MM_MEM_ALIGNED_MASK_512b)&(~MM_MEM_ALIGNED_MASK_512b))
#define MM_MEM_ALIGN_1024b(mem) MM_CAST_PTR((MM_CAST_ADDR(mem)+MM_MEM_ALIGNED_MASK_1024b)&(~MM_MEM_ALIGNED_MASK_1024b))
/* (...) */
#define MM_MEM_ALIGN_1KB(mem)   MM_CAST_PTR((MM_CAST_ADDR(mem)+MM_MEM_ALIGNED_MASK_1KB)&(~MM_MEM_ALIGNED_MASK_1KB))
#define MM_MEM_ALIGN_2KB(mem)   MM_CAST_PTR((MM_CAST_ADDR(mem)+MM_MEM_ALIGNED_MASK_2KB)&(~MM_MEM_ALIGNED_MASK_2KB))
#define MM_MEM_ALIGN_4KB(mem)   MM_CAST_PTR((MM_CAST_ADDR(mem)+MM_MEM_ALIGNED_MASK_4KB)&(~MM_MEM_ALIGNED_MASK_4KB))
#define MM_MEM_ALIGN_2MB(mem)   MM_CAST_PTR((MM_CAST_ADDR(mem)+MM_MEM_ALIGNED_MASK_2MB)&(~MM_MEM_ALIGNED_MASK_2MB))

/*
 * Temporal folder path
 */
#define MM_DEFAULT_TMP_FOLDER ""

char* mm_get_tmp_folder(void);
void mm_set_tmp_folder(char* const tmp_folder_path);

/******************************************************************************
 * UnitMemory
 *   Allocate relative small chunks of memory relying on the regular memory manager,
 *   usually malloc/calloc using a BuddySystem (Helper functions)
 ******************************************************************************/
#define mm_alloc(type) ((type*)mm_malloc_(1,sizeof(type),false,0))
#define mm_malloc(num_bytes) (mm_malloc_(1,num_bytes,false,0))
#define mm_calloc(num_elements,type,clear_mem) ((type*)mm_malloc_(num_elements,sizeof(type),clear_mem,0))
#define mm_malloc_uint64(variable_name,init_value) uint64_t* const variable_name = mm_malloc(sizeof(uint64_t)); *variable_name = init_value;
#define mm_malloc_uint32(variable_name,init_value) uint32_t* const variable_name = mm_malloc(sizeof(uint32_t)); *variable_name = init_value;
#define mm_malloc_uint16(variable_name,init_value) uint16_t* const variable_name = mm_malloc(sizeof(uint16_t)); *variable_name = init_value;
#define mm_malloc_uint8(variable_name,init_value)  uint8_t*  const variable_name = mm_malloc(sizeof(uint8_t));  *variable_name = init_value;
#define mm_malloc_int64(variable_name,init_value)  int64_t* const variable_name = mm_malloc(sizeof(int64_t)); *variable_name = init_value;
#define mm_malloc_int32(variable_name,init_value)  int32_t* const variable_name = mm_malloc(sizeof(int32_t)); *variable_name = init_value;
#define mm_malloc_int16(variable_name,init_value)  int16_t* const variable_name = mm_malloc(sizeof(int16_t)); *variable_name = init_value;
#define mm_malloc_int8(variable_name,init_value)   int8_t*  const variable_name = mm_malloc(sizeof(int8_t));  *variable_name = init_value;
void* mm_malloc_(uint64_t const num_elements,const uint64_t size_element,const bool init_mem,const int init_value);
void* mm_malloc_nothrow(uint64_t const num_elements,const uint64_t size_element,const bool init_mem,const int init_value);
void* mm_realloc(void* mem_addr,const uint64_t num_bytes);
void* mm_realloc_nothrow(void* mem_addr,const uint64_t num_bytes);
void mm_free(void* mem_addr);

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

/*
 * Checkers
 */
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
mm_t* mm_bulk_malloc(const uint64_t num_bytes,const bool init_mem);
mm_t* mm_bulk_mmalloc(const uint64_t num_bytes,const bool use_huge_pages);
mm_t* mm_bulk_mmalloc_temp(const uint64_t num_bytes);
void mm_bulk_free(mm_t* const mem_manager);

mm_t* mm_bulk_mmap_file(char* const file_name,const mm_mode mode,const bool populate_page_tables);
mm_t* mm_bulk_load_file(char* const file_name);
mm_t* mm_bulk_mload_file(char* const file_name);

// Accessors
void* mm_get_mem(mm_t* const mem_manager);
void* mm_get_base_mem(mm_t* const mem_manager);
mm_mode mm_get_mode(mm_t* const mem_manager);
uint64_t mm_get_allocated(mm_t* const mem_manager);
char* mm_get_mfile_name(mm_t* const mem_manager);

// Seek
uint64_t mm_get_current_position(mm_t* const mem_manager);
bool mm_eom(mm_t* const mem_manager);
void mm_seek(mm_t* const mem_manager,const uint64_t byte_position);

void mm_skip_forward(mm_t* const mem_manager,const uint64_t num_bytes);
void mm_skip_backward(mm_t* const mem_manager,const uint64_t num_bytes);
void mm_skip_uint64(mm_t* const mem_manager);
void mm_skip_uint32(mm_t* const mem_manager);
void mm_skip_uint16(mm_t* const mem_manager);
void mm_skip_uint8(mm_t* const mem_manager);
void mm_skip_align(mm_t* const mem_manager,const uint64_t num_bytes);
void mm_skip_align_16(mm_t* const mem_manager);
void mm_skip_align_32(mm_t* const mem_manager);
void mm_skip_align_64(mm_t* const mem_manager);
void mm_skip_align_128(mm_t* const mem_manager);
void mm_skip_align_512(mm_t* const mem_manager);
void mm_skip_align_1024(mm_t* const mem_manager);
void mm_skip_align_4KB(mm_t* const mem_manager);
void mm_skip_align_mempage(mm_t* const mem_manager);

// Read
#define mm_read(mem_manager,var) mm_copy_mem(mem_manager,&var,sizeof(var))
uint64_t mm_read_uint64(mm_t* const mem_manager);
uint32_t mm_read_uint32(mm_t* const mem_manager);
uint16_t mm_read_uint16(mm_t* const mem_manager);
uint8_t mm_read_uint8(mm_t* const mem_manager);
void* mm_read_mem(mm_t* const mem_manager,const uint64_t num_bytes);
void mm_copy_mem(
    mm_t* const mem_manager,
    void* const dst,
    const uint64_t num_bytes);

// Write
#define mm_write(mem_manager,var) mm_write_mem(mem_manager,&var,sizeof(var))
void mm_write_uint64(mm_t* const mem_manager,const uint64_t data);
void mm_write_uint32(mm_t* const mem_manager,const uint32_t data);
void mm_write_uint16(mm_t* const mem_manager,const uint16_t data);
void mm_write_uint8(mm_t* const mem_manager,const uint8_t data);
void mm_write_mem(mm_t* const mem_manager,void* const src,const uint64_t num_bytes);

// Status
int64_t mm_get_page_size(void);
int64_t mm_get_mem_available_virtual(void);
int64_t mm_get_mem_available_cached(void);
int64_t mm_get_mem_available_free(void);
int64_t mm_get_mem_available_total(void);
int64_t mm_get_mem_total(void);

#endif /* MEMORY_MANAGEMENT_H_ */
