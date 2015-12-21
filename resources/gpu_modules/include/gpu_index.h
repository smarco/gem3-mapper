#ifndef GPU_INDEX_H_
#define GPU_INDEX_H_

#include "gpu_commons.h"
#include "gpu_devices.h"

/* Defines related to BWT representation */
#define GPU_FMI_BITS_LOAD_PER_THREAD  128     // accesses of 16 Bytes / thread

#define GPU_FMI_ENTRY_LENGTH        ((GPU_FMI_COUNTERS_PER_ENTRY * GPU_UINT64_LENGTH) \
                                    + (GPU_FMI_ENTRY_SIZE * GPU_FMI_BWT_CHAR_LENGTH))       //   512 bits
#define GPU_FMI_THREADS_PER_ENTRY   (GPU_FMI_ENTRY_LENGTH / GPU_FMI_BITS_LOAD_PER_THREAD)   //   4 threads
#define GPU_FMI_ENTRIES_PER_WARP    (GPU_WARP_SIZE / GPU_FMI_THREADS_PER_ENTRY)
#define GPU_FMI_ENTRIES_PER_BLOCK   (GPU_THREADS_PER_BLOCK / GPU_FMI_THREADS_PER_ENTRY)
#define GPU_FMI_SLICES_PER_THREAD   (GPU_FMI_BITS_LOAD_PER_THREAD / GPU_UINT32_LENGTH)      // slices of 16 Bytes / thread


/*****************************
Internal Objects (General)
*****************************/

typedef struct {
  uint64_t        bwtSize;
  uint64_t        numEntries;
  gpu_fmi_entry_t *h_fmi;
  gpu_fmi_entry_t **d_fmi;
  memory_alloc_t  *memorySpace;
  gpu_module_t    activeModules;
} gpu_index_buffer_t;


/*************************************
Specific types for the Devices (GPUs)
**************************************/

// Just to cast and exchange memory between threads for the shared mem
typedef union {
  uint4 s;
  uint32_t v[GPU_FMI_SLICES_PER_THREAD];
} gpu_fmi_exch_bmp_mem_t;

typedef union {
  uint4 v[GPU_FMI_THREADS_PER_ENTRY];
} gpu_fmi_device_entry_t;

typedef struct {
  uint32_t bitmaps[GPU_FMI_BWT_CHAR_LENGTH];
} gpu_index_bitmap_entry_t;

typedef struct {
  uint64_t counters[GPU_FMI_NUM_COUNTERS];
} gpu_index_counter_entry_t;


/* Functions to initialize the index data on the DEVICE */
gpu_error_t gpu_init_index_dto(gpu_index_buffer_t* const index);
gpu_error_t gpu_init_index(gpu_index_buffer_t **index, const char* const indexRaw,
                           const uint64_t bwtSize, const gpu_index_coding_t indexCoding,
                           const uint32_t numSupportedDevices, const gpu_module_t activeModules);
gpu_error_t gpu_transfer_index_CPU_to_GPUs(gpu_index_buffer_t* const index, gpu_device_info_t** const devices);

/* Transform index functions */
gpu_error_t gpu_transform_index(const char* const indexRaw, gpu_index_buffer_t* const fmi, const gpu_index_coding_t indexCoding);

/* Primitives to build indexes */
gpu_error_t gpu_index_build_FMI(gpu_index_buffer_t* const fmi, gpu_index_bitmap_entry_t* const h_bitmap_BWT,
                                const gpu_index_counter_entry_t* const h_counters_FMI);

/* Stream index functions  */
gpu_error_t gpu_write_index(FILE* fp, const gpu_index_buffer_t* const index);
gpu_error_t gpu_read_index(FILE* fp, gpu_index_buffer_t* const index);

/* Functions to release the index data from the DEVICE & HOST */
gpu_error_t gpu_free_index_host(gpu_index_buffer_t* index);
gpu_error_t gpu_free_unused_index_host(gpu_index_buffer_t *index, gpu_device_info_t** const devices);
gpu_error_t gpu_free_index_device(gpu_index_buffer_t* index, gpu_device_info_t** const devices);
gpu_error_t gpu_free_index(gpu_index_buffer_t **index, gpu_device_info_t** const devices);

#include "gpu_io.h"
#endif /* GPU_INDEX_H_ */
