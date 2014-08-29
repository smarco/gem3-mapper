/*
 * PROJECT: GEMMapper
 * FILE: bpm_align_gpu.c
 * DATE: 06/06/2012
 * AUTHOR(S): Alejandro Chacon <alejandro.chacon@uab.es>
 *            Santiago Marco-Sola <santiagomsola@gmail.com>
 */

/*
 * TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
 *   1.- Common prefix to constants like NUM_SUB_ENTRIES or accessors
 *   2.- Remove duplicates of structures definition. E.g qryEntry_t
 *   3.- Remove redefinitions of structures
 *   4.- Check @bpm_gpu_init::initMyers
 *     4.1 Ref conversion
 *     4.2 ARCHITECTURE_TYPE call
 *   5.- I need to understand the reason to store the pattern (qryEntry_t[]) in such interleaved fashion
 *     5.1 Check adaptor
 *   6.- Different CUDA-cards but, different buffer sizes? (Are buffers homogeneous?)
 *     NO
 *   7.- Lost of performance due to underused buffers
 *   8.- MACRO_VAR to get CUDA_COMPATIBLE
 * TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
 *
 *  A. Autoconfig to detect CUDA_NVCC
 *  B. PREFIXES for API
 *  C. Function to detect CUDA card (to be used at main() to reject --cuda parameter(s))
 *  D. Link Nsight remote run
 *  E. Vtune Amplifier non-commercial
 *
 */

#include "bpm_align_gpu.h"
#include "../resources/myers_gpu/myers-interface.h"

#define BPM_GPU_BUFFER_SIZE BUFFER_SIZE_32M

/*
 * BPM-GPU Setup
 */
GEM_INLINE bpm_gpu_buffer_collection_t* bpm_gpu_init(
    const dna_text_t* const enc_text,const uint32_t num_buffers,
    const int32_t average_query_size,const int32_t candidates_per_query) {
  GEM_CHECK_POSITIVE(average_query_size); // FIXME More here
  GEM_CHECK_POSITIVE(candidates_per_query);
  // Allocate Buffer Collection
  bpm_gpu_buffer_collection_t* const buffer_collection = mm_alloc(bpm_gpu_buffer_collection_t);
  buffer_collection->bpm_gpu_buffers = mm_calloc(num_buffers,bpm_gpu_buffer_t,true);
  buffer_collection->num_buffers = num_buffers;
  // Initialize Myers
  const char* const text = (const char* const) dna_text_get_buffer(enc_text);
  const uint64_t text_length = dna_text_get_length(enc_text);
  initMyers(&buffer_collection->internal_buffers,num_buffers,CONVERT_B_TO_MB(BPM_GPU_BUFFER_SIZE),
      text,GEM,text_length,average_query_size,candidates_per_query,ARCH_SUPPORTED);
  // Initialize Buffers
  uint64_t i;
  for (i=0;i<num_buffers;++i) {
    bpm_gpu_buffer_t* const bpm_gpu_buffer = buffer_collection->bpm_gpu_buffers + i;
    bpm_gpu_buffer->buffer = buffer_collection->internal_buffers[i];
    bpm_gpu_buffer->num_PEQ_entries = 0;
    bpm_gpu_buffer->num_queries = 0;
    bpm_gpu_buffer->num_candidates = 0;
    bpm_gpu_buffer->pattern_id = 0;
    bpm_gpu_buffer->candidates_left = 0;
  }
  // Return
  return buffer_collection;
}
GEM_INLINE void bpm_gpu_destroy(bpm_gpu_buffer_collection_t* const buffer_collection) {
  // Destroy (Myers)
  endMyers(&buffer_collection->internal_buffers);
  // Free HUB
  mm_free(buffer_collection->bpm_gpu_buffers);
  mm_free(buffer_collection);
}
/*
 * Buffer Accessors
 */
GEM_INLINE uint64_t bpm_gpu_buffer_get_num_candidates(bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  return bpm_gpu_buffer->num_candidates;
}
GEM_INLINE uint64_t bpm_gpu_buffer_get_num_queries(bpm_gpu_buffer_t* const bpm_gpu_buffer) {
  return bpm_gpu_buffer->num_queries;
}

#define NUM_SUB_ENTRIES (128/32)  // FIXME FIXME NUM_SUB_ENTRIES => get_limit()
#define BPM_GPU_PATTERN_NUM_SUB_ENTRIES NUM_SUB_ENTRIES
#define BPM_GPU_PATTERN_ENTRY_LENGTH    (NUM_SUB_ENTRIES*UINT32_LENGTH)
#define BPM_GPU_PATTERN_ENTRY_SIZE      (BPM_GPU_PATTERN_ENTRY_LENGTH/UINT8_SIZE)
#define BPM_GPU_PATTERN_ALPHABET_LENGTH PEQ_ALPHABET_SIZE

GEM_INLINE bool bpm_gpu_buffer_fits_in_buffer(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,
    const uint64_t num_patterns,const uint64_t total_pattern_length,const uint64_t total_candidates) {
  // Calculate dimensions
  const uint64_t pattern_PEQ_num_entries = DIV_CEIL(total_pattern_length,BPM_GPU_PATTERN_ENTRY_LENGTH);
  // Get Limits
  const uint64_t max_PEQ_entries = getMaxPEQEntries(bpm_gpu_buffer->buffer);
  const uint64_t max_queries = getMaxQueries(bpm_gpu_buffer->buffer);
  // Check available space in buffer for the pattern
  if (bpm_gpu_buffer->num_queries+num_patterns > max_queries ||
      bpm_gpu_buffer->num_PEQ_entries+pattern_PEQ_num_entries > max_PEQ_entries) {
    // Check if the pattern can fit into an empty buffer
    gem_cond_fatal_error(pattern_PEQ_num_entries > max_PEQ_entries,
        BPM_GPU_MAX_PATTERN_LENGTH,pattern_PEQ_num_entries,max_PEQ_entries);
    return false;
  }
  // Check available space in buffer for the candidates
  const uint64_t max_candidates = getMaxCandidates(bpm_gpu_buffer->buffer);
  if (bpm_gpu_buffer->num_candidates+total_candidates > max_candidates) {
    // Check if the pattern can fit into an empty buffer
    gem_cond_fatal_error(total_candidates > max_candidates,
        BPM_GPU_MAX_CANDIDATES,total_candidates,max_candidates);
    return false;
  }
  // Ok, go on
  return true;
}
GEM_INLINE void bpm_gpu_buffer_put_pattern(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,const bpm_pattern_t* const bpm_pattern) {
  // Calculate PEQ dimensions
  const uint64_t pattern_num_words = bpm_pattern->pattern_num_words;
  const uint64_t pattern_PEQ_length = bpm_pattern->PEQ_length;
  const uint64_t pattern_PEQ_num_entries = DIV_CEIL(pattern_PEQ_length,BPM_GPU_PATTERN_ENTRY_LENGTH);
  /*
   * FIXME
   *  1.- Is max_PEQ_entries nominal space for |alphabet|*max_PEQ_entries == sizeof(qryEntry_t)?
   *  2.- Check ratios of allocated space. i.e. N*max_PEQ_entries <-> M*max_queries <-> Z*max_candidates
   */
  // Insert query metadata
  bpm_gpu_buffer->pattern_id = bpm_gpu_buffer->num_queries;
  (bpm_gpu_buffer->num_queries)++;
  qryInfo_t* const query_info = getPEQInfoBuffer(bpm_gpu_buffer->buffer) + bpm_gpu_buffer->pattern_id;
  query_info->posEntry = bpm_gpu_buffer->num_PEQ_entries;
  query_info->size = pattern_PEQ_num_entries; // FIXME Revise this
  // Insert query pattern
  qryEntry_t* const query_pattern = getPEQBuffer(bpm_gpu_buffer->buffer) + bpm_gpu_buffer->num_PEQ_entries;
  const uint8_t* const pattern_PEQ = bpm_pattern->PEQ;
  uint64_t entry, words8_offset;
  // Iterate over all entries
  for (entry=0,words8_offset=0;entry<pattern_PEQ_num_entries;++entry) {
    qryEntry_t* const query_pattern_entry = query_pattern + entry;
    // Iterate over all sub-entries
    uint64_t enc_char, sub_entry;
    for (enc_char=0;enc_char<BPM_GPU_PATTERN_ALPHABET_LENGTH;++enc_char) {
      for (sub_entry=0;sub_entry<BPM_GPU_PATTERN_NUM_SUB_ENTRIES;++sub_entry) {
        uint8_t* const bitmap = (uint8_t*) &(query_pattern_entry->bitmap[enc_char][sub_entry]);
        // Fill subentry bitmap
        uint64_t subentry_wp;
        for (subentry_wp=0;subentry_wp<4;++subentry_wp) {
          const uint64_t word8_pos = words8_offset + subentry_wp;
          if (word8_pos < pattern_num_words) {
            bitmap[subentry_wp] = pattern_PEQ[BPM_PATTERN_PEQ_IDX(enc_char,word8_pos,pattern_num_words)];
          } else {
            bitmap[subentry_wp] = 0;
          }
        }
      }
    }
    words8_offset+=4; // Next entry
  }
  bpm_gpu_buffer->num_PEQ_entries += pattern_PEQ_num_entries;
}
GEM_INLINE void bpm_gpu_buffer_put_candidate(
    bpm_gpu_buffer_t* const bpm_gpu_buffer,
    const uint64_t candidate_text_position,const uint64_t candidate_length) {
  // Insert candidate
  candInfo_t* const query_candidate = getCandidatesBuffer(bpm_gpu_buffer->buffer) + bpm_gpu_buffer->num_candidates;
  query_candidate->query = bpm_gpu_buffer->pattern_id;
  query_candidate->position = candidate_text_position;
  query_candidate->size = candidate_length;
  ++(bpm_gpu_buffer->num_candidates);
}
