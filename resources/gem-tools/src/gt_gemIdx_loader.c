/*
 * PROJECT: GEM-Tools library
 * FILE: gt_gemIdx_loader.c
 * DATE: 01/02/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_gemIdx_loader.h"

GT_INLINE void gt_gemIdx_read_header(gt_mm* const mm) {
//  char gem_magic_number[6];
//  uint16_t gem_version;
//  // Read magic number
//  gem_magic_number[0] = gt_mm_read_uint8(mm);
//  gem_magic_number[1] = gt_mm_read_uint8(mm);
//  gem_magic_number[2] = gt_mm_read_uint8(mm);
//  gem_magic_number[3] = gt_mm_read_uint8(mm);
//  gem_magic_number[4] = gt_mm_read_uint8(mm);
//  gem_magic_number[5] = gt_mm_read_uint8(mm);
//  // Read version
//  gem_version = gt_mm_read_uint16(mm);
//  gt_debug("MagicNumber[%c%c%c%d%d%d] Version[%d]",
//      gem_magic_number[0],gem_magic_number[1],gem_magic_number[2],
//      gem_magic_number[3],gem_magic_number[4],gem_magic_number[5],gem_version);
  // In short
  gt_mm_skip_uint64(mm); // Read magic number + Read version
}

GT_INLINE void gt_gemIdx_load_archive(
    char* const index_file_name,gt_sequence_archive* const sequence_archive,const bool load_sequences) {
  GT_NULL_CHECK(index_file_name);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  /*
   * Load the file
   */
  gt_mm* const mm = gt_mm_bulk_mmap_file(index_file_name,GT_MM_READ_ONLY,false);
  /*
   * Read Multiarchive
   */
  gt_gemIdx_read_header(mm);
  gt_mm_skip_uint64(mm); // Multiarchive type
  const uint64_t parameters_num = gt_mm_read_uint64(mm); // Parameters Num
  if (parameters_num > 0) gt_mm_skip_forward(mm,parameters_num*8); // Parameters
  gt_mm_skip_uint64(mm); // Archives Num
  /*
   * Read Archive (We only read the first one (the only one ever generated))
   */
  gt_gemIdx_read_header(mm);
  gt_mm_skip_uint64(mm); // Index Type
  gt_mm_skip_uint64(mm); // Filter Type
  gt_mm_skip_uint64(mm); // Complement Type
  /*
   * Read locator
   */
  gt_gemIdx_read_header(mm);
  const uint64_t num_intervals = gt_mm_read_uint64(mm); // Num Intervals
  const uint64_t tags_cumulative_length = gt_mm_read_uint64(mm); // Tags cumulative length
  gem_loc_t* const seq_info = gt_mm_read_mem(mm,num_intervals*sizeof(gem_loc_t));
  char* const tags = gt_mm_read_mem(mm,tags_cumulative_length);
  // Add sequences to the sequence archive
  gt_segmented_sequence* seg_seq = NULL;
  uint64_t i, last_offset=UINT64_MAX, seq_length;
  for (i=0;i<num_intervals;++i) {
    // if (i%10000==0) { gt_debug("Loaded %lu/%lu",i,num_intervals); }
    // gt_debug("#%lu :: %s offset=%ld [bot,top]=[%lu,%lu]",i,
    //    tags+seq_info[i].tag_offset-1,seq_info[i].sequence_offset,seq_info[i].bot,seq_info[i].top);
    if (seq_info[i].tag_offset <= 0) continue; // Skip negative strand sequences
    if (seq_info[i].tag_offset!=last_offset) { // New sequence
      if (last_offset!=UINT64_MAX) { // Close last sequence
        seg_seq->sequence_total_length = seq_length;
      }
      seg_seq = gt_segmented_sequence_new();
      gt_string_set_string(seg_seq->seq_name,tags+seq_info[i].tag_offset-1);
      gt_sequence_archive_add_bed_sequence(sequence_archive,seg_seq);
      last_offset=seq_info[i].tag_offset;
    }
    // Keep sequence length
    seq_length = seq_info[i].sequence_offset+(seq_info[i].top-seq_info[i].bot);
    // Store interval
    gt_vector* const interval_vector =
        gt_sequence_archive_get_bed_intervals_vector_dyn(sequence_archive,gt_string_get_string(seg_seq->seq_name));
    gt_vector_insert(interval_vector,seq_info[i],gem_loc_t);
  }
  if (last_offset!=UINT64_MAX) seg_seq->sequence_total_length = seq_length; // Close last sequence
  // Check if we need to read the sequences themselves
  if (load_sequences) {
    gt_mm_skip_align_64(mm); // Align to 64bits
    /*
     * Read FMI
     */
    gt_gemIdx_read_header(mm);
    gt_mm_skip_uint64(mm); // Index Type
    const uint64_t text_length = gt_mm_read_uint64(mm); // Text Length
    const uint64_t sampling_rate = gt_mm_read_uint64(mm); // Sampling Rate
    const uint64_t counter_bytes = gt_mm_read_uint64(mm); // Counter Bytes
    const uint64_t csa_size = counter_bytes*((text_length+1+(1<<sampling_rate))/(1<<sampling_rate));
    gt_mm_skip_forward(mm,8+2*csa_size); // D+I
    /*
     * Read BWT
     */
    gt_mm_skip_align_64(mm); // Align to 64bits
    gt_gemIdx_read_header(mm);
    const uint64_t bwt_length = gt_mm_read_uint64(mm); // BWT Length
    const uint64_t num_blocks = (bwt_length+256)/256;
    const uint64_t bwt_size = (7*8) + (num_blocks*224);
    gt_mm_skip_forward(mm,bwt_size);
    /*
     * Read BED
     */
    gt_mm_skip_align_64(mm); // Align to 64bits
    gt_gemIdx_read_header(mm);
    const uint64_t bed_length = gt_mm_read_uint64(mm); // BED Length
    const uint64_t bed_size = 8*((bed_length+63)/64)*3;

    /*
     * Dump BED into the sequence archive
     */
    // uint64_t* const ptr_block = gt_mm_read_mem(mm,bed_size);
    sequence_archive->mm = gt_mm_bulk_mmalloc(bed_size,false);
    sequence_archive->bed = gt_mm_get_base_mem(sequence_archive->mm);
    gt_fm_bulk_read_file(index_file_name,sequence_archive->bed,gt_mm_get_current_position(mm),bed_size);
  }
  // Free MM
  gt_mm_free(mm);
}

/*
 * Retrieve sequences from GEMindex
 */
GT_INLINE int64_t gt_gemIdx_get_bed_sequence_string(
  gt_sequence_archive* const sequence_archive,char* const seq_id,
  const uint64_t position,const uint64_t length,gt_string* const string) {
  GT_SEQUENCE_BED_ARCHIVE_CHECK(sequence_archive);
  // Locate BED position
  const int64_t seq_position = position; // Guarantee signed arithmetic
  gt_vector* const intervals_vector = gt_sequence_archive_get_bed_intervals_vector(sequence_archive,seq_id);
  if (intervals_vector==NULL) return GT_GEMIDX_SEQ_NOT_FOUND;
  gem_loc_t* const intervals = gt_vector_get_mem(intervals_vector,gem_loc_t);
  uint64_t inf = 0, sup = gt_vector_get_used(intervals_vector)-1;
  while (sup != inf) {
    const uint64_t mid = (sup+inf)/2;
    if (seq_position > intervals[mid].sequence_offset+(intervals[mid].top-intervals[mid].bot)) {
      inf = mid+1;
    } else {
      sup = mid;
    }
  }
  if (seq_position < intervals[sup].sequence_offset ||
      seq_position >= intervals[sup].sequence_offset+(intervals[sup].top-intervals[sup].bot)) return GT_GEMIDX_INTERVAL_NOT_FOUND;
  // Get BED position
  const uint64_t bed_position = intervals[sup].bot + (seq_position-intervals[sup].sequence_offset);
  // Clean string
  gt_string_clear(string);
  // Decode sequence string
  const uint64_t be_block_pos = (bed_position/GT_CDNA_BLOCK_CHARS)*GT_CDNA_BLOCK_BITMAPS;
  uint64_t be_block_mod = bed_position%GT_CDNA_BLOCK_CHARS;
  uint64_t *ptr_block = sequence_archive->bed + be_block_pos;
  uint64_t chunk_2 = *ptr_block>>be_block_mod; ++ptr_block;
  uint64_t chunk_1 = *ptr_block>>be_block_mod; ++ptr_block;
  uint64_t chunk_0 = *ptr_block>>be_block_mod; ++ptr_block;
  uint64_t i;
  for (i=0;i<length;++i,++be_block_mod) {
    if (gt_expect_false(be_block_mod==GT_CDNA_BLOCK_CHARS)) {
      be_block_mod=0;
      chunk_2 = *ptr_block; ++ptr_block;
      chunk_1 = *ptr_block; ++ptr_block;
      chunk_0 = *ptr_block; ++ptr_block;
      GT_PREFETCH(ptr_block);
    }
    const uint8_t enc_c = (((chunk_0&GT_CDNA_EXTRACT_MASK)) |
                           ((chunk_1&GT_CDNA_EXTRACT_MASK)<<1) |
                           ((chunk_2&GT_CDNA_EXTRACT_MASK)<<2));
    if (gt_expect_false(enc_c > GT_CDNA_ENC_CHAR_N)) {
      break;
    } else {
      gt_string_append_char(string,gt_cdna_decode[enc_c]);
      chunk_0>>=1; chunk_1>>=1; chunk_2>>=1;
    }
  }
  gt_string_append_eos(string);
  return i;
}
