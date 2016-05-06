/*
 * PROJECT: GEMMapper
 * FILE: sa_builder.c
 * DATE: 07/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Implements a data storage structure as to classify SA positions with respect
 *   to it's stating k-mer. SA positions are stored in blocks of a given size
 */

#include "fm_index/sa_builder.h"
#include "data_structures/dna_text.h"
#include "stats/stats_vector.h"

/*
 * Errors
 */
#define GEM_ERROR_SA_BUILDER_LIMIT_MAX_BLOCK_MEM "SA Builder. Maximum SA-bucket (%"PRIu64" GB) larger than max-memory-block allowed (%"PRIu64" GB)"
#define GEM_ERROR_SA_BUILDER_LIMIT_MAX_BLOCK_ADDRESSED "SA Builder. Maximum SA-block (%"PRIu64" GB) cannot be addressed (%"PRIu64" GB)"
#define GEM_ERROR_SA_BUILDER_STORE_BLOCK_NOT_FILLED "SA Builder. K-mer block not filled according to previous counting"
#define GEM_ERROR_SA_BUILDER_SEQUENCE_MIN_LENGTH "SA Builder. Index total length (%"PRIu64") is below minimum threshold (%"PRIu64")"

/*
 * Debug
 */
#define SA_BUILDER_DEBUG_SPLIT_BLOCKS           true
#define SA_BUILDER_DUMP_BWT_NUM_BLOCKS         10000

/*
 * Global
 */
sa_builder_t* global_sa_builder;
uint64_t global_enc_text_length;
const uint8_t* global_enc_text;
uint8_t* global_enc_bwt;
sampled_sa_builder_t* global_sampled_sa;
const uint8_t* ds_shallow_text_limit;
// Stats Ranges
uint64_t kmer_count_range[] = {0,10,100,1000,10000,100000};
uint64_t suffix_cmp_length_range[] = {0,1,2,3,4,5,6,7,8,9,10,100,1000,10000};

/*
 * Constants
 */
#define SA_BUILDER_STORE_SUFFIXES_TICKER_STEP 10000000

/*
 * Prototypes
 */
void sa_builder_record_kmer_count_stats(sa_builder_t* const sa_builder);

/*
 * Setup
 */
sa_builder_t* sa_builder_new(
    char* const name_prefix,
    dna_text_t* const enc_text,
    const uint64_t num_threads,
    const uint64_t max_memory) {
  // Allocate sa_builder
  sa_builder_t* const sa_builder = mm_alloc(sa_builder_t);
  // Prepare the text
  sa_builder->enc_text = enc_text;
  sa_builder->name_prefix = name_prefix;
  // Fill the circular k-mers positions
  const uint64_t text_length = dna_text_get_length(sa_builder->enc_text);
  uint8_t* const enc_text_buffer = dna_text_get_text(sa_builder->enc_text);
  uint64_t i;
  enc_text_buffer[-2] = enc_text_buffer[text_length-2];
  enc_text_buffer[-1] = enc_text_buffer[text_length-1];
  for (i=0;i<SA_BWT_PADDED_LENGTH;++i) {
    enc_text_buffer[text_length+i] = enc_text_buffer[i];
  }
  // K-mer counting/splitting
  sa_builder->num_kmers = pow(SA_BUILDER_ALPHABET_SIZE,SA_BUILDER_KMER_LENGTH);
  sa_builder->kmer_count = mm_calloc(sa_builder->num_kmers,uint64_t,true);
  // Threads
  sa_builder->pthreads = mm_calloc(num_threads,pthread_t,false);
  sa_builder->num_threads = num_threads;
  sa_builder->max_thread_memory = max_memory/num_threads;
  // Stats
  sa_builder->kmer_count_stats = stats_vector_customed_range_new(kmer_count_range,5,100000);
  sa_builder->group_size_stats = stats_vector_customed_range_new(kmer_count_range,5,100000);
  sa_builder->suffix_cmp_length = stats_vector_customed_range_new(suffix_cmp_length_range,13,10000);
  // Misc
  ticker_mutex_enable(&sa_builder->ticker);
  // Return
  return sa_builder;
}
void sa_builder_delete(sa_builder_t* const sa_builder) {
  mm_free(sa_builder->kmer_count);
  mm_free(sa_builder->pthreads);
  gem_unlink(sa_builder->sa_positions_file_name);
  mm_free(sa_builder->sa_positions_file_name);
  fm_close(sa_builder->sa_positions_file);
  mm_free(sa_builder->sa_groups);
  stats_vector_delete(sa_builder->kmer_count_stats);
  stats_vector_delete(sa_builder->group_size_stats);
  stats_vector_delete(sa_builder->suffix_cmp_length);
  mm_free(sa_builder);
}
/*
 * 1.- Count all suffixes
 */
void sa_builder_count_suffixes(
    sa_builder_t* const sa_builder,
    uint64_t* const character_occurrences,
    const bool verbose) {
  // Ticker
  ticker_t ticker;
  ticker_percentage_reset(&ticker,verbose,"Building-BWT::Counting K-mers",1,1,true);
  // Init
  const uint64_t text_length = dna_text_get_length(sa_builder->enc_text);
  gem_cond_fatal_error(text_length < SA_BUILDER_KMER_LENGTH,
      SA_BUILDER_SEQUENCE_MIN_LENGTH,text_length,(uint64_t)SA_BUILDER_KMER_LENGTH);
  const uint8_t* const enc_text = dna_text_get_text(sa_builder->enc_text);
  uint64_t* const kmer_count = sa_builder->kmer_count;
  uint64_t i;
  uint64_t kmer_idx = 0;
  // Fill initial k-mer index
  for (i=0;i<SA_BWT_CYCLIC_LENGTH;++i) {
    const uint8_t enc_char = enc_text[i];
    kmer_idx = (kmer_idx<<DNA_EXT_RANGE_BITS) | enc_char;
    ++(character_occurrences[enc_char]);
  }
  // Count suffixes of all text
  for (;i<text_length;++i) {
    const uint8_t enc_char = enc_text[i];
    kmer_idx = (kmer_idx<<DNA_EXT_RANGE_BITS) | enc_char;
    ++(kmer_count[SA_BUILDER_KMER_MASK_INDEX(kmer_idx)]);
    ++(character_occurrences[enc_char]);
  }
  const uint64_t extended_text_length = text_length + SA_BWT_CYCLIC_LENGTH;
  for (;i<extended_text_length;++i) {
    const uint8_t enc_char = enc_text[i];
    kmer_idx = (kmer_idx<<DNA_EXT_RANGE_BITS) | enc_char;
    ++(kmer_count[SA_BUILDER_KMER_MASK_INDEX(kmer_idx)]);
  }
  // Stats
  sa_builder_record_kmer_count_stats(sa_builder);
  // Finish ticker
  ticker_finish(&ticker);
}
/*
 * 2.- Add all suffixes
 */
uint64_t sa_builder_calculate_max_bucket_size(sa_builder_t* const sa_builder) {
  const uint64_t* const kmer_count = sa_builder->kmer_count;
  uint64_t max_bucket_size = 0, i;
  for (i=0;i<sa_builder->num_kmers;++i) { // Find the largest block
    max_bucket_size = (max_bucket_size < kmer_count[i]) ? kmer_count[i] : max_bucket_size;
  }
  return max_bucket_size*UINT64_SIZE;
}
uint64_t sa_builder_calculate_num_sa_groups(
    sa_builder_t* const sa_builder,
    const uint64_t block_size) {
  const uint64_t* const kmer_count = sa_builder->kmer_count;
  const uint64_t kmers_per_block = block_size/UINT64_SIZE;
  uint64_t num_sa_groups = 0, block_positions_acc = 0, i;
  for (i=0;i<sa_builder->num_kmers;++i) {
    // Check k-mer count (Bucket Type B)
    const uint64_t new_block_positions_acc = block_positions_acc + kmer_count[i];
    if (gem_expect_true(block_positions_acc == 0 || new_block_positions_acc <= kmers_per_block)) {
      block_positions_acc = new_block_positions_acc;
    } else  {
      ++num_sa_groups;
      block_positions_acc = kmer_count[i];
    }
  }
  return (block_positions_acc > 0) ? num_sa_groups+1 : num_sa_groups;
}
void sa_builder_sa_groups_prepare(
    sa_builder_t* const sa_builder,
    const uint64_t block_size) {
  // Setup groups
  uint64_t* const kmer_count = sa_builder->kmer_count;
  const uint64_t kmers_per_block = block_size/UINT64_SIZE;
  uint64_t num_sa_group = 0, block_positions_acc = 0, sa_offset = 0;
  uint64_t i;
  for (i=0;i<sa_builder->num_kmers;++i) {
    // Assign a group
    const uint64_t current_kmer_count = kmer_count[i];
    // Check k-mer count (Bucket Type B)
    const uint64_t new_block_positions_acc = block_positions_acc + current_kmer_count;
    if (gem_expect_true(block_positions_acc == 0 || new_block_positions_acc <= kmers_per_block)) {
      kmer_count[i] = num_sa_group;
      block_positions_acc = new_block_positions_acc;
    } else  {
      sa_builder->sa_groups[num_sa_group].sa_offset = sa_offset;
      sa_builder->sa_groups[num_sa_group].num_sa_positions = block_positions_acc;
      sa_builder->sa_groups[num_sa_group].sa_positions_file = fm_open_file(sa_builder->sa_positions_file_name,FM_WRITE);
      fm_seek(sa_builder->sa_groups[num_sa_group].sa_positions_file,sa_offset*UINT64_SIZE);
      sa_offset += block_positions_acc;
      ++num_sa_group;
      kmer_count[i] = num_sa_group;
      block_positions_acc = current_kmer_count;
    }
  }
  if (block_positions_acc > 0) {
    sa_builder->sa_groups[num_sa_group].sa_offset = sa_offset;
    sa_builder->sa_groups[num_sa_group].num_sa_positions = block_positions_acc;
    sa_builder->sa_groups[num_sa_group].sa_positions_file = fm_open_file(sa_builder->sa_positions_file_name,FM_WRITE);
    fm_seek(sa_builder->sa_groups[num_sa_group].sa_positions_file,sa_offset*UINT64_SIZE);
  }
}
void sa_builder_sa_groups_distribute(sa_builder_t* const sa_builder) {
  // Calculate number of groups per thread
  const uint64_t num_threads = sa_builder->num_threads;
  uint64_t i, num_groups_per_thread = DIV_CEIL(sa_builder->num_sa_groups,num_threads);
  // Split k-mers across threads
  uint64_t current_thread_responsible = 0, num_groups_assign_to_thread = 1;
  for (i=0;i<sa_builder->num_sa_groups;++i) {
    // Assign thread responsible
    sa_builder->sa_groups[i].thread_responsible = current_thread_responsible;
    // Check groups per thread
    if (num_groups_assign_to_thread < num_groups_per_thread) {
      ++num_groups_assign_to_thread;
    } else {
      num_groups_assign_to_thread = 1;
      if (current_thread_responsible+1 < num_threads) {
        ++current_thread_responsible;
        if (current_thread_responsible < num_threads) {
          const uint64_t remaining_groups = sa_builder->num_sa_groups-i;
          const uint64_t remaining_threads = num_threads-current_thread_responsible;
          num_groups_per_thread = DIV_CEIL(remaining_groups,remaining_threads);
        }
      }
    }
  }
}
void sa_builder_sa_groups_cleanup(sa_builder_t* const sa_builder) {
  uint64_t i;
  sa_group_t* const sa_groups = sa_builder->sa_groups;
  for (i=0;i<sa_builder->num_sa_groups;++i) {
    fm_close(sa_groups[i].sa_positions_file);
  }
}
void sa_builder_store_suffixes_prepare(sa_builder_t* const sa_builder) {
  // Calculate maximum bucket size (for sorting)
  const uint64_t max_bucket_size = sa_builder_calculate_max_bucket_size(sa_builder);
  gem_cond_fatal_error(sa_builder->max_thread_memory < max_bucket_size,SA_BUILDER_LIMIT_MAX_BLOCK_MEM,
      CONVERT_B_TO_GB(max_bucket_size),CONVERT_B_TO_GB(sa_builder->max_thread_memory));
  // Calculate SA-block size (for sorting)
  const uint64_t sa_length = dna_text_get_length(sa_builder->enc_text);
  const uint64_t preferred_block_size = (sa_length*UINT64_SIZE)/SA_BUILDER_NUM_WRITTERS;
  const uint64_t block_size = MIN(sa_builder->max_thread_memory,preferred_block_size);
  sa_builder->block_size = block_size;
  sa_builder->max_bucket_size = max_bucket_size;
  // Allocate SA-positions memory
  char* tmp_file_basename = gem_strbasename(sa_builder->name_prefix);
  char* tmp_file_path = gem_strcat(mm_get_tmp_folder(),tmp_file_basename);
  sa_builder->sa_positions_file_name = gem_strcat(tmp_file_path,".sa.tmp");
  sa_builder->sa_positions_file = fm_open_file(sa_builder->sa_positions_file_name,FM_WRITE); // size(sa_builder->sa_positions_file) = sa_length*UINT64_SIZE
  mm_free(tmp_file_basename);
  mm_free(tmp_file_path);
  // Prepare SA-Groups
  const uint64_t num_sa_groups = sa_builder_calculate_num_sa_groups(sa_builder,block_size);
  sa_builder->num_sa_groups = num_sa_groups;
  sa_builder->sa_groups = mm_calloc(num_sa_groups,sa_group_t,true);
  sa_builder_sa_groups_prepare(sa_builder,block_size); // Setup SA-Groups
  sa_builder_sa_groups_distribute(sa_builder);
}
void sa_builder_store_sa_pos(
    sa_builder_t* const sa_builder,
    sa_group_t* const group,
    const uint64_t sa_pos,
    const uint64_t kmer_idx) {
  // Add the suffix (Piggybacking the 2BWT + 6SuffixCache)
  fm_write_uint64(group->sa_positions_file,sa_pos | SA_COMPACTED_TEXT_MASK_PIGGYBACKING(kmer_idx));
}
void* sa_builder_store_suffixes_thread(const uint8_t thread_id) {
  const uint64_t text_length = dna_text_get_length(global_sa_builder->enc_text);
  const uint8_t* const enc_text = dna_text_get_text(global_sa_builder->enc_text);
  uint64_t i, kmer_idx=0, sa_pos=0, tp = 0;
  // Fill k-mer index
  kmer_idx = enc_text[text_length-2];
  kmer_idx = (kmer_idx<<DNA_EXT_RANGE_BITS) | enc_text[text_length-1];
  for (i=0;i<SA_BWT_CYCLIC_LENGTH;++i) {
    kmer_idx = (kmer_idx<<DNA_EXT_RANGE_BITS) | enc_text[i];
  }
  // Count suffixes of all text
  const uint64_t extended_text_length = dna_text_get_length(global_sa_builder->enc_text)+SA_BWT_CYCLIC_LENGTH;
  const uint64_t* const kmer_count = global_sa_builder->kmer_count;
  sa_group_t* const sa_groups = global_sa_builder->sa_groups;
  for (sa_pos=0;i<extended_text_length;++i,++sa_pos,++tp) {
    kmer_idx = (kmer_idx<<DNA_EXT_RANGE_BITS) | enc_text[i];
    sa_group_t* const group = sa_groups + kmer_count[SA_BUILDER_KMER_MASK_INDEX(kmer_idx)];
    if (group->thread_responsible == thread_id) {
      sa_builder_store_sa_pos(global_sa_builder,group,sa_pos,kmer_idx);
    }
    // Ticker update
    if (thread_id==0 && tp==SA_BUILDER_STORE_SUFFIXES_TICKER_STEP) {
      ticker_update(&global_sa_builder->ticker,tp); tp = 0;
    }
  }
  // Return
  return NULL;
}
void sa_builder_store_suffixes(
    sa_builder_t* const sa_builder,
    const bool verbose) {
  /*
   * Prepare storage for Writing
   */
  global_sa_builder = sa_builder; // Set global data
  sa_builder_store_suffixes_prepare(sa_builder); // Begin
  // DEBUG
  sa_builder_display_stats(gem_info_get_stream(),sa_builder,true);
  /*
   * Launch Writing Threads
   */
  const uint64_t ticker_max = dna_text_get_length(sa_builder->enc_text)+SA_BWT_CYCLIC_LENGTH;
  ticker_percentage_reset(&sa_builder->ticker,
      verbose,"Building-BWT::Generating SA-Positions",
      ticker_max,SA_BUILDER_STORE_SUFFIXES_TICKER_STEP,true);
  const uint64_t num_threads = sa_builder->num_threads;
  uint64_t i;
  for (i=0;i<num_threads;++i) {
    // Launch thread
    gem_cond_fatal_error(pthread_create(sa_builder->pthreads+i,0,
        (void* (*)(void*))sa_builder_store_suffixes_thread,(void*)(i)),SYS_THREAD_CREATE);
  }
  for (i=0;i<num_threads;++i) {
    gem_cond_fatal_error(pthread_join(sa_builder->pthreads[i],0),SYS_THREAD_JOIN);
  }
  sa_builder_sa_groups_cleanup(sa_builder); // Flush and free FM handlers
  ticker_finish(&sa_builder->ticker);
}
/*
 * 3.1 Sort all suffixes (Subsidiary sort)
 */
// Basic Sort Suffixes Functions
int sa_builder_suffix_cmp(const uint64_t* const a,const uint64_t* const b) {
  // Eliminate the encoded text (PiggyBack)
  const uint64_t sa_pos_a = SA_POS_MASK_POSITION(*a);
  const void* a1 = global_enc_text+sa_pos_a;
  PREFETCH(a1);
  const uint64_t sa_pos_b = SA_POS_MASK_POSITION(*b);
  const void* a2 = global_enc_text+sa_pos_b;
  PREFETCH(a2);
  // Compare the 6 next characters (PiggyBack)
  const uint64_t cached_text_sa_pos_a = SA_POS_MASK_SUFFIX_CACHE(*a);
  const uint64_t cached_text_sa_pos_b = SA_POS_MASK_SUFFIX_CACHE(*b);
  const uint64_t cmp_cached = cached_text_sa_pos_a - cached_text_sa_pos_b;
  if (cmp_cached) {
    return (int)(cached_text_sa_pos_a>>40)-(int)(cached_text_sa_pos_b>>40);
  }
  // Full compare
  const bool a_lt_b = (sa_pos_a < sa_pos_b);
  const uint64_t cmp_length = (a_lt_b) ? global_enc_text_length-sa_pos_b : global_enc_text_length-sa_pos_a;
  const int cmp = memcmp(a1,a2,cmp_length);
  if (cmp) return cmp;
  // Circular comparison
  if (a_lt_b) {
    return memcmp(global_enc_text+sa_pos_a+cmp_length,global_enc_text,global_enc_text_length-(sa_pos_a+cmp_length));
  } else {
    return memcmp(global_enc_text,global_enc_text+sa_pos_b+cmp_length,global_enc_text_length-(sa_pos_b+cmp_length));
  }
}
/*
 *  Multikey quicksort (from Bentley-Sedgewick)
 *    Stops when text_depth reaches Shallow_depth_limit
 *    that is when we have found that the current set of strings
 *    have @ds_shallow_limit chars in common
 */
// Auxiliary procedures and macro for Bentley-Sedgewick's multikey quicksort
void sa_builder_vector_swap(uint64_t* a,uint64_t* b,uint64_t length) {
  while (length-- > 0) {
    const uint64_t t = *a;
    *a++ = *b;
    *b++ = t;
  }
}
uint8_t sa_builder_get_char(const uint64_t* const sa_position) {
  return SA_POS_MASK_FIRST_CHARACTER(*sa_position);
}
uint64_t sa_builder_word64(const uint8_t* const text,const uint64_t* const sa_position) {
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  #ifdef __MACH__ // OS X
    return (uint64_t)text[SA_POS_MASK_POSITION(*sa_position)]   << 9 |
           (uint64_t)text[SA_POS_MASK_POSITION(*sa_position)+1] << 6 |
           (uint64_t)text[SA_POS_MASK_POSITION(*sa_position)+2] << 3 |
           (uint64_t)text[SA_POS_MASK_POSITION(*sa_position)+3];
  #else
    //  _int64 _bswap64(__int64 x);
    //    Reverses the byte order of x. Swaps 8 bytes; bits 0-7 are swapped with bits 56-63,
    //    bits 8-15 are swapped with bits 48-55, bits 16-23 are swapped with bits 40-47,
    //    and bits 24-31 are swapped with bits 32-39.
    const uint64_t word64 = *((uint64_t*)(text + SA_POS_MASK_POSITION(*sa_position)));
    return __bswap_64(word64);
  #endif
#else
  const uint64_t word64 = *((uint64_t*)(text + SA_POS_MASK_POSITION(*sa_position)));
  return word64;
#endif
}
uint64_t sa_builder_word64_cached(const uint64_t* const sa_position) {
  return SA_POS_MASK_SUFFIX_CACHE(*sa_position);
}
uint64_t* sa_builder_med3(
    uint64_t* const a,
    uint64_t* const b,
    uint64_t* const c) {
  uint8_t va, vb, vc;
  va = sa_builder_get_char(a);
  vb = sa_builder_get_char(b);
  if (va == vb) return a;
  vc = sa_builder_get_char(c);
  if (vc == va || vc == vb) return c;
  return va < vb ?
        (vb < vc ? b : (va < vc ? c : a ) ) :
        (vb > vc ? b : (va < vc ? a : c ) ) ;
}
// Multikey quicksort (from Bentley-Sedgewick)
void sa_builder_ds_shallow_mkq(
    uint64_t* const a,
    const uint64_t n,
    const uint8_t* text) {
#define swap2(a,b) { t = *(a); *(a) = *(b); *(b) = t; }
  uint64_t partval, val;
  uint64_t *pa, *pb, *pc, *pd, *pl, *pm, *pn;
  uint64_t t, r;
  const uint8_t* next_depth;
//  // Small arrays, use insertions sort
//  if (n < DS_MK_QUICKSORT_THRESHOLD) {
//    shallow_inssort_lcp(a,n,text); // FIXME
//    return;
//  }
  /*
   * Choose Pivot & Swap pointers
   */
  bool repeat_pivot = true;
  do {
    // Choose pivot
    pl = a;
    pm = a + (n/2);
    pn = a + (n-1);
    if (n > 30) { // On big arrays use pseudomedian of 9
      const uint64_t d = (n/8);
      pl = sa_builder_med3(pl,pl+d,pl+2*d);
      pm = sa_builder_med3(pm-d,pm,pm+d);
      pn = sa_builder_med3(pn-2*d,pn-d,pn);
    }
    pm = sa_builder_med3(pl,pm,pn);
    swap2(a,pm);
    partval = sa_builder_word64(text,a);
    pa = pb = a + 1;
    pc = pd = a + n-1;
    // Partition
    for (;;) {
      while (pb <= pc && (val=sa_builder_word64(text,pb)) <= partval) {
        if (val == partval) { swap2(pa,pb); pa++; }
        pb++;
      }
      while (pb <= pc && (val=sa_builder_word64(text,pc)) >= partval) {
        if (val == partval) { swap2(pc,pd); pd--; }
        pc--;
      }
      if (pb > pc) break;
      swap2(pb,pc);
      pb++;
      pc--;
    }
    if (pa > pd) {
      // All values were equal to partval: make it simpler
      next_depth = text+4;
      if (next_depth >= ds_shallow_text_limit) {
        // helped_sort(a,n,next_depth-Text); // FIXME
        qsort(a,n,sizeof(uint64_t),(int (*)(const void *,const void *))sa_builder_suffix_cmp);
        return;
      } else {
        text = next_depth;
      }
    } else {
      repeat_pivot = false;
    }
  } while (repeat_pivot);
  /*
   * Partition & Sort
   *   a[] into the values smaller, equal, and larger that partval
   */
  pn = a + n;
  r = MIN(pa-a,pb-pa);     // min{|SA==|,|SA<|}
  sa_builder_vector_swap(a,pb-r,r);
  r = MIN(pd-pc,pn-pd-1);  // min{|SA>|,|SA==|}
  sa_builder_vector_swap(pb,pn-r,r);
  // Sort smaller strings (SA<)
  if ((r = pb-pa) > 1) {
    sa_builder_ds_shallow_mkq(a,r,text);
  }
  // Sort strings starting with partval
  next_depth = text+4;
  if (next_depth < ds_shallow_text_limit) {
    sa_builder_ds_shallow_mkq(a+r,pa-pd+n-1,next_depth);
  } else {
    // helped_sort(a+r,pa-pd+n-1,next_depth-Text); // FIXME
    qsort(a+r,pa-pd+n-1,sizeof(uint64_t),(int (*)(const void *,const void *))sa_builder_suffix_cmp);
  }
  // Sort greater strings (SA>)
  if ((r = pd-pc) > 1) {
    sa_builder_ds_shallow_mkq(a+n-r,r,text);
  }
}
void sa_builder_ds_shallow_mkq_cached(
    uint64_t* const a,
    const uint64_t n,
    const uint8_t* text) {
  uint64_t partval, val;
  uint64_t *pa, *pb, *pc, *pd, *pl, *pm, *pn;
  uint64_t t, r;
  /*
   * Choose Pivot & Swap pointers
   */
  // Choose pivot
  pl = a;
  pm = a + (n/2);
  pn = a + (n-1);
  if (n > 30) { // On big arrays use pseudomedian of 9
    const uint64_t d = (n/8);
    pl = sa_builder_med3(pl,pl+d,pl+2*d);
    pm = sa_builder_med3(pm-d,pm,pm+d);
    pn = sa_builder_med3(pn-2*d,pn-d,pn);
  }
  pm = sa_builder_med3(pl,pm,pn);
  swap2(a,pm);
  partval = sa_builder_word64_cached(a);
  pa = pb = a + 1;
  pc = pd = a + n-1;
  // Partition
  for (;;) {
    while (pb <= pc && (val=sa_builder_word64_cached(pb)) <= partval) {
      if (val == partval) { swap2(pa,pb); pa++; }
      pb++;
    }
    while (pb <= pc && (val=sa_builder_word64_cached(pc)) >= partval) {
      if (val == partval) { swap2(pc,pd); pd--; }
      pc--;
    }
    if (pb > pc) break;
    swap2(pb,pc);
    pb++;
    pc--;
  }
  /*
   * Partition & Sort
   *   a[] into the values smaller, equal, and larger that partval
   */
  pn = a + n;
  r = MIN(pa-a,pb-pa);     // min{|SA==|,|SA<|}
  sa_builder_vector_swap(a,pb-r,r);
  r = MIN(pd-pc,pn-pd-1);  // min{|SA>|,|SA==|}
  sa_builder_vector_swap(pb,pn-r,r);
  // Sort smaller strings (SA<)
  if ((r = pb-pa) > 1) {
    sa_builder_ds_shallow_mkq_cached(a,r,text);
  }
  // Sort strings starting with partval
  sa_builder_ds_shallow_mkq(a+r,pa-pd+n-1,text+4);
  // Sort greater strings (SA>)
  if ((r = pd-pc) > 1) {
    sa_builder_ds_shallow_mkq_cached(a+n-r,r,text);
  }
}
void* sa_builder_sort_suffixes_thread(uint64_t thread_id) {
  gem_thread_register_id(thread_id+1);
  // SA sampling rate
  const uint64_t sa_sampling_rate_pow2 = global_sampled_sa!=NULL ? global_sampled_sa->sa_sampling_rate : 0;
  const uint64_t text_sampling_rate_pow2 = global_sampled_sa!=NULL ? global_sampled_sa->text_sampling_rate : 0;
  // Retrieve SA chunks
  fm_t* const sa_file_reader = global_sa_builder->sa_file_reader[thread_id];
  vector_t* const buffer = vector_new(global_sa_builder->block_size/8,uint64_t);
  uint64_t group_id;
  for (group_id=0;group_id<global_sa_builder->num_sa_groups;++group_id) {
    if (global_sa_builder->sa_groups[group_id].thread_responsible != thread_id) continue;
    // Sort SA-Block
    const sa_group_t* const sa_group = global_sa_builder->sa_groups+group_id;
    vector_reserve(buffer,sa_group->num_sa_positions,false);
    vector_set_used(buffer,sa_group->num_sa_positions);
    fm_seek(sa_file_reader,sa_group->sa_offset*UINT64_SIZE);
    fm_read_mem(sa_file_reader,vector_get_mem(buffer,uint64_t),sa_group->num_sa_positions*UINT64_SIZE);
    sa_builder_ds_shallow_mkq_cached(vector_get_mem(buffer,uint64_t),vector_get_used(buffer),global_enc_text); // FM's ds_sort
    // Store BWT
    const uint64_t* const sa_chunk = vector_get_mem(buffer,uint64_t);
    uint64_t block_position;
    for (block_position=0;block_position<sa_group->num_sa_positions;++block_position) {
      // Store BWT
      const uint64_t sa_idx = sa_group->sa_offset+block_position;
      const uint8_t sa_char = SA_POS_MASK_GET_BWT1(sa_chunk[block_position]);
      global_enc_bwt[sa_idx] = sa_char;
      // Store SA-samples (hybrid)
      if (global_sampled_sa!=NULL) {
        // Store SA-Space sample (Compacted Samples / CPU)
        const uint64_t text_position = SA_POS_MASK_POSITION(sa_chunk[block_position]);
        // Sampling-SA Direct (Text)
        if (MOD_POW2(text_position,text_sampling_rate_pow2) == 0) {
          // Store SA-sample (Hybrid-SA)
          sampled_sa_builder_set_sample(global_sampled_sa,thread_id,sa_idx,text_position);
        } else {
          // Sampling-SA Inverse (SA)
          if (MOD_POW2(sa_idx,sa_sampling_rate_pow2) == 0) {
            // Store SA-sample (Hybrid-SA)
            sampled_sa_builder_set_sample(global_sampled_sa,thread_id,sa_idx,text_position);
          }
        }
        // Store SA-Space sample (Raw sample / GPU)
        if (global_sampled_sa->sa_raw_samples!=NULL) {
          if (MOD_POW2(sa_idx,sa_sampling_rate_pow2) == 0) {
            const uint64_t sa_space_idx = DIV_POW2(sa_idx,sa_sampling_rate_pow2);
            global_sampled_sa->sa_raw_samples[sa_space_idx] = text_position;
          }
        }
      }
    }
    // Ticker update
    ticker_update_mutex(&global_sa_builder->ticker,1);
  }
  // Free
  vector_delete(buffer);
  // Finish
  return NULL;
}
/*
 * 3.3 Sort all suffixes
 */
void sa_builder_sort_suffixes_prepare_groups(sa_builder_t* const sa_builder) {
  const uint64_t num_threads = sa_builder->num_threads;
  sa_builder->sa_file_reader = mm_calloc(num_threads,fm_t*,false);
  uint64_t i;
  for (i=0;i<num_threads;++i) {
    sa_builder->sa_file_reader[i] = fm_open_file(sa_builder->sa_positions_file_name,FM_READ);
  }
}
//#include "libittnotify.h"
void sa_builder_sort_suffixes(
    sa_builder_t* const sa_builder,
    dna_text_t* const enc_bwt,
    sampled_sa_builder_t* const sampled_sa,
    const bool verbose) {
  /*
   * Prepare
   */
  // Store global information to all threads
  global_sa_builder = sa_builder;
  global_enc_text_length = dna_text_get_length(sa_builder->enc_text);
  global_enc_text = dna_text_get_text(sa_builder->enc_text);
  global_enc_bwt = dna_text_get_text(enc_bwt);
  global_sampled_sa = sampled_sa;
  ds_shallow_text_limit = global_enc_text + DS_SHALLOW_LIMIT;
  // Prepare ticket
  ticker_percentage_reset(&sa_builder->ticker,verbose,"Building-BWT::Sorting SA",sa_builder->num_sa_groups,1,true);
  // Prepare sorting groups
  sa_builder_sort_suffixes_prepare_groups(sa_builder);
  /*
   * Sort suffixes
   */
  // Create & Spawn sorting-threads
  const uint64_t num_threads = sa_builder->num_threads;
  uint64_t i;
  for (i=0;i<num_threads;++i) {
    // Launch thread
    gem_cond_fatal_error(pthread_create(sa_builder->pthreads+i,0,
        (void* (*)(void*))sa_builder_sort_suffixes_thread,(void*)(i)),SYS_THREAD_CREATE);
  }
  for (i=0;i<num_threads;++i) {
    gem_cond_fatal_error(pthread_join(sa_builder->pthreads[i],0),SYS_THREAD_JOIN);
  }
  /*
   * Free
   */
  ticker_finish(&sa_builder->ticker);
  for (i=0;i<num_threads;++i) {
    fm_close(sa_builder->sa_file_reader[i]);
  }
  mm_free(sa_builder->sa_file_reader);
//  __itt_resume();
//  __itt_pause();
}
/*
 * Stats
 */
void sa_builder_record_kmer_count_stats(sa_builder_t* const sa_builder) {
  // Allocate
  uint64_t i, max=0;
  for (i=0;i<sa_builder->num_kmers;++i) {
    const uint64_t kmer_count = sa_builder->kmer_count[i];
    if (kmer_count > max) max = kmer_count;
    stats_vector_inc(sa_builder->kmer_count_stats,kmer_count);
  }
  sa_builder->kmer_count_max = max;
}
void sa_builder_display_stats(
    FILE* const stream,
    sa_builder_t* const sa_builder,
    const bool display_groups) {
  tab_fprintf(stream,"[GEM]>SA.Builder.Stats\n");
  tab_fprintf(stream,"  => Text.Length %"PRIu64"\n",dna_text_get_length(sa_builder->enc_text));
  tab_fprintf(stream,"  => Total.Kmers %"PRIu64"\n",sa_builder->num_kmers);
  tab_fprintf(stream,"    => Kmers.distribution\n");
  tab_global_add(4);
  stats_vector_display(stream,sa_builder->kmer_count_stats,false,true,NULL);
  tab_global_subtract(4);
  // Block Stats
  const uint64_t sa_length = dna_text_get_length(sa_builder->enc_text);
  const uint64_t sa_size = sa_length*UINT64_SIZE;
  const uint64_t preferred_block_size = (sa_length*UINT64_SIZE)/SA_BUILDER_NUM_WRITTERS;
  tab_fprintf(stream,"  => Block.File.Size %"PRIu64" MB\n",CONVERT_B_TO_MB(sa_size));
  tab_fprintf(stream,"    => Block.Max %"PRIu64" MB (%2.3f%%)\n",
      CONVERT_B_TO_MB(sa_builder->max_bucket_size),PERCENTAGE(sa_builder->max_bucket_size,sa_size));
  tab_fprintf(stream,"    => Block.Prefered %"PRIu64" MB for %"PRIu64" writers (%2.3f%%)\n",
      CONVERT_B_TO_MB(preferred_block_size),SA_BUILDER_NUM_WRITTERS,PERCENTAGE(preferred_block_size,sa_size));
  tab_fprintf(stream,"    => Block.Size %"PRIu64" MB\n",
      CONVERT_B_TO_MB(sa_builder->block_size),PERCENTAGE(sa_builder->block_size,sa_size));
  uint64_t i;
  // Groups Stats
  tab_fprintf(stream,"  => Num.Groups %"PRIu64" \n",sa_builder->num_sa_groups);
  if (display_groups) {
    for (i=0;i<sa_builder->num_sa_groups;++i) {
      tab_fprintf(stream,"    => Group[%04lu]\tRange=[%8lu,%8lu)\t%"PRIu64" kmers\n",
          i,sa_builder->sa_groups[i].sa_offset,
          sa_builder->sa_groups[i].sa_offset+sa_builder->sa_groups[i].num_sa_positions,
          sa_builder->sa_groups[i].num_sa_positions);
    }
  } else {
    uint64_t max=0;
    for (i=0;i<sa_builder->num_sa_groups;++i) {
      const uint64_t group_size = sa_builder->sa_groups[i].num_sa_positions;
      if (group_size > max) max = group_size;
      stats_vector_inc(sa_builder->group_size_stats,group_size);
    }
    tab_fprintf(stream,"    => Groups.distribution\n");
    tab_global_add(4);
    stats_vector_display(stream,sa_builder->group_size_stats,false,true,NULL);
    tab_global_subtract(4);
  }
  // Calculate load balancing
  uint64_t* const kmers_per_thread = mm_calloc(sa_builder->num_threads,uint64_t,true);
  uint64_t* const groups_per_thread = mm_calloc(sa_builder->num_threads,uint64_t,true);
  for (i=0;i<sa_builder->num_sa_groups;++i) {
    const uint64_t group_size = sa_builder->sa_groups[i].num_sa_positions;
    kmers_per_thread[sa_builder->sa_groups[i].thread_responsible] += group_size;
    groups_per_thread[sa_builder->sa_groups[i].thread_responsible]++;
  }
  tab_fprintf(stream,"  => Load.Balance (NumKmers / NumBlocks / NumSubBuckets)\n");
  for (i=0;i<sa_builder->num_threads;++i) {
    tab_fprintf(stream,"    => Thread[%"PRIu64"] \t %"PRIu64"(%2.3f%%) \t %"PRIu64"(%2.3f%%)\n",i,
        kmers_per_thread[i],PERCENTAGE(kmers_per_thread[i],sa_length),
        groups_per_thread[i],PERCENTAGE(groups_per_thread[i],sa_builder->num_sa_groups));
  }
  mm_free(kmers_per_thread);
  mm_free(groups_per_thread);
  // Flush
  fflush(stream);
}
/*
 * Debug
 */
void sa_builder_debug_print_sa(
    FILE* stream,
    sa_builder_t* const sa_builder,
    const uint64_t sa_position,
    const uint64_t sa_suffix_length) {
  const uint8_t* const enc_text = dna_text_get_text(sa_builder->enc_text);
  const uint64_t enc_text_length = dna_text_get_length(sa_builder->enc_text);
  const uint64_t suffix_pos = SA_POS_MASK_POSITION(sa_position);
  // fprintf(stream,"Suffix=%011lu\t\t",suffix_pos);
  fprintf(stream,"Suffix=%011"PRIu64"\t%c%c\t",suffix_pos,
      dna_decode(SA_POS_MASK_GET_BWT2(sa_position)),
	  dna_decode(SA_POS_MASK_GET_BWT1(sa_position)));
  // Print begin-suffix
  uint64_t num_printed_chars = 0;
  uint64_t i = suffix_pos;
  while (i < enc_text_length) {
    fprintf(stream,"%c",dna_decode(enc_text[i++]));
    if (++num_printed_chars >= sa_suffix_length) {
      fprintf(stream,"\n");
      return;
    }
  }
  i = 0;
  while (i < suffix_pos) {
    fprintf(stream,"%c",dna_decode(enc_text[i++]));
    if (++num_printed_chars >= sa_suffix_length) {
      fprintf(stream,"\n");
      return;
    }
  }
  fprintf(stream,"\n");
}
