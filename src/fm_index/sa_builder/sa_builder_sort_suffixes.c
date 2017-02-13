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
 */

#include "text/dna_text.h"
#include "fm_index/sa_builder/sa_builder.h"
#include "stats/stats_vector.h"

/*
 * Global
 */
sa_builder_t* global_sa_builder;
uint64_t global_enc_text_length;
const uint8_t* global_enc_text;
uint8_t* global_enc_bwt;
sampled_sa_builder_t* global_sampled_sa;
const uint8_t* ds_shallow_text_limit;

/*
 * Subsidiary sort
 */
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
int sa_builder_suffix_n_cmp(
    const uint64_t text_position_a,
    const uint64_t text_position_b,
    const uint64_t n) {
  // Eliminate the encoded text (PiggyBack)
  const void* a1 = global_enc_text+text_position_a;
  PREFETCH(a1);
  const void* a2 = global_enc_text+text_position_b;
  PREFETCH(a2);
  // Compute compare length
  const bool a_lt_b = (text_position_a < text_position_b);
  uint64_t cmp_length = (a_lt_b) ? global_enc_text_length-text_position_b : global_enc_text_length-text_position_a;
  if (n < cmp_length) cmp_length = n;
  // N compare
  return memcmp(a1,a2,cmp_length);
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
  const uint64_t sa_sampling_rate_pow2 = global_sampled_sa!=NULL ?
      global_sampled_sa->sa_sampling_rate : SAMPLING_RATE_NONE;
  const uint64_t text_sampling_rate_pow2 = global_sampled_sa!=NULL ?
      global_sampled_sa->text_sampling_rate : SAMPLING_RATE_NONE;
  const bool sa_sampling_enabled = (sa_sampling_rate_pow2 != SAMPLING_RATE_NONE);
  const bool text_sampling_enabled = (text_sampling_rate_pow2 != SAMPLING_RATE_NONE);
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
        if (text_sampling_enabled && MOD_POW2(text_position,text_sampling_rate_pow2) == 0) {
          // Store SA-sample (Hybrid-SA)
          sampled_sa_builder_set_sample(global_sampled_sa,thread_id,sa_idx,text_position);
        } else {
          // Sampling-SA Inverse (SA)
          if (sa_sampling_enabled && MOD_POW2(sa_idx,sa_sampling_rate_pow2) == 0) {
            // Store SA-sample (Hybrid-SA)
            sampled_sa_builder_set_sample(global_sampled_sa,thread_id,sa_idx,text_position);
          }
        }
        // Store SA-Space sample (Raw sample / GPU)
        if (global_sampled_sa->sa_raw_samples!=NULL) {
          if (sa_sampling_enabled && MOD_POW2(sa_idx,sa_sampling_rate_pow2) == 0) {
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
 * Sort all suffixes
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
