/*
 * PROJECT: GEMMapper
 * FILE: fm_index.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "fm_index.h"

/*
 * FM-Index Model & Version
 */
#define FM_INDEX_MODEL_NO  1004ul

/*
 * Builder
 */
GEM_INLINE void fm_index_write(
    fm_t* const file_manager,dna_text_t* const bwt_text,uint64_t* const character_occurrences,
    sampled_sa_builder_t* const sampled_sa,const bool check,const bool verbose) {
  // Write Header
  const uint64_t text_length = dna_text_get_length(bwt_text);
  const uint64_t proper_length = log2(text_length)/2;
  fm_write_uint64(file_manager,FM_INDEX_MODEL_NO);
  fm_write_uint64(file_manager,text_length);
  fm_write_uint64(file_manager,proper_length);
  // Write Sampled-SA & Free Samples
  sampled_sa_builder_write(file_manager,sampled_sa);
  if (verbose) sampled_sa_builder_print(gem_info_get_stream(),sampled_sa);
  sampled_sa_builder_delete_samples(sampled_sa); // Free Samples (Just the samples)
  // Generate BWT-Bitmap & rank_mtable
  bwt_builder_t* const bwt_builder = bwt_builder_new(bwt_text,character_occurrences,sampled_sa,check,verbose);
  if (verbose) bwt_builder_print(gem_info_get_stream(),bwt_builder);
  sampled_sa_builder_delete(sampled_sa); // Free Sampled-SA
  // Build mrank table
  rank_mtable_t* const rank_mtable = rank_mtable_builder_new(bwt_builder,verbose);
  if (verbose) rank_mtable_print(gem_info_get_stream(),rank_mtable);
  // Write mrank table
  rank_mtable_write(file_manager,rank_mtable);
  rank_mtable_builder_delete(rank_mtable); // Free
  // Write BWT
  bwt_builder_write(file_manager,bwt_builder);
  bwt_builder_delete(bwt_builder); // Free
}
GEM_INLINE void fm_index_reverse_write(
    fm_t* const file_manager,dna_text_t* const bwt_reverse_text,
    uint64_t* const character_occurrences,const bool check,const bool verbose) {
  // Generate BWT-Bitmap Reverse
  bwt_reverse_builder_t* const bwt_reverse_builder =
      bwt_reverse_builder_new(bwt_reverse_text,character_occurrences,check,verbose);
  if (verbose) bwt_reverse_builder_print(gem_info_get_stream(),bwt_reverse_builder);
  // Write BWT
  bwt_reverse_builder_write(file_manager,bwt_reverse_builder);
  bwt_reverse_builder_delete(bwt_reverse_builder); // Free
}
/*
 * Loader
 */
GEM_INLINE fm_index_t* fm_index_read_mem(mm_t* const memory_manager,const bool check) {
  // Allocate handler
  fm_index_t* const fm_index = mm_alloc(fm_index_t);
  // Read Header
  const uint64_t fm_index_model_no = mm_read_uint64(memory_manager);
  gem_cond_fatal_error(fm_index_model_no!=FM_INDEX_MODEL_NO,FM_INDEX_WRONG_MODEL_NO,fm_index_model_no,FM_INDEX_MODEL_NO);
  fm_index->text_length = mm_read_uint64(memory_manager);
  fm_index->proper_length = mm_read_uint64(memory_manager);
  // Load Sampled SA
  fm_index->sampled_sa = sampled_sa_read_mem(memory_manager);
  // Load rank_mtable
  fm_index->rank_table = rank_mtable_read_mem(memory_manager);
  // Load BWT
  fm_index->bwt = bwt_read_mem(memory_manager,check);
  fm_index->bwt_reverse = bwt_reverse_read_mem(memory_manager,check);
  // Return
  return fm_index;
}
GEM_INLINE bool fm_index_check(const fm_index_t* const fm_index,const bool verbose) {
  FM_INDEX_CHECK(fm_index);
  GEM_NOT_IMPLEMENTED(); // TODO Implement
  return true;
}
GEM_INLINE void fm_index_delete(fm_index_t* const fm_index) {
  FM_INDEX_CHECK(fm_index);
  // Delete Sampled SA
  sampled_sa_delete(fm_index->sampled_sa);
  // Delete rank_mtable
  rank_mtable_delete(fm_index->rank_table);
  // Delete BWT
  bwt_delete(fm_index->bwt);
  bwt_reverse_delete(fm_index->bwt_reverse);
  // Free handler
  mm_free(fm_index);
}

/*
 * Accessors
 */
GEM_INLINE uint64_t fm_index_get_length(const fm_index_t* const fm_index) {
  FM_INDEX_CHECK(fm_index);
  return fm_index->text_length;
}
GEM_INLINE double fm_index_get_proper_length(const fm_index_t* const fm_index) {
  FM_INDEX_CHECK(fm_index);
  return fm_index->proper_length;
}
GEM_INLINE uint64_t fm_index_get_size(const fm_index_t* const fm_index) {
  const uint64_t sampled_sa_size = sampled_sa_get_size(fm_index->sampled_sa); // Sampled SuffixArray positions
  const uint64_t bwt_size = bwt_get_size(fm_index->bwt); // BWT structure
  const uint64_t bwt_reverse_size = bwt_reverse_get_size(fm_index->bwt_reverse); // BWT Reverse structure
  const uint64_t rank_table_size = rank_mtable_get_size(fm_index->rank_table); // Memoizated intervals
  return sampled_sa_size + bwt_size + bwt_reverse_size + rank_table_size;
}
/*
 * FM-Index Bidirectional Operators
 */
GEM_INLINE void fm_index_precompute_2erank_forward(
    const fm_index_t* const fm_index,fm_2erank_elms_t* const fm_2erank_elms,const uint64_t position) {
  // Compute all exclusive ranks
  bwt_block_locator_t block_loc;
  bwt_block_elms_t block_elms;
  bwt_reverse_t* const bwt_reverse = fm_index->bwt_reverse;
  bwt_reverse_precompute(bwt_reverse,position,&block_loc,&block_elms);
  // Store eranks
  fm_2erank_elms->pranks[ENC_DNA_CHAR_A] = bwt_reverse_precomputed_erank(bwt_reverse,ENC_DNA_CHAR_A,&block_loc,&block_elms);
  fm_2erank_elms->pranks[ENC_DNA_CHAR_C] = bwt_reverse_precomputed_erank(bwt_reverse,ENC_DNA_CHAR_C,&block_loc,&block_elms);
  fm_2erank_elms->pranks[ENC_DNA_CHAR_G] = bwt_reverse_precomputed_erank(bwt_reverse,ENC_DNA_CHAR_G,&block_loc,&block_elms);
  fm_2erank_elms->pranks[ENC_DNA_CHAR_T] = bwt_reverse_precomputed_erank(bwt_reverse,ENC_DNA_CHAR_T,&block_loc,&block_elms);
  fm_2erank_elms->pranks[ENC_DNA_CHAR_N] = bwt_reverse_precomputed_erank(bwt_reverse,ENC_DNA_CHAR_N,&block_loc,&block_elms);
  fm_2erank_elms->pranks[ENC_DNA_CHAR_SEP] = bwt_reverse_precomputed_erank(bwt_reverse,ENC_DNA_CHAR_SEP,&block_loc,&block_elms);
}
GEM_INLINE void fm_index_precompute_2erank_backward(
    const fm_index_t* const fm_index,fm_2erank_elms_t* const fm_2erank_elms,const uint64_t position) {
  // Compute all exclusive ranks
  bwt_block_locator_t block_loc;
  bwt_block_elms_t block_elms;
  bwt_t* const bwt = fm_index->bwt;
  bwt_precompute(bwt,position,&block_loc,&block_elms);
  // Store eranks
  fm_2erank_elms->pranks[ENC_DNA_CHAR_A] = bwt_precomputed_erank(bwt,ENC_DNA_CHAR_A,&block_loc,&block_elms);
  fm_2erank_elms->pranks[ENC_DNA_CHAR_C] = bwt_precomputed_erank(bwt,ENC_DNA_CHAR_C,&block_loc,&block_elms);
  fm_2erank_elms->pranks[ENC_DNA_CHAR_G] = bwt_precomputed_erank(bwt,ENC_DNA_CHAR_G,&block_loc,&block_elms);
  fm_2erank_elms->pranks[ENC_DNA_CHAR_T] = bwt_precomputed_erank(bwt,ENC_DNA_CHAR_T,&block_loc,&block_elms);
  fm_2erank_elms->pranks[ENC_DNA_CHAR_N] = bwt_precomputed_erank(bwt,ENC_DNA_CHAR_N,&block_loc,&block_elms);
  fm_2erank_elms->pranks[ENC_DNA_CHAR_SEP] = bwt_precomputed_erank(bwt,ENC_DNA_CHAR_SEP,&block_loc,&block_elms);
}
/*
 * Compute the lo-delta/hi-delta difference as to maintain coherence between forward & backward search interval
 * Eg. char_enc == C
 *    +---+ -\
 *    | A |  | > lo-delta
 *    +---+ -/
 *    | C |
 *    +---+ -\
 *    | G |  |
 *    +---+  |
 *    | T |  | > hi-delta
 *    +---+  |
 *    | N |  |
 *    +---+ -/
 */
GEM_INLINE uint64_t fm_2erank_elms_compute_reverse_erank_lo_delta(
    fm_2erank_elms_t* const lo_2erank_elms,fm_2erank_elms_t* const hi_2erank_elms,const uint8_t char_enc) {
  // Switch character
  uint64_t acc = 0;
  switch (char_enc) {
    case ENC_DNA_CHAR_SEP:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_N]-lo_2erank_elms->pranks[ENC_DNA_CHAR_N];
      // no break
    case ENC_DNA_CHAR_N:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_T]-lo_2erank_elms->pranks[ENC_DNA_CHAR_T];
      // no break
    case ENC_DNA_CHAR_T:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_G]-lo_2erank_elms->pranks[ENC_DNA_CHAR_G];
      // no break
    case ENC_DNA_CHAR_G:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_C]-lo_2erank_elms->pranks[ENC_DNA_CHAR_C];
      // no break
    case ENC_DNA_CHAR_C:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_A]-lo_2erank_elms->pranks[ENC_DNA_CHAR_A];
      // no break
    case ENC_DNA_CHAR_A:
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Return
  return acc;
}
GEM_INLINE uint64_t fm_2erank_elms_compute_reverse_erank_hi_delta(
    fm_2erank_elms_t* const lo_2erank_elms,fm_2erank_elms_t* const hi_2erank_elms,const uint8_t char_enc) {
  // Switch character
  uint64_t acc = 0;
  switch (char_enc) {
    case ENC_DNA_CHAR_A:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_C]-lo_2erank_elms->pranks[ENC_DNA_CHAR_C];
      // no break
    case ENC_DNA_CHAR_C:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_G]-lo_2erank_elms->pranks[ENC_DNA_CHAR_G];
      // no break
    case ENC_DNA_CHAR_G:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_T]-lo_2erank_elms->pranks[ENC_DNA_CHAR_T];
      // no break
    case ENC_DNA_CHAR_T:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_N]-lo_2erank_elms->pranks[ENC_DNA_CHAR_N];
      // no break
    case ENC_DNA_CHAR_N:
      acc += hi_2erank_elms->pranks[ENC_DNA_CHAR_SEP]-lo_2erank_elms->pranks[ENC_DNA_CHAR_SEP];
      // no break
    case ENC_DNA_CHAR_SEP:
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
  // Return
  return acc;
}
GEM_INLINE void fm_index_precomputed_2query_forward(
    fm_2interval_t* const fm_2interval,fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms,const uint8_t char_enc) {
  // Assign forward values (eranks computed forward)
  fm_2interval->forward_lo = lo_2erank_elms->pranks[char_enc];
  fm_2interval->forward_hi = hi_2erank_elms->pranks[char_enc];
  // Assign backward values (eranks computed forward)
  fm_2interval->backward_lo += fm_2erank_elms_compute_reverse_erank_lo_delta(lo_2erank_elms,hi_2erank_elms,char_enc);
  fm_2interval->backward_hi -= fm_2erank_elms_compute_reverse_erank_hi_delta(lo_2erank_elms,hi_2erank_elms,char_enc);
}
GEM_INLINE void fm_index_precomputed_2query_backward(
    fm_2interval_t* const fm_2interval,fm_2erank_elms_t* const lo_2erank_elms,
    fm_2erank_elms_t* const hi_2erank_elms,const uint8_t char_enc) {
  // Assign backward values (eranks computed backward)
  fm_2interval->backward_lo = lo_2erank_elms->pranks[char_enc];
  fm_2interval->backward_hi = hi_2erank_elms->pranks[char_enc];
  // Assign forward values (eranks computed backward)
  fm_2interval->forward_lo += fm_2erank_elms_compute_reverse_erank_lo_delta(lo_2erank_elms,hi_2erank_elms,char_enc);
  fm_2interval->forward_hi -= fm_2erank_elms_compute_reverse_erank_hi_delta(lo_2erank_elms,hi_2erank_elms,char_enc);
}
GEM_INLINE void fm_index_2query_forward(
    const fm_index_t* const fm_index,fm_2interval_t* const fm_2interval,const uint8_t char_enc) {
  // Precompute erank
  fm_2erank_elms_t lo_2erank_elms, hi_2erank_elms;
  fm_index_precompute_2erank_forward(fm_index,&lo_2erank_elms,fm_2interval->forward_lo);
  fm_index_precompute_2erank_forward(fm_index,&hi_2erank_elms,fm_2interval->forward_hi);
  // Compute next 2interval
  fm_index_precomputed_2query_forward(fm_2interval,&lo_2erank_elms,&hi_2erank_elms,char_enc);
}
GEM_INLINE void fm_index_2query_backward(
    const fm_index_t* const fm_index,fm_2interval_t* const fm_2interval,const uint8_t char_enc) {
  // Precompute erank
  fm_2erank_elms_t lo_2erank_elms, hi_2erank_elms;
  fm_index_precompute_2erank_backward(fm_index,&lo_2erank_elms,fm_2interval->forward_lo);
  fm_index_precompute_2erank_backward(fm_index,&hi_2erank_elms,fm_2interval->forward_hi);
  // Compute next 2interval
  fm_index_precomputed_2query_backward(fm_2interval,&lo_2erank_elms,&hi_2erank_elms,char_enc);
}
/*
 * FM-Index High-level Operators
 */
// Compute SA[i]
GEM_INLINE uint64_t fm_index_lookup(const fm_index_t* const fm_index,uint64_t bwt_position) {
  const uint64_t bwt_length = fm_index_get_length(fm_index);
  gem_fatal_check(bwt_position>=bwt_length,FM_INDEX_INDEX_OOR,bwt_position,bwt_length);
  const bwt_t* const bwt = fm_index->bwt;
  const sampled_sa_t* const sampled_sa = fm_index->sampled_sa;
  bool is_sampled = false;
  uint64_t dist=0;
  // LF until we find a sampled position
  bwt_position = bwt_LF(bwt,bwt_position,&is_sampled);
  while (!is_sampled) {
    ++dist;
    bwt_position = bwt_LF(bwt,bwt_position,&is_sampled);
  }
  PROF_ADD_COUNTER(GP_FMIDX_LOOKUP_DIST,dist);
  // Recover sampled position & adjust
  return (sampled_sa_get_sample(sampled_sa,bwt_position) + dist) % bwt_length;
}
// Compute SA^(-1)[i]
GEM_INLINE uint64_t fm_index_inverse_lookup(const fm_index_t* const fm_index,const uint64_t text_position) {
  GEM_NOT_IMPLEMENTED(); // TODO Implement
//#ifndef SUPPRESS_CHECKS
//  gem_cond_fatal_error(i<0||i>a->text_length,BWT_INDEX,i);
//#endif
//  register const idx_t refl=a->text_length-i;
//  register const idx_t sam_quo=GEM_DIV_BY_POWER_OF_TWO_64(refl,a->sampling_rate_log);
//  register idx_t sam_rem=GEM_MOD_BY_POWER_OF_TWO_64(refl,a->sampling_rate_log), pos;
//  GEM_UNALIGNED_ACCESS(pos,a->I,a->cntr_bytes,sam_quo);
//  while (sam_rem--)
//    pos=bwt_LF(a->bwt,pos);
//  return pos;
  return 0;
}
// Compute Psi[i]
GEM_INLINE uint64_t fm_index_psi(const fm_index_t* const fm_index,const uint64_t bwt_position) {
  GEM_NOT_IMPLEMENTED(); // TODO Implement
  /*
  if (!i)
    return a->last;
  --i;
  register int c=0;
  while (a->C[c+1]<=i)
    ++c;
  i-=a->C[c];
  register idx_t l=0,m,res=a->bwt->n;
  while (l<res) {
    m=(l+res)/2;
    if (i<fm_rankc(a,c,m))
      res=m;
    else
      l=m+1;
  }
  return res;
  */
  // return fmi_inverse(a,(fmi_lookup(a,i)+1)%a->bwt->n);
  return 0;
}
// Decode fm_index->text[bwt_position..bwt_position+length-1] into @buffer.
GEM_INLINE uint64_t fm_index_decode(
    const fm_index_t* const fm_index,const uint64_t bwt_position,const uint64_t length,char* const buffer) {
  GEM_NOT_IMPLEMENTED(); // TODO Implement
//  if (__builtin_expect(a->bed!=0,1)) {
//    return be_decode(a->bed,i,len,p);
//  } else {
//    register idx_t pos=fmi_inverse(a,i+len),j;
//    for (j=0;j<len;++j) {
//      register const ch_t c=bwt_char(a->bwt,pos);
//      p[len-1-j]=c;
//      if (__builtin_expect(c==CHAR_ENC_SEP||c==CHAR_ENC_EOT,false))
//        p[len-1-j]=0;
//      pos=bwt_LF(a->bwt,pos);
//    }
//    p[len]=0;
//    return len;
//  }
//
  return 0;
}
GEM_INLINE uint64_t fm_index_decode_raw(
    const fm_index_t* const fm_index,const uint64_t bwt_position,const uint64_t length,char* const buffer) {
  GEM_NOT_IMPLEMENTED(); // TODO Implement
//#ifndef SUPPRESS_CHECKS
//  gem_cond_fatal_error(i<0||i>a->text_length,BWT_INDEX,i);
//  gem_cond_fatal_error(len<0||len>a->text_length,BWT_LEN,len);
//  gem_cond_fatal_error(i+len>a->text_length,PARAM_INV);
//#endif
//  if (__builtin_expect(a->bed!=0,1)) {
//    return be_decode_raw(a->bed,i,len,p);
//  } else {
//    register idx_t pos=fmi_inverse(a,i+len),j;
//    for (j=0;j<len;++j) {
//      p[len-1-j] = bwt_char(a->bwt,pos);
//      pos=bwt_LF(a->bwt,pos);
//    }
//    p[len]=0;
//    return len;
//  }
  return 0;
}
// Basic backward search
GEM_INLINE void fm_index_bsearch_pure(
    const fm_index_t* const fm_index,const uint8_t* const key,
    uint64_t key_length,uint64_t* const hi_out,uint64_t* const lo_out) {
  FM_INDEX_CHECK(fm_index);
  GEM_CHECK_NULL(key);
  GEM_CHECK_POSITIVE(key_length);
  GEM_CHECK_NULL(hi_out); GEM_CHECK_NULL(lo_out);
  // Query lookup table
  uint64_t lo=0, hi=fm_index_get_length(fm_index);
  // Continue with ranks against the FM-Index
//  uint64_t j; for (j=0;j<key_length;++j) printf("%c",dna_decode(key[j])); printf("\n");;
  while (key_length > 0 && hi > lo) {
    const uint8_t c = key[--key_length];
    lo = bwt_erank(fm_index->bwt,c,lo);
    hi = bwt_erank(fm_index->bwt,c,hi);
//    printf("%lu  %lu\n",lo,hi);
  }
//  printf("\n");
  // Return results
  *hi_out=hi;
  *lo_out=lo;
}
GEM_INLINE void fm_index_reverse_bsearch_pure(
    const fm_index_t* const fm_index,const uint8_t* const key,
    const uint64_t key_length,uint64_t* const hi_out,uint64_t* const lo_out) {
  FM_INDEX_CHECK(fm_index);
  GEM_CHECK_NULL(key);
  GEM_CHECK_POSITIVE(key_length);
  GEM_CHECK_NULL(hi_out); GEM_CHECK_NULL(lo_out);
  // Init
  fm_2interval_t fm_2interval;
  fm_2interval.backward_lo = 0;
  fm_2interval.backward_hi = fm_index_get_length(fm_index);
  fm_2interval.forward_lo = 0;
  fm_2interval.forward_hi = fm_index_get_length(fm_index);
  // Search using 2interval and eranks against the FM-Index
//  uint64_t j; for (j=0;j<key_length;++j) printf("%c",dna_decode(key[j])); printf("\n");
  uint64_t i = 0;
  while (i < key_length && fm_2interval.forward_lo < fm_2interval.forward_hi) {
    fm_index_2query_forward(fm_index,&fm_2interval,key[i++]);
//    printf("%lu  %lu\n",fm_2interval.backward_lo,fm_2interval.backward_hi);
  }
//  printf("\n");
  // Return results
  *lo_out = fm_2interval.backward_lo;
  *hi_out = fm_2interval.backward_hi;
}
GEM_INLINE void fm_index_bsearch_debug(
    const fm_index_t* const fm_index,
    const uint8_t* const key,uint64_t key_length,
    uint64_t* const hi_out,uint64_t* const lo_out) {
  FM_INDEX_CHECK(fm_index);
  GEM_CHECK_NULL(key);
  GEM_CHECK_POSITIVE(key_length);
  GEM_CHECK_NULL(hi_out); GEM_CHECK_NULL(lo_out);
  // Query lookup table
  uint64_t lo=0, hi=fm_index_get_length(fm_index);
  // Continue with ranks against the FM-Index
//  bwt_block_locator_t bwt_block_locator_lo, bwt_block_locator_hi;
//  bwt_block_elms_t bwt_block_elms_lo, bwt_block_elms_hi;
  while (key_length > 0 && hi > lo) {
    const uint8_t c = key[--key_length];
//    if (bwt_is_same_bucket(lo,hi)) {
//      bwt_precompute_interval(fm_index->bwt,lo,hi,&bwt_block_locator_lo,&bwt_block_elms_lo);
//      bwt_precomputed_erank_interval(fm_index->bwt,c,&lo,&hi,&bwt_block_locator_lo,&bwt_block_elms_lo);
//    } else {
//      bwt_precompute(fm_index->bwt,lo,&bwt_block_locator_lo,&bwt_block_elms_lo);
//      bwt_precompute(fm_index->bwt,hi,&bwt_block_locator_hi,&bwt_block_elms_hi);
//      lo = bwt_precomputed_erank(fm_index->bwt,c,&bwt_block_locator_lo,&bwt_block_elms_lo);
//      hi = bwt_precomputed_erank(fm_index->bwt,c,&bwt_block_locator_hi,&bwt_block_elms_hi);
//    }
    lo = bwt_erank(fm_index->bwt,c,lo);
    hi = bwt_erank(fm_index->bwt,c,hi);
    printf("> %lu\t%lu\n",lo,hi);
  }
  // Return results
  *hi_out=hi;
  *lo_out=lo;
}
GEM_INLINE void fm_index_bsearch(
    const fm_index_t* const fm_index,const uint8_t* const key,
    uint64_t key_length,uint64_t* const hi_out,uint64_t* const lo_out) {
  FM_INDEX_CHECK(fm_index);
  GEM_CHECK_NULL(key);
  GEM_CHECK_POSITIVE(key_length);
  GEM_CHECK_NULL(hi_out); GEM_CHECK_NULL(lo_out);
  uint64_t lo, hi;
  // Rank queries against the lookup table
  const rank_mtable_t* const rank_mtable = fm_index->rank_table;
  rank_mquery_t query;
  rank_mquery_new(&query);
  //rank_mtable_fetch(rank_mtable,&query,&lo,&hi); printf("> %lu\t%lu\n",lo,hi);
  while (key_length > 0 && !rank_mquery_is_exhausted(&query)) {
    const uint8_t c = key[key_length-1];
    if (c >= ENC_DNA_CHAR_N) break;
    rank_mquery_add_char(rank_mtable,&query,c); // Rank query (calculate offsets)
    --key_length;
    //rank_mtable_fetch(rank_mtable,&query,&lo,&hi); printf("> %lu\t%lu\n",lo,hi);
  }
  // Query lookup table
  rank_mtable_fetch(rank_mtable,&query,&lo,&hi);
  // Check zero interval
  if (hi==lo || key_length==0) {
    *hi_out=hi;
    *lo_out=lo;
  } else {
    // Continue with ranks against the FM-Index
    while (key_length > 0 && hi > lo) {
      const uint8_t c = key[--key_length];
      lo = bwt_erank(fm_index->bwt,c,lo);
      hi = bwt_erank(fm_index->bwt,c,hi);
      //printf("> %lu\t%lu\n",lo,hi);
    }
    // Return results
    *hi_out=hi;
    *lo_out=lo;
  }
}
GEM_INLINE uint64_t fm_index_bsearch_continue(
    const fm_index_t* const fm_index,const char* const key,const uint64_t key_length,
    const bool* const allowed_repl,uint64_t last_hi,uint64_t last_lo,uint64_t begin_pos,
    const uint64_t end_pos,uint64_t* const res_hi,uint64_t* const res_lo) {
  // Extends the exact search of a key from a given position and interval
  const bwt_t* const bwt = fm_index->bwt;
  // Continue search
  while (begin_pos > end_pos) {
    if (last_lo==last_hi) {
      *res_lo=last_lo;
      *res_hi=last_hi;
      return end_pos;
    }
    // Query step
    const uint8_t enc_char = key[--begin_pos];
    if (!allowed_repl[enc_char]) {
      ++begin_pos; break;
    }
    last_lo=bwt_erank(bwt,enc_char,last_lo);
    last_hi=bwt_erank(bwt,enc_char,last_hi);
  }
  *res_lo=last_lo;
  *res_hi=last_hi;
  return begin_pos;
}
/*
 * Display
 */
GEM_INLINE void fm_index_print(FILE* const stream,const fm_index_t* const fm_index) {
  // Calculate some figures
  const uint64_t sampled_sa_size = sampled_sa_get_size(fm_index->sampled_sa); // Sampled SuffixArray positions
  const uint64_t rank_table_size = rank_mtable_get_size(fm_index->rank_table); // Memoizated intervals
  const uint64_t bwt_size = bwt_get_size(fm_index->bwt); // BWT structure
  const uint64_t bwt_reverse_size = bwt_reverse_get_size(fm_index->bwt_reverse); // BWT Reverse structure
  const uint64_t fm_index_size = sampled_sa_size+bwt_size+rank_table_size;
  tab_fprintf(stream,"[GEM]>FM.Index\n");
  tab_fprintf(stream,"  => FM.Index.Size  %lu MB (100%%)\n",CONVERT_B_TO_MB(fm_index_size));
  tab_fprintf(stream,"    => Sampled.SA   %lu MB (%2.3f%%)\n",CONVERT_B_TO_MB(sampled_sa_size),PERCENTAGE(sampled_sa_size,fm_index_size));
  tab_fprintf(stream,"    => Rank.mTable  %lu MB (%2.3f%%)\n",CONVERT_B_TO_MB(rank_table_size),PERCENTAGE(rank_table_size,fm_index_size));
  tab_fprintf(stream,"    => BWT          %lu MB (%2.3f%%)\n",CONVERT_B_TO_MB(bwt_size),PERCENTAGE(bwt_size,fm_index_size));
  tab_fprintf(stream,"    => BWT-Reverse  %lu MB (%2.3f%%)\n",CONVERT_B_TO_MB(bwt_reverse_size),PERCENTAGE(bwt_reverse_size,fm_index_size));
  tab_global_inc();
  // Sampled SuffixArray positions
  sampled_sa_print(stream,fm_index->sampled_sa,false);
  // Memoizated intervals
  rank_mtable_print(stream,fm_index->rank_table);
  // BWT structure
  bwt_print(stream,fm_index->bwt);
  // BWT-Reverse structure
  bwt_reverse_print(stream,fm_index->bwt_reverse);
  tab_global_dec();
  // Flush
  fflush(stream);
}

