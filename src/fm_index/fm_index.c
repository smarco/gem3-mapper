/*
 * PROJECT: GEMMapper
 * FILE: fm_index.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "fm_index/fm_index.h"
#include "fm_index/rank_mtable.h"
#include "fm_index/rank_mtable_builder.h"

/*
 * FM-Index Model & Version
 */
#define FM_INDEX_MODEL_NO  1006ull

/*
 * Errors
 */
#define GEM_ERROR_FM_INDEX_WRONG_MODEL_NO "FM-index error. Wrong FM-Index Model %"PRIu64" (Expected model %"PRIu64")"
#define GEM_ERROR_FM_INDEX_INDEX_OOR "FM-index error. Index position %"PRIu64" out-of-range BWT[0,%"PRIu64")"


/*
 * Builder
 */
bwt_builder_t* fm_index_write(
    fm_t* const file_manager,
    const bool fm_index_reverse,
    dna_text_t* const bwt_text,
    uint64_t* const character_occurrences,
    sampled_sa_builder_t* const sampled_sa,
    const bool check,
    const bool verbose) {
  // Write Header
  const uint64_t text_length = dna_text_get_length(bwt_text);
  const uint64_t proper_length = log2(text_length)/2;
  fm_write_uint64(file_manager,FM_INDEX_MODEL_NO);
  fm_write_uint64(file_manager,text_length);
  fm_write_uint64(file_manager,proper_length);
  fm_write_uint64(file_manager,fm_index_reverse);
  // Write Sampled-SA & Free Samples
  sampled_sa_builder_write(file_manager,sampled_sa);
  if (verbose) sampled_sa_builder_print(gem_info_get_stream(),sampled_sa);
  sampled_sa_builder_delete_samples(sampled_sa); // Free Samples (Just the samples)
  // Generate BWT-Bitmap & rank_mtable
  bwt_builder_t* const bwt_builder = bwt_builder_new(bwt_text,character_occurrences,sampled_sa,check,verbose);
  if (verbose) bwt_builder_print(gem_info_get_stream(),bwt_builder);
  // Build mrank table
  rank_mtable_t* const rank_mtable = rank_mtable_builder_new(bwt_builder,verbose);
  if (verbose) rank_mtable_print(gem_info_get_stream(),rank_mtable);
  // Write mrank table
  rank_mtable_builder_write(file_manager,rank_mtable);
  rank_mtable_builder_delete(rank_mtable); // Free
  // Write BWT
  bwt_builder_write(file_manager,bwt_builder);
  // Return BWT
  return bwt_builder;
}
bwt_reverse_builder_t* fm_index_reverse_write(
    fm_t* const file_manager,
    dna_text_t* const bwt_reverse_text,
    uint64_t* const character_occurrences,
    const bool check,
    const bool verbose) {
  // Generate BWT-Bitmap Reverse
  bwt_reverse_builder_t* const bwt_reverse_builder =
      bwt_reverse_builder_new(bwt_reverse_text,character_occurrences,check,verbose);
  if (verbose) bwt_reverse_builder_print(gem_info_get_stream(),bwt_reverse_builder);
  // Build mrank table Reverse
  rank_mtable_t* const rank_mtable = rank_mtable_reverse_builder_new(bwt_reverse_builder,verbose);
  if (verbose) rank_mtable_print(gem_info_get_stream(),rank_mtable);
  // Write mrank table Reverse
  rank_mtable_builder_write(file_manager,rank_mtable);
  rank_mtable_builder_delete(rank_mtable); // Free
  // Write BWT Reverse
  bwt_reverse_builder_write(file_manager,bwt_reverse_builder);
  // Return BWT
  return bwt_reverse_builder;
}
/*
 * Loader
 */
fm_index_t* fm_index_read_mem(mm_t* const memory_manager,const bool check) {
  // Allocate handler
  fm_index_t* const fm_index = mm_alloc(fm_index_t);
  // Read Header
  const uint64_t fm_index_model_no = mm_read_uint64(memory_manager);
  gem_cond_fatal_error(fm_index_model_no!=FM_INDEX_MODEL_NO,FM_INDEX_WRONG_MODEL_NO,fm_index_model_no,(uint64_t)FM_INDEX_MODEL_NO);
  fm_index->text_length = mm_read_uint64(memory_manager);
  fm_index->proper_length = mm_read_uint64(memory_manager);
  fm_index->fm_index_reverse = mm_read_uint64(memory_manager);
  // Load Sampled SA
  fm_index->sampled_sa = sampled_sa_read_mem(memory_manager);
  // Load rank_mtable
  fm_index->rank_table = rank_mtable_read_mem(memory_manager);
  // Load BWT
  fm_index->bwt = bwt_read_mem(memory_manager,check);
  // Reverse Index
  if (fm_index->fm_index_reverse) {
    // Load rank_mtable Reverse
    fm_index->rank_table_reverse = rank_mtable_read_mem(memory_manager);
    // Load BWT Reverse
    fm_index->bwt_reverse = bwt_reverse_read_mem(memory_manager,check);
  } else {
    fm_index->rank_table_reverse = NULL;
    fm_index->bwt_reverse = NULL;
  }
  // Return
  return fm_index;
}
void fm_index_delete(fm_index_t* const fm_index) {
  // Delete Sampled SA
  sampled_sa_delete(fm_index->sampled_sa);
  // Delete rank_mtable
  rank_mtable_delete(fm_index->rank_table);
  if (fm_index->fm_index_reverse) {
    rank_mtable_delete(fm_index->rank_table_reverse);
  }
  // Delete BWT
  bwt_delete(fm_index->bwt);
  if (fm_index->fm_index_reverse) {
    bwt_reverse_delete(fm_index->bwt_reverse);
  }
  // Free handler
  mm_free(fm_index);
}

/*
 * Accessors
 */
uint64_t fm_index_get_length(const fm_index_t* const fm_index) {
  return fm_index->text_length;
}
double fm_index_get_proper_length(const fm_index_t* const fm_index) {
  return fm_index->proper_length;
}
uint64_t fm_index_get_size(const fm_index_t* const fm_index) {
  const uint64_t sampled_sa_size = sampled_sa_get_size(fm_index->sampled_sa); // Sampled SuffixArray positions
  const uint64_t bwt_size = bwt_get_size(fm_index->bwt); // BWT structure
  const uint64_t bwt_reverse_size = (fm_index->fm_index_reverse) ?
      bwt_reverse_get_size(fm_index->bwt_reverse) : 0; // BWT Reverse structure
  const uint64_t rank_table_size = rank_mtable_get_size(fm_index->rank_table); // Memoizated intervals
  const uint64_t rank_table_reverse_size = (fm_index->fm_index_reverse) ?
      rank_mtable_get_size(fm_index->rank_table_reverse) : 0; // Memoizated Reverse intervals
  return sampled_sa_size + bwt_size + bwt_reverse_size + rank_table_size + rank_table_reverse_size;
}
/*
 * FM-Index High-level Operators
 */
uint64_t fm_index_decode(const fm_index_t* const fm_index,uint64_t bwt_position) {
  // Parameters
  const bwt_t* const bwt = fm_index->bwt;
  const uint64_t bwt_length = fm_index_get_length(fm_index);
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
uint64_t fm_index_encode(const fm_index_t* const fm_index,const uint64_t text_position) {
  GEM_NOT_IMPLEMENTED(); // TODO Implement
// // Compute SA^(-1)[i]
//  register const idx_t refl=a->text_length-i;
//  register const idx_t sam_quo=GEM_DIV_BY_POWER_OF_TWO_64(refl,a->sampling_rate_log);
//  register idx_t sam_rem=GEM_MOD_BY_POWER_OF_TWO_64(refl,a->sampling_rate_log), pos;
//  GEM_UNALIGNED_ACCESS(pos,a->I,a->cntr_bytes,sam_quo);
//  while (sam_rem--)
//    pos=bwt_LF(a->bwt,pos);
//  return pos;
  return 0;
}
uint64_t fm_index_psi(const fm_index_t* const fm_index,const uint64_t bwt_position) {
  GEM_NOT_IMPLEMENTED(); // TODO Implement
  /* // Compute Psi[i]
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
void fm_index_retrieve_bwt_sampled(
    const fm_index_t* const fm_index,
    uint64_t bwt_position,
    uint64_t* const sampled_bwt_position,
    uint64_t* const lf_dist) {
  // Parameters
  const bwt_t* const bwt = fm_index->bwt;
  bool is_sampled = false;
  *lf_dist=0;
  // Retrieve BWT-pos sampled LF (until we find a sampled position)
  uint64_t bwt_position_next;
  bwt_position_next = bwt_LF(bwt,bwt_position,&is_sampled);
  while (!is_sampled) {
    ++(*lf_dist);
    bwt_position = bwt_position_next;
    bwt_position_next = bwt_LF(bwt,bwt_position,&is_sampled);
  }
  *sampled_bwt_position = bwt_position;
  PROF_ADD_COUNTER(GP_FMIDX_LOOKUP_DIST,*lf_dist);
}
void fm_index_retrieve_sa_sample(
    const fm_index_t* const fm_index,
    const uint64_t sampled_bwt_position,
    const uint64_t lf_dist,
    uint64_t* const text_position) {
  // Parameters
  const uint64_t bwt_length = fm_index_get_length(fm_index);
  const sampled_sa_t* const sampled_sa = fm_index->sampled_sa;
  // Recover sampled position (SA-Sample) & adjust
  const uint64_t sampling_erank = bwt_sampling_erank(fm_index->bwt,sampled_bwt_position);
  *text_position = (sampled_sa_get_sample(sampled_sa,sampling_erank) + lf_dist) % bwt_length;
}
/*
 * Display
 */
void fm_index_print(
    FILE* const stream,
    const fm_index_t* const fm_index) {
  // Calculate some figures
  const uint64_t sampled_sa_size = sampled_sa_get_size(fm_index->sampled_sa); // Sampled SuffixArray positions
  const uint64_t rank_table_size = rank_mtable_get_size(fm_index->rank_table); // Memoizated intervals
  const uint64_t bwt_size = bwt_get_size(fm_index->bwt); // BWT structure
  const uint64_t bwt_reverse_size = bwt_reverse_get_size(fm_index->bwt_reverse); // BWT Reverse structure
  const uint64_t fm_index_size = sampled_sa_size+bwt_size+rank_table_size;
  tab_fprintf(stream,"[GEM]>FM.Index\n");
  tab_fprintf(stream,"  => FM.Index.Size  %"PRIu64" MB (100%%)\n",CONVERT_B_TO_MB(fm_index_size));
  tab_fprintf(stream,"    => Sampled.SA   %"PRIu64" MB (%2.3f%%)\n",CONVERT_B_TO_MB(sampled_sa_size),PERCENTAGE(sampled_sa_size,fm_index_size));
  tab_fprintf(stream,"    => Rank.mTable  %"PRIu64" MB (%2.3f%%)\n",CONVERT_B_TO_MB(rank_table_size),PERCENTAGE(rank_table_size,fm_index_size));
  tab_fprintf(stream,"    => BWT          %"PRIu64" MB (%2.3f%%)\n",CONVERT_B_TO_MB(bwt_size),PERCENTAGE(bwt_size,fm_index_size));
  tab_fprintf(stream,"    => BWT-Reverse  %"PRIu64" MB (%2.3f%%)\n",CONVERT_B_TO_MB(bwt_reverse_size),PERCENTAGE(bwt_reverse_size,fm_index_size));
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
void fm_index_decode_print_benchmark(
    FILE* const stream,
    const fm_index_t* const fm_index,
    uint64_t bwt_position) {
  // Parameters
  const bwt_t* const bwt = fm_index->bwt;
  const sampled_sa_t* const sampled_sa = fm_index->sampled_sa;
  const uint64_t sampling_mod = sampled_sa_get_sa_sampling_rate(sampled_sa);
  // Profile Data
  const uint64_t bwt_position_base = bwt_position;
  uint8_t char_enc;
  uint64_t num_chars[5] = {0,0,0,0,0};
  // LF loop (until we find a sampled position)
  uint64_t dist=0;
  while (bwt_position % sampling_mod != 0) {
    // LF
    bwt_position = bwt_LF_(bwt,bwt_position,&char_enc);
    ++num_chars[char_enc];
    ++dist;
  }
  // Print benchmark line
  fprintf(stream,"%lu\t%lu\t%lu\t|\t%lu\t%lu\t%lu\t%lu\t%lu\n",
      bwt_position_base,bwt_position,dist,
      num_chars[0],num_chars[1],num_chars[2],num_chars[3],num_chars[4]);
  fflush(stream);
}
