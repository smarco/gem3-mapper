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
 * Errors
 */
#define GEM_ERROR_SA_BUILDER_LIMIT_MAX_BLOCK_MEM "SA Builder. Maximum SA-bucket (%"PRIu64" GB) larger than max-memory-block allowed (%"PRIu64" GB)"

/*
 * Global
 */
extern sa_builder_t* global_sa_builder;

/*
 * Store all suffixes
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
      sa_builder_store_sa_pos(group,sa_pos,kmer_idx);
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
  global_sa_builder = sa_builder; // Set global data
  // Launch Writing Threads
  const uint64_t ticker_max = dna_text_get_length(sa_builder->enc_text)+SA_BWT_CYCLIC_LENGTH;
  ticker_percentage_reset(&sa_builder->ticker,
      verbose,"Building-BWT::Generating SA-Positions",ticker_max,1,true);
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
