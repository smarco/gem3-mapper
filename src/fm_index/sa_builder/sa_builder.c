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
 *   Implements a data storage structure as to classify SA positions with respect
 *   to it's stating k-mer. SA positions are stored in blocks of a given size
 */

#include "text/dna_text.h"
#include "fm_index/sa_builder/sa_builder.h"
#include "stats/stats_vector.h"

/*
 * Debug
 */
#define SA_BUILDER_DEBUG_SPLIT_BLOCKS           true
#define SA_BUILDER_DUMP_BWT_NUM_BLOCKS         10000

/*
 * Global
 */
uint64_t kmer_count_range[] = {0,10,100,1000,10000,100000};
uint64_t suffix_cmp_length_range[] = {0,1,2,3,4,5,6,7,8,9,10,100,1000,10000};

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
 * Stats
 */
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
