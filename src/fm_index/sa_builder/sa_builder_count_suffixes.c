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
#define GEM_ERROR_SA_BUILDER_SEQUENCE_MIN_LENGTH "SA Builder. Index total length (%"PRIu64") is below minimum threshold (%"PRIu64")"

/*
 * Prototypes
 */
void sa_builder_record_kmer_count_stats(sa_builder_t* const sa_builder);

/*
 * Count all suffixes
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
