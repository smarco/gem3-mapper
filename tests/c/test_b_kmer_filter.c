/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2026 by Jacob Schreiber  <jacob.schreiber@umassmed.edu>
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
 * AUTHOR(S): Jacob Schreiber <jacob.schreiber@umassmed.edu> (unit tests, written with Claude Code)
 * DESCRIPTION:
 *   The k-mer counting filter (q-gram lemma, key cut into n tiles): its
 *   bound on the edit distance against hand-worked values and a brute-force
 *   distance. Each tile of t bases holds t-k+1 k-mers and one edit destroys
 *   at most k of them, so a tile whose best text window shares c k-mers
 *   with it needs at least floor((t-k+1-c)/k) edits; the bound sums tiles.
 *   gem-mapper discards a candidate whose bound exceeds the error budget,
 *   so the bound must never exceed a distance within that budget (gem's own
 *   debug check in filtering_region_verify.c says the same).
 *
 *   Internal signatures called (a speed iteration may change them; update
 *   this driver with them):
 *     kmer_counting_compile_nway(kmer_counting_nway_t*,uint8_t* key,uint64_t key_length,
 *         uint64_t num_tiles,uint64_t kmer_length,bool count_pattern_kmers,mm_allocator_t*)
 *     kmer_counting_min_bound_nway(kmer_counting_nway_t*,const uint8_t* text,uint64_t text_length,
 *         uint64_t max_error,mm_allocator_t*)
 *     kmer_counting_destroy(kmer_counting_nway_t*,mm_allocator_t*)
 *     pattern_tiled_compile(pattern_tiled_t*,uint8_t* key,uint64_t key_length,uint64_t max_error,
 *         uint64_t kmer_tiles,uint64_t kmer_length,uint64_t kmer_enabled,mm_allocator_t*)
 *     alignment_verify_edit_kmer(alignment_t*,pattern_tiled_t*,uint8_t* key,uint64_t key_length,
 *         uint8_t* text,uint64_t text_length,uint64_t max_error,uint64_t kmer_tiles,uint64_t kmer_length)
 */

#include "gem3_test.h"
#include "text/dna_text.h"
#include "align/alignment.h"
#include "align/align_kmer_filter_nway.h"
#include "align/pattern/pattern_tiled.h"

#define MAX_LEN 1200

static uint64_t naive_distance(
    const uint8_t* const key,const uint64_t m,
    const uint8_t* const text,const uint64_t n) {
  static uint64_t prev[MAX_LEN+1], curr[MAX_LEN+1];
  uint64_t best = UINT64_MAX, i, j;
  for (i=0;i<=m;++i) prev[i] = i;
  for (j=1;j<=n;++j) {
    curr[0] = 0;
    for (i=1;i<=m;++i) {
      const bool match = key[i-1]==text[j-1] && key[i-1]!=ENC_DNA_CHAR_N;
      uint64_t v = prev[i-1] + (match ? 0 : 1);
      if (prev[i]+1 < v) v = prev[i]+1;
      if (curr[i-1]+1 < v) v = curr[i-1]+1;
      curr[i] = v;
    }
    if (curr[m] < best) best = curr[m];
    memcpy(prev,curr,(m+1)*sizeof(uint64_t));
  }
  return best;
}

static void random_seq(gt_rng_t* const rng,uint8_t* const out,const uint64_t n,const uint64_t n_per_mille) {
  uint64_t i;
  for (i=0;i<n;++i) {
    out[i] = gt_rng_chance(rng,n_per_mille,1000) ? ENC_DNA_CHAR_N : (uint8_t)gt_rng_below(rng,4);
  }
}

static uint64_t embed_with_edits(
    gt_rng_t* const rng,const uint8_t* const key,const uint64_t m,const uint64_t edits,
    const uint64_t left,const uint64_t right,uint8_t* const text) {
  static uint8_t body[MAX_LEN];
  uint64_t b = 0, i, e;
  for (i=0;i<m;++i) body[b++] = key[i];
  for (e=0;e<edits && b>1;++e) {
    const uint64_t pos = gt_rng_below(rng,b);
    switch (gt_rng_below(rng,3)) {
      case 0: body[pos] = (uint8_t)((body[pos]+1+gt_rng_below(rng,3))%4); break;
      case 1:
        if (b+1+left+right < MAX_LEN) {
          memmove(body+pos+1,body+pos,b-pos); body[pos] = (uint8_t)gt_rng_below(rng,4); ++b;
        }
        break;
      default: memmove(body+pos,body+pos+1,b-pos-1); --b; break;
    }
  }
  random_seq(rng,text,left,0);
  memcpy(text+left,body,b);
  random_seq(rng,text+left+b,right,0);
  return left+b+right;
}

static mm_allocator_t* new_allocator(mm_slab_t** const slab) {
  *slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM);
  return mm_allocator_new(*slab);
}

static uint64_t bound_of(
    uint8_t* const key,const uint64_t m,const uint8_t* const text,const uint64_t n,
    const uint64_t tiles,const uint64_t k,const uint64_t max_error,mm_allocator_t* const mm_allocator) {
  kmer_counting_nway_t kc;
  memset(&kc,0,sizeof(kc));
  kmer_counting_compile_nway(&kc,key,m,tiles,k,true,mm_allocator);
  const uint64_t bound = kmer_counting_min_bound_nway(&kc,text,n,max_error,mm_allocator);
  kmer_counting_destroy(&kc,mm_allocator);
  return bound;
}

/*
 * A key with no k-mer in common with the text: every tile contributes
 * floor((t-k+1)/k). An exact copy: 0.
 */
GT_TEST(kmer_bound_hand_examples) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  // Key: (ACGT)x25, 100 bp, no k-mer of a single base; text: 140 A
  static uint8_t key[100], polyA[140], copy[140];
  uint64_t i;
  for (i=0;i<100;++i) key[i] = (uint8_t)(i % 4);
  for (i=0;i<140;++i) polyA[i] = ENC_DNA_CHAR_A;
  static const struct { uint64_t tiles, k, expected; } cases[] = {
    { 1, 6, 15 },          // 95/6
    { 2, 6, 14 },          // tiles 50,50: 45/6 + 45/6
    { 3, 4, 21 },          // tiles 34,34,32: 31/4 + 31/4 + 29/4
    { 4, 5, 16 },          // tiles 25 x4: 21/5 x4
    { 5, 3, 30 },          // tiles 20 x5: 18/3 x5
    { 6, 3, 29 },          // tiles 17 x5 + 15: 15/3 x5 + 13/3
    { 2, 7, 12 },          // 44/7 + 44/7
  };
  uint64_t c;
  for (c=0;c<sizeof(cases)/sizeof(cases[0]);++c) {
    const uint64_t b = bound_of(key,100,polyA,140,cases[c].tiles,cases[c].k,50,mm_allocator);
    GT_CHECK_MSG(b==cases[c].expected,"tiles=%" PRIu64 " k=%" PRIu64 ": bound %" PRIu64 " expected %" PRIu64,
        cases[c].tiles,cases[c].k,b,cases[c].expected);
    // The key itself inside flanks of A
    for (i=0;i<20;++i) copy[i] = ENC_DNA_CHAR_A;
    memcpy(copy+20,key,100);
    for (i=120;i<140;++i) copy[i] = ENC_DNA_CHAR_A;
    GT_CHECK_EQ_UINT(bound_of(key,100,copy,140,cases[c].tiles,cases[c].k,20,mm_allocator),0);
  }
  // A text N counts as A (enc % 4): an all-N text behaves like poly-A
  static uint8_t polyN[140];
  for (i=0;i<140;++i) polyN[i] = ENC_DNA_CHAR_N;
  GT_CHECK_EQ_UINT(bound_of(key,100,polyN,140,2,6,50,mm_allocator),14);
  // Tiles shorter than k: the filter is disabled and never discards
  GT_CHECK_EQ_UINT(bound_of(key,10,polyA,40,4,6,2,mm_allocator),0);
  // k outside 3..7 falls back to 7 (gem warns once per pattern)
  GT_CHECK_EQ_UINT(bound_of(key,100,polyA,140,2,9,50,mm_allocator),12);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * Random keys, edited copies and 1-6 tiles (the 1-, 2-, 3-, 4- and n-tile
 * kernels), k = 3..7: whenever the true distance is within max_error the
 * bound is at most that distance. Candidate texts are at least as long as
 * the key minus 2*max_error, as gem-mapper's are (shorter ones are a fatal
 * error in pattern_tiling_init). Keys whose last tile is shorter than k-1
 * are left out: its k-mer count underflows
 * (known-bugs/B-kmer-filter-short-last-tile.md).
 */
GT_TEST(kmer_bound_never_exceeds_distance) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  gt_rng_t rng = gt_rng_new(707);
  static uint8_t key[MAX_LEN], text[MAX_LEN];
  uint64_t trial, within = 0, positive = 0;
  for (trial=0;trial<2400;++trial) {
    const uint64_t m = gt_rng_range(&rng,30,400);
    const uint64_t tiles = 1 + trial % 6;
    const uint64_t k = 3 + (trial/6) % 5;
    const uint64_t max_error = gt_rng_range(&rng,0,m/6);
    const uint64_t tile_length = DIV_CEIL(m,tiles);
    if (tile_length >= k && m-(tiles-1)*tile_length+1 < k) continue;
    random_seq(&rng,key,m,trial%7==0 ? 10 : 0);
    uint64_t n;
    if (trial % 11 == 10) {
      n = m + 2*max_error;
      random_seq(&rng,text,n,0);
    } else {
      n = embed_with_edits(&rng,key,m,gt_rng_below(&rng,max_error+2),
          gt_rng_range(&rng,1,max_error+1),gt_rng_range(&rng,1,max_error+1),text);
    }
    const uint64_t naive = naive_distance(key,m,text,n);
    const uint64_t bound = bound_of(key,m,text,n,tiles,k,max_error,mm_allocator);
    if (bound > 0) ++positive;
    if (naive <= max_error) {
      ++within;
      GT_CHECK_MSG(bound<=naive,
          "trial %" PRIu64 " (m=%" PRIu64 " n=%" PRIu64 " tiles=%" PRIu64 " k=%" PRIu64 " max_error=%" PRIu64 "): bound %" PRIu64 " > distance %" PRIu64,
          trial,m,n,tiles,k,max_error,bound,naive);
    }
  }
  GT_CHECK(within > 1000 && positive > 200);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * alignment_verify_edit_kmer marks a candidate discarded (bound INF)
 * exactly when the k-mer bound exceeds max_error.
 */
GT_TEST(kmer_verify_discards_above_max_error) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  gt_rng_t rng = gt_rng_new(808);
  static uint8_t key[MAX_LEN], text[MAX_LEN];
  uint64_t trial, discarded = 0, kept = 0;
  for (trial=0;trial<600;++trial) {
    const uint64_t m = gt_rng_range(&rng,40,300);
    const uint64_t tiles = 1 + trial % 4;
    const uint64_t k = 4 + trial % 4;
    const uint64_t max_error = gt_rng_range(&rng,1,m/8);
    random_seq(&rng,key,m,0);
    const uint64_t n = (trial % 3 == 0) ?
        (random_seq(&rng,text,m+2*max_error,0), m+2*max_error) :
        embed_with_edits(&rng,key,m,gt_rng_below(&rng,2*max_error+2),max_error,max_error,text);
    pattern_tiled_t tiled;
    memset(&tiled,0,sizeof(tiled));
    pattern_tiled_compile(&tiled,key,m,max_error,tiles,k,false,mm_allocator);
    alignment_t alignment;
    alignment_tile_t tile;
    alignment.alignment_tiles = &tile;
    alignment_init(&alignment,m,0,n,max_error,1,m);
    alignment_verify_edit_kmer(&alignment,&tiled,key,m,text,n,max_error,tiles,k);
    const uint64_t bound = bound_of(key,m,text,n,tiles,k,max_error,mm_allocator);
    if (bound > max_error) {
      ++discarded;
      GT_CHECK_EQ_UINT(alignment.distance_min_bound,ALIGN_DISTANCE_INF);
    } else {
      ++kept;
      GT_CHECK_EQ_UINT(alignment.distance_min_bound,ALIGN_DISTANCE_UNKNOWN);
    }
    pattern_tiled_destroy(&tiled,mm_allocator);
  }
  GT_CHECK(discarded > 50 && kept > 50);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

GT_CASES(
  GT_CASE(kmer_bound_hand_examples),
  GT_CASE(kmer_bound_never_exceeds_distance),
  GT_CASE(kmer_verify_discards_above_max_error)
)
