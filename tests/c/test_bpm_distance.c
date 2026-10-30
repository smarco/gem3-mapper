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
 *   BPM (Myers' bit-vector) edit distance against a brute-force dynamic
 *   program. The distance is semi-global: the whole pattern against any
 *   substring of the text ending at the reported column; N never matches.
 *
 *   Internal signatures called (a speed iteration may change them; update
 *   this driver with them):
 *     bpm_pattern_compile(bpm_pattern_t*,uint8_t* pattern,uint64_t length,mm_allocator_t*)
 *     bpm_pattern_destroy(bpm_pattern_t*,mm_allocator_t*)
 *     bpm_compute_edit_distance(const bpm_pattern_t*,const uint8_t* text,uint64_t text_length,
 *         uint64_t* match_distance,uint64_t* match_column,uint64_t max_distance,bool quick_abandon)
 *     bpm_compute_edit_distance_raw(bpm_pattern_t*,const uint8_t* text,uint64_t text_length,
 *         uint64_t* position,uint64_t* distance)
 */

#include "gem3_test.h"
#include "text/dna_text.h"
#include "align/alignment.h"
#include "align/align_bpm_pattern.h"
#include "align/align_bpm_distance.h"

#define MAX_LEN 1200

/*
 * Brute force: D[0][j]=0 (free start in the text), D[i][0]=i,
 * D[i][j] = min(D[i-1][j-1]+mismatch, D[i-1][j]+1, D[i][j-1]+1), mismatch
 * when the symbols differ or either is N. Returns min over j>=1 of D[m][j]
 * and the first (0-based) text column reaching it; UINT64_MAX for an empty
 * text.
 */
static uint64_t naive_distance(
    const uint8_t* const pattern,const uint64_t m,
    const uint8_t* const text,const uint64_t n,
    uint64_t* const column) {
  static uint64_t prev[MAX_LEN+1], curr[MAX_LEN+1];
  uint64_t best = UINT64_MAX, i, j;
  *column = UINT64_MAX;
  for (i=0;i<=m;++i) prev[i] = i;
  for (j=1;j<=n;++j) {
    curr[0] = 0;
    for (i=1;i<=m;++i) {
      const bool match = pattern[i-1]==text[j-1] && pattern[i-1]!=ENC_DNA_CHAR_N;
      uint64_t v = prev[i-1] + (match ? 0 : 1);
      if (prev[i]+1 < v) v = prev[i]+1;
      if (curr[i-1]+1 < v) v = curr[i-1]+1;
      curr[i] = v;
    }
    if (curr[m] < best) {
      best = curr[m];
      *column = j-1;
    }
    memcpy(prev,curr,(m+1)*sizeof(uint64_t));
  }
  return best;
}

static void encode(const char* const ascii,uint8_t* const out,uint64_t* const length) {
  uint64_t i;
  for (i=0;ascii[i];++i) out[i] = dna_encode(ascii[i]);
  *length = i;
}

static void random_seq(gt_rng_t* const rng,uint8_t* const out,const uint64_t n,const uint64_t n_per_mille) {
  uint64_t i;
  for (i=0;i<n;++i) {
    out[i] = gt_rng_chance(rng,n_per_mille,1000) ? ENC_DNA_CHAR_N : (uint8_t)gt_rng_below(rng,4);
  }
}

// Copy of the pattern with random edits, between random flanks
static uint64_t embed_with_edits(
    gt_rng_t* const rng,const uint8_t* const pattern,const uint64_t m,
    const uint64_t edits,uint8_t* const text,const uint64_t max_text) {
  uint8_t body[MAX_LEN];
  uint64_t b = 0, i, e;
  for (i=0;i<m && b<MAX_LEN;++i) body[b++] = pattern[i];
  for (e=0;e<edits && b>0;++e) {
    const uint64_t pos = gt_rng_below(rng,b);
    switch (gt_rng_below(rng,3)) {
      case 0: body[pos] = (uint8_t)((body[pos]+1+gt_rng_below(rng,3))%4); break; // substitution
      case 1: if (b < MAX_LEN) { memmove(body+pos+1,body+pos,b-pos); body[pos]=(uint8_t)gt_rng_below(rng,4); ++b; } break;
      default: memmove(body+pos,body+pos+1,b-pos-1); --b; break; // deletion
    }
  }
  const uint64_t room = max_text > b ? max_text-b : 0;
  const uint64_t left = room ? gt_rng_below(rng,room/2+1) : 0;
  const uint64_t right = room ? gt_rng_below(rng,room/2+1) : 0;
  random_seq(rng,text,left,0);
  memcpy(text+left,body,b);
  random_seq(rng,text+left+b,right,0);
  return left+b+right;
}

static mm_allocator_t* new_allocator(mm_slab_t** const slab) {
  *slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM);
  return mm_allocator_new(*slab);
}

static const uint64_t pattern_lengths[] = {
  1, 2, 3, 7, 31, 32, 33, 63, 64, 65, 100, 127, 128, 129, 150, 191, 192, 193, 250, 256, 257, 300, 513
};
#define NUM_PATTERN_LENGTHS (sizeof(pattern_lengths)/sizeof(pattern_lengths[0]))

/*
 * Hand-worked examples
 */
GT_TEST(bpm_hand_examples) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  static const struct {
    const char* pattern;
    const char* text;
    uint64_t distance;
    uint64_t column;
  } cases[] = {
    { "ACGT", "TTACGTTT", 0, 5 },   // exact occurrence ending at text[5]
    { "ACGT", "ACT",      1, 2 },   // one deletion from the pattern
    { "ACGT", "AACCGTT",  1, 5 },   // one insertion in the text (ACCGT)
    { "ACGT", "AGGT",     1, 3 },   // one substitution
    { "AAAA", "CCCC",     4, 0 },   // nothing matches: every column scores 4, first wins
    { "ANA",  "ANA",      1, 2 },   // N never matches, even N against N
    { "GATTACA", "GATTACAGATTACA", 0, 6 }, // the first of two exact occurrences
  };
  uint64_t c;
  for (c=0;c<sizeof(cases)/sizeof(cases[0]);++c) {
    uint8_t p[64], t[64];
    uint64_t m, n, naive_col;
    encode(cases[c].pattern,p,&m);
    encode(cases[c].text,t,&n);
    // The brute force agrees with the hand-worked value
    GT_CHECK_EQ_UINT(naive_distance(p,m,t,n,&naive_col),cases[c].distance);
    GT_CHECK_EQ_UINT(naive_col,cases[c].column);
    bpm_pattern_t bpm;
    bpm_pattern_compile(&bpm,p,m,mm_allocator);
    uint64_t distance = 0, column = 0;
    GT_CHECK(bpm_compute_edit_distance_raw(&bpm,t,n,&column,&distance));
    GT_CHECK_MSG(distance==cases[c].distance && column==cases[c].column,
        "raw %s in %s: (%" PRIu64 ",%" PRIu64 ") expected (%" PRIu64 ",%" PRIu64 ")",
        cases[c].pattern,cases[c].text,distance,column,cases[c].distance,cases[c].column);
    // Cut-off at m (corrected to m-1 inside): found exactly when distance < m
    const bool found = bpm_compute_edit_distance(&bpm,t,n,&distance,&column,m,false);
    GT_CHECK_EQ_INT(found,cases[c].distance < m);
    if (found) {
      GT_CHECK_EQ_UINT(distance,cases[c].distance);
      GT_CHECK_EQ_UINT(column,cases[c].column);
    } else {
      GT_CHECK_EQ_UINT(distance,ALIGN_DISTANCE_INF);
    }
    bpm_pattern_destroy(&bpm,mm_allocator);
  }
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * Raw distance (no cut-off) equals the brute force on random inputs
 */
GT_TEST(bpm_raw_matches_naive) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  gt_rng_t rng = gt_rng_new(11);
  static uint8_t pattern[MAX_LEN], text[MAX_LEN];
  uint64_t trial;
  for (trial=0;trial<2000;++trial) {
    const uint64_t m = pattern_lengths[trial % NUM_PATTERN_LENGTHS];
    random_seq(&rng,pattern,m,trial%3==0 ? 20 : 0);
    uint64_t n;
    if (trial % 4 == 3) {
      n = gt_rng_range(&rng,1,400);
      random_seq(&rng,text,n,10);
    } else {
      n = embed_with_edits(&rng,pattern,m,gt_rng_below(&rng,m/8+2),text,m+gt_rng_below(&rng,120));
      if (n == 0) { text[0] = 0; n = 1; }
    }
    uint64_t naive_col;
    const uint64_t naive = naive_distance(pattern,m,text,n,&naive_col);
    bpm_pattern_t bpm;
    bpm_pattern_compile(&bpm,pattern,m,mm_allocator);
    uint64_t distance = 0, column = 0;
    const bool found = bpm_compute_edit_distance_raw(&bpm,text,n,&column,&distance);
    GT_CHECK_MSG(found && distance==naive && column==naive_col,
        "trial %" PRIu64 " (m=%" PRIu64 ", n=%" PRIu64 "): bpm (%" PRIu64 ",%" PRIu64 ") naive (%" PRIu64 ",%" PRIu64 ")",
        trial,m,n,distance,column,naive,naive_col);
    bpm_pattern_destroy(&bpm,mm_allocator);
  }
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * Cut-off distance: reports a match exactly when the brute-force minimum is
 * within min(max_distance, m-1), with the same distance and first column;
 * with and without quick abandon.
 */
GT_TEST(bpm_cutoff_matches_naive) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  gt_rng_t rng = gt_rng_new(23);
  static uint8_t pattern[MAX_LEN], text[MAX_LEN];
  uint64_t trial;
  for (trial=0;trial<2000;++trial) {
    const uint64_t m = pattern_lengths[trial % NUM_PATTERN_LENGTHS];
    random_seq(&rng,pattern,m,trial%5==0 ? 20 : 0);
    uint64_t n;
    if (trial % 6 == 5) {
      n = gt_rng_range(&rng,1,400);
      random_seq(&rng,text,n,0);
    } else {
      n = embed_with_edits(&rng,pattern,m,gt_rng_below(&rng,m/6+2),text,m+gt_rng_below(&rng,150));
      if (n == 0) { text[0] = 0; n = 1; }
    }
    uint64_t naive_col;
    const uint64_t naive = naive_distance(pattern,m,text,n,&naive_col);
    bpm_pattern_t bpm;
    bpm_pattern_compile(&bpm,pattern,m,mm_allocator);
    const uint64_t thresholds[] = {
      0, naive>0 ? naive-1 : 0, naive, naive+1, naive+7, m/10, m/4, m-1, m, m+5, 1000
    };
    uint64_t k;
    for (k=0;k<sizeof(thresholds)/sizeof(thresholds[0]);++k) {
      const uint64_t max_distance = thresholds[k];
      const uint64_t effective = max_distance >= m ? m-1 : max_distance;
      const bool expect = naive <= effective;
      int qa;
      for (qa=0;qa<2;++qa) {
        uint64_t distance = 0, column = 0;
        const bool found = bpm_compute_edit_distance(&bpm,text,n,&distance,&column,max_distance,qa==1);
        if (expect) {
          GT_CHECK_MSG(found && distance==naive && column==naive_col,
              "trial %" PRIu64 " (m=%" PRIu64 ", n=%" PRIu64 ", max=%" PRIu64 ", quick_abandon=%d): "
              "found=%d (%" PRIu64 ",%" PRIu64 ") naive (%" PRIu64 ",%" PRIu64 ")",
              trial,m,n,max_distance,qa,(int)found,distance,column,naive,naive_col);
        } else {
          GT_CHECK_MSG(!found && distance==ALIGN_DISTANCE_INF,
              "trial %" PRIu64 " (m=%" PRIu64 ", n=%" PRIu64 ", max=%" PRIu64 ", quick_abandon=%d): "
              "found=%d distance %" PRIu64 " but naive %" PRIu64 " > %" PRIu64,
              trial,m,n,max_distance,qa,(int)found,distance,naive,effective);
        }
      }
    }
    bpm_pattern_destroy(&bpm,mm_allocator);
  }
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

GT_CASES(
  GT_CASE(bpm_hand_examples),
  GT_CASE(bpm_raw_matches_naive),
  GT_CASE(bpm_cutoff_matches_naive)
)
