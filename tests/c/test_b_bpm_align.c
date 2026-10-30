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
 *   BPM (Myers' bit-vector) alignment with traceback, BPM pattern tiles and
 *   the tiled edit-distance filter, against a brute-force dynamic program.
 *   The distance is semi-global (whole key, any text substring); N never
 *   matches. The tiled filter must never discard a text within the error
 *   budget: its distance bound is at most the true distance.
 *
 *   Internal signatures called (a speed iteration may change them; update
 *   this driver with them):
 *     bpm_pattern_compile(bpm_pattern_t*,uint8_t* key,uint64_t length,mm_allocator_t*)
 *     bpm_pattern_compile_tiles(bpm_pattern_t*,uint64_t offset_words64,uint64_t tile_length,bpm_pattern_t* tile)
 *     bpm_compute_edit_distance(const bpm_pattern_t*,const uint8_t* text,uint64_t text_length,
 *         uint64_t* distance,uint64_t* column,uint64_t max_distance,bool quick_abandon)
 *     align_bpm_match(bpm_pattern_t*,uint8_t* key,uint8_t* text,uint64_t text_length,
 *         uint64_t max_distance,bool left_gap_alignment,match_alignment_t*,vector_t* cigar,mm_allocator_t*)
 *     pattern_tiled_compile(pattern_tiled_t*,uint8_t* key,uint64_t key_length,uint64_t max_error,
 *         uint64_t kmer_tiles,uint64_t kmer_length,uint64_t kmer_enabled,mm_allocator_t*)
 *     alignment_init(alignment_t*,key_length,text_begin_offset,text_end_offset,max_error,num_tiles,tile_length)
 *     alignment_verify_edit_bpm(alignment_t*,pattern_tiled_t*,uint8_t* key,uint8_t* text,uint64_t max_error)
 *   and reads pattern_tiled_t (num_tiles, tile_length, tiles[].tile_offset,
 *   tile_length, bpm_pattern_tile), alignment_t.distance_min_bound,
 *   match_alignment_t and cigar_element_t.
 */

#include "gem3_test.h"
#include "text/dna_text.h"
#include "align/alignment.h"
#include "align/align_bpm.h"
#include "align/align_bpm_pattern.h"
#include "align/align_bpm_distance.h"
#include "align/pattern/pattern_tiled.h"
#include "matches/matches_cigar.h"

#define MAX_LEN 1800

/*
 * Brute force: the minimum over text end columns of the edit distance of
 * the whole key against a text substring ending there, and the first
 * (0-based) column reaching it.
 */
static uint64_t naive_distance(
    const uint8_t* const key,const uint64_t m,
    const uint8_t* const text,const uint64_t n,
    uint64_t* const column) {
  static uint64_t prev[MAX_LEN+1], curr[MAX_LEN+1];
  uint64_t best = UINT64_MAX, i, j;
  *column = UINT64_MAX;
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
    if (curr[m] < best) { best = curr[m]; *column = j-1; }
    memcpy(prev,curr,(m+1)*sizeof(uint64_t));
  }
  return best;
}

/*
 * Replays a gem CIGAR (cigar_ins consumes text only, cigar_del key only)
 * from text_offset; returns the text bases consumed or -1, and the number
 * of edits it spells.
 */
static int64_t replay_cigar(
    const uint8_t* const key,const uint64_t key_length,
    const uint8_t* const text,const uint64_t text_length,const uint64_t text_offset,
    vector_t* const cigar_vector,const uint64_t cigar_offset,const uint64_t cigar_length,
    uint64_t* const edits) {
  const cigar_element_t* const cigar = vector_get_elm(cigar_vector,cigar_offset,cigar_element_t);
  uint64_t k = 0, t = text_offset, i;
  *edits = 0;
  for (i=0;i<cigar_length;++i) {
    const cigar_element_t* const e = cigar + i;
    int64_t l;
    switch (e->type) {
      case cigar_match:
        if (e->length <= 0 || k+e->length > key_length || t+e->length > text_length) return -1;
        for (l=0;l<e->length;++l) {
          if (key[k+l]!=text[t+l] || key[k+l]==ENC_DNA_CHAR_N) return -1;
        }
        k += e->length; t += e->length;
        break;
      case cigar_mismatch:
        if (k >= key_length || t >= text_length) return -1;
        if (key[k]==text[t] && key[k]!=ENC_DNA_CHAR_N) return -1;
        if (e->mismatch!=text[t]) return -1;
        ++k; ++t; ++(*edits);
        break;
      case cigar_ins:
        if (e->length <= 0 || t+e->length > text_length) return -1;
        t += e->length; *edits += e->length;
        break;
      case cigar_del:
        if (e->length <= 0 || k+e->length > key_length) return -1;
        k += e->length; *edits += e->length;
        break;
      default:
        return -1;
    }
  }
  if (k!=key_length) return -1;
  return (int64_t)(t-text_offset);
}

static void random_seq(gt_rng_t* const rng,uint8_t* const out,const uint64_t n,const uint64_t n_per_mille) {
  uint64_t i;
  for (i=0;i<n;++i) {
    out[i] = gt_rng_chance(rng,n_per_mille,1000) ? ENC_DNA_CHAR_N : (uint8_t)gt_rng_below(rng,4);
  }
}

// Copy of key with random single-base edits between random flanks
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

static void encode(const char* const ascii,uint8_t* const out,uint64_t* const length) {
  uint64_t i;
  for (i=0;ascii[i];++i) out[i] = dna_encode(ascii[i]);
  *length = i;
}

static mm_allocator_t* new_allocator(mm_slab_t** const slab) {
  *slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM);
  return mm_allocator_new(*slab);
}

static void reset_alignment(match_alignment_t* const ma,vector_t* const cigar_vector) {
  memset(ma,0,sizeof(*ma));
  ma->cigar_offset = vector_get_used(cigar_vector);
}

/*
 * Hand-worked BPM alignments
 */
GT_TEST(bpm_align_hand_examples) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  vector_t* const cigar_vector = vector_new(64,cigar_element_t);
  static const struct {
    const char* key;
    const char* text;
    uint64_t distance;
    uint64_t begin;     // text offset of the first aligned base
    uint64_t end;       // text offset after the last aligned base
  } cases[] = {
    { "ACGTACGTAC", "GGACGTACGTACGG", 0, 2, 12 },
    { "ACGTACGTAC", "GGACGTTCGTACGG", 1, 2, 12 },   // one substitution
    { "ACGTACGTAC", "GGACGTCGTACGG",  1, 2, 11 },   // key base missing from the text
    { "ACGTACGTAC", "GGACGTAACGTACGG", 1, 2, 13 },  // extra text base
    { "ACGTNCGTAC", "GGACGTNCGTACGG", 1, 2, 12 },   // N against N is an edit
  };
  uint64_t c;
  for (c=0;c<sizeof(cases)/sizeof(cases[0]);++c) {
    uint8_t key[64], text[64];
    uint64_t m, n, col;
    encode(cases[c].key,key,&m);
    encode(cases[c].text,text,&n);
    GT_CHECK_EQ_UINT(naive_distance(key,m,text,n,&col),cases[c].distance);
    GT_CHECK_EQ_UINT(col+1,cases[c].end);
    bpm_pattern_t bpm;
    bpm_pattern_compile(&bpm,key,m,mm_allocator);
    match_alignment_t ma;
    reset_alignment(&ma,cigar_vector);
    align_bpm_match(&bpm,key,text,n,5,true,&ma,cigar_vector,mm_allocator);
    GT_CHECK_EQ_UINT(ma.score,cases[c].distance);
    GT_CHECK_EQ_UINT(ma.match_position,cases[c].begin);
    GT_CHECK_EQ_UINT(ma.match_position+ma.effective_length,cases[c].end);
    uint64_t edits;
    const int64_t span = replay_cigar(key,m,text,n,ma.match_position,
        cigar_vector,ma.cigar_offset,ma.cigar_length,&edits);
    GT_CHECK_EQ_INT(span,ma.effective_length);
    GT_CHECK_EQ_UINT(edits,cases[c].distance);
    bpm_pattern_destroy(&bpm,mm_allocator);
  }
  vector_delete(cigar_vector);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * Gap placement in a homopolymer: leftmost with left_gap_alignment,
 * rightmost otherwise.
 */
GT_TEST(bpm_align_left_gap_alignment) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  vector_t* const cigar_vector = vector_new(64,cigar_element_t);
  uint8_t key[32], text[32];
  uint64_t m, n;
  encode("ACGATTTTGCAG",key,&m);   // one T more than the text
  encode("ACGATTTGCAG",text,&n);
  bpm_pattern_t bpm;
  bpm_pattern_compile(&bpm,key,m,mm_allocator);
  int left;
  for (left=0;left<2;++left) {
    match_alignment_t ma;
    reset_alignment(&ma,cigar_vector);
    align_bpm_match(&bpm,key,text,n,3,left==1,&ma,cigar_vector,mm_allocator);
    GT_CHECK_EQ_UINT(ma.score,1);
    GT_REQUIRE(ma.cigar_length==3);
    const cigar_element_t* const e = vector_get_elm(cigar_vector,ma.cigar_offset,cigar_element_t);
    GT_CHECK(e[0].type==cigar_match && e[1].type==cigar_del && e[1].length==1 && e[2].type==cigar_match);
    GT_CHECK_EQ_INT(e[0].length,left ? 4 : 7);
  }
  bpm_pattern_destroy(&bpm,mm_allocator);
  vector_delete(cigar_vector);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * Random keys (1-700 bp, some N) and edited copies: the score is the
 * brute-force distance when it is within max_distance (corrected to at most
 * key length - 1), the alignment ends at the first column reaching it, and
 * the CIGAR spells exactly that many edits. Otherwise no alignment.
 */
GT_TEST(bpm_align_matches_naive) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  vector_t* const cigar_vector = vector_new(4096,cigar_element_t);
  gt_rng_t rng = gt_rng_new(404);
  static uint8_t key[MAX_LEN], text[MAX_LEN];
  static const uint64_t lengths[] = { 1, 2, 5, 31, 63, 64, 65, 100, 127, 128, 129, 150, 200, 255, 256, 257, 400, 640 };
  const uint64_t num_lengths = sizeof(lengths)/sizeof(lengths[0]);
  uint64_t trial;
  for (trial=0;trial<1200;++trial) {
    const uint64_t m = lengths[trial % num_lengths];
    random_seq(&rng,key,m,trial%4==0 ? 15 : 0);
    const uint64_t n = (trial % 9 == 8) ?
        (random_seq(&rng,text,m+20,0), m+20) :
        embed_with_edits(&rng,key,m,gt_rng_below(&rng,m/8+2),gt_rng_below(&rng,30),gt_rng_below(&rng,30),text);
    if (n == 0) continue;
    uint64_t col;
    const uint64_t naive = naive_distance(key,m,text,n,&col);
    const uint64_t max_distance = (trial % 3 == 0) ? naive : gt_rng_range(&rng,0,m/4+2);
    const uint64_t effective = (max_distance >= m) ? m-1 : max_distance;
    bpm_pattern_t bpm;
    bpm_pattern_compile(&bpm,key,m,mm_allocator);
    int left;
    for (left=0;left<2;++left) {
      match_alignment_t ma;
      vector_clear(cigar_vector);
      reset_alignment(&ma,cigar_vector);
      align_bpm_match(&bpm,key,text,n,max_distance,left==1,&ma,cigar_vector,mm_allocator);
      if (naive <= effective) {
        GT_CHECK_MSG((uint64_t)ma.score==naive,"trial %" PRIu64 " (m=%" PRIu64 " n=%" PRIu64 " max=%" PRIu64 "): score %d naive %" PRIu64,
            trial,m,n,max_distance,ma.score,naive);
        uint64_t edits;
        const int64_t span = replay_cigar(key,m,text,n,ma.match_position,
            cigar_vector,ma.cigar_offset,ma.cigar_length,&edits);
        GT_CHECK_MSG(span>=0 && span==ma.effective_length && edits==naive,
            "trial %" PRIu64 " (m=%" PRIu64 "): replay span %" PRId64 " (reported %" PRId64 ") edits %" PRIu64 " naive %" PRIu64,
            trial,m,span,ma.effective_length,edits,naive);
        GT_CHECK_MSG(ma.match_position+ma.effective_length==col+1,
            "trial %" PRIu64 ": alignment ends at %" PRIu64 ", first optimal column %" PRIu64,
            trial,ma.match_position+ma.effective_length,col+1);
      } else {
        GT_CHECK_MSG(ma.score==(int32_t)ALIGN_DISTANCE_INF && ma.cigar_length==0,
            "trial %" PRIu64 ": naive %" PRIu64 " > %" PRIu64 " but score %d, %" PRIu64 " CIGAR elements",
            trial,naive,effective,ma.score,ma.cigar_length);
      }
    }
    bpm_pattern_destroy(&bpm,mm_allocator);
  }
  vector_delete(cigar_vector);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * Long keys are cut into BPM tiles. Every tile covers its own stretch of
 * the key (the stretches partition it), and the tile's BPM pattern computes
 * the same cut-off distance as a pattern compiled from that stretch alone.
 */
GT_TEST(bpm_pattern_tiles_equal_substring_patterns) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  gt_rng_t rng = gt_rng_new(505);
  static uint8_t key[MAX_LEN], text[MAX_LEN];
  static const uint64_t lengths[] = { 100, 300, 513, 700, 1024, 1300 };
  static const uint64_t errors[] = { 0, 5, 30, 130, 300 };
  uint64_t li, ei, trial;
  for (li=0;li<sizeof(lengths)/sizeof(lengths[0]);++li) {
    for (ei=0;ei<sizeof(errors)/sizeof(errors[0]);++ei) {
      const uint64_t m = lengths[li];
      const uint64_t max_error = MIN(errors[ei],m);
      random_seq(&rng,key,m,5);
      pattern_tiled_t tiled;
      memset(&tiled,0,sizeof(tiled));
      pattern_tiled_compile(&tiled,key,m,max_error,2,6,false,mm_allocator);
      GT_REQUIRE(tiled.num_tiles>=1);
      uint64_t t, offset = 0;
      for (t=0;t<tiled.num_tiles;++t) {
        pattern_tile_t* const tile = tiled.tiles + t;
        GT_CHECK_EQ_UINT(tile->tile_offset,offset);
        GT_CHECK(tile->tile_length>0);
        offset += tile->tile_length;
        bpm_pattern_t direct;
        bpm_pattern_compile(&direct,key+tile->tile_offset,tile->tile_length,mm_allocator);
        for (trial=0;trial<6;++trial) {
          const uint64_t n = embed_with_edits(&rng,key+tile->tile_offset,tile->tile_length,
              gt_rng_below(&rng,tile->tile_length/10+2),gt_rng_below(&rng,20),gt_rng_below(&rng,20),text);
          const uint64_t cutoff = gt_rng_range(&rng,0,tile->tile_length/5+1);
          uint64_t d1, c1, d2, c2;
          const bool f1 = bpm_compute_edit_distance(&tile->bpm_pattern_tile,text,n,&d1,&c1,cutoff,false);
          const bool f2 = bpm_compute_edit_distance(&direct,text,n,&d2,&c2,cutoff,false);
          GT_CHECK_MSG(f1==f2 && d1==d2 && (!f1 || c1==c2),
              "key %" PRIu64 " bp, tile %" PRIu64 " (%" PRIu64 "+%" PRIu64 "): tile (%d,%" PRIu64 ",%" PRIu64 ") direct (%d,%" PRIu64 ",%" PRIu64 ")",
              m,t,tile->tile_offset,tile->tile_length,(int)f1,d1,c1,(int)f2,d2,c2);
        }
        bpm_pattern_destroy(&direct,mm_allocator);
      }
      GT_CHECK_EQ_UINT(offset,m);
      pattern_tiled_destroy(&tiled,mm_allocator);
    }
  }
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * The tiled BPM filter (alignment_init + alignment_verify_edit_bpm) on a
 * candidate text of key +- max_error. With one tile the bound is the true
 * distance when that is within max_error, otherwise the candidate is
 * discarded. With several tiles each tile gets its share of max_error, so
 * a candidate with clustered errors can be discarded; a kept candidate's
 * bound is at most the true distance and max_error, and an exact copy is
 * kept with bound 0.
 */
GT_TEST(bpm_tiled_filter_bounds) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  gt_rng_t rng = gt_rng_new(606);
  static uint8_t key[MAX_LEN], text[MAX_LEN];
  static const uint64_t lengths[] = { 60, 100, 150, 250, 400, 600, 900, 1200 };
  uint64_t trial, single_tile = 0, multi_tile_kept = 0, discarded = 0;
  for (trial=0;trial<500;++trial) {
    const uint64_t m = lengths[trial % (sizeof(lengths)/sizeof(lengths[0]))];
    const uint64_t max_error = gt_rng_range(&rng,1,m/8);
    random_seq(&rng,key,m,0);
    const uint64_t edits = (trial % 5 == 0) ? 0 : gt_rng_below(&rng,max_error+3);
    const uint64_t n = embed_with_edits(&rng,key,m,edits,
        gt_rng_below(&rng,max_error+1),gt_rng_below(&rng,max_error+1),text);
    uint64_t col;
    const uint64_t naive = naive_distance(key,m,text,n,&col);
    pattern_tiled_t tiled;
    memset(&tiled,0,sizeof(tiled));
    pattern_tiled_compile(&tiled,key,m,max_error,2,6,false,mm_allocator);
    alignment_t alignment;
    alignment_tile_t tiles[16];
    GT_REQUIRE(tiled.num_tiles <= 16);
    alignment.alignment_tiles = tiles;
    alignment_init(&alignment,m,0,n,max_error,tiled.num_tiles,tiled.tile_length);
    alignment_verify_edit_bpm(&alignment,&tiled,key,text,max_error);
    const uint64_t bound = alignment.distance_min_bound;
    if (bound==ALIGN_DISTANCE_INF) ++discarded;
    if (tiled.num_tiles==1) {
      ++single_tile;
      const uint64_t expected = (naive <= max_error) ? naive : ALIGN_DISTANCE_INF;
      GT_CHECK_MSG(bound==expected,"trial %" PRIu64 " (m=%" PRIu64 " max_error=%" PRIu64 "): bound %" PRIu64 " expected %" PRIu64,
          trial,m,max_error,bound,expected);
    } else {
      if (naive==0) GT_CHECK_MSG(bound==0,"trial %" PRIu64 ": exact copy, bound %" PRIu64,trial,bound);
      if (bound!=ALIGN_DISTANCE_INF) {
        ++multi_tile_kept;
        GT_CHECK_MSG(bound<=naive && bound<=max_error,
            "trial %" PRIu64 " (m=%" PRIu64 " max_error=%" PRIu64 " tiles=%" PRIu64 "): kept with bound %" PRIu64 ", distance %" PRIu64,
            trial,m,max_error,tiled.num_tiles,bound,naive);
      }
    }
    pattern_tiled_destroy(&tiled,mm_allocator);
  }
  // The inputs exercise every branch
  GT_CHECK(single_tile > 100 && multi_tile_kept > 50 && discarded > 20);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

GT_CASES(
  GT_CASE(bpm_align_hand_examples),
  GT_CASE(bpm_align_left_gap_alignment),
  GT_CASE(bpm_align_matches_naive),
  GT_CASE(bpm_pattern_tiles_equal_substring_patterns),
  GT_CASE(bpm_tiled_filter_bounds)
)
