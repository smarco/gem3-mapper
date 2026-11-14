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
 *   Smith-Waterman-Gotoh (gap-affine) alignment against a brute-force Gotoh
 *   dynamic program, and the SWG score helpers against hand-worked values.
 *
 *   The key is aligned end to end; the text may be entered (begin_free) and
 *   left (end_free) anywhere. A gap of length k scores gap_open + k*gap_ext
 *   (both negative); N never matches, not even N. gem's CIGAR types:
 *   cigar_ins consumes text only (SAM D), cigar_del consumes key only (SAM I).
 *
 *   Internal signatures called (a speed iteration may change them; update
 *   this driver with them):
 *     align_swg(key,key_length,text,text_length,const swg_penalties_t*,max_bandwidth,
 *         begin_free,end_free,left_gap_alignment,match_alignment_t*,vector_t* cigar,mm_allocator_t*)
 *     align_swg_banded(<same arguments as align_swg>)
 *     align_swg_base(key,key_length,text,text_length,const swg_penalties_t*,left_gap_alignment,
 *         match_alignment_t*,vector_t* cigar,mm_allocator_t*)
 *     align_swg_score_{deletion,insertion}(const swg_penalties_t*,int32_t length)
 *     align_swg_score_mismatch(const swg_penalties_t*), align_swg_score_match(p,int32_t length)
 *     align_swg_score_cigar_element(p,const cigar_element_t*)
 *     align_swg_score_cigar_excluding_clipping(p,vector_t*,cigar_offset,cigar_length)
 *     align_swg_score_compute_max_score_bound(p,edit_distance,key_length)
 *     align_swg_score_compute_max_edit_bound(p,swg_score,key_length)
 *     search_parameters_init(search_parameters_t*), search_instantiate_values(search_parameters_t*,uint64_t)
 *     search_configure_alignment_match_scores(search_parameters_t*,int32_t)
 *     search_configure_alignment_mismatch_scores(search_parameters_t*,int32_t)
 *   and reads these fields: search_parameters_t (mapping_mode, the alignment
 *   fractions and their *_nominal values, swg_penalties, check_type,
 *   mapq_model_se, candidate_verification.verification_strategy),
 *   match_alignment_t, cigar_element_t.
 */

#include "gem3_test.h"
#include "text/dna_text.h"
#include "align/align_swg.h"
#include "align/align_swg_banded.h"
#include "align/align_swg_score.h"
#include "matches/matches_cigar.h"
#include "archive/search/archive_search_se_parameters.h"

#define MAX_LEN 320
#define NEG_INF (-1000000000)

typedef struct {
  int32_t match;     // > 0
  int32_t mismatch;  // penalty > 0
  int32_t open;      // penalty >= 0
  int32_t ext;       // penalty > 0
} penalties_t;

static void set_penalties(swg_penalties_t* const p,const penalties_t* const q) {
  int a, b;
  for (a=0;a<DNA__N_RANGE;++a) {
    for (b=0;b<DNA__N_RANGE;++b) {
      p->matching_score[a][b] = (a==b && a!=ENC_DNA_CHAR_N) ? q->match : -q->mismatch;
    }
  }
  p->generic_match_score = q->match;
  p->generic_mismatch_score = -q->mismatch;
  p->gap_open_score = -q->open;
  p->gap_extension_score = -q->ext;
}

static int32_t sub_score(const penalties_t* const q,const uint8_t a,const uint8_t b) {
  return (a==b && a!=ENC_DNA_CHAR_N) ? q->match : -q->mismatch;
}

/*
 * Brute-force Gotoh: H (any state), E (gap consuming text), F (gap consuming
 * key). The key is global; the text prefix is free with begin_free and the
 * suffix with end_free. Returns the optimal score.
 */
static int32_t gotoh(
    const uint8_t* const key,const int m,const uint8_t* const text,const int n,
    const penalties_t* const q,const bool begin_free,const bool end_free) {
  static int32_t H[MAX_LEN+1][MAX_LEN+1], E[MAX_LEN+1][MAX_LEN+1], F[MAX_LEN+1][MAX_LEN+1];
  const int32_t go = q->open, ge = q->ext;
  int i, j;
  for (i=0;i<=m;++i) {
    for (j=0;j<=n;++j) {
      if (i==0 && j==0) { H[0][0] = 0; E[0][0] = F[0][0] = NEG_INF; continue; }
      int32_t e = NEG_INF, f = NEG_INF, h = NEG_INF;
      if (j>0) e = MAX(H[i][j-1]-go-ge,E[i][j-1]-ge);
      if (i>0) f = MAX(H[i-1][j]-go-ge,F[i-1][j]-ge);
      if (i>0 && j>0) h = H[i-1][j-1] + sub_score(q,key[i-1],text[j-1]);
      if (i==0 && begin_free) h = 0;
      E[i][j] = e; F[i][j] = f;
      H[i][j] = MAX(h,MAX(e,f));
    }
  }
  if (!end_free) return H[m][n];
  int32_t best = NEG_INF;
  for (j=0;j<=n;++j) best = MAX(best,H[m][j]);
  return best;
}

/*
 * Replays a gem CIGAR over key and text (from text_offset): checks each
 * element against the sequences, scores it independently, and returns the
 * number of text bases consumed (or -1 on an inconsistency, with a note).
 */
static int64_t replay_cigar(
    const uint8_t* const key,const uint64_t key_length,
    const uint8_t* const text,const uint64_t text_length,const uint64_t text_offset,
    vector_t* const cigar_vector,const uint64_t cigar_offset,const uint64_t cigar_length,
    const penalties_t* const q,int32_t* const score,uint64_t* const edits) {
  const cigar_element_t* const cigar = vector_get_elm(cigar_vector,cigar_offset,cigar_element_t);
  uint64_t k = 0, t = text_offset, i;
  int32_t s = 0;
  *edits = 0;
  for (i=0;i<cigar_length;++i) {
    const cigar_element_t* const e = cigar + i;
    switch (e->type) {
      case cigar_match: {
        int64_t l;
        if (e->length <= 0 || k+e->length > key_length || t+e->length > text_length) return -1;
        for (l=0;l<e->length;++l) {
          if (key[k+l]!=text[t+l] || key[k]==ENC_DNA_CHAR_N || key[k+l]==ENC_DNA_CHAR_N) {
            gt_note("match element over unequal bases at key %" PRIu64,k+l);
            return -1;
          }
        }
        s += q->match*e->length; k += e->length; t += e->length;
        break;
      }
      case cigar_mismatch:
        if (k >= key_length || t >= text_length) return -1;
        if (key[k]==text[t] && key[k]!=ENC_DNA_CHAR_N) { gt_note("mismatch over equal bases"); return -1; }
        if (e->mismatch!=text[t]) { gt_note("mismatch element records %u, text has %u",e->mismatch,text[t]); return -1; }
        s -= q->mismatch; ++k; ++t; ++(*edits);
        break;
      case cigar_ins:
        if (e->length <= 0 || t+e->length > text_length) return -1;
        s -= q->open + q->ext*e->length; t += e->length; *edits += e->length;
        break;
      case cigar_del:
        if (e->length <= 0 || k+e->length > key_length) return -1;
        s -= q->open + q->ext*e->length; k += e->length; *edits += e->length;
        break;
      default:
        gt_note("unexpected CIGAR element type %d",(int)e->type);
        return -1;
    }
  }
  if (k!=key_length) { gt_note("CIGAR consumes %" PRIu64 " of %" PRIu64 " key bases",k,key_length); return -1; }
  *score = s;
  return (int64_t)(t-text_offset);
}

/* Adjacent elements of one type are merged and no element is empty */
static bool cigar_is_canonical(vector_t* const cigar_vector,const uint64_t offset,const uint64_t length) {
  const cigar_element_t* const cigar = vector_get_elm(cigar_vector,offset,cigar_element_t);
  uint64_t i;
  for (i=0;i<length;++i) {
    if (cigar[i].type!=cigar_mismatch && cigar[i].length<=0) return false;
    if (i>0 && cigar[i].type==cigar[i-1].type && cigar[i].type!=cigar_mismatch) return false;
  }
  return true;
}

static void random_seq(gt_rng_t* const rng,uint8_t* const out,const uint64_t n,const uint64_t n_per_mille) {
  uint64_t i;
  for (i=0;i<n;++i) {
    out[i] = gt_rng_chance(rng,n_per_mille,1000) ? ENC_DNA_CHAR_N : (uint8_t)gt_rng_below(rng,4);
  }
}

// Copy of key with random edits (subs, indels of 1-4 bases) and random flanks
static uint64_t mutate(
    gt_rng_t* const rng,const uint8_t* const key,const uint64_t m,const uint64_t edits,
    const uint64_t left,const uint64_t right,uint8_t* const out) {
  uint8_t body[MAX_LEN];
  uint64_t b = 0, i, e;
  for (i=0;i<m;++i) body[b++] = key[i];
  for (e=0;e<edits && b>0;++e) {
    const uint64_t pos = gt_rng_below(rng,b);
    const uint64_t len = gt_rng_range(rng,1,4);
    switch (gt_rng_below(rng,3)) {
      case 0: body[pos] = (uint8_t)((body[pos]+1+gt_rng_below(rng,3))%4); break;
      case 1:
        if (b+len < MAX_LEN-2*40) {
          memmove(body+pos+len,body+pos,b-pos);
          for (i=0;i<len;++i) body[pos+i] = (uint8_t)gt_rng_below(rng,4);
          b += len;
        }
        break;
      default:
        if (pos+len <= b) { memmove(body+pos,body+pos+len,b-pos-len); b -= len; }
        break;
    }
  }
  random_seq(rng,out,left,0);
  memcpy(out+left,body,b);
  random_seq(rng,out+left+b,right,0);
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
  ma->match_position = 0;
  ma->match_text_offset = 0;
  ma->cigar_offset = vector_get_used(cigar_vector);
  ma->cigar_length = 0;
  ma->effective_length = 0;
  ma->score = 0;
}

static const penalties_t default_penalties = { 1, 4, 6, 1 };

/*
 * Default penalties and the score helpers (README: 1,4,6,1, a k-base gap
 * costs 6 + k)
 */
GT_TEST(swg_score_helpers_hand) {
  search_parameters_t sp;
  memset(&sp,0,sizeof(sp));
  search_parameters_init(&sp);
  swg_penalties_t* const p = &sp.swg_penalties;
  int a, b;
  for (a=0;a<DNA__N_RANGE;++a) {
    for (b=0;b<DNA__N_RANGE;++b) {
      GT_CHECK_EQ_INT(p->matching_score[a][b],(a==b && a!=ENC_DNA_CHAR_N) ? 1 : -4);
    }
  }
  GT_CHECK_EQ_INT(p->gap_open_score,-6);
  GT_CHECK_EQ_INT(p->gap_extension_score,-1);
  GT_CHECK_EQ_INT(align_swg_score_match(p,1),1);
  GT_CHECK_EQ_INT(align_swg_score_match(p,37),37);
  GT_CHECK_EQ_INT(align_swg_score_mismatch(p),-4);
  GT_CHECK_EQ_INT(align_swg_score_deletion(p,1),-7);
  GT_CHECK_EQ_INT(align_swg_score_deletion(p,5),-11);
  GT_CHECK_EQ_INT(align_swg_score_insertion(p,1),-7);
  GT_CHECK_EQ_INT(align_swg_score_insertion(p,3),-9);
  // Elements
  cigar_element_t e;
  memset(&e,0,sizeof(e));
  e.type = cigar_match; e.length = 12;
  GT_CHECK_EQ_INT(align_swg_score_cigar_element(p,&e),12);
  e.type = cigar_mismatch; e.mismatch = ENC_DNA_CHAR_G;
  GT_CHECK_EQ_INT(align_swg_score_cigar_element(p,&e),-4);
  e.type = cigar_ins; e.length = 2;
  GT_CHECK_EQ_INT(align_swg_score_cigar_element(p,&e),-8);
  e.type = cigar_del; e.length = 4;
  GT_CHECK_EQ_INT(align_swg_score_cigar_element(p,&e),-10);
  // A whole CIGAR, with and without the leading/trailing deletions that
  // excluding_clipping skips: 2D 10M X 5M 1I 3M 3D
  vector_t* const cigar_vector = vector_new(16,cigar_element_t);
  uint64_t len = 0;
  matches_cigar_vector_append_deletion(cigar_vector,&len,2,cigar_attr_none);
  matches_cigar_vector_append_match(cigar_vector,&len,10,cigar_attr_none);
  matches_cigar_vector_append_mismatch(cigar_vector,&len,ENC_DNA_CHAR_C,cigar_attr_none);
  matches_cigar_vector_append_match(cigar_vector,&len,5,cigar_attr_none);
  matches_cigar_vector_append_insertion(cigar_vector,&len,1,cigar_attr_none);
  matches_cigar_vector_append_match(cigar_vector,&len,3,cigar_attr_none);
  matches_cigar_vector_append_deletion(cigar_vector,&len,3,cigar_attr_none);
  GT_REQUIRE(len==7);
  GT_CHECK_EQ_INT(align_swg_score_cigar(p,cigar_vector,0,len),-8+10-4+5-7+3-9);
  GT_CHECK_EQ_INT(align_swg_score_cigar_excluding_clipping(p,cigar_vector,0,len),10-4+5-7+3);
  GT_CHECK_EQ_INT(align_swg_score_cigar_excluding_clipping(p,cigar_vector,1,5),10-4+5-7+3);
  GT_CHECK_EQ_INT(align_swg_score_cigar_excluding_clipping(p,cigar_vector,0,0),0);
  // Score bounds (pins the formulas: L-e matches plus the lighter of one
  // e-base gap or e mismatches)
  GT_CHECK_EQ_INT(align_swg_score_compute_max_score_bound(p,0,100),100);
  GT_CHECK_EQ_INT(align_swg_score_compute_max_score_bound(p,1,100),99-4);
  GT_CHECK_EQ_INT(align_swg_score_compute_max_score_bound(p,3,100),97-9);
  GT_CHECK_EQ_INT(align_swg_score_compute_max_score_bound(p,10,150),140-16);
  GT_CHECK_EQ_INT(align_swg_score_compute_max_edit_bound(p,80,100),14);  // ceil(20/4)=5, 20-6=14
  GT_CHECK_EQ_INT(align_swg_score_compute_max_edit_bound(p,97,100),100); // diff 3 < g(1): a gap cannot explain it
  GT_CHECK_EQ_INT(align_swg_score_compute_max_edit_bound(p,50,100),44);
  vector_delete(cigar_vector);
}

/*
 * Defaults stated in `gem-mapper --help=developer` and their instantiation
 * to a read length (a value below 1 is a fraction of the read length, 1 or
 * more an absolute number). --alignment-global-min-score and
 * --alignment-curation-min-end-context are left out: their --help defaults
 * disagree with the code (known-bugs/B-help-defaults.md).
 */
GT_TEST(search_parameters_defaults_and_instantiation) {
  search_parameters_t sp;
  memset(&sp,0,sizeof(sp));
  search_parameters_init(&sp);
  GT_CHECK_EQ_INT(sp.mapping_mode,mapping_adaptive_filtering_fast);
  GT_CHECK(sp.complete_search_error==0.04);
  GT_CHECK(sp.complete_strata_after_best==1.0);
  GT_CHECK(sp.alignment_max_error==0.12);
  GT_CHECK(sp.alignment_max_bandwidth==0.20);
  GT_CHECK(sp.alignment_global_min_identity==0.80);
  GT_CHECK_EQ_INT(sp.alignment_local,local_alignment_if_unmapped);
  GT_CHECK(sp.alignment_local_min_identity==40.0);
  GT_CHECK(sp.alignment_local_min_swg_threshold==20.0);
  GT_CHECK(!sp.alignment_force_full_swg);
  GT_CHECK(sp.alignment_scaffolding_min_matching_length==10.0);
  GT_CHECK(sp.cigar_curation);
  GT_CHECK_EQ_INT(sp.match_alignment_model,match_alignment_model_gap_affine);
  GT_CHECK_EQ_INT(sp.candidate_verification.verification_strategy,verification_chained);
  GT_CHECK_EQ_INT(sp.mapq_model_se,mapq_model_gem);
  GT_CHECK_EQ_INT(sp.check_type,archive_check_nothing);
  GT_CHECK_EQ_INT(sp.swg_penalties.generic_match_score,1);
  GT_CHECK_EQ_INT(sp.swg_penalties.generic_mismatch_score,-4);
  // 100 bp
  search_instantiate_values(&sp,100);
  GT_CHECK_EQ_UINT(sp.complete_search_error_nominal,4);
  GT_CHECK_EQ_UINT(sp.complete_strata_after_best_nominal,1);
  GT_CHECK_EQ_UINT(sp.alignment_max_error_nominal,12);
  GT_CHECK_EQ_UINT(sp.alignment_max_bandwidth_nominal,20);
  GT_CHECK_EQ_UINT(sp.alignment_global_min_identity_nominal,80);
  GT_CHECK_EQ_UINT(sp.alignment_local_min_identity_nominal,40);
  GT_CHECK_EQ_INT(sp.alignment_local_min_swg_threshold_nominal,20);
  GT_CHECK_EQ_UINT(sp.alignment_scaffolding_min_matching_length_nominal,10);
  // 36 bp: fractions round down, absolute values stay
  search_instantiate_values(&sp,36);
  GT_CHECK_EQ_UINT(sp.complete_search_error_nominal,1);
  GT_CHECK_EQ_UINT(sp.alignment_max_error_nominal,4);
  GT_CHECK_EQ_UINT(sp.alignment_max_bandwidth_nominal,7);
  GT_CHECK_EQ_UINT(sp.alignment_global_min_identity_nominal,28);
  GT_CHECK_EQ_UINT(sp.alignment_local_min_identity_nominal,40);
  // Changed values
  sp.alignment_max_error = 7.0;
  sp.complete_search_error = 0.5;
  sp.alignment_global_min_identity = 0.95;
  search_instantiate_values(&sp,250);
  GT_CHECK_EQ_UINT(sp.alignment_max_error_nominal,7);
  GT_CHECK_EQ_UINT(sp.complete_search_error_nominal,125);
  GT_CHECK_EQ_UINT(sp.alignment_global_min_identity_nominal,237);
  // -A sets every match cell and the generic score; -B every mismatch cell
  // (its generic score is known-bugs/B-mismatch-penalty-sign.md)
  search_configure_alignment_match_scores(&sp,3);
  search_configure_alignment_mismatch_scores(&sp,7);
  int a, b;
  for (a=0;a<DNA__N_RANGE;++a) {
    for (b=0;b<DNA__N_RANGE;++b) {
      GT_CHECK_EQ_INT(sp.swg_penalties.matching_score[a][b],(a==b && a!=ENC_DNA_CHAR_N) ? 3 : -7);
    }
  }
  GT_CHECK_EQ_INT(sp.swg_penalties.generic_match_score,3);
}

/*
 * Hand-worked alignments (default penalties)
 */
GT_TEST(swg_hand_examples) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  vector_t* const cigar_vector = vector_new(64,cigar_element_t);
  swg_penalties_t p;
  set_penalties(&p,&default_penalties);
  static const struct {
    const char* key;
    const char* text;
    bool begin_free, end_free;
    int32_t score;
    uint64_t position;   // text offset where the alignment starts
    uint64_t span;       // text bases covered
  } cases[] = {
    { "ACGTACGTAC", "ACGTACGTAC", false, false, 10, 0, 10 },
    { "ACGTACGTAC", "ACGTTCGTAC", false, false, 9-4, 0, 10 },      // one mismatch
    { "ACGTACGTAC", "TTTACGTACGTACTT", true, true, 10, 3, 10 },    // free flanks
    { "ACGTACGTAC", "ACGTAGTAC", false, false, 9-7, 0, 9 },        // key base missing from text (1I)
    { "ACGTAGTAC", "ACGTACGTAC", false, false, 9-7, 0, 10 },       // text base missing from key (1D)
    { "ACGTANGTAC", "ACGTANGTAC", false, false, 9-4, 0, 10 },      // N never matches
    { "AAAACCCCGGGGTTTT", "AAAACCCCTTTT", false, false, 12-10, 0, 12 }, // one 4-base gap (6+4)
  };
  uint64_t c;
  for (c=0;c<sizeof(cases)/sizeof(cases[0]);++c) {
    uint8_t key[64], text[64];
    uint64_t m, n;
    encode(cases[c].key,key,&m);
    encode(cases[c].text,text,&n);
    GT_CHECK_EQ_INT(gotoh(key,(int)m,text,(int)n,&default_penalties,cases[c].begin_free,cases[c].end_free),cases[c].score);
    match_alignment_t ma;
    reset_alignment(&ma,cigar_vector);
    align_swg(key,m,text,n,&p,64,cases[c].begin_free,cases[c].end_free,true,&ma,cigar_vector,mm_allocator);
    GT_CHECK_MSG(ma.score==cases[c].score,"case %" PRIu64 ": score %d expected %d",c,ma.score,cases[c].score);
    GT_CHECK_EQ_UINT(ma.match_position,cases[c].position);
    GT_CHECK_EQ_UINT(ma.effective_length,cases[c].span);
    int32_t s; uint64_t edits;
    const int64_t span = replay_cigar(key,m,text,n,ma.match_position,
        cigar_vector,ma.cigar_offset,ma.cigar_length,&default_penalties,&s,&edits);
    GT_CHECK_EQ_INT(span,cases[c].span);
    GT_CHECK_EQ_INT(s,cases[c].score);
  }
  vector_delete(cigar_vector);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * The trivial shapes align_swg handles without a DP
 */
GT_TEST(swg_trivial_shapes) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  vector_t* const cigar_vector = vector_new(64,cigar_element_t);
  swg_penalties_t p;
  set_penalties(&p,&default_penalties);
  uint8_t key[8], text[8];
  uint64_t m, n;
  match_alignment_t ma;
  const cigar_element_t* e;
  // Empty key, text of 5, global: the whole text is a gap (6+5)
  encode("ACGTA",text,&n);
  reset_alignment(&ma,cigar_vector);
  align_swg(key,0,text,n,&p,10,false,false,true,&ma,cigar_vector,mm_allocator);
  GT_CHECK_EQ_INT(ma.score,-11);
  GT_REQUIRE(ma.cigar_length==1);
  e = vector_get_elm(cigar_vector,ma.cigar_offset,cigar_element_t);
  GT_CHECK(e->type==cigar_ins && e->length==5);
  // Empty key, begin free: no CIGAR, position moves past the text
  reset_alignment(&ma,cigar_vector);
  align_swg(key,0,text,n,&p,10,true,false,true,&ma,cigar_vector,mm_allocator);
  GT_CHECK_EQ_INT(ma.score,0);
  GT_CHECK_EQ_UINT(ma.cigar_length,0);
  GT_CHECK_EQ_UINT(ma.match_position,5);
  // Key of 3, empty text: the whole key is a gap (6+3)
  encode("ACG",key,&m);
  reset_alignment(&ma,cigar_vector);
  align_swg(key,m,text,0,&p,10,false,false,true,&ma,cigar_vector,mm_allocator);
  GT_CHECK_EQ_INT(ma.score,-9);
  GT_REQUIRE(ma.cigar_length==1);
  e = vector_get_elm(cigar_vector,ma.cigar_offset,cigar_element_t);
  GT_CHECK(e->type==cigar_del && e->length==3);
  // One base each: match, mismatch, N against N
  static const char* const pairs[][2] = { {"A","A"}, {"A","C"}, {"N","N"}, {"G","N"} };
  static const int32_t scores[] = { 1, -4, -4, -4 };
  uint64_t i;
  for (i=0;i<4;++i) {
    encode(pairs[i][0],key,&m);
    encode(pairs[i][1],text,&n);
    reset_alignment(&ma,cigar_vector);
    align_swg(key,1,text,1,&p,10,false,false,true,&ma,cigar_vector,mm_allocator);
    GT_CHECK_EQ_INT(ma.score,scores[i]);
    GT_REQUIRE(ma.cigar_length==1);
    e = vector_get_elm(cigar_vector,ma.cigar_offset,cigar_element_t);
    if (i==0) {
      GT_CHECK(e->type==cigar_match && e->length==1);
    } else {
      GT_CHECK(e->type==cigar_mismatch && e->mismatch==text[0]);
    }
  }
  vector_delete(cigar_vector);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * Gap placement in a homopolymer: left_gap_alignment puts the gap at the
 * leftmost equivalent position, otherwise at the rightmost.
 */
GT_TEST(swg_left_gap_alignment) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  vector_t* const cigar_vector = vector_new(64,cigar_element_t);
  swg_penalties_t p;
  set_penalties(&p,&default_penalties);
  uint8_t key[32], text[32];
  uint64_t m, n;
  // Key has one T more than the text: a cigar_del inside TTTT
  encode("ACGATTTTGCAG",key,&m);
  encode("ACGATTTGCAG",text,&n);
  int left;
  for (left=0;left<2;++left) {
    match_alignment_t ma;
    reset_alignment(&ma,cigar_vector);
    align_swg(key,m,text,n,&p,10,false,false,left==1,&ma,cigar_vector,mm_allocator);
    GT_CHECK_EQ_INT(ma.score,11-7);
    GT_REQUIRE(ma.cigar_length==3);
    const cigar_element_t* const e = vector_get_elm(cigar_vector,ma.cigar_offset,cigar_element_t);
    GT_CHECK(e[0].type==cigar_match && e[1].type==cigar_del && e[1].length==1 && e[2].type==cigar_match);
    GT_CHECK_EQ_INT(e[0].length,left ? 4 : 7);
  }
  // Text has one T more than the key: a cigar_ins inside TTTT
  for (left=0;left<2;++left) {
    match_alignment_t ma;
    reset_alignment(&ma,cigar_vector);
    align_swg(text,n,key,m,&p,10,false,false,left==1,&ma,cigar_vector,mm_allocator);
    GT_CHECK_EQ_INT(ma.score,11-7);
    GT_REQUIRE(ma.cigar_length==3);
    const cigar_element_t* const e = vector_get_elm(cigar_vector,ma.cigar_offset,cigar_element_t);
    GT_CHECK(e[0].type==cigar_match && e[1].type==cigar_ins && e[1].length==1 && e[2].type==cigar_match);
    GT_CHECK_EQ_INT(e[0].length,left ? 4 : 7);
  }
  vector_delete(cigar_vector);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * Random keys, edited copies in random flanks and random penalties: with a
 * band at least as wide as both sequences the score is the brute-force
 * optimum, and the CIGAR replays over the sequences to that score.
 */
GT_TEST(swg_wide_band_matches_bruteforce) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  vector_t* const cigar_vector = vector_new(1024,cigar_element_t);
  gt_rng_t rng = gt_rng_new(101);
  static uint8_t key[MAX_LEN], text[MAX_LEN];
  uint64_t trial;
  for (trial=0;trial<1500;++trial) {
    penalties_t q = default_penalties;
    if (trial % 3 != 0) {
      q.match = (int32_t)gt_rng_range(&rng,1,3);
      q.mismatch = (int32_t)gt_rng_range(&rng,1,8);
      q.open = (int32_t)gt_rng_range(&rng,0,12);
      q.ext = (int32_t)gt_rng_range(&rng,1,4);
    }
    swg_penalties_t p;
    set_penalties(&p,&q);
    const uint64_t m = gt_rng_range(&rng,2,150);
    random_seq(&rng,key,m,trial%4==0 ? 30 : 0);
    const bool begin_free = gt_rng_chance(&rng,1,2);
    const bool end_free = gt_rng_chance(&rng,1,2);
    const uint64_t left = begin_free ? gt_rng_below(&rng,25) : 0;
    const uint64_t right = end_free ? gt_rng_below(&rng,25) : 0;
    uint64_t n;
    if (trial % 7 == 6) {
      n = gt_rng_range(&rng,2,150);
      random_seq(&rng,text,n,0);
    } else {
      n = mutate(&rng,key,m,gt_rng_below(&rng,m/10+2),left,right,text);
    }
    if (n < 2) continue;
    const uint64_t band = m + n + 2;
    const int32_t best = gotoh(key,(int)m,text,(int)n,&q,begin_free,end_free);
    int left_gap;
    for (left_gap=0;left_gap<2;++left_gap) {
      match_alignment_t ma;
      vector_clear(cigar_vector);
      reset_alignment(&ma,cigar_vector);
      align_swg(key,m,text,n,&p,band,begin_free,end_free,left_gap==1,&ma,cigar_vector,mm_allocator);
      GT_CHECK_MSG(ma.score==best,
          "trial %" PRIu64 " (m=%" PRIu64 " n=%" PRIu64 " free=%d/%d pen=%d,%d,%d,%d left_gap=%d): score %d, brute force %d",
          trial,m,n,begin_free,end_free,q.match,q.mismatch,q.open,q.ext,left_gap,ma.score,best);
      int32_t s; uint64_t edits;
      const int64_t span = replay_cigar(key,m,text,n,ma.match_position,
          cigar_vector,ma.cigar_offset,ma.cigar_length,&q,&s,&edits);
      GT_CHECK_MSG(span>=0 && s==ma.score && span==ma.effective_length,
          "trial %" PRIu64 ": CIGAR replay span %" PRId64 " score %d; reported span %" PRId64 " score %d",
          trial,span,s,ma.effective_length,ma.score);
      GT_CHECK(cigar_is_canonical(cigar_vector,ma.cigar_offset,ma.cigar_length));
      if (!begin_free) GT_CHECK_EQ_UINT(ma.match_position,0);
      if (!end_free) GT_CHECK_MSG(span<0 || ma.match_position+span==n,
          "trial %" PRIu64 ": end not free but alignment ends at %" PRIu64 " of %" PRIu64,trial,ma.match_position+span,n);
    }
  }
  vector_delete(cigar_vector);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * Narrow bands: the reported score is that of the reported CIGAR and never
 * above the optimum; substitution-only copies stay on the diagonal, so any
 * band finds the optimum for them.
 */
GT_TEST(swg_narrow_band) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  vector_t* const cigar_vector = vector_new(1024,cigar_element_t);
  gt_rng_t rng = gt_rng_new(202);
  static uint8_t key[MAX_LEN], text[MAX_LEN];
  uint64_t trial;
  for (trial=0;trial<1500;++trial) {
    const penalties_t q = default_penalties;
    swg_penalties_t p;
    set_penalties(&p,&q);
    const uint64_t m = gt_rng_range(&rng,20,200);
    random_seq(&rng,key,m,0);
    const uint64_t band = gt_rng_range(&rng,1,12);
    match_alignment_t ma;
    if (trial % 2 == 0) {
      // Substitutions only, same length: optimum on the diagonal
      memcpy(text,key,m);
      uint64_t k;
      const uint64_t subs = gt_rng_below(&rng,m/10+1);
      for (k=0;k<subs;++k) {
        const uint64_t pos = gt_rng_below(&rng,m);
        text[pos] = (uint8_t)((text[pos]+1+gt_rng_below(&rng,3))%4);
      }
      const int32_t best = gotoh(key,(int)m,text,(int)m,&q,false,false);
      vector_clear(cigar_vector);
      reset_alignment(&ma,cigar_vector);
      align_swg_banded(key,m,text,m,&p,band,false,false,true,&ma,cigar_vector,mm_allocator);
      GT_CHECK_MSG(ma.score==best,"trial %" PRIu64 " band %" PRIu64 ": score %d, optimum %d",trial,band,ma.score,best);
    } else {
      const uint64_t n = mutate(&rng,key,m,gt_rng_below(&rng,m/15+2),0,0,text);
      if (n + band <= m || n > m + band) continue; // outside the band: no alignment
      const int32_t best = gotoh(key,(int)m,text,(int)n,&q,false,false);
      vector_clear(cigar_vector);
      reset_alignment(&ma,cigar_vector);
      align_swg_banded(key,m,text,n,&p,band,false,false,true,&ma,cigar_vector,mm_allocator);
      GT_CHECK_MSG(ma.score<=best,"trial %" PRIu64 ": banded score %d above the optimum %d",trial,ma.score,best);
      int32_t s; uint64_t edits;
      const int64_t span = replay_cigar(key,m,text,n,ma.match_position,
          cigar_vector,ma.cigar_offset,ma.cigar_length,&q,&s,&edits);
      GT_CHECK_MSG(span==(int64_t)n && s==ma.score,
          "trial %" PRIu64 " (m=%" PRIu64 " n=%" PRIu64 " band=%" PRIu64 "): replay span %" PRId64 " score %d, reported %d",
          trial,m,n,band,span,s,ma.score);
    }
  }
  // Outside the band in both directions: no alignment
  {
    swg_penalties_t p;
    set_penalties(&p,&default_penalties);
    random_seq(&rng,key,50,0);
    random_seq(&rng,text,80,0);
    match_alignment_t ma;
    reset_alignment(&ma,cigar_vector);
    align_swg_banded(key,50,text,80,&p,10,false,false,true,&ma,cigar_vector,mm_allocator);
    GT_CHECK_EQ_INT(ma.score,SWG_SCORE_MIN);
    reset_alignment(&ma,cigar_vector);
    align_swg_banded(text,80,key,50,&p,10,false,false,true,&ma,cigar_vector,mm_allocator);
    GT_CHECK_EQ_INT(ma.score,SWG_SCORE_MIN);
  }
  vector_delete(cigar_vector);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * align_swg_base (the full DP --check-alignments uses): both text ends free
 */
GT_TEST(swg_base_matches_bruteforce) {
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  vector_t* const cigar_vector = vector_new(1024,cigar_element_t);
  gt_rng_t rng = gt_rng_new(303);
  static uint8_t key[MAX_LEN], text[MAX_LEN];
  uint64_t trial;
  for (trial=0;trial<600;++trial) {
    const penalties_t q = default_penalties;
    swg_penalties_t p;
    set_penalties(&p,&q);
    const uint64_t m = gt_rng_range(&rng,5,120);
    random_seq(&rng,key,m,trial%5==0 ? 30 : 0);
    const uint64_t n = mutate(&rng,key,m,gt_rng_below(&rng,m/10+2),
        gt_rng_below(&rng,30),gt_rng_below(&rng,30),text);
    if (n < 1) continue;
    const int32_t best = gotoh(key,(int)m,text,(int)n,&q,true,true);
    match_alignment_t ma;
    vector_clear(cigar_vector);
    reset_alignment(&ma,cigar_vector);
    align_swg_base(key,m,text,n,&p,true,&ma,cigar_vector,mm_allocator);
    GT_CHECK_MSG(ma.score==best,"trial %" PRIu64 ": score %d, brute force %d",trial,ma.score,best);
    int32_t s; uint64_t edits;
    const int64_t span = replay_cigar(key,m,text,n,ma.match_position,
        cigar_vector,ma.cigar_offset,ma.cigar_length,&q,&s,&edits);
    GT_CHECK_MSG(span==ma.effective_length && s==best,
        "trial %" PRIu64 ": replay span %" PRId64 " score %d",trial,span,s);
  }
  vector_delete(cigar_vector);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

GT_CASES(
  GT_CASE(swg_score_helpers_hand),
  GT_CASE(search_parameters_defaults_and_instantiation),
  GT_CASE(swg_hand_examples),
  GT_CASE(swg_trivial_shapes),
  GT_CASE(swg_left_gap_alignment),
  GT_CASE(swg_wide_band_matches_bruteforce),
  GT_CASE(swg_narrow_band),
  GT_CASE(swg_base_matches_bruteforce)
)
