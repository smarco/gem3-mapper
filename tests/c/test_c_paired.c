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
 *   Paired-end matches: pair orientation, layout and template length,
 *   the pair relation (concordant, discordant, invalid) under the
 *   orientation/layout tables, -l/-L and the estimated template-length
 *   distribution, pairing all matches of both ends, discordant search
 *   modes, the order of paired maps, paired classification, paired
 *   selection (-M) and the paired-end MAPQ models.
 *
 *   Expected template lengths are the SAM span (leftmost to rightmost
 *   mapped base). Two known deviations are not tested here (known-bugs/
 *   C-tlen-sam-spec.md, C-pair-layout-end2-contains-end1.md): RF pairs that
 *   overlap get the overlap length, and a pair whose end 2 contains end 1 is
 *   classified 'overlap' instead of 'contain'.
 *
 *   Internal signatures called (a speed iteration may change them; update
 *   this driver with them):
 *     paired_map_compute_{event_distance,edit_distance,swg_score,error_quality}(const match_trace_t*,const match_trace_t*)
 *     paired_matches_{new,clear,delete,is_mapped,get_max_complete_stratum,get_primary_map,
 *         get_num_maps,get_maps,get_num_discordant_maps,get_discordant_maps,limit_maps,
 *         recompute_metrics,add,add_map,cross_pair,find_pairs,find_discordant_pairs}
 *     (no prototype in a header:) paired_matches_cmp_swg_score(const paired_map_t*,const paired_map_t*),
 *         paired_matches_add_map_insert_sorted, paired_matches_compute_{orientation,layout,
 *         template_length,relation}, paired_matches_find_pairs_{locate,sorted},
 *         matches_add_match_trace_insert_sorted, archive_score_matches_pe_{default,classify},
 *         archive_search_pe_use_recovery_extension
 *     paired_matches_classify, paired_matches_classify_{search_accomplished,subdominant_end}
 *     matches_predictors_compute_pe(matches_predictors_t*,paired_matches_t*)
 *     archive_select_pe_matches(select_parameters_t*,mapper_stats_t*,paired_matches_t*)
 *     archive_score_matches_pe(archive_search_t*,archive_search_t*,paired_matches_t*)
 *     archive_search_pe_use_shortcut_extension(archive_search_t*,archive_search_t*,matches_t*)
 *     search_paired_parameters_init(search_paired_parameters_t*)
 *     mapper_stats_{new,delete,template_init,template_length_sample,template_length_is_reliable,
 *         template_length_get_mean,template_length_get_stddev,template_length_get_expected_{min,max},
 *         template_length_get_num_samples} (area D; used to set up the distribution)
 */

#include "gem3_test.h"
#include "matches/matches.h"
#include "matches/paired_map.h"
#include "matches/paired_matches.h"
#include "matches/classify/matches_classify.h"
#include "matches/classify/paired_matches_classify.h"
#include "matches/classify/matches_predictors.h"
#include "archive/search/archive_search.h"
#include "archive/search/archive_search_pe_parameters.h"
#include "archive/search/archive_search_pe_stages.h"
#include "archive/search/archive_select.h"
#include "archive/score/archive_score_pe.h"
#include "mapper/mapper_stats.h"
#include <math.h>

/* Defined in src/ without a prototype in a header */
int paired_matches_cmp_swg_score(const paired_map_t* const a,const paired_map_t* const b);
void paired_matches_add_map_insert_sorted(paired_matches_t* const paired_matches,paired_map_t* const paired_map);
pair_orientation_t paired_matches_compute_orientation(
    const match_trace_t* const match_trace_end1,const match_trace_t* const match_trace_end2);
pair_layout_t paired_matches_compute_layout(
    const match_trace_t* const match_trace_end1,const match_trace_t* const match_trace_end2);
uint64_t paired_matches_compute_template_length(
    const match_trace_t* const match_trace_end1,const match_trace_t* const match_trace_end2,
    const pair_orientation_t pair_orientation,const pair_layout_t pair_layout);
pair_relation_t paired_matches_compute_relation(
    const search_paired_parameters_t* const parameters,mapper_stats_t* const mapper_stats,
    match_trace_t* const match_trace_end1,match_trace_t* const match_trace_end2,
    pair_orientation_t* const pair_orientation,pair_layout_t* const pair_layout,
    uint64_t* const template_length,double* const template_length_sigma);
uint64_t paired_matches_find_pairs_locate(
    match_trace_t** const match_trace_sorted,const uint64_t num_match_traces,const char* const sequence_name);
void matches_add_match_trace_insert_sorted(matches_t* const matches,match_trace_t* const match_trace);
void archive_score_matches_pe_default(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);
void archive_score_matches_pe_classify(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);
bool archive_search_pe_use_recovery_extension(
    archive_search_t* const archive_search_end1,archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches);

/*
 * Match traces
 */
static char name_a[] = "chrA", name_b[] = "chrB";
static void make_trace(
    match_trace_t* const t,char* const name,const uint64_t position,const uint64_t length,
    const strand_t strand,const int32_t swg,const uint64_t edit) {
  memset(t,0,sizeof(*t));
  t->type = match_type_regular;
  t->sequence_name = name;
  t->text_position = position;
  t->strand = strand;
  t->bs_strand = bs_strand_none;
  t->swg_score = swg;
  t->edit_distance = edit;
  t->event_distance = edit;
  t->error_quality = 30.0f;
  t->match_alignment.effective_length = (int64_t)length;
  t->match_alignment.match_position = position + (name==name_b ? 1000000 : 0);
}
static uint64_t t_begin(const match_trace_t* const t) { return t->text_position; }
static uint64_t t_end(const match_trace_t* const t) { return t->text_position + (uint64_t)t->match_alignment.effective_length; }
static uint64_t ref_span(const match_trace_t* const a,const match_trace_t* const b) {
  return MAX(t_end(a),t_end(b)) - MIN(t_begin(a),t_begin(b));
}
// FF/RR by strands; otherwise FR when the forward end starts at or before the reverse end
static pair_orientation_t ref_orientation(const match_trace_t* const a,const match_trace_t* const b) {
  if (a->strand==Forward && b->strand==Forward) return pair_orientation_FF;
  if (a->strand==Reverse && b->strand==Reverse) return pair_orientation_RR;
  const match_trace_t* const fwd = (a->strand==Forward) ? a : b;
  const match_trace_t* const rev = (a->strand==Forward) ? b : a;
  return (fwd->text_position <= rev->text_position) ? pair_orientation_FR : pair_orientation_RF;
}
// Separate with a gap of at least one base; contain when end 1 covers end 2;
// otherwise overlap (mates that touch count as overlapping)
static pair_layout_t ref_layout(const match_trace_t* const a,const match_trace_t* const b) {
  if (t_end(a) < t_begin(b) || t_end(b) < t_begin(a)) return pair_layout_separate;
  if (t_begin(a) <= t_begin(b) && t_end(b) <= t_end(a)) return pair_layout_contain;
  return pair_layout_overlap;
}
static bool end2_contains_end1_only(const match_trace_t* const a,const match_trace_t* const b) {
  const bool b_in_a = t_begin(a) <= t_begin(b) && t_end(b) <= t_end(a);
  const bool a_in_b = t_begin(b) <= t_begin(a) && t_end(a) <= t_end(b);
  return a_in_b && !b_in_a;
}

GT_TEST(paired_map_scores) {
  match_trace_t a, b;
  make_trace(&a,name_a,100,100,Forward,92,3);
  make_trace(&b,name_a,400,100,Reverse,85,5);
  a.event_distance = 2; b.event_distance = 4;
  a.error_quality = 30.0f; b.error_quality = 21.0f;
  GT_CHECK_EQ_UINT(paired_map_compute_edit_distance(&a,&b),8);
  GT_CHECK_EQ_UINT(paired_map_compute_event_distance(&a,&b),6);
  GT_CHECK_EQ_INT(paired_map_compute_swg_score(&a,&b),177);
  GT_CHECK(fabsf(paired_map_compute_error_quality(&a,&b)-25.5f) < 1e-6f);
}

GT_TEST(orientation_layout_template_length) {
  gt_rng_t rng = gt_rng_new(201);
  uint64_t trial, checked = 0;
  for (trial=0;trial<20000;++trial) {
    match_trace_t a, b;
    make_trace(&a,name_a,gt_rng_below(&rng,400),gt_rng_range(&rng,1,150),
        gt_rng_below(&rng,2) ? Forward : Reverse,0,0);
    make_trace(&b,name_a,gt_rng_below(&rng,400),gt_rng_range(&rng,1,150),
        gt_rng_below(&rng,2) ? Forward : Reverse,0,0);
    const pair_orientation_t o = paired_matches_compute_orientation(&a,&b);
    GT_CHECK_EQ_INT(o,ref_orientation(&a,&b));
    if (end2_contains_end1_only(&a,&b)) continue; // known bug: classified overlap
    const pair_layout_t l = paired_matches_compute_layout(&a,&b);
    GT_CHECK_MSG(l==ref_layout(&a,&b),"[%" PRIu64 ",%" PRIu64 ") [%" PRIu64 ",%" PRIu64 "): layout %d",
        t_begin(&a),t_end(&a),t_begin(&b),t_end(&b),(int)l);
    if (o==pair_orientation_RF && l==pair_layout_overlap) continue; // known bug: overlap length
    const uint64_t tl = paired_matches_compute_template_length(&a,&b,o,l);
    GT_CHECK_MSG(tl==ref_span(&a,&b),"[%" PRIu64 ",%" PRIu64 ")%c [%" PRIu64 ",%" PRIu64 ")%c: TL %" PRIu64 " span %" PRIu64,
        t_begin(&a),t_end(&a),a.strand==Forward?'+':'-',t_begin(&b),t_end(&b),b.strand==Forward?'+':'-',tl,ref_span(&a,&b));
    ++checked;
  }
  GT_CHECK(checked > 10000);
  // Hand-worked: FR 300 bp fragment of two 100 bp mates, either order
  match_trace_t a, b;
  make_trace(&a,name_a,1000,100,Forward,0,0);
  make_trace(&b,name_a,1200,100,Reverse,0,0);
  GT_CHECK_EQ_INT(paired_matches_compute_orientation(&a,&b),pair_orientation_FR);
  GT_CHECK_EQ_INT(paired_matches_compute_orientation(&b,&a),pair_orientation_FR);
  GT_CHECK_EQ_INT(paired_matches_compute_layout(&a,&b),pair_layout_separate);
  GT_CHECK_EQ_UINT(paired_matches_compute_template_length(&a,&b,pair_orientation_FR,pair_layout_separate),300);
  GT_CHECK_EQ_UINT(paired_matches_compute_template_length(&b,&a,pair_orientation_FR,pair_layout_separate),300);
  // Mate-pair (RF): reverse mate on the left
  a.strand = Reverse; b.strand = Forward;
  GT_CHECK_EQ_INT(paired_matches_compute_orientation(&a,&b),pair_orientation_RF);
  GT_CHECK_EQ_UINT(paired_matches_compute_template_length(&a,&b,pair_orientation_RF,pair_layout_separate),300);
  // Mates sharing the same 100 bp: contain, span 100
  make_trace(&b,name_a,1000,100,Forward,0,0);
  GT_CHECK_EQ_INT(paired_matches_compute_layout(&a,&b),pair_layout_contain);
  GT_CHECK_EQ_UINT(paired_matches_compute_template_length(&a,&b,pair_orientation_FR,pair_layout_contain),100);
}

/*
 * Pair relation under the parameter tables, -l/-L and the template-length
 * distribution
 */
static mapper_stats_t* stats_with_prior(void) {
  mapper_stats_t* const stats = mapper_stats_new();
  mapper_stats_template_init(stats,100,0,1000); // gem-mapper's prior (template_length_estimation_*)
  return stats;
}
typedef struct {
  pair_relation_t relation;
  uint64_t tl;
  double sigma;
} ref_relation_t;
static ref_relation_t ref_relation(
    const search_paired_parameters_t* const p,mapper_stats_t* const stats,
    const match_trace_t* const a,const match_trace_t* const b) {
  ref_relation_t r = { pair_relation_concordant, 0, 0.0 };
  if (a->bs_strand!=b->bs_strand) { r.relation = pair_relation_invalid; return r; }
  const pair_relation_t ro = p->pair_orientation[ref_orientation(a,b)];
  const pair_relation_t rl = p->pair_layout[ref_layout(a,b)];
  if (ro==pair_relation_invalid || rl==pair_relation_invalid) { r.relation = pair_relation_invalid; return r; }
  if (ro==pair_relation_discordant || rl==pair_relation_discordant) r.relation = pair_relation_discordant;
  r.tl = ref_span(a,b);
  if (r.tl < p->min_template_length || r.tl > p->max_template_length) { r.relation = pair_relation_invalid; return r; }
  if (mapper_stats_template_length_is_reliable(stats)) {
    const double mean = mapper_stats_template_length_get_mean(stats);
    const double sd = mapper_stats_template_length_get_stddev(stats);
    r.sigma = fabs((double)r.tl-mean)/sd;
    if (r.tl > mapper_stats_template_length_get_expected_max(stats) ||
        r.tl < mapper_stats_template_length_get_expected_min(stats)) r.relation = pair_relation_discordant;
  }
  return r;
}

GT_TEST(search_paired_parameters_defaults) {
  search_paired_parameters_t p;
  memset(&p,0xff,sizeof(p));
  search_paired_parameters_init(&p);
  GT_CHECK(!p.paired_end_search);
  GT_CHECK(!p.paired_end_extension_shortcut);
  GT_CHECK(p.paired_end_extension_recovery);
  GT_CHECK_EQ_INT(p.pair_discordant_search,pair_discordant_search_only_if_no_concordant);
  // --help: --pair-orientation FR, --discordant-pair-orientation RF,
  // --pair-layout separate,overlap, --discordant-pair-layout contain
  GT_CHECK_EQ_INT(p.pair_orientation[pair_orientation_FR],pair_relation_concordant);
  GT_CHECK_EQ_INT(p.pair_orientation[pair_orientation_RF],pair_relation_discordant);
  GT_CHECK_EQ_INT(p.pair_orientation[pair_orientation_FF],pair_relation_invalid);
  GT_CHECK_EQ_INT(p.pair_orientation[pair_orientation_RR],pair_relation_invalid);
  GT_CHECK_EQ_INT(p.pair_layout[pair_layout_separate],pair_relation_concordant);
  GT_CHECK_EQ_INT(p.pair_layout[pair_layout_overlap],pair_relation_concordant);
  GT_CHECK_EQ_INT(p.pair_layout[pair_layout_contain],pair_relation_discordant);
  // --min-template-length disabled, --max-template-length 10000
  GT_CHECK_EQ_UINT(p.min_template_length,0);
  GT_CHECK_EQ_UINT(p.max_template_length,10000);
  GT_CHECK_EQ_UINT(p.template_length_estimation_samples,100);
}

GT_TEST(relation_tables_and_limits) {
  gt_rng_t rng = gt_rng_new(202);
  search_paired_parameters_t p;
  uint64_t trial;
  mapper_stats_t* const reliable = stats_with_prior();
  mapper_stats_t* const unreliable = mapper_stats_new();
  GT_REQUIRE(mapper_stats_template_length_is_reliable(reliable));
  GT_REQUIRE(!mapper_stats_template_length_is_reliable(unreliable));
  for (trial=0;trial<20000;++trial) {
    search_paired_parameters_init(&p);
    uint64_t k;
    if (trial%2) {
      for (k=0;k<4;++k) p.pair_orientation[k] = (pair_relation_t)gt_rng_below(&rng,3);
      for (k=0;k<3;++k) p.pair_layout[k] = (pair_relation_t)gt_rng_below(&rng,3);
    }
    p.min_template_length = gt_rng_chance(&rng,1,2) ? 0 : gt_rng_below(&rng,800);
    p.max_template_length = gt_rng_chance(&rng,1,2) ? 10000 : gt_rng_range(&rng,100,3000);
    mapper_stats_t* const stats = gt_rng_chance(&rng,1,4) ? unreliable : reliable;
    match_trace_t a, b;
    make_trace(&a,name_a,gt_rng_below(&rng,2500),gt_rng_range(&rng,50,150),gt_rng_below(&rng,2) ? Forward : Reverse,0,0);
    make_trace(&b,name_a,gt_rng_below(&rng,2500),gt_rng_range(&rng,50,150),gt_rng_below(&rng,2) ? Forward : Reverse,0,0);
    if (gt_rng_chance(&rng,1,10)) b.bs_strand = bs_strand_G2A;
    if (end2_contains_end1_only(&a,&b)) continue;
    if (ref_orientation(&a,&b)==pair_orientation_RF && ref_layout(&a,&b)==pair_layout_overlap) continue;
    pair_orientation_t o;
    pair_layout_t l;
    uint64_t tl = 0;
    double sigma = -1.0;
    const pair_relation_t rel = paired_matches_compute_relation(&p,stats,&a,&b,&o,&l,&tl,&sigma);
    const ref_relation_t want = ref_relation(&p,stats,&a,&b);
    GT_CHECK_MSG(rel==want.relation,"trial %" PRIu64 ": relation %d expected %d (TL %" PRIu64 ")",trial,(int)rel,(int)want.relation,want.tl);
    if (rel!=pair_relation_invalid) {
      GT_CHECK_EQ_INT(o,ref_orientation(&a,&b));
      GT_CHECK_EQ_INT(l,ref_layout(&a,&b));
      GT_CHECK_EQ_UINT(tl,want.tl);
      GT_CHECK_MSG(fabs(sigma-want.sigma) < 1e-9,"trial %" PRIu64 ": sigma %f expected %f",trial,sigma,want.sigma);
    } else if (want.tl!=0 && a.bs_strand==b.bs_strand) {
      // Outside -l/-L: the deviation is set to the maximum
      GT_CHECK(sigma==MAX_TEMPLATE_LENGTH_SIGMAS);
    }
  }
  // Hand-worked limits with the default tables: FR, 300 bp
  search_paired_parameters_init(&p);
  match_trace_t a, b;
  make_trace(&a,name_a,1000,100,Forward,0,0);
  make_trace(&b,name_a,1200,100,Reverse,0,0);
  pair_orientation_t o; pair_layout_t l; uint64_t tl; double sigma;
  p.max_template_length = 299;
  GT_CHECK_EQ_INT(paired_matches_compute_relation(&p,reliable,&a,&b,&o,&l,&tl,&sigma),pair_relation_invalid);
  p.max_template_length = 300;
  GT_CHECK_EQ_INT(paired_matches_compute_relation(&p,reliable,&a,&b,&o,&l,&tl,&sigma),pair_relation_concordant);
  p.min_template_length = 301;
  GT_CHECK_EQ_INT(paired_matches_compute_relation(&p,reliable,&a,&b,&o,&l,&tl,&sigma),pair_relation_invalid);
  p.min_template_length = 300;
  GT_CHECK_EQ_INT(paired_matches_compute_relation(&p,reliable,&a,&b,&o,&l,&tl,&sigma),pair_relation_concordant);
  // Beyond the prior's expected range (mean 500 + 6 sd) but within -L: discordant
  search_paired_parameters_init(&p);
  make_trace(&b,name_a,2900,100,Reverse,0,0);
  GT_CHECK_EQ_INT(paired_matches_compute_relation(&p,reliable,&a,&b,&o,&l,&tl,&sigma),pair_relation_discordant);
  GT_CHECK_EQ_UINT(tl,2000);
  GT_CHECK_EQ_INT(paired_matches_compute_relation(&p,unreliable,&a,&b,&o,&l,&tl,&sigma),pair_relation_concordant);
  GT_CHECK(sigma==0.0);
  // Different bisulfite strands never pair
  b.bs_strand = bs_strand_C2T;
  GT_CHECK_EQ_INT(paired_matches_compute_relation(&p,unreliable,&a,&b,&o,&l,&tl,&sigma),pair_relation_invalid);
  mapper_stats_delete(reliable);
  mapper_stats_delete(unreliable);
}

/*
 * Pair all matches of both ends: every same-sequence combination is
 * classified; concordant pairs are kept in order, discordant ones apart
 */
static paired_matches_t* new_paired_matches(void) { return paired_matches_new(); }
static void add_end_matches(matches_t* const matches,match_trace_t* const traces,const uint64_t n) {
  uint64_t i;
  for (i=0;i<n;++i) {
    matches_add_match_trace_insert_sorted(matches,traces+i);
    matches_counters_add(matches->counters,traces[i].edit_distance,1);
    matches_metrics_update(&matches->metrics,traces[i].event_distance,traces[i].edit_distance,traces[i].swg_score);
  }
}
static bool map_is(const paired_map_t* const m,const match_trace_t* const a,const match_trace_t* const b) {
  return m->match_trace_end1==a && m->match_trace_end2==b;
}

GT_TEST(find_pairs_all_combinations) {
  gt_rng_t rng = gt_rng_new(203);
  static match_trace_t e1[8], e2[8];
  mm_slab_t* const slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM);
  mm_allocator_t* const mm_allocator = mm_allocator_new(slab);
  mapper_stats_t* const stats = stats_with_prior();
  search_parameters_t sp;
  search_parameters_init(&sp);
  uint64_t trial, i, j;
  for (trial=0;trial<3000;++trial) {
    search_paired_parameters_init(&sp.search_paired_parameters);
    sp.search_paired_parameters.pair_discordant_search = (pair_discordant_search_t)gt_rng_below(&rng,3);
    sp.search_paired_parameters.max_template_length = (trial%3==0) ? gt_rng_range(&rng,200,3000) : 10000;
    sp.search_paired_parameters.min_template_length = (trial%5==0) ? gt_rng_below(&rng,800) : 0;
    // Mates 100 bp at multiples of 500 on two sequences: separate, or identical
    // (contain). Both ends have matches, as when gem-mapper pairs them.
    const uint64_t n1 = gt_rng_range(&rng,1,6), n2 = gt_rng_range(&rng,1,6);
    for (i=0;i<n1;++i) {
      make_trace(e1+i,gt_rng_below(&rng,4) ? name_a : name_b,500*gt_rng_below(&rng,10),100,
          gt_rng_below(&rng,2) ? Forward : Reverse,(int32_t)gt_rng_range(&rng,70,100),gt_rng_below(&rng,6));
      if (gt_rng_chance(&rng,1,10)) e1[i].bs_strand = bs_strand_C2T;
    }
    for (i=0;i<n2;++i) {
      make_trace(e2+i,gt_rng_below(&rng,4) ? name_a : name_b,500*gt_rng_below(&rng,10),100,
          gt_rng_below(&rng,2) ? Forward : Reverse,(int32_t)gt_rng_range(&rng,70,100),gt_rng_below(&rng,6));
      if (gt_rng_chance(&rng,1,10)) e2[i].bs_strand = bs_strand_C2T;
    }
    paired_matches_t* const pm = new_paired_matches();
    add_end_matches(pm->matches_end1,e1,n1);
    add_end_matches(pm->matches_end2,e2,n2);
    paired_matches_find_pairs(pm,&sp,stats,mm_allocator);
    // Expected
    uint64_t want_conc = 0, want_disc = 0;
    for (i=0;i<n1;++i) {
      for (j=0;j<n2;++j) {
        if (strcmp(e1[i].sequence_name,e2[j].sequence_name)!=0) continue;
        const ref_relation_t r = ref_relation(&sp.search_paired_parameters,stats,e1+i,e2+j);
        uint64_t found = 0, k;
        if (r.relation==pair_relation_concordant) {
          ++want_conc;
          for (k=0;k<paired_matches_get_num_maps(pm);++k) {
            paired_map_t* const m = paired_matches_get_maps(pm)[k];
            if (!map_is(m,e1+i,e2+j)) continue;
            ++found;
            GT_CHECK_EQ_UINT(m->template_length,r.tl);
            GT_CHECK(fabs(m->template_length_sigma-r.sigma) < 1e-9);
            GT_CHECK_EQ_INT(m->pair_relation,pair_relation_concordant);
            GT_CHECK_EQ_INT(m->swg_score,e1[i].swg_score+e2[j].swg_score);
            GT_CHECK_EQ_UINT(m->edit_distance,e1[i].edit_distance+e2[j].edit_distance);
            GT_CHECK_EQ_UINT(m->index_position,MIN(e1[i].match_alignment.match_position,e2[j].match_alignment.match_position));
          }
          GT_CHECK_MSG(found==1,"trial %" PRIu64 ": concordant pair (%" PRIu64 ",%" PRIu64 ") found %" PRIu64 " times",trial,i,j,found);
        } else if (r.relation==pair_relation_discordant &&
                   sp.search_paired_parameters.pair_discordant_search!=pair_discordant_search_never) {
          ++want_disc;
          for (k=0;k<paired_matches_get_num_discordant_maps(pm);++k) {
            found += map_is(paired_matches_get_discordant_maps(pm)[k],e1+i,e2+j);
          }
          GT_CHECK_MSG(found==1,"trial %" PRIu64 ": discordant pair (%" PRIu64 ",%" PRIu64 ") found %" PRIu64 " times",trial,i,j,found);
        }
      }
    }
    GT_CHECK_EQ_UINT(paired_matches_get_num_maps(pm),want_conc);
    GT_CHECK_EQ_UINT(paired_matches_get_num_discordant_maps(pm),want_disc);
    GT_CHECK_EQ_UINT(paired_matches_counters_get_total_count(pm),want_conc);
    // Best pairs first: SWG score high to low
    paired_map_t** const maps = paired_matches_get_maps(pm);
    for (i=0;i+1<paired_matches_get_num_maps(pm);++i) {
      GT_CHECK(maps[i]->swg_score >= maps[i+1]->swg_score);
      GT_CHECK(paired_matches_cmp_swg_score(maps[i],maps[i+1]) <= 0);
    }
    // Discordant search: appended after the concordant pairs when asked to
    const uint64_t num_conc = paired_matches_get_num_maps(pm);
    paired_matches_find_discordant_pairs(pm,&sp);
    const pair_discordant_search_t mode = sp.search_paired_parameters.pair_discordant_search;
    const bool added = want_disc>0 && (mode==pair_discordant_search_always ||
        (mode==pair_discordant_search_only_if_no_concordant && num_conc==0));
    GT_CHECK_EQ_UINT(paired_matches_get_num_maps(pm),num_conc+(added ? want_disc : 0));
    GT_CHECK_EQ_UINT(paired_matches_get_num_discordant_maps(pm),added ? 0 : want_disc);
    for (i=0;i<paired_matches_get_num_maps(pm);++i) {
      GT_CHECK_EQ_INT(paired_matches_get_maps(pm)[i]->pair_relation,
          i<num_conc ? pair_relation_concordant : pair_relation_discordant);
    }
    GT_CHECK_EQ_INT(paired_matches_is_mapped(pm),paired_matches_get_num_maps(pm)>0);
    paired_matches_delete(pm);
  }
  // Locate: first trace of a sequence in a genomic-position-sorted list
  match_trace_t s[4];
  match_trace_t* sorted[4];
  make_trace(s+0,name_a,10,100,Forward,0,0);
  make_trace(s+1,name_a,50,100,Forward,0,0);
  make_trace(s+2,name_b,5,100,Forward,0,0);
  make_trace(s+3,name_b,7,100,Forward,0,0);
  for (i=0;i<4;++i) sorted[i] = s+i;
  GT_CHECK_EQ_UINT(paired_matches_find_pairs_locate(sorted,4,"chrA"),0);
  GT_CHECK_EQ_UINT(paired_matches_find_pairs_locate(sorted,4,"chrB"),2);
  GT_CHECK_EQ_UINT(paired_matches_find_pairs_locate(sorted,4,"chrC"),4);
  mapper_stats_delete(stats);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * Ordering of paired maps
 */
static paired_map_t make_map(
    const pair_relation_t rel,const int32_t swg,const double sigma,const uint64_t tl,
    const pair_layout_t layout,const uint64_t index_position) {
  paired_map_t m;
  memset(&m,0,sizeof(m));
  m.pair_relation = rel;
  m.swg_score = swg;
  m.template_length_sigma = sigma;
  m.template_length = tl;
  m.pair_layout = layout;
  m.index_position = index_position;
  return m;
}
GT_TEST(paired_map_order) {
  paired_map_t a, b;
  // Concordant before discordant, whatever the score
  a = make_map(pair_relation_concordant,100,0.5,300,pair_layout_separate,5);
  b = make_map(pair_relation_discordant,190,0.5,300,pair_layout_separate,1);
  GT_CHECK(paired_matches_cmp_swg_score(&a,&b) < 0);
  GT_CHECK(paired_matches_cmp_swg_score(&b,&a) > 0);
  // Higher score first
  b = make_map(pair_relation_concordant,190,0.5,300,pair_layout_separate,1);
  GT_CHECK(paired_matches_cmp_swg_score(&b,&a) < 0);
  // Same score: fewer whole sigmas of template-length deviation first
  a = make_map(pair_relation_concordant,180,2.9,300,pair_layout_separate,5);
  b = make_map(pair_relation_concordant,180,1.2,300,pair_layout_separate,9);
  GT_CHECK(paired_matches_cmp_swg_score(&b,&a) < 0);
  // ... equal whole sigmas, same layout: index position
  a = make_map(pair_relation_concordant,180,1.9,300,pair_layout_separate,5);
  b = make_map(pair_relation_concordant,180,1.2,320,pair_layout_separate,9);
  GT_CHECK(paired_matches_cmp_swg_score(&a,&b) < 0);
  GT_CHECK(paired_matches_cmp_swg_score(&a,&a)==0);
  // Unknown deviation (sigma 0, distribution not reliable): shorter template first
  a = make_map(pair_relation_concordant,180,0.0,400,pair_layout_separate,1);
  b = make_map(pair_relation_concordant,180,0.0,350,pair_layout_separate,9);
  GT_CHECK(paired_matches_cmp_swg_score(&b,&a) < 0);
  // Random maps inserted in any order come out sorted
  gt_rng_t rng = gt_rng_new(204);
  static paired_map_t maps[64];
  uint64_t trial, i;
  for (trial=0;trial<200;++trial) {
    paired_matches_t* const pm = new_paired_matches();
    const uint64_t n = gt_rng_range(&rng,1,40);
    for (i=0;i<n;++i) {
      maps[i] = make_map((pair_relation_t)gt_rng_below(&rng,2),(int32_t)gt_rng_range(&rng,150,200),
          (10+gt_rng_below(&rng,30))/10.0,gt_rng_range(&rng,200,400),pair_layout_separate,gt_rng_below(&rng,1000));
      maps[i].edit_distance = gt_rng_below(&rng,5);
      maps[i].event_distance = maps[i].edit_distance;
      if (trial%2) paired_matches_add_map_insert_sorted(pm,maps+i);
      else paired_matches_add_map(pm,maps+i);
    }
    GT_REQUIRE(paired_matches_get_num_maps(pm)==n);
    paired_map_t** const pmaps = paired_matches_get_maps(pm);
    for (i=0;i+1<n;++i) {
      GT_CHECK(paired_matches_cmp_swg_score(pmaps[i],pmaps[i+1]) <= 0);
      GT_CHECK(pmaps[i]->pair_relation==pmaps[i+1]->pair_relation || pmaps[i]->pair_relation==pair_relation_concordant);
    }
    GT_CHECK(paired_matches_get_primary_map(pm)==pmaps[0]);
    if (trial%2==0) {
      // add_map also counts and updates the metrics
      matches_metrics_t expected;
      matches_metrics_init(&expected);
      for (i=0;i<n;++i) matches_metrics_update(&expected,maps[i].event_distance,maps[i].edit_distance,maps[i].swg_score);
      GT_CHECK_EQ_UINT(paired_matches_counters_get_total_count(pm),n);
      GT_CHECK_EQ_INT(pm->metrics.max1_swg_score,expected.max1_swg_score);
      GT_CHECK_EQ_INT(pm->metrics.max2_swg_score,expected.max2_swg_score);
      GT_CHECK_EQ_UINT(pm->metrics.min1_edit_distance,expected.min1_edit_distance);
      GT_CHECK_EQ_UINT(pm->metrics.min2_edit_distance,expected.min2_edit_distance);
      // Limit keeps the best; recompute follows the kept maps
      const uint64_t keep = gt_rng_range(&rng,1,n);
      paired_matches_limit_maps(pm,keep);
      GT_CHECK_EQ_UINT(paired_matches_get_num_maps(pm),keep);
      paired_matches_recompute_metrics(pm);
      matches_metrics_init(&expected);
      for (i=0;i<keep;++i) matches_metrics_update(&expected,pmaps[i]->event_distance,pmaps[i]->edit_distance,pmaps[i]->swg_score);
      GT_CHECK_EQ_INT(pm->metrics.max1_swg_score,expected.max1_swg_score);
      GT_CHECK_EQ_INT(pm->metrics.max2_swg_score,expected.max2_swg_score);
      GT_CHECK_EQ_UINT(pm->metrics.min1_edit_distance,expected.min1_edit_distance);
      GT_CHECK_EQ_UINT(pm->metrics.min2_edit_distance,expected.min2_edit_distance);
    }
    paired_matches_clear(pm,true);
    GT_CHECK_EQ_UINT(paired_matches_get_num_maps(pm),0);
    GT_CHECK(!paired_matches_is_mapped(pm));
    paired_matches_delete(pm);
  }
}

GT_TEST(add_and_cross_pair) {
  paired_matches_t* const pm = new_paired_matches();
  search_parameters_t sp;
  search_parameters_init(&sp);
  mapper_stats_t* const stats = mapper_stats_new(); // not reliable: sigma 0
  match_trace_t a, b, c;
  make_trace(&a,name_a,1000,100,Forward,95,1);
  make_trace(&b,name_a,1250,100,Reverse,90,2);
  make_trace(&c,name_a,1250,100,Forward,90,2); // FF with a: invalid
  a.match_alignment.match_position = 77;
  b.match_alignment.match_position = 55;
  paired_matches_cross_pair(pm,&sp,stats,&a,&b);
  paired_matches_cross_pair(pm,&sp,stats,&a,&c);
  GT_REQUIRE(paired_matches_get_num_maps(pm)==1);
  GT_CHECK_EQ_UINT(paired_matches_get_num_discordant_maps(pm),0);
  paired_map_t* const m = paired_matches_get_primary_map(pm);
  GT_CHECK(m->match_trace_end1==&a && m->match_trace_end2==&b);
  GT_CHECK_EQ_INT(m->pair_orientation,pair_orientation_FR);
  GT_CHECK_EQ_INT(m->pair_layout,pair_layout_separate);
  GT_CHECK_EQ_UINT(m->template_length,350);
  GT_CHECK(m->template_length_sigma==0.0);
  GT_CHECK_EQ_UINT(m->index_position,55);
  GT_CHECK_EQ_INT(m->swg_score,185);
  GT_CHECK_EQ_UINT(m->edit_distance,3);
  GT_CHECK_EQ_UINT(m->mapq_score,0);
  // RF: discordant, kept apart (unless 'never')
  b.strand = Reverse; a.strand = Forward;
  match_trace_t d;
  make_trace(&d,name_a,700,100,Reverse,90,2); // reverse end left of the forward one: RF
  paired_matches_cross_pair(pm,&sp,stats,&a,&d);
  GT_CHECK_EQ_UINT(paired_matches_get_num_maps(pm),1);
  GT_REQUIRE(paired_matches_get_num_discordant_maps(pm)==1);
  GT_CHECK_EQ_INT(paired_matches_get_discordant_maps(pm)[0]->pair_orientation,pair_orientation_RF);
  GT_CHECK_EQ_UINT(paired_matches_get_discordant_maps(pm)[0]->template_length,400);
  sp.search_paired_parameters.pair_discordant_search = pair_discordant_search_never;
  paired_matches_cross_pair(pm,&sp,stats,&a,&d);
  GT_CHECK_EQ_UINT(paired_matches_get_num_discordant_maps(pm),1);
  // Direct add of a discordant pair
  paired_matches_add(pm,&a,&d,pair_relation_discordant,pair_orientation_RF,pair_layout_separate,400,1.5);
  GT_CHECK_EQ_UINT(paired_matches_get_num_discordant_maps(pm),2);
  GT_CHECK_EQ_UINT(paired_matches_get_num_maps(pm),1);
  // Clearing the pairs keeps the single-end matches unless asked
  add_end_matches(pm->matches_end1,&a,1);
  paired_matches_clear(pm,false);
  GT_CHECK_EQ_UINT(paired_matches_get_num_maps(pm),0);
  GT_CHECK_EQ_UINT(paired_matches_get_num_discordant_maps(pm),0);
  GT_CHECK_EQ_UINT(matches_get_num_match_traces(pm->matches_end1),1);
  paired_matches_clear(pm,true);
  GT_CHECK_EQ_UINT(matches_get_num_match_traces(pm->matches_end1),0);
  mapper_stats_delete(stats);
  paired_matches_delete(pm);
}

/*
 * Max complete stratum of a pair, classification and search depth
 */
GT_TEST(paired_classify) {
  paired_matches_t* const pm = new_paired_matches();
  matches_t* const m1 = pm->matches_end1;
  matches_t* const m2 = pm->matches_end2;
  // MCS: the smaller of the two ends; an extended end counts its own; a
  // limited (exact-match-capped) end counts 0
  m1->max_complete_stratum = 4; m2->max_complete_stratum = 2;
  GT_CHECK_EQ_UINT(paired_matches_get_max_complete_stratum(pm),2);
  m1->matches_extended = true;
  GT_CHECK_EQ_UINT(paired_matches_get_max_complete_stratum(pm),4);
  m1->matches_extended = false;
  m2->limited_exact_matches = true;
  GT_CHECK_EQ_UINT(paired_matches_get_max_complete_stratum(pm),0);
  m2->limited_exact_matches = false;
  m2->matches_extended = true;
  GT_CHECK_EQ_UINT(paired_matches_get_max_complete_stratum(pm),2);
  m2->matches_extended = false;
  // Unmapped
  paired_matches_classify(pm);
  GT_CHECK_EQ_INT(pm->classification.matches_class,matches_class_unmapped);
  GT_CHECK(!paired_matches_classify_search_accomplished(pm));
  // Unique pair (edit 1): wdelta = mcs(2) - 1
  static paired_map_t maps[4];
  maps[0] = make_map(pair_relation_concordant,190,0.5,300,pair_layout_separate,1);
  maps[0].edit_distance = maps[0].event_distance = 1;
  paired_matches_add_map(pm,maps+0);
  paired_matches_classify(pm);
  GT_CHECK_EQ_INT(pm->classification.matches_class,matches_class_unique);
  GT_CHECK_EQ_INT(pm->classification.delta_group,-1);
  GT_CHECK_EQ_INT(pm->classification.wdelta_group,1);
  GT_CHECK(!paired_matches_classify_search_accomplished(pm));
  m1->max_complete_stratum = 5; m2->max_complete_stratum = 5;
  GT_CHECK(paired_matches_classify_search_accomplished(pm)); // wdelta 4
  // Second pair (edit 4): delta 3, wdelta min(5,4)-1
  maps[1] = make_map(pair_relation_concordant,170,0.5,300,pair_layout_separate,2);
  maps[1].edit_distance = maps[1].event_distance = 4;
  paired_matches_add_map(pm,maps+1);
  paired_matches_classify(pm);
  GT_CHECK_EQ_INT(pm->classification.matches_class,matches_class_mmap);
  GT_CHECK_EQ_INT(pm->classification.delta_group,3);
  GT_CHECK_EQ_INT(pm->classification.wdelta_group,3);
  // Subdominant end: an extension can only help when the end's mcs plus the
  // extended match's distance is below the second-best pair's distance (4)
  match_trace_t ext;
  make_trace(&ext,name_a,0,100,Forward,90,2);
  m2->max_complete_stratum = 1;
  GT_CHECK(!paired_matches_classify_subdominant_end(pm,m2,&ext)); // 4 > 1+2
  m2->max_complete_stratum = 2;
  GT_CHECK(paired_matches_classify_subdominant_end(pm,m2,&ext));  // 4 <= 2+2
  // Third pair with the best score: a tie
  maps[2] = make_map(pair_relation_concordant,190,0.5,300,pair_layout_separate,3);
  maps[2].edit_distance = maps[2].event_distance = 1;
  paired_matches_add_map(pm,maps+2);
  paired_matches_classify(pm);
  GT_CHECK_EQ_INT(pm->classification.delta_group,0);
  paired_matches_t* const empty = new_paired_matches();
  GT_CHECK(!paired_matches_classify_subdominant_end(empty,empty->matches_end2,&ext));
  paired_matches_delete(empty);
  paired_matches_delete(pm);
}

/*
 * Paired selection (-M) and template-length sampling of unique concordant pairs
 */
GT_TEST(select_pe_matches) {
  select_parameters_t s;
  select_parameters_init(&s);
  mapper_stats_t* const stats = mapper_stats_new();
  static paired_map_t maps[8];
  uint64_t i;
  // One concordant pair: its template length is sampled
  paired_matches_t* pm = new_paired_matches();
  maps[0] = make_map(pair_relation_concordant,190,0.5,321,pair_layout_separate,1);
  paired_matches_add_map(pm,maps+0);
  archive_select_pe_matches(&s,stats,pm);
  GT_CHECK_EQ_UINT(mapper_stats_template_length_get_num_samples(stats),1);
  GT_CHECK(mapper_stats_template_length_get_mean(stats)==321.0);
  paired_matches_delete(pm);
  // One discordant pair: not sampled
  pm = new_paired_matches();
  maps[0].pair_relation = pair_relation_discordant;
  paired_matches_add_map(pm,maps+0);
  archive_select_pe_matches(&s,stats,pm);
  GT_CHECK_EQ_UINT(mapper_stats_template_length_get_num_samples(stats),1);
  paired_matches_delete(pm);
  // Seven pairs, -M 3: the best three; several pairs are not sampled
  pm = new_paired_matches();
  for (i=0;i<7;++i) {
    maps[i] = make_map(pair_relation_concordant,200-(int32_t)i,0.5,300,pair_layout_separate,i);
    paired_matches_add_map(pm,maps+i);
  }
  s.max_reported_matches = 3;
  archive_select_pe_matches(&s,stats,pm);
  GT_CHECK_EQ_UINT(paired_matches_get_num_maps(pm),3);
  GT_CHECK(paired_matches_get_maps(pm)[2]==maps+2);
  GT_CHECK_EQ_UINT(mapper_stats_template_length_get_num_samples(stats),1);
  paired_matches_delete(pm);
  // Unpaired: each end selected as single-end (-M 2)
  pm = new_paired_matches();
  static match_trace_t t[6];
  for (i=0;i<6;++i) make_trace(t+i,name_a,100*i,100,Forward,100-(int32_t)i,i);
  add_end_matches(pm->matches_end1,t,4);
  add_end_matches(pm->matches_end2,t+4,2);
  s.max_reported_matches = 2;
  s.min_reported_strata_nominal = 0;
  archive_select_pe_matches(&s,stats,pm);
  GT_CHECK_EQ_UINT(matches_get_num_match_traces(pm->matches_end1),2);
  GT_CHECK_EQ_UINT(matches_get_num_match_traces(pm->matches_end2),2);
  paired_matches_delete(pm);
  mapper_stats_delete(stats);
}

/*
 * Paired-end MAPQ models on hand-built pairs
 */
typedef struct {
  archive_search_t* as1;
  archive_search_t* as2;
  paired_matches_t* pm;
  match_trace_t t[4];
} pe_case_t;
static void pe_case_init(
    pe_case_t* const c,const mapq_model_t model,const uint64_t num_maps,
    const pair_relation_t relation,const bool end_tie,const bool pair_tie,const uint64_t mcs,
    const double sub_sigma) {
  c->as1 = calloc(1,sizeof(archive_search_t));
  c->as2 = calloc(1,sizeof(archive_search_t));
  search_parameters_init(&c->as1->search_parameters);
  search_parameters_init(&c->as2->search_parameters);
  c->as1->search_parameters.mapq_model_pe = model;
  c->pm = new_paired_matches();
  make_trace(c->t+0,name_a,1000,100,Forward,98,1);
  make_trace(c->t+1,name_a,1200,100,Reverse,97,1);
  make_trace(c->t+2,name_a,5000,100,Forward,98,1);
  make_trace(c->t+3,name_a,5200,100,Reverse,pair_tie ? 97 : 80,pair_tie ? 1 : 6);
  c->t[0].mapq_score = 40; c->t[1].mapq_score = 21;
  c->pm->matches_end1->classification.delta_group = end_tie ? 0 : 3;
  c->pm->matches_end2->classification.delta_group = 3;
  c->pm->matches_end1->max_complete_stratum = mcs;
  c->pm->matches_end2->max_complete_stratum = mcs;
  c->pm->metrics.read_length = 200;
  c->pm->metrics.swg_match_score = 1;
  c->pm->metrics.proper_length = 9.0;
  uint64_t i;
  for (i=0;i<num_maps;++i) {
    paired_matches_add(c->pm,c->t+2*i,c->t+2*i+1,relation,pair_orientation_FR,pair_layout_separate,
        300,i==0 ? 0.4 : sub_sigma);
  }
  if (relation==pair_relation_discordant) paired_matches_find_discordant_pairs(c->pm,&c->as1->search_parameters);
}
static uint8_t pe_case_score(pe_case_t* const c) {
  archive_score_matches_pe(c->as1,c->as2,c->pm);
  return paired_matches_get_num_maps(c->pm) ? paired_matches_get_primary_map(c->pm)->mapq_score : 255;
}
static void pe_case_destroy(pe_case_t* const c) {
  paired_matches_delete(c->pm);
  free(c->as1);
  free(c->as2);
}
static uint8_t pe_score(
    const mapq_model_t model,const uint64_t num_maps,const pair_relation_t relation,
    const bool end_tie,const bool pair_tie,const uint64_t mcs,const double sub_sigma) {
  pe_case_t c;
  pe_case_init(&c,model,num_maps,relation,end_tie,pair_tie,mcs,sub_sigma);
  const uint8_t mapq = pe_case_score(&c);
  pe_case_destroy(&c);
  return mapq;
}
GT_TEST(score_pe_models) {
  const pair_relation_t conc = pair_relation_concordant, disc = pair_relation_discordant;
  // Discordant primary: 0 in both models
  GT_CHECK_EQ_UINT(pe_score(mapq_model_gem,1,disc,false,false,5,0.0),0);
  GT_CHECK_EQ_UINT(pe_score(mapq_model_classify,1,disc,false,false,5,0.0),0);
  // An end that ties, the pair ties: the second pair's sigma, capped (4 / 30)
  GT_CHECK_EQ_UINT(pe_score(mapq_model_gem,2,conc,true,true,5,2.7),2);
  GT_CHECK_EQ_UINT(pe_score(mapq_model_gem,2,conc,true,true,5,7.0),4);
  GT_CHECK_EQ_UINT(pe_score(mapq_model_classify,2,conc,true,true,5,2.7),12);
  GT_CHECK_EQ_UINT(pe_score(mapq_model_classify,2,conc,true,true,5,45.0),40);
  // An end ties, the pair is unique: 5 + min(wdelta,10) / 50 + min(wdelta,30); wdelta = mcs - 2
  GT_CHECK_EQ_UINT(pe_score(mapq_model_gem,1,conc,true,false,5,0.0),8);
  GT_CHECK_EQ_UINT(pe_score(mapq_model_gem,1,conc,true,false,40,0.0),15);
  GT_CHECK_EQ_UINT(pe_score(mapq_model_gem,1,conc,true,false,1,0.0),5);
  GT_CHECK_EQ_UINT(pe_score(mapq_model_classify,1,conc,true,false,5,0.0),53);
  GT_CHECK_EQ_UINT(pe_score(mapq_model_classify,1,conc,true,false,60,0.0),80);
  // Ends do not tie, the pair ties: 0 / 2
  GT_CHECK_EQ_UINT(pe_score(mapq_model_gem,2,conc,false,true,5,0.0),0);
  GT_CHECK_EQ_UINT(pe_score(mapq_model_classify,2,conc,false,true,5,0.0),2);
  // Unique pair well separated (wdelta 3): 60 / 161
  GT_CHECK_EQ_UINT(pe_score(mapq_model_gem,1,conc,false,false,5,0.0),60);
  GT_CHECK_EQ_UINT(pe_score(mapq_model_classify,1,conc,false,false,5,0.0),161);
  // Unique pair, wdelta <= 1: the mean of the ends' MAPQ (40, 21) / 100 + mean
  GT_CHECK_EQ_UINT(pe_score(mapq_model_gem,1,conc,false,false,3,0.0),30);
  GT_CHECK_EQ_UINT(pe_score(mapq_model_classify,1,conc,false,false,3,0.0),130);
  // Second pair clearly worse (delta 5 > 1, wdelta min(mcs,7)-2): 60
  GT_CHECK_EQ_UINT(pe_score(mapq_model_gem,2,conc,false,false,9,0.0),60);
  // Subdominant pairs get 0 (gem)
  {
    pe_case_t c;
    pe_case_init(&c,mapq_model_gem,2,conc,false,false,9,0.0);
    paired_matches_get_maps(c.pm)[1]->mapq_score = 9;
    pe_case_score(&c);
    GT_CHECK_EQ_UINT(paired_matches_get_maps(c.pm)[1]->mapq_score,0);
    pe_case_destroy(&c);
  }
  // A local end: 1 (classify)
  {
    pe_case_t c;
    pe_case_init(&c,mapq_model_classify,1,conc,false,false,5,0.0);
    c.t[1].type = match_type_local;
    GT_CHECK_EQ_UINT(pe_case_score(&c),1);
    pe_case_destroy(&c);
  }
  // No model, no alignment model, no pairs: nothing changes
  {
    pe_case_t c;
    pe_case_init(&c,mapq_model_none,1,conc,false,false,5,0.0);
    paired_matches_get_primary_map(c.pm)->mapq_score = 77;
    GT_CHECK_EQ_UINT(pe_case_score(&c),77);
    pe_case_destroy(&c);
    pe_case_init(&c,mapq_model_gem,1,conc,false,false,5,0.0);
    c.as1->search_parameters.match_alignment_model = match_alignment_model_none;
    paired_matches_get_primary_map(c.pm)->mapq_score = 77;
    GT_CHECK_EQ_UINT(pe_case_score(&c),77);
    pe_case_destroy(&c);
    pe_case_init(&c,mapq_model_gem,0,conc,false,false,5,0.0);
    GT_CHECK_EQ_UINT(pe_case_score(&c),255);
    archive_score_matches_pe_default(c.as1,c.as2,c.pm);
    archive_score_matches_pe_classify(c.as1,c.as2,c.pm);
    GT_CHECK_EQ_UINT(paired_matches_get_num_maps(c.pm),0);
    pe_case_destroy(&c);
  }
}

GT_TEST(predictors_pe) {
  pe_case_t c;
  pe_case_init(&c,mapq_model_gem,2,pair_relation_concordant,false,false,9,3.0);
  c.pm->metrics.avg_region_length = 18;
  c.pm->metrics.max_region_length = 27;
  c.pm->metrics.kmer_frequency = 0.25;
  matches_predictors_t p;
  matches_predictors_compute_pe(&p,c.pm);
  GT_CHECK_EQ_INT(c.pm->classification.matches_class,matches_class_mmap);
  GT_CHECK(fabs(p.primary_edit_distance_norm-198.0/200.0) < 1e-12);
  GT_CHECK(fabs(p.primary_swg_score_norm-195.0/200.0) < 1e-12);
  GT_CHECK(fabs(p.primary_template_size_sigma_norm-0.04) < 1e-12);
  GT_CHECK(fabs(p.subdominant_edit_distance_norm-193.0/200.0) < 1e-12);
  GT_CHECK(fabs(p.subdominant_swg_score_norm-178.0/200.0) < 1e-12);
  GT_CHECK_EQ_UINT(p.mcs_end1,9);
  GT_CHECK_EQ_UINT(p.mcs_end2,0);
  GT_CHECK_EQ_UINT(p.mapq_end1,40);
  GT_CHECK_EQ_UINT(p.mapq_end2,21);
  GT_CHECK(fabs(p.avg_region_length_norm-2.0) < 1e-12);
  GT_CHECK(fabs(p.max_region_length_norm-3.0) < 1e-12);
  pe_case_destroy(&c);
  // One pair: the unique form
  pe_case_init(&c,mapq_model_gem,1,pair_relation_concordant,false,false,9,0.0);
  matches_predictors_compute_pe(&p,c.pm);
  GT_CHECK_EQ_INT(c.pm->classification.matches_class,matches_class_unique);
  GT_CHECK(p.subdominant_template_size_sigma_norm==1.0);
  GT_CHECK(p.subdominant_swg_score_norm==0.0);
  pe_case_destroy(&c);
}

/*
 * Extension (mate rescue) switches
 */
GT_TEST(pe_extension_switches) {
  archive_search_t* const as1 = calloc(1,sizeof(archive_search_t));
  archive_search_t* const as2 = calloc(1,sizeof(archive_search_t));
  search_parameters_init(&as1->search_parameters);
  search_paired_parameters_init(&as1->search_parameters.search_paired_parameters);
  as1->approximate_search.pattern.key_length = 100;
  as2->approximate_search.pattern.key_length = 100;
  GT_CHECK(archive_search_pe_use_recovery_extension(as1,as2,NULL));
  as2->approximate_search.pattern.key_length = 300;
  GT_CHECK(archive_search_pe_use_recovery_extension(as1,as2,NULL));
  as2->approximate_search.pattern.key_length = 301;   // longer reads are not extended
  GT_CHECK(!archive_search_pe_use_recovery_extension(as1,as2,NULL));
  as2->approximate_search.pattern.key_length = 100;
  as1->approximate_search.pattern.key_length = 0;     // empty end
  GT_CHECK(!archive_search_pe_use_recovery_extension(as1,as2,NULL));
  as1->approximate_search.pattern.key_length = 100;
  as1->search_parameters.search_paired_parameters.paired_end_extension_recovery = false;
  GT_CHECK(!archive_search_pe_use_recovery_extension(as1,as2,NULL));
  // The shortcut extension is disabled in this version
  as1->search_parameters.search_paired_parameters.paired_end_extension_shortcut = true;
  GT_CHECK(!archive_search_pe_use_shortcut_extension(as1,as2,NULL));
  free(as1);
  free(as2);
}

GT_CASES(
  GT_CASE(paired_map_scores),
  GT_CASE(orientation_layout_template_length),
  GT_CASE(search_paired_parameters_defaults),
  GT_CASE(relation_tables_and_limits),
  GT_CASE(find_pairs_all_combinations),
  GT_CASE(paired_map_order),
  GT_CASE(add_and_cross_pair),
  GT_CASE(paired_classify),
  GT_CASE(select_pe_matches),
  GT_CASE(score_pe_models),
  GT_CASE(predictors_pe),
  GT_CASE(pe_extension_switches)
)
