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
 *   The loaded GEM index (archive, locator, archive text, FM-index: BWT,
 *   rank table, sampled SA, bidirectional queries; run-length text) against
 *   naive structures computed in this driver from the multi-FASTA it indexes.
 *
 *   Each case writes a small multi-FASTA to a temporary directory, builds an
 *   index from it in-process with the archive builder (the steps gem-indexer
 *   runs), loads it with archive_read() as gem-mapper does, and compares:
 *   the indexed text (forward contigs with N runs of 50 or more removed,
 *   one separator after each interval, then the reverse complement and a
 *   final separator) with one rebuilt here from the FASTA; the suffix array
 *   (cyclic rotations sorted by encoded character, A<C<G<T<N<separator) with
 *   fm_index_decode(); every BWT character and exclusive rank; backward,
 *   rank-table and bidirectional search intervals with naive counts; and
 *   locator positions with the FASTA coordinates.
 *
 *   Internal signatures called (a speed iteration may change them; update
 *   this driver with them):
 *     archive_builder_new(char* fasta,fm_t* out,char* prefix,archive_type,bool gpu_index,
 *         uint64_t ns_threshold,sampling_rate_t sa_rate,sampling_rate_t text_rate,
 *         uint64_t num_threads,uint64_t max_memory,FILE* info_file)
 *     archive_builder_text_process(archive_builder_t*,bool verbose)
 *     archive_builder_text_generate_run_length(archive_builder_t*,bool verbose)
 *     archive_builder_write_header/_write_locator(archive_builder_t*), locator_builder_delete(locator_builder_t*)
 *     archive_builder_index_build_bwt(archive_builder_t*,bool gpu,bool dump_bwt,bool dump_sa,bool verbose)
 *     archive_builder_write_index(archive_builder_t*,bool gpu,bool check,bool verbose), archive_builder_delete()
 *     archive_read(char* file_name,bool read_text_only), archive_delete(archive_t*)
 *     fm_index_get_length(const fm_index_t*), fm_index_decode(const fm_index_t*,uint64_t)
 *     fm_index_retrieve_bwt_sampled(const fm_index_t*,uint64_t,uint64_t* sampled,uint64_t* lf_dist)
 *     fm_index_retrieve_sa_sample(const fm_index_t*,uint64_t sampled,uint64_t lf_dist,uint64_t* text_position)
 *     fm_index_bsearch / fm_index_bsearch_pure(const fm_index_t*,const uint8_t* key,uint64_t length,uint64_t* hi,uint64_t* lo)
 *     bwt_sbasic_char(bwt,pos), bwt_sbasic_erank(bwt,enc,pos), bwt_sbasic_erank_interval(bwt,enc,lo,hi,&lo,&hi),
 *     bwt_sbasic_prefetch(bwt,pos,bwt_block_locator_t*), bwt_sbasic_prefetched_erank(bwt,enc,pos,locator),
 *     bwt_sbasic_precompute(bwt,pos,locator,bwt_block_elms_t*), bwt_sbasic_precomputed_erank(bwt,enc,locator,elms),
 *     bwt_sbasic_LF(bwt,pos,bool* sampled), bwt_sbasic_prefetched_LF(bwt,pos,bool*,locator),
 *     bwt_sbasic_LF__enc(bwt,pos,uint8_t* enc,bool* sampled), bwt_sbasic_prefetched_LF__enc(...)
 *     rank_mquery_new/rank_mquery_add_char/rank_mquery_is_exhausted, rank_mtable_prefetch, rank_mtable_fetch(mtable,query,&lo,&hi)
 *     fm_index_2query_init, fm_index_2query_forward_query / _backward_query(fm_index,in,out,enc),
 *     fm_index_2query_forward_precompute / _backward_precompute(fm_index,interval,lo_elms,hi_elms),
 *     fm_index_2query_precomputed_forward_query / _backward_query(lo_elms,hi_elms,in,out,enc)
 *     locator_map(const locator_t*,uint64_t index_position,location_t*): for a reverse-strand
 *         index position, location.position is the forward coordinate of that base plus one
 *         (gem-mapper subtracts the match length from it, src/matches/matches.c)
 *     locator_lookup_interval(const locator_t*,uint64_t), locator_interval_get_tag(locator,interval)
 *     locator_inverse_map_position(locator,tag,strand,bs_strand,text_position)
 *     locator_lookup_rl_interval(const locator_t*,uint64_t rl_position)
 *     archive_text_retrieve(archive_text_t*,uint64_t pos,uint64_t len,bool rc,bool rl,text_trace_t*,mm_allocator_t*)
 *     archive_text_get_position_strand(archive_text_t*,uint64_t)
 *     archive_text_rl_position_translate(archive_text_t*,uint64_t rl_position,mm_allocator_t*)
 *     dna_text_get_length(const dna_text_t*), dna_text_retrieve_sequence(const dna_text_t*,uint64_t)
 *     sampled_sa_get_sa_sampling_rate / sampled_sa_get_text_sampling_rate(const sampled_sa_t*)
 *     sampled_rl_get_sample(sampled_rl_t*,uint64_t)
 */

#include "gem3_test.h"
#include <unistd.h>
#include "archive/archive.h"
#include "archive/archive_text.h"
#include "archive/archive_text_rl.h"
#include "archive/sampled_rl.h"
#include "archive/locator.h"
#include "archive/builder/archive_builder.h"
#include "archive/builder/archive_builder_text.h"
#include "archive/builder/archive_builder_index.h"
#include "fm_index/fm_index.h"
#include "fm_index/fm_index_search.h"
#include "fm_index/fm_index_query.h"
#include "fm_index/rank_mtable.h"
#include "fm_index/sampled_sa.h"
#include "fm_index/bwt/bwt.h"
#include "text/dna_text.h"
#include "text/text_trace.h"
#include "filtering/region_profile/region_profile.h"
#include "filtering/region_profile/region_profile_fixed.h"
#include "filtering/region_profile/region_profile_adaptive.h"
#include "filtering/region_profile/region_profile_mem.h"
#include "filtering/region_profile/region_profile_optimum.h"
#include "filtering/region_profile/region_profile_schedule.h"

#define NS_THRESHOLD 50
#define MAX_CONTIGS 8
#define MAX_TEXT 16384

/*
 * Test genome: contigs given as segments so the expected indexed text is
 * known without re-implementing the indexer
 */
typedef struct {
  char name[32];
  char* seq;              // FASTA sequence (case and IUPAC codes as written)
  uint64_t length;
} contig_t;
typedef struct {
  contig_t contigs[MAX_CONTIGS];
  uint64_t num_contigs;
  // Expected indexed text (full: forward + RC + final separator)
  uint8_t text[MAX_TEXT];
  uint64_t text_length;
  uint64_t forward_length;
  // Forward index position -> contig and offset (-1 on separators)
  int64_t pos_contig[MAX_TEXT];
  int64_t pos_offset[MAX_TEXT];
} genome_t;

static char random_base(gt_rng_t* const rng) {
  return "ACGT"[gt_rng_below(rng,4)];
}
static uint8_t enc_of(const char c) {
  switch (c) {
    case 'A': case 'a': return ENC_DNA_CHAR_A;
    case 'C': case 'c': return ENC_DNA_CHAR_C;
    case 'G': case 'g': return ENC_DNA_CHAR_G;
    case 'T': case 't': return ENC_DNA_CHAR_T;
    default: return ENC_DNA_CHAR_N; // N and IUPAC codes
  }
}
static uint8_t enc_complement(const uint8_t e) {
  return (e < 4) ? (uint8_t)(3-e) : e;
}
static void contig_append(contig_t* const contig,const char* const s,const uint64_t n) {
  memcpy(contig->seq+contig->length,s,n);
  contig->length += n;
  contig->seq[contig->length] = '\0';
}
static void contig_append_random(contig_t* const contig,gt_rng_t* const rng,const uint64_t n) {
  uint64_t i;
  for (i=0;i<n;++i) {
    char c = random_base(rng);
    contig_append(contig,&c,1);
  }
}
static void contig_append_char(contig_t* const contig,const char c,const uint64_t n) {
  uint64_t i;
  for (i=0;i<n;++i) contig_append(contig,&c,1);
}
/*
 * Expected forward text from the FASTA (the indexer's documented rules):
 * runs of N/IUPAC shorter than NS_THRESHOLD stay as N; longer runs are
 * dropped and end the current interval; every non-empty interval is
 * followed by one separator. Then the reverse complement of everything but
 * the last separator, a separator, and a final separator.
 */
static void genome_compute_text(genome_t* const g) {
  uint64_t p = 0, c;
  for (c=0;c<g->num_contigs;++c) {
    const contig_t* const contig = g->contigs + c;
    uint64_t i = 0, interval_length = 0;
    while (i < contig->length) {
      if (enc_of(contig->seq[i]) == ENC_DNA_CHAR_N) {
        uint64_t j = i;
        while (j < contig->length && enc_of(contig->seq[j]) == ENC_DNA_CHAR_N) ++j;
        if (j-i < NS_THRESHOLD) {
          for (;i<j;++i) {
            g->text[p] = ENC_DNA_CHAR_N; g->pos_contig[p] = c; g->pos_offset[p] = i; ++p; ++interval_length;
          }
        } else {
          if (interval_length > 0) {
            g->text[p] = ENC_DNA_CHAR_SEP; g->pos_contig[p] = -1; g->pos_offset[p] = -1; ++p;
          }
          interval_length = 0;
          i = j;
        }
      } else {
        g->text[p] = enc_of(contig->seq[i]); g->pos_contig[p] = c; g->pos_offset[p] = i; ++p; ++interval_length;
        ++i;
      }
    }
    if (interval_length > 0) {
      g->text[p] = ENC_DNA_CHAR_SEP; g->pos_contig[p] = -1; g->pos_offset[p] = -1; ++p;
    }
  }
  g->forward_length = p;
  uint64_t k;
  for (k=0;k+1<g->forward_length;++k) {
    g->text[p++] = enc_complement(g->text[g->forward_length-2-k]);
  }
  g->text[p++] = ENC_DNA_CHAR_SEP;
  g->text[p++] = ENC_DNA_CHAR_SEP;
  g->text_length = p;
}
static void genome_free(genome_t* const g) {
  uint64_t c;
  for (c=0;c<g->num_contigs;++c) free(g->contigs[c].seq);
  free(g);
}
static contig_t* genome_add_contig(genome_t* const g,const char* const name) {
  contig_t* const contig = g->contigs + g->num_contigs++;
  snprintf(contig->name,sizeof(contig->name),"%s",name);
  contig->seq = calloc(MAX_TEXT,1);
  contig->length = 0;
  return contig;
}
/*
 * Genome with repeats, short and long N runs (inside contigs, at contig
 * starts and ends), IUPAC codes, soft-masked bases and a very short contig
 */
static genome_t* genome_new(const uint64_t seed) {
  genome_t* const g = calloc(1,sizeof(genome_t));
  gt_rng_t rng = gt_rng_new(seed);
  char repeat[160];
  uint64_t i;
  for (i=0;i<sizeof(repeat);++i) repeat[i] = random_base(&rng);
  // chrA: random, an exact repeat copy, a short N run, IUPAC codes, lowercase
  contig_t* c = genome_add_contig(g,"chrA");
  contig_append_random(c,&rng,300);
  contig_append(c,repeat,sizeof(repeat));
  contig_append_random(c,&rng,200);
  contig_append_char(c,'N',7);            // kept (shorter than the threshold)
  contig_append_random(c,&rng,150);
  contig_append(c,"RYKM",4);              // IUPAC codes: indexed as N
  contig_append_random(c,&rng,120);
  for (i=0;i<90;++i) { char b = (char)tolower(random_base(&rng)); contig_append(c,&b,1); }
  contig_append_random(c,&rng,200);
  contig_append(c,repeat,sizeof(repeat));
  contig_append_random(c,&rng,100);
  // chrB: starts with a long N run, a long N run inside, the repeat reverse-complemented
  c = genome_add_contig(g,"chrB");
  contig_append_char(c,'N',70);
  contig_append_random(c,&rng,250);
  contig_append_char(c,'N',NS_THRESHOLD);  // removed (equal to the threshold)
  contig_append_random(c,&rng,180);
  for (i=0;i<sizeof(repeat);++i) {
    const char b = repeat[sizeof(repeat)-1-i];
    const char rc = (b=='A') ? 'T' : (b=='C') ? 'G' : (b=='G') ? 'C' : 'A';
    contig_append(c,&rc,1);
  }
  contig_append_random(c,&rng,140);
  contig_append_char(c,'N',NS_THRESHOLD-1); // kept
  contig_append_random(c,&rng,60);
  contig_append_char(c,'A',40);            // homopolymer
  contig_append_random(c,&rng,90);
  // chrC: ends with a long N run
  c = genome_add_contig(g,"chrC");
  contig_append_random(c,&rng,330);
  contig_append_char(c,'N',61);
  // chrD: very short
  c = genome_add_contig(g,"chrD");
  contig_append_random(c,&rng,37);
  // chrE: tandem repeat
  c = genome_add_contig(g,"chrE");
  contig_append_random(c,&rng,80);
  for (i=0;i<30;++i) contig_append(c,"CAG",3);
  contig_append_random(c,&rng,80);
  genome_compute_text(g);
  return g;
}

/*
 * Index fixture: FASTA -> in-process index -> archive_read()
 */
typedef struct {
  char dir[512];
  char fasta[600];
  char prefix[600];
  char gem[600];
  archive_t* archive;
} fixture_t;

static void write_fasta(const genome_t* const g,const char* const path) {
  FILE* const f = fopen(path,"w");
  GT_REQUIRE(f != NULL);
  uint64_t c, i;
  for (c=0;c<g->num_contigs;++c) {
    fprintf(f,">%s some description\n",g->contigs[c].name);
    for (i=0;i<g->contigs[c].length;i+=60) {
      const uint64_t n = MIN(60,g->contigs[c].length-i);
      fwrite(g->contigs[c].seq+i,1,n,f);
      fputc('\n',f);
    }
  }
  fclose(f);
}
static void fixture_build(
    fixture_t* const fx,const genome_t* const g,const archive_type type,
    const sampling_rate_t sa_rate,const sampling_rate_t text_rate,const bool run_length) {
  const char* tmp = getenv("TMPDIR");
  snprintf(fx->dir,sizeof(fx->dir),"%s/gem3-test-a-XXXXXX",(tmp && *tmp) ? tmp : "/tmp");
  GT_REQUIRE(mkdtemp(fx->dir) != NULL);
  snprintf(fx->fasta,sizeof(fx->fasta),"%s/genome.fa",fx->dir);
  snprintf(fx->prefix,sizeof(fx->prefix),"%s/index",fx->dir);
  snprintf(fx->gem,sizeof(fx->gem),"%s/index.gem",fx->dir);
  write_fasta(g,fx->fasta);
  // The builder logs a line even when not verbose; keep TAP output clean
  FILE* const log_stream = gem_log_get_stream();
  FILE* const devnull = fopen("/dev/null","w");
  if (devnull) gem_log_set_stream(devnull);
  fm_t* const out = fm_open_file(fx->gem,FM_WRITE);
  archive_builder_t* const builder = archive_builder_new(
      fx->fasta,out,fx->prefix,type,false,NS_THRESHOLD,sa_rate,text_rate,1,BUFFER_SIZE_1G,NULL);
  archive_builder_text_process(builder,false);
  if (run_length) archive_builder_text_generate_run_length(builder,false);
  archive_builder_write_header(builder);
  archive_builder_write_locator(builder);
  locator_builder_delete(builder->locator);
  archive_builder_index_build_bwt(builder,false,false,false,false);
  archive_builder_write_index(builder,false,false,false);
  archive_builder_delete(builder);
  gem_log_set_stream(log_stream);
  if (devnull) fclose(devnull);
  fx->archive = archive_read(fx->gem,false);
}
static void fixture_destroy(fixture_t* const fx) {
  archive_delete(fx->archive);
  unlink(fx->gem);
  unlink(fx->fasta);
  char sa_tmp[700];
  snprintf(sa_tmp,sizeof(sa_tmp),"%s.sa.tmp",fx->prefix);
  unlink(sa_tmp);
  rmdir(fx->dir);
}

/*
 * Naive suffix array of the cyclic text
 */
static const uint8_t* naive_text;
static uint64_t naive_length;
static int rotation_cmp(const void* const a,const void* const b) {
  const uint64_t i = *(const uint64_t*)a, j = *(const uint64_t*)b;
  uint64_t k;
  for (k=0;k<naive_length;++k) {
    const uint8_t ci = naive_text[(i+k)%naive_length], cj = naive_text[(j+k)%naive_length];
    if (ci != cj) return (ci < cj) ? -1 : 1;
  }
  return 0;
}
static uint64_t* naive_sa(const uint8_t* const text,const uint64_t n) {
  uint64_t* const sa = malloc(n*sizeof(uint64_t));
  uint64_t i;
  for (i=0;i<n;++i) sa[i] = i;
  naive_text = text;
  naive_length = n;
  qsort(sa,n,sizeof(uint64_t),rotation_cmp);
  return sa;
}
// Compare the first m characters of rotation i with key (m = key length)
static int prefix_cmp(const uint8_t* const text,const uint64_t n,const uint64_t i,const uint8_t* const key,const uint64_t m) {
  uint64_t k;
  for (k=0;k<m;++k) {
    const uint8_t c = text[(i+k)%n];
    if (c != key[k]) return (c < key[k]) ? -1 : 1;
  }
  return 0;
}
// [lo,hi): rotations that start with key
static void naive_range(
    const uint8_t* const text,const uint64_t n,const uint64_t* const sa,
    const uint8_t* const key,const uint64_t m,uint64_t* const lo,uint64_t* const hi) {
  uint64_t l = 0, h = n;
  while (l < h) { const uint64_t mid = (l+h)/2; if (prefix_cmp(text,n,sa[mid],key,m) < 0) l = mid+1; else h = mid; }
  *lo = l;
  h = n;
  while (l < h) { const uint64_t mid = (l+h)/2; if (prefix_cmp(text,n,sa[mid],key,m) <= 0) l = mid+1; else h = mid; }
  *hi = l;
}

/*
 * Shared state for one built index
 */
typedef struct {
  genome_t* g;
  fixture_t fx;
  uint64_t* sa;
  uint8_t* bwt;      // naive BWT
} indexed_t;
static void indexed_new(indexed_t* const ix,const uint64_t seed,
    const sampling_rate_t sa_rate,const sampling_rate_t text_rate) {
  ix->g = genome_new(seed);
  fixture_build(&ix->fx,ix->g,archive_dna_full,sa_rate,text_rate,false);
  ix->sa = naive_sa(ix->g->text,ix->g->text_length);
  ix->bwt = malloc(ix->g->text_length);
  const uint64_t n = ix->g->text_length;
  uint64_t i;
  for (i=0;i<n;++i) ix->bwt[i] = ix->g->text[(ix->sa[i]+n-1)%n];
}
static void indexed_delete(indexed_t* const ix) {
  fixture_destroy(&ix->fx);
  free(ix->sa);
  free(ix->bwt);
  genome_free(ix->g);
}

/*
 * The indexed text equals the FASTA's, forward and reverse complement
 */
static void check_text(indexed_t* const ix) {
  const genome_t* const g = ix->g;
  archive_t* const archive = ix->fx.archive;
  archive_text_t* const text = archive->text;
  GT_CHECK_EQ_UINT(text->forward_text_length,g->forward_length);
  GT_CHECK_EQ_UINT(fm_index_get_length(archive->fm_index),g->text_length);
  GT_CHECK_EQ_UINT(archive_get_index_length(archive),g->text_length);
  GT_REQUIRE(dna_text_get_length(text->enc_text) >= g->forward_length);
  const uint8_t* const stored = dna_text_retrieve_sequence(text->enc_text,0);
  uint64_t i, mismatches = 0;
  for (i=0;i<g->forward_length;++i) mismatches += (stored[i] != g->text[i]);
  GT_CHECK_EQ_UINT(mismatches,0);
  // archive_text_retrieve: any position and length, either strand
  mm_slab_t* const slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM);
  mm_allocator_t* const mm_allocator = mm_allocator_new(slab);
  gt_rng_t rng = gt_rng_new(5);
  uint64_t trial;
  for (trial=0;trial<4000;++trial) {
    // Ranges stay within one strand's half, as candidate texts do (they lie
    // inside one locator interval)
    const uint64_t length = gt_rng_range(&rng,1,150);
    const bool reverse_half = gt_rng_chance(&rng,1,2);
    const uint64_t position = reverse_half ?
        g->forward_length + gt_rng_below(&rng,g->forward_length-length) :
        gt_rng_below(&rng,g->forward_length-length+1);
    const bool rc = gt_rng_chance(&rng,1,2);
    text_trace_t trace;
    mm_allocator_push_state(mm_allocator);
    archive_text_retrieve(text,position,length,rc,false,&trace,mm_allocator);
    GT_CHECK_EQ_UINT(trace.text_length,length);
    GT_CHECK_EQ_INT(trace.strand,(position < g->forward_length) ? Forward : Reverse);
    GT_CHECK_EQ_INT(archive_text_get_position_strand(text,position),(position < g->forward_length) ? Forward : Reverse);
    uint64_t k, bad = 0;
    for (k=0;k<length;++k) {
      const uint8_t expected = rc ? enc_complement(g->text[position+length-1-k]) : g->text[position+k];
      bad += (trace.text[k] != expected);
    }
    GT_CHECK_MSG(bad==0,"archive_text_retrieve(pos=%" PRIu64 ",len=%" PRIu64 ",rc=%d): %" PRIu64 " wrong bases",
        position,length,(int)rc,bad);
    mm_allocator_pop_state(mm_allocator);
  }
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * BWT characters, exclusive ranks and LF against the naive BWT
 */
static void check_bwt(indexed_t* const ix) {
  const uint64_t n = ix->g->text_length;
  const fm_index_t* const fm_index = ix->fx.archive->fm_index;
  const bwt_t* const bwt = fm_index->bwt;
  uint64_t C[DNA_EXT_RANGE+1] = {0}, occ[DNA_EXT_RANGE] = {0};
  uint64_t i, c;
  for (i=0;i<n;++i) ++C[ix->g->text[i]+1];
  for (c=1;c<=DNA_EXT_RANGE;++c) C[c] += C[c-1];
  uint64_t char_errors = 0, rank_errors = 0, lf_errors = 0, pre_errors = 0;
  for (i=0;i<=n;++i) {
    // Exclusive rank of every character that appears in reads or text
    for (c=0;c<=ENC_DNA_CHAR_N;++c) {
      const uint64_t expected = C[c] + occ[c];
      if (bwt_sbasic_erank(bwt,(uint8_t)c,i) != expected) {
        if (rank_errors++ < 5) gt_note("erank(%" PRIu64 ",%" PRIu64 ")=%" PRIu64 " expected %" PRIu64,
            c,i,bwt_sbasic_erank(bwt,(uint8_t)c,i),expected);
      }
      bwt_block_locator_t loc;
      bwt_block_elms_t elms;
      bwt_sbasic_precompute(bwt,i,&loc,&elms);
      if (bwt_sbasic_precomputed_erank(bwt,(uint8_t)c,&loc,&elms) != expected) ++pre_errors;
      bwt_sbasic_prefetch(bwt,i,&loc);
      if (bwt_sbasic_prefetched_erank(bwt,(uint8_t)c,i,&loc) != expected) ++pre_errors;
    }
    if (i == n) break;
    if (bwt_sbasic_char(bwt,i) != ix->bwt[i]) ++char_errors;
    // LF of an unsampled position: the BWT position of the preceding suffix
    bool sampled;
    uint8_t enc;
    const uint64_t lf = bwt_sbasic_LF__enc(bwt,i,&enc,&sampled);
    bwt_block_locator_t loc;
    bwt_sbasic_prefetch(bwt,i,&loc);
    bool sampled_p;
    const uint64_t lf_p = bwt_sbasic_prefetched_LF(bwt,i,&sampled_p,&loc);
    bool sampled_q;
    const uint64_t lf_q = bwt_sbasic_LF(bwt,i,&sampled_q);
    if (sampled != sampled_p || sampled != sampled_q || lf != lf_p || lf != lf_q) ++lf_errors;
    if (!sampled) {
      if (enc != ix->bwt[i]) ++lf_errors;
      const uint64_t expected = C[ix->bwt[i]] + occ[ix->bwt[i]];
      if (lf != expected || ix->sa[lf] != (ix->sa[i]+n-1)%n) ++lf_errors;
    }
    ++occ[ix->bwt[i]];
  }
  GT_CHECK_EQ_UINT(char_errors,0);
  GT_CHECK_EQ_UINT(rank_errors,0);
  GT_CHECK_EQ_UINT(pre_errors,0);
  GT_CHECK_EQ_UINT(lf_errors,0);
  // Interval ranks
  gt_rng_t rng = gt_rng_new(17);
  uint64_t trial, interval_errors = 0;
  for (trial=0;trial<3000;++trial) {
    const uint64_t hi = gt_rng_range(&rng,1,n);
    const uint64_t lo = hi - gt_rng_range(&rng,0,MIN(hi,40));
    if (!bwt_sbasic_is_same_bucket(lo,hi)) continue;
    const uint8_t enc = (uint8_t)gt_rng_below(&rng,4);
    uint64_t lo_out, hi_out;
    bwt_sbasic_erank_interval(bwt,enc,lo,hi,&lo_out,&hi_out);
    if (lo_out != bwt_sbasic_erank(bwt,enc,lo) || hi_out != bwt_sbasic_erank(bwt,enc,hi)) ++interval_errors;
  }
  GT_CHECK_EQ_UINT(interval_errors,0);
}

/*
 * Decoding every BWT position gives the naive suffix array
 */
static void check_decode(indexed_t* const ix) {
  const uint64_t n = ix->g->text_length;
  const fm_index_t* const fm_index = ix->fx.archive->fm_index;
  uint64_t i, errors = 0, errors2 = 0;
  for (i=0;i<n;++i) {
    const uint64_t pos = fm_index_decode(fm_index,i);
    if (pos != ix->sa[i]) {
      if (errors++ < 5) gt_note("decode(%" PRIu64 ")=%" PRIu64 " expected %" PRIu64,i,pos,ix->sa[i]);
    }
    uint64_t sampled, dist, text_position;
    fm_index_retrieve_bwt_sampled(fm_index,i,&sampled,&dist);
    fm_index_retrieve_sa_sample(fm_index,sampled,dist,&text_position);
    if (text_position != ix->sa[i]) ++errors2;
  }
  GT_CHECK_EQ_UINT(errors,0);
  GT_CHECK_EQ_UINT(errors2,0);
}

/*
 * Search keys: substrings of the text (both strands, some over the short N
 * runs), mutated substrings and random strings
 */
static uint64_t make_key(indexed_t* const ix,gt_rng_t* const rng,uint8_t* const key,const bool acgt_only) {
  const uint64_t n = ix->g->text_length;
  const uint64_t m = gt_rng_range(rng,1,gt_rng_chance(rng,1,4) ? 12 : 64);
  uint64_t i;
  switch (gt_rng_below(rng,4)) {
    case 0: // random
      for (i=0;i<m;++i) key[i] = (uint8_t)gt_rng_below(rng,4);
      return m;
    default: { // substring (maybe mutated), without separators
      uint64_t tries;
      for (tries=0;tries<1000;++tries) {
        const uint64_t p = gt_rng_below(rng,n-m);
        bool ok = true;
        for (i=0;i<m && ok;++i) {
          const uint8_t c = ix->g->text[p+i];
          ok = (c < ENC_DNA_CHAR_SEP) && (!acgt_only || c < ENC_DNA_CHAR_N);
          key[i] = c;
        }
        if (!ok) continue;
        if (gt_rng_chance(rng,1,3)) key[gt_rng_below(rng,m)] = (uint8_t)gt_rng_below(rng,4);
        return m;
      }
      for (i=0;i<m;++i) key[i] = (uint8_t)gt_rng_below(rng,4);
      return m;
    }
  }
}
// Keys that run into the kept N runs (e.g. a read over the 7 N of chrA)
static uint64_t make_key_with_n(indexed_t* const ix,gt_rng_t* const rng,uint8_t* const key) {
  const uint64_t n = ix->g->text_length;
  uint64_t tries;
  for (tries=0;tries<100000;++tries) {
    const uint64_t m = gt_rng_range(rng,2,40);
    const uint64_t p = gt_rng_below(rng,n-m);
    bool has_n = false, ok = true;
    uint64_t i;
    for (i=0;i<m && ok;++i) {
      key[i] = ix->g->text[p+i];
      ok = key[i] < ENC_DNA_CHAR_SEP;
      has_n = has_n || key[i]==ENC_DNA_CHAR_N;
    }
    if (ok && has_n) return m;
  }
  return 0;
}

static void check_bsearch(indexed_t* const ix) {
  const uint64_t n = ix->g->text_length;
  const fm_index_t* const fm_index = ix->fx.archive->fm_index;
  gt_rng_t rng = gt_rng_new(29);
  uint8_t key[128];
  uint64_t trial, errors = 0, nonzero = 0;
  for (trial=0;trial<6000;++trial) {
    const uint64_t m = (trial%10==9) ? make_key_with_n(ix,&rng,key) : make_key(ix,&rng,key,false);
    if (m == 0) continue;
    uint64_t lo, hi, nlo, nhi, plo, phi;
    naive_range(ix->g->text,n,ix->sa,key,m,&nlo,&nhi);
    fm_index_bsearch(fm_index,key,m,&hi,&lo);
    fm_index_bsearch_pure(fm_index,key,m,&phi,&plo);
    // An empty interval may be reported at any position
    const bool ok = (nhi > nlo) ? (lo==nlo && hi==nhi && plo==nlo && phi==nhi) : (hi==lo && phi==plo);
    if (!ok && errors++ < 5) {
      gt_note("bsearch m=%" PRIu64 ": [%" PRIu64 ",%" PRIu64 ") pure [%" PRIu64 ",%" PRIu64 ") naive [%" PRIu64 ",%" PRIu64 ")",
          m,lo,hi,plo,phi,nlo,nhi);
    }
    nonzero += (nhi > nlo);
  }
  GT_CHECK_EQ_UINT(errors,0);
  GT_CHECK(nonzero > 3000); // most keys occur
}

/*
 * Rank table: every k-mer up to 5 and random ones up to the table depth
 */
static void check_rank_mtable(indexed_t* const ix) {
  const uint64_t n = ix->g->text_length;
  const rank_mtable_t* const mtable = ix->fx.archive->fm_index->rank_table;
  uint8_t key[128];
  uint64_t errors = 0, k, code;
  for (k=1;k<=5;++k) {
    for (code=0;code<(1ull<<(2*k));++code) {
      uint64_t i;
      for (i=0;i<k;++i) key[i] = (uint8_t)((code >> (2*i)) & 3);
      rank_mquery_t query;
      rank_mquery_new(&query);
      for (i=k;i-->0;) rank_mquery_add_char(mtable,&query,key[i]);
      GT_CHECK_EQ_UINT(rank_mquery_get_level(&query),k);
      uint64_t lo, hi, nlo, nhi;
      rank_mtable_prefetch(mtable,&query);
      rank_mtable_fetch(mtable,&query,&lo,&hi);
      naive_range(ix->g->text,n,ix->sa,key,k,&nlo,&nhi);
      if ((nhi>nlo) ? (lo!=nlo || hi!=nhi) : (lo!=hi)) {
        if (errors++ < 5) gt_note("mtable k=%" PRIu64 " code=%" PRIu64 ": [%" PRIu64 ",%" PRIu64 ") naive [%" PRIu64 ",%" PRIu64 ")",
            k,code,lo,hi,nlo,nhi);
      }
    }
  }
  gt_rng_t rng = gt_rng_new(31);
  uint64_t trial;
  for (trial=0;trial<3000;++trial) {
    uint64_t i;
    const uint64_t m = make_key(ix,&rng,key,true) % RANK_MTABLE_SEARCH_DEPTH + 1;
    rank_mquery_t query;
    rank_mquery_new(&query);
    for (i=m;i-->0;) {
      GT_CHECK(!rank_mquery_is_exhausted(&query));
      rank_mquery_add_char(mtable,&query,key[i]);
    }
    GT_CHECK_EQ_INT(rank_mquery_is_exhausted(&query),m==RANK_MTABLE_SEARCH_DEPTH);
    uint64_t lo, hi, nlo, nhi;
    rank_mtable_fetch(mtable,&query,&lo,&hi);
    naive_range(ix->g->text,n,ix->sa,key,m,&nlo,&nhi);
    if ((nhi>nlo) ? (lo!=nlo || hi!=nhi) : (lo!=hi)) ++errors;
  }
  GT_CHECK_EQ_UINT(errors,0);
}

/*
 * Bidirectional search: the backward interval is the range of the pattern,
 * the forward interval the range of its reverse complement
 */
static void check_2query(indexed_t* const ix) {
  const uint64_t n = ix->g->text_length;
  const fm_index_t* const fm_index = ix->fx.archive->fm_index;
  gt_rng_t rng = gt_rng_new(37);
  uint8_t seed[128], pattern[256], rc[256];
  uint64_t trial, errors = 0, steps = 0;
  for (trial=0;trial<1500;++trial) {
    const uint64_t m = make_key(ix,&rng,seed,true);
    // Grow the pattern from a random middle position, extending right
    // (forward) or left (backward) at random, until it covers the key
    uint64_t left = gt_rng_below(&rng,m), right = left; // pattern = seed[left,right)
    fm_2interval_t interval;
    fm_index_2query_init(fm_index,&interval);
    GT_CHECK(interval.backward_lo==0 && interval.backward_hi==n && interval.forward_lo==0 && interval.forward_hi==n);
    const bool precomputed = gt_rng_chance(&rng,1,2);
    while (left > 0 || right < m) {
      const bool forward = (left == 0) || (right < m && gt_rng_chance(&rng,1,2));
      fm_2interval_t next;
      if (forward) {
        const uint8_t c = seed[right++];
        if (precomputed) {
          fm_2erank_elms_t lo_elms, hi_elms;
          fm_index_2query_forward_precompute(fm_index,&interval,&lo_elms,&hi_elms);
          fm_index_2query_precomputed_forward_query(&lo_elms,&hi_elms,&interval,&next,c);
        } else {
          fm_index_2query_forward_query(fm_index,&interval,&next,c);
        }
      } else {
        const uint8_t c = seed[--left];
        if (precomputed) {
          fm_2erank_elms_t lo_elms, hi_elms;
          fm_index_2query_backward_precompute(fm_index,&interval,&lo_elms,&hi_elms);
          fm_index_2query_precomputed_backward_query(&lo_elms,&hi_elms,&interval,&next,c);
        } else {
          fm_index_2query_backward_query(fm_index,&interval,&next,c);
        }
      }
      interval = next;
      const uint64_t len = right-left;
      uint64_t i;
      for (i=0;i<len;++i) {
        pattern[i] = seed[left+i];
        rc[len-1-i] = enc_complement(seed[left+i]);
      }
      uint64_t blo, bhi, flo, fhi;
      naive_range(ix->g->text,n,ix->sa,pattern,len,&blo,&bhi);
      naive_range(ix->g->text,n,ix->sa,rc,len,&flo,&fhi);
      ++steps;
      if (bhi == blo) break; // the pattern no longer occurs
      if (interval.backward_lo!=blo || interval.backward_hi!=bhi ||
          interval.forward_lo!=flo || interval.forward_hi!=fhi) {
        if (errors++ < 5) gt_note("2query len=%" PRIu64 " %s: backward [%" PRIu64 ",%" PRIu64 ") naive [%" PRIu64 ",%" PRIu64 "),"
            " forward [%" PRIu64 ",%" PRIu64 ") naive [%" PRIu64 ",%" PRIu64 ")",len,forward ? "fwd" : "bwd",
            interval.backward_lo,interval.backward_hi,blo,bhi,interval.forward_lo,interval.forward_hi,flo,fhi);
        break;
      }
    }
  }
  GT_CHECK_EQ_UINT(errors,0);
  GT_CHECK(steps > 10000);
}

/*
 * Locator: index positions map to FASTA coordinates and back
 */
static void check_locator(indexed_t* const ix) {
  const genome_t* const g = ix->g;
  locator_t* const locator = ix->fx.archive->locator;
  const uint64_t fwd = g->forward_length;
  uint64_t p, errors = 0, inverse_errors = 0, interval_errors = 0;
  for (p=0;p+2<g->text_length;++p) {
    // Forward position p, or reverse position p (base fwd-2-(p-fwd) of the forward text)
    const bool reverse = p >= fwd;
    const uint64_t f = reverse ? fwd-2-(p-fwd) : p;
    if (g->pos_contig[f] < 0) continue; // separator
    const contig_t* const contig = g->contigs + g->pos_contig[f];
    location_t location;
    locator_map(locator,p,&location);
    const int64_t expected = reverse ? g->pos_offset[f]+1 : g->pos_offset[f];
    if (strcmp(location.tag,contig->name)!=0 || location.strand!=(reverse ? Reverse : Forward) ||
        location.position!=expected || location.bs_strand!=bs_strand_none) {
      if (errors++ < 5) gt_note("locator_map(%" PRIu64 ") = %s:%c:%" PRId64 " expected %s:%c:%" PRId64,
          p,location.tag,location.strand==Forward ? '+' : '-',location.position,
          contig->name,reverse ? '-' : '+',expected);
    }
    const locator_interval_t* const interval = locator_lookup_interval(locator,p);
    if (!(interval->begin_position <= p && p < interval->end_position) ||
        strcmp(locator_interval_get_tag(locator,interval),contig->name)!=0) ++interval_errors;
    if (!reverse) {
      const uint64_t back = locator_inverse_map_position(
          locator,(uint8_t*)contig->name,Forward,bs_strand_none,(uint64_t)g->pos_offset[f]);
      if (back != p) ++inverse_errors;
    }
  }
  GT_CHECK_EQ_UINT(errors,0);
  GT_CHECK_EQ_UINT(interval_errors,0);
  GT_CHECK_EQ_UINT(inverse_errors,0);
  GT_CHECK(locator_get_size(locator) > 0);
}


/*
 * Region profiles (seeding): every strategy's regions carry the exact
 * interval of their substring, plus each strategy's own definition
 */
typedef struct {
  indexed_t* ix;
  mm_slab_t* slab;
  mm_allocator_t* mm_allocator;
  uint64_t errors;
} rp_ctx_t;
static uint64_t rp_count(rp_ctx_t* const ctx,const uint8_t* const key,const uint64_t length) {
  uint64_t lo, hi;
  naive_range(ctx->ix->g->text,ctx->ix->g->text_length,ctx->ix->sa,key,length,&lo,&hi);
  return hi-lo;
}
#define RP_FAIL(ctx,...) do { if ((ctx)->errors++ < 8) gt_note(__VA_ARGS__); } while (0)
// Each region lies in the key and its [lo,hi) is the naive interval of key[begin,end)
static void rp_check_intervals(
    rp_ctx_t* const ctx,const region_profile_t* const rp,
    const uint8_t* const key,const uint64_t key_length,const char* const label) {
  uint64_t i;
  for (i=0;i<rp->num_filtering_regions;++i) {
    const region_search_t* const r = rp->filtering_region + i;
    if (!(r->begin < r->end && r->end <= key_length)) {
      RP_FAIL(ctx,"%s: region %" PRIu64 " [%" PRIu64 ",%" PRIu64 ") outside key of %" PRIu64,label,i,r->begin,r->end,key_length);
      continue;
    }
    uint64_t lo, hi;
    naive_range(ctx->ix->g->text,ctx->ix->g->text_length,ctx->ix->sa,key+r->begin,r->end-r->begin,&lo,&hi);
    const bool ok = (hi > lo) ? (r->lo==lo && r->hi==hi) : (r->lo==r->hi);
    if (!ok) RP_FAIL(ctx,"%s: region [%" PRIu64 ",%" PRIu64 ") interval [%" PRIu64 ",%" PRIu64 ") naive [%" PRIu64 ",%" PRIu64 ")",
        label,r->begin,r->end,r->lo,r->hi,lo,hi);
  }
}
// Scheduling: selected regions get @degree, the others none; totals add up
static void rp_check_schedule(
    rp_ctx_t* const ctx,const region_profile_t* const rp,
    const uint64_t threshold,const uint64_t degree,const char* const label) {
  uint64_t i, selected = 0, total = 0;
  for (i=0;i<rp->num_filtering_regions;++i) {
    const region_search_t* const r = rp->filtering_region + i;
    const uint64_t cand = r->hi - r->lo;
    const uint64_t expected = (cand <= threshold) ? degree : REGION_FILTER_NONE;
    if (r->degree != expected) RP_FAIL(ctx,"%s: region %" PRIu64 " (%" PRIu64 " candidates) degree %" PRIu64 " expected %" PRIu64,
        label,i,cand,r->degree,expected);
    if (cand <= threshold) { ++selected; total += cand; }
  }
  if (rp->num_filtered_regions != selected || rp->total_candidates != total) {
    RP_FAIL(ctx,"%s: %" PRIu64 " regions / %" PRIu64 " candidates scheduled, expected %" PRIu64 " / %" PRIu64,
        label,rp->num_filtered_regions,rp->total_candidates,selected,total);
  }
}
// The whole key as one exact region (scheduled before any strategy-specific rule)
static bool rp_is_exact_match(const region_profile_t* const rp,const uint64_t key_length) {
  return rp->num_filtering_regions==1 && rp->filtering_region[0].begin==0 &&
      rp->filtering_region[0].end==key_length;
}
static void rp_fixed(rp_ctx_t* const ctx,region_profile_t* const rp,const uint8_t* const key,const uint64_t L) {
  static const uint64_t params[][4] = { // length, step, error, max_candidates
    {20,20,0,100}, {12,5,1,3}, {30,40,2,1000}, {11,11,0,0},
  };
  uint64_t p;
  for (p=0;p<sizeof(params)/sizeof(params[0]);++p) {
    const uint64_t len = params[p][0], step = params[p][1], err = params[p][2], maxc = params[p][3];
    region_profile_generate_fixed(rp,ctx->ix->fx.archive->fm_index,key,L,len,step,err,maxc);
    // Regions [k*step, k*step+len) for every k that fits
    uint64_t expected = 0, k;
    for (k=0;k*step+len<=L;++k) {
      if (k < rp->num_filtering_regions &&
          (rp->filtering_region[k].begin!=k*step || rp->filtering_region[k].end!=k*step+len)) {
        RP_FAIL(ctx,"fixed(%" PRIu64 ",%" PRIu64 "): region %" PRIu64 " is [%" PRIu64 ",%" PRIu64 ")",len,step,k,
            rp->filtering_region[k].begin,rp->filtering_region[k].end);
      }
      ++expected;
    }
    if (rp->num_filtering_regions != expected) RP_FAIL(ctx,"fixed: %" PRIu64 " regions, expected %" PRIu64,rp->num_filtering_regions,expected);
    rp_check_intervals(ctx,rp,key,L,"fixed");
    if (!rp_is_exact_match(rp,L)) rp_check_schedule(ctx,rp,maxc,err+1,"fixed");
  }
}
static void rp_factors(rp_ctx_t* const ctx,region_profile_t* const rp,const uint8_t* const key,const uint64_t L) {
  static const uint64_t params[][3] = { {5,0,100}, {3,2,10}, {8,1,1000}, {1,0,5} }; // regions, error, max_candidates
  uint64_t p;
  for (p=0;p<sizeof(params)/sizeof(params[0]);++p) {
    const uint64_t num = params[p][0], err = params[p][1], maxc = params[p][2];
    region_profile_generate_factors(rp,ctx->ix->fx.archive->fm_index,key,L,num,err,maxc);
    // num contiguous regions covering the key; lengths differ by at most one, shorter first
    bool ok = rp->num_filtering_regions==num;
    uint64_t i, pos = 0;
    for (i=0;ok && i<num;++i) {
      const region_search_t* const r = rp->filtering_region + i;
      const uint64_t len = r->end - r->begin;
      ok = r->begin==pos && (len==L/num || len==L/num+1) && (i==0 || len >= rp->filtering_region[i-1].end-rp->filtering_region[i-1].begin);
      pos = r->end;
    }
    ok = ok && pos==L;
    if (!ok) RP_FAIL(ctx,"factors(%" PRIu64 "): not an equal partition of %" PRIu64,num,L);
    rp_check_intervals(ctx,rp,key,L,"factors");
    if (!rp_is_exact_match(rp,L)) rp_check_schedule(ctx,rp,maxc,err+1,"factors");
  }
}
static void rp_cks(rp_ctx_t* const ctx,region_profile_t* const rp,const uint8_t* const key,const uint64_t L) {
  static const uint64_t params[][3] = { {3,20,100}, {2,12,5}, {10,11,1000} }; // regions, length, max_candidates
  uint64_t p;
  for (p=0;p<sizeof(params)/sizeof(params[0]);++p) {
    const uint64_t num = params[p][0], len = params[p][1], maxc = params[p][2];
    region_profile_generate_cks(rp,ctx->ix->fx.archive->fm_index,key,L,len,num,maxc);
    // The fixed partition [k*len,(k+1)*len), of which the (up to) num regions
    // with fewest candidates (and at most maxc) are selected
    if (rp->num_filtering_regions != L/len) RP_FAIL(ctx,"CKS: %" PRIu64 " regions",rp->num_filtering_regions);
    uint64_t i, j, selected = 0, passing = 0, total = 0;
    for (i=0;i<rp->num_filtering_regions;++i) {
      const region_search_t* const r = rp->filtering_region + i;
      if (r->begin!=i*len || r->end!=(i+1)*len) RP_FAIL(ctx,"CKS: region %" PRIu64 " [%" PRIu64 ",%" PRIu64 ")",i,r->begin,r->end);
      if (r->hi-r->lo <= maxc) ++passing;
      if (r->degree==REGION_FILTER_DEGREE_ZERO) {
        ++selected; total += r->hi-r->lo;
        if (r->hi-r->lo > maxc) RP_FAIL(ctx,"CKS: selected a region over the threshold");
        for (j=0;j<rp->num_filtering_regions;++j) {
          const region_search_t* const o = rp->filtering_region + j;
          if (o->degree!=REGION_FILTER_DEGREE_ZERO && (o->hi-o->lo) < (r->hi-r->lo)) {
            RP_FAIL(ctx,"CKS: unselected region %" PRIu64 " has fewer candidates than selected %" PRIu64,j,i);
          }
        }
      } else if (r->degree!=REGION_FILTER_NONE) {
        RP_FAIL(ctx,"CKS: degree %" PRIu64,r->degree);
      }
    }
    if (!rp_is_exact_match(rp,L) && (selected != MIN(num,passing) ||
        rp->num_filtered_regions != selected || rp->total_candidates != total)) {
      RP_FAIL(ctx,"CKS(%" PRIu64 "): %" PRIu64 " selected of %" PRIu64 " passing",num,selected,passing);
    }
    rp_check_intervals(ctx,rp,key,L,"CKS");
  }
}
// Minimum total candidates over num non-overlapping windows of length len
static uint64_t ops_brute_force(rp_ctx_t* const ctx,const uint8_t* const key,const uint64_t L,const uint64_t num,const uint64_t len) {
  // best[k][i]: min total of k windows within key[0,i)
  static uint64_t best[16][256];
  uint64_t k, i;
  for (i=0;i<=L;++i) best[0][i] = 0;
  for (k=1;k<=num;++k) {
    for (i=0;i<=L;++i) {
      uint64_t v = (i>0) ? best[k][i-1] : UINT64_MAX;
      if (i >= len && best[k-1][i-len] != UINT64_MAX) {
        const uint64_t take = best[k-1][i-len] + rp_count(ctx,key+i-len,len);
        if (take < v) v = take;
      }
      best[k][i] = v;
    }
  }
  return best[num][L];
}
static void rp_ops(rp_ctx_t* const ctx,region_profile_t* const rp,const uint8_t* const key,const uint64_t L) {
  static const uint64_t params[][3] = { {3,20,100}, {4,12,10}, {2,30,1000} }; // regions, length, max_candidates
  uint64_t p;
  for (p=0;p<sizeof(params)/sizeof(params[0]);++p) {
    const uint64_t num = params[p][0], len = params[p][1], maxc = params[p][2];
    if (L/num < len) continue;
    region_profile_generate_optimum_fixed(rp,ctx->ix->fx.archive->fm_index,key,L,len,num,maxc);
    // num non-overlapping regions of length len with the minimum total of candidates
    uint64_t i, total = 0;
    if (rp->num_filtering_regions != num) RP_FAIL(ctx,"OPS: %" PRIu64 " regions, expected %" PRIu64,rp->num_filtering_regions,num);
    for (i=0;i<rp->num_filtering_regions;++i) {
      const region_search_t* const r = rp->filtering_region + i;
      if (r->end-r->begin != len) RP_FAIL(ctx,"OPS: region of length %" PRIu64,r->end-r->begin);
      if (i>0 && r->end > rp->filtering_region[i-1].begin) RP_FAIL(ctx,"OPS: overlapping or unordered regions");
      total += r->hi-r->lo;
    }
    const uint64_t optimum = ops_brute_force(ctx,key,L,num,len);
    if (total != optimum) RP_FAIL(ctx,"OPS(%" PRIu64 ",%" PRIu64 "): total %" PRIu64 " optimum %" PRIu64,num,len,total,optimum);
    rp_check_intervals(ctx,rp,key,L,"OPS");
    if (!rp_is_exact_match(rp,L)) rp_check_schedule(ctx,rp,maxc,REGION_FILTER_DEGREE_ZERO,"OPS");
  }
}
// Minimum total occurrences over partitions of the key into num contiguous regions
static uint64_t opp_brute_force(rp_ctx_t* const ctx,const uint8_t* const key,const uint64_t L,const uint64_t num) {
  static uint64_t best[16][256];
  uint64_t k, i, j;
  for (i=0;i<=L;++i) best[0][i] = (i==0) ? 0 : UINT64_MAX;
  for (k=1;k<=num;++k) {
    for (i=0;i<=L;++i) {
      best[k][i] = UINT64_MAX;
      for (j=0;j<i;++j) {
        if (best[k-1][j]==UINT64_MAX) continue;
        bool has_n = false;
        uint64_t t;
        for (t=j;t<i;++t) has_n = has_n || key[t]==ENC_DNA_CHAR_N;
        const uint64_t occ = has_n ? 0 : rp_count(ctx,key+j,i-j);
        if (best[k-1][j]+occ < best[k][i]) best[k][i] = best[k-1][j]+occ;
      }
    }
  }
  return best[num][L];
}
static void rp_opp(rp_ctx_t* const ctx,region_profile_t* const rp,const uint8_t* const key,const uint64_t L) {
  static const uint64_t nums[] = {2,4,6};
  uint64_t p;
  for (p=0;p<sizeof(nums)/sizeof(nums[0]);++p) {
    const uint64_t num = nums[p];
    region_profile_generate_optimum_variable(rp,ctx->ix->fx.archive->fm_index,key,L,num);
    // num contiguous regions covering the key with the minimum total of
    // occurrences (each region's lo is 0, hi its count: see the FIXME in
    // region_profile_optimum.c)
    uint64_t i, pos = 0, total = 0;
    if (rp->num_filtering_regions != num) RP_FAIL(ctx,"OPP: %" PRIu64 " regions",rp->num_filtering_regions);
    for (i=0;i<rp->num_filtering_regions;++i) {
      const region_search_t* const r = rp->filtering_region + i;
      if (r->begin != pos || r->end <= r->begin) RP_FAIL(ctx,"OPP: not a partition at region %" PRIu64,i);
      pos = r->end;
      bool has_n = false;
      uint64_t t;
      for (t=r->begin;t<r->end;++t) has_n = has_n || key[t]==ENC_DNA_CHAR_N;
      const uint64_t occ = has_n ? 0 : rp_count(ctx,key+r->begin,r->end-r->begin);
      if (r->hi-r->lo != occ) RP_FAIL(ctx,"OPP: region [%" PRIu64 ",%" PRIu64 ") count %" PRIu64 " naive %" PRIu64,r->begin,r->end,r->hi-r->lo,occ);
      total += occ;
    }
    if (pos != L) RP_FAIL(ctx,"OPP: partition ends at %" PRIu64 " of %" PRIu64,pos,L);
    const uint64_t optimum = opp_brute_force(ctx,key,L,num);
    if (total != optimum || rp->total_candidates != optimum) {
      RP_FAIL(ctx,"OPP(%" PRIu64 "): total %" PRIu64 " (reported %" PRIu64 ") optimum %" PRIu64,num,total,rp->total_candidates,optimum);
    }
  }
}
static void rp_mem(rp_ctx_t* const ctx,region_profile_t* const rp,const uint8_t* const key,const uint64_t L) {
  region_profile_generate_mem(rp,ctx->ix->fx.archive->fm_index,key,L,20);
  // One region per non-N end position (right to left): the longest
  // occurring substring ending there
  uint64_t i, e = L, num = 0;
  for (i=0;i<rp->num_filtering_regions;++i) {
    while (e > 0 && key[e-1]==ENC_DNA_CHAR_N) --e;
    const region_search_t* const r = rp->filtering_region + i;
    if (r->end != e) { RP_FAIL(ctx,"MEM: region %" PRIu64 " ends at %" PRIu64 " expected %" PRIu64,i,r->end,e); break; }
    const bool left_max = r->begin==0 || key[r->begin-1]==ENC_DNA_CHAR_N || rp_count(ctx,key+r->begin-1,r->end-r->begin+1)==0;
    if (!left_max) RP_FAIL(ctx,"MEM: [%" PRIu64 ",%" PRIu64 ") extends left",r->begin,r->end);
    --e; ++num;
  }
  while (e > 0) { if (key[e-1]!=ENC_DNA_CHAR_N) ++num; --e; }
  if (num != rp->num_filtering_regions) RP_FAIL(ctx,"MEM: %" PRIu64 " regions, expected %" PRIu64,rp->num_filtering_regions,num);
  rp_check_intervals(ctx,rp,key,L,"MEM");
  if (!rp_is_exact_match(rp,L)) rp_check_schedule(ctx,rp,20,REGION_FILTER_DEGREE_ZERO,"MEM");
}
static void rp_smem(rp_ctx_t* const ctx,region_profile_t* const rp,const uint8_t* const key,const uint64_t L) {
  region_profile_generate_smem(rp,ctx->ix->fx.archive->fm_index,key,L,20);
  // Every region occurs, cannot be extended to the left, and every key
  // position outside N is covered by some region
  uint64_t i;
  bool* const covered = calloc(L,sizeof(bool));
  for (i=0;i<rp->num_filtering_regions;++i) {
    const region_search_t* const r = rp->filtering_region + i;
    if (r->hi == r->lo) RP_FAIL(ctx,"SMEM: empty region [%" PRIu64 ",%" PRIu64 ")",r->begin,r->end);
    const bool left_max = r->begin==0 || key[r->begin-1]==ENC_DNA_CHAR_N || rp_count(ctx,key+r->begin-1,r->end-r->begin+1)==0;
    if (!left_max) RP_FAIL(ctx,"SMEM: [%" PRIu64 ",%" PRIu64 ") extends left",r->begin,r->end);
    uint64_t t;
    for (t=r->begin;t<r->end && t<L;++t) covered[t] = true;
  }
  for (i=0;i<L;++i) if (key[i]!=ENC_DNA_CHAR_N && !covered[i]) { RP_FAIL(ctx,"SMEM: position %" PRIu64 " uncovered",i); break; }
  free(covered);
  rp_check_intervals(ctx,rp,key,L,"SMEM");
}
static void rp_adaptive(rp_ctx_t* const ctx,region_profile_t* const rp,const uint8_t* const key,const uint64_t L) {
  static const uint64_t params[][3] = { {20,4,1}, {5,2,2}, {100,8,1}, {1,1,1} }; // threshold, steps, dec. factor
  uint64_t p, z;
  for (p=0;p<sizeof(params)/sizeof(params[0]);++p) {
    for (z=0;z<2;++z) {
      region_profile_model_t model;
      region_profile_model_init(&model);
      model.region_th = params[p][0]; model.max_steps = params[p][1]; model.dec_factor = params[p][2];
      region_profile_generate_adaptive(rp,ctx->ix->fx.archive->fm_index,key,L,&model,ALL,z==1);
      // Disjoint regions from right to left, without N, each with at most
      // threshold candidates (or none), unless the whole key is one region
      uint64_t i;
      for (i=0;i<rp->num_filtering_regions;++i) {
        const region_search_t* const r = rp->filtering_region + i;
        if (i>0 && r->end > rp->filtering_region[i-1].begin) RP_FAIL(ctx,"adaptive: overlapping regions");
        uint64_t t;
        for (t=r->begin;t<r->end && t<L;++t) if (key[t]==ENC_DNA_CHAR_N) { RP_FAIL(ctx,"adaptive: region over N"); break; }
        if (!rp_is_exact_match(rp,L) && r->hi-r->lo > model.region_th) {
          RP_FAIL(ctx,"adaptive(th=%" PRIu64 "): region [%" PRIu64 ",%" PRIu64 ") has %" PRIu64 " candidates",
              model.region_th,r->begin,r->end,r->hi-r->lo);
        }
      }
      rp_check_intervals(ctx,rp,key,L,"adaptive");
      if (!rp_is_exact_match(rp,L)) rp_check_schedule(ctx,rp,ALL,REGION_FILTER_DEGREE_ZERO,"adaptive");
      // An exact match of the whole key is reported as one region
      if (rp_count(ctx,key,L) > 0 && rp_count(ctx,key,L) <= model.region_th && !rp_is_exact_match(rp,L)) {
        bool has_n = false; uint64_t t;
        for (t=0;t<L;++t) has_n = has_n || key[t]==ENC_DNA_CHAR_N;
        if (!has_n && rp->num_filtering_regions==0) RP_FAIL(ctx,"adaptive: occurring key without regions");
      }
    }
  }
  // Limited: regions closed at the threshold or at L/num_regions bases
  region_profile_model_t model;
  region_profile_model_init(&model);
  model.num_regions = 4;
  region_profile_generate_adaptive_limited(rp,ctx->ix->fx.archive->fm_index,key,L,&model);
  uint64_t i;
  for (i=0;i<rp->num_filtering_regions;++i) {
    const region_search_t* const r = rp->filtering_region + i;
    if (i>0 && r->end > rp->filtering_region[i-1].begin) RP_FAIL(ctx,"adaptive-limited: overlapping regions");
    if (r->end-r->begin > L/model.num_regions && r->hi-r->lo > model.region_th && !rp_is_exact_match(rp,L)) {
      RP_FAIL(ctx,"adaptive-limited: long region [%" PRIu64 ",%" PRIu64 ") over threshold",r->begin,r->end);
    }
  }
  rp_check_intervals(ctx,rp,key,L,"adaptive-limited");
}
// Rank-table then BWT queries character by character (right to left)
static void rp_query_character(rp_ctx_t* const ctx,const uint8_t* const key,const uint64_t L) {
  fm_index_t* const fm_index = ctx->ix->fx.archive->fm_index;
  const uint64_t mmd = fm_index->rank_table->min_matching_depth;
  rank_mquery_t query;
  rank_mquery_new(&query);
  uint64_t lo = 0, hi = fm_index_get_length(fm_index), k;
  for (k=1;k<=L;++k) {
    const uint8_t c = key[L-k];
    if (c==ENC_DNA_CHAR_N) break;
    region_profile_query_character(fm_index,&query,&lo,&hi,c);
    if (k < mmd) continue; // below the table's minimum depth the interval is not fetched
    const uint64_t count = rp_count(ctx,key+L-k,k);
    uint64_t nlo, nhi;
    naive_range(ctx->ix->g->text,ctx->ix->g->text_length,ctx->ix->sa,key+L-k,k,&nlo,&nhi);
    if ((count>0) ? (lo!=nlo || hi!=nhi) : (lo!=hi)) {
      RP_FAIL(ctx,"query_character: suffix of %" PRIu64 ": [%" PRIu64 ",%" PRIu64 ") naive [%" PRIu64 ",%" PRIu64 ")",k,lo,hi,nlo,nhi);
      break;
    }
    if (count==0) break;
  }
}
// Mean log4 frequency of 11-mers sampled from the end of the key: each
// sample is 11 consecutive bases followed (to its left) by one more base,
// which is skipped; an N restarts the sample
static void rp_kmer_frequency(rp_ctx_t* const ctx,region_profile_t* const rp,const uint8_t* const key,const uint64_t L) {
  region_profile_compute_kmer_frequency(rp,ctx->ix->fx.archive->fm_index,key,L);
  double sum = 0.0;
  uint64_t samples = 0, run = 0;
  int64_t i;
  for (i=(int64_t)L-1;i>=0;--i) {
    if (key[i]==ENC_DNA_CHAR_N) { run = 0; continue; }
    if (run < RANK_MTABLE_SEARCH_DEPTH) { ++run; continue; }
    sum += log((double)rp_count(ctx,key+i+1,RANK_MTABLE_SEARCH_DEPTH));
    ++samples;
    run = 0;
  }
  if (samples == 0) return; // undefined (0/0)
  const double expected = sum/((double)samples*log(4.0));
  if (fabs(rp->kmer_frequency-expected) > 1e-3*(1.0+fabs(expected))) {
    RP_FAIL(ctx,"kmer_frequency %f expected %f",rp->kmer_frequency,expected);
  }
}
static uint64_t rp_make_key(indexed_t* const ix,gt_rng_t* const rng,uint8_t* const key) {
  const genome_t* const g = ix->g;
  const uint64_t L = gt_rng_range(rng,40,150);
  uint64_t tries, i;
  for (tries=0;tries<10000;++tries) {
    const uint64_t p = gt_rng_below(rng,g->text_length-L-2);
    bool ok = true;
    for (i=0;i<L && ok;++i) { key[i] = g->text[p+i]; ok = key[i] < ENC_DNA_CHAR_SEP; }
    if (ok) break;
  }
  // Sequencing errors, sometimes an N, sometimes a random read
  const uint64_t edits = gt_rng_below(rng,6);
  for (i=0;i<edits;++i) key[gt_rng_below(rng,L)] = (uint8_t)gt_rng_below(rng,4);
  if (gt_rng_chance(rng,1,5)) key[gt_rng_below(rng,L)] = ENC_DNA_CHAR_N;
  if (gt_rng_chance(rng,1,10)) for (i=0;i<L;++i) key[i] = (uint8_t)gt_rng_below(rng,4);
  return L;
}
GT_TEST(region_profile_strategies) {
  indexed_t ix;
  indexed_new(&ix,404,SAMPLING_RATE_NONE,SAMPLING_RATE_4);
  rp_ctx_t ctx;
  ctx.ix = &ix;
  ctx.slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM);
  ctx.mm_allocator = mm_allocator_new(ctx.slab);
  ctx.errors = 0;
  gt_rng_t rng = gt_rng_new(41);
  uint8_t key[256];
  uint64_t trial;
  for (trial=0;trial<150;++trial) {
    const uint64_t L = rp_make_key(&ix,&rng,key);
    mm_allocator_push_state(ctx.mm_allocator);
    region_profile_t rp;
    region_profile_init(&rp,L);
    region_profile_inject_mm(&rp,ctx.mm_allocator);
    GT_CHECK_EQ_UINT(rp.pattern_length,L);
    rp_fixed(&ctx,&rp,key,L);
    rp_factors(&ctx,&rp,key,L);
    rp_cks(&ctx,&rp,key,L);
    rp_ops(&ctx,&rp,key,L);
    if (trial % 3 == 0 && L <= 100) rp_opp(&ctx,&rp,key,L);
    rp_mem(&ctx,&rp,key,L);
    rp_smem(&ctx,&rp,key,L);
    rp_adaptive(&ctx,&rp,key,L);
    rp_query_character(&ctx,key,L);
    rp_kmer_frequency(&ctx,&rp,key,L);
    mm_allocator_pop_state(ctx.mm_allocator);
  }
  GT_CHECK_EQ_UINT(ctx.errors,0);
  mm_allocator_delete(ctx.mm_allocator);
  mm_slab_delete(ctx.slab);
  indexed_delete(&ix);
}
// Helpers that reshape a profile: sort by position, fill gaps, merge small
// regions, error limits, exact-match detection
GT_TEST(region_profile_helpers) {
  mm_slab_t* const slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM);
  mm_allocator_t* const mm_allocator = mm_allocator_new(slab);
  region_profile_t rp;
  region_profile_init(&rp,100);
  region_profile_inject_mm(&rp,mm_allocator);
  region_profile_clear(&rp);
  region_profile_allocate_regions(&rp,4);
  // Regions given out of order: [60,80) [10,20) [30,45)
  const uint64_t spans[][3] = { {60,80,1}, {10,20,1}, {30,45,2} };
  uint64_t i;
  for (i=0;i<3;++i) {
    rp.filtering_region[i].begin = spans[i][0];
    rp.filtering_region[i].end = spans[i][1];
    rp.filtering_region[i].degree = spans[i][2];
  }
  rp.num_filtering_regions = 3;
  GT_CHECK(!region_profile_has_exact_matches(&rp));
  region_profile_sort_by_position(&rp);
  GT_CHECK(rp.filtering_region[0].begin==10 && rp.filtering_region[1].begin==30 && rp.filtering_region[2].begin==60);
  // Fill gaps: the gaps become degree-0 regions; the whole key is covered
  uint8_t key[100];
  memset(key,ENC_DNA_CHAR_A,sizeof(key));
  region_profile_fill_gaps(&rp,key,100,0);
  static const uint64_t filled[][3] = { {0,10,0}, {10,20,1}, {20,30,0}, {30,45,2}, {45,60,0}, {60,80,1}, {80,100,0} };
  GT_REQUIRE(rp.num_filtering_regions == 7);
  for (i=0;i<7;++i) {
    GT_CHECK_EQ_UINT(rp.filtering_region[i].begin,filled[i][0]);
    GT_CHECK_EQ_UINT(rp.filtering_region[i].end,filled[i][1]);
    GT_CHECK_EQ_UINT(rp.filtering_region[i].degree,filled[i][2]);
  }
  // Error limits: min is the region's degree, max leaves room for the others
  region_profile_compute_error_limits(&rp,4,6);
  for (i=0;i<7;++i) {
    GT_CHECK_EQ_UINT(rp.filtering_region[i].min,filled[i][2]);
    GT_CHECK_EQ_UINT(rp.filtering_region[i].max,6-(4-filled[i][2]));
  }
  // Merge degree-0 regions shorter than 12 into a neighbour: [0,10) has no
  // region before it and the next one has a degree, so the two merge into
  // a degree-0 [0,20); [20,30) then joins the degree-0 region before it
  region_profile_merge_small_regions(&rp,12);
  static const uint64_t merged[][3] = { {0,30,0}, {30,45,2}, {45,60,0}, {60,80,1}, {80,100,0} };
  GT_REQUIRE(rp.num_filtering_regions == 5);
  for (i=0;i<5;++i) {
    GT_CHECK_EQ_UINT(rp.filtering_region[i].begin,merged[i][0]);
    GT_CHECK_EQ_UINT(rp.filtering_region[i].end,merged[i][1]);
    GT_CHECK_EQ_UINT(rp.filtering_region[i].degree,merged[i][2]);
  }
  // One region over the whole key with candidates: exact matches
  region_profile_clear(&rp);
  region_profile_allocate_regions(&rp,1);
  rp.filtering_region[0].begin = 0; rp.filtering_region[0].end = 100;
  rp.filtering_region[0].lo = 5; rp.filtering_region[0].hi = 7;
  rp.num_filtering_regions = 1;
  GT_CHECK(region_profile_has_exact_matches(&rp));
  rp.filtering_region[0].hi = 5;
  GT_CHECK(!region_profile_has_exact_matches(&rp));
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * Run-length index (gem-indexer -r): the FM-index is built over the text
 * with runs of up to TEXT_RL_MAX_RUN_LENGTH equal characters collapsed;
 * the archive keeps the full text and samples of RL-to-text positions
 */
static uint64_t naive_rl_encode(
    const uint8_t* const text,const uint64_t n,uint8_t* const rl,uint64_t* const run_start,uint32_t* const acc) {
  uint64_t r = 0, i, run = 0;
  for (i=0;i<n;++i) {
    if (i>0 && text[i]==text[i-1] && run < TEXT_RL_MAX_RUN_LENGTH) {
      ++run;
    } else {
      rl[r] = text[i];
      if (run_start) run_start[r] = i;
      ++r;
      run = 1;
    }
    if (acc) acc[r-1] = (uint32_t)(i+1);
  }
  return r;
}
GT_TEST(rl_encode) {
  gt_rng_t rng = gt_rng_new(53);
  static uint8_t text[4096], rl[4096], naive[4096];
  static uint32_t acc[4096], naive_acc[4096];
  // Hand-worked: AAACGGT
  static const uint8_t hand[] = {0,0,0,1,2,2,3};
  uint64_t rl_length;
  archive_text_rl_encode(hand,7,rl,&rl_length,acc);
  GT_CHECK_EQ_UINT(rl_length,4);
  GT_CHECK(rl[0]==0 && rl[1]==1 && rl[2]==2 && rl[3]==3);
  GT_CHECK(acc[0]==3 && acc[1]==4 && acc[2]==6 && acc[3]==7);
  GT_CHECK_EQ_UINT(archive_text_rl_get_run_length(acc,0),3);
  GT_CHECK_EQ_UINT(archive_text_rl_get_run_length(acc,2),2);
  GT_CHECK_EQ_UINT(archive_text_rl_get_decoded_offset_exl(acc,0),0);
  GT_CHECK_EQ_UINT(archive_text_rl_get_decoded_offset_exl(acc,3),6);
  // A run of 12 is split at TEXT_RL_MAX_RUN_LENGTH
  uint8_t twelve[12];
  memset(twelve,0,sizeof(twelve));
  archive_text_rl_encode(twelve,12,rl,&rl_length,acc);
  GT_CHECK_EQ_UINT(rl_length,2);
  GT_CHECK(acc[0]==TEXT_RL_MAX_RUN_LENGTH && acc[1]==12);
  // Random texts with long runs
  uint64_t trial, errors = 0;
  for (trial=0;trial<500;++trial) {
    const uint64_t n = gt_rng_range(&rng,1,4000);
    uint64_t i = 0;
    while (i < n) {
      const uint8_t c = (uint8_t)gt_rng_below(&rng,5);
      uint64_t run = gt_rng_chance(&rng,1,4) ? gt_rng_range(&rng,1,35) : 1;
      while (run-- && i < n) text[i++] = c;
    }
    archive_text_rl_encode(text,n,rl,&rl_length,acc);
    const uint64_t naive_length = naive_rl_encode(text,n,naive,NULL,naive_acc);
    if (rl_length != naive_length || memcmp(rl,naive,naive_length)!=0 ||
        memcmp(acc,naive_acc,naive_length*sizeof(uint32_t))!=0) {
      ++errors; continue;
    }
    for (i=0;i<rl_length;++i) {
      const uint64_t prev = (i>0) ? naive_acc[i-1] : 0;
      if (archive_text_rl_get_run_length(acc,i) != naive_acc[i]-prev) ++errors;
      if (archive_text_rl_get_decoded_offset_exl(acc,i) != prev) ++errors;
    }
  }
  GT_CHECK_EQ_UINT(errors,0);
}
GT_TEST(index_run_length) {
  genome_t* const g = genome_new(505);
  fixture_t fx;
  fixture_build(&fx,g,archive_dna_full,SAMPLING_RATE_NONE,SAMPLING_RATE_4,true);
  archive_t* const archive = fx.archive;
  const uint64_t n = g->text_length;
  uint8_t* const rl = malloc(n);
  uint64_t* const run_start = malloc(n*sizeof(uint64_t));
  const uint64_t rl_n = naive_rl_encode(g->text,n,rl,run_start,NULL);
  GT_CHECK(archive->text->run_length);
  GT_CHECK_EQ_UINT(archive->text->forward_text_length,g->forward_length);
  GT_CHECK_EQ_UINT(fm_index_get_length(archive->fm_index),rl_n);
  // The archive text is the full (not run-length) text
  const uint8_t* const stored = dna_text_retrieve_sequence(archive->text->enc_text,0);
  uint64_t i, errors = 0;
  for (i=0;i<g->forward_length;++i) errors += (stored[i] != g->text[i]);
  GT_CHECK_EQ_UINT(errors,0);
  // Suffix array of the RL text
  uint64_t* const sa = naive_sa(rl,rl_n);
  for (i=0,errors=0;i<rl_n;++i) errors += (fm_index_decode(archive->fm_index,i) != sa[i]);
  GT_CHECK_EQ_UINT(errors,0);
  // Backward search of RL keys
  gt_rng_t rng = gt_rng_new(59);
  uint64_t trial;
  for (trial=0,errors=0;trial<2000;++trial) {
    uint8_t key[64];
    const uint64_t m = gt_rng_range(&rng,1,40);
    const uint64_t p = gt_rng_below(&rng,rl_n-m);
    bool ok = true;
    for (i=0;i<m && ok;++i) { key[i] = rl[p+i]; ok = key[i] < ENC_DNA_CHAR_SEP; }
    if (!ok) continue;
    uint64_t lo, hi, nlo, nhi;
    fm_index_bsearch(archive->fm_index,key,m,&hi,&lo);
    naive_range(rl,rl_n,sa,key,m,&nlo,&nhi);
    if (lo!=nlo || hi!=nhi) ++errors;
  }
  GT_CHECK_EQ_UINT(errors,0);
  // RL positions translate to the text position where their run starts
  mm_slab_t* const slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM);
  mm_allocator_t* const mm_allocator = mm_allocator_new(slab);
  uint64_t translate_errors = 0, sample_errors = 0, locator_errors = 0;
  for (i=0;i<rl_n;++i) {
    if (i % SAMPLED_RL_SAMPLING_RATE == 0) {
      if (sampled_rl_get_sample(archive->text->sampled_rl,i/SAMPLED_RL_SAMPLING_RATE) != run_start[i]) ++sample_errors;
    }
    if (run_start[i] + TEXT_RL_MAX_RUN_LENGTH*SAMPLED_RL_SAMPLING_RATE < g->forward_length) {
      if (archive_text_rl_position_translate(archive->text,i,mm_allocator) != run_start[i]) {
        if (translate_errors++ < 5) gt_note("translate(%" PRIu64 ")=%" PRIu64 " expected %" PRIu64,
            i,archive_text_rl_position_translate(archive->text,i,mm_allocator),run_start[i]);
      }
    }
    // Locator: the RL interval holding a base's run is the base's contig
    const uint64_t t = run_start[i];
    if (t+2 >= n || g->text[t]==ENC_DNA_CHAR_SEP) continue;
    const uint64_t f = (t >= g->forward_length) ? g->forward_length-2-(t-g->forward_length) : t;
    const locator_interval_t* const interval = locator_lookup_rl_interval(archive->locator,i);
    if (!(interval->rl_begin_position <= i && i < interval->rl_end_position) ||
        strcmp(locator_interval_get_tag(archive->locator,interval),g->contigs[g->pos_contig[f]].name)!=0 ||
        interval->strand != ((t >= g->forward_length) ? Reverse : Forward)) ++locator_errors;
  }
  GT_CHECK_EQ_UINT(sample_errors,0);
  GT_CHECK_EQ_UINT(translate_errors,0);
  GT_CHECK_EQ_UINT(locator_errors,0);
  // Retrieving text with its RL encoding
  for (trial=0,errors=0;trial<500;++trial) {
    const uint64_t length = gt_rng_range(&rng,1,150);
    const uint64_t position = gt_rng_below(&rng,g->forward_length-length+1);
    text_trace_t trace;
    mm_allocator_push_state(mm_allocator);
    archive_text_retrieve(archive->text,position,length,false,true,&trace,mm_allocator);
    uint8_t naive[256];
    uint32_t naive_acc[256];
    const uint64_t naive_length = naive_rl_encode(g->text+position,length,naive,NULL,naive_acc);
    if (trace.rl_text_length != naive_length || memcmp(trace.rl_text,naive,naive_length)!=0 ||
        memcmp(trace.rl_runs_acc,naive_acc,naive_length*sizeof(uint32_t))!=0) ++errors;
    mm_allocator_pop_state(mm_allocator);
  }
  GT_CHECK_EQ_UINT(errors,0);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
  free(sa);
  free(rl);
  free(run_start);
  fixture_destroy(&fx);
  genome_free(g);
}

/*
 * Cases
 */
static void check_all(const uint64_t seed,const sampling_rate_t sa_rate,const sampling_rate_t text_rate) {
  indexed_t ix;
  indexed_new(&ix,seed,sa_rate,text_rate);
  const sampled_sa_t* const ssa = ix.fx.archive->fm_index->sampled_sa;
  // The getters return the rate itself (2^rate); no SA-space sampling is a rate above any index length
  if (sa_rate == SAMPLING_RATE_NONE) {
    GT_CHECK(sampled_sa_get_sa_sampling_rate(ssa) > ix.g->text_length);
  } else {
    GT_CHECK_EQ_UINT(sampled_sa_get_sa_sampling_rate(ssa),1ull<<sa_rate);
  }
  GT_CHECK_EQ_UINT(sampled_sa_get_text_sampling_rate(ssa),1ull<<text_rate);
  check_text(&ix);
  check_bwt(&ix);
  check_decode(&ix);
  check_bsearch(&ix);
  check_rank_mtable(&ix);
  check_2query(&ix);
  check_locator(&ix);
  indexed_delete(&ix);
}
GT_TEST(index_default_sampling) {
  // gem-indexer's CPU defaults: no SA-space sampling, text sampling rate 4
  check_all(101,SAMPLING_RATE_NONE,SAMPLING_RATE_4);
}
GT_TEST(index_sa_and_text_sampling) {
  check_all(202,SAMPLING_RATE_8,SAMPLING_RATE_8);
}
GT_TEST(index_sparse_text_sampling) {
  check_all(303,SAMPLING_RATE_NONE,SAMPLING_RATE_64);
}

GT_CASES(
  GT_CASE(index_default_sampling),
  GT_CASE(index_sa_and_text_sampling),
  GT_CASE(index_sparse_text_sampling),
  GT_CASE(region_profile_strategies),
  GT_CASE(region_profile_helpers),
  GT_CASE(rl_encode),
  GT_CASE(index_run_length)
)
