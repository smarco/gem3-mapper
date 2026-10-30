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
 *   Area D: the slab and the stack-like allocator every search uses.
 *
 *   Properties checked (not addresses, segment counts or reuse order):
 *   memory handed out is writable, never overlaps another live block, keeps
 *   its contents until freed, is 8-byte aligned and zeroed when asked;
 *   blocks may be freed in any order; a pushed state, popped, releases what
 *   was allocated after it without touching what was allocated before; the
 *   allocator survives clear and many rounds that span several segments.
 *
 *   Internal signatures called (a speed iteration may change them; update
 *   this driver with them):
 *     mm_slab_t* mm_slab_new_(uint64_t slab_unit_size,uint64_t slab_group_size,uint64_t max_memory)
 *     mm_slab_unit_t* mm_slab_request(mm_slab_t*); void mm_slab_return(mm_slab_t*,mm_slab_unit_t*)
 *     uint64_t mm_slab_get_slab_size(mm_slab_t*); void mm_slab_delete(mm_slab_t*)
 *     mm_allocator_t* mm_allocator_new(mm_slab_t*); mm_allocator_clear/delete(mm_allocator_t*)
 *     void* mm_allocator_allocate(mm_allocator_t*,uint64_t num_bytes,bool zero_mem)
 *       (through mm_allocator_malloc/mm_allocator_calloc)
 *     void mm_allocator_free(mm_allocator_t*,void*)
 *     mm_allocator_push_memory_state/mm_allocator_pop_memory_state(mm_allocator_t*)
 */

#include "gem3_test.h"
#include "system/mm_slab.h"
#include "system/mm_allocator.h"

typedef struct { uint8_t* p; uint64_t n; uint8_t tag; } block_t;

static void fill(block_t* const b) {
  uint64_t i;
  for (i=0;i<b->n;++i) b->p[i] = (uint8_t)(b->tag + i*31);
}
static bool intact(const block_t* const b) {
  uint64_t i;
  for (i=0;i<b->n;++i) if (b->p[i] != (uint8_t)(b->tag + i*31)) return false;
  return true;
}
static bool overlaps(const block_t* const a,const block_t* const b) {
  return a->p < b->p + b->n && b->p < a->p + a->n;
}

GT_TEST(slab_units_are_disjoint_and_reusable) {
  const uint64_t unit = 64*1024;
  mm_slab_t* const slab = mm_slab_new_(unit,4*unit,MM_UNLIMITED_MEM);
  GT_CHECK_EQ_UINT(mm_slab_get_slab_size(slab),unit);
  enum { N = 40 };
  mm_slab_unit_t* units[N];
  block_t blocks[N];
  int i, j;
  for (i=0;i<N;++i) {
    units[i] = mm_slab_request(slab);
    blocks[i].p = units[i]->memory;
    blocks[i].n = unit;
    blocks[i].tag = (uint8_t)i;
    fill(blocks+i);
  }
  for (i=0;i<N;++i) {
    GT_CHECK_MSG(intact(blocks+i),"unit %d overwritten",i);
    for (j=0;j<i;++j) GT_CHECK_MSG(!overlaps(blocks+i,blocks+j),"units %d and %d overlap",i,j);
  }
  for (i=0;i<N;i+=2) mm_slab_return(slab,units[i]);
  for (i=0;i<N;i+=2) {
    units[i] = mm_slab_request(slab);
    blocks[i].p = units[i]->memory;
    blocks[i].tag = (uint8_t)(100+i);
    fill(blocks+i);
  }
  for (i=0;i<N;++i) {
    GT_CHECK(intact(blocks+i));
    for (j=0;j<i;++j) GT_CHECK(!overlaps(blocks+i,blocks+j));
  }
  for (i=0;i<N;++i) mm_slab_return(slab,units[i]);
  mm_slab_delete(slab);
}

GT_TEST(allocator_blocks_are_disjoint_and_kept) {
  mm_slab_t* const slab = mm_slab_new_(BUFFER_SIZE_1M,BUFFER_SIZE_1M,MM_UNLIMITED_MEM);
  mm_allocator_t* const mma = mm_allocator_new(slab);
  gt_rng_t rng = gt_rng_new(51);
  enum { N = 3000 };
  static block_t live[N];
  int n = 0, round, i, j;
  for (round=0;round<6;++round) {
    // Allocate (some zeroed), spanning several 1 MB segments
    while (n < N) {
      block_t b;
      b.n = 1 + gt_rng_below(&rng, gt_rng_below(&rng,10)==0 ? 20000 : 300);
      const bool zero = gt_rng_below(&rng,3)==0;
      b.p = zero ? mm_allocator_calloc(mma,b.n,uint8_t,true) : mm_allocator_malloc(mma,b.n);
      GT_REQUIRE(b.p!=NULL);
      GT_CHECK_MSG(((uintptr_t)b.p % 8)==0,"block of %" PRIu64 " bytes not 8-aligned",b.n);
      if (zero) {
        uint64_t k;
        for (k=0;k<b.n;++k) if (b.p[k]) { GT_CHECK_MSG(false,"calloc byte %" PRIu64 " not zero",k); break; }
      }
      b.tag = (uint8_t)gt_rng_below(&rng,256);
      fill(&b);
      live[n++] = b;
    }
    for (i=0;i<n;++i) GT_REQUIRE_MSG(intact(live+i),"round %d: block %d overwritten",round,i);
    // Sorted-by-address disjointness check (insertion into a sorted copy)
    static block_t sorted[N];
    for (i=0;i<n;++i) {
      j = i;
      while (j>0 && sorted[j-1].p > live[i].p) { sorted[j] = sorted[j-1]; --j; }
      sorted[j] = live[i];
    }
    for (i=1;i<n;++i) GT_REQUIRE_MSG(!overlaps(sorted+i-1,sorted+i),"round %d: blocks overlap",round);
    // Free a random half, in random order
    for (i=n-1;i>0;--i) {
      j = (int)gt_rng_below(&rng,(uint64_t)i+1);
      block_t t = live[i]; live[i] = live[j]; live[j] = t;
    }
    const int keep = n/2;
    for (i=keep;i<n;++i) mm_allocator_free(mma,live[i].p);
    n = keep;
    for (i=0;i<n;++i) GT_REQUIRE_MSG(intact(live+i),"round %d: survivor %d overwritten after frees",round,i);
  }
  for (i=0;i<n;++i) mm_allocator_free(mma,live[i].p);
  // After clear, the allocator serves memory again
  mm_allocator_clear(mma);
  block_t b = { mm_allocator_malloc(mma,1000), 1000, 7 };
  fill(&b);
  GT_CHECK(intact(&b));
  mm_allocator_delete(mma);
  mm_slab_delete(slab);
}

GT_TEST(allocator_push_pop_releases_later_blocks) {
  mm_slab_t* const slab = mm_slab_new_(BUFFER_SIZE_1M,BUFFER_SIZE_1M,MM_UNLIMITED_MEM);
  mm_allocator_t* const mma = mm_allocator_new(slab);
  gt_rng_t rng = gt_rng_new(61);
  enum { K = 200 };
  block_t base[K];
  int i, round;
  for (i=0;i<K;++i) {
    base[i].n = 1 + gt_rng_below(&rng,500);
    base[i].p = mm_allocator_malloc(mma,base[i].n);
    base[i].tag = (uint8_t)i;
    fill(base+i);
  }
  // Many rounds of push, allocate a lot (several segments), pop: the
  // blocks allocated before the push are never touched
  for (round=0;round<40;++round) {
    mm_allocator_push_state(mma);
    static block_t tmp[2000];
    const int m = 200 + (int)gt_rng_below(&rng,1800);
    for (i=0;i<m;++i) {
      tmp[i].n = 1 + gt_rng_below(&rng,4000);
      tmp[i].p = mm_allocator_malloc(mma,tmp[i].n);
      tmp[i].tag = (uint8_t)(round+i);
      fill(tmp+i);
      int k;
      for (k=0;k<K;k+=37) GT_REQUIRE(!overlaps(tmp+i,base+k));
    }
    // A nested state
    mm_allocator_push_state(mma);
    block_t inner = { mm_allocator_malloc(mma,3000), 3000, 99 };
    fill(&inner);
    mm_allocator_pop_state(mma);
    for (i=0;i<m;++i) GT_REQUIRE_MSG(intact(tmp+i),"round %d: block %d lost to a nested pop",round,i);
    mm_allocator_pop_state(mma);
    for (i=0;i<K;++i) GT_REQUIRE_MSG(intact(base+i),"round %d: base block %d overwritten",round,i);
  }
  for (i=0;i<K;++i) mm_allocator_free(mma,base[i].p);
  mm_allocator_delete(mma);
  mm_slab_delete(slab);
}

GT_CASES(
  GT_CASE(slab_units_are_disjoint_and_reusable),
  GT_CASE(allocator_blocks_are_disjoint_and_kept),
  GT_CASE(allocator_push_pop_releases_later_blocks)
)
