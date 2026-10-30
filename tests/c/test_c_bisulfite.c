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
 *   Bisulfite read conversion before searching a bisulfite index: C-to-T
 *   (read 1 of a directional library) or G-to-A (read 2), the end chosen
 *   from the read name (inferred), the mode, or the base composition
 *   (non-stranded: a read with more C than G is taken as G-to-A); and the
 *   restoration of the read as sequenced before output.
 *
 *   Single-end modes '2' and 'interleaved' are not tested at this level:
 *   gem-mapper never passes them to sequence_init (known-bugs/
 *   bisulfite-read-2-ignored.md).
 *
 *   Internal signatures called (a speed iteration may change them; update
 *   this driver with them):
 *     sequence_bisulfite_convert(string_t* input,string_t* backup,const char* BS_table)
 *     sequence_bisulfite_check_cg_depletion_se(string_t*)
 *     sequence_bisulfite_check_cg_depletion_pe(string_t*,string_t*)
 *     sequence_bisulfite_process_se(sequence_t*,bisulfite_read_t)
 *     sequence_bisulfite_process_pe(sequence_t*,sequence_t*,bisulfite_read_t)
 *     sequence_bisulfite_restore_se(sequence_t*), sequence_bisulfite_restore_pe(sequence_t*,sequence_t*)
 *     sequence_init(sequence_t*,bisulfite_read_t,mm_allocator_t*), sequence_set_read (area A)
 *     dna_bisulfite_C2T_table, dna_bisulfite_G2A_table
 *   and reads these fields: sequence_t (read, bs_original_sequence, end_info).
 */

#include "gem3_test.h"
#include "text/sequence.h"
#include "text/sequence_bisulfite.h"
#include "text/dna_text.h"

#define MAX_READ 400

/* Defined in src/text/sequence_bisulfite.c without a prototype in a header */
void sequence_bisulfite_convert(string_t* const input_sequence,string_t* const backup_sequence,const char* const BS_table);
bool sequence_bisulfite_check_cg_depletion_se(string_t* const read_sequence);
bool sequence_bisulfite_check_cg_depletion_pe(string_t* const read_sequence1,string_t* const read_sequence2);

static const char alphabet[] = "ACGTNacgtn";

static void random_read(gt_rng_t* const rng,char* const out,const uint64_t length,const uint64_t c_bias) {
  uint64_t i;
  for (i=0;i<length;++i) {
    // c_bias: 0 none, 1 C-rich, 2 G-rich
    uint64_t k = gt_rng_below(rng,10);
    if (c_bias==1 && gt_rng_chance(rng,1,3)) k = gt_rng_below(rng,2) ? 1 : 6;
    if (c_bias==2 && gt_rng_chance(rng,1,3)) k = gt_rng_below(rng,2) ? 2 : 7;
    out[i] = alphabet[k];
  }
  out[length] = '\0';
}
// Expected conversions (upper-case output; N stays N)
static char c2t(const char c) {
  switch (c) {
    case 'A': case 'a': return 'A';
    case 'C': case 'c': return 'T';
    case 'G': case 'g': return 'G';
    case 'T': case 't': return 'T';
    default: return 'N';
  }
}
static char g2a(const char c) {
  switch (c) {
    case 'A': case 'a': return 'A';
    case 'C': case 'c': return 'C';
    case 'G': case 'g': return 'A';
    case 'T': case 't': return 'T';
    default: return 'N';
  }
}
static void convert(const char* const in,char* const out,char (*f)(const char)) {
  uint64_t i;
  for (i=0;in[i];++i) out[i] = f(in[i]);
  out[i] = '\0';
}
static uint64_t count_ci(const char* const s,const char base) {
  uint64_t i, n = 0;
  for (i=0;s[i];++i) n += (s[i]==base || s[i]==base+('a'-'A'));
  return n;
}
static bool read_equals(sequence_t* const sequence,const char* const expected) {
  const uint64_t n = strlen(expected);
  return sequence_get_length(sequence)==n && strncmp(sequence_get_read(sequence),expected,n)==0;
}
static void set_read(sequence_t* const sequence,const char* const text) {
  sequence_set_read(sequence,(char*)text,strlen(text));
}
static mm_allocator_t* new_allocator(mm_slab_t** const slab) {
  *slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM);
  return mm_allocator_new(*slab);
}

GT_TEST(convert_tables) {
  gt_rng_t rng = gt_rng_new(301);
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  static char read[MAX_READ+1], want[MAX_READ+1];
  uint64_t trial;
  for (trial=0;trial<500;++trial) {
    const uint64_t length = gt_rng_range(&rng,1,MAX_READ);
    random_read(&rng,read,length,0);
    string_t input, backup;
    string_init(&input,10,mm_allocator);
    string_init(&backup,10,mm_allocator);
    string_set_buffer(&input,read,length);
    const bool use_c2t = trial%2==0;
    sequence_bisulfite_convert(&input,&backup,use_c2t ? dna_bisulfite_C2T_table : dna_bisulfite_G2A_table);
    convert(read,want,use_c2t ? c2t : g2a);
    GT_CHECK_EQ_UINT(string_get_length(&input),length);
    GT_CHECK_MSG(strncmp(string_get_buffer(&input),want,length)==0,"trial %" PRIu64 ": converted read differs",trial);
    GT_CHECK_EQ_UINT(string_get_length(&backup),length);
    GT_CHECK_MSG(strncmp(string_get_buffer(&backup),read,length)==0,"trial %" PRIu64 ": backup differs",trial);
  }
  // Hand-worked
  string_t s, b;
  string_init(&s,10,mm_allocator);
  string_init(&b,10,mm_allocator);
  string_set_buffer(&s,"ACGTNacgtn",10);
  sequence_bisulfite_convert(&s,&b,dna_bisulfite_C2T_table);
  GT_CHECK(strncmp(string_get_buffer(&s),"ATGTNATGTN",10)==0);
  string_set_buffer(&s,"ACGTNacgtn",10);
  sequence_bisulfite_convert(&s,&b,dna_bisulfite_G2A_table);
  GT_CHECK(strncmp(string_get_buffer(&s),"ACATNACATN",10)==0);
  GT_CHECK(strncmp(string_get_buffer(&b),"ACGTNacgtn",10)==0);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

GT_TEST(cg_depletion) {
  gt_rng_t rng = gt_rng_new(302);
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  static char r1[MAX_READ+1], r2[MAX_READ+1];
  uint64_t trial;
  for (trial=0;trial<2000;++trial) {
    const uint64_t l1 = gt_rng_range(&rng,1,200), l2 = gt_rng_range(&rng,1,200);
    random_read(&rng,r1,l1,gt_rng_below(&rng,3));
    random_read(&rng,r2,l2,gt_rng_below(&rng,3));
    string_t s1, s2;
    string_init(&s1,10,mm_allocator);
    string_init(&s2,10,mm_allocator);
    string_set_buffer(&s1,r1,l1);
    string_set_buffer(&s2,r2,l2);
    // More C than G: the read looks G-to-A converted
    GT_CHECK_EQ_INT(sequence_bisulfite_check_cg_depletion_se(&s1),count_ci(r1,'C') > count_ci(r1,'G'));
    // Pair: end 1's C plus end 2's G against end 1's G plus end 2's C
    GT_CHECK_EQ_INT(sequence_bisulfite_check_cg_depletion_pe(&s1,&s2),
        count_ci(r1,'C')+count_ci(r2,'G') > count_ci(r1,'G')+count_ci(r2,'C'));
  }
  string_t s;
  string_init(&s,10,mm_allocator);
  string_set_buffer(&s,"ACGT",4);  // tie
  GT_CHECK(!sequence_bisulfite_check_cg_depletion_se(&s));
  string_set_buffer(&s,"cCgT",4);  // case-insensitive: 2 C, 1 G
  GT_CHECK(sequence_bisulfite_check_cg_depletion_se(&s));
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * Single-end: which conversion each mode applies, and the restore
 */
GT_TEST(process_se_modes) {
  gt_rng_t rng = gt_rng_new(303);
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  static char read[MAX_READ+1], want[MAX_READ+1];
  uint64_t trial;
  for (trial=0;trial<2000;++trial) {
    const uint64_t length = gt_rng_range(&rng,1,MAX_READ);
    random_read(&rng,read,length,gt_rng_below(&rng,3));
    static const bisulfite_read_t modes[] = { bisulfite_read_inferred, bisulfite_read_1, bisulfite_non_stranded };
    const bisulfite_read_t mode = modes[trial%3];
    const sequence_end_t end_info = (sequence_end_t)gt_rng_below(&rng,3);
    sequence_t sequence;
    memset(&sequence,0,sizeof(sequence));
    sequence_init(&sequence,mode,mm_allocator);
    set_read(&sequence,read);
    sequence.end_info = end_info;
    sequence_bisulfite_process_se(&sequence,mode);
    // Expected conversion
    sequence_end_t expected_end;
    switch (mode) {
      case bisulfite_read_inferred: expected_end = end_info; break; // from the /1 or /2 suffix
      case bisulfite_read_1: expected_end = paired_end1; break;
      default: expected_end = (count_ci(read,'C') > count_ci(read,'G')) ? paired_end2 : paired_end1; break;
    }
    if (expected_end==paired_end1) convert(read,want,c2t);
    else if (expected_end==paired_end2) convert(read,want,g2a);
    else strcpy(want,read);
    GT_CHECK_MSG(read_equals(&sequence,want),"trial %" PRIu64 " (mode %d, end %d): converted read differs",
        trial,(int)mode,(int)end_info);
    if (mode!=bisulfite_read_1) GT_CHECK_EQ_INT(sequence.bs_sequence_end,expected_end);
    // Restore: the read as sequenced
    sequence_bisulfite_restore_se(&sequence);
    GT_CHECK_MSG(read_equals(&sequence,read),"trial %" PRIu64 ": restored read differs",trial);
  }
  // Mode 2 (as sequence_init sets it up): G-to-A
  sequence_t sequence;
  memset(&sequence,0,sizeof(sequence));
  sequence_init(&sequence,bisulfite_read_2,mm_allocator);
  set_read(&sequence,"ACGTTGCA");
  sequence_bisulfite_process_se(&sequence,bisulfite_read_2);
  GT_CHECK(read_equals(&sequence,"ACATTACA"));
  sequence_bisulfite_restore_se(&sequence);
  GT_CHECK(read_equals(&sequence,"ACGTTGCA"));
  // Inferred without a suffix: no conversion, restore leaves the read alone
  memset(&sequence,0,sizeof(sequence));
  sequence_init(&sequence,bisulfite_read_inferred,mm_allocator);
  set_read(&sequence,"CCCGGG");
  sequence_bisulfite_process_se(&sequence,bisulfite_read_inferred);
  GT_CHECK(read_equals(&sequence,"CCCGGG"));
  GT_CHECK(string_is_null(&sequence.bs_original_sequence));
  sequence_bisulfite_restore_se(&sequence);
  GT_CHECK(read_equals(&sequence,"CCCGGG"));
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * Paired-end: end 1 C-to-T and end 2 G-to-A, swapped in non-stranded mode
 * when the composition says the pair is from the complementary strand
 */
GT_TEST(process_pe_modes) {
  gt_rng_t rng = gt_rng_new(304);
  mm_slab_t* slab;
  mm_allocator_t* const mm_allocator = new_allocator(&slab);
  static char r1[MAX_READ+1], r2[MAX_READ+1], w1[MAX_READ+1], w2[MAX_READ+1];
  uint64_t trial;
  for (trial=0;trial<2000;++trial) {
    const uint64_t l1 = gt_rng_range(&rng,1,200), l2 = gt_rng_range(&rng,1,200);
    random_read(&rng,r1,l1,gt_rng_below(&rng,3));
    random_read(&rng,r2,l2,gt_rng_below(&rng,3));
    const bisulfite_read_t mode = (trial%2) ? bisulfite_non_stranded : bisulfite_read_inferred;
    sequence_t s1, s2;
    memset(&s1,0,sizeof(s1));
    memset(&s2,0,sizeof(s2));
    sequence_init(&s1,mode,mm_allocator);
    sequence_init(&s2,mode,mm_allocator);
    set_read(&s1,r1);
    set_read(&s2,r2);
    sequence_bisulfite_process_pe(&s1,&s2,mode);
    const bool flip = mode==bisulfite_non_stranded &&
        count_ci(r1,'C')+count_ci(r2,'G') > count_ci(r1,'G')+count_ci(r2,'C');
    convert(r1,w1,flip ? g2a : c2t);
    convert(r2,w2,flip ? c2t : g2a);
    GT_CHECK_MSG(read_equals(&s1,w1) && read_equals(&s2,w2),"trial %" PRIu64 " (flip %d): converted pair differs",trial,(int)flip);
    sequence_bisulfite_restore_pe(&s1,&s2);
    GT_CHECK_MSG(read_equals(&s1,r1) && read_equals(&s2,r2),"trial %" PRIu64 ": restored pair differs",trial);
  }
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

GT_CASES(
  GT_CASE(convert_tables),
  GT_CASE(cg_depletion),
  GT_CASE(process_se_modes),
  GT_CASE(process_pe_modes)
)
