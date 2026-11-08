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
 *   Filtering candidates module provides functions to verify filtering-regions
 *   against its corresponding region of text in the index and compute the
 *   distance of the alignment between both
 */

#include "filtering/candidates/filtering_candidates_verify.h"
#include "filtering/candidates/filtering_candidates_process.h"
#include "filtering/region/filtering_region_verify.h"
#include "align/alignment.h"
#include "align/align_bpm_distance.h"

/*
 * Debug
 */
#define DEBUG_FILTERING_CANDIDATES  GEM_DEEP_DEBUG

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Exact-Extend Candidates
 */
void filtering_candidates_extend_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    matches_t* const matches) {
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  filtering_region_t** const regions = filtering_candidates_get_regions(filtering_candidates);
  const uint64_t num_regions = filtering_candidates_get_num_regions(filtering_candidates);
  // Extend all discarded regions
  uint64_t i;
  for (i=0;i<num_regions;++i) {
    // Fetch region and retrieve Text
    filtering_region_t* const filtering_region = regions[i];
    filtering_region_retrieve_text(filtering_region,
        pattern,archive_text,filtering_candidates->mm_allocator);
    // Exact extend and region chain scaffold
    text_trace_t* const text_trace = &filtering_region->text_trace;
    if (filtering_region->match_scaffold.scaffold_type <= scaffold_none) {
      match_scaffold_region_chain(
          &filtering_region->match_scaffold,pattern,
          text_trace,matches,filtering_candidates->mm_allocator);
    }
  }
  // Sort by coverage
  filtering_candidates_sort_regions_by_scaffold_coverage(filtering_candidates);
}
void filtering_candidates_extend_discarded_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    matches_t* const matches) {
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  filtering_region_t** const regions_discarded = filtering_candidates_get_discarded_regions(filtering_candidates);
  const uint64_t num_regions_discarded = filtering_candidates_get_num_discarded_regions(filtering_candidates);
  // Extend all discarded regions
  uint64_t i;
  for (i=0;i<num_regions_discarded;++i) {
    // Fetch region and retrieve Text
    filtering_region_t* const filtering_region = regions_discarded[i];
    filtering_region_retrieve_text(filtering_region,
        pattern,archive_text,filtering_candidates->mm_allocator);
    // Exact extend and region chain scaffold
    text_trace_t* const text_trace = &filtering_region->text_trace;
    if (filtering_region->match_scaffold.scaffold_type <= scaffold_none) {
      match_scaffold_region_chain(
          &filtering_region->match_scaffold,pattern,
          text_trace,matches,filtering_candidates->mm_allocator);
    }
  }
  // Sort by coverage
  filtering_candidates_sort_discarded_by_scaffold_coverage(filtering_candidates);
}
/*
 * Candidate Verification (1-tile patterns of 1-4 words: BPM of several regions at a time)
 *   Per chunk of regions: prepare every region in order (text, alignment, k-mer
 *   filter), run BPM on the regions still undecided together, then store each
 *   result and sort the regions into accepted and discarded in order. BPM allocates
 *   nothing, so the allocations come in the same order as in the per-region loop.
 */
#define FC_VERIFY_BATCH_CHUNK 64
static void filtering_candidates_verify_candidates_batched(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    const bool kmer_count_filter,
    filtering_region_t** const regions_in,
    filtering_region_t** const regions_out,
    filtering_region_t** const regions_discarded,
    const uint64_t num_filtering_regions,
    uint64_t* const num_regions_accepted,
    uint64_t* const num_regions_discarded) {
  // Parameters
  pattern_tiled_t* const pattern_tiled = &pattern->pattern_tiled;
  const uint64_t max_error = pattern->max_effective_filtering_error;
  // Verify by chunks
  const uint8_t* texts[FC_VERIFY_BATCH_CHUNK];
  uint64_t text_lengths[FC_VERIFY_BATCH_CHUNK], pending[FC_VERIFY_BATCH_CHUNK];
  uint64_t distances[FC_VERIFY_BATCH_CHUNK], columns[FC_VERIFY_BATCH_CHUNK];
  uint64_t num_accepted = 0, num_discarded = 0, chunk_begin, n;
  for (chunk_begin=0;chunk_begin<num_filtering_regions;chunk_begin+=FC_VERIFY_BATCH_CHUNK) {
    const uint64_t chunk_end = MIN(chunk_begin+FC_VERIFY_BATCH_CHUNK,num_filtering_regions);
    // Prepare
    uint64_t num_pending = 0;
    for (n=chunk_begin;n<chunk_end;++n) {
      filtering_region_t* const filtering_region = regions_in[n];
      if (filtering_region_verify_prepare(
          filtering_candidates,filtering_region,kmer_count_filter,pattern)) {
        const alignment_tile_t* const alignment_tile = filtering_region->alignment.alignment_tiles;
        texts[num_pending] = filtering_region->text_trace.text_padded + alignment_tile->text_begin_offset;
        text_lengths[num_pending] = alignment_tile->text_end_offset - alignment_tile->text_begin_offset;
        pending[num_pending++] = n;
      }
    }
    // BPM
    bpm_compute_edit_distance_batch(
        &pattern_tiled->tiles->bpm_pattern_tile,texts,text_lengths,
        num_pending,pattern_tiled->tiles->max_error,distances,columns);
    // Store results and sort the regions
    uint64_t p = 0;
    for (n=chunk_begin;n<chunk_end;++n) {
      filtering_region_t* const filtering_region = regions_in[n];
      if (p<num_pending && pending[p]==n) {
        alignment_verify_edit_bpm_1tile_result(
            &filtering_region->alignment,pattern_tiled,distances[p],columns[p],max_error);
        ++p;
      }
      if (filtering_region_verify_check(filtering_region)) {
        regions_out[num_accepted++] = filtering_region;
      } else {
        regions_discarded[num_discarded++] = filtering_region;
      }
    }
  }
  *num_regions_accepted = num_accepted;
  *num_regions_discarded = num_discarded;
}
/*
 * Candidate Verification
 */
#define FILTERING_CANDIDATES_VERIFY_PREFETCH_MAX 64
static void filtering_candidates_verify_prefetch_text(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t** const regions,
    const uint64_t num_regions) {
  // Positions of a run-length text are translated on retrieval
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  if (archive_text->run_length) return;
  const uint8_t* const text = dna_text_get_text(archive_text->enc_text);
  const uint64_t num_prefetched = MIN(num_regions,FILTERING_CANDIDATES_VERIFY_PREFETCH_MAX);
  uint64_t n;
  for (n=0;n<num_prefetched;++n) {
    filtering_region_t* const filtering_region = regions[n];
    if (filtering_region->alignment.distance_min_bound != ALIGN_DISTANCE_UNKNOWN) continue;
    if (filtering_region->text_trace.text != NULL) continue;
    // The window archive_text_retrieve() reads (a reverse-strand one at its forward projection)
    const uint64_t text_position = filtering_region->text_begin_position;
    const uint64_t text_length = filtering_region->text_end_position - text_position;
    const uint64_t position =
        (text_position < archive_text->forward_text_length || archive_text->explicit_complement) ?
        text_position : archive_text_get_projection(archive_text,text_position,text_length);
    const uintptr_t end = (uintptr_t)(text+position+text_length);
    uintptr_t line;
    for (line=(uintptr_t)(text+position)&~(uintptr_t)63;line<end;line+=64) {
      __builtin_prefetch((const void*)line,0,3);
    }
  }
}
uint64_t filtering_candidates_verify_candidates(
    filtering_candidates_t* const filtering_candidates,
    pattern_t* const pattern,
    matches_t* const matches) {
  // Check number of filtering regions
  const uint64_t num_filtering_regions = filtering_candidates_get_num_regions(filtering_candidates);
  if (num_filtering_regions==0) return 0;
  PROFILE_START(GP_FC_VERIFICATION,PROFILE_LEVEL);
  PROFILE_START(GP_FC_VERIFY_CANDIDATES,PROFILE_LEVEL);
  // Parameters
  search_parameters_t* const search_parameters = filtering_candidates->search_parameters;
  const bool candidate_drop_off = search_parameters->candidate_verification.candidate_local_drop_off;
  const verification_strategy_t verification_strategy = search_parameters->candidate_verification.verification_strategy;
  const bool kmer_count_filter = (verification_strategy >= verification_chained);
  // Candidate drop-off filter
  if (candidate_drop_off) {
    filtering_candidates_extend_candidates(filtering_candidates,pattern,matches);
  }
  // Verify all regions
  filtering_region_t** const regions_in = filtering_candidates_get_regions(filtering_candidates);
  filtering_region_t** const regions_out = regions_in;
  filtering_region_t** const regions_discarded =
      filtering_candidates_reserve_discarded_regions(filtering_candidates,num_filtering_regions);
  // Issue every region's text misses at once rather than one region at a time
  filtering_candidates_verify_prefetch_text(filtering_candidates,regions_in,num_filtering_regions);
  uint64_t n, num_regions_accepted = 0, num_regions_discarded = 0;
  const match_alignment_model_t match_alignment_model = search_parameters->match_alignment_model;
  const uint64_t num_words64 = pattern->pattern_tiled.bpm_pattern.pattern_num_words64;
  const bool batched_bpm = !candidate_drop_off &&
      (match_alignment_model==match_alignment_model_levenshtein ||
       match_alignment_model==match_alignment_model_gap_affine ||
       match_alignment_model==match_alignment_model_none) &&
      pattern->key_length > 0 &&
      pattern->pattern_tiled.num_tiles==1 &&
      num_words64 >= 1 && num_words64 <= BPM_BATCH_MAX_WORDS;
  if (batched_bpm) {
    filtering_candidates_verify_candidates_batched(
        filtering_candidates,pattern,kmer_count_filter,regions_in,regions_out,
        regions_discarded,num_filtering_regions,&num_regions_accepted,&num_regions_discarded);
  } else {
    for (n=0;n<num_filtering_regions;++n) {
      // Verify region
      filtering_region_t* const filtering_region = regions_in[n];
      if (filtering_region_verify(filtering_candidates,filtering_region,kmer_count_filter,pattern)) {
        regions_out[num_regions_accepted++] = filtering_region;
      } else {
        // Heuristic filter
        if (candidate_drop_off && num_regions_accepted==0) {
          const float coverage = filtering_region->match_scaffold.scaffolding_coverage;
          const float key_length = pattern->key_length;
          // Check coverage
          if (coverage/key_length < 1.0) {
            for (;n<num_filtering_regions;++n) { // Discard all regions
              regions_discarded[num_regions_discarded++] = regions_in[n];
            }
            break;
          }
        }
        // Discard
        regions_discarded[num_regions_discarded++] = filtering_region;
      }
    }
  }
  // Update Used
  filtering_candidates_set_num_regions(filtering_candidates,num_regions_accepted);
  filtering_candidates_add_num_discarded_regions(filtering_candidates,num_regions_discarded);
  // DEBUG
  gem_cond_debug_block(DEBUG_FILTERING_CANDIDATES) {
    tab_fprintf(gem_log_get_stream(),"[GEM]>Filtering.Candidates (verify_regions)\n");
    tab_global_inc();
    filtering_candidates_print_regions(gem_log_get_stream(),filtering_candidates,false);
    tab_global_dec();
  }
  PROFILE_STOP(GP_FC_VERIFY_CANDIDATES,PROFILE_LEVEL);
  PROFILE_STOP(GP_FC_VERIFICATION,PROFILE_LEVEL);
  return num_regions_accepted;
}
