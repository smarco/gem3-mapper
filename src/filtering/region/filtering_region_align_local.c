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

#include "filtering/region/filtering_region_align_local.h"

/*
 * Profile
 */
#define PROFILE_LEVEL PMED

/*
 * Exclude tiles
 */
void filtering_candidates_align_local_exclude_tiles(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern) {
  // Parameters
  match_scaffold_t* const match_scaffold = &filtering_region->match_scaffold;
  const uint64_t num_alignment_regions = match_scaffold->num_alignment_regions;
  // Check number of tiles
  alignment_t* const alignment = &filtering_region->alignment;
  const uint64_t num_pattern_tiles = pattern->alignment_filters.num_tiles;
  alignment_tile_t* const alignment_tiles = alignment->alignment_tiles;
  if (num_pattern_tiles > 1) {
    // Sort alignment-regions by text-offsets
    match_scaffold_sort_alignment_regions(match_scaffold);
    // Exclude tiles without supporting matching-region
    match_alignment_region_t* match_alignment_region = match_scaffold->alignment_regions;
    uint64_t tile_pos, tile_key_begin = 0, region_pos;
    for (tile_pos=0,region_pos=0;tile_pos<num_pattern_tiles;++tile_pos) {
      alignment_filters_tile_t* const filters_tile = pattern->alignment_filters.tiles + tile_pos;
      const uint64_t tile_key_end = tile_key_begin + filters_tile->tile_length;
      // Skip resolved tiles
      if (alignment_tiles[tile_pos].distance!=ALIGN_DISTANCE_UNKNOWN) continue;
      // Find nearest alignment-region
      while (match_alignment_region!=NULL &&
             match_alignment_region_get_key_end(match_alignment_region) <= tile_key_begin) {
        match_alignment_region = (++region_pos < num_alignment_regions) ? match_alignment_region+1 : NULL;
      }
      // Check overlapping
      if (match_alignment_region==NULL ||
          tile_key_end <= match_alignment_region_get_key_begin(match_alignment_region)) {
        alignment_tiles[tile_pos].distance = ALIGN_DISABLED;
      }
      // Next
      tile_key_begin = tile_key_end;
    }
  }
}
/*
 * Rank local-region
 */
void filtering_candidates_align_local_rank(
    filtering_candidates_t* const filtering_candidates,
    filtering_region_t* const filtering_region,
    pattern_t* const pattern) {
  PROFILE_START(GP_FC_REALIGN_LOCAL_RANK,PROFILE_LEVEL);
  // Parameters
  archive_text_t* const archive_text = filtering_candidates->archive->text;
  text_collection_t* const text_collection = &filtering_candidates->text_collection;
  alignment_tile_t* const alignment_tiles = filtering_region->alignment.alignment_tiles;
  alignment_filters_t* const filters = &pattern->alignment_filters;
  alignment_filters_tile_t* const filters_tiles = filters->tiles;
  const uint64_t num_tiles = filters->num_tiles;
  // Fetch text
  filtering_region_retrieve_text(filtering_region,pattern,archive_text,text_collection);
  text_trace_t* const text_trace = text_collection_get_trace(text_collection,filtering_region->text_trace_offset);
  // Rank based on the number of tiles supporting the alignment
  PROF_ADD_COUNTER(GP_FC_REALIGN_LOCAL_RANK_TILES,num_tiles);
  uint64_t tile_pos, region_rank = 0;
  for (tile_pos=0;tile_pos<num_tiles;++tile_pos) {
    // Scaffold tile
    alignment_tile_t* const alignment_tile = alignment_tiles + tile_pos;
    alignment_filters_tile_t* const filters_tile = filters_tiles + tile_pos;
    if (alignment_tile->distance==ALIGN_DISTANCE_UNKNOWN) {
      // Try k-mer filter
      alignment_tile->distance =
          alignment_verify_levenshtein_kmer_filter(
              alignment_tile,filters_tile,pattern->key,
              text_trace->text_padded,filters->mm_stack);
    }
    // Check tile status
    if (alignment_tile->distance!=ALIGN_DISTANCE_INF) ++region_rank;
  }
  filtering_region->alignment.distance_rank = region_rank;
  PROF_ADD_COUNTER(GP_FC_REALIGN_LOCAL_RANK_TOTAL,region_rank);
  PROF_ADD_COUNTER(GP_FC_REALIGN_LOCAL_RANK_ZERO,region_rank==0?1:0);
  PROFILE_STOP(GP_FC_REALIGN_LOCAL_RANK,PROFILE_LEVEL);
}
