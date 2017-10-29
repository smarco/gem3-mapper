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
 *   Mapper-profile module provides functions to display the
 *   GEM-profiling counters involved in the mapping search
 */

#include "mapper/mapper_profile_search.h"
#include "mapper/mapper_profile_misc.h"

#ifdef GEM_PROFILE /* GEM_PROFILE ENABLED */
/*
 * Region Profile
 */
void mapper_profile_print_region_profile(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Region.Profile\n");
  tab_fprintf(stream,"  --> Num.Profiles             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILES_COMPUTED),NULL,"    ",true);
  tab_fprintf(stream,"  --> Num.Regions              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_NUM_REGIONS),NULL,"    ",true);
  tab_fprintf(stream,"    --> Num.Regions.Filtered   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_NUM_REGIONS_FILTERED),
                       PROF_GET_COUNTER(GP_REGION_PROFILE_NUM_REGIONS),"    ",true);
  tab_fprintf(stream,"  --> Region.length            ");
  SAMPLER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_REGION_LENGTH),NULL,"nt  ");
  tab_fprintf(stream,"  --> Region.candidates        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_REGION_CANDIDATES),NULL,"cand",true);
  tab_fprintf(stream,"  --> Read.candidates          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_TOTAL_CANDIDATES),NULL,"cand",true);
}
/*
 * Candidates Generation
 */
void mapper_profile_print_candidate_generation(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Generate.Candidates\n");
  tab_fprintf(stream,"  => TIME.Generate.Candidates         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_GENERATE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Dynamic.Filtering         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_GENERATE_CANDIDATES_DYNAMIC_FILTERING),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> Generate.Candidates\n");
  tab_fprintf(stream,"    --> Regions.Elegible              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS),
                       PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS),"regions",true);
  tab_fprintf(stream,"      --> Regions.Processed           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_PROCESSED),
                       PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS),"regions",true);
  tab_fprintf(stream,"      --> Regions.Skipped             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_SKIPPED),
                       PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS),"regions",true);
  tab_fprintf(stream,"    --> Candidates.Generated          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),
                       PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"pos",true);
}
/*
 * Candidate Verification
 */
void mapper_profile_print_candidate_verification(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Candidate.Verification\n");
  // Verifying
  tab_fprintf(stream,"  => TIME.Process.Candidates                      ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_PROCESS_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Decode.Positions                      ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_DECODE_POSITIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Compose.Regions                       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_COMPOSE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  => TIME.Verifying                               ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_VERIFICATION),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Verify.Candidates                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_VERIFY_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Verify.Candidate.Region             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_VERIFY_CANDIDATES_REGION),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Kmer.Counting                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_KMER_COUNTER_FILTER),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.BPM.Distance                      ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_DISTANCE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> Filtering.Candidates\n");
  tab_fprintf(stream,"    --> Candidate.Positions                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"positions",true);
  tab_fprintf(stream,"      --> Decode.LF.Dist                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FMIDX_LOOKUP_DIST),NULL,"lf       ",true);
  tab_fprintf(stream,"    --> Candidate.Regions                         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Duplicates.Removed                      ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGIONS_DUPLICATED),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Length                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),NULL,"nt       ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Alignment.Regions     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGION_ALIGNMENT_REGIONS_TOTAL),NULL,"regions  ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Alignment.Coverage    ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGION_ALIGNMENT_COVERAGE),"");
  tab_fprintf(stream,"      --> Alignment.Tiles                         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_TILES),PROF_GET_COUNTER(GP_CANDIDATE_TILES),"tiles    ",true);
  tab_fprintf(stream,"      --> Kmer.Counting\n");
  tab_fprintf(stream,"        --> Kmer.Counting.Candidates              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_KMER_COUNTER_FILTER_CANDIDATES),PROF_GET_COUNTER(GP_FC_KMER_COUNTER_FILTER_CANDIDATES),"         ",true);
  tab_fprintf(stream,"        --> Kmer.Counting.Discarded               ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_KMER_COUNTER_FILTER_DISCARDED),PROF_GET_COUNTER(GP_FC_KMER_COUNTER_FILTER_CANDIDATES),"         ",true);
  tab_fprintf(stream,"        --> Kmer.Counting.Accepted                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_KMER_COUNTER_FILTER_ACCEPTED),PROF_GET_COUNTER(GP_FC_KMER_COUNTER_FILTER_CANDIDATES),"         ",true);
  tab_fprintf(stream,"        --> Kmer.Counting.Accepted.NA             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_KMER_COUNTER_FILTER_NA),PROF_GET_COUNTER(GP_FC_KMER_COUNTER_FILTER_CANDIDATES),"         ",true);
  tab_fprintf(stream,"      --> BPM\n");
  tab_fprintf(stream,"        --> BPM.Candidates                        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_BPM_FILTER_CANDIDATES),PROF_GET_COUNTER(GP_FC_BPM_FILTER_CANDIDATES),"         ",true);
  tab_fprintf(stream,"        --> BPM.Discarded                         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_BPM_FILTER_DISCARDED),PROF_GET_COUNTER(GP_FC_BPM_FILTER_CANDIDATES),"         ",true);
  tab_fprintf(stream,"        --> BPM.Accepted                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_BPM_FILTER_ACCEPTED),PROF_GET_COUNTER(GP_FC_BPM_FILTER_CANDIDATES),"         ",true);
  tab_fprintf(stream,"        --> BPM.Tiles.Computed                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BMP_DISTANCE_NUM_TILES_VERIFIED),PROF_GET_COUNTER(GP_CANDIDATE_TILES),"tiles    ",true);
  tab_fprintf(stream,"          --> BPM.Tile.Tall                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BPM_DISTANCE_KEY_LENGTH),NULL,"nt       ",true);
  tab_fprintf(stream,"          --> BPM.Tile.Length                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BPM_DISTANCE_TEXT_LENGTH),NULL,"nt       ",true);
  tab_fprintf(stream,"          --> BPM.Cells.Computed                  ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BPM_DISTANCE_CELLS),NULL,"cells    ",true);
  tab_fprintf(stream,"          --> BPM.Quick.Abandon                   ");
    COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BPM_DISTANCE_QUICK_ABANDON),PROF_GET_COUNTER(GP_BMP_DISTANCE_NUM_TILES_VERIFIED),"         ",true);
  tab_fprintf(stream,"      --> Alignment.Discarded                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_DISCARDED_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Alignment.Accepted                      ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),"regions  ",true);
}
/*
 * Candidate realign
 */
void mapper_profile_print_candidate_realign(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Candidate.Realign\n");
  tab_fprintf(stream,"  => TIME.Realign                                  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Realign.Exact                          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHES_ALIGN_EXACT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Scaffold.Alignment                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCH_SCAFFOLD_ALIGNMENT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Realign.SWG                            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHES_ALIGN_SWG),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Realign.SWG.Core.Banded              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_SWG_ALIGN_BANDED),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> Realign.Regions\n");
  tab_fprintf(stream,"    --> Candidate.Positions                        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"positions",true);
  tab_fprintf(stream,"    --> Candidate.Regions                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"    --> Accepted.Regions                           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Aligned.Pruned         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_SELECT_PRUNE_HIT),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Aligned.Cached         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_CACHE_SEARCH_HIT),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Scaffold                                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_ALIGNMENT_ADAPTIVE),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"    --> Aligned.Regions                            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ALIGNED_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Aligned.Exact                            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ALIGNED_EXACT),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Aligned.Inexact                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ALIGNED_INEXACT),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Accepted               ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ALIGNED_ACCEPTED),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Discarded              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ALIGNED_DISCARDED),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"        --> Discarded.unscaffolded                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ALIGNED_DISCARDED_SCAFFOLD),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"        --> Discarded.SWG.unaligned                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ALIGNED_DISCARDED_SWG),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"        --> Discarded.SWG.minThreshold             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ALIGNED_DISCARDED_SWG_THRESHOLD),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"        --> Discarded.minMatchingBases             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ALIGNED_DISCARDED_MATCHING_BASES),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"    --> Matches\n");
  tab_fprintf(stream,"      --> Matches.Added                            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCHES_MAPS_ADDED),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Matches.Dups                             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCHES_MAPS_DUP),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"  |> Realign.Length\n");
  tab_fprintf(stream,"    --> Candidate.Regions.Length                   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),"nt       ",true);
  tab_fprintf(stream,"    --> Aligned.Regions.Length                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ALIGNED_REGIONS_LENGTH),PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),"nt       ",true);
  tab_fprintf(stream,"    --> SWG.Banded.Aligned.Length                  ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SWG_ALIGN_BANDED_LENGTH),PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),"nt       ",true);
  tab_fprintf(stream,"    --> SWG.Banded.Aligned.Cells.Computed          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SWG_ALIGN_BANDED_CELLS),NULL,"cells    ",true);
}
void mapper_profile_print_candidate_realign_local(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Candidate.Realign\n");
  tab_fprintf(stream,"  => TIME.Realign.Local                            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_LOCAL_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Rank.Regions                           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_LOCAL_RANK),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> Realign.Local\n");
  tab_fprintf(stream,"    --> Candidate.Regions                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGION_LOCAL),
                       PROF_GET_COUNTER(GP_CANDIDATE_REGION_LOCAL),"regions  ",true);
  tab_fprintf(stream,"    --> Candidate.Rank\n");
  tab_fprintf(stream,"      --> Candidate.Tiles                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_REALIGN_LOCAL_RANK_TILES),NULL,"tiles    ",true);
  tab_fprintf(stream,"      --> Candidate.Tiles.Supporting               ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_REALIGN_LOCAL_RANK_TOTAL),
                       PROF_GET_COUNTER(GP_FC_REALIGN_LOCAL_RANK_TILES),"tiles    ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Rank.Zero              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_REALIGN_LOCAL_RANK_ZERO),
                       PROF_GET_COUNTER(GP_CANDIDATE_REGION_LOCAL),"regions  ",true);
  tab_fprintf(stream,"    --> Candidate.Regions.Aligned.Local            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGION_LOCAL_ALIGNED),
                       PROF_GET_COUNTER(GP_CANDIDATE_REGION_LOCAL),"regions  ",true);
}
/*
 * Candidate Scaffold
 */
void mapper_profile_print_candidate_scaffold(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Candidate.Scaffold\n");
  tab_fprintf(stream,"  => TIME.Scaffold.Alignment                       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCH_SCAFFOLD_ALIGNMENT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Scaffold.Chain                         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCH_SCAFFOLD_CHAIN_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Scaffold.Edit                          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCH_SCAFFOLD_EDIT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Scaffold.Edit.BPM.Matrix             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_COMPUTE_MATRIX),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> Scaffold                                      ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_ALIGNMENT_ADAPTIVE),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"    --> Scaffold.Chain                             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_CHAIN_REGIONS_SCAFFOLDS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Scaffold.Alignment.Regions               ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_ALIGNMENT_REGIONS),NULL,"regions  ",true);
  tab_fprintf(stream,"      --> Scaffold.Alignment.Coverage              ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_ALIGNMENT_COVERAGE),"");
  tab_fprintf(stream,"      --> Scaffold.Chain.Coverage                  ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_CHAIN_REGIONS_COVERAGE),"");
  tab_fprintf(stream,"      --> Scaffold.Edit                            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_SCAFFOLDS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"        --> Scaffold.Edit.Coverage                 ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_COVERAGE),"");
  tab_fprintf(stream,"        --> Scaffold.Edit.Tiles\n");
  tab_fprintf(stream,"          --> Scaffold.Edit.Tiles.Total            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_TOTAL),
                       PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_TOTAL),"tiles    ",true);
  tab_fprintf(stream,"          --> Scaffold.Edit.Tiles.Skipped          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_SKIPPED),
                       PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_TOTAL),"tiles    ",true);
  tab_fprintf(stream,"          --> Scaffold.Edit.Tiles.Align            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_ALIGN),
                       PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_TOTAL),"tiles    ",true);
  tab_fprintf(stream,"          --> Scaffold.Edit.Tiles.DistanceBound    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_DISTANCE_BOUND),NULL,"error    ",true);
  tab_fprintf(stream,"        --> Scaffold.Edit.Cells.Computed           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_CELLS),NULL,"cells    ",true);
}
/*
 * Candidate misc
 */
void mapper_profile_print_candidate_misc(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Candidate.Misc\n");
  tab_fprintf(stream,"  |> Sorting\n");
  tab_fprintf(stream,"    --> Sort.Candidate.by.Position          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_SORT_BY_POSITION),NULL,"         ",true);
  tab_fprintf(stream,"    --> Sort.Candidate.by.Distance          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_SORT_BY_ALIGN_DIST),NULL,"         ",true);
  tab_fprintf(stream,"    --> Sort.Candidate.by.Coverage          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_SORT_BY_COVERAGE),NULL,"         ",true);
}
/*
 * Neighborhood Search
 */
void mapper_profile_print_neighborhood_search(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.NeighborhoodSearch\n");
  tab_fprintf(stream,"  => TIME.NeighborhoodSearch                       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_NEIGHBORHOOD_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.NS.Generation                          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_NS_GENERATION),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.NS.Verification                        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_NS_VERIFICATION),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.NS.Align                               ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_NS_ALIGN),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> NS\n");
  tab_fprintf(stream,"  --> NS.Nodes                                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_NS_NODES),PROF_GET_COUNTER(GP_NS_NODES),"nodes",true);
  tab_fprintf(stream,"    --> NS.Nodes.Success                           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_NS_NODES_SUCCESS),PROF_GET_COUNTER(GP_NS_NODES),"nodes",true);
  tab_fprintf(stream,"    --> NS.Nodes.Fail                              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_NS_NODES_FAIL),PROF_GET_COUNTER(GP_NS_NODES),"nodes",true);
  tab_fprintf(stream,"    --> NS.Nodes.Prevented.MaxError                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_NS_NODES_DP_PREVENT),NULL,"nodes",true);
  tab_fprintf(stream,"  --> NS.Depth                                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_NS_SEARCH_DEPTH),NULL,"     ",true);
  tab_fprintf(stream,"  --> NS.Cells.Computed                            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_NS_DP_CELLS_COMPUTED),NULL,"     ",true);
  tab_fprintf(stream,"  --> NS.Candidates.Branchwise                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_NS_BRANCH_CANDIDATES_GENERATED),
                       PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"     ",true);
//  tab_fprintf(stream,"  --> NS.Candidates.Operation                      ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_NS_OPERATION_CANDIDATES_GENERATED),
//                       PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"     ",true);
  tab_fprintf(stream,"  --> NS.Candidates.TotalSearch                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_NS_SEARCH_CANDIDATES_GENERATED),
                       PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"     ",true);
}
/*
 * Approximate Search
 */
void mapper_profile_print_approximate_search(FILE* const stream) {
  fprintf(stream,    "[GEM]>Profile.Approximate.Search.Stages\n");
  tab_fprintf(stream,"  => TIME.Approximate.Search.Building.Blocks\n");
  tab_fprintf(stream,"    => TIME.Approximate.Search                ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_MAIN),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Neighborhood.Search             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_NEIGHBORHOOD_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Region.Profile                  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_REGION_PROFILE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Generate.Candidates             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_GENERATE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Process.Candidates              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_PROCESS_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Verifying                       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_VERIFICATION),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Realign                         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.NS.Generation                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_NS_GENERATION),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Realign.Local                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_LOCAL_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  => TIME.Approximate.Search.Stages\n");
  tab_fprintf(stream,"    => TIME.Approximate.Search                ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_MAIN),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Filtering.Adaptive              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTERING_ADATIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Neighborhood.Search             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_NEIGHBORHOOD_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Local.Alignment                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_LOCAL_ALIGN),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> Approximate.Search.Stages\n");
  //
  tab_fprintf(stream,"    --> Approximate.Search                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MAPPER_NUM_READS),
      PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"reads",true);
  tab_fprintf(stream,"      --> Filtering.Filtering.Adaptive        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_FILTERING_ADATIVE_CALL),
      PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"reads",true);
  tab_fprintf(stream,"      --> Neighborhood.Search                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_CALL),
      PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"reads",true);
  tab_fprintf(stream,"      --> Local.Alignment                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_LOCAL_ALIGN_CALL),
      PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"reads",true);
  //
  tab_fprintf(stream,"    --> MCS.Approximate.Search                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_MCS),
      PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"mcs  ",true);
  tab_fprintf(stream,"      --> MCS.Filtering.Filtering.Adaptive    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_FILTERING_ADATIVE_MCS),
      PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"mcs  ",true);
  tab_fprintf(stream,"      --> MCS.Neighborhood.Search             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_NEIGHBORHOOD_SEARCH_MCS),
      PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"mcs  ",true);
  //
  tab_fprintf(stream,"    --> Accuracy.Control                      ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_CALLS),
      PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_CALLS),"calls",true);
  tab_fprintf(stream,"      --> Accuracy.Hit                        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_HIT),
      PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_CALLS),"calls",true);
  tab_fprintf(stream,"      --> Accuracy.Miss                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_MISS),
      PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_CALLS),"calls",true);
  tab_fprintf(stream,"        --> Accuracy.Miss.Unmapped            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_UNMAPPED),
      PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_CALLS),"calls",true);
  tab_fprintf(stream,"        --> Accuracy.Miss.Tie                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_TIE),
      PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_CALLS),"calls",true);
  tab_fprintf(stream,"        --> Accuracy.Miss.MMaps               ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_MMAP),
      PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_CALLS),"calls",true);
  tab_fprintf(stream,"        --> Accuracy.Miss.Unique              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_UNIQUE),
      PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_CALLS),"calls",true);
  tab_fprintf(stream,"        --> Accuracy.Miss.Rescued.MaxMatches  ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_MAX_MATCHES),
      PROF_GET_COUNTER(GP_MATCHES_ACCURACY_CASE_CALLS),"calls",true);
}
/*
 * Ranks Profile
 */
void mapper_profile_print_mapper_ranks(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.RANKS\n");
  tab_fprintf(stream,"  =>  RANKS.Mapper               ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_MAPPER_ALL),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"  => RANKS.Archive.Search        ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_ARCHIVE_SEARCH_SE),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"    => RANKS.Region.Profile      ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_REGION_PROFILE),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"    => RANKS.Generate.Candidates ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_AS_GENERATE_CANDIDATES),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"    => RANKS.Process.Candidates  ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_FC_PROCESS_CANDIDATES),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"      => RANKS.Decode.Positions  ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_FC_DECODE_POSITIONS),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"      => RANKS.Compose.Regions   ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_FC_COMPOSE_REGIONS),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"    => RANKS.Verify.Candidates   ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_FC_VERIFICATION),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"      => RANKS.Decode.Positions  ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_FC_DECODE_POSITIONS),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"  => RANKS.Select.Matches        ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_ARCHIVE_SELECT_SE_MATCHES),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
}
/*
 * Approximate Search Profile Summary
 */
void mapper_profile_print_approximate_search_summary(FILE* const stream) {
  // Approximate Search
  mapper_profile_print_approximate_search(stream);
  // Region Profile
  mapper_profile_print_region_profile(stream);
  // Candidate Generation
  mapper_profile_print_candidate_generation(stream);
  // Candidate Verification
  mapper_profile_print_candidate_verification(stream);
  // Candidate Scaffold
  mapper_profile_print_candidate_scaffold(stream);
  // Candidate Realign
  mapper_profile_print_candidate_realign(stream);
  mapper_profile_print_candidate_realign_local(stream);
  // Candidate Misc
  mapper_profile_print_candidate_misc(stream);
  // Neighborhood Search
  mapper_profile_print_neighborhood_search(stream);
  mapper_profile_print_mapper_ranks(stream);
  // I/O
  mapper_profile_print_io(stream);
  // Efficiency Ratios
  mapper_profile_print_mapper_efficiency_ratios(stream);
  // Checks
  mapper_profile_print_checks(stream);
}
#else /* GEM_PROFILE DISABLED */
/*
 * Region Profile
 */
void mapper_profile_print_region_profile(FILE* const stream) {}
/*
 * Candidates Generation
 */
void mapper_profile_print_candidate_generation(FILE* const stream) {}
/*
 * Candidate Verification
 */
void mapper_profile_print_candidate_verification(FILE* const stream) {}
/*
 * Candidate realign
 */
void mapper_profile_print_candidate_realign(FILE* const stream) {}
void mapper_profile_print_candidate_realign_local(FILE* const stream) {}
/*
 * Candidate Scaffold
 */
void mapper_profile_print_candidate_scaffold(FILE* const stream) {}
/*
 * Candidate misc
 */
void mapper_profile_print_candidate_misc(FILE* const stream) {}
/*
 * Neighborhood Search
 */
void mapper_profile_print_neighborhood_search(FILE* const stream) {}
/*
 * Approximate Search
 */
void mapper_profile_print_approximate_search(FILE* const stream) {}
/*
 * Approximate Search Profile Summary
 */
void mapper_profile_print_approximate_search_summary(FILE* const stream) {}
#endif
