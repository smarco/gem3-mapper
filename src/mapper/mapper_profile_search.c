/*
 * PROJECT: GEMMapper
 * FILE: mapper_profile_search.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper/mapper_profile_search.h"
#include "mapper/mapper_profile_misc.h"

#ifdef GEM_PROFILE /* GEM_PROFILE ENABLED */
/*
 * Region Profile
 */
void mapper_profile_print_region_profile_fixed(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Region.Profile {FIXED}\n");
  tab_fprintf(stream,"  --> Num.Profiles             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_FIXED),NULL,"    ",true);
  tab_fprintf(stream,"  --> Num.Regions              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_FIXED_NUM_REGIONS),NULL,"    ",true);
  tab_fprintf(stream,"    --> Num.Regions.Standard   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_FIXED_NUM_REGIONS_STANDARD),
                       PROF_GET_COUNTER(GP_REGION_PROFILE_FIXED_NUM_REGIONS),"    ",true);
  tab_fprintf(stream,"    --> Num.Regions.Unique     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_FIXED_NUM_REGIONS_UNIQUE),
                       PROF_GET_COUNTER(GP_REGION_PROFILE_FIXED_NUM_REGIONS),"    ",true);
  tab_fprintf(stream,"  --> Region.length            ");
  SAMPLER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_FIXED_REGION_LENGTH),NULL,"nt  ");
  tab_fprintf(stream,"  --> Region.candidates        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_FIXED_REGION_CANDIDATES),NULL,"cand",true);
  tab_fprintf(stream,"  --> Read.candidates          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_FIXED_TOTAL_CANDIDATES),NULL,"cand",true);
}
void mapper_profile_print_region_profile_lightweight(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Region.Profile {LIGHTWEIGHT}\n");
  tab_fprintf(stream,"  --> Num.Profiles             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT),NULL,"    ",true);
  tab_fprintf(stream,"  --> Num.Regions              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT_NUM_REGIONS),NULL,"    ",true);
  tab_fprintf(stream,"    --> Num.Regions.Standard   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT_NUM_REGIONS_STANDARD),
                       PROF_GET_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT_NUM_REGIONS),"    ",true);
  tab_fprintf(stream,"    --> Num.Regions.Unique     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT_NUM_REGIONS_UNIQUE),
                       PROF_GET_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT_NUM_REGIONS),"    ",true);
  tab_fprintf(stream,"  --> Region.length            ");
  SAMPLER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT_REGION_LENGTH),NULL,"nt  ");
  tab_fprintf(stream,"  --> Region.candidates        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT_REGION_CANDIDATES),NULL,"cand",true);
  tab_fprintf(stream,"  --> Read.candidates          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_LIGHTWEIGHT_TOTAL_CANDIDATES),NULL,"cand",true);
}
void mapper_profile_print_region_profile_heavyweight(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Region.Profile {HEAVYWEIGHT}\n");
  tab_fprintf(stream,"  --> Num.Profiles             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT),NULL,"    ",true);
  tab_fprintf(stream,"  --> Num.Regions              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT_NUM_REGIONS),NULL,"    ",true);
  tab_fprintf(stream,"    --> Num.Regions.Standard   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT_NUM_REGIONS_STANDARD),
                       PROF_GET_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT_NUM_REGIONS),"    ",true);
  tab_fprintf(stream,"    --> Num.Regions.Unique     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT_NUM_REGIONS_UNIQUE),
                       PROF_GET_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT_NUM_REGIONS),"    ",true);
  tab_fprintf(stream,"  --> Region.length            ");
  SAMPLER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT_REGION_LENGTH),NULL,"nt  ");
  tab_fprintf(stream,"  --> Region.candidates        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT_REGION_CANDIDATES),NULL,"cand",true);
  tab_fprintf(stream,"  --> Read.candidates          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_HEAVYWEIGHT_TOTAL_CANDIDATES),NULL,"cand",true);
}
void mapper_profile_print_region_profile_delimit(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Region.Profile {DELIMIT}\n");
  tab_fprintf(stream,"  --> Num.Profiles             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_DELIMIT),NULL,"    ",true);
  tab_fprintf(stream,"  --> Num.Regions              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_DELIMIT_NUM_REGIONS),NULL,"    ",true);
  tab_fprintf(stream,"    --> Num.Regions.Standard   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_DELIMIT_NUM_REGIONS_STANDARD),
                       PROF_GET_COUNTER(GP_REGION_PROFILE_DELIMIT_NUM_REGIONS),"    ",true);
  tab_fprintf(stream,"    --> Num.Regions.Unique     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_DELIMIT_NUM_REGIONS_UNIQUE),
                       PROF_GET_COUNTER(GP_REGION_PROFILE_DELIMIT_NUM_REGIONS),"    ",true);
  tab_fprintf(stream,"  --> Region.length            ");
  SAMPLER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_DELIMIT_REGION_LENGTH),NULL,"nt  ");
  tab_fprintf(stream,"  --> Region.candidates        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_DELIMIT_REGION_CANDIDATES),NULL,"cand",true);
  tab_fprintf(stream,"  --> Read.candidates          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_DELIMIT_TOTAL_CANDIDATES),NULL,"cand",true);
}
/*
 * Candidates Generation
 */
void mapper_profile_print_candidate_generation(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Generate.Candidates\n");
  tab_fprintf(stream,"  => TIME.Generate.Candidates         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_GENERATE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Generate.Candidates.D2    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_GENERATE_CANDIDATES_SEARCH_D2),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Generate.Candidates.D1    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_GENERATE_CANDIDATES_SEARCH_D1),PROF_GET_TIMER(GP_MAPPER_ALL));
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
  tab_fprintf(stream,"      --> Generate.D2                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D2),
                       PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS),"regions",true);
  tab_fprintf(stream,"        --> Generate.D2.Hit           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D2_HIT),
                       PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS),"regions",true);
  tab_fprintf(stream,"      --> Generate.D1                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D1),
                       PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS),"regions",true);
  tab_fprintf(stream,"        --> Generate.D1.Hit           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D1_HIT),
                       PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS),"regions",true);
  tab_fprintf(stream,"      --> Generate.D0                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D0_HIT),
                       PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS),"regions",true);
  tab_fprintf(stream,"        --> Generate.D0.Hit           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D0_HIT),
                       PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_NUM_ELEGIBLE_REGIONS),"regions",true);
  tab_fprintf(stream,"    --> Candidates.Generated          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),
                       PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"pos",true);
  tab_fprintf(stream,"      --> Candidates.Generated.D2     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D2_HIT_CANDIDATES),
                       PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"pos",true);
  tab_fprintf(stream,"      --> Candidates.Generated.D1     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D1_HIT_CANDIDATES),
                       PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"pos",true);
  tab_fprintf(stream,"      --> Candidates.Generated.D0     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_GENERATE_CANDIDATES_SEARCH_D0_HIT_CANDIDATES),
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
  tab_fprintf(stream,"      --> Candidate.Regions.Matching.Regions      ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGION_MATCHING_REGIONS_TOTAL),NULL,"regions  ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Matching.Coverage     ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGION_MATCHING_COVERAGE),"");
  tab_fprintf(stream,"      --> Kmer.Counting\n");
  tab_fprintf(stream,"        --> Kmer.Counting.Discarded               ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_KMER_COUNTER_FILTER_DISCARDED),PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),"         ",true);
  tab_fprintf(stream,"        --> Kmer.Counting.Accepted                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_KMER_COUNTER_FILTER_ACCEPTED),PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),"         ",true);
  tab_fprintf(stream,"          --> Kmer.Counting.NA                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_KMER_COUNTER_FILTER_NA),PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),"         ",true);
  tab_fprintf(stream,"      --> BPM\n");
  tab_fprintf(stream,"        --> BPM.Tiles                             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BMP_DISTANCE_NUM_TILES),PROF_GET_COUNTER(GP_BMP_DISTANCE_NUM_TILES),"tiles    ",true);
  tab_fprintf(stream,"          --> BPM.Tiles.Computed                  ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BMP_DISTANCE_NUM_TILES_VERIFIED),PROF_GET_COUNTER(GP_BMP_DISTANCE_NUM_TILES),"tiles    ",true);
  tab_fprintf(stream,"          --> BPM.Tile.Tall                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BPM_DISTANCE_KEY_LENGTH),NULL,"nt       ",true);
  tab_fprintf(stream,"          --> BPM.Tile.Length                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BPM_DISTANCE_TEXT_LENGTH),NULL,"nt       ",true);
  tab_fprintf(stream,"        --> BPM.Cells.Computed                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BPM_DISTANCE_CELLS),NULL,"cells    ",true);
  tab_fprintf(stream,"        --> BPM.Discarded                         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_DISCARDED_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),"         ",true);
  tab_fprintf(stream,"          --> BPM.Quick.Abandon                   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BPM_DISTANCE_QUICK_ABANDON),PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),"         ",true);
  tab_fprintf(stream,"        --> BPM.Accepted                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),"         ",true);
  tab_fprintf(stream,"    --> Accepted.Regions                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
}
/*
 * Candidate realign
 */
void mapper_profile_print_candidate_realign(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Candidate.Realign\n");
  tab_fprintf(stream,"  => TIME.Realign                                  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Realign.Cache\n");
  tab_fprintf(stream,"      => TIME.Realign.Cache.Compute.Footpint       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_CACHE_COMPUTE_FOOTPRINT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Realign.Cache.Search                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_CACHE_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Scaffold.Alignment                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCH_SCAFFOLD_ALIGNMENT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Scaffold.Region.Chain                ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCH_SCAFFOLD_CHAIN_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Scaffold.Edit                        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCH_SCAFFOLD_EDIT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Realign.Exact                          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHES_ALIGN_EXACT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Realign.SWG                            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHES_ALIGN_SWG),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Realign.SWG.Core.Banded              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_SWG_ALIGN_BANDED),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  => TIME.Realign.Local                            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_LOCAL_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> Realign.Regions\n");
  tab_fprintf(stream,"    --> Candidate.Positions                        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"positions",true);
  tab_fprintf(stream,"    --> Candidate.Regions                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"    --> Accepted.Regions                           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Scaffold\n");
  tab_fprintf(stream,"        --> Scaffold.Region.Chain                  ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_CHAIN_REGIONS_SCAFFOLDS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"          --> Scaffold.Region.Matching.Regions     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_ALIGNMENT_MATCHING_REGIONS),NULL,"regions  ",true);
  tab_fprintf(stream,"          --> Scaffold.Region.Matching.Coverage    ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_ALIGNMENT_MATCHING_COVERAGE),"");
  tab_fprintf(stream,"          --> Scaffold.Region.Chain.Coverage       ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_CHAIN_REGIONS_COVERAGE),"");
  tab_fprintf(stream,"        --> Scaffold.Edit                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_SCAFFOLDS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"          --> Scaffold.Edit.Coverage               ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_COVERAGE),"");
  tab_fprintf(stream,"          --> Scaffold.Edit.Tiles\n");
  tab_fprintf(stream,"            --> Scaffold.Edit.Tiles.Total          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_TOTAL),
                       PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_TOTAL),"tiles    ",true);
  tab_fprintf(stream,"            --> Scaffold.Edit.Tiles.Skipped        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_SKIPPED),
                       PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_TOTAL),"tiles    ",true);
  tab_fprintf(stream,"            --> Scaffold.Edit.Tiles.Align          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_ALIGN),
                       PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_TILES_TOTAL),"tiles    ",true);
  tab_fprintf(stream,"          --> Scaffold.Edit.Cells.Computed         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCH_SCAFFOLD_EDIT_CELLS),NULL,"cells    ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Aligned.Pruned         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_SELECT_PRUNE_HIT),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Aligned.Cached         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_CACHE_SEARCH_HIT),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"    --> Aligned.Regions                            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ALIGNED_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Aligned.Exact          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ALIGNED_EXACT),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Aligned.Inexact        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ALIGNED_INEXACT),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"  |> Realign.Length\n");
  tab_fprintf(stream,"    --> Candidate.Regions.Length                   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),"nt       ",true);
  tab_fprintf(stream,"    --> Aligned.Regions.Length                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ALIGNED_REGIONS_LENGTH),PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),"nt       ",true);
  tab_fprintf(stream,"    --> SWG.Banded.Aligned.Length                  ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SWG_ALIGN_BANDED_LENGTH),PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),"nt       ",true);
  tab_fprintf(stream,"    --> SWG.Banded.Aligned.Cells.Computed          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SWG_ALIGN_BANDED_CELLS),NULL,"cells    ",true);
  tab_fprintf(stream,"  |> Realign.Local\n");
  tab_fprintf(stream,"    --> Candidate.Regions                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGION_LOCAL),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"    --> Candidate.Regions.Aligned.Local            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGION_LOCAL_ALIGNED),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
}
/*
 * Neighborhood Search
 */
void mapper_profile_print_neighborhood_search(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.NeighborhoodSearch\n");
  tab_fprintf(stream,"  => TIME.NS.BestMatch                         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_NS_BEST_MATCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> NS.Nodes\n");
  tab_fprintf(stream,"    --> Nodes.Explored                         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_NS_NODES_EXPLORED),PROF_GET_COUNTER(GP_NS_NODES_EXPLORED),"nodes",true);
  tab_fprintf(stream,"      --> Nodes.Explored.mTable                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_NS_NODES_EXPLORED_MTABLE),PROF_GET_COUNTER(GP_NS_NODES_EXPLORED),"nodes",true);
  tab_fprintf(stream,"    --> Nodes.Match                            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_NS_NODE_SUCCESS),PROF_GET_COUNTER(GP_NS_NODES_EXPLORED),"nodes",true);
  tab_fprintf(stream,"    --> Nodes.Closed                           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_NS_NODE_CLOSED),PROF_GET_COUNTER(GP_NS_NODES_EXPLORED),"nodes",true);
  tab_fprintf(stream,"      --> Nodes.Failed.Optimization            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_NS_FAILED_OPT),PROF_GET_COUNTER(GP_NS_NODES_EXPLORED),"nodes",true);
//  tab_fprintf(stream,"      --> Nodes.Closed.Depth                   ");
//  SAMPLER_PRINT(stream,&_ns_nodes_closed_depth,PROF_GET_COUNTER(GP_NS_NODES_EXPLORED),"bases");
}
void mapper_profile_print_neighborhood_search_ranks(FILE* const stream) {
  /*TODO*/
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
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_NSEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Region.Profile                  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_REGION_PROFILE_ADAPTIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Generate.Candidates             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_GENERATE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Process.Candidates              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_PROCESS_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Verifying                       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_VERIFICATION),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Realign                         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  => TIME.Approximate.Search.Stages\n");
  tab_fprintf(stream,"    => TIME.Approximate.Search                ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_MAIN),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Filtering.Exact                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTERING_EXACT),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Filtering.Inexact               ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTERING_INEXACT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Filtering.Local                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTERING_LOCAL_ALIGN),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Read.Recovery                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_READ_RECOVERY),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Neighborhood.Search             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_NEIGHBORHOOD_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> Approximate.Search.Stages\n");
  tab_fprintf(stream,"    --> Filtering.Exact.Mapped                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_FILTERING_EXACT_MAPPED),PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"reads",true);
  tab_fprintf(stream,"    --> Filtering.Exact.MCS                   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_FILTERING_EXACT_MCS),PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"mcs  ",true);
  //  tab_fprintf(stream,"    --> Filtering.Inexact.Mapped              ");
  //  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_FILTERING_INEXACT_MAPPED),PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"reads",true);
  //  tab_fprintf(stream,"      --> MCS.Filtering.Inexact               ");
  //  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_FILTERING_INEXACT_MCS),PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"mcs  ",true);
}
void mapper_profile_print_approximate_search_ranks(FILE* const stream) {
  /*TODO*/
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
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_REGION_PROFILE_ADAPTIVE),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
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
void mapper_profile_print_approximate_search_summary(
    FILE* const stream,const bool paired_end,
    const bool cuda_workflow,const bool map_output,
    const uint64_t num_threads) {
  // Approximate Search
  mapper_profile_print_approximate_search(stream);
  // Region Profile
  if (cuda_workflow) {
    mapper_profile_print_region_profile_fixed(stream);
    //mapper_profile_print_region_profile_heavyweight(stream);
    mapper_profile_print_region_profile_lightweight(stream);
  } else {
    mapper_profile_print_region_profile_lightweight(stream);
    // mapper_profile_print_region_profile_delimit(stream);
  }
  // Candidate Generation
  mapper_profile_print_candidate_generation(stream);
  // Candidate Verification
  mapper_profile_print_candidate_verification(stream);
  // Candidate Realign
  mapper_profile_print_candidate_realign(stream);
  // Neighborhood Search
  if (num_threads==1) {
    mapper_profile_print_neighborhood_search(stream);
    mapper_profile_print_mapper_ranks(stream);
  }
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
void mapper_profile_print_region_profile_fixed(FILE* const stream) {}
void mapper_profile_print_region_profile_lightweight(FILE* const stream) {}
void mapper_profile_print_region_profile_heavyweight(FILE* const stream) {}
void mapper_profile_print_region_profile_delimit(FILE* const stream) {}
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
/*
 * Neighborhood Search
 */
void mapper_profile_print_neighborhood_search(FILE* const stream) {}
void mapper_profile_print_neighborhood_search_ranks(FILE* const stream) {}
/*
 * Approximate Search
 */
void mapper_profile_print_approximate_search(FILE* const stream) {}
void mapper_profile_print_approximate_search_ranks(FILE* const stream) {}
/*
 * Approximate Search Profile Summary
 */
void mapper_profile_print_approximate_search_summary(
    FILE* const stream,const bool paired_end,
    const bool cuda_workflow,const bool map_output,
    const uint64_t num_threads) {}
#endif
