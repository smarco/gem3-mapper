/*
 * PROJECT: GEMMapper
 * FILE: mapper_profile.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper_profile.h"
#include "neighborhood_search.h"

#ifndef GEM_NOPROFILE /* GEM_PROFILE ENABLED */

GEM_INLINE void mapper_profile_print_mem_structs(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.MEM\n");
  // index_mem()
  // In buffers (total, used)
  // Out buffers (total, used)
  // mm_stack
  // mm_slabs => mm_pool
  // bpm_buffers
  // filtering_candidates_vectors
  // interval_sets
  // text_collection
  // matches
}
/*
 * General Profile (Efficiency Ratios)
 */
GEM_INLINE void mapper_profile_print_mapper_efficiency_ratios(FILE* const stream) {
  ///* Efficiency Ratios */
  //#define GEM_STATS_SHOW_EFFICIENCY_RATIOS(FILE,GSC_NUM_READS,SC_TIME_TH,GSC_NUM_MAPS)
  //  fprintf(FILE, "Ranks/Read               %10.3f n", (float)GET_RANK_COUNTER(SC_TIME_TH)/(float)GET_COUNTER(GSC_NUM_READS));
  //  fprintf(FILE, "Time/Read                %10.3f ms n", GET_TIME(SC_GEM_PAIR_TH)/(float)GET_COUNTER(GSC_NUM_READS)*1000.0);
  //  fprintf(FILE, "Reads/Timen");
  //  fprintf(FILE, "  --> Reads/s            %10.3f n", (float)GET_COUNTER(GSC_NUM_READS)/GET_TIME(SC_TIME_TH));
  //  fprintf(FILE, "  --> MegaReads/h        %10.3f n", (float)GET_COUNTER(GSC_NUM_READS)/GET_TIME(SC_TIME_TH)*3600.0/1000000.0);
  //  fprintf(FILE, "  --> GigaReads/d        %10.3f n", (float)GET_COUNTER(GSC_NUM_READS)/GET_TIME(SC_TIME_TH)*3600.0*24.0/1000000000.0);
  //  fprintf(FILE, "Ranks/Alignment          %10.3f n", (float)GET_RANK_COUNTER(SC_TIME_TH)/(float)GET_COUNTER(GSC_NUM_MAPS));
  //  fprintf(FILE, "Time/Alignment           %10.3f us n", GET_TIME(SC_TIME_TH)/(float)GET_COUNTER(GSC_NUM_MAPS)*1000000.0)
  //  >>>> fprintf(FILE, "      --> Matches.Duplicates            %10lu (%3.2f %%) n", GET_COUNTER(GSC_ERASED_DUPLICATES), GET_COUNT_PERCENTAGE(GSC_ERASED_DUPLICATES,SC_PAIR_MATCHES))
}
/*
 * I/O
 */
GEM_INLINE void mapper_profile_print_output(FILE* const stream,const bool paired_end,const bool map_output) {
  if (map_output) {
    if (paired_end) {
      fprintf(stream,"=> TIME.Output.MAP.PE               ");
      TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_MAP_PE),PROF_GET_TIMER(GP_MAPPER_ALL));
    } else {
      fprintf(stream,"=> TIME.Output.MAP.SE               ");
      TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_MAP_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
    }
  } else {
    if (paired_end) {
      fprintf(stream,"=> TIME.Output.SAM.PE               ");
      TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_SAM_PE),PROF_GET_TIMER(GP_MAPPER_ALL));
    } else {
      fprintf(stream,"=> TIME.Output.SAM.SE               ");
      TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_SAM_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
    }
  }
}
GEM_INLINE void mapper_profile_print_io(
    FILE* const stream,const bool paired_end,const bool map_output) {
  // High-level I/O
  fprintf(stream,    "[GEM]>Profile.Mapper.IO\n");
  tab_fprintf(stream,"  => TIME.IO.HighLevel\n");
  tab_fprintf(stream,"    => TIME.Load.Index                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_LOAD_INDEX),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Parse.Input.FASTQ          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FASTA_PARSE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Buffer.Input             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      "); mapper_profile_print_output(stream,paired_end,map_output);
  // Low-level I/O
  tab_fprintf(stream,"  => TIME.IO.LowLevel\n");
  tab_fprintf(stream,"    => TIME.BufferedInput.Reload       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.InputFile.ReloadBuffer   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FILL_BUFFER),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        --> Buffer.ReloadBuffer.Read   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BUFFERED_INPUT_BUFFER_SIZE),NULL,"B",true);
  tab_fprintf(stream,"    => TIME.BufferedOutput.Dump        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_OUTPUT_DUMP),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.OutputFile.WriteBuffer   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_WRITE_BUFFER),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        --> Bytes.written              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BYTES_WRITTEN),NULL,"B",true);
  tab_fprintf(stream,"        --> Buffer.requests            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS),NULL,"",false);
  tab_fprintf(stream,"          --> Extensions               ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_EXTENSIONS),PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS),"",false);
  tab_fprintf(stream,"          --> Stalls                   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS),NULL,"",false);
  tab_fprintf(stream,"            --> Stalls.allBusy         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS_BUSY),
      PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS),"",false);
  tab_fprintf(stream,"            --> Stalls.freePreempted   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS_NOT_PRIORITY),
      PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS),"",false);
}
/*
 * Checks
 */
GEM_INLINE void mapper_profile_print_checks(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Checks\n");
  tab_fprintf(stream,"  --> Num.Reads.Checked              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_NUM_READS),PROF_GET_COUNTER(GP_CHECK_NUM_READS),"reads",true);
  tab_fprintf(stream,"  --> Num.Maps.Checked               ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_NUM_MAPS),PROF_GET_COUNTER(GP_CHECK_NUM_MAPS),"maps ",true);
  tab_fprintf(stream,"    --> Num.Maps.Incorrect           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_INCORRECT),PROF_GET_COUNTER(GP_CHECK_NUM_MAPS),"maps ",true);
  tab_fprintf(stream,"    --> Num.Maps.Suboptimal          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_SUBOPTIMAL),PROF_GET_COUNTER(GP_CHECK_NUM_MAPS),"maps ",true);
  tab_fprintf(stream,"      --> Num.Maps.Subdominant       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_SUBOPTIMAL_SUBDOMINANT),PROF_GET_COUNTER(GP_CHECK_NUM_MAPS),"maps ",true);
  tab_fprintf(stream,"      --> Num.Maps.Distance          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_SUBOPTIMAL_DISTANCE),NULL,"     ",true);
  tab_fprintf(stream,"      --> Num.Maps.Score             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_SUBOPTIMAL_SCORE),PROF_GET_COUNTER(GP_CHECK_SUBOPTIMAL_SCORE),"     ",true);
  tab_fprintf(stream,"      --> Num.Maps.ScoreDiff         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_SUBOPTIMAL_DIFF),PROF_GET_COUNTER(GP_CHECK_SUBOPTIMAL_SCORE),"     ",true);
}
/*
 * Region Profile
 */
GEM_INLINE void mapper_profile_print_region_profile_minimal(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Region.Profile {MINIMAL}\n");
  tab_fprintf(stream,"  --> Num.Profiles             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_MINIMAL),NULL,"    ",true);
  tab_fprintf(stream,"  --> Num.Regions              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_MINIMAL_NUM_REGIONS),NULL,"    ",true);
  tab_fprintf(stream,"    --> Num.Regions.Standard   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_MINIMAL_NUM_REGIONS_STANDARD),
                       PROF_GET_COUNTER(GP_REGION_PROFILE_MINIMAL_NUM_REGIONS),"    ",true);
  tab_fprintf(stream,"    --> Num.Regions.Unique     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_MINIMAL_NUM_REGIONS_UNIQUE),
                       PROF_GET_COUNTER(GP_REGION_PROFILE_MINIMAL_NUM_REGIONS),"    ",true);
  tab_fprintf(stream,"  --> Region.length            ");
  SAMPLER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_MINIMAL_REGION_LENGTH),NULL,"nt  ");
  tab_fprintf(stream,"  --> Region.candidates        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_MINIMAL_REGION_CANDIDATES),NULL,"cand",true);
  tab_fprintf(stream,"  --> Read.candidates          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_MINIMAL_TOTAL_CANDIDATES),NULL,"cand",true);
}
GEM_INLINE void mapper_profile_print_region_profile_boost(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Region.Profile {BOOST}\n");
  tab_fprintf(stream,"  --> Num.Profiles             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_BOOST),NULL,"    ",true);
  tab_fprintf(stream,"  --> Num.Regions              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_BOOST_NUM_REGIONS),NULL,"    ",true);
  tab_fprintf(stream,"    --> Num.Regions.Standard   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_BOOST_NUM_REGIONS_STANDARD),
                       PROF_GET_COUNTER(GP_REGION_PROFILE_BOOST_NUM_REGIONS),"    ",true);
  tab_fprintf(stream,"    --> Num.Regions.Unique     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_BOOST_NUM_REGIONS_UNIQUE),
                       PROF_GET_COUNTER(GP_REGION_PROFILE_BOOST_NUM_REGIONS),"    ",true);
  tab_fprintf(stream,"  --> Region.length            ");
  SAMPLER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_BOOST_REGION_LENGTH),NULL,"nt  ");
  tab_fprintf(stream,"  --> Region.candidates        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_BOOST_REGION_CANDIDATES),NULL,"cand",true);
  tab_fprintf(stream,"  --> Read.candidates          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_BOOST_TOTAL_CANDIDATES),NULL,"cand",true);
}
GEM_INLINE void mapper_profile_print_region_profile_delimit(FILE* const stream) {
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
 * Approximate Search
 */
GEM_INLINE void mapper_profile_print_approximate_search(FILE* const stream) {
  fprintf(stream,    "[GEM]>Profile.Approximate.Search.Stages\n");
  tab_fprintf(stream,"  => TIME.Archive.Search.SE.Filtering.Stages  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Filtering.Exact                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTERING_EXACT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Filtering.ExactBoost              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTERING_EXACT_BOOST),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Filtering.Inexact                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTERING_INEXACT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Filtering.Local                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTERING_LOCAL_ALIGNMENTS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Read.Recovery                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_READ_RECOVERY),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> Archive.Search.SE.Filtering.Stages\n");
  tab_fprintf(stream,"    --> Filtering.Exact.Mapped                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_FILTERING_EXACT_MAPPED),PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"reads",true);
  tab_fprintf(stream,"    --> Filtering.ExactBoost.Mapped           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_FILTERING_EXACT_BOOST_MAPPED),PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"reads",true);
  tab_fprintf(stream,"    --> Filtering.Inexact.Mapped              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_FILTERING_INEXACT_MAPPED),PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"reads",true);
  tab_fprintf(stream,"    --> Filtering.Local.Mapped                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_FILTERING_LOCAL_ALIGNMENTS_MAPPED),PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"reads",true);
  tab_fprintf(stream,"    --> MCS.Achieved                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_ADAPTIVE_MCS),PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"mcs  ",true);
  tab_fprintf(stream,"      --> MCS.Filtering.Exact                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTERING_EXACT_MCS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      --> MCS.Filtering.ExactBoost            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTERING_EXACT_BOOST_MCS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      --> MCS.Filtering.Inexact               ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTERING_INEXACT_MCS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      --> MCS.Filtering.Local                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTERING_LOCAL_ALIGNMENTS_MCS),PROF_GET_TIMER(GP_MAPPER_ALL));
}
GEM_INLINE void mapper_profile_print_approximate_search_ranks(FILE* const stream) { /* TODO */ }
/*
 * Archive Search
 */
GEM_INLINE void mapper_profile_print_archive_search_se(FILE* const stream) {
  fprintf(stream,    "[GEM]>Profile.ArchiveSearch.SE\n");
  tab_fprintf(stream,"  => TIME.Archive.Search.SE             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Approximate.Search          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_MAIN),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Neighborhood.Search       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_NEIGHBORHOOD_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Region.Profile            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_REGION_PROFILE_ADAPTIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Generate.Candidates       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_GENERATE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Process.Candidates        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_PROCESS_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Verifying                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_VERIFICATION),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Realign                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Select.Matches              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SELECT_SE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Score.Matches             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SCORE_SE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
}
GEM_INLINE void mapper_profile_print_archive_search_pe(FILE* const stream) {
  fprintf(stream,    "[GEM]>Profile.ArchiveSearch.PE\n");
  tab_fprintf(stream,"  => TIME.Archive.Search.PE                         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    (A) Map.Both.Ends\n");
  tab_fprintf(stream,"        => TIME.Archive.Init                        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_INIT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Archive.Generate.Candidates         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_GENERATE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Archive.Finish.Search               ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Archive.Find.Paired.Matches         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_PAIRED_MATCHES_FIND_PAIRS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Archive.Select.PairedMatches        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SELECT_PE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Score.PairedMatches                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SCORE_PE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    (B) Paired.Filtering\n");
  tab_fprintf(stream,"        => TIME.Archive.Discard.Filtering.Regions   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_DISCARD_FILTERING_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Archive.Verify.Filtering.Regions    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    (C) Map.Extension\n");
  tab_fprintf(stream,"        => TIME.Archive.Extension                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.Filtering.Pair.By.Extension       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_EXTENSION),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.Filtering.Recovery.By.Extension   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_RECOVER_BY_EXTENSION),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.Archive.Extend.Candidates\n");
  tab_fprintf(stream,"            => TIME.Filtering.Extend.Match          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_MATCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"              => TIME.Filtering.Retrieve.Candidates ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_RETRIEVE_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"              => TIME.Filtering.Verify.Candidates   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"              => TIME.Filtering.Realign.Candidates  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    |> Archive.Discard.Filtering.Regions\n");
  tab_fprintf(stream,"      --> Paired.Filtered                         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_PAIRED_FILTERED),
                       PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"       ",true);
  tab_fprintf(stream,"        --> Paired.Filtering.Success              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_PAIRED_FILTERING_SUCCESS),
                       PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"       ",true);
  tab_fprintf(stream,"        --> Discarded.Filtering.Regions           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_DISCARD_FILTERING_REGIONS_NOT_CONCORDANT),
                       PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_DISCARD_FILTERING_REGIONS_TOTAL),"regions",true);
  tab_fprintf(stream,"    |> Archive.Extend.Candidates (Including Recovery)\n");
  tab_fprintf(stream,"      --> Candidate.Regions.Length                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_EXTEND_VERIFY_CANDIDATE_LENGTH),NULL,"nt",true);
  tab_fprintf(stream,"      --> Num.Extended.End1                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END1),
                       PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"ext",true);
  tab_fprintf(stream,"        --> Extended.End1.Success                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END1_SUCCESS),
                       PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"ext",true);
  tab_fprintf(stream,"      --> Num.Extended.End2                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END2),
                       PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"ext",true);
  tab_fprintf(stream,"        --> Extended.End2.Success                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END2_SUCCESS),
                       PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"ext",true);
  tab_fprintf(stream,"      --> Num.Recovered.By.Extension              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_RECOVER_BY_EXTENSION_HIT),
                       PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"ext",true);
  tab_fprintf(stream,"        --> Recovered.By.Extension.End1           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_RECOVER_BY_EXTENSION_END1),
                       PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"ext",true);
  tab_fprintf(stream,"        --> Recovered.By.Extension.End2           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_RECOVER_BY_EXTENSION_END2),
                       PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"ext",true);
}
GEM_INLINE void mapper_profile_print_archive_search_group(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Archive.Search.Group\n");
  // Search groups (BPM-Buffers)
  tab_fprintf(stream,"  |> Archive.Search.Group\n");
  tab_fprintf(stream,"    --> Candidates.processed          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BPM_GPU_BUFFER_NUM_CANDIDATES),NULL,"candidates",true);
  tab_fprintf(stream,"    --> Candidates.length             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BPM_GPU_BUFFER_CANDIDATES_LENGTH),NULL,"nt        ",true);
  tab_fprintf(stream,"    --> Num.Circular.BPMBuffers.Used  ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_GROUP_BUFFERS_USED),NULL,"buffers   ",true);
  tab_fprintf(stream,"      --> Queries                     ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_BPM_GPU_BUFFER_USAGE_QUERIES));
  tab_fprintf(stream,"      --> PeqEntries                  ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_BPM_GPU_BUFFER_USAGE_PEQ_ENTRIES));
  tab_fprintf(stream,"      --> Candidates                  ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_BPM_GPU_BUFFER_USAGE_CANDIDATES));
}
GEM_INLINE void mapper_profile_print_archive_search_cuda(FILE* const stream,const bool map_output) {
  tab_fprintf(stream,"[GEM]>Profile.CUDA.ArchiveSearch\n");
  tab_fprintf(stream,"  => TIME.CUDA.Thread                        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => 0.TIME.CUDA.Buffers.Init              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_INIT),PROF_GET_TIMER(GP_MAPPER_ALL));
  // Generating
  tab_fprintf(stream,"    => 1.TIME.CUDA.Generating                ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD_GENERATING),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Buffer.Input                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Parse.Input.FASTQ              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FASTA_PARSE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Prepare.Sequence       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PREPARE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Generate.Candidates    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_GENERATE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Copy.Candidates        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_COPY_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  // Verifying
  tab_fprintf(stream,"    => 2.TIME.CUDA.Verifying\n");
  tab_fprintf(stream,"      => {GPU}\n");
  tab_fprintf(stream,"        => TIME.CUDA.Send.Delay              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.CUDA.Duty.Cycle              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_CHECK_TIME),NULL);
  tab_fprintf(stream,"        => TIME.CUDA.Retrieve.Delay          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_GROUP_RETRIEVE_CANDIDATES_DELAY),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => {CPU}\n");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD_SELECTING),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Archive.Retrieve.Candidates  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.Retrieve.Candidates        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_RETRIEVE_BPM_BUFFER_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.Realign                    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_BPM_BUFFER_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Archive.Finish.Search        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Archive.Select.Matches       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SELECT_SE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Restart.Unfit.Searches       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD_RESTART_UNFIT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        "); mapper_profile_print_output(stream,false,map_output);
  // Archive Search Groups
  mapper_profile_print_archive_search_group(stream);
}
GEM_INLINE void mapper_profile_print_archive_search_pe_cuda(FILE* const stream,const bool map_output) {
  tab_fprintf(stream,"[GEM]>Profile.CUDA.ArchiveSearch\n");
  tab_fprintf(stream,"  => TIME.CUDA.Thread                        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => 0.TIME.CUDA.Buffers.Init              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_INIT),PROF_GET_TIMER(GP_MAPPER_ALL));
  // Generating
  tab_fprintf(stream,"    => 1.TIME.CUDA.Generating                ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD_GENERATING),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Buffer.Input                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Parse.Input.FASTQ              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FASTA_PARSE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Generate.Candidates    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_GENERATE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Copy.Candidates        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_COPY_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  // Verifying
  tab_fprintf(stream,"    => 2.TIME.CUDA.Verifying                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD_SELECTING),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => {GPU}\n");
  tab_fprintf(stream,"        => TIME.CUDA.Send.Delay              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.CUDA.Duty.Cycle              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_CHECK_TIME),NULL);
  tab_fprintf(stream,"        => TIME.CUDA.Retrieve.Delay          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_GROUP_RETRIEVE_CANDIDATES_DELAY),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => {CPU}\n");
  tab_fprintf(stream,"        => TIME.Archive.Retrieve.Candidates  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.Retrieve.Candidates        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_RETRIEVE_BPM_BUFFER_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.Realign                    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_BPM_BUFFER_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Archive.Finish.Search        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        "); mapper_profile_print_output(stream,false,map_output);
  // Archive Search Groups
  mapper_profile_print_archive_search_group(stream);
}
/*
 * Filtering Verification
 */
GEM_INLINE void mapper_profile_print_verifying(FILE* const stream,const bool verification_profile) {
  tab_fprintf(stream,"[GEM]>Profile.Filtering\n");
  // Verifying
  tab_fprintf(stream,"  => TIME.Process.Candidates                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_PROCESS_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Decode.Positions                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_DECODE_POSITIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Compose.Regions                    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_COMPOSE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  => TIME.Verifying                            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_VERIFICATION),PROF_GET_TIMER(GP_MAPPER_ALL));
  if (verification_profile) {
  tab_fprintf(stream,"    => TIME.Retrieve.Candidate.Regions         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_RETRIEVE_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Verify.Candidate.Regions           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_VERIFY_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  //tab_fprintf(stream,"  => TIME.Kmer.Counting                      ");
  //TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_KMER_COUNTER_FILTER),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.BPM.Align                        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_TILED),PROF_GET_TIMER(GP_MAPPER_ALL));
  }
  tab_fprintf(stream,"  |> Filtering.Candidates\n");
  tab_fprintf(stream,"    --> Candidate.Positions                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"positions",true);
  tab_fprintf(stream,"      --> Duplicates.Removed                      ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS_DUPLICATED),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"positions",true);
  tab_fprintf(stream,"      --> Decode.LF.Dist                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FMIDX_LOOKUP_DIST),NULL,"lf       ",true);
  tab_fprintf(stream,"    --> Candidate.Regions                         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Duplicates.Removed                      ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGIONS_DUPLICATED),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Length                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),NULL,"nt       ",true);
  tab_fprintf(stream,"      --> Candidate.Tiles                         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BMP_TILED_NUM_TILES),PROF_GET_COUNTER(GP_BMP_TILED_NUM_TILES),"tiles    ",true);
  tab_fprintf(stream,"        --> Candidate.Tiles.Verified              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BMP_TILED_NUM_TILES_VERIFIED),PROF_GET_COUNTER(GP_BMP_TILED_NUM_TILES),"tiles    ",true);
  if (verification_profile) {
  tab_fprintf(stream,"        --> BPM.Accepted                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_LEVENSHTEIN_ACCEPTED),PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),"         ",true);
  tab_fprintf(stream,"          --> BPM.Quick-Abandon                   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BPM_QUICK_ABANDON),PROF_GET_COUNTER(GP_BMP_TILED_NUM_TILES_VERIFIED),"         ",true);
  }
  //  tab_fprintf(stream,"      --> FC.Kmer.Counting.Discarded                ");
  //  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_KMER_COUNTER_FILTER_DISCARDED),PROF_GET_COUNTER(GP_FC_NUM_CANDIDATE_REGIONS),"",true);
}
GEM_INLINE void mapper_profile_print_realign(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Realign\n");
  tab_fprintf(stream,"  => TIME.Realign                                  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Realign.Exact                          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHES_ALIGN_EXACT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Chain.Matching.Regions                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHING_REGIONS_CHAIN),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Extend.Matching.Regions                ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHING_REGIONS_EXTEND),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Scaffold.Matching.Regions              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHING_REGIONS_SCAFFOLD),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Realign.SWG                            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHES_ALIGN_SWG),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.SWG.Core.Banded                      ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_SWG_ALIGN_BANDED),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> Realign.Regions\n");
  tab_fprintf(stream,"    --> Candidate.Positions                        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"positions",true);
  tab_fprintf(stream,"    --> Candidate.Regions                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"    --> Accepted.Regions                           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Accepted.Exact         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_EXACT),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Accepted.Inexact       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_INEXACT),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"        --> Matching.Regions.Chain.Success         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCHING_REGIONS_CHAIN_SUCCESS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"          --> Matching.Regions.Coverage            ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_MATCHING_REGIONS_CHAIN_COVERAGE));
  tab_fprintf(stream,"          --> Matching.Regions.Extended.Coverage   ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_MATCHING_REGIONS_EXTEND_COVERAGE));
  tab_fprintf(stream,"        --> Matching.Regions.Scaffolded            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCHING_REGIONS_SCAFFOLDED),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"          --> Matching.Regions.Scaffolded.Coverage ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_MATCHING_REGIONS_SCAFFOLD_COVERAGE));
  tab_fprintf(stream,"  |> Realign.Length\n");
  tab_fprintf(stream,"    --> Candidate.Regions.Length                   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),"nt       ",true);
  tab_fprintf(stream,"    --> Accepted.Regions.Length                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS_LENGTH),PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),"nt       ",true);
  tab_fprintf(stream,"    --> SWG.Banded.Aligned.Length                  ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SWG_ALIGN_BANDED_LENGTH),PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),"nt       ",true);
}
GEM_INLINE void mapper_profile_print_filtering(FILE* const stream,const bool verification_profile) {
  // Verifying
  mapper_profile_print_verifying(stream,verification_profile);
  // Realign
  mapper_profile_print_realign(stream);
}
GEM_INLINE void mapper_profile_print_filtering_verifying_ranks(FILE* const stream) { /*TODO*/ }
/*
 * Neighborhood Search
 */
GEM_INLINE void mapper_profile_print_neighborhood_search(FILE* const stream) {
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
  tab_fprintf(stream,"      --> Nodes.Closed.Depth                   ");
  SAMPLER_PRINT(stream,&_ns_nodes_closed_depth,PROF_GET_COUNTER(GP_NS_NODES_EXPLORED),"bases");
}
GEM_INLINE void mapper_profile_print_neighborhood_search_ranks(FILE* const stream) { /*TODO*/ }
/*
 * Generate Candidates
 */
GEM_INLINE void mapper_profile_print_generate_candidates(FILE* const stream) {
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
 * Ranks Profile
 */
GEM_INLINE void mapper_profile_print_mapper_ranks(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Mapper\n");
  tab_fprintf(stream,"  => RANKS.Mapper                  ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_MAPPER_ALL),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"    => RANKS.Archive.Search        ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_ARCHIVE_SEARCH_SE),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"      => RANKS.Region.Profile      ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_REGION_PROFILE_ADAPTIVE),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"      => RANKS.Generate.Candidates ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_AS_GENERATE_CANDIDATES),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"      => RANKS.Process.Candidates  ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_FC_PROCESS_CANDIDATES),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"        => RANKS.Decode.Positions  ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_FC_DECODE_POSITIONS),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"        => RANKS.Compose.Regions   ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_FC_COMPOSE_REGIONS),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"      => RANKS.Verify.Candidates   ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_FC_VERIFICATION),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"        => RANKS.Decode.Positions  ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_FC_DECODE_POSITIONS),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"    => RANKS.Select.Matches        ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_ARCHIVE_SELECT_SE_MATCHES),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
}
/*
 * Mapper
 */
GEM_INLINE void mapper_profile_print_mapper_commons(
    FILE* const stream,const bool paired_end,const bool map_output,const uint64_t num_threads) {
  // Archive Search SE
  if (paired_end) {
    mapper_profile_print_archive_search_pe(stream);
    mapper_profile_print_archive_search_se(stream);
  } else {
    mapper_profile_print_archive_search_se(stream);
  }
  // Approximate Search
  mapper_profile_print_approximate_search(stream);
  // Generate Candidates
  mapper_profile_print_generate_candidates(stream);
  // Region Profile
  mapper_profile_print_region_profile_minimal(stream);
  mapper_profile_print_region_profile_boost(stream);
  mapper_profile_print_region_profile_delimit(stream);
  // Filtering Verification
  mapper_profile_print_filtering(stream,true);
  if (num_threads==1) {
    // Neighborhood Search
    mapper_profile_print_neighborhood_search(stream);
    // Ranks
    mapper_profile_print_mapper_ranks(stream);
  }
  // I/O
  mapper_profile_print_io(stream,paired_end,map_output);
  // Checks
  mapper_profile_print_checks(stream);
}
GEM_INLINE void mapper_profile_print_mapper_single_end(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {
  // Main
  tab_fprintf(stream,"[GEM]>Profile.Mapper\n");
  tab_fprintf(stream,"  => TIME.Mapper                           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_ALL),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Load.Index                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_LOAD_INDEX),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Parse.Input.FASTQ              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FASTA_PARSE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Archive.Search.SE              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    "); mapper_profile_print_output(stream,false,map_output);
  // Commons
  mapper_profile_print_mapper_commons(stream,false,map_output,num_threads);
}
GEM_INLINE void mapper_profile_print_mapper_single_end_cuda(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {
  // Main
  tab_fprintf(stream,"[GEM]>Profile.Mapper\n");
  tab_fprintf(stream,"  => TIME.Mapper                           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_ALL),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Load.Index                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_LOAD_INDEX),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.CUDA.Thread                    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD),PROF_GET_TIMER(GP_MAPPER_ALL));
  // CUDA
  mapper_profile_print_archive_search_cuda(stream,map_output);
  // Commons
  mapper_profile_print_mapper_commons(stream,false,map_output,num_threads);
}
/*
 * Mapper PE
 */
GEM_INLINE void mapper_profile_print_mapper_paired_end(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {
  // All
  tab_fprintf(stream,"[GEM]>Profile.Mapper\n");
  tab_fprintf(stream,"  => TIME.Mapper                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_ALL),PROF_GET_TIMER(GP_MAPPER_ALL));
  // Commons
  mapper_profile_print_mapper_commons(stream,true,map_output,num_threads);
}
GEM_INLINE void mapper_profile_print_mapper_paired_end_cuda(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {
  // All
  tab_fprintf(stream,"[GEM]>Profile.Mapper\n");
  tab_fprintf(stream,"  => TIME.Mapper                           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_ALL),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Load.Index                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_LOAD_INDEX),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.CUDA.Thread                    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD),PROF_GET_TIMER(GP_MAPPER_ALL));
  // CUDA
  mapper_profile_print_archive_search_pe_cuda(stream,map_output);
  // Commons
  mapper_profile_print_mapper_commons(stream,true,map_output,num_threads);
}
#else /* GEM_PROFILE DISABLED */
GEM_INLINE void mapper_profile_print_mapper_single_end(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {}
GEM_INLINE void mapper_profile_print_mapper_single_end_cuda(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {}
GEM_INLINE void mapper_profile_print_mapper_paired_end(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {}
GEM_INLINE void mapper_profile_print_mapper_paired_end_cuda(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {}
#endif
