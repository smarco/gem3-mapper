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
GEM_INLINE void mapper_profile_print_output(
    FILE* const stream,const bool paired_end,const bool map_output) {
  if (map_output) {
    if (paired_end) {
      tab_fprintf(stream,"  => TIME.Output.MAP.PE                  ");
      TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_MAP_PE),PROF_GET_TIMER(GP_MAPPER_ALL));
    } else {
      tab_fprintf(stream,"  => TIME.Output.MAP.SE                  ");
      TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_MAP_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
    }
  } else {
    if (paired_end) {
      tab_fprintf(stream,"  => TIME.Output.SAM.PE                  ");
      TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_SAM_PE),PROF_GET_TIMER(GP_MAPPER_ALL));
    } else {
      tab_fprintf(stream,"  => TIME.Output.SAM.SE                  ");
      TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_SAM_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
    }
  }
}
GEM_INLINE void mapper_profile_print_io(
    FILE* const stream,const bool paired_end,const bool map_output) {
  // Low-level I/O
  fprintf(stream,    "[GEM]>Profile.Mapper.io.LowLevel\n");
  tab_fprintf(stream,"  => TIME.BufferedInput.Reload       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.InputFile.ReloadBuffer   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FILL_BUFFER),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      --> Buffer.ReloadBuffer.Read   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BUFFERED_INPUT_BUFFER_SIZE),NULL,"B",true);
  tab_fprintf(stream,"  => TIME.BufferedOutput.Dump        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_OUTPUT_DUMP),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.OutputFile.WriteBuffer   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_WRITE_BUFFER),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      --> Bytes.written              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BYTES_WRITTEN),NULL,"B",true);
  tab_fprintf(stream,"      --> Buffer.requests            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS),NULL,"",false);
  tab_fprintf(stream,"        --> Extensions               ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_EXTENSIONS),PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS),"",false);
  tab_fprintf(stream,"        --> Stalls                   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS),NULL,"",false);
  tab_fprintf(stream,"          --> Stalls.allBusy         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS_BUSY),
      PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS),"",false);
  tab_fprintf(stream,"          --> Stalls.freePreempted   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS_NOT_PRIORITY),
      PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS),"",false);
  // High-level I/O
  fprintf(stream,    "[GEM]>Profile.Mapper.io.HighLevel\n");
  tab_fprintf(stream,"  => TIME.Load.Index                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_LOAD_INDEX),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  => TIME.Parse.Input.FASTQ              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FASTA_PARSE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Buffer.Input                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED),PROF_GET_TIMER(GP_MAPPER_ALL));
  mapper_profile_print_output(stream,paired_end,map_output);
}
/*
 * Checks
 */
GEM_INLINE void mapper_profile_print_checks(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Checks\n");
  tab_fprintf(stream,"  --> Num.Reads.Checked              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_NUM_READS),NULL,"reads",true);
  tab_fprintf(stream,"  --> Num.Maps.Checked               ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_NUM_MAPS),NULL,"maps ",true);
  tab_fprintf(stream,"    --> Num.Maps.Incorrect           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_INCORRECT),NULL,"maps ",true);
  tab_fprintf(stream,"    --> Num.Maps.Suboptimal          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_SUBOPTIMAL),NULL,"maps ",true);
  tab_fprintf(stream,"      --> Num.Maps.Distance          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_SUBOPTIMAL_DISTANCE),NULL,"     ",true);
  tab_fprintf(stream,"      --> Num.Maps.Score             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_SUBOPTIMAL_SCORE),NULL,"     ",true);
  tab_fprintf(stream,"      --> Num.Maps.ScoreDiff         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_SUBOPTIMAL_DIFF),NULL,"     ",true);
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
  tab_fprintf(stream,"    => TIME.Filtering.Inexact                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTERING_INEXACT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> Archive.Search.SE.Filtering.Stages\n");
  tab_fprintf(stream,"    --> Filtering.Exact.Mapped                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_FILTERING_EXACT_MAPPED),PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"reads",true);
  tab_fprintf(stream,"    --> Filtering.Inexact.Mapped              ");
    COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_AS_FILTERING_INEXACT_MAPPED),PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"reads",true);
//#define GEM_STATS_SHOW_MISMS_SEARCH_TIME_STATS(FILE)
//  fprintf(FILE, "TIME.Exact.Search               %4.1f (%3.2f %%) n", GET_TIME(SC_EXACT_SEARCH), GET_TIME_PERCENTAGE(SC_EXACT_SEARCH,SC_MISMS_SEARCH));
//  fprintf(FILE, "TIME.Small.Reads                %4.1f (%3.2f %%) n", GET_TIME(SC_SMALL_READS), GET_TIME_PERCENTAGE(SC_SMALL_READS,SC_MISMS_SEARCH));
//  fprintf(FILE, "TIME.Read.Recovery              %4.1f (%3.2f %%) n", GET_TIME(SC_READ_RECOVERY), GET_TIME_PERCENTAGE(SC_READ_RECOVERY,SC_MISMS_SEARCH));
//  fprintf(FILE, "TIME.ATH-filter                 %4.1f (%3.2f %%) n", GET_TIME(SC_ATH), GET_TIME_PERCENTAGE(SC_ATH,SC_MISMS_SEARCH));
//  fprintf(FILE, "  --> ATH.Extract.Profile       %4.1f (%3.2f %%) n", GET_TIME(SC_REGION_PROFILE), GET_TIME_PERCENTAGE(SC_REGION_PROFILE,SC_MISMS_SEARCH));
//  fprintf(FILE, "  --> ATH.ZERO-filter           %4.1f (%3.2f %%) n", GET_TIME(SC_ZERO_FILTER), GET_TIME_PERCENTAGE(SC_ZERO_FILTER,SC_MISMS_SEARCH));
//  fprintf(FILE, "  --> ATH.PMS                   %4.1f (%3.2f %%) n", GET_TIME(SC_PROBING_DELTA), GET_TIME_PERCENTAGE(SC_PROBING_DELTA,SC_MISMS_SEARCH));
//  fprintf(FILE, "  --> ATH.Query.Stage           %4.1f (%3.2f %%) n", GET_TIME(SC_ATH_QUERY), GET_TIME_PERCENTAGE(SC_ATH_QUERY,SC_MISMS_SEARCH));
//  fprintf(FILE, "  --> ATH.Filter.Stage          %4.1f (%3.2f %%) n", GET_TIME(SC_ATH_FILTER), GET_TIME_PERCENTAGE(SC_ATH_FILTER,SC_MISMS_SEARCH));
//  fprintf(FILE, "TIME.PA-filter                  %4.1f (%3.2f %%) n", GET_TIME(SC_PROGRESSIVE), GET_TIME_PERCENTAGE(SC_PROGRESSIVE,SC_MISMS_SEARCH));
//  fprintf(FILE, "  --> PA.Query.Stage            %4.1f (%3.2f %%) n", GET_TIME(SC_PAF_QUERY), GET_TIME_PERCENTAGE(SC_PAF_QUERY,SC_MISMS_SEARCH));
//  fprintf(FILE, "    --> PA.Fails                %4.1f (%3.2f %%) n", GET_TIME(SC_PA_FAIL), GET_TIME_PERCENTAGE(SC_PA_FAIL,SC_MISMS_SEARCH));
//  fprintf(FILE, "      --> PA.Full               %4.1f (%3.2f %%) n", GET_TIME(SC_PA_FULL), GET_TIME_PERCENTAGE(SC_PA_FULL,SC_MISMS_SEARCH));
//  fprintf(FILE, "  --> PA.Filter.Stage           %4.1f (%3.2f %%) n", GET_TIME(SC_PAF_FILTER), GET_TIME_PERCENTAGE(SC_PAF_FILTER,SC_MISMS_SEARCH))
//
//#define GEM_STATS_SHOW_FAST_MAPPING_TIME_STATS(FILE)
//  fprintf(FILE, "TIME.Exact.Search               %4.1f (%3.2f %%) n", GET_TIME(SC_EXACT_SEARCH), GET_TIME_PERCENTAGE(SC_EXACT_SEARCH,SC_MISMS_SEARCH));
//  fprintf(FILE, "TIME.Small.Reads                %4.1f (%3.2f %%) n", GET_TIME(SC_SMALL_READS), GET_TIME_PERCENTAGE(SC_SMALL_READS,SC_MISMS_SEARCH));
//  fprintf(FILE, "TIME.Fast-filter                %4.1f (%3.2f %%) n", GET_TIME(SC_ATH), GET_TIME_PERCENTAGE(SC_ATH,SC_MISMS_SEARCH));
//  fprintf(FILE, "  --> Fast.Extract.Profile      %4.1f (%3.2f %%) n", GET_TIME(SC_REGION_PROFILE), GET_TIME_PERCENTAGE(SC_REGION_PROFILE,SC_MISMS_SEARCH));
//  fprintf(FILE, "  --> Fast.ZERO-filter          %4.1f (%3.2f %%) n", GET_TIME(SC_ZERO_FILTER), GET_TIME_PERCENTAGE(SC_ZERO_FILTER,SC_MISMS_SEARCH));
//  fprintf(FILE, "  --> Fast.Query.Stage          %4.1f (%3.2f %%) n", GET_TIME(SC_FILTER_REGIONS), GET_TIME_PERCENTAGE(SC_FILTER_REGIONS,SC_MISMS_SEARCH));
//  fprintf(FILE, "  --> Fast.Filter.Stage         %4.1f (%3.2f %%) n", GET_TIME(SC_ATH_FILTER), GET_TIME_PERCENTAGE(SC_ATH_FILTER,SC_MISMS_SEARCH))


///* Calls Statistics */
//#define GEM_STATS_SHOW_CALLS_STATS(FILE)
//  fprintf(FILE, "Calls.Mismatched.Search           %10lu  (%3.2f %%) {%1.4f ms/call} n", GET_COUNTER(GSC_ASM_FUNC),
//                                                                                        GET_COUNT_PERCENTAGE(GSC_ASM_FUNC,GSC_ASM_FUNC),
//                                                                                        GET_MS_TIME_PER_CALL(SC_MISMS_SEARCH,GSC_ASM_FUNC));
//  fprintf(FILE, "Calls.Exact.Search                %10lu  (%3.2f %%) {%1.4f ms/call} n", GET_COUNTER(SC_EXACT_SEARCH),
//                                                                                        GET_COUNT_PERCENTAGE(SC_EXACT_SEARCH,GSC_ASM_FUNC),
//                                                                                        GET_MS_TIME_PER_CALL(SC_EXACT_SEARCH,SC_EXACT_SEARCH));
//  fprintf(FILE, "Calls.Small.Reads                 %10lu  (%3.2f %%) {%1.4f ms/call} n", GET_COUNTER(SC_SMALL_READS),
//                                                                                        GET_COUNT_PERCENTAGE(SC_SMALL_READS,GSC_ASM_FUNC),
//                                                                                        GET_MS_TIME_PER_CALL(SC_SMALL_READS,SC_SMALL_READS));
//  fprintf(FILE, "Calls.ATH-filter                  %10lu  (%3.2f %%) {%1.4f ms/call} n", GET_COUNTER(SC_ATH),
//                                                                                        GET_COUNT_PERCENTAGE(SC_ATH,GSC_ASM_FUNC),
//                                                                                        GET_MS_TIME_PER_CALL(SC_ATH,SC_ATH));
//  fprintf(FILE, "  --> ATH.ZERO-filter.Hit         %10lu  (%3.2f %%) (%3.2f filters/call) {%1.4f ms/call} n", GET_COUNTER(SC_ZERO_FILTER),
//                                                                                                             GET_COUNT_PERCENTAGE(SC_ZERO_FILTER,GSC_ASM_FUNC),
//                                                                                                             GET_COUNT_DIV(GSC_ZERO_FILTER_CAND,SC_ZERO_FILTER),
//                                                                                                             GET_MS_TIME_PER_CALL(SC_ZERO_FILTER,SC_ZERO_FILTER));
//  fprintf(FILE, "  --> ATH.PMS                     %10lu  (%3.2f %%) (%3.2f filters/call) {%1.4f ms/call} n", GET_COUNTER(SC_PROBING_DELTA),
//                                                                                                             GET_COUNT_PERCENTAGE(SC_PROBING_DELTA,GSC_ASM_FUNC),
//                                                                                                             GET_COUNT_DIV(GSC_DELTA_PROBE_CAND,SC_PROBING_DELTA),
//                                                                                                             GET_MS_TIME_PER_CALL(SC_PROBING_DELTA,SC_PROBING_DELTA));
//  fprintf(FILE, "    --> ATH.PMS.Hit               %10lu  (%3.2f %%) n", GET_COUNTER(GSC_DELTA_PROBE_HIT),
//                                                                        GET_COUNT_PERCENTAGE(GSC_DELTA_PROBE_HIT,GSC_ASM_FUNC));
//  fprintf(FILE, "  --> ATH.Hit                     %10lu  (%3.2f %%) (%3.2f filters/call) n", GET_COUNTER(GSC_ATH_HIT),
//                                                                                             GET_COUNT_PERCENTAGE(GSC_ATH_HIT,GSC_ASM_FUNC),
//                                                                                             GET_COUNT_DIV(GSC_ATH_FILTER_CAND,GSC_ATH_HIT));
//  fprintf(FILE, "Calls.PA-filter                   %10lu  (%3.2f %%) (%3.2f filters/call) {%1.4f ms/call} n", GET_COUNTER(SC_PROGRESSIVE),
//                                                                                                             GET_COUNT_PERCENTAGE(SC_PROGRESSIVE,GSC_ASM_FUNC),
//                                                                                                             GET_COUNT_DIV(GSC_PAF_FILTER_CAND,SC_PROGRESSIVE),
//                                                                                                             GET_MS_TIME_PER_CALL(SC_PROGRESSIVE,SC_PROGRESSIVE))
//
//#define GEM_STATS_FAST_MAPPING_SCOPE(FILE) {
//  fprintf(FILE, "Fast-mapping scope n");
//  register uint64_t i;
//  for (i=0; i<20; ++i) {
//    fprintf(FILE, "  --> Strata %2lu     %10lu     (%3.2f %%) n", i,
//      GET_COUNTER(GSC_FAST_MAPPING_MCS+i),
//      GET_COUNT_PERCENTAGE(GSC_FAST_MAPPING_MCS+i,GSC_ASM_FUNC));
//  }
//}
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
  tab_fprintf(stream,"      => TIME.Basic.Cases               ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_BASIC),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Small.Reads             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_SMALL_READS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Read.Recovery           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_READ_RECOVERY),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Region.Profile            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_REGION_PROFILE_ADAPTIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Generate.Candidates       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_GENERATE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Process.Candidates        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_PROCESS_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Verifying                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_VERIFICATION),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Select.Matches              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SELECT_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Score.Matches             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SCORE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
}
GEM_INLINE void mapper_profile_print_archive_search_pe(FILE* const stream) {
  fprintf(stream,    "[GEM]>Profile.ArchiveSearch.PE\n");
  tab_fprintf(stream,"  => TIME.Archive.Search.PE                    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Archive.Generate.Candidates        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_GENERATE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Archive.Discard.Filtering.Regions  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_DISCARD_FILTERING_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Archive.Finish.Search              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_FINISH_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL)); // TODO
  tab_fprintf(stream,"    => TIME.Archive.Select.PairedMatches       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SELECT_PE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Archive.Find.Paired.Matches        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_PAIRED_MATCHES_FIND_PAIRS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Archive.Extend.Candidates          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Filtering.Extend.Match           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_MATCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Filtering.Retrieve.Candidates  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_RETRIEVE_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Filtering.Verify.Candidates    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Filtering.Realign.Candidates   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    |> Archive.Discard.Filtering.Regions\n");
  tab_fprintf(stream,"      --> Discarded.Filtering.Regions          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_DISCARD_FILTERING_REGIONS_NOT_CONCORDANT),
                       PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_DISCARD_FILTERING_REGIONS_TOTAL),"regions",true);
  tab_fprintf(stream,"    |> Archive.Extend.Candidates\n");
  tab_fprintf(stream,"      --> Num.Extended.End1                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END1),
                       PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"ext",true);
  tab_fprintf(stream,"        --> Extended.End1.Success              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END1_SUCCESS),
                       PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"ext",true);
  tab_fprintf(stream,"      --> Num.Extended.End2                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END2),
                       PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"ext",true);
  tab_fprintf(stream,"        --> Extended.End1.Success              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_END2_SUCCESS),
                       PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"ext",true);
  tab_fprintf(stream,"      --> Matches.Found                        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_FOUND),NULL,"   ",true);
}
GEM_INLINE void mapper_profile_print_archive_search_group(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Archive.Search.Group\n");
  // Search groups
  tab_fprintf(stream,"  => Archive.Search.Group\n");
  tab_fprintf(stream,"    --> Num.Buffers.Used           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_GROUP_BUFFERS_USED),NULL,"buffers",true);
  // BPM-Buffers
  tab_fprintf(stream,"  => BPM.Buffers\n");
  tab_fprintf(stream,"    --> Queries                    ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_BPM_GPU_BUFFER_USAGE_QUERIES));
  tab_fprintf(stream,"    --> PeqEntries                 ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_BPM_GPU_BUFFER_USAGE_PEQ_ENTRIES));
  tab_fprintf(stream,"    --> Candidates                 ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_BPM_GPU_BUFFER_USAGE_CANDIDATES));
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
  tab_fprintf(stream,"    => TIME.Realign.Accepted.Regions           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
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
  tab_fprintf(stream,"    --> Accepted.Regions                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Accepted.Regions.Length                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS_LENGTH),PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),"regions  ",true);
  }
  //  tab_fprintf(stream,"        --> FC.Regions.Coverage.Extended            ");
  //  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_FC_CANDIDATE_REGIONS_EXT_COVERAGE));
  //  tab_fprintf(stream,"      --> FC.Kmer.Counting.Discarded                ");
  //  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_KMER_COUNTER_FILTER_DISCARDED),PROF_GET_COUNTER(GP_FC_NUM_CANDIDATE_REGIONS),"",true);
}
GEM_INLINE void mapper_profile_print_realign(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Realign\n");
  tab_fprintf(stream,"  => TIME.Realign                                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Realign.Exact                         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHES_ALIGN_EXACT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Realign.Hamming                       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHES_ALIGN_HAMMING),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Realign.Levenshtein                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHES_ALIGN_LEVENSHTEIN),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Realign.SWG                           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHES_ALIGN_SWG),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.SWG.Core.Banded                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_SWG_ALIGN_BANDED),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.SWG.Core.Full                       ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_SWG_ALIGN_FULL),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> Realign\n");
  tab_fprintf(stream,"    --> Candidate.Positions                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"positions",true);
  tab_fprintf(stream,"    --> Candidate.Regions                         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"    --> Candidate.Regions.Accepted                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"      --> Regions.Accepted.Coverage               ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS_COVERAGE));
  tab_fprintf(stream,"      --> Regions.Accepted.Coverage.Extended      ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS_EXT_COVERAGE));
  tab_fprintf(stream,"    --> Candidate.Regions.Chained                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS_CHAINED),PROF_GET_COUNTER(GP_ACCEPTED_REGIONS),"         ",true);
  tab_fprintf(stream,"      --> Candidate.Regions.Length                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),"nt       ",true);
  tab_fprintf(stream,"      --> Accepted.Regions.Length                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS_LENGTH),PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),"nt       ",true);
  tab_fprintf(stream,"      --> Aligned.Regions.SWG.Banded.Length       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SWG_ALIGN_BANDED_LENGTH),PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),"nt       ",true);
//  tab_fprintf(stream,"      --> Aligned.Regions.SWG.Full.Length         ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SWG_ALIGN_FULL_LENGTH),PROF_GET_COUNTER(GP_CANDIDATE_REGION_LENGTH),"nt       ",true);
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
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_ARCHIVE_SELECT_MATCHES),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
}
/*
 * Mapper
 */
GEM_INLINE void mapper_profile_print_mapper_single_end(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {
  // I/O
  mapper_profile_print_io(stream,false,map_output);
  // All
  tab_fprintf(stream,"[GEM]>Profile.Mapper\n");
  tab_fprintf(stream,"  => TIME.Mapper                           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_ALL),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Load.Index                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_LOAD_INDEX),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Parse.Input.FASTQ              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FASTA_PARSE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Archive.Search.SE              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_global_inc();
  mapper_profile_print_output(stream,false,map_output);
  tab_global_dec();
  // Archive Search SE
  mapper_profile_print_archive_search_se(stream);
  // Approximate Search
  mapper_profile_print_approximate_search(stream);
  // Generate Candidates
  mapper_profile_print_generate_candidates(stream);
  // Region Profile
  mapper_profile_print_region_profile_minimal(stream);
  mapper_profile_print_region_profile_delimit(stream);
  // Filtering Verification
  mapper_profile_print_filtering(stream,true);
  if (num_threads==1) {
    // Neighborhood Search
    mapper_profile_print_neighborhood_search(stream);
    // Ranks
    mapper_profile_print_mapper_ranks(stream);
  }
  // Checks
  mapper_profile_print_checks(stream);
}
GEM_INLINE void mapper_profile_print_mapper_single_end_cuda(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {
  // I/O
  mapper_profile_print_io(stream,false,map_output);
  // All
  tab_fprintf(stream,"[GEM]>Profile.Mapper\n");
  tab_fprintf(stream,"  => TIME.Mapper                           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_ALL),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Load.Index                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_LOAD_INDEX),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.CUDA.Thread                    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD),PROF_GET_TIMER(GP_MAPPER_ALL));
  // CUDA
  tab_fprintf(stream,"[GEM]>Profile.CUDA.ArchiveSearch\n");
  tab_fprintf(stream,"  => TIME.CUDA.Thread                        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => 0.TIME.CUDA.Buffers.Init              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_INIT),PROF_GET_TIMER(GP_MAPPER_ALL));
  // Generating
  tab_fprintf(stream,"    => 1.TIME.CUDA.Thread.Generating         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD_GENERATING),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Buffer.Input                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Parse.Input.FASTQ              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FASTA_PARSE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Prepare.Sequence       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PREPARE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Generate.Candidates    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_GENERATE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Copy.Candidates        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_COPY_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  // Verifying
  tab_fprintf(stream,"    => 2.TIME.CUDA.Thread.Verifying {ON GPU}\n");
  tab_fprintf(stream,"      => TIME.CUDA.Send.Delay                ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.CUDA.Duty.Cycle                ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_CHECK_TIME),NULL);
  tab_fprintf(stream,"      => TIME.CUDA.Retrieve.Delay            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES_DELAY),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => 3.TIME.CUDA.Thread.Selecting          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD_SELECTING),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Retrieve.Candidates    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Retrieve.Candidates          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_RETRIEVE_BPM_BUFFER_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Realign                      ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Finish.Search          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_FINISH_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Select.Matches         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SELECT_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_global_inc();tab_global_inc();
  mapper_profile_print_output(stream,false,map_output);
  tab_global_dec();tab_global_dec();
  // Archive Search Groups
  mapper_profile_print_archive_search_group(stream);
  // Archive Search SE
  mapper_profile_print_archive_search_se(stream);
  // Approximate Search
  mapper_profile_print_approximate_search(stream);
  // Generate Candidates
  mapper_profile_print_generate_candidates(stream);
  // Region Profile
  mapper_profile_print_region_profile_minimal(stream);
  mapper_profile_print_region_profile_delimit(stream);
  // Filtering Verification
  mapper_profile_print_filtering(stream,true);
  if (num_threads==1) {
    // Neighborhood Search
    mapper_profile_print_neighborhood_search(stream);
    // Ranks
    mapper_profile_print_mapper_ranks(stream);
  }
  // Checks
  mapper_profile_print_checks(stream);
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
  // I/O
  mapper_profile_print_io(stream,true,map_output);
  // Archive Search
  mapper_profile_print_archive_search_pe(stream);
  // Approximate Search
  mapper_profile_print_approximate_search(stream);
  // Generate Candidates
  mapper_profile_print_generate_candidates(stream);
  // Region Profile
  mapper_profile_print_region_profile_minimal(stream);
  mapper_profile_print_region_profile_delimit(stream);
  // Filtering Verification
  mapper_profile_print_filtering(stream,true);
  if (num_threads==1) {
    // Neighborhood Search
    mapper_profile_print_neighborhood_search(stream);
    // Ranks
    mapper_profile_print_mapper_ranks(stream);
  }
}
#else /* GEM_PROFILE DISABLED */
GEM_INLINE void mapper_profile_print_mapper_single_end(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {}
GEM_INLINE void mapper_profile_print_mapper_single_end_cuda(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {}
GEM_INLINE void mapper_profile_print_mapper_paired_end(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {}
#endif
