/*
 * PROJECT: GEMMapper
 * FILE: mapper_profile.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper_profile.h"

/*
 * System
 */
GEM_INLINE void mapper_profile_print_system_info(FILE* const stream) {
  // TODO
//  host_print();
//  proc_print();
//  mem_print();
//  disk_print();
}
GEM_INLINE void mapper_profile_print_io(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.IO\n");
  // IN
  tab_fprintf(stream,"  => Input\n");
  tab_fprintf(stream,"    => TIME.InputFile.FillBuffer   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FILL_BUFFER),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.BufferedInput.Reload   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      --> Buffer.Reload.Read       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BUFFERED_INPUT_BUFFER_SIZE),NULL,"B",true);
  // OUT
  tab_fprintf(stream,"  => Output\n");
  tab_fprintf(stream,"    => TIME.OutputFile.WriteBuffer  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_WRITE_BUFFER),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      --> Bytes.written             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BYTES_WRITTEN),NULL,"B",true);
  tab_fprintf(stream,"      --> Buffer.requests           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS),NULL,"",false);
  tab_fprintf(stream,"        --> Extensions              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_EXTENSIONS),PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS),"",false);
  tab_fprintf(stream,"        --> Stalls                  ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS),NULL,"",false);
  tab_fprintf(stream,"          --> Stalls.all.busy       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS_BUSY),
      PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS),"",false);
  tab_fprintf(stream,"          --> Stalls.low.prio       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS_NOT_PRIORITY),
      PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS),"",false);
}
/*
 * Global Mapper
 */
GEM_INLINE void mapper_profile_print_mapper_adaptive(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Mapper\n");
  tab_fprintf(stream,"  => TIME.Mapper              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_ALL),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Load.Index        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_LOAD_INDEX),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Archive.Search    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Region.Profile  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_REGION_PROFILE_SOFT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Generating      ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTER_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Verifying       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_VERIFY),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Decoding      ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_DECODE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          --> FM.lookup.dist  ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FMIDX_LOOKUP_DIST),NULL,"lf",true);
  tab_fprintf(stream,"        => TIME.Check         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_CHECK),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Select.Matches    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SELECT_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Output.MAP.SE     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_MAP_SE),PROF_GET_TIMER(GP_MAPPER_ALL));

  ///* Time Statistics */
  //#define GEM_STATS_SHOW_GENERAL_TIME_STATS(FILE,NAME,SC_GEM_APP)
  //  fprintf(FILE, "TIME.GEM.Full.%s           %4.1f (%3.2f %%) n",NAME, GET_TIME(SC_GEM_APP), GET_TIME_PERCENTAGE(SC_GEM_APP,SC_GEM_APP));
  //  fprintf(FILE, "  --> Load Index                %4.1f (%3.2f %%) n", GET_TIME(TSC_GEM_LOAD_INDEX), GET_TIME_PERCENTAGE(TSC_GEM_LOAD_INDEX,SC_GEM_APP));
  //  fprintf(FILE, "  --> In Open                   %4.1f (%3.2f %%) n", GET_TIME(TSC_GEM_IN_OPEN), GET_TIME_PERCENTAGE(TSC_GEM_IN_OPEN,SC_GEM_APP));
  //  fprintf(FILE, "  --> Out Open                  %4.1f (%3.2f %%) n", GET_TIME(TSC_GEM_OUT_OPEN), GET_TIME_PERCENTAGE(TSC_GEM_OUT_OPEN,SC_GEM_APP));
  //  fprintf(FILE, "TIME.Input                      %4.1f (%3.2f %%) n", GET_TIME(TSC_INPUT), GET_TIME_PERCENTAGE(TSC_INPUT,SC_GEM_APP));
  //  fprintf(FILE, "TIME.Output                     %4.1f (%3.2f %%) n", GET_TIME(TSC_OUTPUT), GET_TIME_PERCENTAGE(TSC_OUTPUT,SC_GEM_APP));
  //  fprintf(FILE, "  --> DumpOutput                %4.1f (%3.2f %%) n", GET_TIME(TSC_DUMP_OUTPUT), GET_TIME_PERCENTAGE(TSC_DUMP_OUTPUT,SC_GEM_APP))

  //#define GEM_STATS_SHOW_SE_TIME_STATS(FILE)
  //  fprintf(FILE, "TIME.Mapper.Thread              %4.1f (%3.2f %%) n", GET_TIME(SC_MAPPER_TH), GET_TIME_PERCENTAGE(SC_MAPPER_TH,SC_GEM_MAPPER));
  //  fprintf(FILE, "  --> Mismatched.Search         %4.1f (%3.2f %%) n", GET_TIME(SC_MISMS_SEARCH), GET_TIME_PERCENTAGE(SC_MISMS_SEARCH,SC_GEM_MAPPER));
  //  fprintf(FILE, "  --> Decode.Matches            %4.1f (%3.2f %%) n", GET_TIME(SC_DECODE_MATCHES), GET_TIME_PERCENTAGE(SC_DECODE_MATCHES,SC_GEM_MAPPER))
}
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
}
/*
 * Approximate string search
 */
GEM_INLINE void mapper_profile_print_approximate_search(FILE* const stream) {
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
GEM_INLINE void mapper_profile_print_approximate_search_ranks(FILE* const stream) {
  ///* Rank Statistics */
  //#define GEM_STATS_SHOW_RANK_STATS(FILE)
  //  fprintf(FILE, "Ranks.GEM.Mapper.Full           %16lu     (%3.2f %%) n", GET_RANK_COUNTER(SC_GEM_MAPPER), GET_RANK_PERCENTAGE(SC_GEM_MAPPER,SC_GEM_MAPPER));
  //  fprintf(FILE, "Ranks.Mapper.Thread             %16lu     (%3.2f %%) n", GET_RANK_COUNTER(SC_MAPPER_TH), GET_RANK_PERCENTAGE(SC_MAPPER_TH,SC_GEM_MAPPER));
  //  fprintf(FILE, "Ranks.Mismatched.Search         %16lu     (%3.2f %%) n", GET_RANK_COUNTER(SC_MISMS_SEARCH), GET_RANK_PERCENTAGE(SC_MISMS_SEARCH,SC_GEM_MAPPER));
  //  fprintf(FILE, "Ranks.Decode.Matches            %16lu     (%3.2f %%) n", GET_RANK_COUNTER(SC_DECODE_MATCHES), GET_RANK_PERCENTAGE(SC_DECODE_MATCHES,SC_GEM_MAPPER));
  //  fprintf(FILE, "Ranks.Exact.Search              %16lu     (%3.2f %%) n", GET_RANK_COUNTER(SC_EXACT_SEARCH), GET_RANK_PERCENTAGE(SC_EXACT_SEARCH,SC_MISMS_SEARCH));
  //  fprintf(FILE, "Ranks.Small.Reads               %16lu     (%3.2f %%) n", GET_RANK_COUNTER(SC_SMALL_READS), GET_RANK_PERCENTAGE(SC_SMALL_READS,SC_MISMS_SEARCH));
  //  fprintf(FILE, "Ranks.ATH-filter                %16lu     (%3.2f %%) n", GET_RANK_COUNTER(SC_ATH), GET_RANK_PERCENTAGE(SC_ATH,SC_MISMS_SEARCH));
  //  fprintf(FILE, "  --> ATH.Extract.Profile       %16lu     (%3.2f %%) n", GET_RANK_COUNTER(SC_REGION_PROFILE), GET_RANK_PERCENTAGE(SC_REGION_PROFILE,SC_MISMS_SEARCH));
  //  fprintf(FILE, "  --> ATH.ZERO-filter           %16lu     (%3.2f %%) n", GET_RANK_COUNTER(SC_ZERO_FILTER), GET_RANK_PERCENTAGE(SC_ZERO_FILTER,SC_MISMS_SEARCH));
  //  fprintf(FILE, "  --> ATH.PMS                   %16lu     (%3.2f %%) n", GET_RANK_COUNTER(SC_PROBING_DELTA), GET_RANK_PERCENTAGE(SC_PROBING_DELTA,SC_MISMS_SEARCH));
  //  fprintf(FILE, "  --> ATH.Query.Stage           %16lu     (%3.2f %%) n", GET_RANK_COUNTER(SC_ATH_QUERY), GET_RANK_PERCENTAGE(SC_ATH_QUERY,SC_MISMS_SEARCH));
  //  fprintf(FILE, "  --> ATH.Filter.Stage          %16lu     (%3.2f %%) n", GET_RANK_COUNTER(SC_ATH_FILTER), GET_RANK_PERCENTAGE(SC_ATH_FILTER,SC_MISMS_SEARCH));
  //  fprintf(FILE, "Ranks.PA-filter                 %16lu     (%3.2f %%) n", GET_RANK_COUNTER(SC_PROGRESSIVE), GET_RANK_PERCENTAGE(SC_PROGRESSIVE,SC_MISMS_SEARCH));
  //  fprintf(FILE, "  --> PA.Query.Stage            %16lu     (%3.2f %%) n", GET_RANK_COUNTER(SC_PAF_QUERY), GET_RANK_PERCENTAGE(SC_PAF_QUERY,SC_MISMS_SEARCH));
  //  fprintf(FILE, "  --> PA.Filter.Stage           %16lu     (%3.2f %%) n", GET_RANK_COUNTER(SC_PAF_FILTER), GET_RANK_PERCENTAGE(SC_PAF_FILTER,SC_MISMS_SEARCH))
  //
}
/*
 * Filtering Generating
 */
GEM_INLINE void mapper_profile_print_filtering_generating(FILE* const stream) {
  ///* Events Statistics */
  //#define GEM_STATS_FILTERING_QUERY_REGIONS(FILE)
  //  fprintf(FILE, "CallsATH.filter.regions           %10lu     (%3.2f %%) n", GET_COUNTER(SC_FILTER_REGIONS), GET_COUNT_PERCENTAGE(SC_FILTER_REGIONS,SC_FILTER_REGIONS));
  //  register const uint64_t total_regions = GET_COUNTER(GSC_NUM_HARD_REGIONS)+GET_COUNTER(GSC_NUM_SOFT_REGIONS);
  //  fprintf(FILE, "  --> NumRegions                  %10lu     (%3.2f) n", total_regions, 100.0);
  //  fprintf(FILE, "    --> SoftRegions               %10lu     (%3.2f avg) n", GET_COUNTER(GSC_NUM_HARD_REGIONS), GET_COUNT_DIV(GSC_NUM_HARD_REGIONS,SC_ATH));
  //  fprintf(FILE, "    --> HardRegions               %10lu     (%3.2f avg) n", GET_COUNTER(GSC_NUM_SOFT_REGIONS), GET_COUNT_DIV(GSC_NUM_SOFT_REGIONS,SC_ATH));
  //  fprintf(FILE, "    --> Profiles.Quit             %10lu     (%3.2f %%) n", GET_COUNTER(GSC_QUIT_PROFILE), GET_COUNT_PERCENTAGE(GSC_QUIT_PROFILE,SC_REGION_PROFILE));
  //  fprintf(FILE, "  --> Filter.Reg.Potential        %10lu     (%3.2f %%) n", GET_COUNTER(GSC_CAND_FILTER_REGIONS), GET_COUNT_PERCENTAGE(GSC_CAND_FILTER_REGIONS,GSC_CAND_FILTER_REGIONS));
  //  fprintf(FILE, "    --> Num.Filtered.Reg          %10lu     (%3.2f %%) n", GET_COUNTER(GSC_FILTER_REGIONS), GET_COUNT_PERCENTAGE(GSC_FILTER_REGIONS,GSC_CAND_FILTER_REGIONS));
  //  fprintf(FILE, "    --> Dyn.Saved.Reg             %10lu     (%3.2f %%) n", GET_COUNTER(GSC_SAVED_FILTER_REGIONS), GET_COUNT_PERCENTAGE(GSC_SAVED_FILTER_REGIONS,GSC_CAND_FILTER_REGIONS));
  //  fprintf(FILE, "  --> ATH.Filters.d>=2            %10lu     (%3.2f %%) n", GET_COUNTER(GSC_ATH_D2), GET_COUNT_PERCENTAGE(GSC_ATH_D2,GSC_FILTER_REGIONS));
  //  fprintf(FILE, "    --> ATH.Filters.d2.Hit        %10lu     (%3.2f %%) (%3.2f %%)n", GET_COUNTER(GSC_ATH_D2_HIT), GET_COUNT_PERCENTAGE(GSC_ATH_D2_HIT,GSC_ATH_D2), GET_COUNT_PERCENTAGE(GSC_ATH_D2_HIT,GSC_FILTER_REGIONS));
  //  fprintf(FILE, "  --> ATH.Filters.d==1            %10lu     (%3.2f %%) n", GET_COUNTER(GSC_ATH_D1), GET_COUNT_PERCENTAGE(GSC_ATH_D1,GSC_FILTER_REGIONS));
  //  fprintf(FILE, "    --> ATH.Filters.d1.Hit        %10lu     (%3.2f %%) (%3.2f %%)n", GET_COUNTER(GSC_ATH_D1_HIT), GET_COUNT_PERCENTAGE(GSC_ATH_D1_HIT,GSC_ATH_D1), GET_COUNT_PERCENTAGE(GSC_ATH_D1_HIT,GSC_FILTER_REGIONS));
  //  fprintf(FILE, "  --> ATH.Filters.d==0            %10lu     (%3.2f %%) (%3.2f %%) n", GET_COUNTER(GSC_ATH_D0), GET_COUNT_PERCENTAGE(GSC_ATH_D0,GSC_ATH_D0), GET_COUNT_PERCENTAGE(GSC_ATH_D0,GSC_FILTER_REGIONS))
  //
}
GEM_INLINE void mapper_profile_print_filtering_generating_ranks(FILE* const stream) {

}
/*
 * Filtering Verification
 */
GEM_INLINE void mapper_profile_print_filtering_verifying(FILE* const stream) {
  //#define GEM_STATS_FILTERING_CHECK_REGIONS(FILE)
  //  fprintf(FILE, "Mappings                          %10lu     (%3.2f %%) n", GET_COUNTER(GSC_NUM_MAPS), GET_COUNT_PERCENTAGE(GSC_NUM_MAPS,GSC_NUM_MAPS));
  //  fprintf(FILE, "  --> FilterCandidates            %10lu     (%3.2f %%) n", GET_COUNTER(GSC_CHECKED_MATCHES), GET_COUNT_PERCENTAGE(GSC_CHECKED_MATCHES,GSC_CHECKED_MATCHES));
  //  fprintf(FILE, "    --> HammingHits               %10lu     (%3.2f %%) n", GET_COUNTER(GSC_HAMMING_HIT), GET_COUNT_PERCENTAGE(GSC_HAMMING_HIT,GSC_CHECKED_MATCHES));
  //  fprintf(FILE, "    --> LevenshteinHits           %10lu     (%3.2f %%) n", GET_COUNTER(GSC_LEVENSHTEIN), GET_COUNT_PERCENTAGE(GSC_LEVENSHTEIN,GSC_CHECKED_MATCHES));
  //  fprintf(FILE, "      --> AverageErrors           %10lu     (%3.2f)    n", GET_COUNTER(GSC_ONLINE_ASM_ERRORS), GET_COUNT_DIV(GSC_ONLINE_ASM_ERRORS,GSC_LEVENSHTEIN));
  //  fprintf(FILE, "      --> AverageMisms            %10lu     (%3.2f)    n", GET_COUNTER(GSC_ONLINE_ASM_MISMS), GET_COUNT_DIV(GSC_ONLINE_ASM_MISMS,GSC_LEVENSHTEIN));
  //  fprintf(FILE, "    --> BigSingleIndelChecks      %10lu     (%3.2f %%) n", GET_COUNTER(GSC_CHECK_BIG_INDEL), GET_COUNT_PERCENTAGE(GSC_CHECK_BIG_INDEL,GSC_CHECKED_MATCHES));
  //  fprintf(FILE, "      --> BigSingleIndelHits      %10lu     (%3.2f %%) (%3.2f %%) n", GET_COUNTER(GSC_BIG_INDEL), GET_COUNT_PERCENTAGE(GSC_BIG_INDEL,GSC_CHECK_BIG_INDEL), GET_COUNT_PERCENTAGE(GSC_BIG_INDEL,GSC_CHECKED_MATCHES));
  //  fprintf(FILE, "  --> CallsCheckCandidates        %10lu     (%3.2f %%) n", GET_COUNTER(SC_CHECK_CAND), GET_COUNT_PERCENTAGE(SC_CHECK_CAND,SC_CHECK_CAND));
  //  fprintf(FILE, "    --> UniqueHits                %10lu     (%3.2f %%) n", GET_COUNTER(GSC_UNIQUE_SHORTCUT), GET_COUNT_PERCENTAGE(GSC_UNIQUE_SHORTCUT,SC_CHECK_CAND));
  //  fprintf(FILE, "      --> CheckSaved              %10lu     (%3.2f %%) n", GET_COUNTER(GSC_UNIQUE_SAVED), GET_COUNT_PERCENTAGE(GSC_UNIQUE_SAVED,GSC_CHECKED_MATCHES))

  //#define GEM_STATS_SHOW_FILTERING_TIME_STATS(FILE)
  //  fprintf(FILE, "TIME.CheckCandidates            %4.1f (%3.2f %%) n", GET_TIME(SC_CHECK_CAND), GET_TIME_PERCENTAGE(SC_CHECK_CAND,SC_MISMS_SEARCH));
  //  fprintf(FILE, "  --> CC.Decode                 %4.1f (%3.2f %%) n", GET_TIME(SC_CHECK_DECODE), GET_TIME_PERCENTAGE(SC_CHECK_DECODE,SC_MISMS_SEARCH));
  //  fprintf(FILE, "  --> CC.Hamming                %4.1f (%3.2f %%) n", GET_TIME(TSC_CHECK_HAMMING), GET_TIME_PERCENTAGE(TSC_CHECK_HAMMING,SC_MISMS_SEARCH));
  //  fprintf(FILE, "  --> CC.Levenshtein            %4.1f (%3.2f %%) n", GET_TIME(TSC_CHECK_LEVENSHTEIN), GET_TIME_PERCENTAGE(TSC_CHECK_LEVENSHTEIN,SC_MISMS_SEARCH));
  //  fprintf(FILE, "  --> CC.BigSingleIndel         %4.1f (%3.2f %%) n", GET_TIME(TSC_CHECK_BIG_INDEL), GET_TIME_PERCENTAGE(TSC_CHECK_BIG_INDEL,SC_MISMS_SEARCH))
  //
}
GEM_INLINE void mapper_profile_print_filtering_verifying_ranks(FILE* const stream) {

}
/*
 * Neighborhood Search
 */
GEM_INLINE void mapper_profile_print_neighborhood_search(FILE* const stream) {
  //#define GEM_STATS_TLS_REGIONS(FILE)
  //  float total_tls_time = GET_TIME(TSC_TLS_LOOK_UP)+GET_TIME(TSC_TLS_LOOK_UP_SOLVE)+GET_TIME(TSC_TLS_RANK);
  //  fprintf(FILE, "TLS                         %4.1f n", total_tls_time);
  //  fprintf(FILE, "  --> TLS.lookup            %4.1f     n", GET_TIME(TSC_TLS_LOOK_UP));
  //  fprintf(FILE, "  --> TLS.solve.lookup      %4.1f     (%10lu avg) n", GET_TIME(TSC_TLS_LOOK_UP_SOLVE), GET_COUNTER(GSC_TLS_CALLS)?GET_COUNTER(GSC_TLS_INTERVALS)/GET_COUNTER(GSC_TLS_CALLS):0);
  //  fprintf(FILE, "  --> TLS.rank              %4.1f     n", GET_TIME(TSC_TLS_RANK));
  //
}
GEM_INLINE void mapper_profile_print_neighborhood_search_ranks(FILE* const stream) {

}
/*
 * Archive Search
 */
GEM_INLINE void mapper_profile_print_archive_search(FILE* const stream) {

}
/*
 * Archive Select
 */
GEM_INLINE void mapper_profile_print_archive_select(FILE* const stream) {

}
/*
 * Archive Search-Group (Dispatcher, BMP-Buffers, ...)
 */
GEM_INLINE void mapper_profile_print_archive_search_group(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.IO\n");
  /*
   * Dispatcher
   */
  // Generating
  tab_fprintf(stream,"  => Dispatcher\n");
  tab_fprintf(stream,"    --> SearchGroup.GENERATING\n");
  tab_fprintf(stream,"        --> Requests                   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_GENERATING),NULL,"",false);
  tab_fprintf(stream,"        --> Extensions                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_GENERATING_EXTENSION),
      PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_GENERATING),"",false);
  tab_fprintf(stream,"        --> Stalls                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_GENERATING_STALLS),NULL,"",false);
  tab_fprintf(stream,"          --> Stalls.all.busy          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_GENERATING_STALLS_BUSY),
      PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_GENERATING_STALLS),"",false);
  tab_fprintf(stream,"          --> Stalls.low.prio          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_GENERATING_STALLS_NOT_PRIORITY),
      PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_GENERATING_STALLS),"",false);
  // Selecting
  tab_fprintf(stream,"    --> SearchGroup.SELECTING\n");
  tab_fprintf(stream,"        --> Requests                   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING),NULL,"",false);
  tab_fprintf(stream,"        --> Extensions                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING_EXTENSION),
      PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING),"",false);
  tab_fprintf(stream,"        --> Stalls                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS),NULL,"",false);
  tab_fprintf(stream,"          --> Stalls.idle              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS_IDLE),
      PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS),"",false);
  tab_fprintf(stream,"          --> Stalls.noSingleGroups    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS_NO_SINGLE_GROUPS),
      PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS),"",false);
  tab_fprintf(stream,"          --> Stalls.extensionNotReady ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS_EXTENSION_NOT_READY),
      PROF_GET_COUNTER(GP_SGDISPATCHER_REQUESTS_SELECTING_STALLS),"",false);
  /*
   * BPM-Buffers
   */
  tab_fprintf(stream,"  => BPM.Buffers\n");
  tab_fprintf(stream,"      --> TIME.Send.Buffer             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      --> TIME.Check.Buffer            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_CHECK_TIME),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      --> Buffer.usage\n");
  tab_fprintf(stream,"        --> Queries                    ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_BPM_GPU_BUFFER_USAGE_QUERIES));
  tab_fprintf(stream,"        --> PeqEntries                 ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_BPM_GPU_BUFFER_USAGE_PEQ_ENTRIES));
  tab_fprintf(stream,"        --> Candidates                 ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_BPM_GPU_BUFFER_USAGE_CANDIDATES));
}

/*
 * GT.Stats like
 */

/*
GEM_INLINE void mapper_profile_print_input_stats(FILE* const stream) {
  //#define GEM_STATS_SHOW_INPUT_STATS(FILE)
  //  fprintf(FILE, "  Input.Statsn");
  //  fprintf(FILE, "    --> Reads %lu n", GET_COUNTER(GSC_NUM_READS));
  //  fprintf(FILE, "    --> Reads average length %lu n", GET_COUNTER(GSC_NUM_BASES)/GET_COUNTER(GSC_NUM_READS))
  //
}
GEM_INLINE void mapper_profile_print_mapping_stats(FILE* const stream) {
  /// Mappings
  //#define GEM_STATS_SHOW_MAPPING_STATS(FILE,TITLE,ARE_PAIRED,GSC_NUM_MAPS,GSC_NUM_PAIRS,
  //    GSC_MAPPED,GSC_0_UNIQUE_MAPPED,GSC_1_UNIQUE_MAPPED,GSC_2_UNIQUE_MAPPED)
  //  fprintf(FILE, "MAPPING.STATS.%sn",TITLE);
  //  fprintf(FILE, "  Total.Matches            %10lu     n", GET_COUNTER(GSC_NUM_MAPS));
  //  fprintf(FILE, "  Total.Reads.%s           %10lu     (%3.2f %%)n", ARE_PAIRED?"PE":"SE", GET_COUNTER(GSC_NUM_PAIRS), GET_COUNT_PERCENTAGE(GSC_NUM_PAIRS,GSC_NUM_PAIRS));
  //  fprintf(FILE, "    --> Mapped             %10lu     (%3.2f %%)n", GET_COUNTER(GSC_MAPPED), GET_COUNT_PERCENTAGE(GSC_MAPPED,GSC_NUM_PAIRS));
  //  fprintf(FILE, "    --> L0.UniquelyMapped  %10lu     (%3.2f %%)n", GET_COUNTER(GSC_0_UNIQUE_MAPPED), GET_COUNT_PERCENTAGE(GSC_0_UNIQUE_MAPPED,GSC_NUM_PAIRS));
  //  fprintf(FILE, "    --> L1.UniquelyMapped  %10lu     (%3.2f %%)n", GET_COUNTER(GSC_1_UNIQUE_MAPPED), GET_COUNT_PERCENTAGE(GSC_1_UNIQUE_MAPPED,GSC_NUM_PAIRS));
  //  fprintf(FILE, "    --> L2.UniquelyMapped  %10lu     (%3.2f %%)n", GET_COUNTER(GSC_2_UNIQUE_MAPPED), GET_COUNT_PERCENTAGE(GSC_2_UNIQUE_MAPPED,GSC_NUM_PAIRS))
  //
}

#define GEM_STATS_SHOW_PE_TIME_STATS(FILE)
  fprintf(FILE, "TIME.PE.Thread                  %4.1f (%3.2f %%) n", GET_TIME(SC_GEM_PAIR_TH), GET_TIME_PERCENTAGE(SC_GEM_PAIR_TH,SC_GEM_PAIR));
  fprintf(FILE, "  --> PairMatches               %4.1f (%3.2f %%) n", GET_TIME(SC_PAIR_MATCHES), GET_TIME_PERCENTAGE(SC_PAIR_MATCHES,SC_GEM_PAIR));
  fprintf(FILE, "  --> ExtendMatches             %4.1f (%3.2f %%) n", GET_TIME(SC_EXT__PAIR_MATCHES), GET_TIME_PERCENTAGE(SC_EXT__PAIR_MATCHES,SC_GEM_PAIR));
  fprintf(FILE, "    --> OnlineASM               %4.1f (%3.2f %%) n", GET_TIME(SC_ONLINE_ASM_SEARCH), GET_TIME_PERCENTAGE(SC_ONLINE_ASM_SEARCH,SC_GEM_PAIR));
  fprintf(FILE, "      --> Duplicates            %8lu (%3.2f %%) n", GET_COUNTER(GSC_ERASED_DUPLICATES), GET_COUNT_PERCENTAGE(GSC_ERASED_DUPLICATES,SC_PAIR_MATCHES))
#define GEM_STATS_SHOW_SE__PE_TIME_STATS(FILE)
  fprintf(FILE, "TIME.SE.Thread                  %4.1f (%3.2f %%) n", GET_TIME(SC_MAPPER_TH), GET_TIME_PERCENTAGE(SC_MAPPER_TH,SC_GEM_PAIR));
  fprintf(FILE, "  --> Mismatched.Search         %4.1f (%3.2f %%) n", GET_TIME(SC_MISMS_SEARCH), GET_TIME_PERCENTAGE(SC_MISMS_SEARCH,SC_GEM_PAIR));
  fprintf(FILE, "TIME.PE.Thread                  %4.1f (%3.2f %%) n", GET_TIME(SC_GEM_PAIR_TH), GET_TIME_PERCENTAGE(SC_GEM_PAIR_TH,SC_GEM_PAIR));
  fprintf(FILE, "  --> PairMatches               %4.1f (%3.2f %%) n", GET_TIME(SC_PAIR_MATCHES), GET_TIME_PERCENTAGE(SC_PAIR_MATCHES,SC_GEM_PAIR));
  fprintf(FILE, "  --> ExtendMatches             %4.1f (%3.2f %%) n", GET_TIME(SC_EXT__PAIR_MATCHES), GET_TIME_PERCENTAGE(SC_EXT__PAIR_MATCHES,SC_GEM_PAIR));
  fprintf(FILE, "    --> OnlineASM               %4.1f (%3.2f %%) n", GET_TIME(SC_ONLINE_ASM_SEARCH), GET_TIME_PERCENTAGE(SC_ONLINE_ASM_SEARCH,SC_GEM_PAIR));
  fprintf(FILE, "      --> Duplicates            %10lu (%3.2f %%) n", GET_COUNTER(GSC_ERASED_DUPLICATES), GET_COUNT_PERCENTAGE(GSC_ERASED_DUPLICATES,SC_PAIR_MATCHES))

#define GEM_STATS_SHOW_PE_EVENTS_STATS(FILE)
  fprintf(FILE, "Mapped Ends n");
  fprintf(FILE, "  --> MappedEnd.end1                     %10lu     (%3.2f %%) n", GET_COUNTER(GSC_PAIR_END1_MAP), 2.0*GET_COUNT_PERCENTAGE(GSC_PAIR_END1_MAP,GSC_NUM_READS));
  fprintf(FILE, "  --> MappedEnd.end2                     %10lu     (%3.2f %%) n", GET_COUNTER(GSC_PAIR_END2_MAP), 2.0*GET_COUNT_PERCENTAGE(GSC_PAIR_END2_MAP,GSC_NUM_READS));
  fprintf(FILE, "ExtendedMatches n");
  fprintf(FILE, "  --> ExtendedMatches.hits               %10lu     (%3.2f %%) n", GET_COUNTER(GSC_ONLINE_ASM_HIT), GET_COUNT_PERCENTAGE(GSC_ONLINE_ASM_HIT,GSC_ONLINE_ASM_HIT));
  fprintf(FILE, "    --> ExtendedMatches.hitsAdded        %10lu     (%3.2f %%) (%3.2f %%) n", GET_COUNTER(GSC_EXT_MATCHES_ADDED), GET_COUNT_PERCENTAGE(GSC_EXT_MATCHES_ADDED,
                                                                                               GSC_ONLINE_ASM_HIT), GET_COUNT_PERCENTAGE(GSC_EXT_MATCHES_ADDED,GSC_NUM_MAPS_EXT));
  fprintf(FILE, "  --> ExtendedMatches.averageErrors      %10lu     (%3.2f %%) n", GET_COUNTER(GSC_ONLINE_ASM_ERRORS), GET_COUNT_DIV(GSC_ONLINE_ASM_ERRORS,GSC_ONLINE_ASM_HIT));
  fprintf(FILE, "  --> ExtendedMatches.averageMisms       %10lu     (%3.2f %%) n", GET_COUNTER(GSC_ONLINE_ASM_MISMS), GET_COUNT_DIV(GSC_ONLINE_ASM_MISMS,GSC_ONLINE_ASM_HIT));
  fprintf(FILE, "  --> ExtendedMatches.averageInsSize     %10lu     (%3.2f %%) n", GET_COUNTER(GSC_ONLINE_ASM_DISTANCE), GET_COUNT_DIV(GSC_ONLINE_ASM_DISTANCE,GSC_ONLINE_ASM_HIT))

*/
