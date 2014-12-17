/*
 * PROJECT: GEMMapper
 * FILE: mapper_profile.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper_profile.h"

#ifndef GEM_NOPROFILE /* GEM_PROFILE ENABLED */

/*
 * System
 */
GEM_INLINE void mapper_profile_print_io(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.IO\n");
  // IN
  tab_fprintf(stream,"  => Input\n");
  tab_fprintf(stream,"    => TIME.BufferedInput.Reload     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.InputFile.ReloadBuffer ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FILL_BUFFER),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        --> Buffer.ReloadBuffer.Read ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BUFFERED_INPUT_BUFFER_SIZE),NULL,"B",true);
  // OUT
  tab_fprintf(stream,"  => Output\n");
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
}
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
 * Global Mapper
 */
GEM_INLINE void mapper_profile_print_mapper_adaptive(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Mapper\n");
  tab_fprintf(stream,"  => TIME.Mapper                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_ALL),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Load.Index               ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_LOAD_INDEX),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Parse.Input.FASTQ        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FASTA_PARSE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Buffer.IO              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.BufferedInput.Reload ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.BufferedOutput.Dump  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_OUTPUT_DUMP),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Archive.Search           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Basic.Cases            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_BASIC),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Region.Profile         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_REGION_PROFILE_ADAPTIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Generate.Candidates    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTER_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Verifying              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_VERIFICATION),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Select.Matches         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SELECT_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Output.MAP.SE            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_MAP_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
}
GEM_INLINE void mapper_profile_print_mapper_adaptive_ranks(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Mapper\n");
  tab_fprintf(stream,"  => RANKS.Mapper                  ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_MAPPER_ALL),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"    => RANKS.Archive.Search        ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_ARCHIVE_SEARCH_SE),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"      => RANKS.Region.Profile      ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_REGION_PROFILE_ADAPTIVE),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"      => RANKS.Generate.Candidates ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_AS_FILTER_REGIONS),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"      => RANKS.Verify.Candidates   ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_FC_VERIFICATION),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"        => RANKS.Decode.Positions  ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_FC_DECODE_POSITIONS),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
  tab_fprintf(stream,"    => RANKS.Select.Matches        ");
  COUNTER_PRINT(stream,PROF_GET_RANK(GP_ARCHIVE_SELECT_MATCHES),PROF_GET_RANK(GP_MAPPER_ALL),"ranks",true);
}
GEM_INLINE void mapper_profile_print_mapper_cuda_adaptive(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.CUDAMapper\n");
  // General
  tab_fprintf(stream,"  => TIME.Mapper                               ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_ALL),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Load.Index                         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_LOAD_INDEX),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.CUDA.Init                          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_INIT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.CUDA.Thread                        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.CUDA.Buffers.Init                ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_INIT),PROF_GET_TIMER(GP_MAPPER_ALL));
  // Generating
  tab_fprintf(stream,"      => TIME.CUDA.Thread.Generating           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD_GENERATING),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Parse.Input.FASTQ              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FASTA_PARSE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Buffer.IO                      ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.BufferedInput.Reload         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.BufferedOutput.Dump          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_OUTPUT_DUMP),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Archive.Search                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_GENERATE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.Region.Profile               ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_REGION_PROFILE_ADAPTIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.Generate.Candidates          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_FILTER_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.Process.Candidates           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_PROCESS_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.Decode.Positions             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_DECODE_POSITIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.Compose.Regions              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_COMPOSE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Copy.Candidates                ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_COPY_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  // Verifying
  tab_fprintf(stream,"      => TIME.CUDA.Thread.Verifying\n");
  tab_fprintf(stream,"        => TIME.CUDA.Send.Delay                ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.CUDA.Duty.Cycle                ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_CHECK_TIME),NULL);
  tab_fprintf(stream,"        => TIME.CUDA.Retrieve.Delay            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES_DELAY),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.CUDA.Thread.Selecting            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD_SELECTING),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Retrieve.Realign.Candidates    ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Select.Matches                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SELECT_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Output.MAP.SE                  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_MAP_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
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
 * Region Profile
 */
GEM_INLINE void mapper_profile_print_region_profile_soft(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Region.Profile {SOFT}\n");
  tab_fprintf(stream,"  --> Num.Regions              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_SOFT_NUM_REGIONS),NULL,"    ",true);
  tab_fprintf(stream,"  --> Region.length            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_SOFT_REGION_LENGTH),NULL,"nt  ",true);
  tab_fprintf(stream,"  --> Region.candidates        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_SOFT_REGION_CANDIDATES),NULL,"cand",true);
  tab_fprintf(stream,"  --> Read.candidates          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_REGION_PROFILE_SOFT_TOTAL_CANDIDATES),NULL,"cand",true);
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
GEM_INLINE void mapper_profile_print_filtering_verifying(FILE* const stream,const bool verification_profile) {
  tab_fprintf(stream,"[GEM]>Profile.Filtering\n");
  // Verifying
  tab_fprintf(stream,"  => TIME.Verifying                            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_VERIFICATION),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Process.Candidates                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_PROCESS_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Decode.Positions                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_DECODE_POSITIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Compose.Regions                  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_COMPOSE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  if (verification_profile) {
  tab_fprintf(stream,"    => TIME.Retrieve.Candidate.Regions         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_RETRIEVE_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Verify.Candidate.Regions           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_VERIFY_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"  => TIME.Kmer.Counting                      ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_KMER_COUNTER_FILTER),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.BPM.Align                        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_TILED),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Realign.Accepted.Regions           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  }
  tab_fprintf(stream,"    => Filtering.Candidates\n");
  tab_fprintf(stream,"      --> Candidate.Positions                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"positions",true);
  tab_fprintf(stream,"        --> FM.lookup.dist                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FMIDX_LOOKUP_DIST),NULL,"lf       ",true);
  tab_fprintf(stream,"      --> Candidate.Regions                         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_POSITIONS),"regions  ",true);
  tab_fprintf(stream,"        --> Candidate.Tiles                         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BMP_TILED_NUM_TILES),PROF_GET_COUNTER(GP_BMP_TILED_NUM_TILES),"tiles    ",true);
  tab_fprintf(stream,"          --> Candidate.Tiles.Verified              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BMP_TILED_NUM_TILES_VERIFIED),PROF_GET_COUNTER(GP_BMP_TILED_NUM_TILES),"tiles    ",true);
  if (verification_profile) {
  tab_fprintf(stream,"          --> BPM.Accepted                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_LEVENSHTEIN_ACCEPTED),PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),"         ",true);
  tab_fprintf(stream,"            --> BPM.Quick-Abandon                   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BPM_QUICK_ABANDON),PROF_GET_COUNTER(GP_BMP_TILED_NUM_TILES_VERIFIED),"         ",true);
  }

  //  tab_fprintf(stream,"        --> FC.Regions.Coverage.Extended            ");
  //  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_FC_CANDIDATE_REGIONS_EXT_COVERAGE));
  //  tab_fprintf(stream,"      --> FC.Kmer.Counting.Discarded                ");
  //  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_KMER_COUNTER_FILTER_DISCARDED),PROF_GET_COUNTER(GP_FC_NUM_CANDIDATE_REGIONS),"",true);

  // Realign
  tab_fprintf(stream,"  => TIME.Realign                                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    --> Regions.Accepted                            ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS),PROF_GET_COUNTER(GP_CANDIDATE_REGIONS),"",true);
  tab_fprintf(stream,"      --> Regions.Accepted.Coverage                 ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS_COVERAGE));
  tab_fprintf(stream,"      --> Regions.Accepted.ExtCoverage              ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS_EXT_COVERAGE));
  tab_fprintf(stream,"      --> Regions.Chained                           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ACCEPTED_REGIONS_CHAINED),PROF_GET_COUNTER(GP_ACCEPTED_REGIONS),"",true);
  tab_fprintf(stream,"    => TIME.Realign.Exact                           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHES_ALIGN_EXACT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Realign.Hamming                         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHES_ALIGN_HAMMING),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Realign.Levenshtein                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHES_ALIGN_LEVENSHTEIN),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Realign.SWG                             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MATCHES_ALIGN_SWG),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      --> TIME.Realign.SWG.Cells.Computed           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_MATCHES_ALIGN_SWG_CELLS),
      PROF_GET_COUNTER(GP_MATCHES_ALIGN_SWG_CELLS_POTENTIAL),"cells",true);
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
  tab_fprintf(stream,"[GEM]>Profile.Archive.Search.Group\n");
  tab_fprintf(stream,"  => Archive.Search.Group\n");
  tab_fprintf(stream,"    --> Num.Buffers.Used           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_GROUP_BUFFERS_USED),NULL,"buffers",true);
  /*
   * BPM-Buffers
   */
  tab_fprintf(stream,"  => BPM.Buffers\n");
  tab_fprintf(stream,"    --> Queries                    ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_BPM_GPU_BUFFER_USAGE_QUERIES));
  tab_fprintf(stream,"    --> PeqEntries                 ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_BPM_GPU_BUFFER_USAGE_PEQ_ENTRIES));
  tab_fprintf(stream,"    --> Candidates                 ");
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

#else /* GEM_PROFILE DISABLED */
GEM_INLINE void mapper_profile_print_io(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_mem_structs(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_mapper_adaptive(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_mapper_adaptive_ranks(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_mapper_cuda_adaptive(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_mapper_efficiency_ratios(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_approximate_search(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_approximate_search_ranks(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_region_profile_soft(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_filtering_generating(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_filtering_generating_ranks(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_filtering_verifying(FILE* const stream,const bool verification_profile) {}
GEM_INLINE void mapper_profile_print_filtering_verifying_ranks(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_neighborhood_search(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_neighborhood_search_ranks(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_archive_search(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_archive_select(FILE* const stream) {}
GEM_INLINE void mapper_profile_print_archive_search_group(FILE* const stream) {}
#endif
