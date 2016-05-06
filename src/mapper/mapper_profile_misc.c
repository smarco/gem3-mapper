/*
 * PROJECT: GEMMapper
 * FILE: mapper_profile_misc.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper/mapper_profile_misc.h"

#ifdef GEM_PROFILE /* GEM_PROFILE ENABLED */
/*
 * I/O
 */
void mapper_profile_print_io(FILE* const stream) {
  // High-level I/O
  fprintf(stream,    "[GEM]>Profile.Mapper.IO\n");
  tab_fprintf(stream,"    => TIME.Load.Index                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_LOAD_INDEX),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Parse.Input.FASTQ              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FASTA_PARSE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Buffer.Input                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.BufferedInput.Reload       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          --> Buffer.ReloadBuffer.Read     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_BUFFERED_INPUT_BUFFER_SIZE),NULL,"B",true);
  tab_fprintf(stream,"    => TIME.BufferedOutput.Dump            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_OUTPUT_DUMP),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.OutputFile.WriteBuffer       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_WRITE_BUFFER),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        --> Bytes.written                  ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BYTES_WRITTEN),NULL,"B",true);
  tab_fprintf(stream,"        --> Buffer.requests                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS),NULL,"",false);
  tab_fprintf(stream,"          --> Extensions                   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_EXTENSIONS),PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS),"",false);
  tab_fprintf(stream,"          --> Stalls                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS),NULL,"",false);
  tab_fprintf(stream,"            --> Stalls.allBusy             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS_BUSY),
      PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS),"",false);
  tab_fprintf(stream,"            --> Stalls.freePreempted       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS_NOT_PRIORITY),
      PROF_GET_COUNTER(GP_OUTPUT_BUFFER_REQUESTS_STALLS),"",false);
}
/*
 * Output MAP/SAM
 */
void mapper_profile_print_map_output(FILE* const stream,const bool paired_end) {
  if (paired_end) {
    fprintf(stream,"=> TIME.Output.MAP.PE                  ");
    TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_MAP_PE),PROF_GET_TIMER(GP_MAPPER_ALL));
  } else {
    fprintf(stream,"=> TIME.Output.MAP.SE                  ");
    TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_MAP_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
  }
}
void mapper_profile_print_sam_output(FILE* const stream,const bool paired_end) {
  if (paired_end) {
    fprintf(stream,"=> TIME.Output.SAM.PE                  ");
    TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_SAM_PE),PROF_GET_TIMER(GP_MAPPER_ALL));
  } else {
    fprintf(stream,"=> TIME.Output.SAM.SE                  ");
    TIMER_PRINT(stream,PROF_GET_TIMER(GP_OUTPUT_SAM_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
  }
}
/*
 * Checks
 */
void mapper_profile_print_checks(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.Checks\n");
  tab_fprintf(stream,"  --> Num.Reads.Checked                             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_NUM_READS),PROF_GET_COUNTER(GP_CHECK_NUM_READS),"reads",true);
  tab_fprintf(stream,"  --> Num.Maps.Checked                              ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_NUM_MAPS),PROF_GET_COUNTER(GP_CHECK_NUM_MAPS),"maps ",true);
  tab_fprintf(stream,"    --> Num.Maps.Incorrect                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_INCORRECT),PROF_GET_COUNTER(GP_CHECK_NUM_MAPS),"maps ",true);
  tab_fprintf(stream,"    --> Num.Maps.Suboptimal.Primary                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_PRIMARY_SUBOPTIMAL),PROF_GET_COUNTER(GP_CHECK_NUM_MAPS),"maps ",true);
  tab_fprintf(stream,"      --> Num.Maps.Suboptimal.Primary.Distance      ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_PRIMARY_SUBOPTIMAL_DISTANCE),NULL,"     ",true);
  tab_fprintf(stream,"      --> Num.Maps.Suboptimal.Primary.Score         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_PRIMARY_SUBOPTIMAL_SCORE),PROF_GET_COUNTER(GP_CHECK_PRIMARY_SUBOPTIMAL),"     ",true);
  tab_fprintf(stream,"      --> Num.Maps.Suboptimal.Primary.ScoreDiff     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_PRIMARY_SUBOPTIMAL_DIFF),PROF_GET_COUNTER(GP_CHECK_PRIMARY_SUBOPTIMAL),"     ",true);
  tab_fprintf(stream,"    --> Num.Maps.Suboptimal.Subdominant             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_SUBDOMINANT_SUBOPTIMAL),PROF_GET_COUNTER(GP_CHECK_NUM_MAPS),"maps ",true);
  tab_fprintf(stream,"      --> Num.Maps.Suboptimal.Subdominant.Distance  ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_SUBDOMINANT_SUBOPTIMAL_DISTANCE),NULL,"     ",true);
  tab_fprintf(stream,"      --> Num.Maps.Suboptimal.Subdominant.Score     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_SUBDOMINANT_SUBOPTIMAL_SCORE),PROF_GET_COUNTER(GP_CHECK_SUBDOMINANT_SUBOPTIMAL),"     ",true);
  tab_fprintf(stream,"      --> Num.Maps.Suboptimal.Subdominant.ScoreDiff ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_CHECK_SUBDOMINANT_SUBOPTIMAL_DIFF),PROF_GET_COUNTER(GP_CHECK_SUBDOMINANT_SUBOPTIMAL),"     ",true);
}
/*
 * Efficiency Ratios
 */
void mapper_profile_print_mapper_efficiency_ratios(FILE* const stream) {
  // Efficiency Ratios
  fprintf(stream,    "[GEM]>Profile.Efficiency.Ratios\n");
  tab_fprintf(stream,"  --> Ranks/Read               %10.3f ranks/read\n",
      (float)COUNTER_GET_TOTAL(PROF_GET_RANK(GP_MAPPER_ALL))/
      (float)COUNTER_GET_TOTAL(PROF_GET_COUNTER(GP_MAPPER_NUM_READS)));
  tab_fprintf(stream,"    --> Ranks/Alignment        %10.3f ranks/alg\n",
      (float)COUNTER_GET_TOTAL(PROF_GET_RANK(GP_MAPPER_ALL))/
      (float)COUNTER_GET_TOTAL(PROF_GET_COUNTER(GP_MATCHES_SE_ADD_NUM_MAPS)));
  tab_fprintf(stream,"  --> Time/Read                %10.3f ms\n",
      (float)TIMER_GET_TOTAL_MS(PROF_GET_TIMER(GP_MAPPER_ALL))/
      (float)COUNTER_GET_TOTAL(PROF_GET_COUNTER(GP_MAPPER_NUM_READS)));
  tab_fprintf(stream,"    --> Time/Alignment         %10.3f us\n",
      (float)TIMER_GET_TOTAL_US(PROF_GET_TIMER(GP_MAPPER_ALL))/
      (float)COUNTER_GET_TOTAL(PROF_GET_COUNTER(GP_MATCHES_SE_ADD_NUM_MAPS)));
  tab_fprintf(stream,"  --> Throughput\n");
  const float reads_per_sec =
      (float)COUNTER_GET_TOTAL(PROF_GET_COUNTER(GP_MAPPER_NUM_READS)) /
      (float)TIMER_GET_TOTAL_S(PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  --> Reads/s            %10.3f\n",reads_per_sec);
  tab_fprintf(stream,"  --> MegaReads/h        %10.3f\n",reads_per_sec*3600.0/1000000.0);
  tab_fprintf(stream,"  --> GigaReads/d        %10.3f\n",reads_per_sec*3600.0*24.0/1000000000.0);
}
//
//void mapper_profile_print_mem_structs(FILE* const stream) {
//  tab_fprintf(stream,"[GEM]>Profile.MEM\n");
//  // index_mem()
//  // In buffers (total, used)
//  // Out buffers (total, used)
//  // mm_stack
//  // mm_slabs => mm_pool
//  // bpm_buffers
//  // filtering_candidates_vectors
//  // interval_sets
//  // text_collection
//  // matches
//}
#else /* GEM_PROFILE DISABLED */
/*
 * I/O
 */
void mapper_profile_print_io(FILE* const stream) {}
/*
 * Output MAP/SAM
 */
void mapper_profile_print_map_output(FILE* const stream,const bool paired_end) {}
void mapper_profile_print_sam_output(FILE* const stream,const bool paired_end) {}
/*
 * Checks
 */
void mapper_profile_print_checks(FILE* const stream) {}
/*
 * Efficiency Ratios
 */
void mapper_profile_print_mapper_efficiency_ratios(FILE* const stream) {}
#endif
