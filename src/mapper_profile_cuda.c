/*
 * PROJECT: GEMMapper
 * FILE: mapper_profile_cuda.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper_profile_cuda.h"
#include "mapper_profile_search.h"
#include "mapper_profile_misc.h"


#ifdef GEM_PROFILE /* GEM_PROFILE ENABLED */
/*
 * Stage Buffers
 */
void mapper_profile_print_search_stage_buffers(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.CUDA.Search.Stages.Buffers\n");
  tab_fprintf(stream,"  => TIME.GPU.Buffer.Region.Profile\n");
  tab_fprintf(stream,"    => TIME.Alloc                            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_SEARCH_ALLOC),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Send                             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_SEARCH_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Receive                          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_SEARCH_RECEIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Duty.Cycle                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_SEARCH_DUTY_CYCLE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  => TIME.GPU.Buffer.Decode.Candidates\n");
  tab_fprintf(stream,"    => TIME.Alloc                            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_DECODE_ALLOC),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Send                             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_DECODE_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Receive                          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_DECODE_RECEIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Duty.Cycle                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_DECODE_DUTY_CYCLE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  => TIME.GPU.Buffer.Verify.Candidates\n");
  tab_fprintf(stream,"    => TIME.Alloc                            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_ALIGN_BPM_ALLOC),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Send                             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_ALIGN_BPM_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Receive                          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_ALIGN_BPM_RECEIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Duty.Cycle                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_ALIGN_BPM_DUTY_CYCLE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> CUDA.Search.Stages.Buffers\n");
  tab_fprintf(stream,"    --> CUDA.RegionProfile\n");
  tab_fprintf(stream,"      --> Buffers.Used                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SEARCH_STAGE_REGION_PROFILE_BUFFERS_USED),NULL,"",true);
  tab_fprintf(stream,"      --> Buffer.Searches                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SEARCH_STAGE_REGION_PROFILE_SEARCHES_IN_BUFFER),NULL,"",true);
  tab_fprintf(stream,"      --> Buffer.Queries                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_FMI_SEARCH_NUM_QUERIES),NULL,"",true);
  tab_fprintf(stream,"        --> Buffer.Queries.Usage             ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_FMI_SEARCH_USAGE_QUERIES),"        ");
  tab_fprintf(stream,"    --> CUDA.DecodeCandidates\n");
  tab_fprintf(stream,"      --> Buffers.Used                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SEARCH_STAGE_DECODE_CANDIDATES_BUFFERS_USED),NULL,"buffers ",true);
  tab_fprintf(stream,"      --> Buffer.Searches                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SEARCH_STAGE_DECODE_CANDIDATES_SEARCHES_IN_BUFFER),NULL,"searches",true);
  tab_fprintf(stream,"      --> Buffer.Queries                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_FMI_DECODE_NUM_QUERIES),NULL,"queries ",true);
  tab_fprintf(stream,"        --> Buffer.Queries.Usage             ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_FMI_DECODE_USAGE_CANDIDATES),"        ");
  tab_fprintf(stream,"    --> CUDA.VerifyCandidates\n");
  tab_fprintf(stream,"      --> Buffers.Used                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SEARCH_STAGE_VERIFY_CANDIDATES_BUFFERS_USED),NULL,"buffers ",true);
  tab_fprintf(stream,"      --> Buffer.Searches                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SEARCH_STAGE_VERIFY_CANDIDATES_SEARCHES_IN_BUFFER),NULL,"searches",true);
  tab_fprintf(stream,"      --> Buffer.Queries                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_NUM_QUERIES),NULL,"queries ",true);
  tab_fprintf(stream,"        --> Buffer.Queries.Usage             ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_USAGE_QUERIES),"        ");
  tab_fprintf(stream,"        --> Buffer.Entries.Usage             ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_USAGE_PEQ_ENTRIES),"        ");
  tab_fprintf(stream,"      --> Buffer.Candidates                  ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_NUM_QUERIES),NULL,"queries ",true);
  tab_fprintf(stream,"        --> Candidates.Length                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_CANDIDATE_LENGTH),NULL,"nt      ",true);
  tab_fprintf(stream,"        --> Candidates.per.Query             ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_CANDIDATE_PER_QUERY),NULL,"nt      ",true);
  tab_fprintf(stream,"        --> Buffer.Candidates.Usage          ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_ALIGN_BPM_USAGE_CANDIDATES),"        ");
}
/*
 * Archive Search SE CUDA
 */
void mapper_profile_print_archive_search_se_cuda(FILE* const stream,const bool map_output) {
  // CUDA SE Search
  tab_fprintf(stream,"[GEM]>Profile.CUDA.Search\n");
  tab_fprintf(stream,"  => TIME.GPU.Init                                ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_COLLECTION_INIT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  => TIME.CUDA.Search                             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.CUDA.STAGE.Region.Profile             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_SE_REGION_PROFILE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Buffer.Input                        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Parse.Input.FASTQ                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FASTA_PARSE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Init                        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_INIT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Region.Profile.Generate     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_GENERATE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Region.Profile.Copy         ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_COPY),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.CUDA.STAGE.Decode.Candidates          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_SE_DECODE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Region.Profile.Retrieve     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_RETRIEVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        --> Region.Profile.Unsuccessful           ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ASSW_REGION_PROFILE_UNSUCCESSFUL),PROF_GET_COUNTER(GP_MAPPER_NUM_READS),"profiles ",true);
  tab_fprintf(stream,"      => TIME.Archive.Decode.Candidates.Copy      ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_DECODE_CANDIDATES_COPY),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.CUDA.Verify.Candidates                ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_SE_VERIFY_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Decode.Candidates.Retrieve  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_DECODE_CANDIDATES_RETRIEVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        --> Decode.Candidates.Unsuccessful        ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ASSW_DECODE_CANDIDATES_UNSUCCESSFUL),PROF_GET_COUNTER(GP_ASSW_DECODE_CANDIDATES),"positions",true);
  tab_fprintf(stream,"      => TIME.Archive.Verify.Candidates.Copy      ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES_COPY),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.CUDA.FinishSearch                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_SE_FINISH_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Verify.Candidates.Retrieve  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES_RETRIEVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Finish.Search               ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Archive.Select.Matches              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SELECT_SE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      ");
  if (map_output) {
    mapper_profile_print_map_output(stream,false);
  } else {
    mapper_profile_print_sam_output(stream,false);
  }
  // Search Stage Buffers
  mapper_profile_print_search_stage_buffers(stream);
}
/*
 * Archive Search PE CUDA
 */
void mapper_profile_print_archive_search_pe_cuda(FILE* const stream,const bool map_output) {
  /*
   * TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
   */
  tab_fprintf(stream,"[GEM]>Profile.CUDA.ArchiveSearch\n");
//  tab_fprintf(stream,"  => TIME.CUDA.Thread                        ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"    => 0.TIME.CUDA.Buffers.Init              ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_INIT),PROF_GET_TIMER(GP_MAPPER_ALL));
//  // Generating
//  tab_fprintf(stream,"    => 1.TIME.CUDA.Generating                ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD_GENERATING),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Buffer.Input                   ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Parse.Input.FASTQ              ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FASTA_PARSE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Generate.Candidates    ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_GENERATE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Copy.Candidates        ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_COPY_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
//  // Verifying
//  tab_fprintf(stream,"    => 2.TIME.CUDA.Verifying                 ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD_SELECTING),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => {GPU}\n");
//  tab_fprintf(stream,"        => TIME.CUDA.Send.Delay              ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"        => TIME.CUDA.Duty.Cycle              ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BPM_GPU_BUFFER_CHECK_TIME),NULL);
//  tab_fprintf(stream,"        => TIME.CUDA.Retrieve.Delay          ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_GROUP_RETRIEVE_CANDIDATES_DELAY),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => {CPU}\n");
//  tab_fprintf(stream,"        => TIME.Archive.Retrieve.Candidates  ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_RETRIEVE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"          => TIME.Retrieve.Candidates        ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_RETRIEVE_BPM_BUFFER_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"          => TIME.Realign                    ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_REALIGN_BPM_BUFFER_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"        => TIME.Archive.Finish.Search        ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"        => TIME.Restart.Unfit.Searches       ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_THREAD_RESTART_UNFIT),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"        "); mapper_profile_print_output(stream,false,map_output);
//  // Archive Search Groups
//  mapper_profile_print_archive_search_group(stream);
}
/*
 * Mapper-CUDA SE
 */
void mapper_profile_print_mapper_se_cuda(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {
  // Main
  tab_fprintf(stream,"[GEM]>Profile.Mapper\n");
  tab_fprintf(stream,"  => TIME.Mapper                           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_ALL),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Load.Index                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_LOAD_INDEX),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.GPU.Init                       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_COLLECTION_INIT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.CUDA.SE                        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    ");
  if (map_output) {
    mapper_profile_print_map_output(stream,false);
  } else {
    mapper_profile_print_sam_output(stream,false);
  }
  // CUDA
  mapper_profile_print_archive_search_se_cuda(stream,map_output);
  // Approximate Search
  mapper_profile_print_approximate_search_summary(stream,true,map_output,num_threads);
}
/*
 * Mapper-CUDA PE
 */
void mapper_profile_print_mapper_pe_cuda(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {
  // All
  tab_fprintf(stream,"[GEM]>Profile.Mapper\n");
  tab_fprintf(stream,"  => TIME.Mapper                           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_ALL),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Load.Index                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_LOAD_INDEX),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.GPU.Init                       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_COLLECTION_INIT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.CUDA.PE                        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_PE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    ");
  if (map_output) {
    mapper_profile_print_map_output(stream,false);
  } else {
    mapper_profile_print_sam_output(stream,false);
  }
  // CUDA
  mapper_profile_print_archive_search_pe_cuda(stream,map_output);
  // Approximate Search
  mapper_profile_print_approximate_search_summary(stream,true,map_output,num_threads);
}
#else /* GEM_PROFILE DISABLED */
/*
 * Mapper-CUDA SE
 */
void mapper_profile_print_mapper_se_cuda(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {}
/*
 * Mapper-CUDA PE
 */
void mapper_profile_print_mapper_pe_cuda(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {}
#endif
