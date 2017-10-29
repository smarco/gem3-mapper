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
 *   GEM-profiling counters involved in the CUDA-mapping workflow
 */

#include "mapper/mapper_profile_cuda.h"
#include "mapper/mapper_profile_search.h"
#include "mapper/mapper_profile_misc.h"


#ifdef GEM_PROFILE /* GEM_PROFILE ENABLED */
/*
 * Stage Buffers
 */
void mapper_profile_print_search_stage_buffers(FILE* const stream) {
  tab_fprintf(stream,"[GEM]>Profile.CUDA.Search.Stages.Buffers\n");
  /*
   * Time
   */
  // Region-Profile
  tab_fprintf(stream,"  => (1) TIME.GPU.Buffer.Region.Profile\n");
  tab_fprintf(stream,"    => TIME.Alloc                            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_SEARCH_ALLOC),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Send                             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_SEARCH_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Receive                          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_SEARCH_RECEIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Duty.Cycle                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_SEARCH_DUTY_CYCLE),PROF_GET_TIMER(GP_MAPPER_ALL));
  // Decode
  tab_fprintf(stream,"  => (2) TIME.GPU.Buffer.Decode.Candidates\n");
  tab_fprintf(stream,"    => TIME.Alloc                            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_DECODE_ALLOC),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Send                             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_DECODE_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Receive                          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_DECODE_RECEIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Duty.Cycle                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_DECODE_DUTY_CYCLE),PROF_GET_TIMER(GP_MAPPER_ALL));
  // Kmer
  tab_fprintf(stream,"  => (3) TIME.GPU.Buffer.Kmer\n");
  tab_fprintf(stream,"    => TIME.Alloc                            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_KMER_FILTER_ALLOC),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Send                             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_KMER_FILTER_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Receive                          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_KMER_FILTER_RECEIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Duty.Cycle                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_KMER_FILTER_DUTY_CYCLE),PROF_GET_TIMER(GP_MAPPER_ALL));
  // BPM-distance
  tab_fprintf(stream,"  => (4) TIME.GPU.Buffer.BPM-Distance\n");
  tab_fprintf(stream,"    => TIME.Alloc                            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_BPM_DISTANCE_ALLOC),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Send                             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_BPM_DISTANCE_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Receive                          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_BPM_DISTANCE_RECEIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Duty.Cycle                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_BPM_DISTANCE_DUTY_CYCLE),PROF_GET_TIMER(GP_MAPPER_ALL));
  // BPM-align
  tab_fprintf(stream,"  => (5) TIME.GPU.Buffer.BPM-Align\n");
  tab_fprintf(stream,"    => TIME.Alloc                            ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_BPM_ALIGN_ALLOC),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Send                             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_BPM_ALIGN_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Receive                          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_BPM_ALIGN_RECEIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"      => TIME.Duty.Cycle                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_BPM_ALIGN_DUTY_CYCLE),PROF_GET_TIMER(GP_MAPPER_ALL));
  /*
   * Occupation
   */
  // Region-Profile
  tab_fprintf(stream,"  |> CUDA.Search.Stages.Buffers\n");
  tab_fprintf(stream,"    --> CUDA.RegionProfile\n");
  tab_fprintf(stream,"      --> Buffers.Used                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SSTAGE_REGION_PROFILE_BUFFERS),NULL,"        ",true);
  tab_fprintf(stream,"      --> Buffer.Searches                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SSTAGE_REGION_PROFILE_BUFFER_SEARCHES),NULL,"        ",true);
  tab_fprintf(stream,"      --> Buffer.Queries                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_FMI_SEARCH_NUM_QUERIES),NULL,"        ",true);
  tab_fprintf(stream,"        --> Buffer.Queries.Usage             ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_FMI_SEARCH_USAGE_QUERIES),"        ");
  // Decode
  tab_fprintf(stream,"    --> CUDA.DecodeCandidates\n");
  tab_fprintf(stream,"      --> Buffers.Used                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SSTAGE_DECODE_BUFFERS),NULL,"buffers ",true);
  tab_fprintf(stream,"      --> Buffer.Searches                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SSTAGE_DECODE_BUFFER_SEARCHES),NULL,"searches",true);
  tab_fprintf(stream,"      --> Buffer.Queries                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_FMI_DECODE_NUM_QUERIES),NULL,"queries ",true);
  tab_fprintf(stream,"        --> Buffer.Queries.Usage             ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_FMI_DECODE_USAGE_CANDIDATES),"        ");
  // Kmer-filter
  tab_fprintf(stream,"    --> CUDA.Kmer.filter\n");
  tab_fprintf(stream,"      --> Buffers.Used                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SSTAGE_KMER_FILTER_BUFFERS),NULL,"buffers ",true);
  tab_fprintf(stream,"      --> Buffer.Searches                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SSTAGE_KMER_FILTER_BUFFER_SEARCHES),NULL,"searches",true);
  tab_fprintf(stream,"      --> Buffer.Queries                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_KMER_FILTER_NUM_QUERIES),NULL,"queries ",true);
  tab_fprintf(stream,"        --> Buffer.Queries.Usage             ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_KMER_FILTER_USAGE_QUERIES),"        ");
  tab_fprintf(stream,"      --> Buffer.Candidates\n");
  tab_fprintf(stream,"        --> Candidates.Query                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_KMER_FILTER_CANDIDATES_PER_QUERY),NULL,"cand    ",true);
  tab_fprintf(stream,"        --> Candidates.Length                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_KMER_FILTER_CANDIDATE_LENGTH),NULL,"nt      ",true);
  tab_fprintf(stream,"        --> Buffer.Candidates.Usage          ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_KMER_FILTER_USAGE_CANDIDATES),"        ");
  // BPM-distance
  tab_fprintf(stream,"    --> CUDA.BPM-Distance\n");
  tab_fprintf(stream,"      --> Buffers.Used                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SSTAGE_BPM_DISTANCE_BUFFERS),NULL,"buffers ",true);
  tab_fprintf(stream,"      --> Buffer.Searches                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SSTAGE_BPM_DISTANCE_BUFFER_SEARCHES),NULL,"searches",true);
  tab_fprintf(stream,"      --> Buffer.Queries                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_BPM_DISTANCE_NUM_QUERIES),NULL,"queries ",true);
  tab_fprintf(stream,"        --> Buffer.Queries.Usage             ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_BPM_DISTANCE_USAGE_QUERIES),"        ");
  tab_fprintf(stream,"        --> Buffer.Entries.Usage             ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_BPM_DISTANCE_USAGE_PEQ_ENTRIES),"        ");
  tab_fprintf(stream,"      --> Buffer.Candidates\n");
  tab_fprintf(stream,"        --> Candidates.Query                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_BPM_DISTANCE_CANDIDATES_PER_TILE),NULL,"cand    ",true);
  tab_fprintf(stream,"        --> Candidates.Length                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_BPM_DISTANCE_CANDIDATES_LENGTH),NULL,"nt      ",true);
  tab_fprintf(stream,"        --> Buffer.Candidates.Usage          ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_BPM_DISTANCE_USAGE_CANDIDATES),"        ");
  tab_fprintf(stream,"        --> BPM.Cells.Computed               ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_BPM_DISTANCE_CELLS),NULL,"cells   ",true);
  // BPM-align
  tab_fprintf(stream,"    --> CUDA.BPM-Align\n");
  tab_fprintf(stream,"      --> Buffers.Used                       ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SSTAGE_BPM_ALIGN_BUFFERS),NULL,"buffers ",true);
  tab_fprintf(stream,"      --> Buffer.Searches                    ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_SSTAGE_BPM_ALIGN_BUFFER_SEARCHES),NULL,"searches",true);
  tab_fprintf(stream,"      --> Buffer.Queries                     ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_BPM_ALIGN_NUM_QUERIES),NULL,"queries ",true);
  tab_fprintf(stream,"        --> Buffer.Queries.Usage             ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_BPM_ALIGN_USAGE_QUERIES),"        ");
  tab_fprintf(stream,"      --> Buffer.Candidates\n");
  tab_fprintf(stream,"        --> Candidates.Query                 ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_BPM_ALIGN_CANDIDATE_PER_TILE),NULL,"cand    ",true);
  tab_fprintf(stream,"        --> Candidates.Length                ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_BPM_ALIGN_CANDIDATE_LENGTH),NULL,"nt      ",true);
  tab_fprintf(stream,"        --> Buffer.Candidates.Usage          ");
  PERCENTAGE_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_BPM_ALIGN_USAGE_CANDIDATES),"        ");
  tab_fprintf(stream,"        --> BPM.Cells.Computed               ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_GPU_BUFFER_BPM_ALIGN_CELLS),NULL,"cells   ",true);
}
/*
 * Archive Search SE CUDA
 */
void mapper_profile_print_archive_search_se_cuda(FILE* const stream,const bool map_output) {
//  // CUDA SE Search
//  tab_fprintf(stream,"[GEM]>Profile.CUDA.Search\n");
//  tab_fprintf(stream,"  => TIME.CUDA.Init                                        ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_COLLECTION_INIT),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"  => TIME.CUDA.Search                                      ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"    => (1)\n");
//  tab_fprintf(stream,"    => TIME.CUDA.STAGE.Region.Profile                      ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_SE_REGION_PROFILE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Buffer.Input                                 ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Parse.Input.FASTQ                            ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FASTA_PARSE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Init                                 ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_INIT),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Region.Profile.Generate              ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_GENERATE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Region.Profile.Copy                  ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_COPY),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.CUDA.Region.Profile.Send                     ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_SEARCH_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"    => (2)\n");
//  tab_fprintf(stream,"    => TIME.CUDA.STAGE.Decode.Candidates                   ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_SE_DECODE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.CUDA.Region.Profile.Receive                  ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_SEARCH_RECEIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Region.Profile.Retrieve              ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_REGION_PROFILE_RETRIEVE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"        => TIME.Archive.Region.Profile.CPU.Adaptive        ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ASSW_REGION_PROFILE_UNSUCCESSFUL),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Decode.Candidates.Copy               ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_DECODE_CANDIDATES_COPY),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.CUDA.Decode.Candidates.Send                  ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_DECODE_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"    => (3)\n");
//  tab_fprintf(stream,"    => TIME.CUDA.Verify.Candidates                         ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_SE_VERIFY_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.CUDA.Decode.Candidates.Receive               ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_FMI_DECODE_RECEIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Decode.Candidates.Retrieve           ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_DECODE_CANDIDATES_RETRIEVE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"        => TIME.Decode.Candidates.Buffered                 ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_DECODE_CANDIDATES_BUFFERED),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"          => TIME.Decode.Candidates.Buffered.Recomputed    ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"        => TIME.Compose.Candidates.Buffered                ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_COMPOSE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Verify.Candidates.Copy               ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES_COPY),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.CUDA.Verify.Candidates.Send                  ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_ALIGN_BPM_SEND),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"    => (4)\n");
//  tab_fprintf(stream,"    => TIME.CUDA.FinishSearch                              ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_SE_FINISH_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.CUDA.Verify.Candidates.Receive               ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_ALIGN_BPM_RECEIVE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Verify.Candidates.Retrieve           ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_VERIFY_CANDIDATES_RETRIEVE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"        => TIME.Filtering.Verify.Candidates.Retrieve       ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_VERIFY_CANDIDATES_BUFFERED),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Finish.Search                        ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE_FINISH_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Select.Matches                       ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SELECT_SE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Score.Matches                        ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SCORE_SE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      ");
//  if (map_output) {
//    mapper_profile_print_map_output(stream,false);
//  } else {
//    mapper_profile_print_sam_output(stream,false);
//  }
//  tab_fprintf(stream,"  |> (A) CUDA.Region.Profile\n");
//  tab_fprintf(stream,"  |> (B) CUDA.Decode.Candidates\n");
//  tab_fprintf(stream,"    --> Decode.Position.Copied                       ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ASSW_DECODE_CANDIDATES_COPIED),NULL,"positions",true);
//  tab_fprintf(stream,"    --> Decode.Position.Retrieved                    ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ASSW_DECODE_CANDIDATES_RETRIVED),NULL,"positions",true);
//  tab_fprintf(stream,"    --> Decode.Candidates.Buffered                   ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_DECODE_POSITIONS),
//                       PROF_GET_COUNTER(GP_FC_DECODE_POSITIONS),"positions",true);
//  tab_fprintf(stream,"    --> Decode.Candidates.Buffered.Recomputed        ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL_TOTAL),
//                       PROF_GET_COUNTER(GP_FC_DECODE_POSITIONS),"positions",true);
//  tab_fprintf(stream,"  |> (C) CUDA.Verify.Candidates\n");
//  tab_fprintf(stream,"    --> Verify.Candidates.Copied                     ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ASSW_VERIFY_CANDIDATES_TILES_COPIED),NULL,"tiles    ",true);
//  tab_fprintf(stream,"    --> Verify.Candidates.Retrieved                  ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ASSW_VERIFY_CANDIDATES_TILES_RETRIVED),NULL,"tiles    ",true);
//  tab_fprintf(stream,"    --> Verify.Tile.Distance.Diff                    ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_VERIFY_CANDIDATES_BUFFERED_DDIFF),NULL,"distance ",true);
//  // Search Stage Buffers
//  mapper_profile_print_search_stage_buffers(stream);
}
/*
 * Archive Search PE CUDA
 */
void mapper_profile_print_archive_search_pe_cuda(FILE* const stream,const bool map_output) {
//  // CUDA SE Search
//  tab_fprintf(stream,"[GEM]>Profile.CUDA.Search\n");
//  tab_fprintf(stream,"  => TIME.GPU.Init                                   ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_GPU_BUFFER_COLLECTION_INIT),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"  => TIME.CUDA.Search                                ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_PE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"    => (1)\n");
//  tab_fprintf(stream,"    => TIME.CUDA.STAGE.Region.Profile                ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_PE_REGION_PROFILE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Buffer.Input                           ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_BUFFERED_INPUT_RELOAD__DUMP_ATTACHED),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Parse.Input.FASTQ                      ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FASTA_PARSE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Init                           ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_INIT),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Region.Profile.Generate        ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_GENERATE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Region.Profile.Copy            ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_COPY),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"    => (2)\n");
//  tab_fprintf(stream,"    => TIME.CUDA.STAGE.Decode.Candidates             ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_PE_DECODE_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Region.Profile.Retrieve        ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_REGION_PROFILE_RETRIEVE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"        => TIME.Archive.Region.Profile.Unsuccessful  ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ASSW_REGION_PROFILE_UNSUCCESSFUL),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Decode.Candidates.Copy         ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_DECODE_CANDIDATES_COPY),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"    => (3)\n");
//  tab_fprintf(stream,"    => TIME.CUDA.Verify.Candidates                   ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_PE_VERIFY_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Decode.Candidates.Retrieve     ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_DECODE_CANDIDATES_RETRIEVE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"        => TIME.Decode.Candidates.Buffered           ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_DECODE_CANDIDATES_BUFFERED),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"          => TIME.Decode.Candidates.Buffered.Fail    ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"          --> Decode.Candidates.Buffered             ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_DECODE_POSITIONS),
//                       PROF_GET_COUNTER(GP_FC_DECODE_POSITIONS),"positions",true);
//  tab_fprintf(stream,"            --> Decode.Candidates.Buffered.Fail      ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_DECODE_CANDIDATES_BUFFERED_UNSUCCESSFUL_TOTAL),
//                       PROF_GET_COUNTER(GP_FC_DECODE_POSITIONS),"positions",true);
//  tab_fprintf(stream,"      => TIME.Archive.Verify.Candidates.Copy         ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_VERIFY_CANDIDATES_COPY),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"    => (4)\n");
//  tab_fprintf(stream,"    => TIME.CUDA.FinishSearch                        ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_CUDA_PE_FINISH_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Verify.Candidates.Retrieve     ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_VERIFY_CANDIDATES_RETRIEVE),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"        => TIME.Retrieve.Verify.Candidates           ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_VERIFY_CANDIDATES_BUFFERED),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"            --> Retrieve.Verify.Tile.Distance.Dif    ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_VERIFY_CANDIDATES_BUFFERED_DDIFF),NULL,"distance ",true);
//  tab_fprintf(stream,"      => TIME.Archive.Finish.Search                  ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_FINISH_SEARCH),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"        => TIME.Archive.Find.Paired.Matches          ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_FIND_PAIRS),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"        => TIME.Archive.Extend.Candidates            ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"          => TIME.Extend.Recovery                    ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"          => TIME.Filtering.Extend.Candidates\n");
//  tab_fprintf(stream,"            => TIME.Filtering.Extend.Match           ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_MATCH),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"              => TIME.Filtering.Retrieve.Candidates  ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_RETRIEVE_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"              => TIME.Filtering.Verify.Candidates    ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"              => TIME.Filtering.Realign.Candidates   ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"          --> Candidate.Regions.Length               ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_EXTEND_VERIFY_CANDIDATES_LENGTH),NULL,"nt",true);
//  tab_fprintf(stream,"      => TIME.Archive.Select.Matches                 ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SELECT_PE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"      => TIME.Archive.Score.Matches                  ");
//  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SCORE_PE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
//  tab_fprintf(stream,"  |> Archive.Extend.Candidates\n");
//  tab_fprintf(stream,"    --> Extensions                                   ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_TOTAL),
//      PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_TOTAL),"reads  ",true);
//  tab_fprintf(stream,"      --> Extension.Matches                          ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_NUM_MATCHES),NULL,"matches",true);
//  tab_fprintf(stream,"      --> Extension.Length                           ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_EXTEND_VERIFY_CANDIDATES_LENGTH),NULL,"nt     ",true);
//  tab_fprintf(stream,"      --> Extension.Regions.Found                    ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_EXTEND_VERIFY_CANDIDATES_FOUND),NULL,"regions",true);
//  tab_fprintf(stream,"    --> Extensions.Shortcut                          ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT_TOTAL),
//      PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_TOTAL),"reads  ",true);
//  tab_fprintf(stream,"      --> Extensions.Shortcut.Success                ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT_SUCCESS),
//      PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_TOTAL),"reads  ",true);
//  tab_fprintf(stream,"    --> Extensions.Recovery                          ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY_TOTAL),
//      PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_TOTAL),"reads  ",true);
//  tab_fprintf(stream,"      --> Extensions.Recovery.Success                ");
//  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY_SUCCESS),
//      PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_TOTAL),"reads  ",true);
//  tab_fprintf(stream,"      ");
//  if (map_output) {
//    mapper_profile_print_map_output(stream,false);
//  } else {
//    mapper_profile_print_sam_output(stream,false);
//  }
//  // Search Stage Buffers
//  mapper_profile_print_search_stage_buffers(stream);
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
  mapper_profile_print_approximate_search_summary(stream);
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
  mapper_profile_print_approximate_search_summary(stream);
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
