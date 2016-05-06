/*
 * PROJECT: GEMMapper
 * FILE: mapper_profile.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "mapper/mapper_profile.h"
#include "mapper/mapper_profile_search.h"
#include "mapper/mapper_profile_misc.h"

#ifdef GEM_PROFILE /* GEM_PROFILE ENABLED */
/*
 * Archive Search
 */
void mapper_profile_print_archive_search_se(FILE* const stream) {
  fprintf(stream,    "[GEM]>Profile.ArchiveSearch.SE\n");
  tab_fprintf(stream,"  => TIME.Archive.Search.SE             ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_SE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Approximate.Search          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_AS_MAIN),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Select.Matches              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SELECT_SE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Score.Matches               ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SCORE_SE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
}
void mapper_profile_print_archive_search_pe(FILE* const stream) {
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
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_FIND_PAIRS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Archive.Select.PairedMatches        ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SELECT_PE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Score.PairedMatches                 ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SCORE_PE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    (B) Map.Extension\n");
  tab_fprintf(stream,"        => TIME.Archive.Extend.Candidates           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.Extend.Shortcut                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.Extend.Recovery                   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"          => TIME.Filtering.Extend.Candidates\n");
  tab_fprintf(stream,"            => TIME.Filtering.Extend.Match          ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_MATCH),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"              => TIME.Filtering.Retrieve.Candidates ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_RETRIEVE_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"              => TIME.Filtering.Verify.Candidates   ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_VERIFY_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"              => TIME.Filtering.Realign.Candidates  ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_FC_EXTEND_REALIGN_CANDIDATE_REGIONS),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    (C) Map.Selection\n");
  tab_fprintf(stream,"        => TIME.Select.Matches                      ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SELECT_PE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"        => TIME.Score.Matches                       ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SCORE_PE_MATCHES),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"  |> Archive.Extend.Candidates\n");
  tab_fprintf(stream,"    --> Extensions                                  ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_TOTAL),
      PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_TOTAL),"reads  ",true);
  tab_fprintf(stream,"      --> Extension.Matches                         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_NUM_MATCHES),NULL,"matches",true);
  tab_fprintf(stream,"      --> Extension.Length                          ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_EXTEND_VERIFY_CANDIDATES_LENGTH),NULL,"nt     ",true);
  tab_fprintf(stream,"      --> Extension.Regions.Found                   ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_FC_EXTEND_VERIFY_CANDIDATES_FOUND),NULL,"regions",true);
  tab_fprintf(stream,"    --> Extensions.Shortcut                         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT_TOTAL),
      PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_TOTAL),"reads  ",true);
  tab_fprintf(stream,"      --> Extensions.Shortcut.Success               ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_SHORTCUT_SUCCESS),
      PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_TOTAL),"reads  ",true);
  tab_fprintf(stream,"    --> Extensions.Recovery                         ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY_TOTAL),
      PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_TOTAL),"reads  ",true);
  tab_fprintf(stream,"      --> Extensions.Recovery.Success               ");
  COUNTER_PRINT(stream,PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTENSION_RECOVERY_SUCCESS),
      PROF_GET_COUNTER(GP_ARCHIVE_SEARCH_PE_EXTEND_CANDIDATES_TOTAL),"reads  ",true);
}
/*
 * Mapper SE
 */
void mapper_profile_print_mapper_se(
    FILE* const stream,const bool map_output,
    const uint64_t num_threads) {
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
  tab_fprintf(stream,"    ");
  if (map_output) {
    mapper_profile_print_map_output(stream,false);
  } else {
    mapper_profile_print_sam_output(stream,false);
  }
  // Archive Search SE
  mapper_profile_print_archive_search_se(stream);
  // Approximate Search
  mapper_profile_print_approximate_search_summary(stream,false,false,map_output,num_threads);
}
/*
 * Mapper PE
 */
void mapper_profile_print_mapper_pe(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {
  // All
  tab_fprintf(stream,"[GEM]>Profile.Mapper\n");
  tab_fprintf(stream,"  => TIME.Mapper                           ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_ALL),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Load.Index                     ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_MAPPER_LOAD_INDEX),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Parse.Input.FASTQ              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_INPUT_FASTA_PARSE_SEQUENCE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    => TIME.Archive.Search.PE              ");
  TIMER_PRINT(stream,PROF_GET_TIMER(GP_ARCHIVE_SEARCH_PE),PROF_GET_TIMER(GP_MAPPER_ALL));
  tab_fprintf(stream,"    ");
  if (map_output) {
    mapper_profile_print_map_output(stream,false);
  } else {
    mapper_profile_print_sam_output(stream,false);
  }
  // Archive Search SE
  mapper_profile_print_archive_search_pe(stream);
  // Approximate Search
  mapper_profile_print_approximate_search_summary(stream,true,false,map_output,num_threads);
}
#else /* GEM_PROFILE DISABLED */
void mapper_profile_print_mapper_se(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {}
void mapper_profile_print_mapper_pe(
    FILE* const stream,const bool map_output,const uint64_t num_threads) {}
#endif
