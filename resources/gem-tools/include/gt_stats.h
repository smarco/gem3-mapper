/*
 * PROJECT: GEM-Tools library
 * FILE: gt_stats.h
 * DATE: 10/12/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_STATS_H_
#define GT_STATS_H_

#include "gt_commons.h"
#include "gt_compact_dna_string.h"
#include "gt_alignment_utils.h"
#include "gt_template_utils.h"

/*
 * Range Definition
 */
#define GT_STATS_MMAP_RANGE_0 0
#define GT_STATS_MMAP_RANGE_1 1
#define GT_STATS_MMAP_RANGE_5 2
#define GT_STATS_MMAP_RANGE_10 3
#define GT_STATS_MMAP_RANGE_50 4
#define GT_STATS_MMAP_RANGE_100 5
#define GT_STATS_MMAP_RANGE_500 6
#define GT_STATS_MMAP_RANGE_1000 7
#define GT_STATS_MMAP_RANGE_BEHOND 8
#define GT_STATS_MMAP_RANGE 9

#define GT_STATS_LENGTH_RANGE_5 0
#define GT_STATS_LENGTH_RANGE_40 1
#define GT_STATS_LENGTH_RANGE_80 2
#define GT_STATS_LENGTH_RANGE_100 3
#define GT_STATS_LENGTH_RANGE_150 4
#define GT_STATS_LENGTH_RANGE_300 5
#define GT_STATS_LENGTH_RANGE_800 6
#define GT_STATS_LENGTH_RANGE_1000 7
#define GT_STATS_LENGTH_RANGE_2000 8
#define GT_STATS_LENGTH_RANGE_5000 9
#define GT_STATS_LENGTH_RANGE_BEHOND 10
#define GT_STATS_LENGTH_RANGE 11

#define GT_STATS_LENGTH__MMAP_RANGE (GT_STATS_LENGTH_RANGE*GT_STATS_MMAP_RANGE)
#define GT_STATS_LENGTH__QUAL_SCORE_RANGE (GT_STATS_LENGTH_RANGE*GT_STATS_QUAL_SCORE_RANGE)
#define GT_STATS_QUAL_SCORE__MMAP_RANGE (GT_STATS_QUAL_SCORE_RANGE*GT_STATS_MMAP_RANGE)

#define GT_STATS_INSS_MIN  (-1000)
#define GT_STATS_INSS_MAX   50010
#define GT_STATS_INSS_STEP  10
#define GT_STATS_INSS_RANGE ((GT_STATS_INSS_MAX-GT_STATS_INSS_MIN+(GT_STATS_INSS_STEP-1))/GT_STATS_INSS_STEP)
#define GT_STATS_INSS_GET_BUCKET(INSS) (((INSS)-GT_STATS_INSS_MIN)/GT_STATS_INSS_STEP)

#define GT_STATS_MISMS_RANGE_0 0
#define GT_STATS_MISMS_RANGE_1 1
#define GT_STATS_MISMS_RANGE_2 2
#define GT_STATS_MISMS_RANGE_3 3
#define GT_STATS_MISMS_RANGE_4 4
#define GT_STATS_MISMS_RANGE_5 5
#define GT_STATS_MISMS_RANGE_6 6
#define GT_STATS_MISMS_RANGE_7 7
#define GT_STATS_MISMS_RANGE_8 8
#define GT_STATS_MISMS_RANGE_9 9
#define GT_STATS_MISMS_RANGE_10 10
#define GT_STATS_MISMS_RANGE_20 11
#define GT_STATS_MISMS_RANGE_50 12
#define GT_STATS_MISMS_RANGE_BEHOND 13
#define GT_STATS_MISMS_RANGE 14

#define GT_STATS_UNIQ_RANGE_0 0
#define GT_STATS_UNIQ_RANGE_1 1
#define GT_STATS_UNIQ_RANGE_2 2
#define GT_STATS_UNIQ_RANGE_3 3
#define GT_STATS_UNIQ_RANGE_10 4
#define GT_STATS_UNIQ_RANGE_50 5
#define GT_STATS_UNIQ_RANGE_100 6
#define GT_STATS_UNIQ_RANGE_500 7
#define GT_STATS_UNIQ_RANGE_BEHOND 8
#define GT_STATS_UNIQ_RANGE_X 9
#define GT_STATS_UNIQ_RANGE 10

#define GT_STATS_DELTA_EDIT_RANGE_0    0
#define GT_STATS_DELTA_EDIT_RANGE_1    1
#define GT_STATS_DELTA_EDIT_RANGE_2    2
#define GT_STATS_DELTA_EDIT_RANGE_3    3
#define GT_STATS_DELTA_EDIT_RANGE_4    4
#define GT_STATS_DELTA_EDIT_RANGE_5    5
#define GT_STATS_DELTA_EDIT_RANGE_6    6
#define GT_STATS_DELTA_EDIT_RANGE_7    7
#define GT_STATS_DELTA_EDIT_RANGE_8    8
#define GT_STATS_DELTA_EDIT_RANGE_9    9
#define GT_STATS_DELTA_EDIT_RANGE_10  10
#define GT_STATS_DELTA_EDIT_RANGE_BEHOND       11
#define GT_STATS_DELTA_EDIT_RANGE_INCONSISTENT 12
#define GT_STATS_DELTA_EDIT_RANGE_UNIQUE       13
#define GT_STATS_DELTA_EDIT_RANGE_UNMAPPED     14
#define GT_STATS_DELTA_EDIT_RANGE              15

#define GT_STATS_DELTA_SWG_RANGE_0    0
#define GT_STATS_DELTA_SWG_RANGE_1    1
#define GT_STATS_DELTA_SWG_RANGE_2    2
#define GT_STATS_DELTA_SWG_RANGE_3    3
#define GT_STATS_DELTA_SWG_RANGE_4    4
#define GT_STATS_DELTA_SWG_RANGE_5    5
#define GT_STATS_DELTA_SWG_RANGE_10   6
#define GT_STATS_DELTA_SWG_RANGE_20   7
#define GT_STATS_DELTA_SWG_RANGE_BEHOND       8
#define GT_STATS_DELTA_SWG_RANGE_INCONSISTENT 9
#define GT_STATS_DELTA_SWG_RANGE_UNIQUE       10
#define GT_STATS_DELTA_SWG_RANGE_UNMAPPED     11
#define GT_STATS_DELTA_SWG_RANGE              12

#define GT_STATS_CLASSES_TIE_PERFECT 0
#define GT_STATS_CLASSES_TIE         1
#define GT_STATS_CLASSES_MMAP_D1     2
#define GT_STATS_CLASSES_MMAP        3
#define GT_STATS_CLASSES_UNIQUE      4
#define GT_STATS_CLASSES_UNMAPPED    5
#define GT_STATS_CLASSES_RANGE       6

#define GT_STATS_MAPQ_BEHOND 61
#define GT_STATS_MAPQ_RANGE  62

#define GT_STATS_LARGE_READ_POS_RANGE 1000
#define GT_STATS_SHORT_READ_POS_RANGE 100

#define GT_STATS_MISMS_BASE_A 0
#define GT_STATS_MISMS_BASE_C 1
#define GT_STATS_MISMS_BASE_G 2
#define GT_STATS_MISMS_BASE_T 3
#define GT_STATS_MISMS_BASE_N 4
#define GT_STATS_MISMS_BASE_RANGE 5

#define GT_STATS_NT_BASE_A 0
#define GT_STATS_NT_BASE_C 1
#define GT_STATS_NT_BASE_G 2
#define GT_STATS_NT_BASE_T 3
#define GT_STATS_NT_BASE_N 4
#define GT_STATS_NT_BASE_RANGE 5

#define GT_STATS_QUAL_SCORE_RANGE 256

#define GT_STATS_NUM_JUNCTION_1      0
#define GT_STATS_NUM_JUNCTION_2      1
#define GT_STATS_NUM_JUNCTION_3      2
#define GT_STATS_NUM_JUNCTION_BEHOND 3
#define GT_STATS_NUM_JUNCTION_RANGE  4

#define GT_STATS_LEN_JUNCTION_100    0
#define GT_STATS_LEN_JUNCTION_1000   1
#define GT_STATS_LEN_JUNCTION_5000   2
#define GT_STATS_LEN_JUNCTION_10000  3
#define GT_STATS_LEN_JUNCTION_50000  4
#define GT_STATS_LEN_JUNCTION_BEHOND 5
#define GT_STATS_LEN_JUNCTION_RANGE  6

#define GT_STATS_MISMS_1_CONTEXT_RANGE ((GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE)*GT_STATS_MISMS_BASE_RANGE)
// DELTEME TODO #define GT_STATS_MISMS_2_CONTEXT_RANGE ((GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_BASE_RANGE)*GT_STATS_MISMS_BASE_RANGE)
#define GT_STATS_INDEL_TRANSITION_1_RANGE ((GT_STATS_MISMS_BASE_RANGE))  /* TODO */
#define GT_STATS_INDEL_TRANSITION_2_RANGE ((GT_STATS_MISMS_BASE_RANGE))  /* TODO */
#define GT_STATS_INDEL_TRANSITION_3_RANGE ((GT_STATS_MISMS_BASE_RANGE))  /* TODO */
#define GT_STATS_INDEL_TRANSITION_4_RANGE ((GT_STATS_MISMS_BASE_RANGE))  /* TODO */
#define GT_STATS_INDEL_1_CONTEXT           (GT_STATS_MISMS_BASE_RANGE)   /* TODO */
#define GT_STATS_INDEL_2_CONTEXT           (GT_STATS_MISMS_BASE_RANGE)   /* TODO */

#define GT_STATS_DIVERSITY_RANGE_0  0
#define GT_STATS_DIVERSITY_RANGE_1  1
#define GT_STATS_DIVERSITY_RANGE_2  2
#define GT_STATS_DIVERSITY_RANGE_3  3
#define GT_STATS_DIVERSITY_RANGE_4  4
#define GT_STATS_DIVERSITY_RANGE_5  5
#define GT_STATS_DIVERSITY_RANGE_10 6
#define GT_STATS_DIVERSITY_RANGE_20 7
#define GT_STATS_DIVERSITY_RANGE_50 8
#define GT_STATS_DIVERSITY_RANGE_BEHOND 9
#define GT_STATS_DIVERSITY_RANGE 10

#define GT_STATS_DOMINANT_RANGE_10  0
#define GT_STATS_DOMINANT_RANGE_25  1
#define GT_STATS_DOMINANT_RANGE_50  2
#define GT_STATS_DOMINANT_RANGE_75  3
#define GT_STATS_DOMINANT_RANGE_90  4
#define GT_STATS_DOMINANT_RANGE_95  5
#define GT_STATS_DOMINANT_RANGE_100 6
#define GT_STATS_DOMINANT_RANGE 7

#define GT_STATS_DIVERSITY_DOMINANT_RANGE (GT_STATS_DIVERSITY_RANGE*GT_STATS_DOMINANT_RANGE)

/*
 * Stats Data Structures
 */
typedef struct {
  // Mismatch/Indel Profile
  uint64_t *mismatches;        /* GT_STATS_MISMS_RANGE */
  uint64_t *levenshtein;       /* GT_STATS_MISMS_RANGE */
  uint64_t *insertion_length;  /* GT_STATS_MISMS_RANGE */
  uint64_t *deletion_length;   /* GT_STATS_MISMS_RANGE */
  uint64_t *errors_events;     /* GT_STATS_MISMS_RANGE */
  uint64_t total_mismatches;
  uint64_t total_levenshtein;
  uint64_t total_indel_length;
  uint64_t total_errors_events;
  uint64_t *error_position;    /* GT_STATS_LARGE_READ_POS_RANGE */
  // Trim/Mapping stats
  uint64_t total_bases;
  uint64_t total_bases_matching;
  uint64_t total_bases_trimmed;
  // Local/Global Alignment
  uint64_t total_local_maps;
  uint64_t total_global_maps;
  // Strandness combinations
  uint64_t single_strand_f;
  uint64_t single_strand_r;
  uint64_t pair_strand_rf;
  uint64_t pair_strand_fr;
  uint64_t pair_strand_ff;
  uint64_t pair_strand_rr;
  // Insert Size Distribution
  uint64_t *inss;               /* GT_STATS_INSS_RANGE */
  // Mismatch/Errors bases
  uint64_t *misms_transition;   /* GT_STATS_MISMS_BASE_RANGE*GT_STATS_MISMS_RANGE */
  uint64_t *qual_score_misms;   /* GT_STATS_QUAL_SCORE_RANGE */
  uint64_t *misms_1context;     /* GT_STATS_MISMS_1_CONTEXT_RANGE */
  uint64_t *indel_transition_1; /* GT_STATS_INDEL_TRANSITION_1_RANGE */
  uint64_t *indel_transition_2; /* GT_STATS_INDEL_TRANSITION_2_RANGE */
  uint64_t *indel_transition_3; /* GT_STATS_INDEL_TRANSITION_3_RANGE */
  uint64_t *indel_transition_4; /* GT_STATS_INDEL_TRANSITION_4_RANGE */
  uint64_t *indel_1context;     /* GT_STATS_INDEL_1_CONTEXT */
  uint64_t *indel_2context;     /* GT_STATS_INDEL_2_CONTEXT */
  uint64_t *qual_score_errors;  /* GT_STATS_QUAL_SCORE_RANGE */
} gt_maps_profile;

typedef struct {
  // General SM
  /*
   * @num_mapped_with_splitmaps, @num_mapped_only_splitmaps
   *   At the level of alignments... not matter if SE or PE or how many maps
   */
  uint64_t num_mapped_with_splitmaps;
  uint64_t num_mapped_only_splitmaps;
  /*
   * @total_splitmaps :: Number of blocks with SM.
   *   Eg 'chr1:+:12345:10>20*90::chr1:+:12445:100' => 1
   *   Eg 'chr1:+:12345:10>20*90::chr1:+:12445:10>20*90' => 2
   *   Eg 'chr1:+:12345:10>20*10>20*80::chr1:+:12445:10>20*90' => 2
   */
  uint64_t total_splitmaps;
  uint64_t total_junctions;
  uint64_t *num_junctions;     /* GT_STATS_NUM_JUNCTION_RANGE */
  uint64_t *length_junctions;  /* GT_STATS_LEN_JUNCTION_RANGE */
  uint64_t *junction_position; /* GT_STATS_SHORT_READ_POS_RANGE */
  // Paired SM combinations
  uint64_t pe_sm_sm;
  uint64_t pe_sm_rm;
  uint64_t pe_rm_rm;
} gt_splitmaps_profile;

typedef struct {
  // Diversity
  uint64_t* local_diversity;           /* GT_STATS_DIVERSITY_RANGE */
  uint64_t* local_dominant;            /* GT_STATS_DOMINANT_RANGE */
  uint64_t* local_diversity__dominant; /* GT_STATS_DIVERSITY_DOMINANT_RANGE */
  uint64_t global_diversity;
  // Quimeras
  uint64_t num_map_quimeras;
  uint64_t num_pair_quimeras;
  // Aux
  gt_shash* _local_diversity_hash;    // Auxiliary hash for diversity
  gt_shash* _global_diversity_hash;   // Auxiliary hash for diversity
} gt_population_profile;

typedef struct {
  // Length stats
  uint64_t min_length;
  uint64_t max_length;
  uint64_t total_bases; // WRT to the read
  uint64_t total_bases_aligned; // WRT to reads mapped
  uint64_t mapped_min_length;
  uint64_t mapped_max_length;
  uint64_t *length;            /* GT_STATS_LENGTH_RANGE */
  uint64_t *length_mapped;     /* GT_STATS_LENGTH_RANGE */
  uint64_t *length__mmap;      /* GT_STATS_LENGTH__MMAP_RANGE */
  uint64_t *length__quality;   /* GT_STATS_LENGTH__QUAL_SCORE_RANGE */
  uint64_t *avg_quality;       /* GT_STATS_QUAL_SCORE_RANGE */
  uint64_t *mmap__avg_quality; /* GT_STATS_QUAL_SCORE__MMAP_RANGE */
  // Nucleotide counting
  uint64_t *nt_counting;       /* GT_STATS_MISMS_BASE_RANGE */
  // Mapped/Maps
  uint64_t num_blocks;         // SE => 1 block. PE => 2 blocks
  uint64_t num_alignments;     // Number of alignments (independently of the type{SE,PE} or the syntax used in the file)
  uint64_t num_maps;
  uint64_t num_mapped;         // number of mapped reads for SE and number of mapped properly paired for templates
  // SE/PE read and mapping counts
  uint64_t num_mapped_reads;   // number of total mapped reads, counts 2 for a mapped paired template and checks unpaired pairs
  // Maps Distribution
  uint64_t *mmap;              /* GT_STATS_MMAP_RANGE */
  uint64_t *uniq;              /* GT_STATS_UNIQ_RANGE */
  uint64_t *delta_edit;        /* GT_STATS_DELTA_EDIT_RANGE */
  uint64_t *delta_swg;         /* GT_STATS_DELTA_SWG_RANGE */
  uint64_t *classes;           /* GT_STATS_CLASSES_RANGE */
  uint64_t *mapq;              /* GT_STATS_MAPQ_RANGE */
  // Error Profile
  gt_maps_profile *maps_profile;
  // Split maps info
  gt_splitmaps_profile* splitmaps_profile;
  // Population profile
  gt_population_profile* population_profile;
} gt_stats;

typedef struct {
  /* Analysis */
  bool first_map;
  bool nucleotide_stats;
  bool maps_profile;
  bool indel_profile; // TODO
  bool splitmap_profile;
  bool population_profile;
  /* Control Flags */
  bool use_map_counters; /* If possible, use counters instead of decoded matches */
} gt_stats_analysis;
#define GT_STATS_ANALYSIS_DEFAULT() \
  { \
    /* Analysis */ \
    .first_map=false, \
    .nucleotide_stats=true, \
    .maps_profile=true, \
    .indel_profile=false, \
    .splitmap_profile=true, \
    .population_profile=true, \
    /* Control Flags */ \
    .use_map_counters = true \
  }

/*
 * STATS Profile
 */
GT_INLINE gt_stats* gt_stats_new();
GT_INLINE void gt_stats_clear(gt_stats *stats);
GT_INLINE void gt_stats_delete(gt_stats *stats);

/*
 * STATS Merge
 */
void gt_stats_merge(gt_stats** const stats,const uint64_t stats_array_size);

/*
 * Calculate stats
 *   NOTE: @seq_archive==NULL if no indel_profile is requested (default)
 */
GT_INLINE void gt_stats_calculate_template_stats(
    gt_stats* const stats,gt_template* const template,gt_sequence_archive* seq_archive,gt_stats_analysis* const stats_analysis);

/*
 * STATS Report Output Printers
 */
GT_INLINE void gt_stats_print_mmap_distribution(FILE* stream,uint64_t* const mmap,const uint64_t num_alignments,const uint64_t num_mapped);
GT_INLINE void gt_stats_print_uniq_distribution(FILE* stream,uint64_t* const uniq,const uint64_t num_alignments);
GT_INLINE void gt_stats_print_inss_distribution(FILE* stream,uint64_t* const inss,const uint64_t num_maps,const bool full_resolution);
GT_INLINE void gt_stats_print_error_event_distribution(FILE* stream,uint64_t* const error,const uint64_t num_maps);
GT_INLINE void gt_stats_print_read_event_positions(FILE* stream,uint64_t* const pos_error,uint64_t const num_errors,uint64_t const max_length);
GT_INLINE void gt_stats_print_num_junctions_distribution(FILE* stream,uint64_t* const num_junctions,uint64_t const total);
GT_INLINE void gt_stats_print_length_junctions_distribution(FILE* stream,uint64_t* const length_junctions,uint64_t const total_junctions);
GT_INLINE void gt_stats_print_junction_position_distribution(FILE* stream,uint64_t* const junction_position,uint64_t const total_junctions,uint64_t const max_length);
GT_INLINE void gt_stats_print_qualities_error_distribution(FILE* stream,uint64_t* const qualities_error,uint64_t const total_error);
GT_INLINE void gt_stats_print_misms_transition_table(FILE* stream,uint64_t* const misms_trans,uint64_t const total_misms);
GT_INLINE void gt_stats_print_misms_transition_table_1context(FILE* stream,uint64_t* const misms_trans,uint64_t const total_misms);

GT_INLINE void gt_stats_print_population_stats(FILE* stream,gt_stats* const stats,const uint64_t num_reads,const bool paired_end);
GT_INLINE void gt_stats_print_split_maps_stats(FILE* stream,gt_stats* const stats,const bool paired_end);
GT_INLINE void gt_stats_print_maps_stats(FILE* stream, gt_stats* const stats,const uint64_t num_reads,const bool paired_end);
GT_INLINE void gt_stats_print_general_stats(FILE* stream,gt_stats* const stats,const uint64_t num_reads,const bool paired_end);

/*
 * Error Messages
 */
#define GT_ERROR_VSTATS_INVALID_MIN_MAX "Invalid step range for stats vector, min_value <= max_value"

#endif /* GT_STATS_H_ */
