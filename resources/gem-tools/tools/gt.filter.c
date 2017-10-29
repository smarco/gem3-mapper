/*
 * PROJECT: GEM-Tools library
 * FILE: gt.filter.c
 * DATE: 02/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Application to filter {MAP,SAM,FASTQ} files and output the filtered result
 */
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "gem_tools.h"

#define GT_FILTER_FLOAT_NO_VALUE (-1.0)

#define gt_filter_cond_fatal_error_msg(condition,error_msg,args...) \
  gt_cond_fatal_error_msg(condition,error_msg ". File '%s', line %"PRIu64"\n",##args, \
      parameters.name_input_file,__buffered_input->current_line_num-1)
#define gt_filter_fatal_error_msg(error_msg,args...) \
  gt_fatal_error_msg(error_msg ". File '%s', line %"PRIu64"\n",##args, \
      parameters.name_input_file,__buffered_input->current_line_num-1)

typedef struct {
  uint64_t min;
  uint64_t max;
} gt_filter_quality_range;

typedef struct {
  /* I/O */
  char* name_input_file;
  bool mmap_input;
  gt_compression_t input_compression;
  bool paired_end;
  char* name_output_file;
  gt_file_format output_format;
  gt_compression_t output_compression;
  bool no_output;
  char* name_discarded_output_file;
  bool discarded_output;
  gt_file_format discarded_output_format;
  char* name_reference_file;
  char* name_gem_index_file;
  /* Filter Read/Qualities */
  bool hard_trim;
  uint64_t left_trim;
  uint64_t right_trim;
  bool restore_trim;
  bool uniform_read;
  bool uniform_read_strict;
  bool qualities_to_offset_33;
  bool qualities_to_offset_64;
  bool remove_qualities;
  bool add_qualities;
  /* Filter Template/Alignments */
  bool mapped;
  bool unmapped;
  int64_t unique_level;
  float min_length;
  float max_length;
  int64_t min_maps;
  int64_t max_maps;
  float max_strata_after_map;
  /* RNA Seq to recalculate counters */
  bool no_split_maps;
  bool only_split_maps;
  uint64_t min_intron_length;
  uint64_t min_block_length;
  /* Filter SE-Maps */
  bool first_map;
  bool keep_first_map;
  bool matches_pruning;
  uint64_t max_decoded_matches;
  uint64_t min_decoded_strata;
  uint64_t max_output_matches;
  uint64_t max_input_matches;
  bool make_counters;
  bool only_unmapped;
  bool only_mapped;
  float min_event_distance;
  float max_event_distance;
  float min_levenshtein_distance;
  float max_levenshtein_distance;
  gt_vector* map_ids;
  bool filter_by_strand_se;
  bool allow_strand_r;
  bool allow_strand_f;
  gt_vector* quality_score_ranges; /* (gt_filter_quality_range) */
  /* Filter PE-Maps */
  int64_t max_inss;
  int64_t min_inss;
  bool filter_by_strand_pe;
  bool allow_strand_rf;
  bool allow_strand_fr;
  bool allow_strand_ff;
  bool allow_strand_rr;
  /* Filter-Realign */
  bool mismatch_recovery;
  bool realign_hamming;
  bool realign_levenshtein;
  /* Checking/Report */
  bool check;
  bool check_format;
  gt_file_format check_file_format;
  /* Hidden */
  bool special_functionality;
  bool error_plot; // Print error distribution (depreciated)
  bool insert_size_plot; // Print insert size distribution (depreciated)
  bool show_sequence_list; // Display sequence list in the GEMindex/.fa...
  bool display_pretty; // Display pretty printed map(s)
  bool group_reads; // Group previously split reads
  bool sample_read; // Sample the read in chunks (annotated by chunk group)
  float split_chunk_size;
  float split_step_size;
  float split_left_trim;
  float split_right_trim;
  float split_min_remainder;
  /* Misc */
  uint64_t num_threads;
  bool verbose;
  /* Control flags */
  bool perform_dna_map_filter; // Any DNA-filtering criteria activated
  bool perform_rna_map_filter; // Any RNA-filtering criteria activated
  bool load_index;
} gt_filter_args;

gt_filter_args parameters = {
    /* I/O */
    .name_input_file=NULL,
    .mmap_input=false,
    .input_compression=GT_COMPRESSION_NONE,
    .paired_end=false,
    .name_output_file=NULL,
    .output_format=FILE_FORMAT_UNKNOWN,
    .output_compression=GT_COMPRESSION_NONE,
    .no_output=false,
    .name_discarded_output_file=NULL,
    .discarded_output = false,
    .discarded_output_format=FILE_FORMAT_UNKNOWN,
    .name_reference_file=NULL,
    .name_gem_index_file=NULL,
    /* Filter Read/Qualities */
    .hard_trim=false,
    .left_trim=0,
    .right_trim=0,
    .restore_trim=false,
    .uniform_read=false,
    .uniform_read_strict=false,
    .qualities_to_offset_33=false,
    .qualities_to_offset_64=false,
    .remove_qualities=false,
    .add_qualities=false,
    /* Filter Template/Alignments */
    .mapped=false,
    .unmapped=false,
    .unique_level=-1,
    .min_length=-1.0,
    .max_length=-1.0,
    .min_maps=-1,
    .max_strata_after_map=-1.0,
    .max_maps=-1,
    /* RNA Seq */
    .no_split_maps=false,
    .only_split_maps=false,
    .min_intron_length=0,
    .min_block_length=0,
    /* Filter SE-Maps */
    .first_map=false,
    .keep_first_map=false,
    .matches_pruning=false,
    .max_decoded_matches=GT_ALL,
    .min_decoded_strata=0,
    .max_output_matches=GT_ALL,
    .max_input_matches=GT_ALL,
    .make_counters=false,
    .only_unmapped=false,
    .only_mapped=false,
    .min_event_distance=GT_FILTER_FLOAT_NO_VALUE,
    .max_event_distance=GT_FILTER_FLOAT_NO_VALUE,
    .min_levenshtein_distance=GT_FILTER_FLOAT_NO_VALUE,
    .max_levenshtein_distance=GT_FILTER_FLOAT_NO_VALUE,
    .map_ids=NULL,
    .filter_by_strand_se=false,
    .allow_strand_r=false,
    .allow_strand_f=false,
    .quality_score_ranges = NULL,
    /* Filter PE-Maps */
    .max_inss=INT64_MAX,
    .min_inss=INT64_MIN,
    .filter_by_strand_pe=false,
    .allow_strand_rf=false,
    .allow_strand_fr=false,
    .allow_strand_ff=false,
    .allow_strand_rr=false,
    /* Filter-Realign */
    .mismatch_recovery=false,
    .realign_hamming=false,
    .realign_levenshtein=false,
    /* Checking/Report */
    .check = false,
    .check_format = false,
    /* Hidden */
    .special_functionality = false,
    .error_plot = false,
    .insert_size_plot = false,
    .show_sequence_list = false,
    .display_pretty = false,
    .group_reads = false,
    .sample_read = false,
    .split_chunk_size = -1.0,
    .split_step_size = -1.0,
    .split_left_trim = -1.0,
    .split_right_trim = -1.0,
    .split_min_remainder = 0.0,
    /* Misc */
    .num_threads=1,
    .verbose=false,
    /* Control flags */
    .perform_dna_map_filter=false,
    .perform_rna_map_filter=false,
    .load_index=false
};
/*
 * I/O Helpers
 */
gt_input_file* gt_filter_open_input_file(char* const name_input_file) {
  switch (parameters.input_compression) {
    case GT_COMPRESSION_NONE:
      return ((name_input_file==NULL) ?
          gt_input_stream_open(stdin) :
          gt_input_file_open(name_input_file,parameters.mmap_input));
      break;
    case GT_COMPRESSION_GZIP:
      return ((name_input_file==NULL) ?
          gt_input_gzip_stream_open(stdin) :
          gt_input_gzip_stream_open(gt_open_FILE(name_input_file,"r")));
      break;
    case GT_COMPRESSION_BZIP:
      return ((name_input_file==NULL) ?
          gt_input_bzip_stream_open(stdin) :
          gt_input_bzip_stream_open(gt_open_FILE(name_input_file,"r")));
      break;
    default:
      GT_INVALID_CASE();
      break;
  }
  return NULL;
}
gt_output_file* gt_filter_open_output_file(char* const name_output_file) {
  switch (parameters.output_compression) {
    case GT_COMPRESSION_NONE:
      return ((name_output_file==NULL) ?
          gt_output_stream_new(stdout,SORTED_FILE) :
          gt_output_file_new(name_output_file,SORTED_FILE));
      break;
    case GT_COMPRESSION_GZIP:
      return ((name_output_file==NULL) ?
          gt_output_gzip_stream_new(stdout,SORTED_FILE) :
          gt_output_gzip_stream_new(gt_open_FILE(name_output_file,"w"),SORTED_FILE));
      break;
    case GT_COMPRESSION_BZIP:
      return ((name_output_file==NULL) ?
          gt_output_bzip_stream_new(stdout,SORTED_FILE) :
          gt_output_bzip_stream_new(gt_open_FILE(name_output_file,"w"),SORTED_FILE));
      break;
    default:
      GT_INVALID_CASE();
      break;
  }
  return NULL;
}
/*
 * Helper to get num maps correctly also for unpaired 
 * mapped pairs
 */
GT_INLINE uint64_t gt_filter_get_num_maps(gt_template* template){
  GT_TEMPLATE_IF_SE_ALINGMENT(template) {
    return gt_template_get_num_mmaps(template);
  } else {
    if (!gt_template_is_mapped(template)) {
      GT_TEMPLATE_REDUCE_BOTH_ENDS(template,alignment_end1,alignment_end2);
      return gt_alignment_get_num_maps(alignment_end1) + gt_alignment_get_num_maps(alignment_end2);
    } else {
      return gt_template_get_num_mmaps(template);
    }
  }
}
/*
 * Checking/(Re)Aligning/MismsRecovery
 */
GT_INLINE void gt_filter_mismatch_recovery_maps(
    char* const name_input_file,const uint64_t current_line_num,
    gt_template* const template,gt_sequence_archive* const sequence_archive) {
  // Unfolded as to report errors in the recovery
  gt_status error_code;
  uint64_t alignment_pos = 0;
  GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
    uint64_t map_pos = 0;
    GT_ALIGNMENT_ITERATE(alignment,map) {
      if ((error_code=gt_map_recover_mismatches_sa(map,alignment->read,sequence_archive))) {
        gt_error_msg("Unrecoverable Alignment '%s':%"PRIu64"\n\tREAD::'"PRIgts"':%"PRIu64":%"PRIu64" ",
            name_input_file,current_line_num,PRIgts_content(template->tag),alignment_pos,map_pos);
        gt_output_map_fprint_map_pretty_sa(stderr,map,alignment->read,sequence_archive);
      }
      ++map_pos;
    }
    gt_alignment_recalculate_counters(alignment);
    ++alignment_pos;
  }
  if (gt_template_get_num_blocks(template)>1) gt_template_recalculate_counters(template);
}
GT_INLINE bool gt_filter_check_maps(
    char* const name_input_file,const uint64_t current_line_num,
    gt_template* const template,gt_sequence_archive* const sequence_archive,
    uint64_t* const total_algs_checked,uint64_t* const total_algs_correct,
    uint64_t* const total_maps_checked,uint64_t* const total_maps_correct) {
  bool alignment_correct=true;
  gt_status error_code;
  uint64_t alignment_pos = 0;
  GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
    uint64_t map_pos = 0;
    GT_ALIGNMENT_ITERATE(alignment,map) {
      if ((error_code=gt_map_check_alignment_sa(map,alignment->read,sequence_archive))) {
        gt_error_msg("Wrong Alignment '%s':%"PRIu64"\n\tREAD::'"PRIgts"':%"PRIu64":%"PRIu64" ",
            name_input_file,current_line_num,PRIgts_content(template->tag),alignment_pos,map_pos);
        gt_output_map_fprint_map_pretty_sa(stderr,map,alignment->read,sequence_archive);
        alignment_correct = false;
      } else {
        ++(*total_maps_correct);
      }
      ++(*total_maps_checked);
      ++map_pos;
    }
    ++alignment_pos;
  }
  ++(*total_algs_checked);
  if (alignment_correct) {
    ++(*total_algs_correct);
    return true;
  } else {
    return false;
  }
}
/*
 * Filtering MAPs functions
 */
void gt_filter_delete_map_ids(gt_vector* filter_map_ids) {
  // Free vector
  if (filter_map_ids!=NULL) {
    GT_VECTOR_ITERATE(filter_map_ids,map_id,pos,gt_string*) {
      gt_string_delete(*map_id);
    }
    gt_vector_delete(filter_map_ids);
  }
}
GT_INLINE bool gt_filter_is_sequence_name_allowed(gt_string* const seq_name) {
  GT_VECTOR_ITERATE(parameters.map_ids,map_id,pos,gt_string*) {
    if (gt_string_equals(seq_name,*map_id)) return true;
  }
  return false;
}
GT_INLINE bool gt_filter_is_quality_value_allowed(const uint64_t quality_score) {
  GT_VECTOR_ITERATE(parameters.quality_score_ranges,quality_range,pos,gt_filter_quality_range) {
    if (quality_score >= quality_range->min && quality_score <= quality_range->max) return true;
  }
  return false;
}
GT_INLINE void gt_filter_prune_matches(gt_template* const template) {
  uint64_t max_num_matches = GT_ALL;
  if (parameters.max_decoded_matches!=GT_ALL || parameters.min_decoded_strata!=0) {
    uint64_t max_strata;
    gt_counters_calculate_num_maps(gt_template_get_counters_vector(template),
        parameters.min_decoded_strata,parameters.max_decoded_matches,&max_strata,&max_num_matches);
  }
  if (parameters.max_output_matches!=GT_ALL) {
    max_num_matches = GT_MIN(max_num_matches,parameters.max_output_matches);
  }
  // Reduce matches
  if (max_num_matches < GT_ALL) {
    gt_template_reduce_mmaps(template,max_num_matches);
  }
}
void gt_alignment_dna_filter(gt_alignment* const alignment_dst,gt_alignment* const alignment_src,const gt_file_format file_format) {
  const uint64_t first_matching_distance = gt_counters_get_min_matching_strata(gt_alignment_get_counters_vector(alignment_src)) - 1;
  // Reduction by unique level (can be calculated beforehand)
  bool pick_only_first_map = false;
  /*
   * (1) Pre-filtering steps
   */
  gt_map* first_map = NULL;
  if (parameters.keep_first_map && gt_alignment_get_num_maps(alignment_src)>0) {
    first_map = gt_map_copy(gt_alignment_get_map(alignment_src,0));
  }
  /*
   * (2) Filtering of maps
   */
  GT_ALIGNMENT_ITERATE(alignment_src,map) {
    // Check sequence name
    if (parameters.map_ids!=NULL) {
      if (!gt_filter_is_sequence_name_allowed(map->seq_name)) continue;
    }
    // Filter strata beyond first mapping
    const int64_t current_stratum = gt_map_get_global_distance(map);
    if (parameters.max_strata_after_map >= 0.0 &&
        (current_stratum-first_matching_distance) > gt_alignment_get_read_proportion(alignment_src,parameters.max_strata_after_map)) break;
    // Check strata
    if (parameters.min_event_distance != GT_FILTER_FLOAT_NO_VALUE || parameters.max_event_distance != GT_FILTER_FLOAT_NO_VALUE) {
      const uint64_t total_distance = gt_map_get_global_distance(map);
      if (parameters.min_event_distance != GT_FILTER_FLOAT_NO_VALUE) {
        if (total_distance < gt_alignment_get_read_proportion(alignment_src,parameters.min_event_distance)) continue;
      }
      if (parameters.max_event_distance != GT_FILTER_FLOAT_NO_VALUE) {
        if (total_distance > gt_alignment_get_read_proportion(alignment_src,parameters.max_event_distance)) continue;
      }
    }
    // Check levenshtein distance
    if (parameters.min_levenshtein_distance != GT_FILTER_FLOAT_NO_VALUE || parameters.max_levenshtein_distance != GT_FILTER_FLOAT_NO_VALUE) {
      const uint64_t total_distance = gt_map_get_global_levenshtein_distance(map);
      if (parameters.min_levenshtein_distance != GT_FILTER_FLOAT_NO_VALUE) {
        if (total_distance < gt_alignment_get_read_proportion(alignment_src,parameters.min_levenshtein_distance)) continue;
      }
      if (parameters.max_levenshtein_distance != GT_FILTER_FLOAT_NO_VALUE) {
        if (total_distance > gt_alignment_get_read_proportion(alignment_src,parameters.max_levenshtein_distance)) continue;
      }
    }
    // Filter strand
    if (parameters.filter_by_strand_se) {
      if (map->strand==FORWARD && !parameters.allow_strand_f) continue;
      if (map->strand==REVERSE && !parameters.allow_strand_r) continue;
    }
    // Filter quality scores
    if (parameters.quality_score_ranges!=NULL) {
      if (!gt_filter_is_quality_value_allowed((file_format==SAM) ? map->phred_score : map->gt_score)) continue;
    }
    /*
     * Insert the map
     */
    gt_alignment_insert_map(alignment_dst,gt_map_copy(map));
    // Skip the rest if first map is enabled
    if (parameters.first_map || pick_only_first_map) break;
  }
  /*
   * (3) Post-filtering steps
   */
  if (parameters.keep_first_map) {
    if (gt_alignment_get_num_maps(alignment_dst)==0) {
      gt_alignment_insert_map(alignment_dst,first_map);
    } else {
      gt_map_delete(first_map);
    }
  }
}
void gt_template_dna_filter(gt_template* const template_dst,gt_template* const template_src,const gt_file_format file_format) {
  /*
   * Filtering workflow
   *   (1) Pre-filtering steps
   *   (2) Filtering of maps (taking them into account individually)
   *   (3) Reduction of all maps (taking them into account as a whole)
   *   (4) Post-filtering steps
   */
  GT_TEMPLATE_IF_SE_ALINGMENT(template_src) {
    GT_TEMPLATE_REDUCTION(template_src,alignment_src);
    GT_TEMPLATE_REDUCTION(template_dst,alignment_dst);
    gt_alignment_dna_filter(alignment_dst,alignment_src,file_format);
  } else {
    if (!gt_template_is_mapped(template_src)) {
      GT_TEMPLATE_REDUCE_BOTH_ENDS(template_src,alignment_src_end1,alignment_src_end2);
      GT_TEMPLATE_REDUCE_BOTH_ENDS(template_dst,alignment_dst_end1,alignment_dst_end2);
      gt_alignment_dna_filter(alignment_dst_end1,alignment_src_end1,file_format);
      gt_alignment_dna_filter(alignment_dst_end2,alignment_src_end2,file_format);
    } else {
      const uint64_t first_matching_distance = gt_counters_get_min_matching_strata(gt_template_get_counters_vector(template_src))-1;
      // Reduction by unique level (can be calculated beforehand)
      bool pick_only_first_map = false;
      /*
       * (1) Pre-filtering steps
       */
      gt_map** first_mmap = NULL;
      gt_mmap_attributes first_mmap_attributes = {0, 0, 0};
      if (parameters.keep_first_map && gt_template_get_num_mmaps(template_src)>0) {
        gt_mmap* const mmap = gt_template_get_mmap(template_src,0);
        first_mmap = gt_mmap_array_copy(mmap->mmap,gt_template_get_num_blocks(template_src));
        first_mmap_attributes = mmap->attributes;
      }
      /*
       * (2) Filtering of maps
       */
      GT_TEMPLATE_ITERATE_MMAP__ATTR(template_src,mmap,mmap_attributes) {
        const int64_t current_stratum = gt_map_get_global_distance(mmap[0]) + gt_map_get_global_distance(mmap[1]);
        if (parameters.max_strata_after_map >= 0.0 &&
            (current_stratum-first_matching_distance) > gt_template_get_read_proportion(template_src,parameters.max_strata_after_map)) break;
        // Check sequence name
        if (parameters.map_ids!=NULL) {
          if (!gt_filter_is_sequence_name_allowed(mmap[0]->seq_name)) continue;
          if (!gt_filter_is_sequence_name_allowed(mmap[1]->seq_name)) continue;
        }
        // Check strata
        if (parameters.min_event_distance != GT_FILTER_FLOAT_NO_VALUE || parameters.max_event_distance != GT_FILTER_FLOAT_NO_VALUE) {
          const int64_t total_distance = gt_map_get_global_distance(mmap[0]) + gt_map_get_global_distance(mmap[1]);
          if (parameters.min_event_distance != GT_FILTER_FLOAT_NO_VALUE) {
            if (total_distance < gt_template_get_read_proportion(template_src,parameters.min_event_distance)) continue;
          }
          if (parameters.max_event_distance != GT_FILTER_FLOAT_NO_VALUE) {
            if (total_distance > gt_template_get_read_proportion(template_src,parameters.max_event_distance)) continue;
          }
        }
        // Check levenshtein distance
        if (parameters.min_levenshtein_distance != GT_FILTER_FLOAT_NO_VALUE || parameters.max_levenshtein_distance != GT_FILTER_FLOAT_NO_VALUE) {
          const int64_t total_distance = gt_map_get_global_levenshtein_distance(mmap[0])+gt_map_get_global_levenshtein_distance(mmap[1]);
          if (parameters.min_levenshtein_distance != GT_FILTER_FLOAT_NO_VALUE) {
            if (total_distance < gt_template_get_read_proportion(template_src,parameters.min_levenshtein_distance)) continue;
          }
          if (parameters.max_levenshtein_distance != GT_FILTER_FLOAT_NO_VALUE) {
            if (total_distance > gt_template_get_read_proportion(template_src,parameters.max_levenshtein_distance)) continue;
          }
        }
        // Check inss
        if (parameters.min_inss > INT64_MIN || parameters.max_inss < INT64_MAX) {
          gt_status error_code;
          const int64_t inss = gt_template_get_insert_size(mmap,&error_code,0,0);
          if (parameters.min_inss > inss || inss > parameters.max_inss) continue;
        }
        // Check strandness
        if (parameters.filter_by_strand_se) {
          if (!parameters.allow_strand_f && (mmap[0]->strand==FORWARD || mmap[1]->strand==FORWARD)) continue;
          if (!parameters.allow_strand_r && (mmap[0]->strand==REVERSE || mmap[1]->strand==REVERSE)) continue;
        }
        if (parameters.filter_by_strand_pe) {
          if (mmap[0]->strand==FORWARD && mmap[1]->strand==REVERSE && !parameters.allow_strand_fr) continue;
          if (mmap[0]->strand==REVERSE && mmap[1]->strand==FORWARD && !parameters.allow_strand_rf) continue;
          if (mmap[0]->strand==FORWARD && mmap[1]->strand==FORWARD && !parameters.allow_strand_ff) continue;
          if (mmap[0]->strand==REVERSE && mmap[1]->strand==REVERSE && !parameters.allow_strand_rr) continue;
        }
        // Filter quality scores
        if (parameters.quality_score_ranges!=NULL) {
          if (!gt_filter_is_quality_value_allowed((file_format==SAM) ? mmap_attributes->phred_score : mmap_attributes->gt_score)) continue;
        }
        /*
         * Insert the map
         */
        gt_map** mmap_copy = gt_mmap_array_copy(mmap,__mmap_num_blocks);
        gt_template_insert_mmap(template_dst,mmap_copy,mmap_attributes);
        free(mmap_copy);
        // Skip the rest if first map is enabled
        if (parameters.first_map || pick_only_first_map) break;
      }
      /*
       * (3) Post-filtering steps
       */
      if (parameters.keep_first_map) {
        if (gt_template_get_num_mmaps(template_dst)==0) {
          gt_template_insert_mmap(template_dst,first_mmap,&first_mmap_attributes);
        }
        free(first_mmap);
      }
    }
  }
}
void gt_alignment_rna_filter(gt_alignment* const alignment_dst,gt_alignment* const alignment_src,const gt_file_format file_format) {
  GT_ALIGNMENT_ITERATE(alignment_src,map) {
    // Check sequence name
    if (parameters.map_ids!=NULL) {
      if (!gt_filter_is_sequence_name_allowed(map->seq_name)) continue;
    }
    // Check SM contained
    const uint64_t num_blocks = gt_map_get_num_blocks(map);
    if (parameters.no_split_maps && num_blocks>1) continue;
    if (parameters.only_split_maps && num_blocks==1) continue;
    // Filter intron length
    if (parameters.min_intron_length > 0) {
      if (gt_map_get_num_blocks(map) > 1) {
        if(gt_map_get_min_intron_length(map) < parameters.min_intron_length){
          continue;
        }
      }
    }
    // Filter block length
    if (parameters.min_block_length > 0) {
      if (gt_map_get_num_blocks(map) > 1) {
        if (gt_map_get_min_block_length(map) < parameters.min_block_length) continue;
      }
    }
    // Insert the map
    gt_alignment_insert_map(alignment_dst,gt_map_copy(map));
    // Skip the rest if best
    if (parameters.first_map) return;
  }
}
void gt_template_rna_filter(gt_template* const template_dst,gt_template* const template_src,const gt_file_format file_format) {
  GT_TEMPLATE_IF_SE_ALINGMENT(template_src) {
    GT_TEMPLATE_REDUCTION(template_src,alignment_src);
    GT_TEMPLATE_REDUCTION(template_dst,alignment_dst);
    /*
     * SE
     */
    gt_alignment_rna_filter(alignment_dst,alignment_src,file_format);
  } else {
    /*
     * PE
     */
    if (!gt_template_is_mapped(template_src)) {
      GT_TEMPLATE_REDUCE_BOTH_ENDS(template_src,alignment_src_end1,alignment_src_end2);
      GT_TEMPLATE_REDUCE_BOTH_ENDS(template_dst,alignment_dst_end1,alignment_dst_end2);
      gt_alignment_rna_filter(alignment_dst_end1,alignment_src_end1,file_format);
      gt_alignment_rna_filter(alignment_dst_end2,alignment_src_end2,file_format);
    } else {
      const uint64_t num_blocks = gt_template_get_num_blocks(template_src);
      GT_TEMPLATE_ITERATE_MMAP__ATTR(template_src,mmap,mmap_attributes) {
        // Check SM contained and get minimum intron length
        uint64_t has_sm = false;
        uint64_t min_intron_length = UINT64_MAX, min_block_length = UINT64_MAX;
        if (parameters.no_split_maps || parameters.only_split_maps || parameters.min_intron_length >= 0) {
          GT_MMAP_ITERATE(mmap,map,end_p) {
            if (gt_map_get_num_blocks(map) > 1) {
              const uint64_t mil = gt_map_get_min_intron_length(map);
              const uint64_t mbl = gt_map_get_min_block_length(map);
              has_sm = true;
              if (mil >= 0 && mil < min_intron_length) min_intron_length = mil;
              if (mbl >= 0 && mbl < min_block_length) min_block_length = mbl;
            }
          }
        }
        if (parameters.no_split_maps && has_sm) continue;
        if (parameters.only_split_maps && !has_sm) continue;
        // Filter intron length
        if (parameters.min_intron_length > 0 && min_intron_length != UINT64_MAX){
          if(min_intron_length < parameters.min_intron_length){
            continue;
          }
        }
        // Filter block length
        if (parameters.min_block_length > 0 && min_block_length != UINT64_MAX){
          if(min_block_length < parameters.min_block_length) continue;
        }
        // Add the mmap
        gt_map** mmap_copy = gt_mmap_array_copy(mmap,num_blocks);
        gt_template_insert_mmap(template_dst,mmap_copy,mmap_attributes);
        gt_free(mmap_copy);
        // Skip the rest if best
        if (parameters.first_map) return;
      }
    }
  }
}
GT_INLINE bool gt_filter_apply_filters(
    const gt_file_format file_format,const uint64_t line_no,
    gt_sequence_archive* const sequence_archive,gt_template* const template) {
  /*
   * Process Read/Qualities // TODO: move out of filter (this is processing)
   */
  const uint64_t has_qualities = gt_template_has_qualities(template);
  if (parameters.remove_qualities && has_qualities) {
    GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
      gt_string_clear(alignment->qualities);
    }
  } else if (parameters.add_qualities && !has_qualities) {
    GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
      const uint64_t read_length = gt_alignment_get_read_length(alignment);
      gt_string_resize(alignment->qualities,read_length+1);
      gt_string_set_length(alignment->qualities,read_length);
      GT_STRING_ITERATE(alignment->qualities,buffer,i) {
        buffer[i]='~';
      }
    }
  }
  if (parameters.uniform_read) {
    if (parameters.uniform_read_strict) {
      GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
        gt_dna_read_uniform_strict_content(alignment->read,alignment->qualities);
      }
    } else {
      GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
        gt_dna_read_uniform_content(alignment->read,alignment->qualities);
      }
    }
  }
  if (has_qualities) {
    if (parameters.qualities_to_offset_33) {
      GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
        gt_qualities_adapt_from_offset64_to_offset33(alignment->qualities);
      }
    }
    if (parameters.qualities_to_offset_64) {
      GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
        gt_qualities_adapt_from_offset33_to_offset64(alignment->qualities);
      }
    }
  }
  /*
   * Template/Alignment Filter
   */
  // Consider mapped/unmapped
  const bool is_mapped = gt_template_is_mapped(template);
  if (parameters.mapped && !is_mapped) return false;
  if (parameters.unmapped && is_mapped) return false;
  // Unique based filtering
  if (parameters.unique_level>=0.0 && is_mapped) {
    if (parameters.unique_level > gt_template_get_uniq_degree(template)) return false;
  }
  // Filter by read length
  if (parameters.min_length>=0.0 || parameters.max_length>=0.0) {
    GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
      const uint64_t read_length = gt_alignment_get_read_length(alignment);
      if (parameters.min_length>=0.0) {
        const uint64_t min_length = gt_alignment_get_read_proportion(alignment,parameters.min_length);
        if (read_length < min_length) return false;
      }
      if (parameters.max_length>=0.0) {
        const uint64_t max_length = gt_alignment_get_read_proportion(alignment,parameters.max_length);
        if (read_length > max_length) return false;
      }
    }
  }
  // Filter by number of maps
  if (parameters.min_maps>=0 || parameters.max_maps>=0) {
    const uint64_t num_maps = gt_template_get_num_mmaps(template);
    if (parameters.min_maps>=0 && num_maps<parameters.min_maps) return false;
    if (parameters.max_maps>=0 && num_maps>parameters.max_maps) return false;
  }
  /*
   * MAP Filter
   */
  // Trim
  if (parameters.hard_trim) {
    gt_template_hard_trim(template,parameters.left_trim,parameters.right_trim);
    gt_template_recalculate_counters(template);
  } else if (parameters.restore_trim) {
    gt_template_restore_trim(template);
    gt_template_recalculate_counters(template);
  }
  // (Re)Align
  if (parameters.realign_levenshtein) {
    gt_template_realign_levenshtein(template,sequence_archive);
  } else if (parameters.realign_hamming) {
    gt_template_realign_hamming(template,sequence_archive);
  } else if (parameters.mismatch_recovery) {
    gt_filter_mismatch_recovery_maps(parameters.name_input_file,line_no,template,sequence_archive);
  }
  // Map DNA-filtering
  if (parameters.perform_dna_map_filter) {
    gt_template *template_filtered = gt_template_dup(template,false,false);
    gt_template_dna_filter(template_filtered,template,file_format);
    gt_template_swap(template,template_filtered);
    gt_template_delete(template_filtered);
  }
  // Map RNA-filtering
  if (parameters.perform_rna_map_filter) {
    gt_template *template_filtered = gt_template_dup(template,false,false);
    gt_template_rna_filter(template_filtered,template,file_format);
    gt_template_swap(template,template_filtered);
    gt_template_delete(template_filtered);
  }
  // Map pruning
  if (parameters.matches_pruning) gt_filter_prune_matches(template);
  // Make counters
  if (parameters.make_counters) {
    gt_template_recalculate_counters(template);
  }
  // Ok, go on
  return true;
}
GT_INLINE void gt_filter__print(
    const gt_file_format file_format,const uint64_t line_no,
    gt_sequence_archive* const sequence_archive,gt_template* const template,
    uint64_t* const total_algs_checked,uint64_t* const total_algs_correct,
    uint64_t* const total_maps_checked,uint64_t* const total_maps_correct,
    gt_buffered_output_file* const buffered_output,gt_generic_printer_attributes* const generic_printer_attributes,
    gt_buffered_output_file* const buffered_discarded_output,gt_generic_printer_attributes* const discarded_output_attributes) {
  bool discaded = false;
  /*
   * Apply Filters
   */
  if (!gt_filter_apply_filters(file_format,line_no,sequence_archive,template)) discaded = true;
  if (parameters.uniform_read) { // Check zero-length reads
    GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
      if (gt_alignment_get_read_length(alignment)==0) return;
    }
  }
  /*
   * Check
   */
  if (!discaded && parameters.check) {
    if (!gt_filter_check_maps(parameters.name_input_file,line_no,
        template,sequence_archive,total_algs_checked,total_algs_correct,total_maps_checked,total_maps_correct)) discaded = true;
  }
  /*
   * Print template
   */
  if (!parameters.no_output && !discaded) {
    if (gt_output_generic_bofprint_template(buffered_output,template,generic_printer_attributes)) {
      gt_error_msg("Fatal error outputting read '"PRIgts"'(InputLine:%"PRIu64")\n",
          PRIgts_content(gt_template_get_string_tag(template)),line_no);
    }
  } else if (discaded && buffered_discarded_output!=NULL) {
    if (gt_output_generic_bofprint_template(buffered_discarded_output,template,discarded_output_attributes)) {
      gt_error_msg("Fatal error outputting read '"PRIgts"'(InputLine:%"PRIu64")\n",
          PRIgts_content(gt_template_get_string_tag(template)),line_no);
    }
  }
}
/*
 * Special funcionality
 */
GT_INLINE void gt_filter_sample_read_print_fastq(
    gt_buffered_output_file* const buffered_output,gt_string* const tag,gt_string* const read,gt_string* const qualities,
    const bool print_segmented_read_info,const uint64_t segment_id,const uint64_t total_segments,
    const uint64_t left_trim,const uint64_t right_trim,const uint64_t chunk_size) {
  gt_bofprintf(buffered_output,"@"PRIgts,PRIgts_content(tag));
  if (print_segmented_read_info) gt_output_bofprint_segmented_read_info(buffered_output,segment_id,total_segments); // Segmented Read
  if (left_trim > 0) {
    gt_bofprintf(buffered_output," lt:Z:%"PRIu64":"PRIgts":"PRIgts,left_trim,
        PRIgts_range_content(read,0,left_trim),
        PRIgts_range_content(qualities,0,left_trim)); // Left-trim
  }
  if (right_trim > 0) {
    gt_bofprintf(buffered_output," rt:Z:%"PRIu64":"PRIgts":"PRIgts,right_trim,
        PRIgts_range_content(read,left_trim+chunk_size,right_trim),
        PRIgts_range_content(qualities,left_trim+chunk_size,right_trim)); // Right-trim
  }
  // Print READ + QUALITIES (trimmed)
  gt_bofprintf(buffered_output,"\n"PRIgts"\n+\n"PRIgts"\n",
      PRIgts_trimmed_content(read,left_trim,right_trim),
      PRIgts_trimmed_content(qualities,left_trim,right_trim));
}
GT_INLINE void gt_filter_sample_read_print_fasta(
    gt_buffered_output_file* const buffered_output,gt_string* const tag,gt_string* const read,
    const bool print_segmented_read_info,const uint64_t segment_id,const uint64_t total_segments,
    const uint64_t left_trim,const uint64_t right_trim,const uint64_t chunk_size) {
  gt_bofprintf(buffered_output,">"PRIgts,PRIgts_content(tag));
  if (print_segmented_read_info) gt_output_bofprint_segmented_read_info(buffered_output,segment_id,total_segments); // Segmented Read
  if (left_trim > 0) {
    gt_bofprintf(buffered_output," lt:Z:%"PRIu64":"PRIgts,left_trim,
        PRIgts_range_content(read,0,left_trim)); // Left-trim
  }
  if (right_trim > 0) {
    gt_bofprintf(buffered_output," rt:Z:%"PRIu64":"PRIgts,right_trim,
        PRIgts_range_content(read,left_trim+chunk_size,right_trim)); // Right-trim
  }
  // Print READ (trimmed)
  gt_bofprintf(buffered_output,"\n"PRIgts"\n",
      PRIgts_trimmed_content(read,left_trim,right_trim));
}
GT_INLINE void gt_filter_group_reads() {
  // Open file IN/OUT
  gt_input_file* input_file = gt_filter_open_input_file(parameters.name_input_file);
  gt_output_file* output_file = gt_filter_open_output_file(parameters.name_output_file);
  // Prepare out-printers
  if (parameters.output_format==FILE_FORMAT_UNKNOWN) parameters.output_format = input_file->file_format; // Select output format
  gt_generic_printer_attributes* const generic_printer_attributes = gt_generic_printer_attributes_new(parameters.output_format);
  // SegmentedRead aux variables
  gt_template* const group_template = gt_template_new();
  uint64_t total_segments = 0, last_segment_id = 0;
  GT_BEGIN_READING_WRITING_LOOP(input_file,output_file,parameters.paired_end,buffered_output,template) {
    // Get group attribute
    gt_segmented_read_info* const segmented_read_info = gt_attributes_get_segmented_read_info(template->attributes);
    if (segmented_read_info==NULL) {
      gt_filter_cond_fatal_error_msg(total_segments!=last_segment_id,
          "Expected SegmentedRead Info => lastRead(%"PRIu64"/%"PRIu64")",last_segment_id,total_segments);
      gt_template_restore_trim(template); // If any
      GT_TEMPLATE_ITERATE_ALIGNMENT(group_template,alignment) {
        gt_attributes_remove(alignment->attributes,GT_ATTR_ID_SEGMENTED_READ_INFO); // If any
      }
      gt_output_generic_bofprint_template(buffered_output,template,generic_printer_attributes); // Print it, as it is
    } else {
      // First, undo the trim
      gt_template_restore_trim(template);
      // Tackle the group merging
      if (last_segment_id==total_segments) {
        /*
         * New group
         */
        gt_filter_cond_fatal_error_msg(segmented_read_info->total_segments==0 || segmented_read_info->segment_id!=1,
            "Wrong SegmentedRead Info (Zero reads in group or not properly sorted)");
        gt_template_clear(group_template,true);
        gt_template_copy(group_template,template,true,true);
        total_segments = segmented_read_info->total_segments;
        last_segment_id = segmented_read_info->segment_id;
      } else if (segmented_read_info->segment_id==last_segment_id+1 && segmented_read_info->segment_id <= total_segments) {
        /*
         * Old group (Keep merging)
         */
        gt_filter_cond_fatal_error_msg(!gt_string_equals(template->tag,group_template->tag),
            "Wrong TAG in Segmented Reads Sequence ('"PRIgts"'/'"PRIgts"')",PRIgts_content(group_template->tag),PRIgts_content(template->tag));
        gt_template_merge_template_mmaps(group_template,template);
        last_segment_id = segmented_read_info->segment_id;
        if (last_segment_id==total_segments) { // Close group
          GT_TEMPLATE_ITERATE_ALIGNMENT(group_template,alignment) {
            gt_attributes_remove(alignment->attributes,GT_ATTR_ID_SEGMENTED_READ_INFO); // If any
          }
          gt_output_generic_bofprint_template(buffered_output,group_template,generic_printer_attributes);
        }
      } else {
        gt_filter_fatal_error_msg("Wrong SegmentedRead Info => Expected(%"PRIu64"/%"PRIu64")::Found(%"PRIu64"/%"PRIu64").",
            segmented_read_info->segment_id,segmented_read_info->total_segments,last_segment_id,total_segments);
      }
    }
  } GT_END_READING_WRITING_LOOP(input_file,output_file,template);
  // Check proper end of merging groups
  gt_filter_cond_fatal_error_msg(total_segments!=last_segment_id,
      "Expected SegmentedRead Info => lastRead(%"PRIu64"/%"PRIu64")",last_segment_id,total_segments);
  // Clean
  gt_template_delete(group_template);
  gt_generic_printer_attributes_delete(generic_printer_attributes);
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
}
GT_INLINE void gt_filter_sample_read() {
  // Open file IN/OUT
  gt_input_file* input_file = gt_filter_open_input_file(parameters.name_input_file);
  gt_output_file* output_file = gt_filter_open_output_file(parameters.name_output_file);
  // Parallel I/O
#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(parameters.num_threads)
#endif
  {
    GT_BEGIN_READING_WRITING_LOOP(input_file,output_file,parameters.paired_end,buffered_output,template) {
        GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
          // Calculate the chunks
          const uint64_t read_length = gt_alignment_get_read_length(alignment);
          const uint64_t split_chunk_size = gt_get_integer_proportion(parameters.split_chunk_size,read_length);
          const uint64_t split_min_remainder = gt_get_integer_proportion(parameters.split_min_remainder,read_length);
          // Check boundaries
          if (split_chunk_size >= read_length || split_chunk_size <= split_min_remainder) {
            if (gt_alignment_has_qualities(alignment)) {
              gt_filter_sample_read_print_fastq(buffered_output,alignment->tag,alignment->read,alignment->qualities,false,1,1,0,0,read_length); // FASTQ
            } else {
              gt_filter_sample_read_print_fasta(buffered_output,alignment->tag,alignment->read,false,1,1,0,0,read_length); // FASTA
            }
            continue;
          }
          uint64_t split_step_size = gt_get_integer_proportion(parameters.split_step_size,read_length);
          if (split_step_size==0) split_step_size=1;
          const uint64_t split_left_trim = gt_get_integer_proportion(parameters.split_left_trim,read_length);
          const uint64_t split_right_trim = gt_get_integer_proportion(parameters.split_right_trim,read_length);
          const uint64_t full_chunks = ((read_length-split_left_trim-split_right_trim-split_chunk_size)/split_step_size)+1;
          uint64_t total_chunks = full_chunks;
          uint64_t left_trim=split_left_trim, right_trim=read_length-split_left_trim-split_chunk_size;
          // Check last chunk (remainder)
          const uint64_t last_left_trim = left_trim+(split_step_size*full_chunks);
          const uint64_t remainder_chunk = read_length-split_right_trim-last_left_trim;
          bool print_remainder_chunk = false;
          if (remainder_chunk > 0 && split_min_remainder > 0 &&
              remainder_chunk < split_chunk_size && remainder_chunk >= split_min_remainder) {
            print_remainder_chunk = true; ++total_chunks;
          }
          uint64_t i;
          for (i=0;i<full_chunks;++i,left_trim+=split_step_size,right_trim-=split_step_size) {
            if (gt_alignment_has_qualities(alignment)) {
              gt_filter_sample_read_print_fastq(
                  buffered_output,alignment->tag,alignment->read,alignment->qualities,true,
                  i+1,total_chunks,left_trim,right_trim,split_chunk_size); // FASTQ
            } else {
              gt_filter_sample_read_print_fasta(
                  buffered_output,alignment->tag,alignment->read,true,
                  i+1,total_chunks,left_trim,right_trim,split_chunk_size); // FASTA
            }
          }
          // Print last chunk (remainder)
          if (print_remainder_chunk) {
            if (gt_alignment_has_qualities(alignment)) {
              gt_filter_sample_read_print_fastq(
                  buffered_output,alignment->tag,alignment->read,alignment->qualities,true,
                  total_chunks,total_chunks,last_left_trim,split_right_trim,remainder_chunk); // FASTQ
            } else {
              gt_filter_sample_read_print_fasta(
                  buffered_output,alignment->tag,alignment->read,true,
                  total_chunks,total_chunks,last_left_trim,split_right_trim,remainder_chunk); // FASTA
            }
          }
        }
    } GT_END_READING_WRITING_LOOP(input_file,output_file,template);
  }
  // Clean
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
}
GT_INLINE void gt_filter_print_insert_size_distribution() {
  // Open file IN/OUT
  gt_input_file* input_file = gt_filter_open_input_file(parameters.name_input_file);
  gt_output_file* output_file = gt_filter_open_output_file(parameters.name_output_file);
  // Parallel I/O
#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(parameters.num_threads)
#endif
  {
    GT_BEGIN_READING_WRITING_LOOP(input_file,output_file,parameters.paired_end,buffered_output,template) {
      // Print insert size
      if (gt_template_get_num_blocks(template)!=2) continue;
      GT_TEMPLATE_ITERATE_(template,mmap) {
        gt_status error_code;
        gt_bofprintf(buffered_output,"%"PRIu64"\n",gt_template_get_insert_size(mmap,&error_code,0,0));
        if (parameters.first_map) break;
      }
    } GT_END_READING_WRITING_LOOP(input_file,output_file,template);
  }
  // Clean
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
}
GT_INLINE void gt_filter_print_error_distribution() {
  // Open file IN/OUT
  gt_input_file* input_file = gt_filter_open_input_file(parameters.name_input_file);
  gt_output_file* output_file = gt_filter_open_output_file(parameters.name_output_file);
  // Parallel I/O
#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(parameters.num_threads)
#endif
  {
    GT_BEGIN_READING_WRITING_LOOP(input_file,output_file,parameters.paired_end,buffered_output,template) {
      // Print levenshtein distance of the maps
      if (parameters.first_map)  {
        uint64_t best_distance = UINT64_MAX;
        GT_TEMPLATE_ITERATE_(template,mmap) {
          const uint64_t dist = gt_map_get_global_levenshtein_distance(*mmap);
          if (dist < best_distance) best_distance = dist;
        }
        if (best_distance < UINT64_MAX) gt_bofprintf(buffered_output,"%"PRIu64"\n",best_distance);
      } else {
        GT_TEMPLATE_ITERATE_(template,mmap) {
          gt_bofprintf(buffered_output,"%"PRIu64"\n",gt_map_get_global_levenshtein_distance(*mmap));
        }
      }
    } GT_END_READING_WRITING_LOOP(input_file,output_file,template);
  }
  // Clean
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
}
/*
 * Handler for opening an archive (GEMIndex/MULTIFastaFile)
 */
gt_sequence_archive* gt_filter_open_sequence_archive(const bool load_sequences) {
  gt_sequence_archive* sequence_archive = NULL;
  gt_log("Loading reference file ...");
  if (parameters.name_gem_index_file!=NULL) { // Load GEM-IDX
    sequence_archive = gt_sequence_archive_new(GT_BED_ARCHIVE);
    gt_gemIdx_load_archive(parameters.name_gem_index_file,sequence_archive,load_sequences);
  } else {
    gt_input_file* const reference_file = gt_input_file_open(parameters.name_reference_file,false);
    sequence_archive = gt_sequence_archive_new(GT_CDNA_ARCHIVE);
    if (gt_input_multifasta_parser_get_archive(reference_file,sequence_archive)!=GT_IFP_OK) {
      gt_fatal_error_msg("Error parsing reference file '%s'\n",parameters.name_reference_file);
    }
    gt_input_file_close(reference_file);
  }
  gt_log("Done.");
  return sequence_archive;
}
GT_INLINE void gt_filter_display_sequence_list(){
  // Show sequence archive summary
  gt_sequence_archive* sequence_archive = gt_filter_open_sequence_archive(false);
  gt_sequence_archive_iterator sequence_archive_it;
  gt_sequence_archive_new_iterator(sequence_archive,&sequence_archive_it);
  gt_segmented_sequence* seq;
  while ((seq=gt_sequence_archive_iterator_next(&sequence_archive_it))) {
    fprintf(stdout,"%s\t%"PRIu64"\n",seq->seq_name->buffer,seq->sequence_total_length);
  }
}
/*
 * I/O Filtering Loop
 */
#define GT_FILTER_CHECK_PARSING_ERROR(FORMAT) \
  ++record_num; \
  if (error_code!=GT_INPUT_STATUS_OK) { \
    gt_error_msg("[#%"PRIu64"]Fatal error parsing "FORMAT"file '%s', line %"PRIu64"\n", \
        record_num,parameters.name_input_file,buffered_input->current_line_num-1); \
    continue; \
  }
void gt_filter_read__write() {
  // Open in file
  gt_input_file* input_file = gt_filter_open_input_file(parameters.name_input_file);
  gt_output_file* output_file, *dicarded_output_file;

  // Open out file
  if (!parameters.no_output) {
    output_file = gt_filter_open_output_file(parameters.name_output_file);
    if (parameters.discarded_output) {
      if (gt_streq(parameters.name_discarded_output_file,"stdout")) {
        dicarded_output_file = gt_output_stream_new(stdout,SORTED_FILE);
      } else if (gt_streq(parameters.name_discarded_output_file,"stderr")) {
        dicarded_output_file = gt_output_stream_new(stderr,SORTED_FILE);
      } else {
        dicarded_output_file = gt_output_file_new(parameters.name_discarded_output_file,SORTED_FILE);
      }
    }
  }

  // Open reference file
  gt_sequence_archive* sequence_archive = NULL;
  if (parameters.load_index) {
    sequence_archive = gt_filter_open_sequence_archive(true);
  }

  // Parallel reading+process
  uint64_t total_algs_checked=0, total_algs_correct=0, total_maps_checked=0, total_maps_correct=0;
#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(parameters.num_threads) reduction(+:total_algs_checked,total_algs_correct,total_maps_checked,total_maps_correct)
#endif
  {
    // Prepare IN/OUT buffers & printers
    gt_status error_code;
    gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input_file);
    gt_buffered_output_file *buffered_output = NULL, *buffered_discarded_output = NULL;
    if (!parameters.no_output) {
      buffered_output = gt_buffered_output_file_new(output_file);
      gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_output);
      if (parameters.discarded_output) {
        buffered_discarded_output = gt_buffered_output_file_new(dicarded_output_file);
        gt_buffered_input_file_attach_buffered_output(buffered_input,buffered_discarded_output);
      }
    }
    // Prepare IN/OUT parser/printer attributes
    gt_generic_printer_attributes *generic_printer_attributes=NULL, *discarded_output_attributes=NULL;
    if (parameters.output_format==FILE_FORMAT_UNKNOWN) parameters.output_format = input_file->file_format; // Select output format
    generic_printer_attributes = gt_generic_printer_attributes_new(parameters.output_format);
    if (parameters.discarded_output) {
      gt_file_format output_format = input_file->file_format;
      if (parameters.discarded_output_format!=FILE_FORMAT_UNKNOWN) output_format=parameters.discarded_output_format;
      discarded_output_attributes = gt_generic_printer_attributes_new(output_format);
    }
    /*
     * READ + PROCCESS Loop
     */
    uint64_t record_num = 0;
    gt_template* template = gt_template_new();
    if (parameters.check_format && parameters.check_file_format==FASTA) {
      /*
       * FASTA I/O loop
       */
      while ((error_code=gt_input_fasta_parser_get_template(buffered_input,template,parameters.paired_end))) {
        GT_FILTER_CHECK_PARSING_ERROR("FASTA ");
        // Apply all filters and print
        gt_filter__print(input_file->file_format,buffered_input->current_line_num-1,sequence_archive,template,
            &total_algs_checked,&total_algs_correct,&total_maps_checked,&total_maps_correct,
            buffered_output,generic_printer_attributes,buffered_discarded_output,discarded_output_attributes);
      }
    } else if (parameters.check_format && parameters.check_file_format==MAP) {
      /*
       * MAP I/O loop
       */
      gt_map_parser_attributes* const attr = gt_input_map_parser_attributes_new(parameters.paired_end);
      while ((error_code=gt_input_map_parser_get_template(buffered_input,template,attr))) {
        GT_FILTER_CHECK_PARSING_ERROR("MAP ");
        // Apply all filters and print
        gt_filter__print(input_file->file_format,buffered_input->current_line_num-1,sequence_archive,template,
            &total_algs_checked,&total_algs_correct,&total_maps_checked,&total_maps_correct,
            buffered_output,generic_printer_attributes,buffered_discarded_output,discarded_output_attributes);
      }
      gt_input_map_parser_attributes_delete(attr);
    } else if (parameters.check_format && parameters.check_file_format==SAM) {
      /*
       * SAM I/O loop
       */
      gt_sam_parser_attributes* const attr = gt_input_sam_parser_attributes_new();
      while ((error_code=gt_input_sam_parser_get_template(buffered_input,template,attr))) {
        GT_FILTER_CHECK_PARSING_ERROR("SAM ");
        // Apply all filters and print
        gt_filter__print(input_file->file_format,buffered_input->current_line_num-1,sequence_archive,template,
            &total_algs_checked,&total_algs_correct,&total_maps_checked,&total_maps_correct,
            buffered_output,generic_printer_attributes,buffered_discarded_output,discarded_output_attributes);
      }
      gt_input_sam_parser_attributes_delete(attr);
    } else {
      /*
       * Generic I/O loop
       */
      gt_generic_parser_attributes* generic_parser_attributes = gt_input_generic_parser_attributes_new(parameters.paired_end);
      gt_input_map_parser_attributes_set_max_parsed_maps(generic_parser_attributes->map_parser_attributes,parameters.max_input_matches); // Limit max-matches
      while ((error_code=gt_input_generic_parser_get_template(buffered_input,template,generic_parser_attributes))) {
        GT_FILTER_CHECK_PARSING_ERROR("");
        // Apply all filters and print
        gt_filter__print(input_file->file_format,buffered_input->current_line_num-1,sequence_archive,template,
            &total_algs_checked,&total_algs_correct,&total_maps_checked,&total_maps_correct,
            buffered_output,generic_printer_attributes,buffered_discarded_output,discarded_output_attributes);
      }
      gt_input_generic_parser_attributes_delete(generic_parser_attributes);
    }
    // Clean
    gt_template_delete(template);
    gt_buffered_input_file_close(buffered_input);
    gt_generic_printer_attributes_delete(generic_printer_attributes);
    if (!parameters.no_output) {
      gt_buffered_output_file_close(buffered_output);
      if (parameters.discarded_output) gt_buffered_output_file_close(buffered_discarded_output);
    }
  }
  /*
   * Print check report
   */
  if (parameters.check) {
    gt_log("Checked %lu alignments. Total.Correct %lu (%2.3f %%). Total.Maps.Correct %lu (%2.3f %%)",
        total_algs_checked,total_algs_correct,GT_GET_PERCENTAGE(total_algs_correct,total_algs_checked),
        total_maps_correct,GT_GET_PERCENTAGE(total_maps_correct,total_maps_checked));
  }
  // Release archive & Clean
  if (sequence_archive) gt_sequence_archive_delete(sequence_archive);
  gt_filter_delete_map_ids(parameters.map_ids);
  if (parameters.quality_score_ranges!=NULL) gt_vector_delete(parameters.quality_score_ranges);
  gt_input_file_close(input_file);
  if (!parameters.no_output) {
    gt_output_file_close(output_file);
    if (parameters.discarded_output)  gt_output_file_close(dicarded_output_file);
  }
}
/*
 * Argument Parsing
 */
void gt_filter_get_coma_separated_arguments_long(char* const parameters_list,const uint64_t num_params,...) {
  uint64_t num_params_parsed = 0;
  // Start va_args
  va_list v_args;
  va_start(v_args,num_params);
  // Start parsing
  char *opt = strtok(parameters_list,",");
  while (opt!=NULL && num_params_parsed<num_params) {
    uint64_t* const uint64_arg = va_arg(v_args,uint64_t*);
    *uint64_arg = atoll(opt);
    opt = strtok(NULL,",");
  }
  // End va_args
  va_end(v_args);
}
GT_INLINE uint64_t gt_filter_get_coma_separated_arguments_float(char* const parameters_list,const uint64_t num_params,...) {
  uint64_t num_params_parsed = 0;
  // Start va_args
  va_list v_args;
  va_start(v_args,num_params);
  // Start parsing
  char *opt = strtok(parameters_list,",");
  while (opt!=NULL && num_params_parsed<num_params) {
    float* const float_arg = va_arg(v_args,float*);
    *float_arg = atof(opt);
    opt = strtok(NULL,",");
    ++num_params_parsed;
  }
  // End va_args
  va_end(v_args);
  return num_params_parsed;
}
void gt_filter_get_discarded_output_arguments(char* const optarg) {
  // Start parsing
  char *opt = strtok(optarg,",");
  parameters.name_discarded_output_file = opt;
  opt = strtok(NULL,","); // Next
  if (opt!=NULL) {
    if (gt_streq(opt,"FASTA")) {
      parameters.discarded_output_format = FASTA;
    } else if (gt_streq(opt,"MAP")) {
      parameters.discarded_output_format = MAP;
    } else if (gt_streq(opt,"SAM")) {
      parameters.discarded_output_format = SAM;
    } else {
      gt_fatal_error_msg("Output format '%s' not recognized",opt);
    }
  }
}
void gt_filter_get_argument_pair_strandness(char* const strandness_opt) {
  char *opt;
  opt = strtok(strandness_opt,",");
  while (opt!=NULL) {
    if (gt_streq(opt,"FR")) {
      parameters.allow_strand_fr = true;
    } else if (gt_streq(opt,"RF")) {
      parameters.allow_strand_rf = true;
    } else if (gt_streq(opt,"FF")) {
      parameters.allow_strand_ff = true;
    } else if (gt_streq(opt,"RR")) {
      parameters.allow_strand_rr = true;
    } else {
      gt_fatal_error_msg("Strandedness option not recognized '%s'\n",opt);
    }
    opt = strtok(NULL,","); // Reload
  }
  parameters.filter_by_strand_pe = true;
}
void gt_filter_get_argument_map_id(char* const maps_ids) {
  // Allocate vector
  parameters.map_ids = gt_vector_new(20,sizeof(gt_string*));
  // Add all the valid map Ids (sequence names)
  char *opt;
  opt = strtok(maps_ids,",");
  while (opt!=NULL) {
    // Get id
    gt_string* map_id = gt_string_new(0);
    gt_string_set_string(map_id,opt);
    // Add to the vector
    gt_vector_insert(parameters.map_ids,map_id,gt_string*);
    // Next
    opt = strtok(NULL,","); // Reload
  }
}
void parse_arguments(int argc,char** argv) {
  struct option* gt_filter_getopt = gt_options_adaptor_getopt(gt_filter_options);
  gt_string* const gt_filter_short_getopt = gt_options_adaptor_getopt_short(gt_filter_options);
  int option, option_index;
  while (true) {
    // Get option &  Select case
    if ((option=getopt_long(argc,argv,
        gt_string_get_string(gt_filter_short_getopt),gt_filter_getopt,&option_index))==-1) break;
    switch (option) {
    /* I/O */
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 'o':
      parameters.name_output_file = optarg;
      if (gt_streq(optarg,"null")) parameters.no_output = true;
      break;
    case 'r':
      parameters.name_reference_file = optarg;
      break;
    case 'I':
      parameters.name_gem_index_file = optarg;
      break;
    case 200:
      parameters.mmap_input = true;
      break;
    case 'z': // gzip-input
      parameters.input_compression = GT_COMPRESSION_GZIP;
      break;
    case 'j': // bzip-input
      parameters.input_compression = GT_COMPRESSION_BZIP;
      break;
    case 'p':
      parameters.paired_end = true;
      break;
    case 201: // output-format
      if (gt_streq(optarg,"FASTA")) {
        parameters.output_format = FASTA;
      } else if (gt_streq(optarg,"MAP")) {
        parameters.output_format = MAP;
      } else if (gt_streq(optarg,"SAM")) {
        parameters.output_format = SAM;
      } else {
        gt_fatal_error_msg("Output format '%s' not recognized",optarg);
      }
      break;
    case 202: // output-compression
      if (gt_streq(optarg,"NONE")) {
        parameters.output_compression = GT_COMPRESSION_NONE;
      } else if (gt_streq(optarg,"GZIP")) {
        parameters.output_compression = GT_COMPRESSION_GZIP;
      } else if (gt_streq(optarg,"BZIP")) {
        parameters.output_compression = GT_COMPRESSION_BZIP;
      } else {
        gt_fatal_error_msg("Output compression '%s' not recognized",optarg);
      }
      break;
    case 203: // discarded-output
      parameters.discarded_output = true;
      gt_filter_get_discarded_output_arguments(optarg);
      break;
    case 204: // no-output
      parameters.no_output = true;
      break;
    /* Filter Read/Qualities */
    case 300: // hard-trim
      parameters.hard_trim = true;
      gt_filter_get_coma_separated_arguments_long(optarg,2,&(parameters.left_trim),&(parameters.right_trim));
      break;
    case 301: // quality-trim
      gt_fatal_error(NOT_IMPLEMENTED);
      break;
    case 302: // restore-trim
      parameters.restore_trim = true;
      break;
    case 303: // uniform-read
      parameters.uniform_read = true;
      if (optarg && gt_streq(optarg,"strict")) parameters.uniform_read_strict = true;
      break;
    case 304: // qualities-to-offset-33
      parameters.qualities_to_offset_33 = true;
      break;
    case 305: // qualities-to-offset-64
      parameters.qualities_to_offset_64 = true;
      break;
    case 306: // remove-qualities
      parameters.remove_qualities = true;
      break;
    case 307: // add-qualities
      parameters.add_qualities = true;
      break;
    /* Filter Template/Alignments */
    case 400:
      parameters.mapped = true;
      break;
    case 401:
      parameters.unmapped = true;
      break;
    case 402:
      parameters.unique_level = atoll(optarg);
      break;
    case 403:
      parameters.min_length = atof(optarg);
      break;
    case 404:
      parameters.max_length = atof(optarg);
      break;
    case 405:
      parameters.min_maps = atof(optarg);
      break;
    case 406:
      parameters.max_maps = atof(optarg);
      break;
    /* Filter Maps */
    case 500: // first-map
      parameters.perform_dna_map_filter = true;
      parameters.first_map = true;
      break;
    case 'k': // keep-first-map
      parameters.keep_first_map = true;
      break;
    case 'd': // max-decoded-matches
      parameters.matches_pruning = true;
      parameters.max_decoded_matches = atoll(optarg);
      break;
    case 'D': // min-decoded-strata
      parameters.matches_pruning = true;
      parameters.min_decoded_strata = atoll(optarg);
      break;
    case 501: // max-output-matches
      parameters.matches_pruning = true;
      parameters.max_output_matches = atoll(optarg);
      break;
    case 502: // max-input-matches
      parameters.max_input_matches = atoll(optarg);
      break;
    case 503: // max-strata-after-map
      parameters.perform_dna_map_filter = true;
      parameters.max_strata_after_map = atof(optarg);
      break;
    case 504: // make-counters
      parameters.make_counters = true;
      break;
    case 505: // min-strata
      parameters.perform_dna_map_filter = true;
      parameters.min_event_distance = atof(optarg);
      break;
    case 506: // max-strata
      parameters.perform_dna_map_filter = true;
      parameters.max_event_distance = atof(optarg);
      break;
    case 507: // min-levenshtein-error
      parameters.perform_dna_map_filter = true;
      parameters.min_levenshtein_distance = atof(optarg);
      break;
    case 508: // max-levenshtein-error
      parameters.perform_dna_map_filter = true;
      parameters.max_levenshtein_distance = atof(optarg);
      break;
    case 509: // map-id
      parameters.perform_dna_map_filter = true;
      gt_filter_get_argument_map_id(optarg);
      break;
    case 510: // strandedness
      parameters.perform_dna_map_filter = true;
      parameters.filter_by_strand_se = true;
      if (gt_streq(optarg,"F")) {
        parameters.allow_strand_f = true;
      } else if (gt_streq(optarg,"R")) {
        parameters.allow_strand_r = true;
      } else {
        gt_fatal_error_msg("Strand '%s' not recognized {'F','R'}",optarg);
      }
      break;
    case 511: // filter-quality
      parameters.perform_dna_map_filter = true;
      gt_filter_quality_range qrange;
      gt_filter_get_coma_separated_arguments_long(optarg,2,&(qrange.min),&(qrange.max));
      // Add it to the vector of ranges
      if (parameters.quality_score_ranges==NULL) {
        parameters.quality_score_ranges = gt_vector_new(4,sizeof(gt_filter_quality_range));
      }
      gt_vector_insert(parameters.quality_score_ranges,qrange,gt_filter_quality_range);
      break;
    /* Filter RNA-Maps */
    case 600: // no-split-maps
      parameters.no_split_maps = true;
      parameters.perform_rna_map_filter = true;
      break;
    case 601: // only-split-maps
      parameters.only_split_maps = true;
      parameters.perform_rna_map_filter = true;
      break;
    case 603: // min-intron-length
      parameters.min_intron_length = atol(optarg);
      parameters.perform_rna_map_filter = true;
      break;
    case 604: // min-block-length
      parameters.min_block_length = atol(optarg);
      parameters.perform_rna_map_filter = true;
      break;
    /* Filter PE-Maps */
    case 700: // pair-strandness
      parameters.perform_dna_map_filter = true;
      gt_filter_get_argument_pair_strandness(optarg);
      break;
    case 701: // min-inss
      parameters.perform_dna_map_filter = true;
      parameters.min_inss = atoll(optarg);
      break;
    case 702: // max-inss
      parameters.perform_dna_map_filter = true;
      parameters.max_inss = atoll(optarg);
      break;
    /* Realign/Check */
    case 800: // mismatch-recovery
      parameters.load_index = true;
      parameters.mismatch_recovery = true;
      break;
    case 801: // hamming-realign
      parameters.load_index = true;
      parameters.realign_hamming = true;
      break;
    case 802: // levenshtein-realign
      parameters.load_index = true;
      parameters.realign_levenshtein = true;
      break;
    /* Checking/Report */
    case 'c': // check
      parameters.load_index = true;
      parameters.check = true;
      break;
    case 'C': // check-only
      parameters.load_index = true;
      parameters.check = true;
      parameters.no_output = true;
      break;
    case 803: // check-format
      parameters.check_format = true;
      if (gt_streq(optarg,"FASTA")) {
        parameters.check_file_format = FASTA;
      } else if (gt_streq(optarg,"MAP")) {
        parameters.check_file_format = MAP;
      } else if (gt_streq(optarg,"SAM")) {
        parameters.check_file_format = SAM;
      } else {
        gt_fatal_error_msg("Check format '%s' not recognized",optarg);
      }
      break;
      /* Split/Grouping */
    case 900: // split-read
      gt_fatal_error(NOT_IMPLEMENTED);
      break;
    case 901: // sample-read
      parameters.special_functionality = true;
      parameters.sample_read = true;
      gt_cond_fatal_error_msg(gt_filter_get_coma_separated_arguments_float(optarg,5,
          &(parameters.split_chunk_size),&(parameters.split_step_size),
          &(parameters.split_left_trim),&(parameters.split_right_trim),
          &(parameters.split_min_remainder))<4,
          "Too few parameters provided to option --split-read");
      break;
    case 902: // group-read-chunks
      parameters.special_functionality = true;
      parameters.group_reads = true;
      break;
    /* Display/Information */
    case 1000:
      parameters.special_functionality = true;
      parameters.error_plot = true;
      break;
    case 1001:
      parameters.special_functionality = true;
      parameters.insert_size_plot = true;
      break;
    case 1002:
      parameters.special_functionality = true;
      parameters.load_index = true;
      parameters.show_sequence_list = true;
      break;
    case 1003:
      parameters.special_functionality = true;
      parameters.load_index = true;
      parameters.display_pretty = true;
      break;
    /* Misc */
    case 't': // threads
#ifdef HAVE_OPENMP
      parameters.num_threads = atol(optarg);
#endif
      gt_cond_fatal_error_msg(parameters.num_threads > GT_MAX_OUTPUT_BUFFERS,
          "Excessive number of threads (maximum %"PRId32")",GT_MAX_OUTPUT_BUFFERS);
      break;
    case 'v': // verbose
      parameters.verbose = true;
      break;
    case 'h': // help
      fprintf(stderr, "USE: ./gt.filter [ARGS]...\n");
      gt_options_fprint_menu(stderr,gt_filter_options,gt_filter_groups,false,false);
      exit(1);
    case 'H': // full-help
      fprintf(stderr, "USE: ./gt.filter [ARGS]...\n");
      gt_options_fprint_menu(stderr,gt_filter_options,gt_filter_groups,false,true);
      exit(1);
    case 'J':
      gt_options_fprint_json_menu(stderr,gt_filter_options,gt_filter_groups,true,false);
      exit(1);
      break;
    case '?':
    default:
      gt_fatal_error_msg("Option not recognized");
    }
  }
  /*
   * Parameters check
   */
  if (parameters.load_index && parameters.name_reference_file==NULL && parameters.name_gem_index_file==NULL) {
    gt_fatal_error_msg("Reference file required");
  }
  // Free
  gt_string_delete(gt_filter_short_getopt);
}
/*
 * Main
 */
int main(int argc,char** argv) {
  // GT error handler
  gt_handle_error_signals();

  // Parsing command-line options
  parse_arguments(argc,argv);

  /*
   * Select functionality
   */
  if (parameters.show_sequence_list) {
    gt_filter_display_sequence_list();
  } else if (parameters.group_reads) {
    gt_filter_group_reads();
  } else if (parameters.sample_read) {
    gt_filter_sample_read();

  // Depreciated
  } else if (parameters.error_plot) {
    gt_filter_print_insert_size_distribution();
  } else if (parameters.insert_size_plot) {
    gt_filter_print_error_distribution();
  // Depreciated

  } else {
    gt_filter_read__write(); // Filter !!
  }

  return 0;
}

