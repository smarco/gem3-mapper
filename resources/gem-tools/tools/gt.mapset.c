/*
 * PROJECT: GEM-Tools library
 * FILE: gt.mapset.c
 * DATE: 08/11/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Utility to perform set operations {UNION,INTERSECTION,DIFFERENCE} over alignment files {MAP,SAM}
 */

#include <getopt.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "gem_tools.h"

/*
 * Constants
 */
#define GT_MAPSET_OPERATIONS "union,intersection,difference,compare,join,merge-map,display-compact"

/*
 * Data structures
 */
typedef enum { GT_MAP_SET_UNKNOWN, GT_MERGE_MAP,
               GT_MAP_SET_INTERSECTION, GT_MAP_SET_UNION, GT_MAP_SET_DIFFERENCE,
               GT_MAP_SET_JOIN, GT_MAP_SET_COMPARE, GT_DISPLAY_COMPACT_MAP,
               GT_MAP_SET_SPECIFICITY_PROFILE } gt_operation;
typedef struct {
  gt_operation operation;
  char* name_input_file_1;
  char* name_input_file_2;
  char* name_output_file;
  bool mmap_input;
  bool paired_end;
  bool files_contain_same_reads;
  double eq_threshold;
  uint8_t mapq_min;
  uint8_t mapq_max;
  bool strict;
  bool verbose;
  uint64_t num_threads;
} gt_stats_args;
gt_stats_args parameters = {
    .operation=GT_MAP_SET_UNKNOWN,
    .name_input_file_1=NULL,
    .name_input_file_2=NULL,
    .name_output_file=NULL,
    .mmap_input=false,
    .paired_end=false,
    .files_contain_same_reads=false,
    .eq_threshold=0.2,
    .mapq_min=0,
    .mapq_max=255,
    .strict=false,
    .verbose=false,
    .num_threads=1
};
/*
 * Global variables
 */
uint64_t current_read_length;
/*
 * Compare procedures
 */
int64_t gt_mapset_map_cmp(gt_map* const map_1,gt_map* const map_2) {
  const uint64_t eq_threshold = (parameters.eq_threshold <= 1.0) ?
      parameters.eq_threshold*current_read_length: parameters.eq_threshold;
  return parameters.strict ? gt_map_cmp(map_1,map_2) : gt_map_range_cmp(map_1,map_2,eq_threshold);
}
int64_t gt_mapset_mmap_cmp(gt_map** const map_1,gt_map** const map_2,const uint64_t num_maps) {
  const uint64_t eq_threshold = (parameters.eq_threshold <= 1.0) ?
      parameters.eq_threshold*current_read_length: parameters.eq_threshold;
  const int64_t cmp = parameters.strict ?
    gt_mmap_cmp(map_1,map_2,num_maps) :
    gt_mmap_range_cmp(map_1,map_2,num_maps,eq_threshold);
  if (num_maps==1 || cmp==0) return cmp;
  // Trick for when I screw it up generating datasets (inverting maps order flag)
  gt_map* map_2_inv[2] = {map_2[1], map_2[0]};
  return parameters.strict ?
      gt_mmap_cmp(map_1,(gt_map** const )map_2_inv,num_maps) :
      gt_mmap_range_cmp(map_1,(gt_map** const )map_2_inv,num_maps,eq_threshold);
}
/*
 * I/O Sync
 */
GT_INLINE gt_status gt_mapset_read_template_sync(
    gt_buffered_input_file* const buffered_input_master,gt_buffered_input_file* const buffered_input_slave,
    gt_buffered_output_file* const buffered_output,gt_template* const template_master,gt_template* const template_slave,
    const gt_operation operation) {
  // Read master
  gt_status error_code_master, error_code_slave;
  gt_output_map_attributes* output_attributes = gt_output_map_attributes_new();
  gt_generic_parser_attributes* generic_parser_attr = gt_input_generic_parser_attributes_new(parameters.paired_end);
  if ((error_code_master=gt_input_generic_parser_get_template(
      buffered_input_master,template_master,generic_parser_attr))==GT_INPUT_STATUS_FAIL) {
    gt_fatal_error_msg("Fatal error parsing file <<Master>>");
  }
  // Read slave
  if ((error_code_slave=gt_input_generic_parser_get_template(
      buffered_input_slave,template_slave,generic_parser_attr))==GT_INPUT_STATUS_FAIL) {
    gt_fatal_error_msg("Fatal error parsing file <<Slave>>");
  }
  // Check EOF conditions
  if (error_code_master==GT_INPUT_STATUS_EOF) {
    if (error_code_slave!=GT_INPUT_STATUS_EOF) {
      gt_fatal_error_msg("<<Slave>> contains more/different reads from <<Master>>");
    }
    return GT_INPUT_STATUS_EOF;
  } else if (error_code_slave==GT_INPUT_STATUS_EOF) { // Slave exhausted. Dump master & return EOF
    do {
      if (error_code_master==GT_INPUT_STATUS_FAIL) gt_fatal_error_msg("Fatal error parsing file <<Master>>");
      if (operation==GT_MAP_SET_UNION || operation==GT_MAP_SET_DIFFERENCE) {
        gt_output_map_bofprint_template(buffered_output,template_master,output_attributes);
      }
    } while ((error_code_master=gt_input_generic_parser_get_template(
                buffered_input_master,template_master,generic_parser_attr)));
    return GT_INPUT_STATUS_EOF;
  }
  // Synch loop
  while (gt_string_cmp(gt_template_get_string_tag(template_master),
      gt_template_get_string_tag(template_slave))) {
    // Print non correlative master's template
    if (operation==GT_MAP_SET_UNION || operation==GT_MAP_SET_DIFFERENCE) {
      gt_output_map_bofprint_template(buffered_output,template_master,output_attributes);
    }
    // Fetch next master's template
    if ((error_code_master=gt_input_generic_parser_get_template(
        buffered_input_master,template_master,generic_parser_attr))!=GT_INPUT_STATUS_OK) {
      gt_fatal_error_msg("<<Slave>> contains more/different reads from <<Master>>");
    }
  }
  return GT_INPUT_STATUS_OK;
}
GT_INLINE gt_status gt_mapset_read_template_get_commom_map(
    gt_buffered_input_file* const buffered_input_master,gt_buffered_input_file* const buffered_input_slave,
    gt_template* const template_master,gt_template* const template_slave) {
  gt_status error_code_master, error_code_slave;
  gt_generic_parser_attributes* generic_parser_attr = gt_input_generic_parser_attributes_new(parameters.paired_end);
  // Read master
  if ((error_code_master=gt_input_generic_parser_get_template(
      buffered_input_master,template_master,generic_parser_attr))==GT_INPUT_STATUS_FAIL) {
    gt_fatal_error_msg("Fatal error parsing file <<Master>>");
  }
  if (error_code_master==GT_INPUT_STATUS_EOF) return GT_INPUT_STATUS_EOF;
  // Read slave
  if ((error_code_slave=gt_input_generic_parser_get_template(
      buffered_input_slave,template_slave,generic_parser_attr))==GT_INPUT_STATUS_FAIL) {
    gt_fatal_error_msg("Fatal error parsing file <<Slave>>");
  }
  if (error_code_slave==GT_INPUT_STATUS_EOF) { // Check EOF conditions
    gt_fatal_error_msg("<<Slave>> is not contained in master <<Master>> (looking for '"PRIgts"')",
        PRIgts_content(gt_template_get_string_tag(template_master)));
  }
  // Synch loop
  while (gt_string_cmp(gt_template_get_string_tag(template_master),gt_template_get_string_tag(template_slave))) {
    // Fetch next slave's template
    if ((error_code_master=gt_input_generic_parser_get_template(
        buffered_input_slave,template_slave,generic_parser_attr))!=GT_INPUT_STATUS_OK) {
      gt_fatal_error_msg("<<Slave>> is not contained in master <<Master>> (looking for '"PRIgts"')",
          PRIgts_content(gt_template_get_string_tag(template_master)));
    }
  }
  return GT_INPUT_STATUS_OK;
}
/*
 * Set operations
 */
void gt_mapset_perform_set_operations() {
  // File IN/OUT
  gt_input_file* input_file_1 = gt_input_file_open(parameters.name_input_file_1,parameters.mmap_input);
  gt_input_file* input_file_2 = (parameters.name_input_file_2==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file_2,parameters.mmap_input);
  if (parameters.name_input_file_2==NULL) GT_SWAP(input_file_1,input_file_2);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);

  // Buffered I/O
  gt_buffered_input_file* buffered_input_1 = gt_buffered_input_file_new(input_file_1);
  gt_buffered_input_file* buffered_input_2 = gt_buffered_input_file_new(input_file_2);
  gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
  gt_buffered_input_file_attach_buffered_output(buffered_input_1,buffered_output);

  // Template I/O (synch)
  gt_template *template_1 = gt_template_new();
  gt_template *template_2 = gt_template_new();
  gt_output_map_attributes* output_attributes = gt_output_map_attributes_new();
  while (gt_mapset_read_template_sync(buffered_input_1,buffered_input_2,
      buffered_output,template_1,template_2,parameters.operation)) {
    // Record current read length
    current_read_length = gt_template_get_total_length(template_1);
    // Apply operation
    gt_template *ptemplate;
    switch (parameters.operation) {
      case GT_MAP_SET_UNION:
        ptemplate=gt_template_union_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_1,template_2);
        break;
      case GT_MAP_SET_INTERSECTION:
        ptemplate=gt_template_intersect_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_1,template_2);
//        // FIXME:: SNAP PE BUG
//        if (gt_template_get_num_mmaps(template_2)==0) {
//          gt_alignment* const tmp = template_2->alignment_end1;
//          template_2->alignment_end1 =  template_2->alignment_end2;
//          template_2->alignment_end2 = tmp;
//          GT_VECTOR_ITERATE(template_2->mmaps,mmap,n,gt_mmap) {
//            gt_map* tmp_map = mmap->mmap[0];
//            mmap->mmap[0] = mmap->mmap[1];
//            mmap->mmap[1] = tmp_map;
//          }
//          gt_template_delete(ptemplate);
//          ptemplate=gt_template_intersect_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_1,template_2);
//        }
//        // FIXME:: SNAP PE BUG
        break;
      case GT_MAP_SET_DIFFERENCE:
        ptemplate=gt_template_subtract_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_1,template_2);
        break;
      default:
        gt_fatal_error(SELECTION_NOT_VALID);
        break;
    }
    // Print template
    gt_output_map_bofprint_template(buffered_output,ptemplate,output_attributes);
    // Delete template
    gt_template_delete(ptemplate);
  }

  // Clean
  gt_template_delete(template_1);
  gt_template_delete(template_2);
  gt_buffered_input_file_close(buffered_input_1);
  gt_buffered_input_file_close(buffered_input_2);
  gt_buffered_output_file_close(buffered_output);
  gt_input_file_close(input_file_1);
  gt_input_file_close(input_file_2);
  gt_output_file_close(output_file);
}
/*
 * Compare operations
 */
void gt_mapset_perform_cmp_operations() {
  // File IN/OUT
  gt_input_file* input_file_1 = gt_input_file_open(parameters.name_input_file_1,parameters.mmap_input);
  gt_input_file* input_file_2 = (parameters.name_input_file_2==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file_2,parameters.mmap_input);
  if (parameters.name_input_file_2==NULL) GT_SWAP(input_file_1,input_file_2);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);

  // Buffered I/O
  gt_buffered_input_file* buffered_input_1 = gt_buffered_input_file_new(input_file_1);
  gt_buffered_input_file* buffered_input_2 = gt_buffered_input_file_new(input_file_2);
  gt_buffered_output_file* buffered_output = gt_buffered_output_file_new(output_file);
  gt_buffered_input_file_attach_buffered_output(buffered_input_1,buffered_output);

  // Print Banner
  if (parameters.operation == GT_MAP_SET_JOIN) {
    gt_bofprintf(buffered_output,"Tag\tCOUNTER(Master)\tCOUNTER(Slave)\tMAP(Master)\tMAP(Slave)\n");
  } else if (parameters.operation == GT_MAP_SET_COMPARE) {
    gt_bofprintf(buffered_output,"TAG\t"
          "COUNTER(Master-Slave)\tCOUNTER(Slave-Master)\tCOUNTER(Slave I Master)\t"
          "MAP(Master-Slave)\tMAP(Slave-Master)\tMAP(Slave I Master)\n");
  }

  // Template I/O (synch)
  gt_template *template_1 = gt_template_new();
  gt_template *template_2 = gt_template_new();
  gt_output_map_attributes* output_map_attributes = gt_output_map_attributes_new();
  while (gt_mapset_read_template_get_commom_map(buffered_input_1,buffered_input_2,template_1,template_2)) {
    // Record current read length
    current_read_length = gt_template_get_total_length(template_1);
    // Apply operation
    switch (parameters.operation) {
      case GT_MAP_SET_JOIN:
        // Print Master's TAG+Counters+Maps
        gt_output_map_bofprint_tag(buffered_output,template_1->tag,template_1->attributes,output_map_attributes);
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_1),
            template_1->attributes,output_map_attributes); // Master's Counters
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_2),
            template_1->attributes,output_map_attributes); // Slave's Counters
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_template_maps(buffered_output,template_1,output_map_attributes); // Master's Maps
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_template_maps(buffered_output,template_2,output_map_attributes); // Slave's Maps
        gt_bofprintf(buffered_output,"\n");
        break;
      case GT_MAP_SET_COMPARE: {
        // Perform simple cmp operations
        gt_template *template_master_minus_slave=gt_template_subtract_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_1,template_2);
        gt_template *template_slave_minus_master=gt_template_subtract_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_2,template_1);
        gt_template *template_intersection=gt_template_intersect_template_mmaps_fx(gt_mapset_mmap_cmp,gt_mapset_map_cmp,template_1,template_2);
        /*
         * Print results :: (TAG (Master-Slave){COUNTER MAPS} (Slave-Master){COUNTER MAPS} (Intersection){COUNTER MAPS})
         */
        gt_output_map_bofprint_tag(buffered_output,template_1->tag,template_1->attributes,output_map_attributes);
        // Counters
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_master_minus_slave),
            template_master_minus_slave->attributes,output_map_attributes); // (Master-Slave){COUNTER}
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_slave_minus_master),
            template_slave_minus_master->attributes,output_map_attributes); // (Slave-Master){COUNTER}
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_intersection),
            template_intersection->attributes,output_map_attributes); // (Intersection){COUNTER}
        // Maps
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_template_maps(buffered_output,template_master_minus_slave,output_map_attributes); // (Master-Slave){COUNTER}
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_template_maps(buffered_output,template_slave_minus_master,output_map_attributes); // (Slave-Master){COUNTER}
        gt_bofprintf(buffered_output,"\t");
        gt_output_map_bofprint_template_maps(buffered_output,template_intersection,output_map_attributes); // (Intersection){COUNTER}
        gt_bofprintf(buffered_output,"\n");
        // Delete templates
        gt_template_delete(template_master_minus_slave);
        gt_template_delete(template_slave_minus_master);
        gt_template_delete(template_intersection);
        }
        break;
      default:
        gt_fatal_error(SELECTION_NOT_VALID);
        break;
    }
  }
  // Clean
  gt_template_delete(template_1);
  gt_template_delete(template_2);
  gt_buffered_input_file_close(buffered_input_1);
  gt_buffered_input_file_close(buffered_input_2);
  gt_buffered_output_file_close(buffered_output);
  gt_input_file_close(input_file_1);
  gt_input_file_close(input_file_2);
  gt_output_file_close(output_file);
}
/*
 * Specificity Profile
 */
#define GT_CLASSES_TIE_PERFECT 0
#define GT_CLASSES_TIE         1
#define GT_CLASSES_MMAP_D1     2
#define GT_CLASSES_MMAP_D2     3
#define GT_CLASSES_MMAP_D3     4
#define GT_CLASSES_MMAP_D4     5
#define GT_CLASSES_MMAP_D5     6
#define GT_CLASSES_MMAP        7
#define GT_CLASSES_UNIQUE      8
#define GT_CLASSES_UNMAPPED    9
#define GT_CLASSES_RANGE      10
typedef struct {
  // Total
  uint64_t num_reads;
  uint64_t num_mapped;
  uint64_t num_maps;
  uint64_t num_signal_maps;
  uint64_t num_noise_maps;
  // Reads Hit
  uint64_t num_reads_hit_first;
  uint64_t num_reads_hit_first_ties;
  uint64_t num_reads_hit_any;
  // Reads Miss
  uint64_t num_reads_miss_ties;
  uint64_t num_reads_miss_2tie;
  uint64_t num_reads_miss_3tie;
  uint64_t num_reads_miss_4tie;
  uint64_t num_reads_miss_moretie;
  uint64_t num_reads_miss_unmapped;
  uint64_t num_reads_miss_misscored;
  uint64_t num_reads_miss_missing;
  // Scoring goodness
  uint64_t num_reads_score_sensitive;
  uint64_t num_reads_score_sensitive_first_hit;
  uint64_t num_reads_score_sensitive_event_hit;
  uint64_t num_reads_score_sensitive_edit_hits;
  // Maps Noise
  uint64_t num_noise_maps_ties;
  uint64_t num_noise_maps_misscored;
  uint64_t num_noise_maps_missing;
  // ROC
  uint64_t mapq_true[UINT8_MAX];
  uint64_t mapq_false[UINT8_MAX];
  // Classes
  uint64_t matches_classes_true[GT_CLASSES_RANGE];
  uint64_t matches_classes_false[GT_CLASSES_RANGE];
} specificity_profile_t;
GT_INLINE specificity_profile_t* gt_specificity_profile_new() {
  // Alloc
  specificity_profile_t* const specificity_profile = gt_alloc(specificity_profile_t);
  // Total
  specificity_profile->num_reads = 0;
  specificity_profile->num_mapped = 0;
  specificity_profile->num_maps = 0;
  specificity_profile->num_reads_hit_first_ties = 0;
  specificity_profile->num_signal_maps = 0;
  specificity_profile->num_noise_maps = 0;
  // Reads Hit
  specificity_profile->num_reads_hit_first = 0;
  specificity_profile->num_reads_hit_any = 0;
  // Reads Miss
  specificity_profile->num_reads_miss_ties = 0;
  specificity_profile->num_reads_miss_2tie = 0;
  specificity_profile->num_reads_miss_3tie = 0;
  specificity_profile->num_reads_miss_4tie = 0;
  specificity_profile->num_reads_miss_moretie = 0;
  specificity_profile->num_reads_miss_unmapped = 0;
  specificity_profile->num_reads_miss_misscored = 0;
  specificity_profile->num_reads_miss_missing = 0;
  // Scoring goodness
  specificity_profile->num_reads_score_sensitive = 0;
  specificity_profile->num_reads_score_sensitive_first_hit = 0;
  specificity_profile->num_reads_score_sensitive_event_hit = 0;
  specificity_profile->num_reads_score_sensitive_edit_hits = 0;
  // Maps Noise
  specificity_profile->num_noise_maps_ties = 0;
  specificity_profile->num_noise_maps_misscored = 0;
  specificity_profile->num_noise_maps_missing = 0;
  // ROC
  memset(specificity_profile->mapq_true,0,UINT8_MAX*sizeof(uint64_t));
  memset(specificity_profile->mapq_false,0,UINT8_MAX*sizeof(uint64_t));
  // Classes
  memset(specificity_profile->matches_classes_true,0,GT_CLASSES_RANGE*sizeof(uint64_t));
  memset(specificity_profile->matches_classes_false,0,GT_CLASSES_RANGE*sizeof(uint64_t));
  // Return
  return specificity_profile;
}
GT_INLINE void gt_specificity_profile_delete(specificity_profile_t* const specificity_profile) {
  gt_free(specificity_profile);
}
GT_INLINE void gt_specificity_profile_print(specificity_profile_t* const specificity_profile) {
  // Data
  const uint64_t num_reads = specificity_profile->num_reads;
  const uint64_t num_mapped = specificity_profile->num_mapped;
  const uint64_t num_reads_hit_first_ties = specificity_profile->num_reads_hit_first_ties;
  const uint64_t num_maps = specificity_profile->num_maps;
  const uint64_t num_signal_maps = specificity_profile->num_signal_maps;
  const uint64_t num_noise_maps = specificity_profile->num_noise_maps;
  // Reads Hit
  const uint64_t num_reads_hit_first = specificity_profile->num_reads_hit_first;
  const uint64_t num_reads_hit_any = specificity_profile->num_reads_hit_any;
  // Reads Miss
  const uint64_t num_reads_miss_ties = specificity_profile->num_reads_miss_ties;
  const uint64_t num_reads_miss_2tie = specificity_profile->num_reads_miss_2tie;
  const uint64_t num_reads_miss_3tie = specificity_profile->num_reads_miss_3tie;
  const uint64_t num_reads_miss_4tie = specificity_profile->num_reads_miss_4tie;
  const uint64_t num_reads_miss_moretie = specificity_profile->num_reads_miss_moretie;
  const uint64_t num_reads_miss_unmapped = specificity_profile->num_reads_miss_unmapped;
  const uint64_t num_reads_miss_misscored = specificity_profile->num_reads_miss_misscored;
  const uint64_t num_reads_miss_missing = specificity_profile->num_reads_miss_missing;
  // Scoring goodness
  const uint64_t num_reads_score_sensitive = specificity_profile->num_reads_score_sensitive;
  const uint64_t num_reads_score_sensitive_first_hit = specificity_profile->num_reads_score_sensitive_first_hit;
  const uint64_t num_reads_score_sensitive_event_hit = specificity_profile->num_reads_score_sensitive_event_hit;
  const uint64_t num_reads_score_sensitive_edit_hits = specificity_profile->num_reads_score_sensitive_edit_hits;
  // Maps Noise
  const uint64_t num_noise_maps_ties = specificity_profile->num_noise_maps_ties;
  const uint64_t num_noise_maps_misscored = specificity_profile->num_noise_maps_misscored;
  const uint64_t num_noise_maps_missing = specificity_profile->num_noise_maps_missing;
  // Print
  fprintf(stderr,"[GEM]>Specificity.Profile\n");
  /*
   * READS
   */
  fprintf(stderr,"  => Num.Reads                         %lu\n",num_reads);
  fprintf(stderr,"    => Num.Read.TP (Hit)  \n");
  fprintf(stderr,"      => Num.Reads.First.TP            %lu (%2.4f%%) (TP is first among reported maps)\n",num_reads_hit_first,100.0*(double)num_reads_hit_first/(double)num_reads);
  fprintf(stderr,"        => Num.Reads.First.TP.Tie      %lu (%2.4f%%) (TP is first though being a tie)\n",num_reads_hit_first_ties,100.0*(double)num_reads_hit_first_ties/(double)num_reads);
  fprintf(stderr,"      => Num.Reads.Any.TP              %lu (%2.4f%%) (TP is among reported maps)\n",num_reads_hit_any,100.0*(double)num_reads_hit_any/(double)num_reads);
  fprintf(stderr,"      => Num.Mapped                    %lu (%2.4f%%) \n",num_mapped,100.0*(double)num_mapped/(double)num_reads);
  fprintf(stderr,"    => Num.Read.FP (Fail)  \n");
  fprintf(stderr,"      => Num.Unmaped            %lu (%2.4f%%)\n",num_reads_miss_unmapped,100.0*(double)num_reads_miss_unmapped/(double)num_reads);
  fprintf(stderr,"      => Num.Missing            %lu (%2.4f%%)\n",num_reads_miss_missing,100.0*(double)num_reads_miss_missing/(double)num_reads);
  fprintf(stderr,"      => Num.Ties               %lu (%2.4f%%)\n",num_reads_miss_ties,100.0*(double)num_reads_miss_ties/(double)num_reads);
  fprintf(stderr,"        => Num.2Ties                 %lu (%2.4f%%)\n",num_reads_miss_2tie,100.0*(double)num_reads_miss_2tie/(double)num_reads_miss_ties);
  fprintf(stderr,"        => Num.3Ties                 %lu (%2.4f%%)\n",num_reads_miss_3tie,100.0*(double)num_reads_miss_3tie/(double)num_reads_miss_ties);
  fprintf(stderr,"        => Num.4Ties                 %lu (%2.4f%%)\n",num_reads_miss_4tie,100.0*(double)num_reads_miss_4tie/(double)num_reads_miss_ties);
  fprintf(stderr,"        => Num.>4Ties                %lu (%2.4f%%)\n",num_reads_miss_moretie,100.0*(double)num_reads_miss_moretie/(double)num_reads_miss_ties);
  fprintf(stderr,"      => Num.MisScored          %lu (%2.4f%%) (not-tie & subdominant)\n",num_reads_miss_misscored,100.0*(double)num_reads_miss_misscored/(double)num_reads);
  fprintf(stderr,"    => Num.Read.Sensitive.To.Score         %lu (%2.4f%%) (MMap, TP found & not-tie)\n",
      num_reads_score_sensitive,100.0*(double)num_reads_score_sensitive/(double)num_reads);
  fprintf(stderr,"      => Score.First.Hit                   %lu (%2.4f%%) (TP is first)\n",
      num_reads_score_sensitive_first_hit,100.0*(double)num_reads_score_sensitive_first_hit/(double)num_reads);
  fprintf(stderr,"      => Score.Event.Subdominant           %lu (%2.4f%%) (TP is subdominant with better event-distance)\n",
      num_reads_score_sensitive_event_hit,100.0*(double)num_reads_score_sensitive_event_hit/(double)num_reads);
  fprintf(stderr,"      => Score.Edit.Subdominant            %lu (%2.4f%%) (TP is subdominant with better edit-distance)\n",
      num_reads_score_sensitive_edit_hits,100.0*(double)num_reads_score_sensitive_edit_hits/(double)num_reads);
  /*
   * MAPS
   */
  fprintf(stderr,"  => Num.Maps                           %lu\n",num_maps);
  fprintf(stderr,"    => Num.Maps.TP                      %lu (%2.4f%%)\n",
      num_signal_maps,100.0*(double)num_signal_maps/(double)num_maps);
  fprintf(stderr,"    => Num.Maps.FP                      %lu (%2.4f%%)\n",
      num_noise_maps,100.0*(double)num_noise_maps/(double)num_maps);
  fprintf(stderr,"      => Num.Maps.Missing                   %lu (%2.4f%%)  (Surplus caused by TP missing)\n",
      num_noise_maps_missing,100.0*(double)num_noise_maps_missing/(double)num_maps);
  fprintf(stderr,"      => Num.Maps.Ties                      %lu (%2.4f%%)  (Surplus caused by TP tie)\n",
      num_noise_maps_ties,100.0*(double)num_noise_maps_ties/(double)num_maps);
  fprintf(stderr,"      => Num.Maps.MisScored                 %lu (%2.4f%%)  (Surplus caused by TP not-tie & subdominant)\n",
      num_noise_maps_misscored,100.0*(double)num_noise_maps_misscored/(double)num_maps);
#define GT_PRINT_CLASS_ACCURACY(RANGE) \
    specificity_profile->matches_classes_true[RANGE],specificity_profile->matches_classes_false[RANGE], \
    100.0*(double)specificity_profile->matches_classes_false[RANGE]/ \
     ((double)specificity_profile->matches_classes_true[RANGE]+(double)specificity_profile->matches_classes_false[RANGE])
  fprintf(stderr,"  => Classes (TP/FP)\n");
  fprintf(stderr,"    [Unique]      => %lu,%lu \t(FPR = %2.4f %%)\n",GT_PRINT_CLASS_ACCURACY(GT_CLASSES_UNIQUE));
  fprintf(stderr,"    [MMap-D>5]    => %lu,%lu \t(FPR = %2.4f %%)\n",GT_PRINT_CLASS_ACCURACY(GT_CLASSES_MMAP));
  fprintf(stderr,"    [MMap-D5]     => %lu,%lu \t(FPR = %2.4f %%)\n",GT_PRINT_CLASS_ACCURACY(GT_CLASSES_MMAP_D5));
  fprintf(stderr,"    [MMap-D4]     => %lu,%lu \t(FPR = %2.4f %%)\n",GT_PRINT_CLASS_ACCURACY(GT_CLASSES_MMAP_D4));
  fprintf(stderr,"    [MMap-D3]     => %lu,%lu \t(FPR = %2.4f %%)\n",GT_PRINT_CLASS_ACCURACY(GT_CLASSES_MMAP_D3));
  fprintf(stderr,"    [MMap-D2]     => %lu,%lu \t(FPR = %2.4f %%)\n",GT_PRINT_CLASS_ACCURACY(GT_CLASSES_MMAP_D2));
  fprintf(stderr,"    [MMap-D1]     => %lu,%lu \t(FPR = %2.4f %%)\n",GT_PRINT_CLASS_ACCURACY(GT_CLASSES_MMAP_D1));
  fprintf(stderr,"    [Tie]         => %lu,%lu \t(FPR = %2.4f %%)\n",GT_PRINT_CLASS_ACCURACY(GT_CLASSES_TIE));
  fprintf(stderr,"    [Tie.Perfect] => %lu,%lu \t(FPR = %2.4f %%)\n",GT_PRINT_CLASS_ACCURACY(GT_CLASSES_TIE_PERFECT));
  fprintf(stderr,"  => ROC (TP/FP)\n");
  int64_t p;
  int64_t acc_true=0, acc_false=0;
  for (p=UINT8_MAX-1;p>=0;--p) {
    acc_true+=specificity_profile->mapq_true[p];
    acc_false+=specificity_profile->mapq_false[p];
    //fprintf(stderr,"%lu,%lu,%lu\n",p,acc_true,acc_false);
    fprintf(stderr,"%lu,%lu,%lu\n",p,specificity_profile->mapq_true[p],specificity_profile->mapq_false[p]);
  }
}
GT_INLINE void gt_mapset_specificity_profile_print_fail(
    gt_buffered_output_file* const buffered_output,const gt_file_format file_format,
    gt_template* const template_master,gt_template* const template_slave,
    gt_output_map_attributes* const output_map_attributes) {
  // Tag
  gt_output_map_bofprint_tag(buffered_output,gt_template_get_string_tag(template_slave),template_slave->attributes,output_map_attributes);
  // Counters
  gt_bofprintf(buffered_output,"\t");
  gt_attributes* const attributes = gt_template_get_num_blocks(template_slave)==2 ?
      template_slave->attributes : template_slave->alignment_end1->attributes;
  gt_output_map_bofprint_counters(buffered_output,gt_template_get_counters_vector(template_slave),attributes,output_map_attributes);
  // Maps
  gt_bofprintf(buffered_output,"\t");
  gt_output_map_bofprint_template_maps(buffered_output,template_master,output_map_attributes);
  gt_bofprintf(buffered_output,"\t");
  uint64_t pos = 0;
  GT_TEMPLATE_ITERATE_MMAP__ATTR(template_slave,mmap_slave,mmap_slave_attr) {
//    const uint8_t mapq_score = (__mmap_slave_num_blocks==1) ?
//        (file_format==MAP ? mmap_slave[0]->gt_score : mmap_slave[0]->phred_score) :
//        (file_format==MAP ? mmap_slave_attr->gt_score : mmap_slave_attr->phred_score);
    if (pos++>0) gt_bofprintf(buffered_output,",");
    GT_MMAP_ITERATE(mmap_slave,map,end_position) {
      if (end_position>0) gt_bofprintf(buffered_output,":::");
      gt_output_map_bofprint_map(buffered_output,map,output_map_attributes);
    }
    if (__mmap_slave_num_blocks==2) {
      gt_bofprintf(buffered_output,":::%"PRIu64"(%ld)",mmap_slave_attr->phred_score,
          gt_map_get_observed_template_size(mmap_slave[0],mmap_slave[1]));
    }
//    else {
//      gt_bofprintf(buffered_output,":::%",);
//    } // TODO Quality average
  }
  gt_bofprintf(buffered_output,"\n");
}
GT_INLINE uint64_t gt_mapset_specificity_profile_compute_num_ties(
    const gt_file_format file_format,gt_template* const template_slave,gt_map** const mmap_match) {
  if (mmap_match==NULL) return 0;
  uint64_t ties_found = 0;
  GT_TEMPLATE_ITERATE_MMAP__ATTR(template_slave,mmap_slave,mmap_slave_attr) {
//    const uint8_t mapq_score = (__mmap_slave_num_blocks==1) ?
//        (file_format==MAP ? mmap_slave[0]->gt_score : mmap_slave[0]->phred_score) :
//        (file_format==MAP ? mmap_slave_attr->gt_score : mmap_slave_attr->phred_score);
    // if (mapq_score < parameters.mapq_min || mapq_score > parameters.mapq_max) continue;
    if (__mmap_slave_num_blocks==1) {
      if (gt_map_cmp_cigar(mmap_slave[0],mmap_match[0])==0) ++ties_found;
    } else {
      if (gt_map_cmp_cigar(mmap_slave[0],mmap_match[0])==0 &&
          gt_map_cmp_cigar(mmap_slave[1],mmap_match[1])==0 /* TODO */) ++ties_found;
    }
  }
  return ties_found;
}
GT_INLINE uint64_t gt_specificity_profile_compute_class(gt_template* const template) {
  const uint64_t num_mmaps = gt_template_get_num_mmaps(template);
  if (num_mmaps==0) {
    return GT_CLASSES_UNMAPPED;
  } else if (num_mmaps==1) {
    return GT_CLASSES_UNIQUE;
  } else {
    uint64_t delta, edit_primary, edit_secondary;
    if (gt_template_get_num_blocks(template)==2) {
      gt_mmap* const primary = gt_template_get_mmap(template,0);
      gt_mmap* const secondary = gt_template_get_mmap(template,1);
      edit_primary =
          gt_map_get_levenshtein_distance(primary->mmap[0]) +
          gt_map_get_levenshtein_distance(primary->mmap[1]);
      edit_secondary =
          gt_map_get_levenshtein_distance(secondary->mmap[0]) +
          gt_map_get_levenshtein_distance(secondary->mmap[1]);
    } else {
      gt_alignment* const alignment = gt_template_get_block(template,0);
      gt_map* const primary = gt_alignment_get_map(alignment,0);
      gt_map* const secondary = gt_alignment_get_map(alignment,1);
      edit_primary = gt_map_get_levenshtein_distance(primary);
      edit_secondary = gt_map_get_levenshtein_distance(secondary);
    }
    if (edit_primary > edit_secondary) {
      return GT_CLASSES_TIE;
    } else {
      delta = edit_secondary - edit_primary;
      if (delta==0) {
        if (gt_template_get_num_blocks(template)==2) {
          gt_mmap* const primary = gt_template_get_mmap(template,0);
          gt_mmap* const secondary = gt_template_get_mmap(template,1);
          if (gt_map_cmp_cigar(primary->mmap[0],secondary->mmap[0])==0 &&
              gt_map_cmp_cigar(primary->mmap[1],secondary->mmap[1])==0) {
            return GT_CLASSES_TIE_PERFECT;
          } else {
            return GT_CLASSES_TIE;
          }
        } else {
          gt_alignment* const alignment = gt_template_get_block(template,0);
          gt_map* const primary = gt_alignment_get_map(alignment,0);
          gt_map* const secondary = gt_alignment_get_map(alignment,1);
          if (gt_map_cmp_cigar(primary,secondary)==0) {
            return GT_CLASSES_TIE_PERFECT;
          } else {
            return GT_CLASSES_TIE;
          }
        }
      } else if (delta==1) {
        return GT_CLASSES_MMAP_D1;
      } else if (delta==2) {
        return GT_CLASSES_MMAP_D2;
      } else if (delta==3) {
        return GT_CLASSES_MMAP_D3;
      } else if (delta==4) {
        return GT_CLASSES_MMAP_D4;
      } else if (delta==5) {
        return GT_CLASSES_MMAP_D5;
      } else {
        return GT_CLASSES_MMAP;
      }
    }
  }
}
GT_INLINE void gt_specificity_profile_account(
    specificity_profile_t* const specificity_profile,const gt_file_format file_format,
    gt_template* const template_master,gt_template* const template_slave,
    gt_buffered_output_file* const buffered_output_ties,gt_buffered_output_file* const buffered_output_missing,
    gt_buffered_output_file* const buffered_output_unmapped,gt_buffered_output_file* const buffered_output_misscored,
    gt_buffered_output_file* const buffered_output_first_tp,gt_buffered_output_file* const buffered_output_sub_tp,
    gt_output_map_attributes* const output_map_attributes) {
  // Compute match-class
  const uint64_t match_class = gt_specificity_profile_compute_class(template_slave);
  // Calculate specificity
  bool first_tp = false, sub_tp = false, any_tp = false;
  uint64_t least_map_distance = UINT64_MAX, least_map_edit_distance = UINT64_MAX;
  uint64_t position = 0, best_pos = 0, maps_skipped = 0;
  uint8_t max_mapq_score = 0;
  GT_TEMPLATE_ITERATE_MMAP__ATTR(template_slave,mmap_slave,mmap_slave_attr) {
    uint8_t mapq_score = (__mmap_slave_num_blocks==1) ?
        (file_format==MAP ? mmap_slave[0]->gt_score : mmap_slave[0]->phred_score) :
        (file_format==MAP ? mmap_slave_attr->gt_score : mmap_slave_attr->phred_score);
    if (__mmap_slave_num_blocks==2 && mmap_slave[0]->phred_score!=mmap_slave[1]->phred_score) {
      mapq_score = GT_MIN(mmap_slave[0]->phred_score,mmap_slave[1]->phred_score);
    }
    if (position==0 && (mapq_score < parameters.mapq_min || mapq_score > parameters.mapq_max)) {
      maps_skipped = gt_template_get_num_mmaps(template_slave);
      break;
    }
    const bool contained = gt_template_is_mmap_contained_fx(gt_mapset_mmap_cmp,template_master,mmap_slave);
    // Display First TP/FP (0/1)
    if (position==0) {
      // gt_bofprintf(buffered_output_tp,"%s\t",template_master->tag->buffer);
      if (contained) {
        gt_bofprintf(buffered_output_first_tp,"1\n");
      } else {
        gt_bofprintf(buffered_output_first_tp,"0\n");
      }
    }
    if (contained && any_tp) continue; // Skip TP duplicates due to repetitions or ambiguity with eq-th
    max_mapq_score = (mapq_score > max_mapq_score) ? mapq_score : max_mapq_score;
    const uint64_t mmap_edit_distance = gt_mmap_get_global_levenshtein_distance(mmap_slave,__mmap_slave_num_blocks);
    least_map_edit_distance = (mmap_edit_distance < least_map_edit_distance) ? mmap_edit_distance : least_map_edit_distance;
    const uint64_t mmap_distance = gt_mmap_get_global_distance(mmap_slave,__mmap_slave_num_blocks);
    least_map_distance = (mmap_distance < least_map_distance) ? mmap_distance : least_map_distance;
    if (contained) {
      best_pos = position;
      if (position==0) {
        first_tp = true;
      } else {
        sub_tp = true;
      }
      any_tp = true;
      specificity_profile->mapq_true[mapq_score]++;
      if (position==0) specificity_profile->matches_classes_true[match_class]++;
      ++specificity_profile->num_signal_maps;
    } else {
      specificity_profile->mapq_false[mapq_score]++;
      if (position==0) specificity_profile->matches_classes_false[match_class]++;
      ++specificity_profile->num_noise_maps;
    }
    ++position;
    ++specificity_profile->num_maps;
  }
  // Compute ties & other metrics
  const uint64_t num_blocks = gt_template_get_num_blocks(template_slave);
  const uint64_t template_slave_num_maps = gt_template_get_num_mmaps(template_slave) - maps_skipped;
  uint64_t ties_found = 0;
  if (num_blocks==1) {
    gt_map* mmap_slave_first = (gt_alignment_get_num_maps(template_slave->alignment_end1) > 0) ?
        gt_alignment_get_map(template_slave->alignment_end1,0) : NULL;
    ties_found = gt_mapset_specificity_profile_compute_num_ties(file_format,template_slave,&mmap_slave_first);
  } else {
    gt_map** mmap_slave_first = (gt_template_get_num_mmaps(template_slave) > 0) ?
        gt_template_get_mmap(template_slave,0)->mmap : NULL;
    ties_found = gt_mapset_specificity_profile_compute_num_ties(file_format,template_slave,mmap_slave_first);
  }
  /*
   * READ level
   */
  ++specificity_profile->num_reads;
  if (template_slave_num_maps > 0) ++specificity_profile->num_mapped;
  if (first_tp) {
    ++specificity_profile->num_reads_hit_first;
    if (ties_found > 1) {
      ++specificity_profile->num_reads_hit_first_ties;
    }
  } else {
    if (ties_found > 1) { // Identify ties
      ++specificity_profile->num_reads_miss_ties;
      if (ties_found==2) ++specificity_profile->num_reads_miss_2tie;
      else if (ties_found==3) ++specificity_profile->num_reads_miss_3tie;
      else if (ties_found==4) ++specificity_profile->num_reads_miss_4tie;
      else ++specificity_profile->num_reads_miss_moretie;
      gt_mapset_specificity_profile_print_fail(
          buffered_output_ties,file_format,template_master,template_slave,output_map_attributes);
    } else {
      if (any_tp) {
        ++specificity_profile->num_reads_miss_misscored;
        gt_mapset_specificity_profile_print_fail(
            buffered_output_misscored,file_format,template_master,template_slave,output_map_attributes);
      } else {
        if (template_slave_num_maps > 0) {
          ++specificity_profile->num_reads_miss_missing;
          gt_mapset_specificity_profile_print_fail(
              buffered_output_missing,file_format,template_master,template_slave,output_map_attributes);
        } else {
          ++specificity_profile->num_reads_miss_unmapped;
          gt_mapset_specificity_profile_print_fail(
              buffered_output_unmapped,file_format,template_master,template_slave,output_map_attributes);
        }
      }
    }
  }
  // Display SUB TP/FP (0/1)
  if (sub_tp) {
    gt_mapset_specificity_profile_print_fail(
        buffered_output_sub_tp,file_format,template_master,template_slave,output_map_attributes);
  }
  /*
   * MAP level
   */
  if (any_tp) {
    ++specificity_profile->num_reads_hit_any;
    if (ties_found > 1) {
      specificity_profile->num_noise_maps_ties += template_slave_num_maps - 1;
    } else if (!first_tp) {
      specificity_profile->num_noise_maps_misscored += template_slave_num_maps - 1;
    }
  } else {
    if (template_slave_num_maps > 0) {
      specificity_profile->num_noise_maps_missing += template_slave_num_maps;
    }
  }
  /*
   * Score eligible
   */
  const bool is_score_elegible = any_tp && ties_found <= 1 && template_slave_num_maps > 1;
  if (is_score_elegible) {
    ++specificity_profile->num_reads_score_sensitive;
    if (num_blocks==1) {
      gt_map* const best_map = gt_alignment_get_map(template_slave->alignment_end1,0);
      if (best_pos == 0) {
        ++specificity_profile->num_reads_score_sensitive_first_hit;
      } else {
        if (least_map_distance < gt_map_get_distance(best_map)) {
          ++specificity_profile->num_reads_score_sensitive_event_hit;
        }
        if (least_map_edit_distance < gt_map_get_levenshtein_distance(best_map)) {
          ++specificity_profile->num_reads_score_sensitive_edit_hits;
        }
      }
    } else {
      // TODO
    }
  }
}
void gt_mapset_specificity_profile() {
  // File IN/OUT
  gt_input_file* input_file_1 = gt_input_file_open(parameters.name_input_file_1,parameters.mmap_input);
  gt_input_file* input_file_2 = (parameters.name_input_file_2==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file_2,parameters.mmap_input);
  if (parameters.name_input_file_2==NULL) GT_SWAP(input_file_1,input_file_2);

  // Output files
  char* const output_file_name_prefix = (parameters.name_output_file==NULL) ? "noname" : parameters.name_output_file;
  char output_file_name[1000];
  sprintf(output_file_name,"sp.%s.ties.rmap",output_file_name_prefix);
  gt_output_file* const output_file_ties = gt_output_file_new(output_file_name,SORTED_FILE);
  sprintf(output_file_name,"sp.%s.missing.rmap",output_file_name_prefix);
  gt_output_file* const output_file_missing = gt_output_file_new(output_file_name,SORTED_FILE);
  sprintf(output_file_name,"sp.%s.unmapped.rmap",output_file_name_prefix);
  gt_output_file* const output_file_unmapped = gt_output_file_new(output_file_name,SORTED_FILE);
  sprintf(output_file_name,"sp.%s.misscored.rmap",output_file_name_prefix);
  gt_output_file* const output_file_misscored = gt_output_file_new(output_file_name,SORTED_FILE);
  sprintf(output_file_name,"sp.%s.tp.rmap",output_file_name_prefix);
  gt_output_file* const output_file_first_tp = gt_output_file_new(output_file_name,SORTED_FILE);
  sprintf(output_file_name,"sp.%s.subtp.rmap",output_file_name_prefix);
  gt_output_file* const output_file_sub_tp = gt_output_file_new(output_file_name,SORTED_FILE);

  // Buffered I/O
  gt_buffered_input_file* buffered_input_1 = gt_buffered_input_file_new(input_file_1);
  gt_buffered_input_file* buffered_input_2 = gt_buffered_input_file_new(input_file_2);
  gt_buffered_output_file* buffered_output_ties = gt_buffered_output_file_new(output_file_ties);
  gt_buffered_output_file* buffered_output_missing = gt_buffered_output_file_new(output_file_missing);
  gt_buffered_output_file* buffered_output_unmapped = gt_buffered_output_file_new(output_file_unmapped);
  gt_buffered_output_file* buffered_output_misscored = gt_buffered_output_file_new(output_file_misscored);
  gt_buffered_output_file* buffered_output_first_tp = gt_buffered_output_file_new(output_file_first_tp);
  gt_buffered_output_file* buffered_output_sub_tp = gt_buffered_output_file_new(output_file_sub_tp);
  gt_buffered_input_file_attach_buffered_output(buffered_input_1,buffered_output_ties);
  gt_buffered_input_file_attach_buffered_output(buffered_input_1,buffered_output_missing);
  gt_buffered_input_file_attach_buffered_output(buffered_input_1,buffered_output_unmapped);
  gt_buffered_input_file_attach_buffered_output(buffered_input_1,buffered_output_misscored);
  gt_buffered_input_file_attach_buffered_output(buffered_input_1,buffered_output_first_tp);
  gt_buffered_input_file_attach_buffered_output(buffered_input_1,buffered_output_sub_tp);

  // Signal/Noise
  specificity_profile_t* specificity_profile = gt_specificity_profile_new();

  // Template I/O (synch)
  gt_template *template_master = gt_template_new();
  gt_template *template_slave = gt_template_new();
  gt_output_map_attributes* output_map_attributes = gt_output_map_attributes_new();
  while (gt_mapset_read_template_get_commom_map(buffered_input_1,buffered_input_2,template_master,template_slave)) {
//    if (strcmp("Sim.Illumina.l75.0000000005",template_master->tag->buffer)==0) {
//      printf("Hola\n");
//    }
    // Record current read length
    current_read_length = gt_template_get_total_length(template_master)/gt_template_get_num_blocks(template_master);
    // Record Stats
    gt_specificity_profile_account(
        specificity_profile,input_file_2->file_format,template_master,template_slave,
        buffered_output_ties,buffered_output_missing,buffered_output_unmapped,
        buffered_output_misscored,buffered_output_first_tp,
        buffered_output_sub_tp,output_map_attributes);
  }

  // Print out results
  gt_specificity_profile_print(specificity_profile);

  // Clean
  gt_specificity_profile_delete(specificity_profile);
  gt_template_delete(template_master);
  gt_template_delete(template_slave);
  gt_buffered_input_file_close(buffered_input_1);
  gt_buffered_input_file_close(buffered_input_2);
  gt_buffered_output_file_close(buffered_output_missing);
  gt_buffered_output_file_close(buffered_output_unmapped);
  gt_buffered_output_file_close(buffered_output_misscored);
  gt_buffered_output_file_close(buffered_output_ties);
  gt_buffered_output_file_close(buffered_output_first_tp);
  gt_buffered_output_file_close(buffered_output_sub_tp);
  gt_input_file_close(input_file_1);
  gt_input_file_close(input_file_2);
  gt_output_file_close(output_file_missing);
  gt_output_file_close(output_file_unmapped);
  gt_output_file_close(output_file_misscored);
  gt_output_file_close(output_file_ties);
  gt_output_file_close(output_file_first_tp);
  gt_output_file_close(output_file_sub_tp);
}
/*
 * Merge-Map
 */
void gt_mapset_perform_merge_map() {
  // Open file IN/OUT
  gt_input_file* input_file_1 = gt_input_file_open(parameters.name_input_file_1,parameters.mmap_input);
  gt_input_file* input_file_2 = (parameters.name_input_file_2==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file_2,parameters.mmap_input);
  if (parameters.name_input_file_2==NULL) GT_SWAP(input_file_1,input_file_2);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);

  // Mutex
  pthread_mutex_t input_mutex = PTHREAD_MUTEX_INITIALIZER;

  // Parallel reading+process
#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(parameters.num_threads)
#endif
  {
    if (parameters.files_contain_same_reads) {
      gt_merge_synch_map_files(&input_mutex,parameters.paired_end,output_file,input_file_1,input_file_2);
    } else {
      gt_merge_unsynch_map_files(&input_mutex,input_file_1,input_file_2,parameters.paired_end,output_file);
    }
  }

  // Clean
  gt_input_file_close(input_file_1);
  gt_input_file_close(input_file_2);
  gt_output_file_close(output_file);
}
/*
 * Display Compact
 */
void gt_mapset_display_compact_map() {
  // Open file IN/OUT
  gt_input_file* input_file = (parameters.name_input_file_1==NULL) ?
      gt_input_stream_open(stdin) : gt_input_file_open(parameters.name_input_file_1,parameters.mmap_input);
  gt_output_file* output_file = (parameters.name_output_file==NULL) ?
      gt_output_stream_new(stdout,SORTED_FILE) : gt_output_file_new(parameters.name_output_file,SORTED_FILE);
#ifdef HAVE_OPENMP
  #pragma omp parallel num_threads(parameters.num_threads)
#endif
  {
    gt_output_map_attributes* const output_map_attributes = gt_output_map_attributes_new();
    output_map_attributes->compact = true;
    GT_BEGIN_READING_WRITING_LOOP(input_file,output_file,parameters.paired_end,buffered_output,template) {
      GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
        // Print compact summary
        gt_bofprintf(buffered_output,"End1::"PRIgts"[%"PRIu64"]\t",
            PRIgts_content(alignment->tag),gt_string_get_length(alignment->read));
        gt_output_map_bofprint_counters(buffered_output,alignment->counters,alignment->attributes,output_map_attributes);
        gt_bofprintf(buffered_output,"\t");
        uint64_t printed = 0;
        GT_ALIGNMENT_ITERATE(alignment,map) {
          if (printed>0) {
            gt_bofprintf(buffered_output,","PRIgts,PRIgts_content(map->seq_name));
          } else {
            gt_bofprintf(buffered_output,PRIgts,PRIgts_content(map->seq_name));
          }
          ++printed;
        }
        gt_bofprintf(buffered_output,"\n");
      }
    } GT_END_READING_WRITING_LOOP(input_file,output_file,template);
    // Clean
    gt_output_map_attributes_delete(output_map_attributes);
  }

  // Clean
  gt_input_file_close(input_file);
  gt_output_file_close(output_file);
}
/*
 * Main Parse
 */
void gt_filter_parse_operation(char* const string_operation) {
  if (gt_streq(string_operation,"INTERSECCTION") || gt_streq(string_operation,"Intersection") || gt_streq(string_operation,"intersection")) {
    parameters.operation = GT_MAP_SET_INTERSECTION;
  } else if (gt_streq(string_operation,"UNION") || gt_streq(string_operation,"Union") || gt_streq(string_operation,"union")) {
    parameters.operation = GT_MAP_SET_UNION;
  } else if (gt_streq(string_operation,"DIFFERENCE") || gt_streq(string_operation,"Difference") || gt_streq(string_operation,"difference")) {
    parameters.operation = GT_MAP_SET_DIFFERENCE;
  } else if (gt_streq(string_operation,"COMPARE") || gt_streq(string_operation,"Compare") || gt_streq(string_operation,"compare")) {
    parameters.operation = GT_MAP_SET_COMPARE;
  } else if (gt_streq(string_operation,"JOIN") || gt_streq(string_operation,"Join") || gt_streq(string_operation,"join")) {
    parameters.operation = GT_MAP_SET_JOIN;
  } else if (gt_streq(string_operation,"MERGE-MAP") || gt_streq(string_operation,"Merge-map") || gt_streq(string_operation,"merge-map")) {
    parameters.operation = GT_MERGE_MAP;
  } else if (gt_streq(string_operation,"DISPLAY-COMPACT") || gt_streq(string_operation,"Display-compact") || gt_streq(string_operation,"display-compact")) {
    parameters.operation = GT_DISPLAY_COMPACT_MAP;
  } else if (gt_streq(string_operation,"SPECIFICITY-PROFILE") || gt_streq(string_operation,"Specificity-profile") || gt_streq(string_operation,"specificity-profile")) {
    parameters.operation = GT_MAP_SET_SPECIFICITY_PROFILE;
  } else {
    if (string_operation[0]=='I' || string_operation[0]=='i') {
      fprintf(stderr,"\tAssuming 'Intersection' ...\n");
      parameters.operation = GT_MAP_SET_INTERSECTION;
    } else if (string_operation[0]=='U' || string_operation[0]=='u') {
      fprintf(stderr,"\tAssuming 'Union' ...\n");
      parameters.operation = GT_MAP_SET_UNION;
    } else if (string_operation[0]=='D' || string_operation[0]=='d') {
      fprintf(stderr,"\tAssuming 'Difference' ...\n");
      parameters.operation = GT_MAP_SET_DIFFERENCE;
    } else if (string_operation[0]=='C' || string_operation[0]=='c') {
      fprintf(stderr,"\tAssuming 'Compare' ...\n");
      parameters.operation = GT_MAP_SET_COMPARE;
    } else if (string_operation[0]=='P' || string_operation[0]=='p') {
      fprintf(stderr,"\tAssuming 'Join' ...\n");
      parameters.operation = GT_MAP_SET_JOIN;
    } else if (string_operation[0]=='M' || string_operation[0]=='m') {
      fprintf(stderr,"\tAssuming 'Merge-map' ...\n");
      parameters.operation = GT_MERGE_MAP;
    } else {
      gt_fatal_error_msg("Unknown operation '%s' in {"GT_MAPSET_OPERATIONS"}",string_operation);
    }
  }
}
void parse_arguments(int argc,char** argv) {
  struct option* gt_mapset_getopt = gt_options_adaptor_getopt(gt_mapset_options);
  gt_string* const gt_mapset_short_getopt = gt_options_adaptor_getopt_short(gt_mapset_options);
  int option, option_index;
  while (true) {
    // Get option & Select case
    if ((option=getopt_long(argc,argv,
        gt_string_get_string(gt_mapset_short_getopt),gt_mapset_getopt,&option_index))==-1) break;
    switch (option) {
    /* Operations */
    case 'O':
      gt_filter_parse_operation(optarg);
      break;
    /* I/O */
    case 300:
      parameters.name_input_file_1 = optarg;
      break;
    case 301:
      parameters.name_input_file_2 = optarg;
      break;
    case 'p':
      parameters.paired_end = true;
      break;
    case 302:
      parameters.mmap_input = true;
      gt_fatal_error(NOT_IMPLEMENTED);
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    /* Compare Function */
    case 's': // files-with-same-reads
      parameters.files_contain_same_reads = true;
      break;
    case 400: // eq-th
      parameters.eq_threshold = atof(optarg);
      break;
    case 401: // strict
      parameters.strict = true;
      break;
    case 'q':
      parameters.mapq_min = atoi(optarg);
      break;
    case 'Q':
      parameters.mapq_max = atoi(optarg);
      break;
    /* Misc */
    case 'v':
      parameters.verbose = true;
      break;
    case 't':
#ifdef HAVE_OPENMP
      parameters.num_threads = atol(optarg);
#endif
      break;
    case 'h':
      fprintf(stderr, "USE: ./gt.mapset [OPERATION] [ARGS]...\n");
      gt_options_fprint_menu(stderr,gt_mapset_options,gt_mapset_groups,true,false);
      exit(1);
    case 'J':
      gt_options_fprint_json_menu(stderr,gt_mapset_options,gt_mapset_groups,true,false);
      exit(1);
      break;
    case '?':
    default:
      gt_fatal_error_msg("Option not recognized");
    }
  }
  // Check parameters
  if (parameters.operation==GT_MAP_SET_UNKNOWN) {
    gt_fatal_error_msg("Please specify operation {"GT_MAPSET_OPERATIONS"}");
  }
  if (parameters.operation!=GT_DISPLAY_COMPACT_MAP && !parameters.name_input_file_1) {
    gt_fatal_error_msg("Input file 1 required (--i1)\n");
  }
  // Free
  gt_string_delete(gt_mapset_short_getopt);
}
/*
 * Main
 */
int main(int argc,char** argv) {
  // GT error handler
  gt_handle_error_signals();
  // Parsing command-line options
  parse_arguments(argc,argv);
  // Do it !
  if (parameters.operation==GT_MERGE_MAP) {
    gt_mapset_perform_merge_map();
  } else if (parameters.operation==GT_DISPLAY_COMPACT_MAP) {
    gt_mapset_display_compact_map();
  } else if (parameters.operation==GT_MAP_SET_INTERSECTION ||
      parameters.operation==GT_MAP_SET_UNION ||
      parameters.operation==GT_MAP_SET_DIFFERENCE) {
    gt_mapset_perform_set_operations();
  } else if (parameters.operation==GT_MAP_SET_COMPARE) {
    gt_mapset_perform_cmp_operations();
  } else if (parameters.operation==GT_MAP_SET_SPECIFICITY_PROFILE) {
    gt_mapset_specificity_profile();
  } else {
    gt_fatal_error_msg("Invalid map-set operation");
  }
  return 0;
}


