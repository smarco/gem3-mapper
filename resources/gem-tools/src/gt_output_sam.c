/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_sam.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_output_sam.h"

/*
 * Constants
 */
#define GT_OUTPUT_SAM_FORMAT_VERSION "1.4"

/*
 * Output SAM Attributes
 */
/* Setup */
GT_INLINE gt_output_sam_attributes* gt_output_sam_attributes_new() {
  gt_output_sam_attributes* const attributes = gt_alloc(gt_output_sam_attributes);
  /* Optional fields */
  attributes->sam_attributes=NULL;
  attributes->attribute_func_params=NULL;
  /* Reset defaults */
  gt_output_sam_attributes_clear(attributes);
  return attributes;
}
GT_INLINE void gt_output_sam_attributes_delete(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  if (attributes->sam_attributes!=NULL) gt_sam_attributes_delete(attributes->sam_attributes);
  if (attributes->attribute_func_params!=NULL) gt_sam_attribute_func_params_delete(attributes->attribute_func_params);
  gt_free(attributes);
}
GT_INLINE void gt_output_sam_attributes_clear(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  /* Format */
  attributes->format = GT_SAM;
  /* Read/Qualities */
  attributes->always_output_read__qualities = true;
  attributes->qualities_offset = GT_QUALS_OFFSET_33;
  /* Maps */
  attributes->max_printable_maps = UINT64_MAX;
  attributes->compact_format = true;
  /* Mismatch/CIGAR string */
  attributes->print_mismatches = false;
  /* SAM Optional Fields */
  attributes->print_optional_fields = true;
  if (attributes->sam_attributes!=NULL) {
    gt_sam_attributes_clear(attributes->sam_attributes);
  } else {
    attributes->sam_attributes = gt_sam_attributes_new();
  }
  if (attributes->attribute_func_params!=NULL) {
    gt_sam_attribute_func_params_clear(attributes->attribute_func_params);
  } else {
    attributes->attribute_func_params = gt_sam_attribute_func_params_new();
  }
}

/* Format */
GT_INLINE void gt_output_sam_attributes_set_format(gt_output_sam_attributes* const attributes,gt_output_sam_format_t const format) {
  GT_NULL_CHECK(attributes);
  attributes->format = format;
}
/* Read/Qualities */
GT_INLINE void gt_output_sam_attributes_dump_read__qualities_once(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  attributes->always_output_read__qualities = false;
}
GT_INLINE void gt_output_sam_attributes_always_dump_read__qualities(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  attributes->always_output_read__qualities = true;
}
GT_INLINE void gt_output_sam_attributes_set_qualities_offset(gt_output_sam_attributes* const attributes,gt_qualities_offset_t const qualities_offset) {
  attributes->qualities_offset = GT_QUALS_OFFSET_33;
}
/* Maps */
GT_INLINE void gt_output_sam_attributes_set_max_printable_maps(gt_output_sam_attributes* const attributes,const uint64_t max_printable_maps) {
  GT_NULL_CHECK(attributes);
  attributes->max_printable_maps = max_printable_maps;
}
GT_INLINE void gt_output_sam_attributes_set_compact_format(gt_output_sam_attributes* const attributes,const bool compact_format) {
  GT_NULL_CHECK(attributes);
  attributes->compact_format = compact_format;
}
/* CIGAR/Mismatch string */
GT_INLINE void gt_output_sam_attributes_set_print_mismatches(gt_output_sam_attributes* const attributes,const bool print_mismatches) {
  GT_NULL_CHECK(attributes);
  attributes->print_mismatches = print_mismatches;
}
/* SAM Optional Fields */
GT_INLINE void gt_output_sam_attributes_set_print_optional_fields(gt_output_sam_attributes* const attributes,const bool print_optional_fields) {
  GT_NULL_CHECK(attributes);
  attributes->print_optional_fields = print_optional_fields;
}
GT_INLINE void gt_output_sam_attributes_set_reference_sequence_archive(gt_output_sam_attributes* const attributes,gt_sequence_archive* const reference_sequence_archive) {
  GT_NULL_CHECK(attributes);
  attributes->attribute_func_params->sequence_archive = reference_sequence_archive;
}
GT_INLINE gt_sam_attributes* gt_output_sam_attributes_get_sam_attributes(gt_output_sam_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  return attributes->sam_attributes;
}

/*
 * // TODO replace with sample
 * SAM Headers
 * @HD  VN:1.0  SO:unsorted
 * @HD  VN:1.3
 *
 * @SQ  SN:chr10   LN:135534747  AS:hg19_ncbi37  SP:human
 * @SQ  SN:chr10   LN:135534747
 * @SQ  SN:chr11   LN:135006516
 *
 * @RG  ID:NOID   PG:tmap  SM:NOSM
 * @RG  ID:0      PG:GEM   PL:ILLUMINA  SM:0
 *
 * @PG  ID:tmap  CL:map4 -f /home/user/references/hsapiens_v37.fa -r /home/user/miseq_S/H.Sapiens.1M.S.MiSeq.low.l250_1.fastq -i fastq -s /home/user/miseq_S/map41.H.Sapiens.1M.S.MiSeq.low.l250_1.fastq.sam   VN:3.0.1
 * @PG  ID:dvtgm PN:stampy     VN:1.0.17_(r1481)  CL:-g /home/devel/user/hsapiens_v37 -h /home/devel/user/hsapiens_v37 --bwaoptions=/home/user/hsapiens_v37.fa --substitutionrate=0.08 --maxbasequal 90 -M /home/user/miseq_S/H.Sapiens.1M.S.MiSeq.low.l250_1.fastq
 * @PG  ID:GEM   PN:gem-2-sam  VN:1.414
 *
 * @CO  TM:Fri, 30 Nov 2012 14:14:13 CET        WD:/home/user/benchmark/DEF/miseq_S/CMD      HN:cn38.bullx   UN:user
 * @CO  BWAVersion: 0.6.1-r104
 *
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS sam_headers
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_headers_sh,gt_sam_headers* const sam_headers);
GT_INLINE gt_status gt_output_sam_gprint_headers_sh(gt_generic_printer* const gprinter,gt_sam_headers* const sam_headers) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  // Print all @HD line (Header line)
  if (sam_headers==NULL || gt_string_is_null(sam_headers->header)) {
    gt_gprintf(gprinter,"@HD\tVN:"GT_OUTPUT_SAM_FORMAT_VERSION"\n");
  } else {
    gt_gprintf(gprinter,"@HD\t"PRIgts"\n",PRIgts_content(sam_headers->header));
  }
  // Print all @SQ lines (Reference sequence dictionary)
  if (sam_headers==NULL || sam_headers->sequence_archive!=NULL) {
    gt_output_sam_gprint_headers_sa(gprinter,sam_headers->sequence_archive);
  }
  // Print all @RG lines (Read group)
  if (sam_headers==NULL || gt_vector_is_empty(sam_headers->read_group)) {
    gt_gprintf(gprinter,"@RG\tID:0\tPG:GTools\tSM:0\n");
  } else {
    GT_VECTOR_ITERATE(sam_headers->read_group,rg_line,line_num,gt_string*) {
      gt_gprintf(gprinter,"@RG\t"PRIgts"\n",PRIgts_content(*rg_line));
    }
  }
  // Print all @PG lines (Program)
  if (sam_headers==NULL || gt_vector_is_empty(sam_headers->program)) {
    gt_gprintf(gprinter,"@PG\tID:GToolsLib\tPN:gt_output_sam\tVN:"GT_VERSION"\n");
  } else {
    GT_VECTOR_ITERATE(sam_headers->program,prog_line,line_num,gt_string*) {
      gt_gprintf(gprinter,"@PG\t"PRIgts"\n",PRIgts_content(*prog_line));
    }
  }
  // Print all @CO lines (Comments)
  if (sam_headers==NULL || gt_vector_is_empty(sam_headers->comments)) {
    // Print Current Date
    gt_gprintf(gprinter,"@CO\tTM:");
    time_t current_time=time(0);
    struct tm local_time;
    localtime_r(&current_time,&local_time);
    gt_gprintf(gprinter,"%4d/%d/%d::%02d:%02d:%02d CET",
        1900+local_time.tm_year,local_time.tm_mon+1,local_time.tm_mday,
        local_time.tm_hour,local_time.tm_min,local_time.tm_sec);
    // Print GT banner
    gt_gprintf(gprinter,"\tGTools v"GT_VERSION" "GT_GIT_URL"\n");
  } else {
    GT_VECTOR_ITERATE(sam_headers->comments,comment,line_num,gt_string*) {
      gt_gprintf(gprinter,"@CO\t"PRIgts"\n",PRIgts_content(*comment));
    }
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS sequence_archive
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_headers_sa,gt_sequence_archive* const sequence_archive);
GT_INLINE gt_status gt_output_sam_gprint_headers_sa(gt_generic_printer* const gprinter,gt_sequence_archive* const sequence_archive) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Print all @SQ lines (Reference sequence dictionary)
  //   @SQ  SN:chr10 LN:135534747 AS:hg19_ncbi37  SP:human
  gt_sequence_archive_iterator sequence_archive_it;
  gt_sequence_archive_new_iterator(sequence_archive,&sequence_archive_it);
  gt_segmented_sequence* seq;
  while ((seq=gt_sequence_archive_iterator_next(&sequence_archive_it))) {
    // @SQ  SN:chr10 LN:135534747
    gt_gprintf(gprinter,"@SQ\tSN:"PRIgts"\tLN:%"PRIu64"\n",PRIgts_content(seq->seq_name),seq->sequence_total_length);
  }
  return 0;
}
/*
 * SAM QNAME (Tag)
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_qname,gt_string* const tag);
GT_INLINE gt_status gt_output_sam_gprint_qname(gt_generic_printer* const gprinter,gt_string* const tag) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  // Print the plain tag (no read pair info, nor extra tag nor nothing)
  char* const tag_buffer = gt_string_get_string(tag);
  int i;
  for (i=0;i<gt_string_get_length(tag);++i) {
    if (tag_buffer[i]==SPACE) break;
  }
  gt_gprintf(gprinter,PRIgts,i,gt_string_get_string(tag));
  return 0;
}
/*
 * SAM Flag
 *
 * 0x1   (Bit 0)  => Template having multiple segments in sequencing
 *                   (The read was part of a pair during sequencing) [read paired]
 * 0x2   (Bit 1)  => Each segment properly aligned according to the aligner
 *                   (The read is mapped in a pair) [read mapped in proper pair]
 * 0x4   (Bit 2)  => Segment unmapped. The query sequence is unmapped [read unmapped]
 * 0x8   (Bit 3)  => Next segment in the template unmapped. The mate is unmapped [mate unmapped]
 * 0x10  (Bit 4)  => SEQ being reverse complemented. Strand of query (0=forward 1=reverse) [read reverse strand]
 * 0x20  (Bit 5)  => SEQ of the next segment in the template being reversed [mate reverse strand]
 * 0x40  (Bit 6)  => The first segment in the template [first in pair]
 * 0x80  (Bit 7)  => The last segment in the template [second in pair]
 * 0x100 (Bit 8)  => Secondary alignment [not primary alignment]
 * 0x200 (Bit 9)  => Not passing quality controls [read fails platform/vendor quality checks]
 * 0x400 (Bit 10) => PCR or optical duplicate [read is PCR or optical duplicate]
 *
 * - RULE1:: Bit 0x4 is the only reliable place to tell whether the segment is unmapped. If 0x4 is set,
 *     no assumptions can be made about RNAME, POS, CIGAR, MAPQ, bits 0x2, 0x10 and 0x100
 *     and the bit 0x20 of the next segment in the template.
 * - RULE2:: If 0x40 and 0x80 are both set, the segment is part of a linear template, but it is neither
 *     the first nor the last segment. If both 0x40 and 0x80 are unset, the index of the segment
 *     in the template is unknown. This may happen for a non-linear template or the index is
 *     lost in data processing.
 * - RULE3:: Bit 0x100 marks the alignment not to be used in certain analyses when the tools in use
 *     are aware of this bit.
 * - RULE4:: If 0x1 is unset, no assumptions can be made about 0x2, 0x8, 0x20, 0x40 and 0x80.
*/
GT_INLINE uint16_t gt_output_sam_calculate_flag_se_map(
    gt_map* const map,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate) {
  return (gt_expect_true(map!=NULL)) ?
      gt_output_sam_calculate_flag_se(true,map->strand,secondary_alignment,not_passing_QC,PCR_duplicate): // Mapped
      gt_output_sam_calculate_flag_se(false,FORWARD,secondary_alignment,not_passing_QC,PCR_duplicate); // Unmapped
}
GT_INLINE uint16_t gt_output_sam_calculate_flag_pe_map(
    gt_map* const map,gt_map* const mate,const bool is_map_first_in_pair,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate) {
  return gt_output_sam_calculate_flag_pe(
      map!=NULL && mate!=NULL, /* read_paired */
      map!=NULL,               /* read_mapped */
      mate!=NULL,              /* mate_strand */
      (map!=NULL) ? map->strand : FORWARD,   /* read_strand */
      (mate!=NULL) ? mate->strand : FORWARD, /* mate_strand */
      is_map_first_in_pair,    /* first_in_pair */
      !is_map_first_in_pair,   /* last_in_pair */
      secondary_alignment,not_passing_QC,PCR_duplicate);
}
GT_INLINE uint16_t gt_output_sam_calculate_flag_pe(
    const bool read_paired,const bool read_mapped,const bool mate_mapped,
    const gt_strand read_strand,const gt_strand mate_strand,
    const bool first_in_pair,const bool last_in_pair,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate) {
  /* 0x1 */
  uint16_t sam_flag = GT_SAM_FLAG_MULTIPLE_SEGMENTS;
  /* 0x4  */
  if (!read_mapped) { // (**RULE1)
    sam_flag |= GT_SAM_FLAG_UNMAPPED;
    /* 0x8  */
    if (!mate_mapped) {
      sam_flag |= GT_SAM_FLAG_NEXT_UNMAPPED;
    } else {
      /* 0x20 */if (mate_strand==REVERSE) sam_flag |= GT_SAM_FLAG_NEXT_REVERSE_COMPLEMENT;
    }
  } else {
    /* 0x10 */  if (read_strand==REVERSE) sam_flag |= GT_SAM_FLAG_REVERSE_COMPLEMENT;
    /* 0x100 */ if (secondary_alignment) sam_flag |= GT_SAM_FLAG_SECONDARY_ALIGNMENT;
    /* 0x8 */
    if (!mate_mapped) {
      sam_flag |= GT_SAM_FLAG_NEXT_UNMAPPED;
    } else {
      /* 0x2 Each segment properly aligned can only take place if both ends are mapped */
      if (read_paired) sam_flag |= GT_SAM_FLAG_PROPERLY_ALIGNED;
      /* 0x20 */
      if (mate_strand==REVERSE) sam_flag |= GT_SAM_FLAG_NEXT_REVERSE_COMPLEMENT;
    }
  }
  /* 0x40 */  if (first_in_pair)  sam_flag |= GT_SAM_FLAG_FIRST_SEGMENT;
  /* 0x80 */  if (last_in_pair)   sam_flag |= GT_SAM_FLAG_LAST_SEGMENT;
  /* 0x200 */ if (not_passing_QC) sam_flag |= GT_SAM_FLAG_NOT_PASSING_QC;
  /* 0x400 */ if (PCR_duplicate)  sam_flag |= GT_SAM_FLAG_PCR_OR_OPTICAL_DUPLICATE;
  return sam_flag;
}
GT_INLINE uint16_t gt_output_sam_calculate_flag_se(
    const bool read_mapped,const gt_strand read_strand,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate) {
  uint16_t sam_flag = 0; // (**RULE4)
  /* 0x4  */
  if (!read_mapped) { // (**RULE1)
    sam_flag |= GT_SAM_FLAG_UNMAPPED;
  } else {
    /* 0x10 */  if (read_strand==REVERSE) sam_flag |= GT_SAM_FLAG_REVERSE_COMPLEMENT;
    /* 0x100 */ if (secondary_alignment) sam_flag |= GT_SAM_FLAG_SECONDARY_ALIGNMENT;
  }
  /* 0x200 */
  if (not_passing_QC) sam_flag |= GT_SAM_FLAG_NOT_PASSING_QC;
  /* 0x400 */
  if (PCR_duplicate) sam_flag |= GT_SAM_FLAG_PCR_OR_OPTICAL_DUPLICATE;
  return sam_flag;
}
GT_INLINE uint16_t gt_output_sam_calculate_flag(
    const bool paired_end,const bool read_paired,
    const bool read_mapped,const bool mate_mapped,
    const gt_strand read_strand,const gt_strand mate_strand,
    const bool first_in_pair,const bool last_in_pair,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate) {
  return (paired_end) ? /* 0x1 */
      gt_output_sam_calculate_flag_pe(
          read_paired,read_mapped,mate_mapped,
          read_strand,mate_strand,first_in_pair,last_in_pair,
          secondary_alignment,not_passing_QC,PCR_duplicate): // PE
      gt_output_sam_calculate_flag_se(read_mapped,read_strand,
          secondary_alignment,not_passing_QC,PCR_duplicate); // SE
}
/*
 * SAM CIGAR
 */
#define GT_OUTPUT_SAM_CIGAR_FORWARD_MATCH() \
  if (misms_pos!=centinel) { \
    gt_gprintf(gprinter,"%"PRIu64"%c",misms_pos-centinel,(attributes->print_mismatches)?'=':'M'); \
    centinel = misms_pos; \
  }
#define GT_OUTPUT_SAM_CIGAR_REVERSE_MATCH() \
  if (misms_pos!=centinel) { \
    gt_gprintf(gprinter,"%"PRIu64"%c",centinel-misms_pos,(attributes->print_mismatches)?'=':'M'); \
    centinel = misms_pos; \
  }
GT_INLINE gt_status gt_output_sam_gprint_map_block_cigar_reverse(gt_generic_printer* const gprinter,gt_map* const map,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  // Map auxiliary variables
  const uint64_t map_length = gt_map_get_base_length(map);
  int64_t centinel = map_length;
  // Mismatches auxiliary variables
  const uint64_t num_misms = gt_map_get_num_misms(map);
  uint64_t misms_n = num_misms;
  gt_misms* misms;
  while (misms_n > 0) {
    misms = gt_map_get_misms(map,misms_n-1);
    const uint64_t misms_pos = gt_misms_get_position(misms);
    switch (misms->misms_type) {
      case MISMS:
        if (attributes->print_mismatches) {
          GT_OUTPUT_SAM_CIGAR_REVERSE_MATCH();
          gt_gprintf(gprinter,"1X");
          --centinel;
        }
        break;
      case INS: // SAM Deletion
        GT_OUTPUT_SAM_CIGAR_REVERSE_MATCH();
        gt_gprintf(gprinter,"%"PRIu64"D",gt_misms_get_size(misms));
        break;
      case DEL: // SAM Insertion
        centinel-=gt_misms_get_size(misms);
        GT_OUTPUT_SAM_CIGAR_REVERSE_MATCH();
        gt_gprintf(gprinter,"%"PRIu64"I",gt_misms_get_size(misms));
        break;
      default:
        gt_error(SELECTION_NOT_VALID);
        return GT_SOE_PRINTING_MISM_STRING;
        break;
    }
    --misms_n;
  }
  if (centinel >= 0) gt_gprintf(gprinter,"%"PRIu64"M",centinel);
  return 0;
}
GT_INLINE gt_status gt_output_sam_gprint_map_block_cigar_forward(gt_generic_printer* const gprinter,gt_map* const map,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  const uint64_t map_length = gt_map_get_base_length(map);
  uint64_t centinel = 0;
  GT_MISMS_ITERATE(map,misms) {
    const uint64_t misms_pos = gt_misms_get_position(misms);
    switch (misms->misms_type) {
      case MISMS:
        if (attributes->print_mismatches) {
          GT_OUTPUT_SAM_CIGAR_FORWARD_MATCH();
          gt_gprintf(gprinter,"1X");
          ++centinel;
        }
        break;
      case INS: // SAM Deletion
        GT_OUTPUT_SAM_CIGAR_FORWARD_MATCH();
        gt_gprintf(gprinter,"%"PRIu64"D",gt_misms_get_size(misms));
        break;
      case DEL: // SAM Insertion
        GT_OUTPUT_SAM_CIGAR_FORWARD_MATCH();
        gt_gprintf(gprinter,"%"PRIu64"I",gt_misms_get_size(misms));
        centinel+=gt_misms_get_size(misms);
        break;
      default:
        gt_error(SELECTION_NOT_VALID);
        return GT_SOE_PRINTING_MISM_STRING;
        break;
    }
  }
  if (centinel < map_length) gt_gprintf(gprinter,"%"PRIu64"M",map_length-centinel);
  return 0;
}
GT_INLINE gt_status gt_output_sam_gprint_map_block_cigar(
    gt_generic_printer* const gprinter,gt_map* const map_block,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map_block);
  gt_status error_code = 0;
  if (gt_map_get_strand(map_block)==REVERSE) {
    // Check following map blocks
    gt_map* const next_map_block = gt_map_get_next_block(map_block);
    if (next_map_block!=NULL && GT_MAP_IS_SAME_SEGMENT(map_block,next_map_block)) { // SplitMap (Otherwise is a quimera)
      error_code = gt_output_sam_gprint_map_block_cigar(gprinter,next_map_block,attributes);
      gt_gprintf(gprinter,"%"PRIu64"N",gt_map_get_junction_size(map_block));
    }
    // Print CIGAR for current map block
    gt_output_sam_gprint_map_block_cigar_reverse(gprinter,map_block,attributes);
  } else {
    // Print CIGAR for current map block
    gt_output_sam_gprint_map_block_cigar_forward(gprinter,map_block,attributes);
    // Check following map blocks
    gt_map* const next_map_block = gt_map_get_next_block(map_block);
    if (next_map_block!=NULL && GT_MAP_IS_SAME_SEGMENT(map_block,next_map_block)) { // SplitMap (Otherwise is a quimera)
      gt_gprintf(gprinter,"%"PRIu64"N",gt_map_get_junction_size(map_block));
      error_code = gt_output_sam_gprint_map_block_cigar(gprinter,next_map_block,attributes);
    }
  }
  return error_code;
}
GT_INLINE gt_status gt_output_sam_gprint_map_cigar(
    gt_generic_printer* const gprinter,gt_map* const map_segment,gt_output_sam_attributes* const attributes,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map_segment);
  gt_status error_code = 0;
  // Check strandness
  if (gt_map_get_strand(map_segment)==FORWARD) {
    if (hard_left_trim_read>0) gt_gprintf(gprinter,"%"PRIu64"H",hard_left_trim_read);
    error_code=gt_output_sam_gprint_map_block_cigar(gprinter,map_segment,attributes);
    if (hard_right_trim_read>0) gt_gprintf(gprinter,"%"PRIu64"H",hard_right_trim_read);
  } else {
    if (hard_right_trim_read>0) gt_gprintf(gprinter,"%"PRIu64"H",hard_right_trim_read);
    error_code=gt_output_sam_gprint_map_block_cigar(gprinter,map_segment,attributes);
    if (hard_left_trim_read>0) gt_gprintf(gprinter,"%"PRIu64"H",hard_left_trim_read);
  }
  return error_code;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map_segment,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_cigar,gt_map* const map_segment,gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_cigar(gt_generic_printer* const gprinter,gt_map* const map_segment,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map_segment);
  return gt_output_sam_gprint_map_cigar(gprinter,map_segment,attributes,0,0);
}
/*
 * SAM CORE fields
 *   (QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL). No EOL is printed
 *   Don't handle quimeras (just print one record out of the first map segment)
 */
GT_INLINE void gt_output_sam_gprint_map_placeholder_xa(gt_generic_printer* const gprinter,gt_map_placeholder* const map_ph,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_NULL_CHECK(map_ph);
  if (map_ph->map==NULL) {
    gt_gprintf(gprinter,";"); return;
  }
  /*
   * XA maps (chr12,+91022,101M,0)
   */
  gt_gprintf(gprinter,PRIgts",%c%lu,",
      PRIgts_content(map_ph->map->seq_name),
      (map_ph->map->strand==FORWARD)?'+':'-',
      gt_map_get_global_coordinate(map_ph->map)); // Print the map
  gt_output_sam_gprint_map_cigar(gprinter,map_ph->map,attributes,map_ph->hard_trim_left,map_ph->hard_trim_right);
  gt_gprintf(gprinter,",%"PRIu64";",gt_map_get_levenshtein_distance(map_ph->map));
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities,map,position,phred_score, \
    hard_left_trim_read,hard_right_trim_read,secondary_alignment,not_passing_QC,PCR_duplicate,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_core_fields_se,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map* const map,const uint64_t position,const uint8_t phred_score,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_core_fields_se(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map* const map,const uint64_t position,const uint8_t phred_score,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  // (1) Print QNAME
  gt_output_sam_gprint_qname(gprinter,tag);
  // (2) Print FLAG
  gt_gprintf(gprinter,"\t%"PRId16,gt_output_sam_calculate_flag_se_map(map,secondary_alignment,not_passing_QC,PCR_duplicate));
  // Is mapped?
  if (gt_expect_true(map!=NULL)) {
    // (3) Print RNAME
    // (4) Print POS
    // (5) Print MAPQ
    gt_gprintf(gprinter,"\t"PRIgts"\t%"PRIu64"\t%"PRIu8"\t",PRIgts_content(map->seq_name),position,phred_score);
    // (6) Print CIGAR
    gt_output_sam_gprint_map_cigar(gprinter,map,attributes,hard_left_trim_read,hard_right_trim_read);
  } else {
    // (3) Print RNAME
    // (4) Print POS
    // (5) Print MAPQ
    // (6) Print CIGAR
    gt_gprintf(gprinter,"\t*\t0\t255\t*");
  }
  //  (7) Print RNEXT
  //  (8) Print PNEXT
  //  (9) Print TLEN
  // (10) Print SEQ
  // (11) Print QUAL
  if (!gt_string_is_null(read) && !gt_string_is_null(qualities)) {
    gt_gprintf(gprinter,"\t*\t0\t0\t"PRIgts"\t"PRIgts,
        PRIgts_trimmed_content(read,hard_left_trim_read,hard_right_trim_read),
        PRIgts_trimmed_content(qualities,hard_left_trim_read,hard_right_trim_read));
  } else if (!gt_string_is_null(read)) {
    gt_gprintf(gprinter,"\t*\t0\t0\t"PRIgts"\t*",PRIgts_trimmed_content(read,hard_left_trim_read,hard_right_trim_read));
  } else if (!gt_string_is_null(qualities)) {
    gt_gprintf(gprinter,"\t*\t0\t0\t*\t"PRIgts,PRIgts_trimmed_content(qualities,hard_left_trim_read,hard_right_trim_read));
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities, \
    map,position,phred_score,mate,mate_position,template_length, \
    hard_left_trim_read,hard_right_trim_read,is_map_first_in_pair,secondary_alignment,not_passing_QC,PCR_duplicate,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_core_fields_pe,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map* const map,const uint64_t position,const uint8_t phred_score,
    gt_map* const mate,const uint64_t mate_position,const int64_t template_length,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool is_map_first_in_pair,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_core_fields_pe(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map* const map,const uint64_t position,const uint8_t phred_score,
    gt_map* const mate,const uint64_t mate_position,const int64_t template_length,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool is_map_first_in_pair,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  // (1) Print QNAME
  gt_output_sam_gprint_qname(gprinter,tag);
  // (2) Print FLAG
  gt_gprintf(gprinter,"\t%"PRId16,gt_output_sam_calculate_flag_pe_map(
      map,mate,is_map_first_in_pair,secondary_alignment,not_passing_QC,PCR_duplicate));
  // (3) Print RNAME
  // (4) Print POS
  // (5) Print MAPQ
  // (6) Print CIGAR
  if (map!=NULL) {
    gt_gprintf(gprinter,"\t"PRIgts"\t%"PRIu64"\t%"PRIu8"\t",PRIgts_content(map->seq_name),position,phred_score);
    gt_output_sam_gprint_map_cigar(gprinter,map,attributes,hard_left_trim_read,hard_right_trim_read); // CIGAR
  } else {
    gt_gprintf(gprinter,"\t*\t0\t255\t*");
  }
  // (7) Print RNEXT
  // (8) Print PNEXT
  // (9) Print TLEN
  if (mate!=NULL) {
    if (map!=NULL && !gt_string_equals(map->seq_name,mate->seq_name)) {
      gt_gprintf(gprinter,"\t"PRIgts"\t%"PRIu64"\t%"PRId64,PRIgts_content(mate->seq_name),mate_position,template_length);
    } else {
      gt_gprintf(gprinter,"\t=\t%"PRIu64"\t%"PRId64,mate_position,template_length);
    }
  } else {
    gt_gprintf(gprinter,"\t*\t0\t0");
  }
  // (10) Print SEQ
  // (11) Print QUAL
  if (!gt_string_is_null(read) && !gt_string_is_null(qualities)) {
    gt_gprintf(gprinter,"\t"PRIgts"\t"PRIgts,
        PRIgts_trimmed_content(read,hard_left_trim_read,hard_right_trim_read),
        PRIgts_trimmed_content(qualities,hard_left_trim_read,hard_right_trim_read));
  } else if (!gt_string_is_null(read)) {
    gt_gprintf(gprinter,"\t"PRIgts"\t*",PRIgts_trimmed_content(read,hard_left_trim_read,hard_right_trim_read));
  } else if (!gt_string_is_null(qualities)) {
    gt_gprintf(gprinter,"\t*\t"PRIgts,PRIgts_trimmed_content(qualities,hard_left_trim_read,hard_right_trim_read));
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities,map_segment,hard_left_trim_read,hard_right_trim_read,secondary_alignment,not_passing_QC,PCR_duplicate,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_map_core_fields_se,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_map* const map_segment,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_map_core_fields_se(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_map* const map_segment,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  return gt_output_sam_gprint_core_fields_se(gprinter,tag,read,qualities,
      map_segment,
      (map_segment!=NULL) ? gt_map_get_global_coordinate(map_segment) : 0,
      (map_segment!=NULL) ? gt_map_get_phred_score(map_segment) : GT_MAP_NO_PHRED_SCORE,
      hard_left_trim_read,hard_right_trim_read,secondary_alignment,not_passing_QC,PCR_duplicate,attributes);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities,map_segment,mate_segment,mmap_attributes, \
    hard_left_trim_read,hard_right_trim_read,is_map_first_in_pair,secondary_alignment,not_passing_QC,PCR_duplicate,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_map_core_fields_pe,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map* const map_segment,gt_map* const mate_segment,gt_mmap_attributes* const mmap_attributes,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool is_map_first_in_pair,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_map_core_fields_pe(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map* const map_segment,gt_map* const mate_segment,gt_mmap_attributes* const mmap_attributes,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool is_map_first_in_pair,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  return gt_output_sam_gprint_core_fields_pe(gprinter,tag,read,qualities,
      map_segment,
      (map_segment!=NULL) ? gt_map_get_global_coordinate(map_segment) : 0,
      (map_segment!=NULL) ? gt_map_get_phred_score(map_segment) : GT_MAP_NO_PHRED_SCORE,
      mate_segment,
      (mate_segment!=NULL) ? gt_map_get_global_coordinate(mate_segment) : 0,
      (map_segment!=NULL && mate_segment!=NULL) ? gt_map_get_observed_template_size(map_segment,mate_segment) : 0,
      hard_left_trim_read,hard_right_trim_read,is_map_first_in_pair,secondary_alignment,not_passing_QC,PCR_duplicate,attributes);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,read,qualities,map_placeholder,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_map_placeholder,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map_placeholder* const map_placeholder,gt_output_sam_attributes* const output_attributes);
GT_INLINE gt_status gt_output_sam_gprint_map_placeholder(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map_placeholder* const map_placeholder,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  if (map_placeholder->type==GT_MAP_PLACEHOLDER) {
    return gt_output_sam_gprint_map_core_fields_se(gprinter,tag,read,qualities,map_placeholder->map,
        map_placeholder->hard_trim_left,map_placeholder->hard_trim_right,
        map_placeholder->secondary_alignment,map_placeholder->not_passing_QC,map_placeholder->PCR_duplicate,
        output_attributes);
  } else {
    return gt_output_sam_gprint_map_core_fields_pe(gprinter,tag,read,qualities,
        map_placeholder->map,map_placeholder->paired_end.mate,map_placeholder->paired_end.mmap_attributes,
        map_placeholder->hard_trim_left,map_placeholder->hard_trim_right,
        map_placeholder->paired_end.paired_end_position==0,
        map_placeholder->secondary_alignment,map_placeholder->not_passing_QC,map_placeholder->PCR_duplicate,output_attributes);
  }
}
/*
 * SAM Optional fields
 *   - SAM Attributes is a shash of gt_sam_attribute (gt_sam_data_attributes.h)
 *   - SAM Attributes (gt_sam_data_attributes.h) can be either a value(int,double,string)
 *       or a function -> f(gt_sam_attribute_func_params* params) returning a value(int,double,string)
 *   - gt_output_sam_print_optional_fields_values() prints all the values contained in @sam_attributes
 *     gt_output_sam_print_optional_fields() prints all attributes.
 *       Those relying on a function, are generating calling that function with @gt_sam_attribute_func_params
 *       as argument (some fields can be NULL, so the attribute function must be ready to deal with that)
 */
GT_INLINE gt_sam_attributes* gt_output_sam_select_sam_attributes(gt_map_placeholder* const map_ph) {
  GT_NULL_CHECK(map_ph);
  // Look into the map
  gt_sam_attributes* attributes = NULL;
  if (map_ph->map!=NULL) {
    attributes = gt_attributes_get_sam_attributes(map_ph->map->attributes);
    if (attributes!=NULL) return attributes;
  }
  // Look into the alignment
  const gt_alignment* alignment = NULL;
  if (map_ph->type==GT_MAP_PLACEHOLDER) {
    alignment = map_ph->single_end.alignment;
  } else if (map_ph->paired_end.template != NULL) {
    alignment = gt_template_get_block(map_ph->paired_end.template,map_ph->paired_end.paired_end_position);
  }
  if (alignment!=NULL) {
    attributes = gt_attributes_get_sam_attributes(alignment->attributes);
    if (attributes!=NULL) return attributes;
  }
  return NULL;
}

#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS sam_attribute,attribute_func_params
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_sam_attribute,
    gt_sam_attribute* const sam_attribute,gt_sam_attribute_func_params* const attribute_func_params);
GT_INLINE gt_status gt_output_sam_gprint_sam_attribute(gt_generic_printer* const gprinter,
    gt_sam_attribute* const sam_attribute,gt_sam_attribute_func_params* const attribute_func_params) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_NULL_CHECK(sam_attribute);
  switch (sam_attribute->attribute_type) {
    // Values
    case SAM_ATTR_INT_VALUE:
      gt_gprintf(gprinter,"%c%c:%c:%ld",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,sam_attribute->i_value);
      break;
    case SAM_ATTR_FLOAT_VALUE:
      gt_gprintf(gprinter,"%c%c:%c:%3.2E",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,sam_attribute->f_value);
      break;
    case SAM_ATTR_STRING_VALUE:
      gt_gprintf(gprinter,"%c%c:%c:"PRIgts,sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,PRIgts_content(sam_attribute->s_value));
      break;
    // Functions
    case SAM_ATTR_INT_FUNC:
      if (sam_attribute->i_func(attribute_func_params)==0) { // Generate i-value
        gt_gprintf(gprinter,"%c%c:%c:%ld",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,attribute_func_params->return_i);
      }
      break;
    case SAM_ATTR_FLOAT_FUNC:
      if (sam_attribute->f_func(attribute_func_params)==0) { // Generate f-value
        gt_gprintf(gprinter,"%c%c:%c:%3.2E",sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,attribute_func_params->return_f);
      }
      break;
    case SAM_ATTR_STRING_FUNC:
      if (sam_attribute->s_func(attribute_func_params)==0) { // Generate s-value
        gt_gprintf(gprinter,"%c%c:%c:"PRIgts,sam_attribute->tag[0],sam_attribute->tag[1],sam_attribute->type_id,PRIgts_content(attribute_func_params->return_s));
      }
      break;
    default:
      GT_INVALID_CASE();
      break;
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS sam_attributes,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_optional_fields_values,
    gt_sam_attributes* sam_attributes,gt_output_sam_attributes* const output_attributes);
GT_INLINE gt_status gt_output_sam_gprint_optional_fields_values(gt_generic_printer* const gprinter,
    gt_sam_attributes* sam_attributes,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  if (!output_attributes->print_optional_fields) return 0;
  if (sam_attributes==NULL) sam_attributes = output_attributes->sam_attributes;
  if (sam_attributes!=NULL) {
    GT_SAM_ATTRIBUTES_CHECK(sam_attributes);
    GT_SAM_ATTRIBUTES_BEGIN_ITERATE(sam_attributes,sam_attribute) {
      switch (sam_attribute->attribute_type) {
        case SAM_ATTR_INT_VALUE:
        case SAM_ATTR_FLOAT_VALUE:
        case SAM_ATTR_STRING_VALUE:
          gt_gprintf(gprinter,"\t");
          gt_output_sam_gprint_sam_attribute(gprinter,sam_attribute,NULL);
          break;
        default:
          break;
      }
    } GT_SAM_ATTRIBUTES_END_ITERATE;
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS sam_attributes,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_optional_fields,
    gt_sam_attributes* sam_attributes,gt_output_sam_attributes* const output_attributes);
GT_INLINE gt_status gt_output_sam_gprint_optional_fields(gt_generic_printer* const gprinter,
    gt_sam_attributes* sam_attributes,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  if (!output_attributes->print_optional_fields) return 0;
  if (sam_attributes!=NULL) {
    GT_SAM_ATTRIBUTES_CHECK(sam_attributes);
    GT_SAM_ATTRIBUTES_BEGIN_ITERATE(sam_attributes,sam_attribute) {
      gt_gprintf(gprinter,"\t");
      gt_output_sam_gprint_sam_attribute(gprinter,sam_attribute,output_attributes->attribute_func_params);
    } GT_SAM_ATTRIBUTES_END_ITERATE;
  }
  if (output_attributes->sam_attributes!=NULL) { // TODO: some short of OPTFIELD preference
    GT_SAM_ATTRIBUTES_CHECK(output_attributes->sam_attributes);
    GT_SAM_ATTRIBUTES_BEGIN_ITERATE(output_attributes->sam_attributes,sam_attribute) {
      gt_gprintf(gprinter,"\t");
      gt_output_sam_gprint_sam_attribute(gprinter,sam_attribute,output_attributes->attribute_func_params);
    } GT_SAM_ATTRIBUTES_END_ITERATE;
  }
  return 0;
}
/*
 * SAM Placeholders Printers
 */
GT_INLINE gt_status gt_output_sam_gprint_map_placeholder_vector_se_compact_xa_list(gt_generic_printer* const gprinter,
    gt_vector* const map_placeholder,const uint64_t primary_position,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_VECTOR_CHECK(map_placeholder);
  GT_NULL_CHECK(attributes);
  if (attributes->max_printable_maps == 0) return 0;
  if (gt_vector_get_used(map_placeholder) > 1) {
    gt_gprintf(gprinter,"\tXA:Z:");
    GT_VECTOR_ITERATE(map_placeholder,map_ph,map_placeholder_position,gt_map_placeholder) {
      // Filter PH
      if (map_ph->type!=GT_MAP_PLACEHOLDER ||
          map_placeholder_position==primary_position) continue;
      // Print the map (XA)
      gt_output_sam_gprint_map_placeholder_xa(gprinter,map_ph,attributes);
    }
  }
  return 0;
}
GT_INLINE gt_status gt_output_sam_gprint_map_placeholder_vector_pe_compact_xa_list(gt_generic_printer* const gprinter,
    gt_vector* const map_placeholder,const uint64_t primary_position,const uint64_t end_position,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_VECTOR_CHECK(map_placeholder);
  GT_NULL_CHECK(attributes);
  if (attributes->max_printable_maps == 0) return 0;
  if (gt_vector_get_used(map_placeholder) > 2) {
    gt_gprintf(gprinter,"\tXA:Z:");
    GT_VECTOR_ITERATE(map_placeholder,map_ph,map_placeholder_position,gt_map_placeholder) {
      // Filter PH
      if (map_ph->type==GT_MAP_PLACEHOLDER ||
          map_ph->paired_end.paired_end_position!=end_position ||
          map_placeholder_position==primary_position) continue;
      // Print the map (XA)
      gt_output_sam_gprint_map_placeholder_xa(gprinter,map_ph,attributes);
    }
  }
  return 0;
}
GT_INLINE gt_status gt_output_sam_gprint_map_placeholder_se_compact(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_vector* const map_placeholder_vector,const uint64_t primary_position,
    gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_NULL_CHECK(tag);
  GT_VECTOR_CHECK(map_placeholder_vector);
  GT_NULL_CHECK(attributes);
  gt_status error_code = 0;
  // Produce RC of the mapping's read/qualities
  gt_string* const read_f = read;
  gt_string* const qualities_f = qualities;
  gt_string *read_rc = NULL, *qualities_r = NULL;
  // Get primary map
  gt_cond_error(primary_position>=gt_vector_get_used(map_placeholder_vector),OUTPUT_SAM_NO_PRIMARY_ALG);
  gt_map_placeholder* const primary_map_ph = gt_vector_get_elm(map_placeholder_vector,primary_position,gt_map_placeholder);
  gt_map* const primary_map = primary_map_ph->map;
  // Print primary MAP
  if (primary_map==NULL || gt_map_get_strand(primary_map)==FORWARD) {
    error_code |= gt_output_sam_gprint_map_placeholder(gprinter,tag,read_f,qualities_f,primary_map_ph,attributes);
  } else {
    if (gt_expect_false(read_rc==NULL)) { // Check RC
      read_rc = gt_dna_string_reverse_complement_dup(read_f);
      qualities_r = gt_string_reverse_dup(qualities_f);
    }
    error_code |= gt_output_sam_gprint_map_placeholder(gprinter,tag,read_rc,qualities_r,primary_map_ph,attributes);
  }
  // Print XA:Z field
  gt_output_sam_gprint_map_placeholder_vector_se_compact_xa_list(gprinter,map_placeholder_vector,primary_position,attributes);
  // Print Optional Fields
  gt_sam_attribute_func_params_set_alignment_info(attributes->attribute_func_params,primary_map_ph); // Set func params for OF
  gt_output_sam_gprint_optional_fields(gprinter,gt_output_sam_select_sam_attributes(primary_map_ph),attributes);
  gt_gprintf(gprinter,"\n");
  // Free
  if (read_rc!=NULL) {
    gt_string_delete(read_rc);
    gt_string_delete(qualities_r);
  }
  return error_code;
}
GT_INLINE gt_status gt_output_sam_gprint_map_placeholder_pe_compact(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_vector* const map_placeholder_vector,const uint64_t primary_position,const uint64_t end_position,
    gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_NULL_CHECK(tag);
  GT_VECTOR_CHECK(map_placeholder_vector);
  GT_NULL_CHECK(attributes);
  gt_status error_code = 0;
  // Produce RC of the mapping's read/qualities
  gt_string* const read_f = read;
  gt_string* const qualities_f = qualities;
  gt_string* read_rc = NULL;
  gt_string* qualities_r = NULL;
  // Get primary map
  gt_cond_error(primary_position>=gt_vector_get_used(map_placeholder_vector),OUTPUT_SAM_NO_PRIMARY_ALG);
  gt_map_placeholder* const primary_map_ph = gt_vector_get_elm(map_placeholder_vector,primary_position,gt_map_placeholder);
  gt_map* const primary_map = primary_map_ph->map;
  // Print primary MAP
  if (primary_map==NULL || gt_map_get_strand(primary_map)==FORWARD) {
    error_code |= gt_output_sam_gprint_map_placeholder(gprinter,tag,read_f,qualities_f,primary_map_ph,attributes);
  } else {
    if (gt_expect_false(read_rc==NULL)) { // Check RC
      read_rc = gt_dna_string_reverse_complement_dup(read_f);
      qualities_r = gt_string_reverse_dup(qualities_f);
    }
    error_code |= gt_output_sam_gprint_map_placeholder(gprinter,tag,read_rc,qualities_r,primary_map_ph,attributes);
  }
  // Print XA:Z field
  gt_output_sam_gprint_map_placeholder_vector_pe_compact_xa_list(gprinter,map_placeholder_vector,primary_position,end_position,attributes);
  // Print Optional Fields
  gt_sam_attribute_func_params_set_alignment_info(attributes->attribute_func_params,primary_map_ph); // Set func params for OF
  gt_output_sam_gprint_optional_fields(gprinter,gt_output_sam_select_sam_attributes(primary_map_ph),attributes);
  gt_gprintf(gprinter,"\n");
  // Free
  if (read_rc!=NULL) {
    gt_string_delete(read_rc);
    gt_string_delete(qualities_r);
  }
  return error_code;
}
GT_INLINE gt_status gt_output_sam_gprint_map_placeholder_vector_se(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_vector* const map_placeholder_vector,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_VECTOR_CHECK(map_placeholder_vector);
  GT_NULL_CHECK(attributes);
  gt_status error_code = 0;
  // Produce RC of the mapping's read/qualities
  gt_string *read_f = read, *qualities_f = qualities;
  gt_string *read_rc=NULL, *qualities_r=NULL;
  // Iterate over all placeholders
  GT_VECTOR_ITERATE(map_placeholder_vector,map_ph,map_ph_it,gt_map_placeholder) {
    if (map_ph->type!=GT_MAP_PLACEHOLDER) continue;
    // Print MAP
    if (map_ph->map==NULL || gt_map_get_strand(map_ph->map)==FORWARD) {
      error_code |= gt_output_sam_gprint_map_placeholder(gprinter,tag,read_f,qualities_f,map_ph,attributes);
    } else {
      if (gt_expect_false(read_rc==NULL && read_f!=NULL)) { // Check RC
        read_rc = gt_dna_string_reverse_complement_dup(read_f);
        qualities_r = gt_string_reverse_dup(qualities_f);
      }
      error_code |= gt_output_sam_gprint_map_placeholder(gprinter,tag,read_rc,qualities_r,map_ph,attributes);
    }
    // Print Optional Fields
    gt_sam_attribute_func_params_set_alignment_info(attributes->attribute_func_params,map_ph); // Set func params for OF
    gt_output_sam_gprint_optional_fields(gprinter,gt_output_sam_select_sam_attributes(map_ph),attributes);
    gt_gprintf(gprinter,"\n");
    // Nullify read & qualities
    if (gt_expect_false(!attributes->always_output_read__qualities && map_ph_it>0)) {
      read_f = NULL; qualities_f = NULL;
      if (read_rc!=NULL) {
        gt_string_delete(read_rc); read_rc = NULL;
        gt_string_delete(qualities_r); qualities_r = NULL;
      }
    }
  }
  // Free
  if (read_rc!=NULL) {
    gt_string_delete(read_rc);
    gt_string_delete(qualities_r);
  }
  return error_code;
}
GT_INLINE gt_status gt_output_sam_gprint_map_placeholder_vector_pe(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_string* const read_end1,gt_string* const read_end2,
    gt_string* const qualities_end1,gt_string* const qualities_end2,
    gt_vector* const map_placeholder_vector,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_VECTOR_CHECK(map_placeholder_vector);
  GT_NULL_CHECK(attributes);
  gt_status error_code = 0;
  // Produce RC of the mapping's read/qualities
  gt_string *read_f_end1 = read_end1, *read_f_end2 = read_end2;
  gt_string *qualities_f_end1 = qualities_end1, *qualities_f_end2 = qualities_end2;
  gt_string *read_rc_end1=NULL, *qualities_r_end1=NULL, *read_rc_end2=NULL, *qualities_r_end2=NULL;
  // Iterate over all placeholders
  GT_VECTOR_ITERATE(map_placeholder_vector,map_ph,map_ph_it,gt_map_placeholder) {
    if (map_ph->type==GT_MAP_PLACEHOLDER) continue;
    // Print MAP
    if (map_ph->paired_end.paired_end_position==0) {
      if (map_ph->map==NULL || gt_map_get_strand(map_ph->map)==FORWARD) {
        error_code |= gt_output_sam_gprint_map_placeholder(gprinter,tag,read_f_end1,qualities_f_end1,map_ph,attributes);
      } else {
        if (gt_expect_false(read_rc_end1==NULL)) { // Check RC
          read_rc_end1 = gt_dna_string_reverse_complement_dup(read_f_end1);
          qualities_r_end1 = gt_string_reverse_dup(qualities_f_end1);
        }
        error_code |= gt_output_sam_gprint_map_placeholder(gprinter,tag,read_rc_end1,qualities_r_end1,map_ph,attributes);
      }
    } else {
      if (map_ph->map==NULL || gt_map_get_strand(map_ph->map)==FORWARD) {
        error_code |= gt_output_sam_gprint_map_placeholder(gprinter,tag,read_f_end2,qualities_f_end2,map_ph,attributes);
      } else {
        if (gt_expect_false(read_rc_end2==NULL)) { // Check RC
          read_rc_end2 = gt_dna_string_reverse_complement_dup(read_f_end2);
          qualities_r_end2 = gt_string_reverse_dup(qualities_f_end2);
        }
        error_code |= gt_output_sam_gprint_map_placeholder(gprinter,tag,read_rc_end2,qualities_r_end2,map_ph,attributes);
      }
    }
    // Print Optional Fields
    gt_sam_attribute_func_params_set_alignment_info(attributes->attribute_func_params,map_ph); // Set func params for OF
    gt_output_sam_gprint_optional_fields(gprinter,gt_output_sam_select_sam_attributes(map_ph),attributes);
    gt_gprintf(gprinter,"\n");
    // Nullify read & qualities
    if (gt_expect_false(!attributes->always_output_read__qualities && map_ph_it>0)) {
      read_f_end1 = NULL; qualities_f_end1 = NULL;
      read_f_end2 = NULL; qualities_f_end2 = NULL;
      if (read_rc_end1!=NULL) {
        gt_string_delete(read_rc_end1); gt_string_delete(qualities_r_end1);
        read_rc_end1 = NULL; qualities_r_end1 = NULL;
      }
      if (read_rc_end2!=NULL) {
        gt_string_delete(read_rc_end2); gt_string_delete(qualities_r_end2);
        read_rc_end2 = NULL; qualities_r_end2 = NULL;
      }
    }
  }
  // Free
  if (read_rc_end1!=NULL) {
    gt_string_delete(read_rc_end1); gt_string_delete(qualities_r_end1);
  }
  if (read_rc_end2!=NULL) {
    gt_string_delete(read_rc_end2); gt_string_delete(qualities_r_end2);
  }
  return error_code;
}
GT_INLINE gt_status gt_output_sam_gprint_alignment_map_placeholder_vector(gt_generic_printer* const gprinter,
    gt_alignment* const alignment,gt_vector* const map_placeholder_vector,const uint64_t primary_position,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  GT_VECTOR_CHECK(map_placeholder_vector);
  GT_NULL_CHECK(attributes);
  // Check qualities
  gt_status error_code = 0;
  gt_string* qualities = alignment->qualities;
  if (attributes->qualities_offset == GT_QUALS_OFFSET_64) {
    qualities = gt_qualities_dup__adapt_offset64_to_offset33(alignment->qualities);
  }
  if (attributes->compact_format) {
    // Print all placeholders (compact)
    error_code = gt_output_sam_gprint_map_placeholder_se_compact(gprinter,
        alignment->tag,alignment->read,qualities,map_placeholder_vector,primary_position,attributes);
  } else {
    // Print all placeholders (explicit)
    error_code = gt_output_sam_gprint_map_placeholder_vector_se(gprinter,
        alignment->tag,alignment->read,qualities,map_placeholder_vector,attributes);
  }
  if (attributes->qualities_offset == GT_QUALS_OFFSET_64) gt_string_delete(qualities);
  return error_code;
}
GT_INLINE gt_status gt_output_sam_gprint_template_map_placeholder_vector(gt_generic_printer* const gprinter,
    gt_template* const template,gt_vector* const map_placeholder_vector,
    const uint64_t primary_position_end1,const uint64_t primary_position_end2,gt_output_sam_attributes* const attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  GT_VECTOR_CHECK(map_placeholder_vector);
  GT_NULL_CHECK(attributes);
  gt_status error_code = 0;
  gt_alignment* const alignment_end1 = gt_template_get_end1(template);
  gt_alignment* const alignment_end2 = gt_template_get_end2(template);
  gt_string* qualities_end1 = alignment_end1->qualities;
  gt_string* qualities_end2 = alignment_end2->qualities;
  if (attributes->qualities_offset == GT_QUALS_OFFSET_64) {
    qualities_end1 = gt_qualities_dup__adapt_offset64_to_offset33(alignment_end1->qualities);
    qualities_end2 = gt_qualities_dup__adapt_offset64_to_offset33(alignment_end2->qualities);
  }
  if (attributes->compact_format) {
    // Print End/1
    error_code|=gt_output_sam_gprint_map_placeholder_pe_compact(gprinter,template->tag,alignment_end1->read,alignment_end1->qualities,
        map_placeholder_vector,primary_position_end1,0,attributes);
    // Print End/2
    error_code|=gt_output_sam_gprint_map_placeholder_pe_compact(gprinter,template->tag,alignment_end2->read,alignment_end2->qualities,
        map_placeholder_vector,primary_position_end2,1,attributes);
  } else {
    // Print all placeholders
    error_code|=gt_output_sam_gprint_map_placeholder_vector_pe(gprinter,
        template->tag,alignment_end1->read,alignment_end2->read,
        alignment_end1->qualities,alignment_end2->qualities,map_placeholder_vector,attributes);
  }
  if (attributes->qualities_offset == GT_QUALS_OFFSET_64) {
    gt_string_delete(qualities_end1); gt_string_delete(qualities_end2);
  }
  return error_code;
}
/*
 * SAM High-level MMap/Map Printers
 */
#define GT_OUTPUT_SAM_SETUP_READ__QUALITIES(alignment,map,read_f,qualities_f,read_rc,qualities_r,rc_complement) \
  /* Produce RC of the mapping's read/qualities */ \
  rc_complement = map!=NULL && gt_map_get_strand(map)==REVERSE; \
  if (alignment!=NULL) { \
    if (alignment->read!=NULL) { \
      read_f = alignment->read; \
      if (rc_complement) read_rc = gt_dna_string_reverse_complement_dup(read_f); \
    } \
    if (alignment->qualities!=NULL) { \
      qualities_f = alignment->qualities; \
      if (rc_complement) qualities_r = gt_string_reverse_dup(qualities_f); \
    } \
  }
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,map_end1,map_end2,mmap_attributes,secondary_alignment,not_passing_QC,PCR_duplicate,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_mmap,
    gt_template* const template,gt_map* const map_end1,gt_map* const map_end2,gt_mmap_attributes* const mmap_attributes,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,gt_output_sam_attributes* const output_attributes);
GT_INLINE gt_status gt_output_sam_gprint_mmap(gt_generic_printer* const gprinter,
    gt_template* const template,gt_map* const map_end1,gt_map* const map_end2,gt_mmap_attributes* const mmap_attributes,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(output_attributes);
  gt_status error_code;
  // Create placeholder template
  gt_map_placeholder ph;
  gt_map_placeholder_set_sam_fields(&ph,not_passing_QC,PCR_duplicate,0,0);
  ph.paired_end.template = template;
  // Allocate ph-vector
  gt_vector* const map_placeholder = gt_vector_new(10,sizeof(gt_map_placeholder));
  // Add mmap
  uint64_t primary_position_end1 = UINT64_MAX, primary_position_end2 = UINT64_MAX;
  gt_map_placeholder *best_ph_end1 = NULL, *best_ph_end2 = NULL;
  gt_map_placeholder_add_mmap(map_end1,map_end2,gt_template_get_end1(template)->read,0,
      map_placeholder,true,gt_map_placeholder_cmp_map_phred_scores,true,best_ph_end1,&primary_position_end1,&ph);
  gt_map_placeholder_add_mmap(map_end2,map_end1,gt_template_get_end2(template)->read,1,
      map_placeholder,true,gt_map_placeholder_cmp_map_phred_scores,true,best_ph_end2,&primary_position_end2,&ph);
  // Print maps !!
  error_code = gt_output_sam_gprint_template_map_placeholder_vector(gprinter,template,
      map_placeholder,primary_position_end1,primary_position_end2,output_attributes);
  // Free
  gt_vector_delete(map_placeholder);
  return error_code;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,alignment,map,secondary_alignment,not_passing_QC,PCR_duplicate,output_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_map,
    gt_template* const template,gt_alignment* const alignment,gt_map* const map,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,gt_output_sam_attributes* const output_attributes);
GT_INLINE gt_status gt_output_sam_gprint_map(gt_generic_printer* const gprinter,
    gt_template* const template,gt_alignment* const alignment,gt_map* const map,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(output_attributes);
  gt_status error_code;
  // Create placeholder template
  gt_map_placeholder ph;
  gt_map_placeholder_set_sam_fields(&ph,not_passing_QC,PCR_duplicate,0,0);
  ph.single_end.template = template;
  ph.single_end.alignment = alignment;
  // Allocate ph-vector
  gt_vector* const map_placeholder = gt_vector_new(4,sizeof(gt_map_placeholder));
  // Add map
  uint64_t primary_position = UINT64_MAX;
  gt_map_placeholder *best_ph = NULL;
  gt_map_placeholder_add_map(map,alignment->read,map_placeholder,true,gt_map_placeholder_cmp_map_phred_scores,true,best_ph,&primary_position,&ph);
  // Print maps !!
  error_code = gt_output_sam_gprint_alignment_map_placeholder_vector(gprinter,
      alignment,map_placeholder,primary_position,output_attributes);
  // Free
  gt_vector_delete(map_placeholder);
  return error_code;
}
/*
 * SAM High-level Template/Alignment Printers
 *   - Optional fields are generated from the first SAM-Attributes object found in the following order:
 *       1.- map->attributes{GT_ATTR_ID_SAM_FLAGS} / mmap_attributes->attributes{GT_ATTR_ID_SAM_FLAGS}
 *       2.- @output_attributes->sam_attributes
 */
GT_INLINE gt_status gt_output_sam_gprint_alignment_(gt_generic_printer* const gprinter,
    gt_alignment* const alignment,gt_map_placeholder* const ph,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(output_attributes);
  gt_status error_code;
  // Create ph-vector with all maps
  const uint64_t num_maps = gt_alignment_get_num_maps(alignment);
  uint64_t primary_position;
  gt_vector* const map_placeholder = gt_vector_new(num_maps,sizeof(gt_map_placeholder));
  gt_map_placeholder_build_from_alignment(alignment,map_placeholder,true,output_attributes->max_printable_maps,
      gt_map_placeholder_cmp_map_phred_scores,&primary_position,ph);
  gt_attributes_add(output_attributes->attribute_func_params->attributes,GT_ATTR_ID_SAM_TAG_NH,&gt_vector_get_used(map_placeholder),uint64_t); // TAG_NH
  // Print maps !!
  error_code = gt_output_sam_gprint_alignment_map_placeholder_vector(gprinter,
      alignment,map_placeholder,primary_position,output_attributes);
  // Free
  gt_vector_delete(map_placeholder);
  return error_code;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_alignment,gt_alignment* const alignment,gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_alignment(gt_generic_printer* const gprinter,gt_alignment* const alignment,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(output_attributes);
  // Get passingQC and PCRDuplicate flags
  bool passing_QC = true, PCR_duplicate = false, *aux;
  aux = gt_attributes_get(alignment->attributes,GT_ATTR_ID_SAM_PASSING_QC);
  if (aux!=NULL) passing_QC = *aux;
  aux = gt_attributes_get(alignment->attributes,GT_ATTR_ID_SAM_PCR_DUPLICATE);
  if (aux!=NULL) PCR_duplicate = *aux;
  // Fill Ph template
  gt_map_placeholder ph;
  gt_map_placeholder_set_sam_fields(&ph,!passing_QC,PCR_duplicate,0,0);
  return gt_output_sam_gprint_alignment_(gprinter,alignment,&ph,output_attributes);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_sam,print_template,gt_template* const template,gt_output_sam_attributes* const attributes);
GT_INLINE gt_status gt_output_sam_gprint_template(gt_generic_printer* const gprinter,gt_template* const template,gt_output_sam_attributes* const output_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(output_attributes);
  // Get passingQC and PCRDuplicate flags
  bool passing_QC = true, PCR_duplicate = false, *aux;
  aux = gt_attributes_get(template->attributes,GT_ATTR_ID_SAM_PASSING_QC);
  if (aux!=NULL) passing_QC = *aux;
  aux = gt_attributes_get(template->attributes,GT_ATTR_ID_SAM_PCR_DUPLICATE);
  if (aux!=NULL) PCR_duplicate = *aux;
  gt_map_placeholder ph;
  gt_map_placeholder_set_sam_fields(&ph,!passing_QC,PCR_duplicate,0,0);
  // Handle reduction to alignment
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    ph.single_end.template = template;
    return gt_output_sam_gprint_alignment_(gprinter,alignment,&ph,output_attributes);
  } GT_TEMPLATE_END_REDUCTION;
  gt_status error_code;
  // Create ph-vector with all mmaps
  const uint64_t num_maps = gt_template_get_num_mmaps(template);
  uint64_t primary_position_end1, primary_position_end2;
  gt_vector* const map_placeholder = gt_vector_new(num_maps,sizeof(gt_map_placeholder));
  gt_map_placeholder_build_from_template(template,map_placeholder,true,true,output_attributes->max_printable_maps,
      gt_map_placeholder_cmp_map_phred_scores,&primary_position_end1,&primary_position_end2,&ph);
  gt_attributes_add(output_attributes->attribute_func_params->attributes,GT_ATTR_ID_SAM_TAG_NH,&gt_vector_get_used(map_placeholder),uint64_t); // TAG_NH
  // Print maps !!
  error_code = gt_output_sam_gprint_template_map_placeholder_vector(gprinter,template,
      map_placeholder,primary_position_end1,primary_position_end2,output_attributes);
  // Free
  gt_vector_delete(map_placeholder);
  return error_code;
}


