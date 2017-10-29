/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_sam.h
 * DATE: 01/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_OUTPUT_SAM_H_
#define GT_OUTPUT_SAM_H_

#include "gt_essentials.h"
#include "gt_dna_string.h"
#include "gt_dna_read.h"
#include "gt_alignment_utils.h"
#include "gt_template_utils.h"
#include "gt_sam_attributes.h"
#include "gt_map_score.h"
#include "gt_output_buffer.h"
#include "gt_buffered_output_file.h"
#include "gt_generic_printer.h"

/*
 * Error Codes
 */
#define GT_SOE_PRINTING_MISM_STRING 10

/*
 * SAM Output attributes
 */
typedef enum { GT_SAM, GT_BAM } gt_output_sam_format_t;
typedef struct {
  /* Format */
  gt_output_sam_format_t format; // TODO
  /* Read/Qualities */
  bool always_output_read__qualities;
  gt_qualities_offset_t qualities_offset;
  /* Maps */
  uint64_t max_printable_maps;
  bool compact_format; // Compact map representation in SAM via XA field
  /* CIGAR/Mismatch string */
  bool print_mismatches;
  /* SAM Optional Fields */
  bool print_optional_fields;
  gt_sam_attributes* sam_attributes; // Optional fields stored as sam_attributes
  gt_sam_attribute_func_params* attribute_func_params; // Parameters provided to generate functional attributes
} gt_output_sam_attributes;
/*
 * BAM record
 */
typedef struct {
  int32_t block_size;  // Length of the remainder of the alignment record
  int32_t refID; // Reference sequence ID, 1 <= refID < n ref; -1 for a read without a mapping position.
  int32_t pos; // pos 0-based leftmost coordinate (= POS  1)
  uint32_t bin_mq_nl; // bin<<16|MAPQ<<8|l read name;
                      //   bin is computed by the reg2bin() function in Section 4.3;
                      //   l read name is the length of read name below (= length(QNAME) + 1).
  uint32_t flag_mq_nl; // FLAG<<16|n cigar op; n cigar op is the number of operations in CIGAR.
  int32_t l_seq; // Length of SEQ
  int32_t next_refID; // Ref-ID of the next segment (1 <= mate refID < n_ref)
  int32_t next_pos; // 0-based leftmost pos of the next segment (= PNEXT - 1)
  int32_t tlen; // Template length (= TLEN)
  gt_string* read_name; // Read name, NULL terminated (QNAME plus a tailing `\0')
  uint32_t* cigar; // op_len<<4|op. `MIDNSHP=X'!`012345678'
  uint8_t* seq; // 4-bit encoded read: `=ACMGRSVTWYHKDBN'! [0; 15];
                //   Other characters mapped to `N';
                //   High nybble first (1st base in the highest 4-bit of the 1st byte)
  char* qual; // Phred base quality (a sequence of 0xFF if absent)
  gt_string* optional_fields; // Optional Fields String. List of auxiliary data (until the end of the alignment block)
} gt_bam_record;

/* Setup */
GT_INLINE gt_output_sam_attributes* gt_output_sam_attributes_new();
GT_INLINE void gt_output_sam_attributes_delete(gt_output_sam_attributes* const attributes);
GT_INLINE void gt_output_sam_attributes_clear(gt_output_sam_attributes* const attributes);
/* Format */
GT_INLINE void gt_output_sam_attributes_set_format(gt_output_sam_attributes* const attributes,gt_output_sam_format_t const format);
/* Read/Qualities */
GT_INLINE void gt_output_sam_attributes_dump_read__qualities_once(gt_output_sam_attributes* const attributes);
GT_INLINE void gt_output_sam_attributes_always_dump_read__qualities(gt_output_sam_attributes* const attributes);
GT_INLINE void gt_output_sam_attributes_set_qualities_offset(gt_output_sam_attributes* const attributes,gt_qualities_offset_t const qualities_offset);
/* Maps */
GT_INLINE void gt_output_sam_attributes_set_max_printable_maps(gt_output_sam_attributes* const attributes,const uint64_t max_printable_maps);
GT_INLINE void gt_output_sam_attributes_set_compact_format(gt_output_sam_attributes* const attributes,const bool compact_format);
/* CIGAR/Mismatch string */
GT_INLINE void gt_output_sam_attributes_set_print_mismatches(gt_output_sam_attributes* const attributes,const bool print_mismatches);
/* SAM Optional Fields */
GT_INLINE void gt_output_sam_attributes_set_print_optional_fields(gt_output_sam_attributes* const attributes,const bool print_optional_fields);
GT_INLINE void gt_output_sam_attributes_set_reference_sequence_archive(gt_output_sam_attributes* const attributes,gt_sequence_archive* const reference_sequence_archive);
GT_INLINE gt_sam_attributes* gt_output_sam_attributes_get_sam_attributes(gt_output_sam_attributes* const attributes);

/*
 * SAM Headers
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_headers_sh,gt_sam_headers* const sam_headers);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_headers_sa,gt_sequence_archive* const sequence_archive);
/*
 * SAM QNAME (Tag)
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_qname,gt_string* const tag);
/*
 * SAM Flag
 *   Beware of the SAM flags, they might cause severe mind injuries...
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
 * - Bit 0x4 is the only reliable place to tell whether the segment is unmapped. If 0x4 is set,
 *     no assumptions can be made about RNAME, POS, CIGAR, MAPQ, bits 0x2, 0x10 and 0x100
 *     and the bit 0x20 of the next segment in the template.
 * - If 0x40 and 0x80 are both set, the segment is part of a linear template, but it is neither
 *     the first nor the last segment. If both 0x40 and 0x80 are unset, the index of the segment
 *     in the template is unknown. This may happen for a non-linear template or the index is
 *     lost in data processing.
 * - Bit 0x100 marks the alignment not to be used in certain analysis when the tools in use
 *     are aware of this bit.
 * - If 0x1 is unset, no assumptions can be made about 0x2, 0x8, 0x20, 0x40 and 0x80.
*/
GT_INLINE uint16_t gt_output_sam_calculate_flag(
    const bool paired_end,const bool read_paired,
    const bool read_mapped,const bool mate_mapped,
    const gt_strand read_strand,const gt_strand mate_strand,
    const bool first_in_pair,const bool last_in_pair,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate);
GT_INLINE uint16_t gt_output_sam_calculate_flag_pe(
    const bool read_paired,const bool read_mapped,const bool mate_mapped,
    const gt_strand read_strand,const gt_strand mate_strand,
    const bool first_in_pair,const bool last_in_pair,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate);
GT_INLINE uint16_t gt_output_sam_calculate_flag_se(
    const bool read_mapped,const gt_strand read_strand,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate);
GT_INLINE uint16_t gt_output_sam_calculate_flag_se_map(
    gt_map* const map,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate);
GT_INLINE uint16_t gt_output_sam_calculate_flag_pe_map(
    gt_map* const map,gt_map* const mate,const bool is_map_first_in_pair,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate);
/*
 * SAM CIGAR
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_cigar,gt_map* const map_segment,gt_output_sam_attributes* const attributes);
/*
 * SAM CORE fields
 *   (QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL). No EOL is printed
 *   Don't handle quimeras (just print one record out of the first map segment)
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_core_fields_se,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map* const map,const uint64_t position,const uint8_t phred_score,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_core_fields_pe,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map* const map,const uint64_t position,const uint8_t phred_score,
    gt_map* const mate,const uint64_t mate_position,const int64_t template_length,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool is_map_first_in_pair,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_map_core_fields_se,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,gt_map* const map_segment,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_map_core_fields_pe,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map* const map_segment,gt_map* const mate_segment,gt_mmap_attributes* const mmap_attributes,
    const uint64_t hard_left_trim_read,const uint64_t hard_right_trim_read,
    const bool is_map_first_in_pair,const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const output_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_map_placeholder,
    gt_string* const tag,gt_string* const read,gt_string* const qualities,
    gt_map_placeholder* const map_placeholder,gt_output_sam_attributes* const output_attributes);
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
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_sam_attribute,gt_sam_attribute* const sam_attribute,gt_sam_attribute_func_params* const attribute_func_params);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_optional_fields_values,gt_sam_attributes* const sam_attributes,gt_output_sam_attributes* const output_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_optional_fields,gt_sam_attributes* const sam_attributes,gt_output_sam_attributes* const output_attributes);
/*
 * SAM High-level MMap/Map Printers
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_mmap,
    gt_template* const template,gt_map* const map_end1,gt_map* const map_end2,gt_mmap_attributes* const mmap_attributes,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const output_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_map,
    gt_template* const template,gt_alignment* const alignment,gt_map* const map,
    const bool secondary_alignment,const bool not_passing_QC,const bool PCR_duplicate,
    gt_output_sam_attributes* const output_attributes);
/*
 * SAM High-level Template/Alignment Printers
 *   - Optional fields are generated from the first SAM-Attributes object found in the following order:
 *       1.- map->attributes{GT_ATTR_ID_SAM_FLAGS} / mmap_attributes->attributes{GT_ATTR_ID_SAM_FLAGS}
 *       2.- @output_attributes->sam_attributes
 */
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_alignment,gt_alignment* const alignment,gt_output_sam_attributes* const output_attributes);
GT_GENERIC_PRINTER_PROTOTYPE(gt_output_sam,print_template,gt_template* const template,gt_output_sam_attributes* const output_attributes);

/*
 * Error Messages
 */
#define GT_ERROR_OUTPUT_SAM_NO_PRIMARY_ALG "Output SAM. No primary alignment specified"

#endif /* GT_OUTPUT_SAM_H_ */
