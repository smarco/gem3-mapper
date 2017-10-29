/*
 * PROJECT: GEM-Tools library
 * FILE: gt_sam_data_attributes.h
 * DATE: 15/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides support to SAM format data structures
 */

#ifndef GT_SAM_DATA_ATTRIBUTES_H_
#define GT_SAM_DATA_ATTRIBUTES_H_

#include "gt_essentials.h"
#include "gt_generic_printer.h"

#include "gt_map.h"
#include "gt_alignment.h"
#include "gt_template_utils.h"
#include "gt_sequence_archive.h"

/*
 * Checkers
 */
#define GT_SAM_HEADERS_CHECK(sam_headers) \
  GT_STRING_CHECK(sam_headers->header); \
  GT_VECTOR_CHECK(sam_headers->read_group); \
  GT_VECTOR_CHECK(sam_headers->program); \
  GT_VECTOR_CHECK(sam_headers->comments)

/*
 * SAM FLAGS
 */
#define GT_SAM_FLAG_MULTIPLE_SEGMENTS 0x1
#define GT_SAM_FLAG_PROPERLY_ALIGNED 0x2
#define GT_SAM_FLAG_UNMAPPED 0x4
#define GT_SAM_FLAG_NEXT_UNMAPPED 0x8
#define GT_SAM_FLAG_REVERSE_COMPLEMENT 0x10
#define GT_SAM_FLAG_NEXT_REVERSE_COMPLEMENT 0x20
#define GT_SAM_FLAG_FIRST_SEGMENT 0x40
#define GT_SAM_FLAG_LAST_SEGMENT 0x80
#define GT_SAM_FLAG_SECONDARY_ALIGNMENT 0x100
#define GT_SAM_FLAG_NOT_PASSING_QC 0x200
#define GT_SAM_FLAG_PCR_OR_OPTICAL_DUPLICATE 0x400
/*
 * SAM File specifics Attribute (SAM Headers)
 */
#define GT_SAM_HEADER_OK   0
#define GT_SAM_HEADER_TAG_INVALID   (-1)
#define GT_SAM_HEADER_VALUE_INVALID (-2)
typedef struct {
  gt_string* header; // @HD
  gt_sequence_archive* sequence_archive; // @SQ
  gt_vector* read_group; /* // @RG (gt_string*) */
  gt_vector* program; // @PG /* (gt_string*) */
  gt_vector* comments; // @ CO /* (gt_string*) */
} gt_sam_headers; // SAM Headers

GT_INLINE gt_sam_headers* gt_sam_header_new(void);
GT_INLINE void gt_sam_header_clear(gt_sam_headers* const sam_headers);
GT_INLINE void gt_sam_header_delete(gt_sam_headers* const sam_headers);

GT_INLINE void gt_sam_header_set_sequence_archive(gt_sam_headers* const sam_headers,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_sam_header_set_header_record(gt_sam_headers* const sam_headers,gt_string* const header_line);
GT_INLINE void gt_sam_header_add_read_group_record(gt_sam_headers* const sam_headers,gt_string* const read_group_record);
GT_INLINE void gt_sam_header_add_program_record(gt_sam_headers* const sam_headers,gt_string* const program_record);
GT_INLINE void gt_sam_header_add_comment(gt_sam_headers* const sam_headers,gt_string* const comment);
/*
 * SAM Optional Fields
 *   - SAM Attributes(optional fields) are just a hash of @gt_sam_attribute
 *     embedded into the general attributes(@gt_shash) of any object(@template,@alignment,@map,...)
 */
typedef enum { SAM_ATTR_INT_VALUE, SAM_ATTR_FLOAT_VALUE, SAM_ATTR_STRING_VALUE,
               SAM_ATTR_INT_FUNC,  SAM_ATTR_FLOAT_FUNC,  SAM_ATTR_STRING_FUNC } gt_sam_attribute_t;
typedef gt_shash gt_sam_attributes;
typedef struct {
  /* Return Values
   *   Depending on the function type, the proper field will be returned/output
   *   as the value for the TAG (function)
   */
  gt_string* return_s; // (In case of returning a gt_string, use this buffer)
  float return_f;
  int32_t return_i;
  /*
   * Alignment record info
   *   Contains all info related to the current SAM record
   *   as to fetch the template/alignment/map info needed
   */
  gt_sequence_archive *sequence_archive; /* Sequence Archive */
  gt_map_placeholder *alignment_info; /* Template/Alignment/Map/etc */
  /*
   * Placeholder Vector
   *   Contains all map/mmap that are going to be printed; pre-processed).
   *   When printing a all maps from a template/alignment a placeholder-vector
   *   is built so to simplify the process
   */
  gt_vector* mmap_placeholder; // gt_vector<gt_map_placeholder>
  /* Attributes */
  gt_attributes* attributes;
} gt_sam_attribute_func_params;
typedef struct {
  char tag[2];
  gt_sam_attribute_t attribute_type;
  char type_id;
  union {
    /* Values */
    int32_t i_value;
    float f_value;
    gt_string* s_value;
    /* Functions */
    gt_status (*i_func)(gt_sam_attribute_func_params*);
    gt_status (*f_func)(gt_sam_attribute_func_params*);
    gt_status (*s_func)(gt_sam_attribute_func_params*);
  };
} gt_sam_attribute;
/*
 * SAM Optional Fields
 *   - SAM Attributes(optional fields) are just a hash of @gt_sam_attribute
 *     embedded into the general attributes(@gt_shash) of any object(@template,@alignment,@map,...)
 */
#define GT_SAM_ATTRIBUTES_CHECK(sam_attributes) GT_HASH_CHECK((gt_shash*)(sam_attributes))
GT_INLINE gt_sam_attributes* gt_sam_attributes_new();
GT_INLINE void gt_sam_attributes_clear(gt_sam_attributes* const sam_attributes);
GT_INLINE void gt_sam_attributes_delete(gt_sam_attributes* const sam_attributes);
GT_INLINE gt_sam_attribute* gt_sam_attributes_get_attribute(gt_sam_attributes* const sam_attributes,char* const tag);
GT_INLINE void gt_sam_attributes_add_attribute(gt_sam_attributes* const sam_attributes,gt_sam_attribute* const sam_attribute);
// General Attributes API
GT_INLINE gt_sam_attributes* gt_attributes_get_sam_attributes(gt_attributes* const general_attributes);
GT_INLINE gt_sam_attributes* gt_attributes_get_sam_attributes_dyn(gt_attributes* const general_attributes);
GT_INLINE void gt_attributes_delete_sam_attributes(gt_attributes* const general_attributes);
GT_INLINE void gt_attributes_clear_sam_attributes(gt_attributes* const general_attributes);
GT_INLINE bool gt_attributes_has_sam_attributes(gt_attributes* const general_attributes);
GT_INLINE gt_sam_attribute* gt_attributes_get_sam_attribute(gt_attributes* const general_attributes,char* const tag);
// Iterator over SAM attributes
#define GT_SAM_ATTRIBUTES_BEGIN_ITERATE(sam_attributes,attribute) \
  GT_SHASH_BEGIN_ELEMENT_ITERATE(sam_attributes,attribute,gt_sam_attribute)
#define GT_SAM_ATTRIBUTES_END_ITERATE GT_SHASH_END_ITERATE
#define GT_ATTRIBUTES_SAM_ATTRIBUTES_BEGIN_ITERATE(general_attributes,attribute) \
  gt_sam_attributes* const sam_attributes_from_##general_attributes = gt_attribute_get_sam_attributes_dyn(general_attributes,GT_ATTR_ID_SAM_ATTRIBUTES);  \
  GT_SAM_ATTRIBUTES_BEGIN_ITERATE(sam_attributes,attribute)
#define GT_ATTRIBUTES_SAM_ATTRIBUTES_END_ITERATE GT_SHASH_END_ITERATE
/*
 * SAM Attribute Setup
 */
GT_INLINE void gt_sam_attribute_set_ivalue(gt_sam_attribute* const sam_attribute,const char* const tag,const char type_id,const int32_t value);
GT_INLINE void gt_sam_attribute_set_fvalue(gt_sam_attribute* const sam_attribute,const char* const tag,const char type_id,const float value);
GT_INLINE void gt_sam_attribute_set_svalue(gt_sam_attribute* const sam_attribute,const char* const tag,const char type_id,gt_string* const string);
GT_INLINE void gt_sam_attribute_set_ifunc(gt_sam_attribute* const sam_attribute,const char* const tag,const char type_id,gt_status (*i_func)(gt_sam_attribute_func_params*));
GT_INLINE void gt_sam_attribute_set_ffunc(gt_sam_attribute* const sam_attribute,const char* const tag,const char type_id,gt_status (*f_func)(gt_sam_attribute_func_params*));
GT_INLINE void gt_sam_attribute_set_sfunc(gt_sam_attribute* const sam_attribute,const char* const tag,const char type_id,gt_status (*s_func)(gt_sam_attribute_func_params*));
/*
 * SAM Attributes Add
 */
GT_INLINE void gt_sam_attributes_add_ivalue(gt_sam_attributes* const sam_attributes,const char* const tag,const char type_id,const int32_t value);
GT_INLINE void gt_sam_attributes_add_fvalue(gt_sam_attributes* const sam_attributes,const char* const tag,const char type_id,const float value);
GT_INLINE void gt_sam_attributes_add_svalue(gt_sam_attributes* const sam_attributes,const char* const tag,const char type_id,gt_string* const string);
GT_INLINE void gt_sam_attributes_add_ifunc(gt_sam_attributes* const sam_attributes,const char* const tag,const char type_id,gt_status (*i_func)(gt_sam_attribute_func_params*));
GT_INLINE void gt_sam_attributes_add_ffunc(gt_sam_attributes* const sam_attributes,const char* const tag,const char type_id,gt_status (*f_func)(gt_sam_attribute_func_params*));
GT_INLINE void gt_sam_attributes_add_sfunc(gt_sam_attributes* const sam_attributes,const char* const tag,const char type_id,gt_status (*s_func)(gt_sam_attribute_func_params*));
/*
 * Functional Attribute (ifunc,ffunc,sfunc,pfunc) parameters
 */
GT_INLINE gt_sam_attribute_func_params* gt_sam_attribute_func_params_new();
GT_INLINE void gt_sam_attribute_func_params_delete(gt_sam_attribute_func_params* const func_params);
GT_INLINE void gt_sam_attribute_func_params_clear(gt_sam_attribute_func_params* const func_params);
GT_INLINE void gt_sam_attribute_func_params_set_sequence_archive(
    gt_sam_attribute_func_params* const func_params,gt_sequence_archive* const sequence_archive);
GT_INLINE void gt_sam_attribute_func_params_set_sam_flags(
    gt_sam_attribute_func_params* const func_params,const bool not_passing_QC,const bool PCR_duplicate);

GT_INLINE gt_map_placeholder* gt_sam_attribute_func_params_get_alignment_info(gt_sam_attribute_func_params* const func_params);
GT_INLINE void gt_sam_attribute_func_params_set_alignment_info(gt_sam_attribute_func_params* const func_params,gt_map_placeholder* const map_placeholder);

GT_INLINE void gt_sam_attribute_func_params_set_se(
    gt_sam_attribute_func_params* const func_params,gt_template* const template,gt_alignment* const alignment,gt_map* const map);
GT_INLINE void gt_sam_attribute_func_params_set_pe(
    gt_sam_attribute_func_params* const func_params,gt_template* const template,
    uint64_t paired_end_position,gt_map* const map,gt_map* const mate,gt_mmap_attributes* const mmap_attributes);

/*
 * Functional Internal Data
 *   Storage for calculation to do once for all template/alignmet/map SAM output
 *   So that fields like NH are calculated just once, printed multiple times and erased at the end
 */

/*
 * GT-library PRE-Implemented Functional Attributes
 *   SAM specification predefined
 */

//  AM  i  The smallest template-independent mapping quality of segments in the rest // FIXME

//  AS  i  Alignment score generated by aligner
// TODO

//  BC  Z  Barcode sequence
//  BQ  Z  Offset to base alignment quality (BAQ), of the same length as the read sequence.
//         At the i-th read base, BAQi = Qi - (BQi - 64) where Qi is the i-th base quality. // FIXME
//  CC  Z  Reference name of the next hit; "=" for the same chromosome
//  CM  i  Edit distance between the color sequence and the color reference (see also NM)
//  CP  i  Leftmost coordinate of the next hit
//  CQ  Z  Color read quality on the original strand of the read. Same encoding as QUAL; same length as CS.
//  CS  Z  Color read sequence on the original strand of the read. The primer base must be included.
//  E2  Z  The 2nd most likely base calls. Same encoding and same length as QUAL.
//  FI  i  The index of segment in the template.
//  FS  Z  Segment suffix.
//  FZ  B,S Flow signal intensities on the original strand of the read, stored as (uint16 t) round(value * 100.0).
//  LB  Z  Library. Value to be consistent with the header RG-LB tag if @RG is present.
//  H0  i  Number of perfect hits
//  H1  i  Number of 1-difference hits (see also NM)
//  H2  i  Number of 2-difference hits
//  HI  i  Query hit index, indicating the alignment record is the i-th one stored in SAM
//  IH  i  Number of stored alignments in SAM that contains the query in the current record
//  MD  Z  String for mismatching positions. Regex : [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
//  MQ  i  Mapping quality of the mate/next segment
GT_INLINE void gt_sam_attributes_add_tag_MQ(gt_sam_attributes* const sam_attributes);

//  NH  i  Number of reported alignments that contains the query in the current record
GT_INLINE void gt_sam_attributes_add_tag_NH(gt_sam_attributes* const sam_attributes);

//  NM  i  Edit distance to the reference, including ambiguous bases but excluding clipping
GT_INLINE void gt_sam_attributes_add_tag_NM(gt_sam_attributes* const sam_attributes);

//  OQ  Z  Original base quality (usually before recalibration). Same encoding as QUAL.
//  OP  i  Original mapping position (usually before realignment)
//  OC  Z  Original CIGAR (usually before realignment)
//  PG  Z  Program. Value matches the header PG-ID tag if @PG is present.
//  PQ  i  Phred likelihood of the template, conditional on both the mapping being correct
GT_INLINE void gt_sam_attributes_add_tag_PQ(gt_sam_attributes* const sam_attributes);
//  PU  Z  Platform unit. Value to be consistent with the header RG-PU tag if @RG is present.
//  Q2  Z  Phred quality of the mate/next segment. Same encoding as QUAL.
//  R2  Z  Sequence of the mate/next segment in the template.

//  RG  Z  Read group. Value matches the header RG-ID tag if @RG is present in the header.
GT_INLINE void gt_sam_attributes_add_tag_RG(gt_sam_attributes* const sam_attributes,gt_string* const read_group);

//  SM  i  Template-independent mapping quality
//  TC  i  The number of segments in the template.
//  U2  Z  Phred probability of the 2nd call being wrong conditional on the best being wrong. The same encoding as QUAL.
//  UQ  i  Phred likelihood of the segment, conditional on the mapping being correct
GT_INLINE void gt_sam_attributes_add_tag_UQ(gt_sam_attributes* const sam_attributes);

/*
 * GT-library PRE-Implemented Functional Attributes
 *   X?  ?  Reserved fields for end users (together with Y? and Z?)
 */

//  X0  i  Number of best hits
//  X1  i  Number of suboptimal hits found by BWA
//  XN  i  Number of ambiguous bases in the reference
//  XM  i  Number of mismatches in the alignment
//  XO  i  Number of gap opens
//  XG  i  Number of gap extentions
/*
 *  XT  A  Type: Unique/Repeat/N/Mate-sw
 * i.e.
 *   XT:A:U => Unique alignment
 *   XT:A:R => Repeat
 *   XT:A:N => Not mapped
 *   XT:A:M => Mate-sw (Read is fixed due to paired end rescue)
 */
GT_INLINE void gt_sam_attributes_add_tag_XT(gt_sam_attributes* const sam_attributes);

/*
 *  XA  Z  Alternative hits; format: (chr,pos,CIGAR,NM;)*
 *    => Implemented by means of setting ...
 *    =>   gt_output_sam_attributes_compact_format_on(@gt_output_sam_attributes)
 *    => at gt_output_sam.h
 */

//  XS  ?  Suboptimal alignment score
//  XF  ?  Support from forward/reverse alignment
//  XE  i  Number of supporting seeds

//  cs  Z  Casava TAG (if any)
GT_INLINE void gt_sam_attributes_add_tag_cs(gt_sam_attributes* const sam_attributes);
//  md  Z  GEM CIGAR String
GT_INLINE void gt_sam_attributes_add_tag_md(gt_sam_attributes* const sam_attributes);
//  XS  A  XS directionality information
GT_INLINE void gt_sam_attributes_add_tag_XS(gt_sam_attributes* const sam_attributes);

#endif /* GT_SAM_DATA_ATTRIBUTES_H_ */
