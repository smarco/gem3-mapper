/*
 * PROJECT: GEM-Tools library
 * FILE: gt_sam_data_attributes.c
 * DATE: 15/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides support to SAM format data structures
 */


#include "gt_sam_attributes.h"

/*
 * SAM File specifics Attribute (SAM Headers)
 */
#define GT_ATTR_SAM_INIT_ELEMENTS 2
GT_INLINE gt_sam_headers* gt_sam_header_new(void) {
  gt_sam_headers* sam_headers = gt_alloc(gt_sam_headers);
  sam_headers->header = gt_string_new(50); // @HD
  sam_headers->read_group = gt_vector_new(GT_ATTR_SAM_INIT_ELEMENTS,sizeof(gt_string*)); // @RG
  sam_headers->program = gt_vector_new(GT_ATTR_SAM_INIT_ELEMENTS,sizeof(gt_string*)); // @PG
  sam_headers->comments = gt_vector_new(GT_ATTR_SAM_INIT_ELEMENTS,sizeof(gt_string*)); // @ CO
  sam_headers->sequence_archive = NULL; // @SQ
  return sam_headers;
}
GT_INLINE void gt_sam_header_clear(gt_sam_headers* const sam_headers) {
  GT_SAM_HEADERS_CHECK(sam_headers);
  // Header
  gt_string_clear(sam_headers->header);
  // Read group
  GT_VECTOR_ITERATE(sam_headers->read_group,record_h,nh,gt_string*) { gt_string_delete(*record_h); }
  gt_vector_clear(sam_headers->read_group);
  // Program
  GT_VECTOR_ITERATE(sam_headers->program,record_p,np,gt_string*) { gt_string_delete(*record_p); }
  gt_vector_clear(sam_headers->program);
  // Comments
  GT_VECTOR_ITERATE(sam_headers->comments,comment,nc,gt_string*) { gt_string_delete(*comment); }
  gt_vector_clear(sam_headers->comments);
  // Seq Archive
  sam_headers->sequence_archive=NULL;
}
GT_INLINE void gt_sam_header_delete(gt_sam_headers* const sam_headers) {
  GT_SAM_HEADERS_CHECK(sam_headers);
  // Clear
  gt_sam_header_clear(sam_headers);
  // Delete
  gt_string_delete(sam_headers->header);
  gt_vector_delete(sam_headers->read_group);
  gt_vector_delete(sam_headers->program);
  gt_vector_delete(sam_headers->comments);
}

GT_INLINE void gt_sam_header_set_sequence_archive(gt_sam_headers* const sam_headers,gt_sequence_archive* const sequence_archive) {
  GT_SAM_HEADERS_CHECK(sam_headers);
  sam_headers->sequence_archive = sequence_archive;
}
GT_INLINE void gt_sam_header_set_header_record(gt_sam_headers* const sam_headers,gt_string* const header_line) {
  GT_SAM_HEADERS_CHECK(sam_headers);
  GT_STRING_CHECK(header_line);
  gt_string_copy(sam_headers->header,header_line);
}
GT_INLINE void gt_sam_header_add_read_group_record(gt_sam_headers* const sam_headers,gt_string* const read_group_record) {
  GT_SAM_HEADERS_CHECK(sam_headers);
  GT_STRING_CHECK(read_group_record);
  gt_vector_insert(sam_headers->read_group,read_group_record,gt_string*);
}
GT_INLINE void gt_sam_header_add_program_record(gt_sam_headers* const sam_headers,gt_string* const program_record) {
  GT_SAM_HEADERS_CHECK(sam_headers);
  GT_STRING_CHECK(program_record);
  gt_vector_insert(sam_headers->program,program_record,gt_string*);
}
GT_INLINE void gt_sam_header_add_comment(gt_sam_headers* const sam_headers,gt_string* const comment) {
  GT_SAM_HEADERS_CHECK(sam_headers);
  GT_STRING_CHECK(comment);
  gt_vector_insert(sam_headers->comments,comment,gt_string*);
}
/*
 * SAM Optional Fields
 *   - SAM Attributes(optional fields) are just a hash of @gt_sam_attribute
 *     embedded into the general attributes(@gt_shash) of any object(@template,@alignment,@map,...)
 */
GT_INLINE gt_sam_attributes* gt_sam_attributes_new() {
  return gt_shash_new();
}
GT_INLINE void gt_sam_attributes_clear(gt_sam_attributes* const sam_attributes) {
  GT_SAM_ATTRIBUTES_CHECK(sam_attributes);
  gt_shash_clear(sam_attributes,true);
}
GT_INLINE void gt_sam_attributes_delete(gt_sam_attributes* const sam_attributes) {
  GT_SAM_ATTRIBUTES_CHECK(sam_attributes);
  gt_shash_delete(sam_attributes,true);
}
GT_INLINE gt_sam_attribute* gt_sam_attributes_get_attribute(gt_sam_attributes* const sam_attributes,char* const tag) {
  GT_SAM_ATTRIBUTES_CHECK(sam_attributes);
  return gt_shash_get_element(sam_attributes,tag);
}
GT_INLINE void gt_sam_attributes_add_attribute(gt_sam_attributes* const sam_attributes,gt_sam_attribute* const sam_attribute) {
  GT_SAM_ATTRIBUTES_CHECK(sam_attributes);
  GT_NULL_CHECK(sam_attribute);
  GT_NULL_CHECK(sam_attribute->tag);
  gt_shash_insert(sam_attributes,sam_attribute->tag,sam_attribute,gt_sam_attribute);
}
/*
 * General Attributes API
 */
GT_INLINE gt_sam_attributes* gt_attributes_get_sam_attributes(gt_attributes* const general_attributes) {
  if (general_attributes==NULL) return NULL;
  return gt_attributes_get(general_attributes,GT_ATTR_ID_SAM_ATTRIBUTES);
}
GT_INLINE gt_sam_attributes* gt_attributes_get_sam_attributes_dyn(gt_attributes* const general_attributes) {
  GT_ATTRIBUTES_CHECK(general_attributes);
  gt_sam_attributes* sam_attributes = gt_attributes_get_sam_attributes(general_attributes);
  if (sam_attributes==NULL) {
    sam_attributes = gt_sam_attributes_new();
    gt_attributes_add_object(general_attributes,GT_ATTR_ID_SAM_ATTRIBUTES,
        sam_attributes,(void*(*)())gt_shash_dup,(void(*)())gt_shash_destroy);
  }
  return sam_attributes;
}
GT_INLINE void gt_attributes_delete_sam_attributes(gt_attributes* const general_attributes) {
  GT_ATTRIBUTES_CHECK(general_attributes);
  gt_attributes_remove(general_attributes,GT_ATTR_ID_SAM_ATTRIBUTES);
}
GT_INLINE void gt_attributes_clear_sam_attributes(gt_attributes* const general_attributes) {
  GT_ATTRIBUTES_CHECK(general_attributes);
  gt_sam_attributes* sam_attributes = gt_attributes_get_sam_attributes(general_attributes);
  if (sam_attributes!=NULL) {
    gt_sam_attributes_clear(sam_attributes);
  }
}
GT_INLINE bool gt_attributes_has_sam_attributes(gt_attributes* const general_attributes) {
  GT_ATTRIBUTES_CHECK(general_attributes);
  return gt_attributes_get_sam_attributes(general_attributes)!=NULL;
}
GT_INLINE gt_sam_attribute* gt_attributes_get_sam_attribute(gt_attributes* const general_attributes,char* const tag) {
  GT_ATTRIBUTES_CHECK(general_attributes);
  gt_sam_attributes* sam_attributes = gt_attributes_get_sam_attributes(general_attributes);
  return (sam_attributes!=NULL) ? gt_sam_attributes_get_attribute(sam_attributes,tag) : NULL;
}
/*
 * SAM Attribute Setup
 */
#define GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag_src) \
  sam_attribute->tag[0]=tag_src[0]; \
  sam_attribute->tag[1]=tag_src[1]
#define GT_ATTRIBUTE_SAM_CMP_TAG(sam_attribute_ptr,tag_src) \
  (sam_attribute_ptr->tag[0]==tag_src[0] && sam_attribute_ptr->tag[1]==tag_src[1])
GT_INLINE void gt_sam_attribute_set_ivalue(gt_sam_attribute* const sam_attribute,const char* const tag,const char type_id,const int32_t value) {
  GT_NULL_CHECK(sam_attribute); // TODO: type checking of i,Z,etc
  GT_NULL_CHECK(tag);
  GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag);
  sam_attribute->attribute_type = SAM_ATTR_INT_VALUE;
  sam_attribute->type_id = type_id;
  sam_attribute->i_value = value;
}
GT_INLINE void gt_sam_attribute_set_fvalue(gt_sam_attribute* const sam_attribute,const char* const tag,const char type_id,const float value) {
  GT_NULL_CHECK(sam_attribute); // TODO: type checking of i,Z,etc
  GT_NULL_CHECK(tag);
  GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag);
  sam_attribute->attribute_type = SAM_ATTR_FLOAT_VALUE;
  sam_attribute->type_id = type_id;
  sam_attribute->f_value = value;
}
GT_INLINE void gt_sam_attribute_set_svalue(gt_sam_attribute* const sam_attribute,const char* const tag,const char type_id,gt_string* const string) {
  GT_NULL_CHECK(sam_attribute); // TODO: type checking of i,Z,etc
  GT_NULL_CHECK(tag);
  GT_STRING_CHECK(string);
  GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag);
  sam_attribute->attribute_type = SAM_ATTR_STRING_VALUE;
  sam_attribute->type_id = type_id;
  sam_attribute->s_value = string;
}
GT_INLINE void gt_sam_attribute_set_ifunc(gt_sam_attribute* const sam_attribute,const char* const tag,const char type_id,gt_status (*i_func)(gt_sam_attribute_func_params*)) {
  GT_NULL_CHECK(sam_attribute); // TODO: type checking of i,Z,etc
  GT_NULL_CHECK(tag);
  GT_NULL_CHECK(i_func);
  GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag);
  sam_attribute->attribute_type = SAM_ATTR_INT_FUNC;
  sam_attribute->type_id = type_id;
  sam_attribute->i_func = i_func;
}
GT_INLINE void gt_sam_attribute_set_ffunc(gt_sam_attribute* const sam_attribute,const char* const tag,const char type_id,gt_status (*f_func)(gt_sam_attribute_func_params*)) {
  GT_NULL_CHECK(sam_attribute); // TODO: type checking of i,Z,etc
  GT_NULL_CHECK(tag);
  GT_NULL_CHECK(f_func);
  GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag);
  sam_attribute->attribute_type = SAM_ATTR_FLOAT_FUNC;
  sam_attribute->type_id = type_id;
  sam_attribute->f_func = f_func;
}
GT_INLINE void gt_sam_attribute_set_sfunc(gt_sam_attribute* const sam_attribute,const char* const tag,const char type_id,gt_status (*s_func)(gt_sam_attribute_func_params*)) {
  GT_NULL_CHECK(sam_attribute); // TODO: type checking of i,Z,etc
  GT_NULL_CHECK(tag);
  GT_NULL_CHECK(s_func);
  GT_ATTRIBUTE_SAM_COPY_TAG(sam_attribute,tag);
  sam_attribute->attribute_type = SAM_ATTR_STRING_FUNC;
  sam_attribute->type_id = type_id;
  sam_attribute->s_func = s_func;
}
/*
 * SAM Attributes Add
 */
GT_INLINE void gt_sam_attributes_add_ivalue(gt_sam_attributes* const sam_attributes,const char* const tag,const char type_id,const int32_t value) {
  GT_SAM_ATTRIBUTES_CHECK(sam_attributes);
  GT_NULL_CHECK(tag);
  gt_sam_attribute* const sam_attribute = gt_alloc(gt_sam_attribute);
  gt_sam_attribute_set_ivalue(sam_attribute,tag,type_id,value);
  gt_sam_attributes_add_attribute(sam_attributes,sam_attribute);
}
GT_INLINE void gt_sam_attributes_add_fvalue(gt_sam_attributes* const sam_attributes,const char* const tag,const char type_id,const float value){
  GT_SAM_ATTRIBUTES_CHECK(sam_attributes);
  GT_NULL_CHECK(tag);
  gt_sam_attribute* const sam_attribute = gt_alloc(gt_sam_attribute);
  gt_sam_attribute_set_fvalue(sam_attribute,tag,type_id,value);
  gt_sam_attributes_add_attribute(sam_attributes,sam_attribute);
}
GT_INLINE void gt_sam_attributes_add_svalue(gt_sam_attributes* const sam_attributes,const char* const tag,const char type_id,gt_string* const string){
  GT_SAM_ATTRIBUTES_CHECK(sam_attributes);
  GT_NULL_CHECK(tag);
  GT_STRING_CHECK(string);
  gt_sam_attribute* const sam_attribute = gt_alloc(gt_sam_attribute);
  gt_sam_attribute_set_svalue(sam_attribute,tag,type_id,string);
  gt_sam_attributes_add_attribute(sam_attributes,sam_attribute);
}
GT_INLINE void gt_sam_attributes_add_ifunc(gt_sam_attributes* const sam_attributes,const char* const tag,const char type_id,gt_status (*i_func)(gt_sam_attribute_func_params*)){
  GT_SAM_ATTRIBUTES_CHECK(sam_attributes);
  GT_NULL_CHECK(tag);
  GT_NULL_CHECK(i_func);
  gt_sam_attribute* const sam_attribute = gt_alloc(gt_sam_attribute);
  gt_sam_attribute_set_ifunc(sam_attribute,tag,type_id,i_func);
  gt_sam_attributes_add_attribute(sam_attributes,sam_attribute);
}
GT_INLINE void gt_sam_attributes_add_ffunc(gt_sam_attributes* const sam_attributes,const char* const tag,const char type_id,gt_status (*f_func)(gt_sam_attribute_func_params*)){
  GT_SAM_ATTRIBUTES_CHECK(sam_attributes);
  GT_NULL_CHECK(tag);
  GT_NULL_CHECK(f_func);
  gt_sam_attribute* const sam_attribute = gt_alloc(gt_sam_attribute);
  gt_sam_attribute_set_ffunc(sam_attribute,tag,type_id,f_func);
  gt_sam_attributes_add_attribute(sam_attributes,sam_attribute);
}
GT_INLINE void gt_sam_attributes_add_sfunc(gt_sam_attributes* const sam_attributes,const char* const tag,const char type_id,gt_status (*s_func)(gt_sam_attribute_func_params*)){
  GT_SAM_ATTRIBUTES_CHECK(sam_attributes);
  GT_NULL_CHECK(tag);
  GT_NULL_CHECK(s_func);
  gt_sam_attribute* const sam_attribute = gt_alloc(gt_sam_attribute);
  gt_sam_attribute_set_sfunc(sam_attribute,tag,type_id,s_func);
  gt_sam_attributes_add_attribute(sam_attributes,sam_attribute);
}
/*
 * Functional Attribute (ifunc,ffunc,sfunc) parameters
 */
#define GT_SAM_ATTR_FUNC_PARAMS_RETURN_S_INIT_LENGTH 20
GT_INLINE gt_sam_attribute_func_params* gt_sam_attribute_func_params_new() {
  gt_sam_attribute_func_params* const func_params = gt_alloc(gt_sam_attribute_func_params);
  /* String (gt_string) buffer */
  func_params->return_s = gt_string_new(GT_SAM_ATTR_FUNC_PARAMS_RETURN_S_INIT_LENGTH);
  /* Attributes */
  func_params->attributes = gt_attributes_new();
  /* Reset defaults */
  gt_sam_attribute_func_params_clear(func_params);
  return func_params;
}
GT_INLINE void gt_sam_attribute_func_params_delete(gt_sam_attribute_func_params* const func_params) {
  GT_NULL_CHECK(func_params);
  /* String (gt_string) buffer */
  if (func_params->return_s!=NULL) gt_string_delete(func_params->return_s);
  /* Attributes */
  gt_attributes_delete(func_params->attributes);
  gt_free(func_params);
}
GT_INLINE void gt_sam_attribute_func_params_clear(gt_sam_attribute_func_params* const func_params) {
  GT_NULL_CHECK(func_params);
  /* Return Values */
  if (func_params->return_s!=NULL) gt_string_clear(func_params->return_s);
  /* Sequence Archive */
  func_params->sequence_archive = NULL;
  /* Template/Alignment source */
  func_params->alignment_info=NULL;
  /* Placeholder Vector */
  func_params->mmap_placeholder = NULL;
}
GT_INLINE void gt_sam_attribute_func_params_set_sequence_archive(
    gt_sam_attribute_func_params* const func_params,gt_sequence_archive* const sequence_archive) {
  GT_NULL_CHECK(func_params);
  func_params->sequence_archive = sequence_archive;
}
GT_INLINE void gt_sam_attribute_func_params_set_sam_flags(
    gt_sam_attribute_func_params* const func_params,const bool not_passing_QC,const bool PCR_duplicate) {
  GT_NULL_CHECK(func_params);
  GT_NULL_CHECK(func_params->alignment_info);
  func_params->alignment_info->not_passing_QC = not_passing_QC;
  func_params->alignment_info->PCR_duplicate = PCR_duplicate;
}
GT_INLINE gt_map_placeholder* gt_sam_attribute_func_params_get_alignment_info(gt_sam_attribute_func_params* const func_params) {
  GT_NULL_CHECK(func_params);
  return func_params->alignment_info;
}
GT_INLINE void gt_sam_attribute_func_params_set_alignment_info(gt_sam_attribute_func_params* const func_params,gt_map_placeholder* const map_placeholder) {
  GT_NULL_CHECK(func_params);
  func_params->alignment_info = map_placeholder;
}
GT_INLINE void gt_sam_attribute_func_params_set_se(
    gt_sam_attribute_func_params* const func_params,gt_template* const template,gt_alignment* const alignment,gt_map* const map) {
  GT_NULL_CHECK(func_params);
  GT_NULL_CHECK(func_params->alignment_info);
  /* Type */
  func_params->alignment_info->type = GT_MAP_PLACEHOLDER;
  /* Map */
  func_params->alignment_info->map = map;
  /* Template/Alignment source */
  func_params->alignment_info->single_end.template = template;
  func_params->alignment_info->single_end.alignment = alignment;
}
GT_INLINE void gt_sam_attribute_func_params_set_pe(
    gt_sam_attribute_func_params* const func_params,gt_template* const template,
    uint64_t paired_end_position,gt_map* const map,gt_map* const mate,gt_mmap_attributes* const mmap_attributes) {
  GT_NULL_CHECK(func_params);
  GT_NULL_CHECK(func_params->alignment_info);
  /* Type */
  func_params->alignment_info->type = (map!=NULL && mate!=NULL) ? GT_MMAP_PLACEHOLDER_PAIRED : GT_MMAP_PLACEHOLDER_UNPAIRED;
  /* Map */
  func_params->alignment_info->map = map;
  /* Template/Alignment source */
  func_params->alignment_info->paired_end.template = template;
  func_params->alignment_info->paired_end.paired_end_position = paired_end_position;
  func_params->alignment_info->paired_end.mate = mate;
  func_params->alignment_info->paired_end.mmap_attributes = mmap_attributes;
}
/*
 * GT-library PRE-Implemented Functional Attributes
 *   SAM specification predefined
 */

//  NH  i  Number of reported alignments that contains the query in the current record
GT_INLINE gt_status gt_sam_attribute_generate_NH(gt_sam_attribute_func_params* func_params) {
  const int32_t* const nh_value = gt_attributes_get(func_params->attributes,GT_ATTR_ID_SAM_TAG_NH);
  if (nh_value==NULL) return -1;
  func_params->return_i = *nh_value;
  return 0;
}
GT_INLINE void gt_sam_attributes_add_tag_NH(gt_sam_attributes* const sam_attributes) {
  gt_sam_attributes_add_ifunc(sam_attributes,"NH",'i',gt_sam_attribute_generate_NH);
}

//  NM  i  Edit distance to the reference, including ambiguous bases but excluding clipping
GT_INLINE gt_status gt_sam_attribute_generate_NM(gt_sam_attribute_func_params* func_params) {
  gt_map* const map = func_params->alignment_info->map;
  if (map == NULL) return -1; // Don't print NM field
  // Get edit distance ( levenshtein distance )
  const int32_t edit_distance = gt_map_get_segment_levenshtein_distance(map);
  // Excluding clipping
  const int32_t edit_distance_no_clipping = edit_distance - gt_map_get_left_trim_length(map) - gt_map_get_right_trim_length(map);
  // Set ivalue
  func_params->return_i = edit_distance_no_clipping;
  // Return OK
  return 0;
}
GT_INLINE void gt_sam_attributes_add_tag_NM(gt_sam_attributes* const sam_attributes) {
  gt_sam_attributes_add_ifunc(sam_attributes,"NM",'i',gt_sam_attribute_generate_NM);
}

//  RG  Z  Read group. Value matches the header RG-ID tag if @RG is present in the header.
GT_INLINE void gt_sam_attributes_add_tag_RG(gt_sam_attributes* const sam_attributes,gt_string* const read_group) {
  gt_sam_attributes_add_svalue(sam_attributes,"RG",'Z',read_group);
}

/*
 * GT-library PRE-Implemented Functional Attributes
 *   X?  ?  Reserved fields for end users (together with Y? and Z?)
 */

/*
 *  XT  A  Type: Unique/Repeat/N/Mate-sw
 * i.e.
 *   XT:A:U => Unique alignment
 *   XT:A:R => Repeat
 *   XT:A:N => Not mapped
 *   XT:A:M => Mate-sw (Read is fixed due to paired end rescue)
 */
typedef enum { GT_XT_UNIQUE, GT_XT_REPEAT, GT_XT_UNMAPPED, GT_XT_MATE_SW } gt_sam_xt_value;
GT_INLINE gt_status gt_sam_attribute_generate_XT(gt_sam_attribute_func_params* func_params) {
  char* xt_char_value_attr = gt_attributes_get(func_params->attributes,GT_ATTR_ID_SAM_TAG_XT);
  char xt_char_value;
  if (xt_char_value_attr==NULL) {
    gt_sam_xt_value xt_value;
    if (func_params->alignment_info->map==NULL) { // Unmapped
      xt_value = GT_XT_UNMAPPED;
    } else if (func_params->alignment_info->type==GT_MAP_PLACEHOLDER) {
      if (func_params->alignment_info->single_end.template!=NULL) {
        const int64_t uniq_degree = gt_template_get_uniq_degree(func_params->alignment_info->single_end.template);
        xt_value = (uniq_degree!=GT_NO_STRATA) ? GT_XT_UNIQUE : GT_XT_REPEAT;
      } else if (func_params->alignment_info->single_end.alignment!=NULL) {
        const int64_t uniq_degree = gt_alignment_get_uniq_degree(func_params->alignment_info->single_end.alignment);
        xt_value = (uniq_degree!=GT_NO_STRATA) ? GT_XT_UNIQUE : GT_XT_REPEAT;
      } else {
        return -1;
      }
    } else { // GT_MMAP_PLACEHOLDER_PAIRED, GT_MMAP_PLACEHOLDER_UNPAIRED
      if (func_params->alignment_info->paired_end.template!=NULL) {
        const int64_t uniq_degree = gt_template_get_uniq_degree(func_params->alignment_info->paired_end.template);
        xt_value = (uniq_degree!=GT_NO_STRATA) ? GT_XT_UNIQUE : GT_XT_REPEAT;
      } else {
        return -1;
      }
    }
    // Set proper value to return
    switch (xt_value) {
      case GT_XT_UNIQUE:   xt_char_value = 'U'; break;
      case GT_XT_REPEAT:   xt_char_value = 'R'; break;
      case GT_XT_UNMAPPED: xt_char_value = 'N'; break;
      case GT_XT_MATE_SW:  xt_char_value = 'M'; break;
      default: return -1; break;
    }
    // Save as Functional Internal Data (let's save computations)
    xt_char_value_attr = &xt_char_value;
    gt_attributes_add(func_params->attributes,GT_ATTR_ID_SAM_TAG_XT,&xt_char_value,char);
  }
  // Return value
  gt_string_clear(func_params->return_s);
  gt_string_append_char(func_params->return_s,*xt_char_value_attr);
  gt_string_append_eos(func_params->return_s);
  return 0;
}
GT_INLINE void gt_sam_attributes_add_tag_XT(gt_sam_attributes* const sam_attributes) {
  gt_sam_attributes_add_sfunc(sam_attributes,"XT",'A',gt_sam_attribute_generate_XT);
}
//  cs  Z  Casava TAG (if any)
GT_INLINE gt_status gt_sam_attribute_generate_cs(gt_sam_attribute_func_params* func_params) {
  gt_string* casava_string = NULL;
  /*
   * Fetch the CASAVA string from the attributes
   *   {GT_MMAP_PLACEHOLDER_PAIRED=0, GT_MMAP_PLACEHOLDER_UNPAIRED=1, GT_MAP_PLACEHOLDER=2}
   */
  if (func_params->alignment_info->type==GT_MAP_PLACEHOLDER) {
    if (func_params->alignment_info->single_end.alignment!=NULL) {
      casava_string = gt_attributes_get(func_params->alignment_info->single_end.alignment->attributes,GT_ATTR_ID_TAG_CASAVA);
    }
  } else {
    if (func_params->alignment_info->paired_end.template!=NULL) {
      casava_string = gt_attributes_get(gt_template_get_block(
              func_params->alignment_info->paired_end.template,
              func_params->alignment_info->paired_end.paired_end_position)->attributes,GT_ATTR_ID_TAG_CASAVA);
    }
  }
  if (casava_string!=NULL) {
    // Copy the CASAVA string to the return_s buffer (as to pass the return string value)
    gt_string_copy(func_params->return_s,casava_string);
    return 0; // Ok, go ahead and print it
  } else {
    return -1; // I don't have such info (don't print this field)
  }
}
GT_INLINE void gt_sam_attributes_add_tag_cs(gt_sam_attributes* const sam_attributes) {
  gt_sam_attributes_add_sfunc(sam_attributes,"cs",'Z',gt_sam_attribute_generate_cs);
}
//  md  Z  GEM CIGAR String
GT_INLINE gt_status gt_sam_attribute_generate_md(gt_sam_attribute_func_params* func_params) {
  if (func_params->alignment_info->map == NULL) return -1; // Don't print anything
  // Print GEM CIGAR String into the buffer
  gt_output_map_sprint_mismatch_string(func_params->return_s,func_params->alignment_info->map,NULL);
  return 0; // OK
}
GT_INLINE void gt_sam_attributes_add_tag_md(gt_sam_attributes* const sam_attributes) {
  gt_sam_attributes_add_sfunc(sam_attributes,"md",'Z',gt_sam_attribute_generate_md);
}

//  XS  A  +/- directionality infomration for split reads
GT_INLINE gt_status gt_sam_attribute_generate_XS(gt_sam_attribute_func_params* func_params) {
  // do not print this for non split maps
  if (func_params->alignment_info->map==NULL) { // Unmapped
    return -1;
  } else if (func_params->alignment_info->type==GT_MAP_PLACEHOLDER) {
    gt_map* map = func_params->alignment_info->map;
    if(!gt_map_has_next_block(map)){
      return -1; // no split map
    }else{
      bool forward = gt_map_get_strand(map) == FORWARD;
      // is the junction in the overlap ?
      uint64_t junctions_start = gt_map_get_end_mapping_position(forward ? map: gt_map_get_next_block(map));
      uint64_t junctions_end = gt_map_get_begin_mapping_position(forward ? gt_map_get_next_block(map): map);
      gt_sequence_archive* index = func_params->sequence_archive;
      gt_string* j_start = gt_string_new(8);
      gt_string* j_end = gt_string_new(8);
      gt_sequence_archive_get_sequence_string(index, gt_map_get_seq_name(map), gt_map_get_strand(map), junctions_start+1, 5, j_start);
      gt_sequence_archive_get_sequence_string(index, gt_map_get_seq_name(map), gt_map_get_strand(map), junctions_end-1-5, 5, j_end);


      gt_string_delete(j_start);
      gt_string_delete(j_end);
    }
  }
//    // Set proper value to return
//    switch (xt_value) {
//      case GT_XT_UNIQUE:   xt_char_value = 'U'; break;
//      case GT_XT_REPEAT:   xt_char_value = 'R'; break;
//      case GT_XT_UNMAPPED: xt_char_value = 'N'; break;
//      case GT_XT_MATE_SW:  xt_char_value = 'M'; break;
//      default: return -1; break;
//    }
//    // Save as Functional Internal Data (let's save computations)
//    xt_char_value_attr = &xt_char_value;
//    gt_attributes_add(func_params->attributes,GT_ATTR_ID_SAM_TAG_XT,&xt_char_value,char);

  return 0; // OK
}
GT_INLINE void gt_sam_attributes_add_tag_XS(gt_sam_attributes* const sam_attributes) {
  gt_sam_attributes_add_sfunc(sam_attributes,"XS",'A',gt_sam_attribute_generate_XS);
}


//  MQ  i  MAPQ score for mate if paired end alignment
GT_INLINE gt_status gt_sam_attribute_generate_MQ(gt_sam_attribute_func_params* func_params) {
  // only for paired mmaps
  if (func_params->alignment_info->map==NULL) { // Unmapped
    return -1;
  } else if (func_params->alignment_info->type==GT_MMAP_PLACEHOLDER_PAIRED) {
  	func_params->return_i=func_params->alignment_info->paired_end.mate->phred_score;
  	return 0;
  }
  return -1; // NOK
}

GT_INLINE void gt_sam_attributes_add_tag_MQ(gt_sam_attributes* const sam_attributes) {
  gt_sam_attributes_add_ifunc(sam_attributes,"MQ",'i',gt_sam_attribute_generate_MQ);
}


//  UQ  i  PHRED encoded likelihood for all read ends
GT_INLINE gt_status gt_sam_attribute_generate_UQ(gt_sam_attribute_func_params* func_params) {
  // only for paired mmaps
  if (func_params->alignment_info->map==NULL) { // Unmapped
    return -1;
  } else if (func_params->alignment_info->type==GT_MMAP_PLACEHOLDER_PAIRED) {
  	uint64_t sc=func_params->alignment_info->paired_end.mmap_attributes->gt_score;
  	func_params->return_i=(sc&0xffff)+((sc>>16)&0xffff);
  	return 0;
  }
  return -1; // NOK
}

GT_INLINE void gt_sam_attributes_add_tag_UQ(gt_sam_attributes* const sam_attributes) {
  gt_sam_attributes_add_ifunc(sam_attributes,"UQ",'i',gt_sam_attribute_generate_UQ);
}

//  PQ  i  PHRED encoded likelihood for template (includes insert size likelihood term)
GT_INLINE gt_status gt_sam_attribute_generate_PQ(gt_sam_attribute_func_params* func_params) {
  // only for paired mmaps
  if (func_params->alignment_info->map==NULL) { // Unmapped
    return -1;
  } else if (func_params->alignment_info->type==GT_MMAP_PLACEHOLDER_PAIRED) {
  	uint64_t sc=func_params->alignment_info->paired_end.mmap_attributes->gt_score;
  	func_params->return_i=(sc&0xffff)+((sc>>16)&0xffff)+((sc>>32)&0xff);
  	return 0;
  }
  return -1; // OK
}

GT_INLINE void gt_sam_attributes_add_tag_PQ(gt_sam_attributes* const sam_attributes) {
  gt_sam_attributes_add_ifunc(sam_attributes,"PQ",'i',gt_sam_attribute_generate_PQ);
}
