/*
 * PROJECT: GEM-Tools library
 * FILE: gt_attributes.c
 * DATE: 15/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides support to data attributes
 */


#include "gt_attributes.h"
#include "gt_sam_attributes.h"

/*
 * General Attribute accessors
 */
GT_INLINE gt_attributes* gt_attributes_new(void) {
  return gt_shash_new();
}
GT_INLINE void gt_attributes_clear(gt_attributes* const attributes) {
  GT_ATTRIBUTES_CHECK(attributes);
  // Check SAM attributes
  gt_sam_attributes* const sam_attributes = gt_attributes_get_sam_attributes(attributes);
  if (sam_attributes) gt_sam_attributes_clear(sam_attributes);
  // Clear general attributes
  gt_shash_clear(attributes,true);
}
GT_INLINE void gt_attributes_delete(gt_attributes* const attributes) {
  GT_ATTRIBUTES_CHECK(attributes);
  // Check SAM attributes
  gt_attributes_delete_sam_attributes(attributes);
  // Delete general attributes
  gt_shash_delete(attributes,true);
}
GT_INLINE void* gt_attributes_get(gt_attributes* const attributes,char* const attribute_id) {
  GT_ATTRIBUTES_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  return gt_shash_get(attributes,attribute_id,void);
}
GT_INLINE bool gt_attributes_is_contained(gt_attributes* const attributes,char* const attribute_id) {
  GT_ATTRIBUTES_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  return gt_shash_is_contained(attributes,attribute_id);
}
GT_INLINE void gt_attributes_add_string(
    gt_attributes* const attributes,char* const attribute_id,gt_string* const attribute_string) {
  GT_ATTRIBUTES_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  GT_STRING_CHECK(attribute_string);
  // Insert attribute
  gt_shash_insert_string(attributes,attribute_id,attribute_string);
}
GT_INLINE void gt_attributes_add_primitive(
    gt_attributes* const attributes,char* const attribute_id,void* const attribute,const size_t element_size) {
  GT_ATTRIBUTES_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  GT_NULL_CHECK(attribute);
  GT_ZERO_CHECK(element_size);
  // We do a copy of the element as to handle it ourselves from here
  void* attribute_cp = gt_malloc(element_size); // Allocate attribute
  memcpy(attribute_cp,attribute,element_size); // Copy attribute
  // Insert attribute
  gt_shash_insert_primitive(attributes,attribute_id,attribute_cp,element_size);
}
GT_INLINE void gt_attributes_add_object(
    gt_attributes* const attributes,char* const attribute_id,
    void* const attribute,void* (*attribute_dup_fx)(),void (*attribute_free_fx)()) {
  GT_ATTRIBUTES_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  GT_NULL_CHECK(attribute);
  GT_NULL_CHECK(attribute_dup_fx);
  GT_NULL_CHECK(attribute_free_fx);
  // Insert attribute
  gt_shash_insert_object(attributes,attribute_id,attribute,attribute_dup_fx,attribute_free_fx);
}
GT_INLINE void gt_attributes_remove(gt_attributes* const attributes,char* const attribute_id) {
  GT_ATTRIBUTES_CHECK(attributes);
  gt_shash_remove(attributes,attribute_id,true);
}
GT_INLINE gt_attributes* gt_attributes_dup(gt_attributes* const attributes) {
  GT_ATTRIBUTES_CHECK(attributes);
  return gt_shash_dup(attributes);
}
GT_INLINE void gt_attributes_copy(gt_attributes* const attributes_dst,gt_attributes* const attributes_src) {
  GT_ATTRIBUTES_CHECK(attributes_dst);
  GT_ATTRIBUTES_CHECK(attributes_src);
  gt_shash_copy(attributes_dst,attributes_src);
}
/*
 * Specific Attributes Handlers
 */
GT_INLINE gt_read_trim* gt_attributes_get_left_trim(gt_attributes* const attributes) {
  GT_ATTRIBUTES_CHECK(attributes);
  return gt_attributes_get(attributes,GT_ATTR_ID_LEFT_TRIM);
}
GT_INLINE gt_read_trim* gt_attributes_get_right_trim(gt_attributes* const attributes) {
  GT_ATTRIBUTES_CHECK(attributes);
  return gt_attributes_get(attributes,GT_ATTR_ID_RIGHT_TRIM);
}
GT_INLINE void gt_attributes_annotate_left_trim(gt_attributes* const attributes,gt_read_trim* const left_trim) {
  GT_ATTRIBUTES_CHECK(attributes);
  gt_read_trim* const previous_left_trim = gt_attributes_get(attributes,GT_ATTR_ID_LEFT_TRIM);
  if (previous_left_trim==NULL) {
    gt_attributes_add(attributes,GT_ATTR_ID_LEFT_TRIM,left_trim,gt_read_trim); // Set LEFT-Trim
  } else {
    previous_left_trim->length += left_trim->length;
    // Read
    if (left_trim->trimmed_read!=NULL) {
      gt_string_right_append_gt_string(previous_left_trim->trimmed_read,left_trim->trimmed_read);
      gt_string_delete(left_trim->trimmed_read);
    }
    // Qualities
    if (left_trim->trimmed_qualities!=NULL) {
      gt_string_right_append_gt_string(previous_left_trim->trimmed_qualities,left_trim->trimmed_qualities);
      gt_string_delete(left_trim->trimmed_qualities);
    }
  }
}
GT_INLINE void gt_attributes_annotate_right_trim(gt_attributes* const attributes,gt_read_trim* const right_trim) {
  GT_ATTRIBUTES_CHECK(attributes);
  gt_read_trim* const previous_right_trim = gt_attributes_get(attributes,GT_ATTR_ID_RIGHT_TRIM);
  if (previous_right_trim==NULL) {
    gt_attributes_add(attributes,GT_ATTR_ID_RIGHT_TRIM,right_trim,gt_read_trim); // Set RIGHT-Trim
  } else {
    previous_right_trim->length += right_trim->length;
    // Read
    if (right_trim->trimmed_read!=NULL) {
      gt_string_left_append_gt_string(previous_right_trim->trimmed_read,right_trim->trimmed_read);
      gt_string_delete(right_trim->trimmed_read);
    }
    // Qualities
    if (right_trim->trimmed_qualities!=NULL) {
      gt_string_left_append_gt_string(previous_right_trim->trimmed_qualities,right_trim->trimmed_qualities);
      gt_string_delete(right_trim->trimmed_qualities);
    }
  }
}
GT_INLINE gt_segmented_read_info* gt_attributes_get_segmented_read_info(gt_attributes* const attributes) {
  GT_ATTRIBUTES_CHECK(attributes);
  return gt_attributes_get(attributes,GT_ATTR_ID_SEGMENTED_READ_INFO);
}
/*
 * SAM Attributes Handlers
 */
GT_INLINE void gt_attributes_add_sam_ivalue(gt_attributes* const attributes,const char* const tag,const char type_id,const int32_t value) {
  GT_ATTRIBUTES_CHECK(attributes);
  GT_NULL_CHECK(tag);
  gt_sam_attributes_add_ivalue(gt_attributes_get_sam_attributes_dyn(attributes),tag,type_id,value);
}
GT_INLINE void gt_attributes_add_sam_fvalue(gt_attributes* const attributes,const char* const tag,const char type_id,const float value) {
  GT_ATTRIBUTES_CHECK(attributes);
  GT_NULL_CHECK(tag);
  gt_sam_attributes_add_fvalue(gt_attributes_get_sam_attributes_dyn(attributes),tag,type_id,value);
}
GT_INLINE void gt_attributes_add_sam_svalue(gt_attributes* const attributes,const char* const tag,const char type_id,gt_string* const string) {
  GT_ATTRIBUTES_CHECK(attributes);
  GT_NULL_CHECK(tag);
  gt_sam_attributes_add_svalue(gt_attributes_get_sam_attributes_dyn(attributes),tag,type_id,string);
}

