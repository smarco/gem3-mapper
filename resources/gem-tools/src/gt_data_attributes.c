/*
 * PROJECT: GEM-Tools library
 * FILE: gt_data_attributes.c
 * DATE: 15/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides support to data attributes (from general attributes to specific to certain file formats)
 */


#include "gt_data_attributes.h"
#include "gt_sam_data_attributes.h"

/*
 * General Attribute accessors
 */
GT_INLINE gt_shash* gt_attribute_new(void) {
  return gt_shash_new();
}
GT_INLINE void gt_attribute_clear(gt_shash* const attributes) {
  GT_HASH_CHECK(attributes);
  // Free string SAM attributes
  if (gt_attribute_has_sam_attributes(attributes)) gt_attribute_sam_clear_attributes(attributes);
  // Actual clear of the hash
  gt_shash_clear(attributes,true);
}
GT_INLINE void gt_attribute_delete(gt_shash* const attributes) {
  GT_HASH_CHECK(attributes);
  // Free string attributes
  if (gt_attribute_has_sam_attributes(attributes)) gt_attribute_sam_delete_attributes(attributes);
  // Actual free of the hash
  gt_shash_delete(attributes,true);
}
GT_INLINE void* gt_attribute_get(gt_shash* const attributes,char* const attribute_id) {
  GT_HASH_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  return gt_shash_get(attributes,attribute_id,void);
}
GT_INLINE void gt_attribute_add_string(
    gt_shash* const attributes,char* const attribute_id,gt_string* const attribute_string) {
  GT_HASH_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  GT_STRING_CHECK(attribute_string);
  // Insert attribute
  gt_shash_insert_string(attributes,attribute_id,attribute_string);
}
GT_INLINE void gt_attribute_add_primitive(
    gt_shash* const attributes,char* const attribute_id,void* const attribute,const size_t element_size) {
  GT_HASH_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  GT_NULL_CHECK(attribute);
  GT_ZERO_CHECK(element_size);
  // We do a copy of the element as to handle it ourselves from here
  register void* attribute_cp = gt_malloc(element_size); // Allocate attribute
  memcpy(attribute_cp,attribute,element_size); // Copy attribute
  // Insert attribute
  gt_shash_insert_primitive(attributes,attribute_id,attribute_cp,element_size);
}
GT_INLINE void gt_attribute_add_object(
    gt_shash* const attributes,char* const attribute_id,
    void* const attribute,void* (*attribute_dup_fx)(),void (*attribute_free_fx)()) {
  GT_HASH_CHECK(attributes);
  GT_NULL_CHECK(attribute_id);
  GT_NULL_CHECK(attribute);
  GT_NULL_CHECK(attribute_dup_fx);
  GT_NULL_CHECK(attribute_free_fx);
  // Insert attribute
  gt_shash_insert_object(attributes,attribute_id,attribute,attribute_dup_fx,attribute_free_fx);
}

