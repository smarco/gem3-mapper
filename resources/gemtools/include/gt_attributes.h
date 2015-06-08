/*
 * PROJECT: GEM-Tools library
 * FILE: gt_data_attributes.h
 * DATE: 15/01/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#ifndef GT_ATTRIBUTES_H_
#define GT_ATTRIBUTES_H_

#include "gt_essentials.h"

/*
 * Attributes IDs
 */
#define GT_ATTR_ID_MAX_COMPLETE_STRATA "MCS"
#define GT_ATTR_ID_NOT_UNIQUE "NOT-UNIQUE"

#define GT_ATTR_ID_TAG_PAIR   "pair"   // (int64_t)
#define GT_ATTR_ID_TAG_CASAVA "casava" // (gt_string)
#define GT_ATTR_ID_TAG_EXTRA  "extra"  // (gt_string)

#define GT_ATTR_ID_LEFT_TRIM  "LTrim"  // (gt_read_trim)
#define GT_ATTR_ID_RIGHT_TRIM "RTrim"  // (gt_read_trim)

#define GT_ATTR_ID_SEGMENTED_READ_INFO "SegmentedReadInfo" // (gt_segmented_read_info)

#define GT_ATTR_ID_SAM_FLAGS "SAM_FLAGS"
#define GT_ATTR_ID_SAM_PRIMARY_ALIGNMENT "SAM_PRIMARY_MAP"
#define GT_ATTR_ID_SAM_PASSING_QC "SAM_PASSING_QC"
#define GT_ATTR_ID_SAM_PCR_DUPLICATE "SAM_PCR_DUPLICATE"
#define GT_ATTR_ID_SAM_ATTRIBUTES "SAM_ATTR"

#define GT_ATTR_ID_SAM_TAG_NH "SAM_NH"
#define GT_ATTR_ID_SAM_TAG_XT "SAM_XT"

/*
 * Attribute Constants
 */
#define GT_PAIR_SE   0
#define GT_PAIR_PE_1 1
#define GT_PAIR_PE_2 2
#define GT_PAIR_BOTH 3

/*
 * Attributes Type
 */
typedef gt_shash gt_attributes;

/*
 * Checkers
 */
#define GT_ATTRIBUTES_CHECK(attributes) GT_HASH_CHECK((gt_shash*)attributes)

/*
 * General Attributes
 */
GT_INLINE gt_attributes* gt_attributes_new(void);
GT_INLINE void gt_attributes_clear(gt_attributes* const attributes);
GT_INLINE void gt_attributes_delete(gt_attributes* const attributes);

GT_INLINE void* gt_attributes_get(gt_attributes* const attributes,char* const attribute_id);
GT_INLINE bool gt_attributes_is_contained(gt_attributes* const attributes,char* const attribute_id);

GT_INLINE void gt_attributes_add_string(
    gt_attributes* const attributes,char* const attribute_id,gt_string* const attribute_string);
GT_INLINE void gt_attributes_add_primitive(
    gt_attributes* const attributes,char* const attribute_id,void* const attribute,const size_t element_size);
GT_INLINE void gt_attributes_add_object(
    gt_attributes* const attributes,char* const attribute_id,
    void* const attribute,void* (*attribute_dup_fx)(),void (*attribute_free_fx)());

#define gt_attributes_add(attributes,attribute_id,attribute,element_type) \
    gt_attributes_add_primitive(attributes,attribute_id,(void*)attribute,sizeof(element_type)) /* use typeof() */

GT_INLINE void gt_attributes_remove(gt_attributes* const attributes,char* const attribute_id);

GT_INLINE gt_attributes* gt_attributes_dup(gt_attributes* const attributes);
GT_INLINE void gt_attributes_copy(gt_attributes* const attributes_dst,gt_attributes* const attributes_src);

/*
 * Specific Attributes
 */

// MAP-format File Attribute
typedef struct {
  bool contains_qualities;
} gt_map_file_format;

// FASTQ/FASTA/MULTIFASTA File Attribute
typedef enum { F_FASTA, F_FASTQ, F_MULTI_FASTA } gt_file_fasta_format;
typedef struct {
  gt_file_fasta_format fasta_format;
} gt_fasta_file_format;

// Read Trim Attribute
typedef enum { GT_TRIM_LEFT, GT_TRIM_RIGHT } gt_read_trim_t;
typedef struct {
  gt_string* trimmed_read;
  gt_string* trimmed_qualities;
  uint64_t length;
} gt_read_trim;

// Segmented Read Attribute
typedef struct {
  uint64_t segment_id;
  uint64_t total_segments;
} gt_segmented_read_info;

/*
 * Specific Attributes Handlers
 */
GT_INLINE gt_read_trim* gt_attributes_get_left_trim(gt_attributes* const attributes);
GT_INLINE gt_read_trim* gt_attributes_get_right_trim(gt_attributes* const attributes);
GT_INLINE void gt_attributes_annotate_left_trim(gt_attributes* const attributes,gt_read_trim* const left_trim);
GT_INLINE void gt_attributes_annotate_right_trim(gt_attributes* const attributes,gt_read_trim* const right_trim);

GT_INLINE gt_segmented_read_info* gt_attributes_get_segmented_read_info(gt_attributes* const attributes);

/*
 * SAM Attributes Handlers
 */
GT_INLINE void gt_attributes_add_sam_ivalue(gt_attributes* const attributes,const char* const tag,const char type_id,const int32_t value);
GT_INLINE void gt_attributes_add_sam_fvalue(gt_attributes* const attributes,const char* const tag,const char type_id,const float value);
GT_INLINE void gt_attributes_add_sam_svalue(gt_attributes* const attributes,const char* const tag,const char type_id,gt_string* const string);

#endif /* GT_ATTRIBUTES_H_ */
