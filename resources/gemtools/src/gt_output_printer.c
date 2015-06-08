/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_printer.c
 * DATE: 01/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "gt_output_printer.h"

/*
 * Print Segmented Read Attribute
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS segment_id,total_segments
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output,print_segmented_read_info,const uint64_t segment_id,const uint64_t total_segments);
GT_INLINE gt_status gt_output_gprint_segmented_read_info(gt_generic_printer* const gprinter,const uint64_t segment_id,const uint64_t total_segments) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  gt_gprintf(gprinter," sr:Z:%"PRIu64":%"PRIu64,segment_id,total_segments);
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS segmented_read_info
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output,print_segmented_read,gt_segmented_read_info* const segmented_read_info);
GT_INLINE gt_status gt_output_gprint_segmented_read(gt_generic_printer* const gprinter,gt_segmented_read_info* const segmented_read_info) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_NULL_CHECK(segmented_read_info);
  return gt_output_gprint_segmented_read_info(gprinter,segmented_read_info->segment_id,segmented_read_info->total_segments);
}
/*
 * Print Trim Attributes
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS trim_type,trimmed_read,trimmed_qualities,length
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output,print_trim,
    const gt_read_trim_t trim_type,gt_string* const trimmed_read,gt_string* const trimmed_qualities,const uint64_t length);
GT_INLINE gt_status gt_output_gprint_trim(gt_generic_printer* const gprinter,
    const gt_read_trim_t trim_type,gt_string* const trimmed_read,gt_string* const trimmed_qualities,const uint64_t length) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  const bool has_read = (trimmed_read!=NULL);
  const bool has_qualities = (trimmed_qualities!=NULL);
  if (trim_type==GT_TRIM_RIGHT) {
    if (has_read && has_qualities) {
      gt_gprintf(gprinter," rt:Z:%"PRIu64":"PRIgts":"PRIgts,length,PRIgts_content(trimmed_read),PRIgts_content(trimmed_qualities));
    } else if (!has_read && has_qualities) {
      gt_gprintf(gprinter," rt:Z:%"PRIu64"::"PRIgts,length,PRIgts_content(trimmed_qualities));
    } else if (has_read && !has_qualities) {
      gt_gprintf(gprinter," rt:Z:%"PRIu64":"PRIgts,length,PRIgts_content(trimmed_read));
    } else { // !has_read && !has_qualities
      gt_gprintf(gprinter," rt:Z:%"PRIu64,length);
    }
  } else { // trim_type==GT_TRIM_LEFT
    if (has_read && has_qualities) {
      gt_gprintf(gprinter," lt:Z:%"PRIu64":"PRIgts":"PRIgts,length,PRIgts_content(trimmed_read),PRIgts_content(trimmed_qualities));
    } else if (!has_read && has_qualities) {
      gt_gprintf(gprinter," lt:Z:%"PRIu64"::"PRIgts,length,PRIgts_content(trimmed_qualities));
    } else if (has_read && !has_qualities) {
      gt_gprintf(gprinter," lt:Z:%"PRIu64":"PRIgts,length,PRIgts_content(trimmed_read));
    } else { // !has_read && !has_qualities
      gt_gprintf(gprinter," lt:Z:%"PRIu64,length);
    }
  }
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS right_trim
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output,print_right_trim,gt_read_trim* const right_trim);
GT_INLINE gt_status gt_output_gprint_right_trim(gt_generic_printer* const gprinter,gt_read_trim* const right_trim) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_NULL_CHECK(right_trim);
  return gt_output_gprint_trim(gprinter,GT_TRIM_RIGHT,right_trim->trimmed_read,right_trim->trimmed_qualities,right_trim->length);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS left_trim
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output,print_left_trim,gt_read_trim* const left_trim);
GT_INLINE gt_status gt_output_gprint_left_trim(gt_generic_printer* const gprinter,gt_read_trim* const left_trim) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_NULL_CHECK(left_trim);
  return gt_output_gprint_trim(gprinter,GT_TRIM_LEFT,left_trim->trimmed_read,left_trim->trimmed_qualities,left_trim->length);
}
/*
 * Prints all TAG attributes (trims,sg,casava,extras,...)
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS attributes,print_casava_flags,print_extra_tag_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output,print_tag_attributes,
    gt_attributes* const attributes,const bool print_casava_flags,const bool print_extra_tag_attributes);
GT_INLINE gt_status gt_output_gprint_tag_attributes(gt_generic_printer* const gprinter,
    gt_attributes* const attributes,const bool print_casava_flags,const bool print_extra_tag_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ATTRIBUTES_CHECK(attributes);
  // Check if we have CASAVA 1.8 attributes
  if (print_casava_flags && gt_attributes_is_contained(attributes,GT_ATTR_ID_TAG_CASAVA)) {
    gt_gprintf(gprinter," "PRIgts,PRIgts_content(gt_attributes_get(attributes,GT_ATTR_ID_TAG_CASAVA)));
  } else {
    // Append /1 /2 if paired
    if (gt_attributes_is_contained(attributes,GT_ATTR_ID_TAG_PAIR)) {
      int64_t p = *((int64_t*)gt_attributes_get(attributes,GT_ATTR_ID_TAG_PAIR));
      if (p > 0) gt_gprintf(gprinter,"/%"PRId64,p);
    }
  }
  // Group info
  gt_segmented_read_info* const segmented_read_info = gt_attributes_get_segmented_read_info(attributes);
  if (segmented_read_info!=NULL) gt_output_gprint_segmented_read(gprinter,segmented_read_info);
  // Print RIGHT-Trim
  gt_read_trim* const right_trim = gt_attributes_get_right_trim(attributes);
  if (right_trim!=NULL) gt_output_gprint_right_trim(gprinter,right_trim);
  // Print LEFT-Trim
  gt_read_trim* const left_trim = gt_attributes_get_left_trim(attributes);
  if (left_trim!=NULL) gt_output_gprint_left_trim(gprinter,left_trim);
  // Print additional info (extra tag)
  if (print_extra_tag_attributes && gt_attributes_is_contained(attributes,GT_ATTR_ID_TAG_EXTRA)) {
    gt_gprintf(gprinter," "PRIgts,PRIgts_content(gt_attributes_get(attributes,GT_ATTR_ID_TAG_EXTRA)));
  }
  return 0;
}
