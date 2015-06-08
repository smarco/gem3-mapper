/*
 * PROJECT: GEM-Tools library
 * FILE: gt_output_map.c
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "gt_output_map.h"
#include "gt_input_map_parser.h"

#define GT_OUTPUT_MAP_COMPACT_COUNTERS_ZEROS_TH 5

/*
 * MAP Printer attributes
 */
GT_INLINE gt_output_map_attributes* gt_output_map_attributes_new() {
	gt_output_map_attributes* attributes = gt_alloc(gt_output_map_attributes);
	gt_output_map_attributes_reset_defaults(attributes);
	return attributes;
}
GT_INLINE void gt_output_map_attributes_delete(gt_output_map_attributes* attributes) {
  GT_NULL_CHECK(attributes);
  gt_free(attributes);
}
GT_INLINE void gt_output_map_attributes_reset_defaults(gt_output_map_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  /* Tag */
	attributes->print_extra = true;
	attributes->print_casava = true;
	/* Counters */
	attributes->compact = false;
	/* Maps */
	attributes->print_scores = true;
	attributes->hex_print_scores = false;
	attributes->max_printable_maps = GT_ALL;
}
/* Tag */
GT_INLINE bool gt_output_map_attributes_is_print_casava(gt_output_map_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  return attributes->print_casava;
}
GT_INLINE void gt_output_map_attributes_set_print_casava(gt_output_map_attributes* const attributes,const bool print_casava) {
  GT_NULL_CHECK(attributes);
  attributes->print_casava = print_casava;
}
GT_INLINE bool gt_output_map_attributes_is_print_extra(gt_output_map_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
  return attributes->print_extra;
}
GT_INLINE void gt_output_map_attributes_set_print_extra(gt_output_map_attributes* const attributes,const bool print_extra) {
  GT_NULL_CHECK(attributes);
  attributes->print_extra = print_extra;
}
/* Maps */
GT_INLINE bool gt_output_map_attributes_is_print_scores(gt_output_map_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
	return attributes->print_scores;
}
GT_INLINE void gt_output_map_attributes_set_print_scores(gt_output_map_attributes* const attributes,const bool print_scores) {
  GT_NULL_CHECK(attributes);
	attributes->print_scores = print_scores;
}
GT_INLINE uint64_t gt_output_map_attributes_get_max_printable_maps(gt_output_map_attributes* const attributes) {
  GT_NULL_CHECK(attributes);
	return attributes->max_printable_maps;
}
GT_INLINE void gt_output_map_attributes_set_max_printable_maps(gt_output_map_attributes* const attributes,const uint64_t max_printable_maps) {
  GT_NULL_CHECK(attributes);
	attributes->max_printable_maps = max_printable_maps;
}
/*
 * TAG building block printers
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS tag,attributes,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_tag,
    gt_string* const tag,gt_attributes* const attributes,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_tag(gt_generic_printer* const gprinter,
    gt_string* const tag,gt_attributes* const attributes,gt_output_map_attributes* const output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_STRING_CHECK(tag);
  GT_ATTRIBUTES_CHECK(attributes);
  // Print the TAG itself
  gt_gprintf(gprinter,PRIgts,PRIgts_content(tag));
  // Print TAG Attributes
  gt_output_gprint_tag_attributes(gprinter,attributes,
      gt_output_map_attributes_is_print_casava(output_map_attributes),
      gt_output_map_attributes_is_print_extra(output_map_attributes));
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_template_tag,gt_template* const template);
GT_INLINE gt_status gt_output_map_gprint_template_tag(
    gt_generic_printer* const gprinter,gt_template* const template) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_output_map_gprint_alignment_tag(gprinter,alignment);
  } GT_TEMPLATE_END_REDUCTION;
  gt_output_map_attributes output_map_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();
  return gt_output_map_gprint_tag(gprinter,template->tag,template->attributes,&output_map_attributes);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_alignment_tag,gt_alignment* const alignment);
GT_INLINE gt_status gt_output_map_gprint_alignment_tag(
    gt_generic_printer* const gprinter,gt_alignment* const alignment) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  gt_output_map_attributes output_map_attributes = GT_OUTPUT_MAP_ATTR_DEFAULT();
  return gt_output_map_gprint_tag(gprinter,alignment->tag,alignment->attributes,&output_map_attributes);
}
GT_INLINE void gt_output_map_gprint_template_reads(
    gt_generic_printer* const gprinter,gt_template* const template,gt_output_map_attributes* output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  GT_OUTPUT_MAP_CHECK_ATTRIBUTES(output_map_attributes);
  // Print READ(s)
  const uint64_t num_blocks = gt_template_get_num_blocks(template);
  uint64_t i = 0;
  gt_gprintf(gprinter,"%s",gt_alignment_get_read(gt_template_get_block(template,i)));
  while (++i<num_blocks) {
    gt_gprintf(gprinter," %s",gt_alignment_get_read(gt_template_get_block(template,i)));
  }
}
GT_INLINE void gt_output_map_gprint_template_qualities(
    gt_generic_printer* const gprinter,gt_template* const template,gt_output_map_attributes* output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  GT_OUTPUT_MAP_CHECK_ATTRIBUTES(output_map_attributes);
  // Print QUALITY
  uint64_t i = 0;
  GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
    if (gt_alignment_has_qualities(alignment)) {
      if (i > 0) {
        gt_gprintf(gprinter," %s",gt_alignment_get_qualities(alignment));
      } else {
        gt_gprintf(gprinter,"%s",gt_alignment_get_qualities(alignment));
      }
    } else if (i > 0) {
      gt_gprintf(gprinter,"");
    }
    ++i;
  }
}
/*
 * Internal MAP printers (take parameters as to control flow/format options)
 */
GT_INLINE gt_status gt_output_map_gprint_mismatch_string_(
    gt_generic_printer* const gprinter,gt_map* const map,gt_output_map_attributes* const output_map_attributes,
    const bool begin_trim,const bool end_trim) {
  GT_NULL_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(output_map_attributes);
  const uint64_t map_length = gt_map_get_base_length(map);
  gt_status error_code = 0;
  uint64_t centinel = 0;
  GT_MISMS_ITERATE(map,misms) {
    const uint64_t misms_pos = gt_misms_get_position(misms);
    if (misms_pos!=centinel) {
      gt_gprintf(gprinter,"%"PRIu64,misms_pos-centinel);
      centinel = misms_pos;
    }
    switch (gt_misms_get_type(misms)) {
      case MISMS:
        gt_gprintf(gprinter,"%c",gt_misms_get_base(misms));
        centinel=misms_pos+1;
        break;
      case INS:
        gt_gprintf(gprinter,">%"PRIu64"+",gt_misms_get_size(misms));
        break;
      case DEL: {
        const uint64_t init_centinel = centinel;
        centinel+=gt_misms_get_size(misms);
        if (gt_expect_false((init_centinel==0 && begin_trim) || (centinel==map_length && end_trim))) { // Trim
          gt_gprintf(gprinter,"(%"PRIu64")",gt_misms_get_size(misms));
        } else {
          gt_gprintf(gprinter,">%"PRIu64"-",gt_misms_get_size(misms));
        }
        break;
      }
      default:
        gt_error(SELECTION_NOT_VALID);
        error_code = GT_MOE_ERROR_PRINTING_MISM_STRING;
        break;
    }
  }
  if (centinel < map_length) {
    gt_gprintf(gprinter,"%"PRIu64,map_length-centinel);
  }
  return error_code;
}
GT_INLINE gt_status gt_output_map_gprint_map_block_(
    gt_generic_printer* const gprinter,gt_map* const map,gt_output_map_attributes* const output_map_attributes,
    const bool begin_trim,const bool end_trim) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(output_map_attributes);
  /*
   * FORMAT => chr11:-:51590050:(5)43T46A9>24*
   */
  // Print sequence name
  gt_gprintf(gprinter,PRIgts,PRIgts_content(gt_map_get_string_seq_name(map)));
  // Print strand
  gt_gprintf(gprinter,GT_MAP_SEP_S"%c",(gt_map_get_strand(map)==FORWARD)?GT_MAP_STRAND_FORWARD_SYMBOL:GT_MAP_STRAND_REVERSE_SYMBOL);
  // Print position
  gt_gprintf(gprinter,GT_MAP_SEP_S"%"PRIu64 GT_MAP_SEP_S,gt_map_get_global_coordinate(map));
  // Print CIGAR
  return gt_output_map_gprint_mismatch_string_(gprinter,map,output_map_attributes,begin_trim,end_trim);
}
GT_INLINE gt_status gt_output_map_gprint_map_(
    gt_generic_printer* const gprinter,gt_map* const map,
    gt_output_map_attributes* const output_map_attributes,const bool print_scores,const bool begin_trim,const bool end_trim) {
  GT_NULL_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(output_map_attributes);
  /*
   * FORMAT => chr11:-:51590050:(5)43T46A9>24*
   */
  gt_status error_code = 0;
  // Print sequence name
  gt_gprintf(gprinter,PRIgts,PRIgts_content(gt_map_get_string_seq_name(map)));
  // Print strand
  gt_gprintf(gprinter,GT_MAP_SEP_S"%c",(gt_map_get_strand(map)==FORWARD)?GT_MAP_STRAND_FORWARD_SYMBOL:GT_MAP_STRAND_REVERSE_SYMBOL);
  // Print position
  gt_gprintf(gprinter,GT_MAP_SEP_S"%"PRIu64 GT_MAP_SEP_S,gt_map_get_global_coordinate(map));
  // Print mismatch string (compact it)
  gt_map* map_it = map;
  gt_map* next_map = NULL;
  bool cigar_pending = true;
  while (cigar_pending) {
    const bool has_next_block = gt_map_has_next_block(map_it);
    error_code|=gt_output_map_gprint_mismatch_string_(gprinter,map_it,output_map_attributes,next_map==NULL,!has_next_block);
    if (has_next_block) {
      next_map = gt_map_get_next_block(map_it);
      const gt_junction_t junction = gt_map_get_junction(map_it);
      if (GT_MAP_IS_SAME_SEGMENT(map_it,next_map) && junction!=QUIMERA) {
        cigar_pending = true;
        switch (junction) {
          case SPLICE:
            gt_gprintf(gprinter,">""%"PRId64"*",gt_map_get_junction_size(map_it));
            break;
          case POSITIVE_SKIP:
            gt_gprintf(gprinter,">""%"PRId64"+",gt_map_get_junction_size(map_it));
            break;
          case NEGATIVE_SKIP:
            gt_gprintf(gprinter,">""%"PRId64"-",gt_map_get_junction_size(map_it));
            break;
          case NO_JUNCTION:
          default:
            error_code=GT_MOE_ERROR_PRINTING_MAP_BLOCKS;
            gt_error(SELECTION_NOT_VALID);
            break;
        }
        map_it = next_map;
      } else {
        cigar_pending = false;
      }
    } else {
      cigar_pending = false;
    }
  }
  // Print quimeras, split-maps across chromosomes, ...
  if (gt_map_has_next_block(map_it)) {
    gt_gprintf(gprinter,GT_MAP_TEMPLATE_SEP);
    error_code|=gt_output_map_gprint_map_(gprinter,next_map,output_map_attributes,false,true,true);
  }
  // Print attributes (scores)
  if (print_scores && gt_map_get_score(map)!=GT_MAP_NO_GT_SCORE) {
  	if (output_map_attributes->hex_print_scores) {
      gt_gprintf(gprinter,GT_MAP_TEMPLATE_SCORE"0x%"PRIx64,gt_map_get_score(map));
  	} else {
  	  gt_gprintf(gprinter,GT_MAP_TEMPLATE_SCORE"%"PRIu64,gt_map_get_score(map));
  	}
  }
  return error_code;
}
GT_INLINE gt_status gt_output_map_gprint_counters_(
    gt_generic_printer* const gprinter,gt_vector* const counters,gt_output_map_attributes* const output_map_attributes,
    const uint64_t max_complete_strata,const bool not_unique_flag) {
  const uint64_t num_counters = gt_vector_get_used(counters);
  uint64_t i;
  // Not unique
  if (not_unique_flag) {
    gt_gprintf(gprinter,GT_MAP_COUNTS_NOT_UNIQUE_S);
    return 0;
  }
  // No counters
  if (num_counters==0) {
    gt_gprintf(gprinter,"0");
    return 0;
  }
  // Print all counters
  for (i=0;i<num_counters;) {
    if (i>0) gt_gprintf(gprinter,"%c",gt_expect_false(i==max_complete_strata)?GT_MAP_MCS:GT_MAP_COUNTS_SEP);
    const uint64_t counter = *gt_vector_get_elm(counters,i,uint64_t);
    if (gt_expect_false(output_map_attributes->compact && counter==0)) {
      uint64_t j=i+1;
      while (j<num_counters && *gt_vector_get_elm(counters,j,uint64_t)==0) ++j;
      if (gt_expect_false((j-i)>=GT_OUTPUT_MAP_COMPACT_COUNTERS_ZEROS_TH)) {
        gt_gprintf(gprinter,"0" GT_MAP_COUNTS_TIMES_S "%"PRIu64,(j-i)); i=j;
      } else {
        gt_gprintf(gprinter,"0"); ++i;
      }
    } else {
      gt_gprintf(gprinter,"%"PRIu64,counter); ++i;
    }
  }
  // MCS (zeros)
  if (max_complete_strata < UINT64_MAX) {
    for (;i<max_complete_strata;++i) {
      if (i>0) {
        gt_gprintf(gprinter,"%c0",GT_MAP_COUNTS_SEP);
      } else {
        gt_gprintf(gprinter,"0");
      }
    }
  }
  return 0;
}
/*
 * MAP building block printers
 *   - If @gt_output_map_attributes==NULL then defaults are applied
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_mismatch_string,gt_map* const map,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_mismatch_string(
    gt_generic_printer* const gprinter,gt_map* const map,gt_output_map_attributes* output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_OUTPUT_MAP_CHECK_ATTRIBUTES(output_map_attributes);
  // Print mismatch string (compact it)
  gt_status error_code = 0;
  gt_map* map_it = map;
  gt_map* next_map = NULL;
  bool cigar_pending = true;
  while (cigar_pending) {
    const bool has_next_block = gt_map_has_next_block(map_it);
    error_code|=gt_output_map_gprint_mismatch_string_(gprinter,map_it,output_map_attributes,next_map==NULL,!has_next_block);
    if (has_next_block) {
      next_map = gt_map_get_next_block(map_it);
      if ((cigar_pending=GT_MAP_IS_SAME_SEGMENT(map_it,next_map))) {
        switch (gt_map_get_junction(map_it)) {
          case SPLICE:
            gt_gprintf(gprinter,">""%"PRIu64"*",gt_map_get_junction_size(map_it));
            break;
          case POSITIVE_SKIP:
            gt_gprintf(gprinter,">""%"PRIu64"+",gt_map_get_junction_size(map_it));
            break;
          case NEGATIVE_SKIP:
            gt_gprintf(gprinter,">""%"PRIu64"-",gt_map_get_junction_size(map_it));
            break;
          case NO_JUNCTION:
          default:
            error_code=GT_MOE_ERROR_PRINTING_MAP_BLOCKS;
            gt_error(SELECTION_NOT_VALID);
            break;
        }
        map_it = next_map;
      }
    } else {
      cigar_pending = false;
    }
  }
  return error_code;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS counters,attributes,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_counters,
    gt_vector* const counters,gt_attributes* const attributes,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_counters(gt_generic_printer* const gprinter,gt_vector* const counters,
    gt_attributes* const attributes,gt_output_map_attributes* output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_VECTOR_CHECK(counters);
  GT_OUTPUT_MAP_CHECK_ATTRIBUTES(output_map_attributes);
  if (attributes!=NULL) {
    uint64_t* const mcs_ptr = (uint64_t*)gt_attributes_get(attributes,GT_ATTR_ID_MAX_COMPLETE_STRATA);
    bool* const not_unique_flag = (bool*)gt_attributes_get(attributes,GT_ATTR_ID_NOT_UNIQUE);
    return gt_output_map_gprint_counters_(gprinter,counters,output_map_attributes,
        ((mcs_ptr!=NULL) ? *mcs_ptr : UINT64_MAX),
        ((not_unique_flag!=NULL) ? *not_unique_flag : false));
  } else {
    return gt_output_map_gprint_counters_(gprinter,counters,output_map_attributes,UINT64_MAX,false);
  }
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_map_block,gt_map* const map,gt_output_map_attributes* output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_map_block(
    gt_generic_printer* const gprinter,gt_map* const map,gt_output_map_attributes* output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_OUTPUT_MAP_CHECK_ATTRIBUTES(output_map_attributes);
  return gt_output_map_gprint_map_block_(gprinter,map,output_map_attributes,true,true);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_map,gt_map* const map,gt_output_map_attributes* output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_map(
    gt_generic_printer* const gprinter,gt_map* const map,gt_output_map_attributes* output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_OUTPUT_MAP_CHECK_ATTRIBUTES(output_map_attributes);
  return gt_output_map_gprint_map_(gprinter,map,output_map_attributes,output_map_attributes->print_scores,true,true);
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS maps,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_map_list,gt_vector* const maps,gt_output_map_attributes* output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_map_list(
    gt_generic_printer* const gprinter,gt_vector* const maps,gt_output_map_attributes* output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_VECTOR_CHECK(maps);
  GT_OUTPUT_MAP_CHECK_ATTRIBUTES(output_map_attributes);
  gt_status error_code = 0;
  GT_VECTOR_ITERATE(maps,map,map_pos,gt_map*) {
    error_code |= gt_output_map_gprint_map_(gprinter,*map,output_map_attributes,output_map_attributes->print_scores,true,true);
  }
  return error_code;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS mmap_placeholder,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_map_placeholder,gt_vector* const mmap_placeholder,gt_output_map_attributes* output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_map_placeholder(gt_generic_printer* const gprinter,
    gt_vector* const mmap_placeholder,gt_output_map_attributes* output_map_attributes) {
  return 0;
  // TODO
//  GT_NULL_CHECK(gprinter);
//  GT_VECTOR_CHECK(mmap_placeholder);
//  GT_OUTPUT_MAP_CHECK_ATTRIBUTES(output_map_attributes);
//  // Determine the number of blocks
//  const uint64_t num_blocks = (template!=NULL) ? gt_template_get_num_blocks(template) : 1;
//  gt_status error_code = 0;
//  if (gt_expect_false(gt_vector_get_used(mmap_placeholder)==0 || output_map_attributes->max_printable_maps==0)) {
//    gt_gprintf(gprinter,GT_MAP_NONE_S);
//  } else {
//    GT_VECTOR_ITERATE(mmap_placeholder,mmap_placeholder,maps_printed,gt_map_placeholder) {
//      if (maps_printed>=output_map_attributes->max_printable_maps) return error_code;
//      if (maps_printed>0) gt_gprintf(gprinter,GT_MAP_NEXT_S);
//      /*
//       * Print MAP
//       */
//      if (mmap_placeholder->type==GT_MAP_PLACEHOLDER || mmap_placeholder->type==GT_MMAP_PLACEHOLDER_UNPAIRED) {
//        gt_alignment* alignment_block = (template!=NULL) ? gt_template_get_block(template,mmap_placeholder->end_position) : alignment;
//        GT_ALIGNMENT_CHECK(alignment_block);
//        // Print preamble
//        uint64_t pos = 0;
//        while (pos < mmap_placeholder->end_position) { gt_gprintf(gprinter,GT_MAP_TEMPLATE_SEP); ++pos; }
//        // Print map
//        gt_map* const map = gt_alignment_get_map(alignment,mmap_placeholder->mmap_position);
//        error_code|=gt_output_map_gprint_map_(gprinter,map,output_map_attributes,output_map_attributes->print_scores,true,true);
//        ++pos;
//        // Print postamble
//        while (pos < num_blocks) { gt_gprintf(gprinter,GT_MAP_TEMPLATE_SEP); ++pos; }
//      } else { // GT_MMAP_PLACEHOLDER_PAIRED
//        GT_TEMPLATE_CHECK(template);
//        gt_mmap_attributes mmap_attributes;
//        gt_map** const mmap = gt_template_get_mmap(template,mmap_placeholder->mmap_position,&mmap_attributes);
//        GT_MMAP_ITERATE_ENDS(mmap,num_blocks,map,end_position) {
//          if (end_position>0) gt_gprintf(gprinter,GT_MAP_TEMPLATE_SEP);
//          error_code|=gt_output_map_gprint_map_(gprinter,map,output_map_attributes,false,true,true);
//        }
//        if (output_map_attributes->print_scores && mmap_attributes.gt_score!=GT_MAP_NO_GT_SCORE) {
//          gt_gprintf(gprinter,GT_MAP_TEMPLATE_SCORE"%"PRIu64,mmap_attributes.gt_score);
//        }
//      }
//    }
//  }
//  return error_code;

}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_template_maps,gt_template* const template,gt_output_map_attributes* output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_template_maps(
    gt_generic_printer* const gprinter,gt_template* const template,gt_output_map_attributes* output_map_attributes) {
  GT_NULL_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  GT_OUTPUT_MAP_CHECK_ATTRIBUTES(output_map_attributes);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_output_map_gprint_alignment_maps(gprinter,alignment,output_map_attributes);
  } GT_TEMPLATE_END_REDUCTION;
  gt_status error_code = 0;
  if (gt_expect_false(gt_template_get_num_mmaps(template)==0 || output_map_attributes->max_printable_maps==0)) {
    gt_gprintf(gprinter,GT_MAP_NONE_S);
  } else {
    const uint64_t num_maps = gt_template_get_num_mmaps(template);
    uint64_t strata = 0, pending_maps = 0, total_maps_printed = 0;
    while (gt_template_get_next_matching_strata(template,strata,&strata,&pending_maps)) {
      GT_TEMPLATE_ITERATE_MMAP__ATTR(template,map_array,map_array_attr) {
        if (map_array_attr->distance!=strata) continue;
        // Print mmap
        --pending_maps;
        if ((total_maps_printed++)>0) gt_gprintf(gprinter,GT_MAP_NEXT_S);
        GT_MMAP_ITERATE(map_array,map,end_position) {
          if (end_position>0) gt_gprintf(gprinter,GT_MAP_TEMPLATE_SEP);
          if (map!=NULL) error_code|=gt_output_map_gprint_map_(gprinter,map,output_map_attributes,false,true,true);
        }
        // Print scores
        if (output_map_attributes->print_scores && map_array_attr!=NULL && map_array_attr->gt_score!=GT_MAP_NO_GT_SCORE) {
        	if(output_map_attributes->hex_print_scores)
            gt_gprintf(gprinter,GT_MAP_TEMPLATE_SCORE"0x%"PRIx64,map_array_attr->gt_score);
        	else gt_gprintf(gprinter,GT_MAP_TEMPLATE_SCORE"%"PRIu64,map_array_attr->gt_score);
        }
        if (total_maps_printed>=output_map_attributes->max_printable_maps || total_maps_printed>=num_maps) return error_code;
        if (pending_maps==0) break;
      }
      if (pending_maps>0) {
        gt_error(TEMPLATE_INCONSISTENT_COUNTERS);
        return GT_MOE_INCONSISTENT_COUNTERS;
      }
      ++strata;
    }
    if (gt_expect_false(total_maps_printed!=num_maps)) {
      gt_error(TEMPLATE_INCONSISTENT_COUNTERS);
      return GT_MOE_INCONSISTENT_COUNTERS;
    }
  }
  return error_code;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_alignment_maps,gt_alignment* const alignment,gt_output_map_attributes* output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_alignment_maps(
    gt_generic_printer* const gprinter,gt_alignment* const alignment,gt_output_map_attributes* output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  GT_OUTPUT_MAP_CHECK_ATTRIBUTES(output_map_attributes);
  gt_status error_code = 0;
  if (gt_expect_false(gt_alignment_get_num_maps(alignment)==0 || output_map_attributes->max_printable_maps==0)) {
    gt_gprintf(gprinter,GT_MAP_NONE_S);
  } else {
    const uint64_t num_maps = gt_alignment_get_num_maps(alignment);
    uint64_t strata = 0, pending_maps = 0, total_maps_printed = 0;
    while (gt_alignment_get_next_matching_strata(alignment,strata,&strata,&pending_maps)) {
      GT_ALIGNMENT_ITERATE(alignment,map) {
        if (gt_map_get_global_distance(map)!=strata) continue;
        // Print map
        --pending_maps;
        if ((total_maps_printed++)>0) gt_gprintf(gprinter,GT_MAP_NEXT_S);
        error_code|=gt_output_map_gprint_map_(gprinter,map,output_map_attributes,output_map_attributes->print_scores,true,true);
        if (total_maps_printed>=output_map_attributes->max_printable_maps || total_maps_printed>=num_maps) return 0;
        if (pending_maps==0) break;
      }
      if (pending_maps > 0) {
        error_code = GT_MOE_INCONSISTENT_COUNTERS;
        gt_error(ALIGNMENT_INCONSISTENT_COUNTERS);
      }
      ++strata;
    }
  }
  return error_code;
}
/*
 * High-level MAP Printers {Alignment/Template}
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_template,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_template(
    gt_generic_printer* const gprinter,gt_template* const template,gt_output_map_attributes* const output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  GT_NULL_CHECK(output_map_attributes);
  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
    return gt_output_map_gprint_alignment(gprinter,alignment,output_map_attributes);
  } GT_TEMPLATE_END_REDUCTION;
  gt_status error_code = 0;
  // Print TAG
  error_code|=gt_output_map_gprint_tag(gprinter,template->tag,template->attributes,output_map_attributes);
  // Print READ(s)
  gt_gprintf(gprinter,"\t");
  gt_output_map_gprint_template_reads(gprinter,template,output_map_attributes);
  // Print QUALITY
  gt_gprintf(gprinter,"\t");
  gt_output_map_gprint_template_qualities(gprinter,template,output_map_attributes);
  // Print COUNTERS
  gt_gprintf(gprinter,"\t");
  error_code|=gt_output_map_gprint_counters_(gprinter,gt_template_get_counters_vector(template),
      output_map_attributes,gt_template_get_mcs(template),gt_template_get_not_unique_flag(template));
  // Print MAPS
  gt_gprintf(gprinter,"\t");
  error_code|=gt_output_map_gprint_template_maps(gprinter,template,output_map_attributes);
  gt_gprintf(gprinter,"\n");
  return error_code;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS alignment,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_alignment,gt_alignment* const alignment,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_alignment(
    gt_generic_printer* const gprinter,gt_alignment* const alignment,gt_output_map_attributes* const output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_ALIGNMENT_CHECK(alignment);
  GT_NULL_CHECK(output_map_attributes);
  gt_status error_code = 0;
  // Print TAG
  error_code|=gt_output_map_gprint_tag(gprinter,alignment->tag,alignment->attributes,output_map_attributes);
  // Print READ(s)
  gt_gprintf(gprinter,"\t%s",gt_alignment_get_read(alignment));
  // Print QUALITY
  if (gt_alignment_has_qualities(alignment)) {
    gt_gprintf(gprinter,"\t%s",gt_alignment_get_qualities(alignment));
  } else {
    gt_gprintf(gprinter,"\t");
  }
  // Print COUNTERS
  gt_gprintf(gprinter,"\t");
  error_code|=gt_output_map_gprint_counters_(gprinter,gt_alignment_get_counters_vector(alignment),
        output_map_attributes,gt_alignment_get_mcs(alignment),gt_alignment_get_not_unique_flag(alignment));
  // Print MAPS
  gt_gprintf(gprinter,"\t");
  error_code|=gt_output_map_gprint_alignment_maps(gprinter,alignment,output_map_attributes);
  gt_gprintf(gprinter,"\n");
  return error_code;
}
/*
 * GEM printer
 *   If the template is paired generates 1 PairedEnd MAP line
 *   If the template is unpaired generates 2 SingleEnd MAP lines
 */
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,output_map_attributes
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_gem_template,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
GT_INLINE gt_status gt_output_map_gprint_gem_template(
    gt_generic_printer* const gprinter,gt_template* const template,gt_output_map_attributes* const output_map_attributes) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_TEMPLATE_CHECK(template);
  if (gt_template_get_num_mmaps(template)>0) {
    return gt_output_map_gprint_template(gprinter,template,output_map_attributes);
  } else {
    gt_status error_code = 0;
    GT_TEMPLATE_ITERATE_ALIGNMENT(template,alignment) {
      if ((error_code|=gt_output_map_gprint_alignment(gprinter,alignment,output_map_attributes))) return error_code;
    }
    return error_code;
  }
}

// TODO
//   this is the code we should aim to...
//#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
//#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS template,output_map_attributes
//GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_template_all,gt_template* const template,gt_output_map_attributes* const output_map_attributes);
//GT_INLINE gt_status gt_output_map_gprint_template_all(
//    gt_generic_printer* const gprinter,gt_template* const template,gt_output_map_attributes* const output_map_attributes) {
//  GT_GENERIC_PRINTER_CHECK(gprinter);
//  GT_TEMPLATE_CHECK(template);
//  GT_TEMPLATE_IF_REDUCES_TO_ALINGMENT(template,alignment) {
//    return gt_output_map_gprint_alignment(gprinter,alignment,output_map_attributes);
//  } GT_TEMPLATE_END_REDUCTION;
//  gt_status error_code = 0;
//  // Print TAG
//  error_code|=gt_output_map_gprint_tag(gprinter,template->tag,template->attributes,output_map_attributes);
//  // Print READ(s)
//  gt_gprintf(gprinter,"\t");
//  gt_output_map_gprint_template_reads(gprinter,template,output_map_attributes);
//  // Print QUALITY
//  gt_gprintf(gprinter,"\t");
//  gt_output_map_gprint_template_qualities(gprinter,template,output_map_attributes);
//  // Build combined structure of paired+unpaired
//  gt_vector* const mmaps_sort = gt_vector_new(gt_template_get_num_mmaps(template),sizeof(gt_map_placeholder));
//  gt_vector* const combined_counters = gt_vector_new(gt_vector_get_used(gt_template_get_counters_vector(template)),sizeof(uint64_t));
//  gt_map_placeholder_create_from_template(template,mmaps_sort,GT_DISTANCE,true,combined_counters);
//  // Print COUNTERS
//  gt_gprintf(gprinter,"\t");
//  error_code|=gt_output_map_gprint_counters_(gprinter,combined_counters,
//      output_map_attributes,gt_template_get_mcs(template),gt_template_get_not_unique_flag(template));
//  // Print MAPS
//  gt_gprintf(gprinter,"\t");
//  error_code|=gt_output_map_gprint_maps_placeholder(gprinter,template,NULL,mmaps_sort,output_map_attributes);
//  gt_gprintf(gprinter,"\n");
//  // Free
//  gt_vector_delete(mmaps_sort);
//  gt_vector_delete(combined_counters);
//  return error_code;
//}


/*
 * Misc. Handy printers
 */
GT_INLINE gt_status gt_output_map_gprint_mismatch_summary_(gt_generic_printer* const gprinter,gt_map* const map,const uint64_t block_num) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  gt_gprintf(gprinter,"{Block %lu}\n",block_num);
  gt_gprintf(gprinter,"  --> GEM-Distance %lu\n",gt_map_get_distance(map));
  gt_gprintf(gprinter,"  --> Levenshtein-Distance %lu\n",gt_map_get_levenshtein_distance(map));
  gt_gprintf(gprinter,"  --> Mismatches %lu\n",gt_map_get_num_mismatch_bases(map));
  gt_gprintf(gprinter,"  --> Indels %lu\n",gt_map_get_num_indels(map));
  gt_gprintf(gprinter,"    --> Deletions %lu\n",gt_map_get_num_deletions(map));
  gt_gprintf(gprinter,"    --> Insertions %lu\n",gt_map_get_num_insertions(map));
  gt_gprintf(gprinter,"  --> Mismatch.list::",block_num);
  GT_MISMS_ITERATE(map,misms) {
    switch (misms->misms_type) {
      case MISMS:
        gt_gprintf(gprinter,"  MIS.%02lu.%c",misms->position,misms->base);
        break;
      case DEL:
        gt_gprintf(gprinter,"  DEL.%02lu<-%02lu>",misms->position,misms->size);
        break;
      case INS:
        gt_gprintf(gprinter,"  INS.%02lu<+%02lu>",misms->position,misms->size);
        break;
    }
  }
  gt_gprintf(gprinter,"\n");
  return 0;
}
#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_mismatch_summary,gt_map* const map);
GT_INLINE gt_status gt_output_map_gprint_mismatch_summary(gt_generic_printer* const gprinter,gt_map* const map) {
  gt_gprintf(gprinter,"[Mismatch/Indel Summary]\n");
  gt_gprintf(gprinter,"  --> Total.Blocks %lu\n",gt_map_get_num_blocks(map));
  uint64_t block_pos = 0;
  GT_MAP_ITERATE(map,map_block) {
    gt_output_map_gprint_mismatch_summary_(gprinter,map_block,block_pos);
    ++block_pos;
  }
  gt_gprintf(gprinter,"  --> Stringify \t ");
  gt_output_map_gprint_map(gprinter,map,NULL);
  gt_gprintf(gprinter,"\n");
  return 0;
}

#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map,pattern,pattern_length,sequence,sequence_length
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_map_block_pretty,gt_map* const map,
    char* const pattern,const uint64_t pattern_length,char* const sequence,const uint64_t sequence_length);
GT_INLINE gt_status gt_output_map_gprint_map_block_pretty(
    gt_generic_printer* const gprinter,gt_map* const map,
    char* const pattern,const uint64_t pattern_length,
    char* const sequence,const uint64_t sequence_length) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_NULL_CHECK(pattern); GT_NULL_CHECK(sequence);
  /*
   * Print the map block (in short)
   */
  gt_gprintf(gprinter,"#{");
  gt_output_map_gprint_map_block(gprinter,map,NULL);
  gt_gprintf(gprinter,"}\n");
  /*
   * Construct the schemes
   */
  const uint64_t approx_bound_length = pattern_length+sequence_length;
  gt_string* const pattern_scheme = gt_string_new(approx_bound_length);
  gt_string* const sequence_scheme = gt_string_new(approx_bound_length);
  gt_string* const relation_scheme = gt_string_new(approx_bound_length);
  uint64_t pattern_centinel=0, sequence_centinel=0;
  // Init misms
  uint64_t misms_offset=0, i;
  const uint64_t num_misms = gt_map_get_num_misms(map);
  gt_misms* misms;
  GT_MAP_CHECK__RELOAD_MISMS_PTR(map,misms_offset,misms,num_misms);
  // Some printing artifacts
  for (i=0;i<4;++i) {
    gt_string_append_char(sequence_scheme,'-');
    gt_string_append_char(relation_scheme,' ');
    gt_string_append_char(pattern_scheme,' ');
  }
  // Traverse the sequence
  bool exception = false, go_on = true;
  while (go_on && (pattern_centinel<pattern_length || sequence_centinel<sequence_length)) {
    if (misms!=NULL && misms->position==pattern_centinel) { // Misms
      switch (misms->misms_type) {
        case MISMS:
          if (pattern_centinel>=pattern_length || sequence_centinel>=sequence_length) {
            exception=true; go_on=false; break;
          }
          gt_string_append_char(sequence_scheme,sequence[sequence_centinel]);
          gt_string_append_char(pattern_scheme,pattern[pattern_centinel]);
          if (pattern[pattern_centinel]==sequence[sequence_centinel]) {
            gt_string_append_char(relation_scheme,'O'); exception=true;
          } else if (misms->base!=sequence[sequence_centinel]) {
            gt_string_append_char(relation_scheme,'#'); exception=true;
          } else {
            gt_string_append_char(relation_scheme,'X');
          }
          ++pattern_centinel; ++sequence_centinel;
          break;
        case INS:
          if (sequence_centinel+misms->size>sequence_length) {
            exception=true; go_on=false; break;
          }
          for (i=0;i<misms->size;++i) {
            gt_string_append_char(sequence_scheme,sequence[sequence_centinel++]);
            gt_string_append_char(relation_scheme,'-');
            gt_string_append_char(pattern_scheme,' ');
          }
          break;
        case DEL:
          if (pattern_centinel+misms->size>pattern_length) {
            exception=true; go_on=false; break;
          }
          for (i=0;i<misms->size;++i) {
            gt_string_append_char(sequence_scheme,' ');
            gt_string_append_char(relation_scheme,'-');
            gt_string_append_char(pattern_scheme,pattern[pattern_centinel++]);
          }
          break;
      }
      ++misms_offset;
      GT_MAP_CHECK__RELOAD_MISMS_PTR(map,misms_offset,misms,num_misms);
    } else { // Match
      if (pattern_centinel>=pattern_length || sequence_centinel>=sequence_length) {
        exception=true; break;
      }
      gt_string_append_char(sequence_scheme,sequence[sequence_centinel]);
      gt_string_append_char(pattern_scheme,pattern[pattern_centinel]);
      if (pattern[pattern_centinel]!=sequence[sequence_centinel]) {
        gt_string_append_char(relation_scheme,'*');
      } else {
        gt_string_append_char(relation_scheme,'|');
      }
      ++pattern_centinel;
      ++sequence_centinel;
    }
  }
  // Fill the rest (in case of exceptions)
  while (pattern_centinel<pattern_length || sequence_centinel<sequence_length) {
    if (sequence_centinel<sequence_length) {
      gt_string_append_char(sequence_scheme,sequence[sequence_centinel]);
    } else {
      gt_string_append_char(sequence_scheme,'?');
    }
    gt_string_append_char(relation_scheme,'!');
    if (pattern_centinel<pattern_length) {
      gt_string_append_char(pattern_scheme,pattern[pattern_centinel]);
    } else {
      gt_string_append_char(pattern_scheme,'?');
    }
    ++pattern_centinel; ++sequence_centinel;
  }
  // Some printing artifacts
  for (i=0;i<4;++i) {
    gt_string_append_char(sequence_scheme,'-');
    gt_string_append_char(relation_scheme,' ');
    gt_string_append_char(pattern_scheme,' ');
  }
  // Append EOS
  gt_string_append_eos(sequence_scheme);
  gt_string_append_eos(relation_scheme);
  gt_string_append_eos(pattern_scheme);
  /*
   * Print the alignment (pretty)
   */
  gt_gprintf(gprinter,PRIgts"\n",PRIgts_content(sequence_scheme));
  gt_gprintf(gprinter,PRIgts"\n",PRIgts_content(relation_scheme));
  gt_gprintf(gprinter,PRIgts"\n",PRIgts_content(pattern_scheme));
  // Free
  gt_string_delete(sequence_scheme);
  gt_string_delete(relation_scheme);
  gt_string_delete(pattern_scheme);
  return exception?1:0;
}

GT_INLINE gt_status gt_output_map_gprint_pretty_map_block_sa(
    gt_generic_printer* const gprinter,gt_map* const map,gt_string* const pattern,
    gt_sequence_archive* const sequence_archive) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Retrieve the sequence
  const uint64_t sequence_length = gt_map_get_length(map);
  gt_string* const sequence = gt_string_new(sequence_length+1);
  gt_status error_code;
  if ((error_code=gt_sequence_archive_retrieve_sequence_chunk(sequence_archive,
      gt_map_get_seq_name(map),gt_map_get_strand(map),gt_map_get_position(map),
      sequence_length,0,sequence))) return error_code;
  // Check Alignment
  return gt_output_map_gprint_map_block_pretty(gprinter,map,
      gt_string_get_string(pattern),gt_string_get_length(pattern),
      gt_string_get_string(sequence),gt_string_get_length(sequence));
}

#undef GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS
#define GT_GENERIC_PRINTER_DELEGATE_CALL_PARAMS map,pattern,sequence_archive
GT_GENERIC_PRINTER_IMPLEMENTATION(gt_output_map,print_map_pretty_sa,
    gt_map* const map,gt_string* const pattern,gt_sequence_archive* const sequence_archive);
GT_INLINE gt_status gt_output_map_gprint_map_pretty_sa(
    gt_generic_printer* const gprinter,gt_map* const map,
    gt_string* const pattern,gt_sequence_archive* const sequence_archive) {
  GT_GENERIC_PRINTER_CHECK(gprinter);
  GT_MAP_CHECK(map);
  GT_STRING_CHECK(pattern);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  // Print map (short)
  const uint64_t num_blocks = gt_map_get_num_blocks(map);
  if (num_blocks > 1) {
    gt_gprintf(gprinter,"SM [ ");
    gt_output_map_gprint_map(gprinter,map,NULL);
    gt_gprintf(gprinter," ] TotalBlocks=%lu\n",gt_map_get_num_blocks(map));
  } else {
    gt_gprintf(gprinter,"RM ");
  }
  // Print all map blocks (Handle SMs)
  gt_status error_code;
  if (num_blocks==1) { // Single-Block
    return gt_output_map_gprint_pretty_map_block_sa(gprinter,map,pattern,sequence_archive);
  } else { // Slit-Maps
    gt_string* read_chunk = gt_string_new(0);
    uint64_t offset = 0, block_num = 0;
    GT_MAP_ITERATE(map,map_block) {
      gt_string_set_nstring(read_chunk,gt_string_get_string(pattern)+offset,gt_map_get_base_length(map_block));
      // Print the mapBlock number
      gt_gprintf(gprinter,"#%lu/%lu",block_num,num_blocks);
      // Print the map block pretty
      if ((error_code=gt_output_map_gprint_pretty_map_block_sa(gprinter,map_block,read_chunk,sequence_archive))) {
        gt_string_delete(read_chunk);
        return error_code;
      }
      offset += gt_map_get_base_length(map_block);
      ++block_num;
    }
    gt_string_delete(read_chunk);
    return 0;
  }
}
