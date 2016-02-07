/*
 * PROJECT: GEMMapper
 * FILE: filtering_region_verify.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef FILTERING_REGION_VERIFY_H_
#define FILTERING_REGION_VERIFY_H_

#include "filtering_region.h"
#include "archive_search_parameters.h"
#include "archive_text.h"
#include "pattern.h"

/*
 * Region Verification
 */
bool filtering_region_verify(
    filtering_region_t* const filtering_region,archive_text_t* const archive_text,
    text_collection_t* const text_collection,pattern_t* const pattern,
    search_parameters_t* const search_parameters,mm_stack_t* const mm_stack);
uint64_t filtering_region_verify_multiple_hits(
    vector_t* const filtering_regions,filtering_region_t* const filtering_region,
    const text_collection_t* const text_collection,search_parameters_t* const search_parameters,
    const pattern_t* const pattern);

uint64_t filtering_region_verify_extension(
    vector_t* const filtering_regions,vector_t* const verified_regions,
    const text_collection_t* const text_collection,const uint64_t text_trace_offset,
    const uint64_t index_position,search_parameters_t* const search_parameters,
    const pattern_t* const pattern);

#endif /* FILTERING_REGION_VERIFY_H_ */
