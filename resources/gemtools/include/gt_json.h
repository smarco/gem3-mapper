/*
 * PROJECT: GEM-Tools library
 * FILE: gt_json.h
 * DATE: 01/06/2012
 * DESCRIPTION: Helper functions for json output
 */

#ifndef GT_JSON_H_
#define GT_JSON_H_

#include "gt_commons.h"
#include "gt_shash.h"
#include "json.h"

/*
 * Creates a one level json object like:
 *   {
 *     "key": value
 *   }
 * You have to call this with <key>,<value>... pairs in the
 * argument list, i.e:
 *   gt_json_int_named_tuple(2, "A", 1, "B", 2);
 */
GT_INLINE JsonNode* gt_json_int_named_tuple(const uint64_t num_elements,...);

/*
 * Convert the given data into a json array
 */
GT_INLINE JsonNode* gt_json_int_array(const uint64_t start, const uint64_t len, uint64_t* const data);

/*
 * Convert a counter hash from string->uint64_t to a json object
 */
GT_INLINE JsonNode* gt_json_int_hash(gt_shash* const data);


#endif /* GT_JSON_H_ */
