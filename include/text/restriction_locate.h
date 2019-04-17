/*
 * restriction_locate.h
 *
 *  Created on: 25 Feb 2019
 *      Author: heath
 */

#ifndef RESTRICTION_LOCATE_H_
#define RESTRICTION_LOCATE_H_

#include "utils/essentials.h"
#include "archive/locator.h"

typedef struct {
	vector_t *restriction_site_offsets;
	uint32_t *restriction_site_block_index;
	const locator_interval_t *start_interval;  // In forward strand only (C2T strand for bisufite indexes)
	const locator_interval_t *end_interval;
	// uint64_t start_index;  // Index to interval
	// uint64_t end_index;
	uint64_t num_blocks;
	bool completed;  // Ready to print
} restriction_site_sequence_locator_t;

typedef struct {
	uint64_t forward_text_length;
	restriction_site_sequence_locator_t *restriction_site_sequence_locator;
	uint64_t num_tags;
	bool bisulfite_index;
}	restriction_site_locator_t;


#endif /* RESTRICTION_LOCATE_H_ */
