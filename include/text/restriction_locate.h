/*
 * restriction_locate.h
 *
 *  Created on: 25 Feb 2019
 *      Author: heath
 */

#ifndef RESTRICTION_LOCATE_H_
#define RESTRICTION_LOCATE_H_

#include "utils/essentials.h"

typedef struct {
	vector_t *restriction_site_offsets;
	vector_t *restriction_site_block_index;
	uint64_t *sequence_offsets;
	uint64_t forward_text_length;
	bool bisulfite_index;
}	restriction_site_locator_t;


#endif /* RESTRICTION_LOCATE_H_ */
