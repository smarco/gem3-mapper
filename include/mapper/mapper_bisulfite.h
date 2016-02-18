/*
 * PROJECT: GEMMapper
 * FILE: mapper_bisulfite.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef MAPPER_BISULFITE_H_
#define MAPPER_BISULFITE_H_

#include "mapper/mapper.h"

/*
 * SE Bisulfite Mapper Thread
 */
void* mapper_SE_bisulfite_thread(mapper_search_t* const mapper_search);

/*
 * PE Bisulfite Mapper Thread
 */
void* mapper_PE_bisulfite_thread(mapper_search_t* const mapper_search);

#endif /* MAPPER_BISULFITE_H_ */
