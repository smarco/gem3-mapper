/*
 * PROJECT: GEMMapper
 * FILE: gruntime.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef GRUNTIME_H_
#define GRUNTIME_H_

#include "system/commons.h"

/*
 * Version
 */
#define GEM_CORE_VERSION 3.1
#define GEM_CORE_VERSION_STR "3.1"

/*
 * GEM Runtime
 */
void gruntime_init(
    const uint64_t num_threads,
    const uint64_t max_memory,
    char* const tmp_folder);
void gruntime_destroy();

#endif /* GRUNTIME_H_ */
