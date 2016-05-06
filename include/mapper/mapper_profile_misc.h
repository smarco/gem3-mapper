/*
 * PROJECT: GEMMapper
 * FILE: mapper_profile_misc.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:  // TODO
 */

#ifndef MAPPER_PROFILE_MISC_H_
#define MAPPER_PROFILE_MISC_H_

#include "utils/essentials.h"
#include "mapper/mapper_profile_counters.h"

/*
 * I/O
 */
void mapper_profile_print_io(FILE* const stream);

/*
 * Output MAP/SAM
 */
void mapper_profile_print_map_output(FILE* const stream,const bool paired_end);
void mapper_profile_print_sam_output(FILE* const stream,const bool paired_end);

/*
 * Checks
 */
void mapper_profile_print_checks(FILE* const stream);

/*
 * Efficiency Ratios
 */
void mapper_profile_print_mapper_efficiency_ratios(FILE* const stream);

#endif /* MAPPER_PROFILE_MISC_H_ */
