/*
 * PROJECT: GEMMapper
 * FILE: bwt.c
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Provides basic routines to encode a DNA text into a Burrows-Wheeler transform
 *              using a compact bit representation and counter buckets as to enhance Occ/rank queries
 */

#include "bwt.h"
#include "profiler.h"

/*
 * Text-sampling
 */
#ifdef SAMPLING_SA_DIRECT
#include "bwt_s1_bm64_sampled.c"
#endif

/*
 * SA-sampling
 */
#ifdef SAMPLING_SA_INVERSE
#include "bwt_s1_bm64.c"
#endif
