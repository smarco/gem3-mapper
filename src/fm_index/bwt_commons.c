/*
 * PROJECT: GEMMapper
 * FILE: bwt_commons.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#include "fm_index/bwt_commons.h"

/*
 * Profile
 */
uint64_t _bwt_ranks = 0; // Bwt rank counter

/*
 * XOR table to mask bitmap depending based on the character (enc)
 */
const int64_t xor_table_1[] = {-1ll, -1ll, -1ll, -1ll, -0ll, -0ll, -0ll, -0ll};
const int64_t xor_table_2[] = {-1ll, -1ll, -0ll, -0ll, -1ll, -1ll, -0ll, -0ll};
const int64_t xor_table_3[] = {-1ll, -0ll, -1ll, -0ll, -1ll, -0ll, -1ll, -0ll};
