/*
 * PROJECT: GEMMapper
 * FILE: hash.h
 * DATE: 2/09/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple Hash Adaptor {String,Integer}
 */

#ifndef HASH_H_
#define HASH_H_

#include "commons.h"

/*
 * Checkers
 */
#define HASH_CHECK(hash) GEM_CHECK_NULL(hash)

/*
 * KeyType specific Hashs
 */
#include "ihash.h"
#include "shash.h"

#endif /* HASH_H_ */
