/*
 * PROJECT: GEMMapper
 * FILE: archive_text_retrieve.h
 * DATE: 06/06/2013
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 */

#ifndef ARCHIVE_TEXT_RETRIEVE_H_
#define ARCHIVE_TEXT_RETRIEVE_H_

#include "essentials.h"

#include "text_collection.h"
#include "archive.h"

/*
 * Position Locator
 */
GEM_INLINE void archive_text_append_starting_position(
    archive_t* const archive,const uint64_t index_position,const uint64_t right_offset,
    svector_iterator_t* const positions_iterator);

/*
 * Text Retriever
 */
GEM_INLINE void archive_text_retrieve(
    archive_t* const archive,const uint64_t index_position,const uint64_t length,
    text_collection_t* const text_collection);


#endif /* ARCHIVE_TEXT_RETRIEVE_H_ */
