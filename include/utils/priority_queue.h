/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Simple Priority Queue
 */

#ifndef PRIORITY_QUEUE_H_
#define PRIORITY_QUEUE_H_

#include "system/commons.h"
#include "utils/vector.h"

/*
 * Priority Queue
 */
typedef struct {
  int64_t priority;
  void* element;
} pqueue_element_t;
typedef struct {
  vector_t* buffer;
  uint64_t num_elements;
} pqueue_t;


/*
 * Setup
 */
pqueue_t* pqueue_new(uint64_t num_initial_elements);
void pqueue_clear(pqueue_t* const pqueue);
void pqueue_delete(pqueue_t* const pqueue);

/*
 * Accessors
 */
#define pqueue_get_num_elements(pqueue) ((pqueue)->num_elements-1)
#define pqueue_is_empty(pqueue) (pqueue_get_num_elements(pqueue)==0)

uint64_t pqueue_top_priority(pqueue_t* const pqueue);
void* pqueue_top_priority_element(pqueue_t* const pqueue);

#define pqueue_push(pqueue,element,priority) pqueue_push_(pqueue,(void* const)element,priority)
#define pqueue_pop(pqueue,type) ((type*)pqueue_pop_(pqueue))
void pqueue_push_(pqueue_t* const pqueue,void* const element,const uint64_t priority);
void* pqueue_pop_(pqueue_t* const pqueue);


#endif /* PRIORITY_QUEUE_H_ */
