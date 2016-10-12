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

#include "utils/priority_queue.h"
#include "system/mm.h"

/*
 * Constants
 */

#define PQUEUE_MIN_INITIAL_ELEMENTS 11

/*
 * Setup
 */
pqueue_t* pqueue_new(uint64_t num_initial_elements) {
  // Alloc
  pqueue_t* pqueue = mm_alloc(pqueue_t);
  // Init
  if (num_initial_elements<PQUEUE_MIN_INITIAL_ELEMENTS) num_initial_elements = PQUEUE_MIN_INITIAL_ELEMENTS;
  pqueue->buffer = vector_new(num_initial_elements,pqueue_element_t);
  vector_set_used(pqueue->buffer,1); // First element never used
  pqueue->num_elements = 1;
  // Return
  return pqueue;
}
void pqueue_clear(pqueue_t* const pqueue) {
  pqueue->num_elements = 1;
}
void pqueue_delete(pqueue_t* const pqueue) {
  // Free buffer
  vector_delete(pqueue->buffer);
  // Free handler
  mm_free(pqueue);
}
/*
 * Accessors
 */
uint64_t pqueue_top_priority(pqueue_t* const pqueue) {
  // Empty case
  if (pqueue->num_elements == 1) return UINT64_MAX;
  pqueue_element_t* const heap = vector_get_mem(pqueue->buffer,pqueue_element_t);
  return heap[1].priority;
}
void* pqueue_top_priority_element(pqueue_t* const pqueue) {
  // Empty case
  if (pqueue->num_elements == 1) return NULL;
  pqueue_element_t* const heap = vector_get_mem(pqueue->buffer,pqueue_element_t);
  return heap[1].element;
}
void pqueue_push_(pqueue_t* const pqueue,void* const element,const uint64_t priority) {
  // Reserve 1 more
  vector_reserve_additional(pqueue->buffer,1);
  vector_inc_used(pqueue->buffer);
  pqueue_element_t* const heap = vector_get_mem(pqueue->buffer,pqueue_element_t);
  // Append at the end & Preserve heap condition (heapify)
  uint64_t n, m;
  n = pqueue->num_elements++;
  while ((m = n/2) && priority < heap[m].priority) {
    heap[n] = heap[m];
    n = m;
  }
  heap[n].element = element;
  heap[n].priority = priority;
}
void* pqueue_pop_(pqueue_t* const pqueue) {
  // Empty case
  if (pqueue->num_elements == 1) return NULL;
  // Shrink by one
  --pqueue->num_elements;
  vector_dec_used(pqueue->buffer);
  // Keep head
  pqueue_element_t* const heap = vector_get_mem(pqueue->buffer,pqueue_element_t);
  void* const top = heap[1].element;
  // Last element goes to the top & Preserve heap condition (heapify)
  const pqueue_element_t* const last = heap + pqueue->num_elements;
  uint64_t n, m;
  n = 1;
  while ((m = 2*n) < pqueue->num_elements) {
    // Left or right ?
    if (m+1 < pqueue->num_elements && heap[m].priority > heap[m+1].priority) m++;
    // Is this the place to be ?
    if (last->priority <= heap[m].priority) break;
    heap[n] = heap[m];
    n = m;
  }
  heap[n] = *last;
  // Return top
  return top;
}

