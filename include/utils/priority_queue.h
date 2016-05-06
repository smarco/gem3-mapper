/*
 * PROJECT: GEMMapper
 * FILE: priority_queue.h
 * DATE: 06/06/2012
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
