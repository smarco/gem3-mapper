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
 * DESCRIPTION:
 *   General & naive vector sorting generic-routine generator
 *   Define macros VECTOR_SORT_NAME, VECTOR_SORT_TYPE & VECTOR_SORT_CMP to instantiate
 *   each sorting routine for a specific datatype
 *     #define VECTOR_SORT_NAME                 filtering_positions
 *     #define VECTOR_SORT_TYPE                 filtering_position_t
 *     #define VECTOR_SORT_CMP(a,b)             (a->text_end_position - b->text_end_position)
 *     #include "vector_sort.h"
 */

#include "system/commons.h"
#include "utils/vector.h"

#ifndef VECTOR_SORT_H_
#define SORT_QUICKSORT_DEPTH 64
#endif

/*
 * Mandatory type/names defines
 */
#define ___(name,type)                   name ## _ ## type
#define __(name,type)                    ___(name,type)
#define _(name)                          __(name,VECTOR_SORT_NAME)

/*
 * Check sorted
 */
bool _(buffer_sort_check)(
    VECTOR_SORT_TYPE* const vector_mem,
    const int num_positions) {
  int i;
  for (i=1;i<num_positions;++i) {
    if (VECTOR_SORT_CMP(vector_mem + (i-1),vector_mem + i) > 0) {
      fprintf(stderr,"VectorSort: Check failed (not sorted)\n");
      return false;
    }
  }
  return true;
}
/*
 * Selectsort (for small arrays)
 */
void _(vector_sort_selection)(
    VECTOR_SORT_TYPE* const vector_mem,
    const int num_positions) {
  const int top = num_positions-1;
  int i,j;
  for (j=0;j<top;j++) {
    // Find minimum
    int min = j;
    for (i=j+1;i<=top;i++) {
      if (VECTOR_SORT_CMP(vector_mem+i,vector_mem+min) < 0) min = i;
    }
    // Swap
    if (min != j) {
      SWAP(vector_mem[j],vector_mem[min]);
    }
  }
}
/*
 * HeapSort (for deep sorting process)
 */
void _(vector_sort_heapify)(
  VECTOR_SORT_TYPE* const vector_mem,
  const int root,const int bottom) {
  int max_elements = 2*root + 1;
  // Find the biggest child
  if (max_elements < bottom) {
    const int otherChild = max_elements + 1;
    max_elements = (VECTOR_SORT_CMP(vector_mem+otherChild,vector_mem+max_elements) > 0) ? otherChild : max_elements;
  } else if(max_elements > bottom) {
    return;
  }
  // Detect proper ordering
  if (VECTOR_SORT_CMP(vector_mem+root,vector_mem+max_elements) >= 0) return;
  // Swap
  SWAP(vector_mem[root],vector_mem[max_elements]);
  // Heapify
  _(vector_sort_heapify)(vector_mem,max_elements,bottom);
}
void _(vector_sort_heapsort)(
    VECTOR_SORT_TYPE* const vector_mem,
    const int num_positions) {
  // Sort by enforcing heap properties
  int i;
  for (i=num_positions/2;i>=0;i--) {
    _(vector_sort_heapify)(vector_mem,i,num_positions-1); // Heapify
  }
  for (i=num_positions-1;i>=1;i--) {
    SWAP(vector_mem[0],vector_mem[i]); // Swap
    _(vector_sort_heapify)(vector_mem,0,i-1); // Heapify
  }
}
/*
 * Quicksort (for first level sort arrays)
 */
int _(vector_sort_quicksort_partition)(
    VECTOR_SORT_TYPE* const vector_mem,
    const int lo,const int hi) {
  // pivot => hi
  int i=lo-1, j;
  for (j=lo;j<hi;++j) {
    if (VECTOR_SORT_CMP(vector_mem+j,vector_mem+hi) <= 0) {
      i++;
      SWAP(vector_mem[i],vector_mem[j]);
    }
  }
  SWAP(vector_mem[i+1],vector_mem[hi]);
  return i + 1;
}
/*
 * General Sort
 */
void _(buffer_sort)(VECTOR_SORT_TYPE* const vector_mem,const int num_positions) {
  // Create a callstack (queue all the vector elements)
  int quicksort_callstack[SORT_QUICKSORT_DEPTH];
  int top = 2;
  quicksort_callstack[0] = 0;
  quicksort_callstack[1] = num_positions-1;
  // Solve all pending sort-task in the callstack
  while (top > 0) {
    // Pop next vector range
    const int hi = quicksort_callstack[--top];
    const int lo = quicksort_callstack[--top];
    if (hi-lo <= 10) {
      _(vector_sort_selection)(vector_mem+lo,hi-lo+1);
    } else {
      // Set pivot element at its correct position in sorted array
      const int pivot = _(vector_sort_quicksort_partition)(vector_mem,lo,hi);
      // Left side of pivot
      if (pivot - 1 > lo) {
        if (top < SORT_QUICKSORT_DEPTH) {
          quicksort_callstack[top++] = lo;
          quicksort_callstack[top++] = pivot - 1;
        } else {
          _(vector_sort_heapsort)(vector_mem+lo,pivot-lo);
        }
      }
      // Right side of pivot
      if (pivot + 1 < hi) {
        if (top < SORT_QUICKSORT_DEPTH) {
          quicksort_callstack[top++] = pivot + 1;
          quicksort_callstack[top++] = hi;
        } else {
          _(vector_sort_heapsort)(vector_mem+pivot+1,hi-pivot);
        }
      }
    }
  }
#ifdef GEM_DEBUG
  // Check sorting
  _(buffer_sort_check)(vector_mem,num_positions);
#endif
}
void _(vector_sort)(vector_t* const vector) {
  // Parameters
  VECTOR_SORT_TYPE* const vector_mem = vector_get_mem(vector,VECTOR_SORT_TYPE);
  const int num_positions = vector_get_used(vector);
  // Sort
  _(buffer_sort)(vector_mem,num_positions);
}


/*
 * Undefs
 */
#undef VECTOR_SORT_NAME
#undef VECTOR_SORT_TYPE
#undef VECTOR_SORT_CMP
#undef ___
#undef __
#undef _
