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
 *   Region-Profile module provides functions to generate an
 *   fixed/static key partition (with a fixed length size)
 */

#include "filtering/region_profile/region_profile_split.h"
#include "filtering/region_profile/region_profile_schedule.h"
#include "fm_index/fm_index_search.h"

/*
 * Display
 */
void region_profile_split_print_fdec(
    FILE* const stream,
    const uint64_t index,
    const double global_fdec,
    const double* const fdec,
    const uint64_t vector_length,
    const uint64_t offset);
void region_profile_split_print_occ(
    FILE* const stream,
    const uint64_t* const occ,
    const uint64_t vector_length,
    const uint64_t offset);

/*
 * Compute vector of dec factors
 */
void region_profile_split_compute_fdec(
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const bool* const allowed_enc,
    uint64_t* const occ,
    double* const fdec) {
  // Init
  const uint64_t text_length = fm_index_get_length(fm_index);
  uint64_t lo = 0, hi = text_length;
  uint64_t further_occ = text_length, further_pos = key_length;
  int64_t i;
  for (i=key_length-1;i>=0;--i) {
    // Fetch character
    const uint8_t enc_char = key[i];
    if (!allowed_enc[enc_char]) break;
    // Query index
    lo = bwt_sbasic_erank(fm_index->bwt,enc_char,lo);
    hi = bwt_sbasic_erank(fm_index->bwt,enc_char,hi);
    if (hi==lo) break;
    // Store further
    if (hi-lo > 5) {
      further_occ = hi-lo;
      further_pos = i;
    }
  }
  // Return average
  const double dec_prop = (double)text_length/(double)further_occ;
  const double dec_steps = key_length-further_pos;
  *occ = hi-lo;
  *fdec = pow(dec_prop,1./dec_steps);
}
void region_profile_split_fdec(
    FILE* const stream,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const bool* const allowed_enc,
    const uint64_t region_length,
    mm_allocator_t* const mm_allocator) {
  // Allocate fdec vector
  mm_allocator_push_state(mm_allocator);
//  double* const fdec = mm_allocator_calloc(mm_allocator,key_length,double,true);
//  double global_fdec;
//  // Compute global fdec
//  uint64_t length;
//  for (length=region_length;length<=key_length;++length) {
//    uint64_t offset = key_length-length;
//    global_fdec = region_profile_split_compute_fdec(fdec,fm_index,key+offset,length,allowed_enc);
//    region_profile_split_print_fdec(stream,1000+offset,global_fdec,fdec,length,offset);
//  }
//  // Compute partial fdec
//  int64_t position = key_length - region_length;
//  while (position > 0) {
//    const uint64_t begin_pos = position;
//    const uint64_t end_pos = position + region_length;
//    // Compute local fdec
//    global_fdec = region_profile_split_compute_fdec(fdec,
//        fm_index,key+begin_pos,end_pos-begin_pos,allowed_enc);
//    region_profile_split_print_fdec(stream,
//        position,global_fdec,fdec,end_pos-begin_pos,begin_pos);
//    // Next
//    position -= 1;
//  }
  // Free
  mm_allocator_pop_state(mm_allocator);
}
/*
 * Compute splitters
 */
void region_profile_splitters(
    region_profile_t* const region_profile,
    fm_index_t* const fm_index,
    const uint8_t* const key,
    const uint64_t key_length,
    const bool* const allowed_enc,
    const uint64_t sample_length,
    const uint64_t num_regions,
    mm_allocator_t* const mm_allocator) {
  // Init
  region_profile_clear(region_profile);
  region_profile_allocate_regions(region_profile,key_length); // Allocate
  // Allocate fdec vector
  mm_allocator_push_state(mm_allocator);
  const int64_t fdec_length = key_length - sample_length;
  double* const fdec = mm_allocator_calloc(mm_allocator,fdec_length,double,true);
  // Compute fdec vector
  uint64_t occ=UINT64_MAX, prev_occ=UINT64_MAX;
  double accumulated = 0.0;
  int64_t position;
  for (position=0;position+sample_length<=key_length;++position) {
    // Compute Occ
    region_profile_split_compute_fdec(
        fm_index,key+position,sample_length,
        allowed_enc,&occ,&fdec[position]);
    // Platteau condition
    if (occ>500 && occ==prev_occ) {
      fdec[position] = 0.0;
    }
    // Next
    accumulated += fdec[position];
    prev_occ = occ;
  }
  // Init filtering regions
  const uint64_t text_length = fm_index_get_length(fm_index);
  region_search_t* filtering_region = region_profile->filtering_region;
  region_profile->num_filtering_regions = 0;
  filtering_region->begin = 0;
  // Init splitters
  const double split_fdec = accumulated/num_regions;
  // uint64_t h = sample_length; //sample_length/2;
  accumulated = 0.0;
  // Split accumulated
  position = 0;
  fprintf(stderr,"Splitter.acc=%2.3f\n[%02"PRIu64"] ",split_fdec,region_profile->num_filtering_regions);
  while (position < fdec_length) {
    fprintf(stderr,"%2.3f ",fdec[position]);
    accumulated += fdec[position];
    if (accumulated >= split_fdec) {
      // Configure region
      const uint64_t length = position-filtering_region->begin;
      const uint64_t end = length >= sample_length ? position+1 : filtering_region->begin+sample_length;
      filtering_region->end = end;
      fprintf(stderr,"= %2.3f (expected=%2.3f)\n",
          accumulated,(double)text_length/(double)accumulated);
      // Next
      ++(region_profile->num_filtering_regions);
      ++filtering_region;
      filtering_region->begin = end;
      accumulated -= split_fdec;
      fprintf(stderr,"[%02"PRIu64"] ",region_profile->num_filtering_regions);
    }
    ++position;
  }
  fprintf(stderr,"= %2.3f (expected=%2.3f)\n",
      accumulated,(double)text_length/(double)accumulated);
  filtering_region->end = key_length;
  ++(region_profile->num_filtering_regions);
  // Query region profile
  region_profile_query_regions(region_profile,fm_index,key);
  // Schedule to filter all
  region_profile_schedule_exact(region_profile,ALL);
  // Free
  mm_allocator_pop_state(mm_allocator);
}
/*
 * Display
 */
void region_profile_split_print_fdec(
    FILE* const stream,
    const uint64_t index,
    const double global_fdec,
    const double* const fdec,
    const uint64_t vector_length,
    const uint64_t offset) {
  // Show fdec
  uint64_t i;
  fprintf(stream,"[#%04"PRIu64"](g=%05.2f) ",index,global_fdec);
  for (i=0;i<offset;++i) fprintf(stream,"      ");
  for (i=0;i<vector_length;++i) {
    fprintf(stream,"%05.2f ",fdec[i] >= 100.0 ? 99.0 : fdec[i]);
  }
  fprintf(stream,"\n");
}
void region_profile_split_print_occ(
    FILE* const stream,
    const uint64_t* const occ,
    const uint64_t vector_length,
    const uint64_t offset) {
  // Show fdec
  uint64_t i;
  fprintf(stream,"[Occ] ");
  for (i=0;i<offset;++i) fprintf(stream,"      ");
  for (i=0;i<vector_length;++i) {
    fprintf(stream,"%06"PRIu64" ",occ[i]);
  }
  fprintf(stream,"\n");
}



//double region_profile_split_compute_fdec(
//    double* const fdec,
//    fm_index_t* const fm_index,
//    const uint8_t* const key,
//    const uint64_t key_length,
//    const bool* const allowed_enc) {
//  // Init
//  const uint64_t text_length = fm_index_get_length(fm_index);
//  uint64_t lo = 0, hi = text_length;
//  uint64_t further_occ = text_length, further_pos = text_length;
//  int64_t i;
//  for (i=key_length-1;i>=0;--i) {
//    uint64_t next_hi, next_lo;
//    // Fetch character
//    const uint8_t enc_char = key[i];
//    if (!allowed_enc[enc_char]) break;
//    // Query index
//    next_lo = bwt_sbasic_erank(fm_index->bwt,enc_char,lo);
//    next_hi = bwt_sbasic_erank(fm_index->bwt,enc_char,hi);
//    if (next_hi==next_lo) break;
//    // Compute fdec
//    fdec[i] = (double)(hi-lo)/(double)(next_hi-next_lo);
//    lo = next_lo;
//    hi = next_hi;
//    // Store further
//    further_occ = next_hi-next_lo;
//    further_pos = i;
//  }
//  // Complete remaining with 1.0 (const)
//  for (;i>=0;--i) fdec[i] = 1.0;
//  // Return average
//  const double dec_prop = (double)text_length/(double)further_occ;
//  const double dec_steps = key_length-further_pos;
//  return pow(dec_prop,1./dec_steps);
//}
///*
// * Compute splitters (geometric mean)
// */
//void region_profile_splitters(
//    region_profile_t* const region_profile,
//    fm_index_t* const fm_index,
//    const uint8_t* const key,
//    const uint64_t key_length,
//    const bool* const allowed_enc,
//    const uint64_t sample_length,
//    const uint64_t num_regions,
//    mm_allocator_t* const mm_allocator) {
//  // Init
//  region_profile_clear(region_profile);
//  region_profile_allocate_regions(region_profile,key_length); // Allocate
//  // Allocate fdec vector
//  mm_allocator_push_state(mm_allocator);
//  const int64_t fdec_length = key_length - sample_length;
//  double* const fdec = mm_allocator_calloc(mm_allocator,fdec_length,double,true);
//  // Compute fdec vector
//  uint64_t occ=UINT64_MAX, prev_occ=UINT64_MAX;
//  double accumulated = 1.0;
//  int64_t position;
//  for (position=0;position+sample_length<=key_length;++position) {
//    // Compute Occ
//    region_profile_split_compute_fdec(
//        fm_index,key+position,sample_length,
//        allowed_enc,&occ,&fdec[position]);
//    // Platteau condition
//    if (occ>500000 && occ==prev_occ) fdec[position] = 1.0;
//    // Next
//    accumulated *= fdec[position];
//    prev_occ = occ;
//  }
//  // Init filtering regions
//  const uint64_t text_length = fm_index_get_length(fm_index);
//  region_search_t* filtering_region = region_profile->filtering_region;
//  region_profile->num_filtering_regions = 0;
//  filtering_region->begin = 0;
//  // Init splitters
//  const double split_fdec = pow(accumulated,1./(double)num_regions);
//  uint64_t h = sample_length; //sample_length/2;
//  accumulated = 1.0;
//  // Split accumulated
//  position = 0;
//  fprintf(stderr,"Splitter.acc=%2.3f\n[%02lu] ",split_fdec,region_profile->num_filtering_regions);
//  while (position < fdec_length) {
//    fprintf(stderr,"%2.3f ",fdec[position]);
//    accumulated *= fdec[position];
//    if (accumulated >= split_fdec) {
//      // Configure region
//      filtering_region->end = position+h;
//      fprintf(stderr,"= %2.3f (expected=%2.3f)\n",
//          accumulated,(double)text_length/(double)accumulated);
//      // Next
//      ++(region_profile->num_filtering_regions);
//      ++filtering_region;
//      filtering_region->begin = position+h;
//      accumulated = 1.0;
//      fprintf(stderr,"[%02lu] ",region_profile->num_filtering_regions);
//    }
//    ++position;
//  }
//  fprintf(stderr,"= %2.3f (expected=%2.3f)\n",
//      accumulated,(double)text_length/(double)accumulated);
//  filtering_region->end = key_length;
//  ++(region_profile->num_filtering_regions);
//  // Query region profile
//  region_profile_query_regions(region_profile,fm_index,key);
//  // Schedule to filter all
//  region_profile_schedule_exact_all(region_profile);
//  // Free
//  mm_allocator_pop_state(mm_allocator);
//}
