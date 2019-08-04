/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_FMI_SEARCH_C_
#define GPU_FMI_SEARCH_C_

#include "../include/gpu_fmi_core.h"

GPU_INLINE __device__ void gpu_fmi_query_decompose(const char currentBase, uint32_t* const bit0, uint32_t* const bit1, bool* const foundN)
{
  // Decomposing base of the seed in a representative bits
  (* bit0)   =  currentBase & 0x1L;
  (* bit1)   = (currentBase & 0x2L) >> 1;
  (* foundN) = (currentBase & 0x4L) >> 2;
}

GPU_INLINE __device__ void gpu_fmi_query_reverse_lookup(const uint64_t* const query, const uint32_t idBase, const uint32_t querySize,
                                                        char* const globalBase, uint64_t* const globalInfoBase, uint32_t* const globalInfoQuery)
{
  //Query meta information to reduce the main memory requests
  uint64_t infoBase  = (* globalInfoBase);  // Packet 8 bases (cached)
  uint32_t infoQuery = (* globalInfoQuery); // Position cached
  //Address the cached block request (hash function)
  const uint32_t idEntryQuery = (querySize - idBase - 1) / GPU_FMI_BASES_PER_QUERY_ENTRY;
  const uint32_t idIntraQuery = (querySize - idBase - 1) % GPU_FMI_BASES_PER_QUERY_ENTRY;
  uint32_t shiftBits, base;
  //Make the query request if data was not cached
  if((infoQuery != idEntryQuery) || (infoQuery == GPU_UINT32_ONES)){
    infoBase  = LDG(query + idEntryQuery);
    infoQuery = idEntryQuery;
  }
  //logical base extraction from cached request
  shiftBits = idIntraQuery * GPU_FMI_BASE_QUERY_LENGTH;
  base      = (infoBase >> shiftBits) & 0xFF;
  //return the requested data
  (* globalBase)      = base;      //requested base
  (* globalInfoBase)  = infoBase;  //contain query cached data
  (* globalInfoQuery) = infoQuery; //contain the query position cached data
}

GPU_INLINE __device__ void advance_step_LF_mapping(const gpu_fmi_device_entry_t* const fmi, const uint32_t bit0, const uint32_t bit1,
                                                   uint64_t* const L, uint64_t* const R, gpu_fmi_exch_bmp_mem_t * const seedExchBMP)
{
  // Recovering the thread configuration (intra-warp group setup)
  const uint32_t globalThreadIdx     = gpu_get_thread_idx();
  const uint32_t localWarpThreadIdx  = globalThreadIdx    % GPU_WARP_SIZE;
  const uint32_t localEntryIdx       = localWarpThreadIdx / GPU_FMI_THREADS_PER_ENTRY;
  const uint32_t localEntryThreadIdx = localWarpThreadIdx % GPU_FMI_THREADS_PER_ENTRY;
  // Communicate along the FMI entry group threads the L o R interval
  uint64_t       interval       = (localEntryIdx % GPU_FMI_ENTRIES_PER_QUERY) ? (* R) : (* L);
  const uint64_t entryIdx       =  interval / GPU_FMI_ENTRY_SIZE;
  const uint32_t bitmapPosition =  interval % GPU_FMI_ENTRY_SIZE;
  // Loading FM-index entry in thread cooperative way
  const uint32_t missedEntry   = (entryIdx % GPU_FMI_ALTERNATE_COUNTERS != bit1) ? 1 : 0;
  const uint64_t entryIdxFixed = (localEntryThreadIdx == 0) ? entryIdx + missedEntry : entryIdx;
  uint4 loadEntry              = fmi[entryIdxFixed].v[localEntryThreadIdx];
  // Compute LF-Mapping (th0 of each group contain the result)
  interval = LF_Mapping(loadEntry, seedExchBMP, missedEntry, localEntryThreadIdx, bitmapPosition, bit1, bit0);
  // Update interval & communicate the th0 (L, R) to the rest of group threads
  const uint32_t lane = GPU_SELECT_OFFSET(localWarpThreadIdx, GPU_FMI_THREADS_PER_QUERY);
  (* L) = shfl_64(interval, lane);
  (* R) = shfl_64(interval, lane + GPU_FMI_THREADS_PER_ENTRY);
}

void __global__ gpu_fmi_asearch_kernel(const gpu_fmi_device_entry_t* const fmi, const uint64_t bwtSize,
                                       const char* const queries, const uint2* const queryInfo, const uint32_t numQueries,
                                       uint2* const regions, ulonglong2* const regIntervals, uint2* const regOffset,
                                       const uint32_t maxExtraSteps, const uint32_t occThreshold, const uint32_t occShrinkFactor,
                                       const uint32_t maxRegionsFactor)
{
  //Thread group-scheduling initializations
  const uint32_t globalThreadIdx     = gpu_get_thread_idx();
  const uint32_t localWarpThreadIdx  = globalThreadIdx % GPU_WARP_SIZE;
  const uint32_t idQuery             = globalThreadIdx / GPU_FMI_THREADS_PER_QUERY;
  //Just threads with asigned work will execute the search
  if (idQuery < numQueries){
    //Setting thread buffer affinity
    const uint32_t    querySize      = queryInfo[idQuery].y;
    const char* const query          = queries + queryInfo[idQuery].x + querySize - 1; //end of query
    ulonglong2* const regionInterval = regIntervals + regions[idQuery].x;
    uint2* const      regionOffset   = regOffset    + regions[idQuery].x;
    //Shared memory space dedicated for the internal FMI entry thread-communications
    __shared__ gpu_fmi_exch_bmp_mem_t   exchBMP[GPU_FMI_ENTRIES_PER_BLOCK];
               gpu_fmi_exch_bmp_mem_t * const seedExchBMP = &exchBMP[threadIdx.x / GPU_FMI_THREADS_PER_ENTRY];
    //Declarations relative to the search
    const uint32_t maxRegions = GPU_MAX(GPU_DIV_CEIL(querySize, maxRegionsFactor), GPU_FMI_MIN_REGIONS);
          uint32_t idBase = 0, idRegion = 0, initBase = 0, endBase = 0;
          uint64_t L = 0, R = bwtSize, occ = R - L;
          bool     foundN;
    //Extracts and locates each seed
    while ((idBase < querySize) && (idRegion < maxRegions)){
      //Query input initializations
      uint32_t bit0, bit1;
      foundN = false;
      //Search initializations
      L = 0; R = bwtSize; occ = R - L;
      idBase = initBase = endBase;
      //Searching for the next seed
      while((occ > occThreshold) && (idBase < querySize) && !foundN){
        // Gathering the base of the seed
        gpu_fmi_query_decompose(LDG(query - idBase), &bit0, &bit1, &foundN);
        idBase++;
        // Advance step FMI reducing the interval search 
        if(!foundN) advance_step_LF_mapping(fmi, bit0, bit1, &L, &R, seedExchBMP);
        occ = R - L;
      }
      //Evaluate current seed (discard or continue exploration)
      endBase = idBase;
      if(occ <= occThreshold){
        uint64_t endL = L, endR = R;
        if(!foundN){
          // Extension initialization
          uint32_t shrinkOccThreshold = occ >> occShrinkFactor, idStep = 0;
          //Last steps extension (exploration for consecutive 4 bases)
          while((idBase < querySize) && (occ != 0) && (idStep < maxExtraSteps) && !foundN){
            // Gathering the base of the seed
            gpu_fmi_query_decompose(LDG(query - idBase), &bit0, &bit1, &foundN);
            idBase++; idStep++;
            // Advance step FMI reducing the interval search 
            if(!foundN) advance_step_LF_mapping(fmi, bit0, bit1, &L, &R, seedExchBMP);
            // Update seed information
            occ = R - L;
            if((occ < shrinkOccThreshold) && (occ != 0) && !foundN){
              endL = L; endR = R;
              endBase = idBase;
              shrinkOccThreshold = occ;
            }
            shrinkOccThreshold >>= occShrinkFactor;
          }
        }
        // Save extracted region (SA intervals + Query position)
        if((localWarpThreadIdx % GPU_FMI_THREADS_PER_QUERY) == 0){
            regionInterval[idRegion] = make_ulonglong2(endL, endR);
            regionOffset[idRegion]   = make_uint2(querySize - endBase, querySize - initBase);
        }
        idRegion++;
      }
    }
    // Save region profile info (number of extracted regions)
    if((localWarpThreadIdx % GPU_FMI_THREADS_PER_QUERY) == 0){
      if((idRegion == 0) && (querySize == (endBase - initBase))){
        // Save extracted region (SA intervals + Query position)
        regionInterval[idRegion] = make_ulonglong2(L, R);
        regionOffset[idRegion]   = make_uint2(querySize - endBase, querySize - initBase);
        idRegion++;
      }
      regions[idQuery].y = idRegion;
    }
  }
}

GPU_INLINE __device__ void gpu_fmi_table_get_positions(const uint32_t idLevel, const uint32_t idTableLeft, const uint2* const offsetsTable,
                                                       uint32_t* const idGlobalL, uint32_t* const idGlobalR)
{
  // Parameter initialization for the table indexes
  uint32_t idL, idR, idTableRight = idTableLeft >> 2;
  uint2 offset;
  // Gathering the intervals delimiting the table entries for the corresponding level
  offset = LDG(&offsetsTable[idLevel]);
  // Gathering the table entries for L & R (specialized table layout)
  idL = offset.x + idTableLeft; idR = offset.x + idTableLeft + 1;
  // Seeds starting by T, uses the right specialized table layout
  if((idTableLeft & GPU_FMI_TABLE_KEY_MASK) == (GPU_FMI_TABLE_ALPHABET_SIZE - 1))
    idR = offset.y + idTableRight;
  // Return the postition table entry intervals
  (* idGlobalL) = idL; (* idGlobalR) = idR;
}

GPU_INLINE __device__ void gpu_fmi_table_linked_lookup(const uint64_t* const query, const uint32_t querySize, uint64_t* const infoBase, uint32_t* const infoQuery,
                                                       const uint64_t* const fmiTable, const uint2* const offsetsTable, const uint32_t maxLevels,
                                                       uint32_t* const globalBase, uint64_t* const globalL, uint64_t* const globalR)
{
  uint32_t idTable = 0, idLevel = 0, idBase = (* globalBase);
  uint64_t L = (* globalL), R = (* globalR);
  bool     foundN = false;
  // Skipping the first n table levels where n = (OCC > OCC_THRESHOLD)
  while((idLevel < maxLevels - 1) && (idBase < querySize) && !foundN){
    uint32_t bit0, bit1, indexBase;
    char base;
    // Gathering the base of the seed
    gpu_fmi_query_reverse_lookup(query, idBase, querySize, &base, infoBase, infoQuery);
    gpu_fmi_query_decompose(base, &bit0, &bit1, &foundN);
    if(!foundN){
      // Creating the hash key to obtain the corresponding FMI interval
      indexBase = bit0 | (bit1 << 1);
      idTable |= (indexBase << (idLevel << 1));
      idLevel++; idBase++;
    }
  }
  // Query the FMI table if seed no not start with N
  if(idLevel){
    uint32_t idL, idR, idLevelRestored;
    gpu_fmi_table_get_positions(idLevel, idTable, offsetsTable, &idL, &idR);
    L = fmiTable[idL]; R = fmiTable[idR];
    // Gather the truly interval
    idLevelRestored = (L & GPU_FMI_TABLE_LINK_MASK) >> GPU_FMI_TABLE_FIELD_LENGTH;
    if(idLevelRestored != idLevel){
      uint32_t idL, idR;
      // Restoring the truly hash table key
      idTable &= ~(GPU_UINT32_ONES << (idLevelRestored << 1));
      gpu_fmi_table_get_positions(idLevelRestored, idTable, offsetsTable, &idL, &idR);
      L = fmiTable[idL]; R = fmiTable[idR];
      idBase -= (idLevel - idLevelRestored);
    }
  }
  // Updating the search parameters
  (* globalL)      = L & GPU_FMI_TABLE_FIELD_MASK;
  (* globalR)      = R & GPU_FMI_TABLE_FIELD_MASK;
  (* globalBase)   = idBase;
}

GPU_INLINE __device__ void gpu_fmi_table_lookup(const uint64_t* const query, const uint32_t querySize, uint32_t* const globalBase,
                                                uint64_t* const infoBase, uint32_t* const infoQuery,
                                                const uint64_t* const fmiTable, const uint2* const offsetsTable,
                                                const uint32_t skipLevels, const uint32_t maxLevels, const uint32_t occThreshold,
                                                uint64_t* const globalL, uint64_t* const globalR, bool* const globalFoundN)
{
  uint32_t idTableLeft = 0, idTableRight = 0, idLevel = 0, idBase = (* globalBase);
  uint64_t L = (* globalL), R = (* globalR), occ = (* globalR) - (* globalL);
  bool     foundN = (* globalFoundN), rightTableProcessing = false;
  uint2    offset;
  // Seeds starting by T, uses the right specialized table layout
  if(idBase < querySize){
    char base;
    gpu_fmi_query_reverse_lookup(query, idBase, querySize, &base, infoBase, infoQuery);
    if(base == (GPU_FMI_TABLE_ALPHABET_SIZE - 1)) rightTableProcessing = true;
  }
  // Skipping the first n table levels where n = (OCC > OCC_THRESHOLD)
  while((idLevel < skipLevels) && (idBase < querySize) && !foundN){
    uint32_t bit0, bit1, indexBase;
    char base;
    // Gathering the base of the seed
    gpu_fmi_query_reverse_lookup(query, idBase, querySize, &base, infoBase, infoQuery);
    gpu_fmi_query_decompose(base, &bit0, &bit1, &foundN);
    idBase++;
    if(!foundN){
      // Creating the hash key to obtain the corresponding FMI interval
      indexBase = bit0 | (bit1 << 1);
      idTableLeft |= (indexBase << (idLevel << 1));
      idTableRight = idTableLeft >> 2;
      idLevel++;
    }
  }
  // Query the FMI table if seed no not start with N
  if(!foundN && idLevel){
    uint32_t idL, idR;
    // Gathering the intervals delimiting the table entries for the corresponding level
    offset = LDG(&offsetsTable[idLevel]);
    // Gathering the table entries for L & R (specialized table layout)
    idL = offset.x + idTableLeft; idR = offset.x + idTableLeft + 1;
    if(rightTableProcessing) idR = offset.y + idTableRight;
    L = fmiTable[idL] & GPU_FMI_TABLE_FIELD_MASK; R = fmiTable[idR] & GPU_FMI_TABLE_FIELD_MASK;
    occ = R - L;
  }
  // Extending the FMI table query to fit in the adatative search requirements
  while((idLevel < maxLevels-1) && (occ > occThreshold) && (idBase < querySize) && !foundN){
    uint32_t bit0, bit1, idL, idR, indexBase;
    char base;
    // Gathering the base of the seed
    gpu_fmi_query_reverse_lookup(query, idBase, querySize, &base, infoBase, infoQuery);
    gpu_fmi_query_decompose(base, &bit0, &bit1, &foundN);
    idBase++;
    // Look-up the table (do not contain Ns intervals)
    if(!foundN){
      // Creating the hash key to obtain the corresponding FMI interval
      indexBase = bit0 | (bit1 << 1);
      idTableLeft |= (indexBase << (idLevel << 1));
      idTableRight = idTableLeft >> 2;
      idLevel++;
      // Gathering the intervals delimiting the table entries for the corresponding level
      offset = LDG(&offsetsTable[idLevel]);
      // Gathering the table entries for L & R (specialized table layout)
      idL = offset.x + idTableLeft; idR = offset.x + idTableLeft + 1;
      if(rightTableProcessing) idR = offset.y + idTableRight;
      L = fmiTable[idL] & GPU_FMI_TABLE_FIELD_MASK; R = fmiTable[idR] & GPU_FMI_TABLE_FIELD_MASK;
      occ = R - L;
    }
  }
  // Updating the search parameters
  (* globalL)      = L;
  (* globalR)      = R;
  (* globalFoundN) = foundN;
  (* globalBase)   = idBase;
}

void __global__ gpu_fmi_asearch_table_linked_kernel(const gpu_fmi_device_entry_t* const fmi, const uint64_t bwtSize,
                                                    const uint64_t* const fmiTable, const uint2* const offsetsTable, const uint32_t maxLevels,
                                                    const char* const queries, const uint2* const queryInfo, const uint32_t numQueries,
                                                    uint2* const regions, ulonglong2* const regIntervals, uint2* const regOffset,
                                                    const uint32_t maxExtraSteps, const uint32_t occThreshold, const uint32_t occShrinkFactor,
                                                    const uint32_t maxRegionsFactor)
{
  // Thread group-scheduling initializations
  const uint32_t globalThreadIdx     = gpu_get_thread_idx();
  const uint32_t localWarpThreadIdx  = globalThreadIdx % GPU_WARP_SIZE;
  const uint32_t idQuery             = globalThreadIdx / GPU_FMI_THREADS_PER_QUERY;
  // Sanity check (threads with asigned work will be executed)
  if (idQuery < numQueries){
    // Setting thread buffer affinity
    const uint32_t        querySize      = queryInfo[idQuery].y; // 1st query position
    const uint64_t* const query          = (uint64_t*) (queries + queryInfo[idQuery].x);
    ulonglong2* const     regionInterval = regIntervals + regions[idQuery].x;
    uint2* const          regionOffset   = regOffset    + regions[idQuery].x;
    // Shared memory space dedicated for the internal FMI entry thread-communications
    __shared__ gpu_fmi_exch_bmp_mem_t   exchBMP[GPU_FMI_ENTRIES_PER_BLOCK];
               gpu_fmi_exch_bmp_mem_t * const seedExchBMP = &exchBMP[threadIdx.x / GPU_FMI_THREADS_PER_ENTRY];
    //Search variable declaration
    const uint32_t maxRegions = GPU_MAX(GPU_DIV_CEIL(querySize, maxRegionsFactor), GPU_FMI_MIN_REGIONS);
          uint32_t infoQuery = GPU_UINT32_ONES, idBase = 0, idRegion = 0, initBase = 0, endBase = 0;
          uint64_t infoBase = 0, L = 0, R = bwtSize, occ = R - L;
          bool     foundN;
          char     base;
    // Extracts and locates each seed
    while ((idBase < querySize) && (idRegion < maxRegions)){
      // Query input initializations
      uint32_t bit0, bit1;
      foundN = false;
      // Search initializations
      L = 0; R = bwtSize;
      idBase = initBase = endBase;
      // LUT initializations
      gpu_fmi_table_linked_lookup(query, querySize, &infoBase, &infoQuery, fmiTable, offsetsTable, maxLevels, &idBase, &L, &R);
      occ = R - L;
      //Searching for the next seed
      while((occ > occThreshold) && (idBase < querySize) && !foundN){
        // Gathering the base of the seed
        gpu_fmi_query_reverse_lookup(query, idBase, querySize, &base, &infoBase, &infoQuery);
        gpu_fmi_query_decompose(base, &bit0, &bit1, &foundN);
        idBase++;
        // Advance step FMI reducing the interval search
        if(!foundN) advance_step_LF_mapping(fmi, bit0, bit1, &L, &R, seedExchBMP);
        occ = R - L;
      }
      //Evaluate current seed (discard or continue exploration)
      endBase = idBase;
      if(occ <= occThreshold){
        uint64_t endL = L, endR = R;
        if(!foundN){
          // Extension initialization
          uint32_t shrinkOccThreshold = occ >> occShrinkFactor, idStep = 0;
          //Last steps extension (exploration for consecutive 4 bases)
          while((idBase < querySize) && (occ != 0) && (idStep < maxExtraSteps) && !foundN){
            // Gathering the base of the seed
            gpu_fmi_query_reverse_lookup(query, idBase, querySize, &base, &infoBase, &infoQuery);
            gpu_fmi_query_decompose(base, &bit0, &bit1, &foundN);
            idBase++; idStep++;
            // Advance step FMI reducing the interval search
            if(!foundN) advance_step_LF_mapping(fmi, bit0, bit1, &L, &R, seedExchBMP);
            // Update seed information
            occ = R - L;
            if((occ < shrinkOccThreshold) && (occ != 0) && !foundN){
              endL = L; endR = R;
              endBase = idBase;
              shrinkOccThreshold = occ;
            }
            shrinkOccThreshold >>= occShrinkFactor;
          }
        }
        // Save extracted region (SA intervals + Query position)
        if((localWarpThreadIdx % GPU_FMI_THREADS_PER_QUERY) == 0){
            regionInterval[idRegion] = make_ulonglong2(endL, endR);
            regionOffset[idRegion]   = make_uint2(querySize - endBase, querySize - initBase);
        }
        idRegion++;
      }
    }
    // Save region profile info (number of extracted regions)
    if((localWarpThreadIdx % GPU_FMI_THREADS_PER_QUERY) == 0){
      if((idRegion == 0) && (querySize == (endBase - initBase))){
        // Save extracted region (SA intervals + Query position)
        regionInterval[idRegion] = make_ulonglong2(L, R);
        regionOffset[idRegion]   = make_uint2(querySize - endBase, querySize - initBase);
        idRegion++;
      }
      regions[idQuery].y = idRegion;
    }
  }
}

void __global__ gpu_fmi_asearch_table_kernel(const gpu_fmi_device_entry_t* const fmi, const uint64_t bwtSize,
                                             const uint64_t* const fmiTable, const uint2* const offsetsTable, const uint32_t skipLevels, const uint32_t maxLevels,
                                             const char* const queries, const uint2* const queryInfo, const uint32_t numQueries,
                                             uint2* const regions, ulonglong2* const regIntervals, uint2* const regOffset,
                                             const uint32_t maxExtraSteps, const uint32_t occThreshold, const uint32_t occShrinkFactor,
                                             const uint32_t maxRegionsFactor)
{
  // Thread group-scheduling initializations
  const uint32_t globalThreadIdx     = gpu_get_thread_idx();
  const uint32_t localWarpThreadIdx  = globalThreadIdx % GPU_WARP_SIZE;
  const uint32_t idQuery             = globalThreadIdx / GPU_FMI_THREADS_PER_QUERY;
  // Sanity check (threads with asigned work will be executed)
  if (idQuery < numQueries){
    // Setting thread buffer affinity
    const uint32_t        querySize      = queryInfo[idQuery].y; // 1st query position
    const uint64_t* const query          = (uint64_t*) (queries + queryInfo[idQuery].x);
    ulonglong2* const     regionInterval = regIntervals + regions[idQuery].x;
    uint2* const          regionOffset   = regOffset    + regions[idQuery].x;
    // Shared memory space dedicated for the internal FMI entry thread-communications
    __shared__ gpu_fmi_exch_bmp_mem_t   exchBMP[GPU_FMI_ENTRIES_PER_BLOCK];
               gpu_fmi_exch_bmp_mem_t * const seedExchBMP = &exchBMP[threadIdx.x / GPU_FMI_THREADS_PER_ENTRY];
    //Search variable declaration
    const uint32_t maxRegions = GPU_MAX(GPU_DIV_CEIL(querySize, maxRegionsFactor), GPU_FMI_MIN_REGIONS);
          uint32_t infoQuery = GPU_UINT32_ONES, idBase = 0, idRegion = 0, initBase = 0, endBase = 0;
          uint64_t infoBase = 0, L = 0, R = bwtSize, occ = R - L;
          bool     foundN;
          char     base;
    // Extracts and locates each seed
    while ((idBase < querySize) && (idRegion < maxRegions)){
      // Query input initializations
      uint32_t bit0, bit1;
      foundN = false;
      // Search initializations
      L = 0; R = bwtSize;
      idBase = initBase = endBase;
      // LUT initializations
      gpu_fmi_table_lookup(query, querySize, &idBase, &infoBase, &infoQuery,
                           fmiTable, offsetsTable, skipLevels, maxLevels,
                           occThreshold, &L, &R, &foundN);
      occ = R - L;
      //Searching for the next seed
      while((occ > occThreshold) && (idBase < querySize) && !foundN){
        // Gathering the base of the seed
        gpu_fmi_query_reverse_lookup(query, idBase, querySize, &base, &infoBase, &infoQuery);
        gpu_fmi_query_decompose(base, &bit0, &bit1, &foundN);
        idBase++;
        // Advance step FMI reducing the interval search
        if(!foundN) advance_step_LF_mapping(fmi, bit0, bit1, &L, &R, seedExchBMP);
        occ = R - L;
      }
      //Evaluate current seed (discard or continue exploration)
      endBase = idBase;
      if(occ <= occThreshold){
        uint64_t endL = L, endR = R;
        if(!foundN){
          // Extension initialization
          uint32_t shrinkOccThreshold = occ >> occShrinkFactor, idStep = 0;
          //Last steps extension (exploration for consecutive 4 bases)
          while((idBase < querySize) && (occ != 0) && (idStep < maxExtraSteps) && !foundN){
            // Gathering the base of the seed
            gpu_fmi_query_reverse_lookup(query, idBase, querySize, &base, &infoBase, &infoQuery);
            gpu_fmi_query_decompose(base, &bit0, &bit1, &foundN);
            idBase++; idStep++;
            // Advance step FMI reducing the interval search
            if(!foundN) advance_step_LF_mapping(fmi, bit0, bit1, &L, &R, seedExchBMP);
            // Update seed information
            occ = R - L;
            if((occ < shrinkOccThreshold) && (occ != 0) && !foundN){
              endL = L; endR = R;
              endBase = idBase;
              shrinkOccThreshold = occ;
            }
            shrinkOccThreshold >>= occShrinkFactor;
          }
        }
        // Save extracted region (SA intervals + Query position)
        if((localWarpThreadIdx % GPU_FMI_THREADS_PER_QUERY) == 0){
            regionInterval[idRegion] = make_ulonglong2(endL, endR);
            regionOffset[idRegion]   = make_uint2(querySize - endBase, querySize - initBase);
        }
        idRegion++;
      }
    }
    // Save region profile info (number of extracted regions)
    if((localWarpThreadIdx % GPU_FMI_THREADS_PER_QUERY) == 0){
      if((idRegion == 0) && (querySize == (endBase - initBase))){
        // Save extracted region (SA intervals + Query position)
        regionInterval[idRegion] = make_ulonglong2(L, R);
        regionOffset[idRegion]   = make_uint2(querySize - endBase, querySize - initBase);
        idRegion++;
      }
      regions[idQuery].y = idRegion;
    }
  }
}

extern "C"
gpu_error_t gpu_fmi_asearch_launch_kernel(const gpu_fmi_device_entry_t* const d_fmi, const uint64_t bwtSize, const gpu_fmi_table_format_t formatTable,
                                          const uint64_t* const d_fmiTable, const uint2* const d_offsetsTable, const uint32_t skipLevels, const uint32_t maxLevels,
                                          const char* const d_queries, const uint2* const d_queryInfo, const uint32_t numQueries,
                                          uint2* const d_regions, ulonglong2* const d_regIntervals, uint2* const d_regOffset)
{
  //GPU thread configuration
  const uint32_t threads = 128;
  const uint32_t blocks  = GPU_DIV_CEIL(numQueries * GPU_FMI_THREADS_PER_QUERY, threads);
  const uint32_t nreps   = 10;
  //Search configuration
  const uint32_t occShrinkFactor  =  GPU_FMI_D;
  const uint32_t occThreshold     =  GPU_FMI_OCC_THRESHOLD;
  const uint32_t maxExtraSteps    =  GPU_FMI_STEPS;
  const uint32_t maxRegionsFactor =  GPU_FMI_RATIO_REGIONS;

  float elapsed_time_ms = 0.0f;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

    for(uint32_t iteration = 0; iteration < nreps; ++iteration){
      switch(formatTable){
        case GPU_FMI_TABLE_DISABLED:
          gpu_fmi_asearch_kernel<<<blocks,threads>>>(d_fmi, bwtSize, d_queries, d_queryInfo, numQueries,
                                                     d_regions, d_regIntervals, d_regOffset,
                                                     maxExtraSteps, occThreshold, occShrinkFactor, maxRegionsFactor);
          break;
        case GPU_FMI_TABLE_MULTILEVEL:
          gpu_fmi_asearch_table_kernel<<<blocks,threads>>>(d_fmi, bwtSize, d_fmiTable, d_offsetsTable, skipLevels, maxLevels,
                                                           d_queries, d_queryInfo, numQueries,
                                                           d_regions, d_regIntervals, d_regOffset,
                                                           maxExtraSteps, occThreshold, occShrinkFactor, maxRegionsFactor);
        break;
        case GPU_FMI_TABLE_MULTILEVEL_LINKED:
          gpu_fmi_asearch_table_linked_kernel<<<blocks,threads>>>(d_fmi, bwtSize, d_fmiTable, d_offsetsTable, maxLevels,
                                                                  d_queries, d_queryInfo, numQueries,
                                                                  d_regions, d_regIntervals, d_regOffset,
                                                                  maxExtraSteps, occThreshold, occShrinkFactor, maxRegionsFactor);
        break;
      }
    }

  cudaEventRecord(stop, 0);
  cudaThreadSynchronize();
  cudaEventElapsedTime(&elapsed_time_ms, start, stop);
  elapsed_time_ms /= nreps;

  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  printf("\t Time Kernel GPU:  %8.2f ms\n", elapsed_time_ms);

  return(SUCCESS);
}

extern "C"
gpu_error_t gpu_fmi_asearch_process_buffer(gpu_buffer_t* const mBuff)
{
  // Getting internal buffers
  const gpu_index_buffer_t* const           index            =  mBuff->index;
  const gpu_fmi_table_t* const              table            = &mBuff->index->fmi.table;
  const gpu_fmi_asearch_queries_buffer_t*   queries          = &mBuff->data.asearch.queries;
  const gpu_fmi_asearch_regions_buffer_t*   regions          = &mBuff->data.asearch.regions;

  // Getting buffer sizes
  const uint32_t                            numQueries       =  mBuff->data.asearch.queries.numQueries;
  const uint32_t                            numMaxQueries    =  mBuff->data.asearch.numMaxQueries;
  const uint32_t                            numBases         =  mBuff->data.asearch.queries.numBases;
  const uint32_t                            numMaxBases      =  mBuff->data.asearch.numMaxBases;
  // Getting device information
  const cudaStream_t                        idStream         =  mBuff->listStreams[mBuff->idStream];
  const uint32_t                            idSupDev         =  mBuff->idSupportedDevice;
  const gpu_device_info_t* const            device           =  mBuff->device[idSupDev];
  // Search configuration
  const uint32_t                            occShrinkFactor  =  mBuff->data.asearch.occShrinkFactor;
  const uint32_t                            occThreshold     =  mBuff->data.asearch.occMinThreshold;
  const uint32_t                            maxExtraSteps    =  mBuff->data.asearch.extraSteps;
  const uint32_t                            maxRegionsFactor =  mBuff->data.asearch.maxRegionsFactor;

  dim3 blocksPerGrid, threadsPerBlock;
  const uint32_t numThreads = numQueries * GPU_FMI_THREADS_PER_QUERY;
  gpu_device_kernel_thread_configuration(device, numThreads, &blocksPerGrid, &threadsPerBlock);
  // Sanity-check (checks buffer overflowing)
  if((numQueries > numMaxQueries) || (numBases > numMaxBases))
    return(E_OVERFLOWING_BUFFER);

  switch(table->formatTableLUT){
    case GPU_FMI_TABLE_DISABLED:
      gpu_fmi_asearch_kernel<<<blocksPerGrid, threadsPerBlock, 0, idStream>>>((gpu_fmi_device_entry_t*) index->fmi.d_fmi[idSupDev], index->fmi.bwtSize,
                                                                              (char*) queries->d_queries, (uint2*) queries->d_queryInfo, numQueries,
                                                                              (uint2*) queries->d_regions, (ulonglong2*) regions->d_intervals, (uint2*) regions->d_regionsOffsets,
                                                                              maxExtraSteps, occThreshold, occShrinkFactor, maxRegionsFactor);

    break;
    case GPU_FMI_TABLE_MULTILEVEL:
      gpu_fmi_asearch_table_kernel<<<blocksPerGrid, threadsPerBlock, 0, idStream>>>((gpu_fmi_device_entry_t*) index->fmi.d_fmi[idSupDev], index->fmi.bwtSize,
                                                                                    (uint64_t*) table->d_fmiTableLUT[idSupDev], (uint2*) table->d_offsetsTableLUT[idSupDev], table->skipLevelsTableLUT, table->maxLevelsTableLUT,
                                                                                    (char*) queries->d_queries, (uint2*) queries->d_queryInfo, numQueries,
                                                                                    (uint2*) queries->d_regions, (ulonglong2*) regions->d_intervals, (uint2*) regions->d_regionsOffsets,
                                                                                    maxExtraSteps, occThreshold, occShrinkFactor, maxRegionsFactor);
    break;
    case GPU_FMI_TABLE_MULTILEVEL_LINKED:
      gpu_fmi_asearch_table_linked_kernel<<<blocksPerGrid, threadsPerBlock, 0, idStream>>>((gpu_fmi_device_entry_t*) index->fmi.d_fmi[idSupDev], index->fmi.bwtSize,
                                                                                           (uint64_t*) table->d_fmiTableLUT[idSupDev], (uint2*) table->d_offsetsTableLUT[idSupDev], table->maxLevelsTableLUT,
                                                                                           (char*) queries->d_queries, (uint2*) queries->d_queryInfo, numQueries,
                                                                                           (uint2*) queries->d_regions, (ulonglong2*) regions->d_intervals, (uint2*) regions->d_regionsOffsets,
                                                                                           maxExtraSteps, occThreshold, occShrinkFactor, maxRegionsFactor);
    break;
  }
  return(SUCCESS);
}

#endif /* GPU_FMI_SEARCH_CU_ */
