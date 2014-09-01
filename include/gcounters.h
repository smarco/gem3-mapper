/*
 * PROJECT: GEMMapper
 * FILE: gcounters.h
 * DATE: 06/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: GEM Profile Counters
 *   Basically, this counters/timers focus on the measurement of:
 *     - Time spent in certain functions/blocks-code (and number of calls/executions)
 *     - Functional profile of execution flow branches (counting)
 */

#ifndef GP_COUNTERS_H_
#define GP_COUNTERS_H_

/*
 * Approximate search
 */
#define GP_AS_MAIN              0
#define GP_AS_READ_RECOVERY     1
#define GP_AS_EXACT_SEARCH      2

#define GP_REGION_PROFILE_SOFT 30

#define GP_AS_FILTER_REGIONS                       60
#define GP_AS_FILTER_REGIONS_NUM_ELEGIBLE_REGIONS  61
#define GP_AS_FILTER_REGIONS_SEARCH_D2             62
#define GP_AS_FILTER_REGIONS_SEARCH_D2_HIT         63
#define GP_AS_FILTER_REGIONS_SEARCH_D1             64
#define GP_AS_FILTER_REGIONS_SEARCH_D1_HIT         65
#define GP_AS_FILTER_REGIONS_SEARCH_D0             66
#define GP_AS_FILTER_REGIONS_PROCESSED             67
#define GP_AS_FILTER_REGIONS_SKIPPED               68


#define GP_
#define GP_
#define GP_
#define GP_
#define GP_

/*
 * Region Profile
 */
#define GP_REGION_PROFILE_ADAPTIVE      80
#define GP_REGION_PROFILE_QUIT_PROFILE  81

//
//
//#define GP_READ_RECOVERY 666
//
//#define GP_ATH_QUERY 666
//
//#define GC_FILTER_REGIONS 666
//#define GC_CAND_FILTER_REGIONS 666
//
//#define SC_FILTER_REGIONS 666
//#define GSC_CAND_FILTER_REGIONS 666
//#define GSC_ATH_D2 666
//#define GSC_ATH_D2_HIT 666
//#define GSC_ATH_D1 666
//#define GSC_ATH_D1_HIT 666
//#define GSC_ATH_D0 666
//#define GSC_ATH_FILTER_CAND 666
//#define SC_FILTER_REGIONS 666
//#define SC_ATH_FILTER 666
//#define GSC_SAVED_FILTER_REGIONS 666
//#define GSC_FILTER_REGIONS 666
//
//#define GP_NUM_SOFT_REGIONS 666
//#define GSC_ATH_HIT 666
//
//#define SC_EXACT_SEARCH 666
//#define SC_SMALL_READS 666
//#define SC_ATH 666
//
//#define SC_REGION_PROFILE 666
//#define SC_PROBING_DELTA 666


#endif /* GP_COUNTERS_H_ */
