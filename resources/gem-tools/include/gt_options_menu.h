/*
 * PROJECT: GEM-Tools library
 * FILE: gt_options_menu.h
 * DATE: 01/06/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#ifndef GT_OPTIONS_MENU_H_
#define GT_OPTIONS_MENU_H_

// Essentials
#include "gt_essentials.h"

/*
 * Options (Tools Menu)
 */
typedef enum { GT_OPT_NO_ARGUMENT=no_argument, GT_OPT_REQUIRED=required_argument, GT_OPT_OPTIONAL=optional_argument } gt_option_t;
typedef enum { GT_OPT_NONE, GT_OPT_INT, GT_OPT_FLOAT, GT_OPT_CHAR, GT_OPT_STRING, GT_OPT_BOOL } gt_option_argument_t;
typedef struct {
  int option_id;       // Integer ID or short character option
  char* long_option;   // Long string option
  gt_option_t option_type;            // Option type
  gt_option_argument_t argument_type; // Type of the argument
  uint64_t group_id;   // Label of the group it belongs to (zero if none)
  bool active;         // Enable/Disable option
  char* command_info;  // Extra command line syntax info
  char* description;   // Brief description
} gt_option;

extern gt_option gt_filter_options[];
extern char* gt_filter_groups[];

extern gt_option gt_stats_options[];
extern char* gt_stats_groups[];

extern gt_option gt_mapset_options[];
extern char* gt_mapset_groups[];

extern gt_option gt_map2sam_options[];
extern char* gt_map2sam_groups[];

GT_INLINE uint64_t gt_options_get_num_options(const gt_option* const options);
GT_INLINE struct option* gt_options_adaptor_getopt(const gt_option* const options);
GT_INLINE gt_string* gt_options_adaptor_getopt_short(const gt_option* const options);
GT_INLINE void gt_options_fprint_menu(
    FILE* const stream,const gt_option* const options,char* gt_filter_groups[],
    const bool print_description,const bool print_inactive);
GT_INLINE void gt_options_fprint_json_menu(
    FILE* const stream,const gt_option* const options,char* groups[],
    const bool print_description,const bool print_inactive);

#define GT_OPTIONS_ITERATE_BEGIN(tool_options_prefix,option_id) \
  struct option* __getopt = gt_options_adaptor_getopt(gt_##tool_options_prefix##_options); \
  gt_string* const __short_getopt = gt_options_adaptor_getopt_short(gt_##tool_options_prefix##_options); \
  int option_id, option_index; \
  while (true) { \
    /* Get option &  Select case */ \
    if ((option_id=getopt_long(argc,argv,gt_string_get_string(__short_getopt),__getopt,&option_index))==-1) break; \
    switch (option_id)

#define GT_OPTIONS_ITERATE_END \
    gt_string_delete(__short_getopt)

#endif /* GT_OPTIONS_MENU_H_ */
