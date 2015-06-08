/*
 * PROJECT: GEM-Tools library
 * FILE: gt_utest_core.c
 * DATE: 03/10/2012
 * DESCRIPTION: // TODO
 */

#include "gt_test.h"

// Include Suites
#include "gt_suite_alignment.c"
#include "gt_suite_template_utils.c"
//#include "gt_suite_template.c"

int main(void) {
  SRunner *sr = srunner_create(gt_alignment_suite());
  srunner_add_suite (sr, gt_template_utils_suite());
  
  // add logging to xml
  srunner_set_xml(sr, "reports/check-test-core.xml");

  // Run the suites
  srunner_run_all(sr,CK_NORMAL);
  uint32_t num_test_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  // Return unit tests
  return (num_test_failed==0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
