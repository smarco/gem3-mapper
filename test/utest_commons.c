/*
 * PROJECT: GEMMapper
 * FILE: utest_commons.c
 * DATE: 03/10/2012
 * DESCRIPTION: // TODO
 */

#include "gem_utest.h"

// Include Suites
#include "suite_stats_vector.c"

int main(void) {
  SRunner *sr = srunner_create(suite_stats_vector());
  //srunner_add_suite(sr,suite_stats_vector());

  // Run the suites
  srunner_run_all(sr,CK_NORMAL);
  uint32_t num_test_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  // Return unit tests
  return (num_test_failed==0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
