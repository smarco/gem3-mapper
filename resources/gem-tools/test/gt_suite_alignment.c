/*
 * PROJECT: GEM-Tools library
 * FILE: gt_suite_alignment.c
 * DATE: 03/10/2012
 * DESCRIPTION: // TODO
 */

#include "gt_test.h"

void gt_alignment_setup(void) {
  // TODO
}

void gt_alignment_teardown(void) {
  // TODO
}

START_TEST(gt_test_alignment_accessors)
{
  // fail_unless(money_amount (five_dollars) == 5,"Amount not set correctly on creation");
  // TODO
  fail_unless(true,"Failed cause zeros");
}
END_TEST

Suite *gt_alignment_suite(void) {
  Suite *s = suite_create("gt_alignment");

  /* Core test case */
  TCase *tc_core = tcase_create("Core");
  tcase_add_checked_fixture(tc_core,gt_alignment_setup,gt_alignment_teardown);
  tcase_add_test(tc_core,gt_test_alignment_accessors);
  // tcase_add_test(tc_core,...);
  suite_add_tcase(s,tc_core);

  return s;
}
