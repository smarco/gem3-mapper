/*
 * PROJECT: GEM-Tools library
 * FILE: gt_suite_alignment.c
 * DATE: 03/10/2012
 * DESCRIPTION: // TODO
 */

#include "gt_test.h"

gt_ihash* ihash;

void gt_ihash_setup(void) {
  ihash = gt_ihash_new();
}

void gt_ihash_teardown(void) {
  gt_ihash_delete(ihash,true);
}

START_TEST(gt_test_ihash_basic_insertions)
{
  // Insert integer
  uint64_t* integer = gt_alloc(uint64_t);
  *integer = 9000;
  gt_ihash_insert(ihash,200,integer,uint64_t);
  fail_unless(*gt_ihash_get(ihash,200,uint64_t)==9000,"Failed inserting integer into ihash");
  // Insert string
  char* string = gt_malloc(200);
  strcpy(string,"HELLO WORLD");
  gt_ihash_insert(ihash,201,string,char*);
  fail_unless(strcmp(gt_ihash_get(ihash,201,char),"HELLO WORLD")==0,"Failed inserting string into ihash");
  // Insert enum
  typedef enum { VALUE_1, VALUE_2, VALUE_3 } custom_enum;
  custom_enum* enumeration = gt_alloc(custom_enum);
  *enumeration = VALUE_2;
  gt_ihash_insert(ihash,202,enumeration,custom_enum);
  fail_unless(*gt_ihash_get(ihash,202,custom_enum)==VALUE_2,"Failed inserting enum into ihash");
  // Insert custom structure
  typedef struct {
    uint64_t field1;
    uint64_t field2;
    uint64_t field3;
  } custom_type;
  custom_type* type = gt_alloc(custom_type);
  type->field2 = 2;
  gt_ihash_insert(ihash,203,type,custom_type);
  fail_unless(gt_ihash_get(ihash,203,custom_type)->field2==2,"Failed inserting custom type into ihash");
  // Check all again
  fail_unless(*gt_ihash_get(ihash,200,uint64_t)==9000,"Failed inserting integer into ihash");
  fail_unless(strcmp(gt_ihash_get(ihash,201,char),"HELLO WORLD")==0,"Failed inserting string into ihash");
  fail_unless(*gt_ihash_get(ihash,202,custom_enum)==VALUE_2,"Failed inserting enum into ihash");
  fail_unless(gt_ihash_get(ihash,203,custom_type)->field2==2,"Failed inserting custom type into ihash");
}
END_TEST

START_TEST(gt_test_ihash_neg_key)
{
  // Insert integer
  uint64_t* integer = gt_alloc(uint64_t);
  *integer = 9000;
  gt_ihash_insert(ihash,-100,integer,uint64_t);
  fail_unless(*gt_ihash_get(ihash,-100,uint64_t)==9000,"Failed inserting integer into ihash (negative key)");
}
END_TEST

Suite *gt_ihash_suite(void) {
  Suite *s = suite_create("gt_ihash");

  /* Core test case */
  TCase *tc_core = tcase_create("ihash accessors");
  tcase_add_checked_fixture(tc_core,gt_ihash_setup,gt_ihash_teardown);
  tcase_add_test(tc_core,gt_test_ihash_basic_insertions);
  tcase_add_test(tc_core,gt_test_ihash_neg_key);
  suite_add_tcase(s,tc_core);

  return s;
}
