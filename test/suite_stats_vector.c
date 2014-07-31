/*
 * PROJECT: GEMMapper
 * FILE: suite_stats_vector.c
 * DATE: 03/10/2012
 * DESCRIPTION: // TODO
 */

#include "gem_utest.h"

//
//START_TEST(gt_test_ihash_basic_insertions)
//{
//  // Insert integer
//  register uint64_t* integer = malloc(sizeof(uint64_t));
//  *integer = 9000;
//  gt_ihash_insert(ihash,200,integer,uint64_t);
//  fail_unless(*gt_ihash_get(ihash,200,uint64_t)==9000,"Failed inserting integer into ihash");
//  // Insert string
//  register char* string = malloc(200);
//  strcpy(string,"HELLO WORLD");
//  gt_ihash_insert(ihash,201,string,char*);
//  fail_unless(strcmp(gt_ihash_get(ihash,201,char),"HELLO WORLD")==0,"Failed inserting string into ihash");
//  // Insert enum
//  typedef enum { VALUE_1, VALUE_2, VALUE_3 } custom_enum;
//  register custom_enum* enumeration = malloc(sizeof(custom_enum));
//  *enumeration = VALUE_2;
//  gt_ihash_insert(ihash,202,enumeration,custom_enum);
//  fail_unless(*gt_ihash_get(ihash,202,custom_enum)==VALUE_2,"Failed inserting enum into ihash");
//  // Insert custom structure
//  typedef struct {
//    uint64_t field1;
//    uint64_t field2;
//    uint64_t field3;
//  } custom_type;
//  register custom_type* type = malloc(sizeof(custom_type));
//  type->field2 = 2;
//  gt_ihash_insert(ihash,203,type,custom_type);
//  fail_unless(gt_ihash_get(ihash,203,custom_type)->field2==2,"Failed inserting custom type into ihash");
//  // Check all again
//  fail_unless(*gt_ihash_get(ihash,200,uint64_t)==9000,"Failed inserting integer into ihash");
//  fail_unless(strcmp(gt_ihash_get(ihash,201,char),"HELLO WORLD")==0,"Failed inserting string into ihash");
//  fail_unless(*gt_ihash_get(ihash,202,custom_enum)==VALUE_2,"Failed inserting enum into ihash");
//  fail_unless(gt_ihash_get(ihash,203,custom_type)->field2==2,"Failed inserting custom type into ihash");
//}
//END_TEST
//

START_TEST(stats_vector_raw_inc)
{
  stats_vector_t* const stats = stats_vector_raw_new(5,10);
  stats_vector_inc(stats,0);
  stats_vector_inc(stats,2);
  stats_vector_inc(stats,4);
  stats_vector_inc(stats,6);
  stats_vector_inc(stats,8);
  stats_vector_print_raw(stderr,stats);
}
END_TEST


Suite *suite_stats_vector(void) {
  Suite *s = suite_create("stats_vector");

  /* Core test case */
  TCase *tc_core = tcase_create("stats_vector");
  //tcase_add_checked_fixture(tc_core,stats_vector_setup,stats_vector_teardown);
  tcase_add_test(tc_core,stats_vector_raw_inc);
  suite_add_tcase(s,tc_core);

  return s;
}
