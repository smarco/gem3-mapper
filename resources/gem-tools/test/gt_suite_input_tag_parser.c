/*
 * PROJECT: GEM-Tools library
 * FILE: gt_suite_input_map_parser.c
 * DATE: 03/10/2012
 * DESCRIPTION: // TODO
 */

#include "gt_test.h"

gt_string* tag;
gt_string* expected;
gt_string* expected_casava;
gt_string* expected_extra;
gt_attributes* attributes;
gt_template* template;
gt_output_map_attributes* output_attributes;
gt_output_fasta_attributes* fastq_attributes;

void gt_input_tag_parser_setup(void) {
	tag = gt_string_new(1024);
	expected = gt_string_new(1024);
	expected_casava = gt_string_new(1024);
	expected_extra = gt_string_new(1024);
	attributes = gt_attributes_new();
	template = gt_template_new();
	output_attributes = gt_output_map_attributes_new();
	fastq_attributes = gt_output_fasta_attributes_new();
}

void gt_input_tag_parser_teardown(void) {
	gt_string_delete(tag);
	gt_string_delete(expected);
	gt_string_delete(expected_casava);
	gt_string_delete(expected_extra);
	gt_attributes_delete(attributes);
	gt_template_delete(template);
	gt_output_map_attributes_delete(output_attributes);
	gt_output_fasta_attributes_delete(fastq_attributes);
}


START_TEST(gt_test_basic_tag_parsing)
{
	char* input[1];

	// basic
	input[0] = "mytag";
	gt_string_set_string(expected, "mytag");
	char* tag_begin = *input;
	int pair = 1;

	fail_unless(gt_input_parse_tag((const char** const)input, tag, attributes) == GT_STATUS_OK, "Basic tag not parsed");
	fail_unless(gt_string_cmp(tag, expected) == 0, "Tag not parsed correctly");
	pair = *((int*)gt_attributes_get(attributes,GT_ATTR_ID_TAG_PAIR));
	fail_unless(pair == 0, "Pair information not parsed");

	// pair info /1
	gt_string_clear(tag);
	input[0] = "mytag/1";
	gt_string_set_string(expected, "mytag");
	tag_begin = *input;
	fail_unless(gt_input_parse_tag((const char** const)input, tag, attributes) == GT_STATUS_OK, "Basic tag not parsed");
	fail_unless(gt_string_cmp(tag, expected) == 0, "Tag not parsed correctly");
	pair = *((int*)gt_attributes_get(attributes, GT_ATTR_ID_TAG_PAIR));
	fail_unless(pair == 1, "Pair information not parsed");

	// pair info /2
	gt_string_clear(tag);
	input[0] = "mytag/2";
	gt_string_set_string(expected, "mytag");
	tag_begin = *input;
	fail_unless(gt_input_parse_tag((const char** const)input, tag, attributes) == GT_STATUS_OK, "Basic tag not parsed");
	fail_unless(gt_string_cmp(tag, expected) == 0, "Tag not parsed correctly");
	pair = *((int*)gt_attributes_get(attributes, GT_ATTR_ID_TAG_PAIR));
	fail_unless(pair == 2, "Pair information not parsed");

	gt_string_clear(tag);
	gt_string_clear(expected_casava);
	gt_string_clear(expected_extra);
	input[0] = "mytag/2 ABC123 XYF";
	gt_string_set_string(expected, "mytag");
	gt_string_set_string(expected_extra, "ABC123 XYF");
	tag_begin = *input;
	fail_unless(gt_input_parse_tag((const char** const)input, tag, attributes) == GT_STATUS_OK, "Basic tag not parsed");
	fail_unless(gt_string_cmp(tag, expected) == 0, "Tag not parsed correctly");
	fail_unless(*((int64_t*)gt_attributes_get(attributes, GT_ATTR_ID_TAG_PAIR)) == 2, "Pair information not parsed");
	fail_unless(gt_string_equals(expected_extra, gt_attributes_get(attributes,GT_ATTR_ID_TAG_EXTRA)), "Extra string not extracted");
}
END_TEST

START_TEST(gt_test_casava_tag_parsing)
{
	char* input[1];

	// basic
	input[0] = "mytag 1:Y:18:ATCACG";
	gt_string_set_string(expected, "mytag");
	gt_string_set_string(expected_casava, "1:Y:18:ATCACG");
	char* tag_begin = *input;
	int pair = 0;

	fail_unless(gt_input_parse_tag((const char** const)input, tag, attributes) == GT_STATUS_OK, "Basic tag not parsed");
	//printf("PARSED TAG "PRIgts"\n",PRIgts_content(tag));
	fail_unless(gt_string_cmp(tag, expected) == 0, "Tag not parsed correctly");
	fail_unless(*((int64_t*)gt_attributes_get(attributes, GT_ATTR_ID_TAG_PAIR)) == 1, "Pair information not parsed, should be 1");
	fail_unless(gt_string_cmp(expected_casava,(gt_string*)gt_attributes_get(attributes, GT_ATTR_ID_TAG_CASAVA)) == 0, "Casava string not extracted");

	// pair info /1
	gt_string_clear(tag);
	gt_string_clear(expected_casava);
	gt_string_clear(expected_extra);
	input[0] = "mytag 2:Y:18:ATCACG B T AAA CCC ### ###";
	gt_string_set_string(expected, "mytag");
	gt_string_set_string(expected_casava, "2:Y:18:ATCACG");
	gt_string_set_string(expected_extra, "B T AAA CCC ### ###");
	tag_begin = *input;
	fail_unless(gt_input_parse_tag((const char** const)input, tag, attributes) == GT_STATUS_OK, "Basic tag not parsed");
	fail_unless(gt_string_cmp(tag, expected) == 0, "Tag not parsed correctly");
	fail_unless(*((int64_t*)gt_attributes_get(attributes, GT_ATTR_ID_TAG_PAIR)) == 2, "Pair information not parsed, should be 2");
	fail_unless(gt_string_cmp(expected_casava,(gt_string*)gt_attributes_get(attributes, GT_ATTR_ID_TAG_CASAVA)) == 0, "Casava string not extracted");
	fail_unless(gt_string_cmp(expected_extra,(gt_string*)gt_attributes_get(attributes, GT_ATTR_ID_TAG_EXTRA)) == 0, "Extra string not extracted");
}
END_TEST

START_TEST(gt_test_casava_tag_parsing_extended)
{
	char* input[1];

	// basic
	input[0] = "mytag X:Y:18:ATCACG/1";
	gt_string_set_string(expected, "mytag_X:Y:18:ATCACG");
	gt_string_set_string(expected_casava, "");
	register char* tag_begin = *input;
	int pair = 0;

	fail_unless(gt_input_parse_tag((const char** const)input, tag, attributes) == GT_STATUS_OK, "Basic tag not parsed");
	fail_unless(gt_string_cmp(tag, expected) == 0, "Tag not parsed correctly");
	fail_unless(*gt_shash_get(attributes, GT_ATTR_ID_TAG_PAIR, int64_t) == 1, "Pair information not parsed, should be 1");
	
	gt_string_clear(tag);
	gt_string_clear(expected_casava);
	gt_string_clear(expected_extra);
	pair = 0;

	input[0] = "@SRR384920.1 HWI-ST382_0049:1:1:1217:1879/1";
	gt_string_set_string(expected, "@SRR384920.1_HWI-ST382_0049:1:1:1217:1879");
	gt_string_set_string(expected_casava, "");
	tag_begin = *input;
	pair = 0;

	fail_unless(gt_input_parse_tag((const char** const)input, tag, attributes) == GT_STATUS_OK, "Basic tag not parsed");
	fail_unless(gt_string_cmp(tag, expected) == 0, "Tag not parsed correctly");
	fail_unless(*gt_shash_get(attributes, GT_ATTR_ID_TAG_PAIR, int64_t) == 1, "Pair information not parsed, should be 1");
	
	
	
	gt_string_clear(tag);
	gt_string_clear(expected_casava);
	gt_string_clear(expected_extra);
}
END_TEST


START_TEST(gt_test_tag_parsing_generic_parser_single_paired)
{
	gt_input_file* input = gt_input_file_open("testdata/single_paired.map", false);
	gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input);

	gt_generic_parser_attributes* attr = gt_input_generic_parser_attributes_new(false);

	gt_string_set_string(tag, "myid");
	gt_alignment* alignment;

	// check first template, one alignment first pair
	fail_unless(gt_input_generic_parser_get_template(buffered_input, template, attr) == GT_STATUS_OK, "Failed to read input");
	fail_unless(gt_string_cmp(template->tag, tag) == 0, "Tag is not myid");
	fail_unless(*((int64_t*)gt_attributes_get(template->attributes, GT_ATTR_ID_TAG_PAIR)) == 1, "Pair information not parsed, should be 1");
	fail_unless(gt_template_get_num_blocks(template) == 1, "Found more than 1 block");
	alignment = gt_template_get_block(template, 0);
	fail_unless(gt_string_cmp(alignment->tag, tag) == 0, "Alignment tag is not myid");
	fail_unless(*((int64_t*)gt_attributes_get(alignment->attributes, GT_ATTR_ID_TAG_PAIR)) == 1, "Alignment pair information not 1");


	// check second template, one alignment first pair
	fail_unless(gt_input_generic_parser_get_template(buffered_input, template, attr) == GT_STATUS_OK, "Failed to read input");
	fail_unless(gt_string_cmp(template->tag, tag) == 0, "Tag is not myid");
	fail_unless(*((int64_t*)gt_attributes_get(template->attributes, GT_ATTR_ID_TAG_PAIR)) == 2, "Pair information not parsed, should be 2");
	fail_unless(gt_template_get_num_blocks(template) == 1, "Found more than 1 block");
	alignment = gt_template_get_block(template, 0);
	fail_unless(gt_string_cmp(alignment->tag, tag) == 0, "Alignment tag is not myid");
	fail_unless(*((int64_t*)gt_attributes_get(alignment->attributes, GT_ATTR_ID_TAG_PAIR)) == 2, "Alignment pair information not 2");

	gt_buffered_input_file_close(buffered_input);
	gt_input_file_close(input);
}
END_TEST

START_TEST(gt_test_tag_parsing_generic_parser_single_paired_map_output)
{
	gt_input_file* input = gt_input_file_open("testdata/single_paired.map", false);
	gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input);
	gt_generic_parser_attributes* attr = gt_input_generic_parser_attributes_new(false);
	// check first template, one alignment first pair
	fail_unless(gt_input_generic_parser_get_template(buffered_input, template, attr) == GT_STATUS_OK, "Failed to read input");
	gt_output_map_sprint_template(expected, template, output_attributes);
	// convert to string
	gt_string_set_string(tag, "myid/1\tACGT\t####\t1\tchr1:+:10:4\n");
	fail_unless(gt_string_cmp(tag, expected) == 0, "Not the right output: '%s'\n", gt_string_get_string(expected));
	gt_buffered_input_file_close(buffered_input);
	gt_input_file_close(input);
}
END_TEST

START_TEST(gt_test_tag_parsing_generic_parser_single_paired_map_output_casava_additional)
{
	gt_input_file* input = gt_input_file_open("testdata/single_paired_casava_additional.map", false);
	gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input);
	gt_generic_parser_attributes* attr = gt_input_generic_parser_attributes_new(false);
	// check first template, one alignment first pair
	fail_unless(gt_input_generic_parser_get_template(buffered_input, template, attr) == GT_STATUS_OK, "Failed to read input");
	gt_output_map_sprint_template(expected, template, output_attributes);

	// convert to string
	gt_string_set_string(tag, "myid 1:Y:18:ATCACG B T AAA CCC ### ###\tACGT\t####\t1\tchr1:+:10:4\n");
	fail_unless(gt_string_cmp(tag, expected) == 0, "Not the right output: '%s'\n", gt_string_get_string(expected));
	gt_buffered_input_file_close(buffered_input);
	gt_input_file_close(input);
}
END_TEST

START_TEST(gt_test_tag_parsing_generic_parser_single_paired_map_output_casava_additional_no_casava)
{
	gt_input_file* input = gt_input_file_open("testdata/single_paired_casava_additional.map", false);
	gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input);
	gt_generic_parser_attributes* attr = gt_input_generic_parser_attributes_new(false);
	gt_output_map_attributes_set_print_casava(output_attributes, false);
	// check first template, one alignment first pair
	fail_unless(gt_input_generic_parser_get_template(buffered_input, template, attr) == GT_STATUS_OK, "Failed to read input");
	gt_output_map_sprint_template(expected, template, output_attributes);

	// convert to string
	gt_string_set_string(tag, "myid/1 B T AAA CCC ### ###\tACGT\t####\t1\tchr1:+:10:4\n");
	fail_unless(gt_string_cmp(tag, expected) == 0, "Not the right output: '%s'\n", gt_string_get_string(expected));
	gt_buffered_input_file_close(buffered_input);
	gt_input_file_close(input);
}
END_TEST


START_TEST(gt_test_tag_parsing_generic_parser_single_paired_map_output_casava_additional_no_casava_no_extra)
{
	gt_input_file* input = gt_input_file_open("testdata/single_paired_casava_additional.map", false);
	gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input);
	gt_generic_parser_attributes* attr = gt_input_generic_parser_attributes_new(false);
	gt_output_map_attributes_set_print_casava(output_attributes, false);
	gt_output_map_attributes_set_print_extra(output_attributes, false);
	// check first template, one alignment first pair
	fail_unless(gt_input_generic_parser_get_template(buffered_input, template, attr) == GT_STATUS_OK, "Failed to read input");
	gt_output_map_sprint_template(expected, template, output_attributes);

	// convert to string
	gt_string_set_string(tag, "myid/1\tACGT\t####\t1\tchr1:+:10:4\n");
	fail_unless(gt_string_cmp(tag, expected) == 0, "Not the right output: '%s'\n", gt_string_get_string(expected));
	gt_buffered_input_file_close(buffered_input);
	gt_input_file_close(input);
}
END_TEST

START_TEST(gt_test_tag_parsing_generic_parser_single_paired_map_output_casava_additional_no_casava_no_extra_fastq)
{
	gt_input_file* input = gt_input_file_open("testdata/single_paired_casava_additional.fastq", false);
	gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input);
	gt_generic_parser_attributes* attr = gt_input_generic_parser_attributes_new(false);

	gt_output_fasta_attributes_set_print_casava(fastq_attributes, false);
	gt_output_fasta_attributes_set_print_extra(fastq_attributes, false);
	// check first template, one alignment first pair
	fail_unless(gt_input_generic_parser_get_template(buffered_input, template, attr) == GT_STATUS_OK, "Failed to read input");
	gt_output_fasta_sprint_template(expected, template, fastq_attributes);

	// convert to string
	gt_string_set_string(tag, "@myid/1\nACGT\n+\n####\n");
	fail_unless(gt_string_cmp(tag, expected) == 0, "Not the right output: '%s'\n", gt_string_get_string(expected));
	gt_buffered_input_file_close(buffered_input);
	gt_input_file_close(input);
}
END_TEST

START_TEST(gt_test_tag_parsing_generic_parser_single_paired_map_output_casava_additional_fasta)
{
	gt_input_file* input = gt_input_file_open("testdata/single_paired_casava_additional.fastq", false);
	gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input);
	gt_generic_parser_attributes* attr = gt_input_generic_parser_attributes_new(false);

	// check first template, one alignment first pair
	fail_unless(gt_input_generic_parser_get_template(buffered_input, template, attr) == GT_STATUS_OK, "Failed to read input");
	gt_output_fasta_sprint_template(expected, template, fastq_attributes);

	// convert to string
	gt_string_set_string(tag, "@myid 1:Y:18:ATCACG B T AAA CCC ### ###\nACGT\n+\n####\n");
	fail_unless(gt_string_cmp(tag, expected) == 0, "Not the right output: '%s'\n", gt_string_get_string(expected));
	gt_buffered_input_file_close(buffered_input);
	gt_input_file_close(input);
}
END_TEST

START_TEST(gt_test_tag_parsing_generic_parser_single_paired_casava_additional)
{
	gt_input_file* input = gt_input_file_open("testdata/single_paired_casava_additional.map", false);
	gt_buffered_input_file* buffered_input = gt_buffered_input_file_new(input);

	gt_generic_parser_attributes* attr = gt_input_generic_parser_attributes_new(false);

	gt_string_set_string(tag, "myid");
	gt_string_set_string(expected_casava, "1:Y:18:ATCACG");
	gt_string_set_string(expected_extra, "B T AAA CCC ### ###");

	gt_alignment* alignment;

	// check first template, one alignment first pair
	fail_unless(gt_input_generic_parser_get_template(buffered_input, template, attr) == GT_STATUS_OK, "Failed to read input");
	fail_unless(gt_string_cmp(template->tag, tag) == 0, "Tag is not myid");
	fail_unless(gt_string_cmp(gt_attributes_get(template->attributes, GT_ATTR_ID_TAG_CASAVA), expected_casava) == 0, "Casava string does not match");
	fail_unless(gt_string_cmp(gt_attributes_get(template->attributes, GT_ATTR_ID_TAG_EXTRA), expected_extra) == 0, "Extra String does not match");
	fail_unless(*((int64_t*)gt_attributes_get(template->attributes, GT_ATTR_ID_TAG_PAIR)) == 1, "Pair information not parsed, should be 1");
	fail_unless(gt_template_get_num_blocks(template) == 1, "Found more than 1 block");
	alignment = gt_template_get_block(template, 0);
	fail_unless(gt_string_cmp(alignment->tag, tag) == 0, "Alignment tag is not myid");
	fail_unless(*((int64_t*)gt_attributes_get(alignment->attributes, GT_ATTR_ID_TAG_PAIR)) == 1, "Alignment pair information not 1");
	fail_unless(gt_string_cmp(gt_attributes_get(alignment->attributes, GT_ATTR_ID_TAG_CASAVA), expected_casava) == 0, "Casava string does not match");
	fail_unless(gt_string_cmp(gt_attributes_get(alignment->attributes, GT_ATTR_ID_TAG_EXTRA), expected_extra) == 0, "Extra String does not match");


	gt_string_set_string(expected_casava, "2:Y:18:ATCACG");
	gt_string_set_string(expected_extra, "B T AAA TTT ### ###");
	// check second template, one alignment first pair
	fail_unless(gt_input_generic_parser_get_template(buffered_input, template, attr) == GT_STATUS_OK, "Failed to read input");
	fail_unless(gt_string_cmp(template->tag, tag) == 0, "Tag is not myid");
	fail_unless(gt_string_cmp(gt_attributes_get(template->attributes, GT_ATTR_ID_TAG_CASAVA), expected_casava) == 0, "Casava string does not match");
	fail_unless(gt_string_cmp(gt_attributes_get(template->attributes, GT_ATTR_ID_TAG_EXTRA), expected_extra) == 0, "Extra String does not match");
	fail_unless(*((int64_t*)gt_attributes_get(template->attributes, GT_ATTR_ID_TAG_PAIR)) == 2, "Pair information not parsed, should be 2");
	fail_unless(gt_template_get_num_blocks(template) == 1, "Found more than 1 block");
	alignment = gt_template_get_block(template, 0);
	fail_unless(gt_string_cmp(alignment->tag, tag) == 0, "Alignment tag is not myid");
	fail_unless(*((int64_t*)gt_attributes_get(alignment->attributes, GT_ATTR_ID_TAG_PAIR)) == 2, "Alignment pair information not 2");
	fail_unless(gt_string_cmp(gt_attributes_get(alignment->attributes, GT_ATTR_ID_TAG_CASAVA), expected_casava) == 0, "Casava string does not match");
	fail_unless(gt_string_cmp(gt_attributes_get(alignment->attributes, GT_ATTR_ID_TAG_EXTRA), expected_extra) == 0, "Extra String does not match");

	gt_buffered_input_file_close(buffered_input);
	gt_input_file_close(input);
}
END_TEST



Suite *gt_input_tag_parser_suite(void) {
  Suite *s = suite_create("gt_input_parser");

  /* String parsers test case */
  TCase *tc_tag_string_parser = tcase_create("TAG parser. String parsers");
  tcase_add_checked_fixture(tc_tag_string_parser,gt_input_tag_parser_setup,gt_input_tag_parser_teardown);
  tcase_add_test(tc_tag_string_parser,gt_test_basic_tag_parsing);
  tcase_add_test(tc_tag_string_parser,gt_test_casava_tag_parsing);
  tcase_add_test(tc_tag_string_parser,gt_test_casava_tag_parsing_extended);
  tcase_add_test(tc_tag_string_parser,gt_test_tag_parsing_generic_parser_single_paired);
  tcase_add_test(tc_tag_string_parser,gt_test_tag_parsing_generic_parser_single_paired_casava_additional);
  tcase_add_test(tc_tag_string_parser,gt_test_tag_parsing_generic_parser_single_paired_map_output);
  tcase_add_test(tc_tag_string_parser,gt_test_tag_parsing_generic_parser_single_paired_map_output_casava_additional);
  tcase_add_test(tc_tag_string_parser,gt_test_tag_parsing_generic_parser_single_paired_map_output_casava_additional_no_casava);
  tcase_add_test(tc_tag_string_parser,gt_test_tag_parsing_generic_parser_single_paired_map_output_casava_additional_no_casava_no_extra);
  tcase_add_test(tc_tag_string_parser,gt_test_tag_parsing_generic_parser_single_paired_map_output_casava_additional_no_casava_no_extra_fastq);
  tcase_add_test(tc_tag_string_parser,gt_test_tag_parsing_generic_parser_single_paired_map_output_casava_additional_fasta);

  suite_add_tcase(s,tc_tag_string_parser);

  return s;
}
