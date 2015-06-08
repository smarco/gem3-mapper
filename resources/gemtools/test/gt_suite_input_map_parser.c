/*
 * PROJECT: GEM-Tools library
 * FILE: gt_suite_input_map_parser.c
 * DATE: 03/10/2012
 * DESCRIPTION: // TODO
 */

#include "gt_test.h"

gt_map** map;
gt_vector* map_list;
gt_alignment* alignment;
gt_template* template;
gt_output_map_attributes* output_attributes;

void gt_input_map_parser_setup(void) {
  map = gt_malloc_(2, sizeof(gt_map), false, false);
  map[0] = gt_map_new();
  map[1] = gt_map_new();
  gt_map_set_base_length(map[0],300);
  gt_map_set_base_length(map[1],300);
  map_list = gt_vector_new(10,sizeof(gt_map*));
  alignment = gt_alignment_new();
  template = gt_template_new();
  output_attributes = gt_output_map_attributes_new();
}

void gt_input_map_parser_teardown(void) {
  gt_map_delete(map[0]);
  gt_map_delete(map[1]);
  free(map);
  gt_output_map_attributes_delete(output_attributes);
}

START_TEST(gt_test_imp_string_map) // TODO
{
  /*
   * Parse old map
   */
  fail_unless(gt_input_map_parse_map("chr7:F127708134G27<+5>30A34<-5>40T88",map, NULL)==0,"Failed parsing old map");
  fail_unless(gt_strcmp(gt_map_get_seq_name(map[0]),"chr7")==0,"Failed parsing old map. Seq Name");
  fail_unless(gt_map_get_strand(map[0])==FORWARD,"Failed parsing old map. Strand");

  /*
   * Parse new map
   */
  fail_unless(gt_input_map_parse_map("chr12:+:9570521:6C2>1+1>1-3T>1-2T12>4-8T23T8",map, NULL)==0,"Failed parsing old split-map");

  /*
   * Parsing single maps
   */
  // Parse old map
  fail_unless(gt_input_map_parse_map("chr7:F127708134G27<+5>30A34<-5>40T88",map, NULL)==0);
  // Parse new SE-map
  fail_unless(gt_input_map_parse_map("chr12:+:9570521:6C2>1+1>1-3T>1-2T12>4-8T23T8",map, NULL)==0);
  // Parse new PE-map
  fail_unless(gt_input_map_parse_map("chr15:-:102516634:66G9::chr15:+:102516358:66>1+10:::7936",map, NULL)==0);

  /*
   * Parsing list of maps
   */
  // Parse multiple old split-map
  // FIXME: GT-77 - check the test
//  fail_unless(gt_input_map_parse_map_list("[26]=chr7:R1203797~chr7:R1203108",map_list, NULL)==0);
//  fail_unless(gt_input_map_parse_map_list("[31;35]=chr16:R[2503415;2503411]~chr16:R2503271",map_list, NULL)==0);
//  fail_unless(gt_input_map_parse_map_list("[30;34]=chr10:F74776624~chr10:F[74790025;74790029]",map_list, NULL)==0);
//  fail_unless(gt_input_map_parse_map_list("[23-50]=chr1:F[188862944-188868041]~chr19:F53208292",map_list, NULL)==0);
//  fail_unless(gt_input_map_parse_map_list("[70-71]=chr1:F188862944~chr19:F[53208292-53208293]",map_list, NULL)==0);
//  fail_unless(gt_input_map_parse_map_list("[26]=chr7:R1203797~chr7:R1203108",map_list, NULL)==0);
//  fail_unless(gt_input_map_parse_map_list(
//      "chrM:F6598<+1>16@0/0,[23]=chr6:R31322884~chr6:R31237276,"
//      "[70-71]=chr1:F188862944~chr19:F[53208292-53208293]",map_list, NULL)==0);
//  // Parse multiple old SE-map
//  fail_unless(gt_input_map_parse_map_list("chr1:F8926499@0/0,chr12:R7027116G39A42@77/2",map_list, NULL)==0);
//  // Parse multiple new SE-map
//  fail_unless(gt_input_map_parse_map_list(
//      "chrX:-:155255234:1T36A37,chrY:-:59358240:1T36A37:200,"
//      "chr15:-:102516664:1>1-28>5+8A37,chr16:+:64108:3>1-30>1+1>4+3A37,"
//      "chr9:+:14540:3>1-34A33A2>1-",map_list, NULL)==0);
//  // Parse multiple new PE-map
//  fail_unless(gt_input_map_parse_map_list(
//      "chr15:-:102516742:(3)3GCA67::chr15:+:102516611:76:::7936,"
//      "chr16:+:64114:1>1-26>1+1>4+47::chr16:-:64196:68A6T:::12224,"
//      "chr1:+:16731:(5)35>92*16(20),chrY:-:59355959:(5)6G4G24>1-3A1CA1AAA1>1-1(20),"
//      "chrX:-:155252953:(5)6G4G24>1-3A1CA1AAA1>1-1(20)",map_list, NULL)==0);

  /*
   * Parsing counters
   */

  /*
   * Parsing Alignments
   */
  fail_unless(gt_input_map_parse_alignment(
      "C0LMTACXX120523:8:2209:19417:37092/1\t"
      "ACTCCAGTCACTCCAGCAAATCTCGGATGCCGTCTTCTGCTTGAACGAAACCAGAACTGTGTGGAGAACAGCTTAA\t"
      "7=7AAAA7<722<C+AA;=C?<3A3+2<):@)?###########################################\t"
      "0+0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:1:0:1:0:1:0:1:1\t"
      "chr10:+:28549954:(5)17A1AA1GC1A1ATGGG1C1CCA1TTC2T1A1TT1(20):65472,"
      "chr13:+:43816250:(5)10A7ACTC2AGA1A2TCCAG3C1CTCTT1TTGC(20):12224,"
      "chr20:+:21833059:(5)5T12GCTC1AAAAG3TGCTGA1C1ACCT2AGTGT(20):1984,"
      "chr4:+:94518781:(5)12G6CTTAT1TCAAAA1CCAGAAGT1CC2GACACT(20):1984,"
      "chr14:-:74094161:(5)17AAACCCAA1AAGGAGGT1A1TGGAACCC1A1CGT(20):1984",alignment)==0);
  fail_unless(gt_input_map_parse_alignment(
      "C0LMTACXX120523:8:2209:19417:37092/1\t"
      "ACTCCAGTCACTCCAGCAAATCTCGGATGCCGTCTTCTGCTTGAACGAAACCAGAACTGTGTGGAGAACAGCTTAA\t"
      "7=7AAAA7<722<C+AA;=C?<3A3+2<):@)?###########################################\t"
      "0+0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:1:0:1:0:1:0:1:1\t"
      "-",alignment)==0);

  /*
   * Parsing Templates
   */
  fail_unless(gt_input_map_parse_template(
      "C0LMTACXX120523:8:2209:19417:37092/1\t"
      "ACTCCAGTCACTCCAGCAAATCTCGGATGCCGTCTTCTGCTTGAACGAAACCAGAACTGTGTGGAGAACAGCTTAA\t"
      "7=7AAAA7<722<C+AA;=C?<3A3+2<):@)?###########################################\t"
      "0+0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:1:0:1:0:1:0:1:1\t"
      "chr10:+:28549954:(5)17A1AA1GC1A1ATGGG1C1CCA1TTC2T1A1TT1(20),"
      "chr13:+:43816250:(5)10A7ACTC2AGA1A2TCCAG3C1CTCTT1TTGC(20),"
      "chr20:+:21833059:(5)5T12GCTC1AAAAG3TGCTGA1C1ACCT2AGTGT(20),"
      "chr4:+:94518781:(5)12G6CTTAT1TCAAAA1CCAGAAGT1CC2GACACT(20),"
      "chr14:-:74094161:(5)17AAACCCAA1AAGGAGGT1A1TGGAACCC1A1CGT(20)",template)==0);

  /*
   * Check mcs delegation for single end
   */
  fail_unless(gt_template_get_mcs(template) == 1);

    /*
   * Parsing Templates
   */
  fail_unless(gt_input_map_parse_template(
      "A/1\t"
      "AAA\t"
      "###\t"
      "0\t"
      "-",template)==0);

}
END_TEST

Suite *gt_input_map_parser_suite(void) {
  Suite *s = suite_create("gt_input_map_parser");

  /* String parsers test case */
  TCase *tc_map_string_parser = tcase_create("MAP parser. String parsers");
  tcase_add_checked_fixture(tc_map_string_parser,gt_input_map_parser_setup,gt_input_map_parser_teardown);
  tcase_add_test(tc_map_string_parser,gt_test_imp_string_map);
  suite_add_tcase(s,tc_map_string_parser);

  return s;
}

//  gt_map* map = gt_map_new();
//
//  map->base_length = 25;
//  gt_map_realign_levenshtein(
//      map,
//      "AAAAAAAAAAAAAAAAAAAAAAAAA",map->base_length,
//      "AAAA",4,false);
//
//  map->base_length = 25;
//  gt_map_realign_levenshtein(
//      map,
//      "AAAAAAAAAAAAGAAAAAAAAAAAA",map->base_length,
//      "AAAATAAAAAAAAAAAAAAAACAAA",25,false);
//
//  map->base_length = 4;
//  gt_map_realign_levenshtein(
//      map,
//      "AAAA",map->base_length,
//      "AAAAAAAAAAAAAAAAAAAAAAAAA",25,false);
//
//  map->base_length = 4;
//  gt_map_realign_levenshtein(
//      map,
//      "CCCC",map->base_length,
//      "AAAAAAAAAACCCCAAAAAAAAAAA",25,false);
//
//  map->base_length = 19;
//  gt_map_realign_levenshtein(
//      map,
//      "AAAAAAAAAAAAAAAAAAA",map->base_length,
//      "AAAAACCAAAACCCAAAAAAACAAA",25,false);
//
//  map->base_length = 25;
//  gt_map_realign_levenshtein(
//      map,
//      "AAAAAAAAAAAAAAAAAAAAAAAAA",map->base_length,
//      "AAAATAAAAAAAAAAAAAAAAAAAAC",26,false);
//
//  map->base_length = 4;
//  gt_map_realign_levenshtein(
//      map,
//      "CCCC",map->base_length,
//      "AAAAAAAAAACCCCAAAAAAAAAAA",25,true);

//  map->base_length = 363;
//  gt_map_realign_levenshtein(
//      map,
//      "TAATTGCTATATCCCTCAAACATCCTTTACCCTGAAATCCCTTCTAATCCATCCTCTGCCACTGCTTCCAGATTATTCTCTCTGAAATCAAGTCTAATCATGTCACT"
//      "TTTTAGCTTAAAATACTTCAATGGCACTCCATAGTTAACCAGACAGGAAGAAAGTAAAGCATACGGTCAAGAGTCCTGGCTCTAGAGTGAGACTGCCTGGGTTCAAA"
//      "ATCCTAGTATGACAGTTAATAAATCTTAATACCTGTGTGAACTTGGGAGGATGACTTCACTTCTCCTTTGCCTTCAGTTGCTTATCTAAATGAGTTAATGTAATGTA"
//      "AAGCACATGCCACACTGAAGTACTTTAATCAATATTAGCTGTTATTGTAAGTTCAAGTTTTGTAGTTAAATT",map->base_length,
//      "TAATTGCTATATCCCTCAAACATCCTTTACCCTGAATCCCTTCTAATCCATCCTCTGCACTGCTTCCAGATTATTCTC"
//      "TCTGAAAATCAAGTCTAATCATGTCACTTTTTAGCTTAAAATACTTCAATGGCACTCCATAGTTAACCAGACAGGAAG"
//      "AAAGTAAAGCATACGGTCAAGAGTCCTGGCTCTAGAGTGAGACTGCCTGGGTTCAAATCCTAGTATGACAGTTAATAA"
//      "ATCTTAATACCTGTGTGAACTTGGGAGGATGACTTCACTTCTCCTTTGCCTCAGTTGCTTTATCTAAATGAGTTAATG"
//      "TATGTAAAGCACATGCCACACTGAAGTACTTTAATCAATATTAGCTGTTATTGTAAGTTCAAGTTTTGTAGTTTAAAT"
//      "TCCTTAAGAAAACTCCCAAAAAACAGACGTCATATCATGATCTTGCCCCTTTCTACTACTTATGAACCTCCCCAAAGCTAT",471,true);

//  map->base_length = 363;
//  gt_map_realign_levenshtein(
//      map,
//      "TTAGATTGGGTTGGCTGGTATGCATGAAAATGACAGACCACTATAATTTTCCTTACAAAGAAAAATCTATGCAGTTGGA"
//      "TGGTTTCTTTTAAAAATGAACTAATTTTGATTATTTGCTAACTTTCCCAGCTTTATTGTCCAGAAACAATAGTCCTTGG"
//      "AATAAAGAAAATGTCAAAGAGTAAAACAAGCCAGCGCATTTAAATTGTAAGAATTATTTTTAAAAAATAAGATTGGACT"
//      "GGACTGCAATTTTATAACTGAACCACATTTAATTCTATCTTGCATGGGGTCACTGCACAACATGATTTGAGTTCTCCTT"
//      "AGAGCTTTCCCATCTCTTCCTAGGAGGCTGAAGAATTTATGGAGACA",map->base_length,
//      "CCACTTTGATTGGGTTGGCTGGTAATGCATTGAAATGACAGACACTATAATTTTCCTTACAAAGAAAAATCTATGCAGTTG"
//      "GATGGTTTCTTTAAAAATGAACTAATTTTGATTATTTGCTAACTTTCCCAGTTCTTATGTCAGAAACAATAGTCCTTGAAT"
//      "AAAGAAAATGTCAAAGAGTAAAAACAAGCCAGCGCATTTAAATTGTAGAATTATTTTTAAAAAATAAGATTGGACTGGACT"
//      "GCAATTTTATAACTGAACACATTTTAATTCTATCTTGCATGGGGTCACTGCACACATGATTTGAGTTCTCCTTAGAGCTTT"
//      "CCCATCTCTTCCTAGGAGGCTGAAGAATTTATGGAGACAAAATAGGCAAGGACATTTCTTAGAGAATATGCAATGCAGTTC"
//      "ATCCAGAATGACATCTTAGAGGTTATTTTG",435,true);
