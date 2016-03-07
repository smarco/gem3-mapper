/*
 * PROJECT: GEMMapper
 * FILE: gem-constructor.c
 * DATE:5/12/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Loads/test things
 */

#include "utils/essentials.h"
#include "utils/packed_integer_array.h"
#include "utils/priority_queue.h"
#include "utils/segmented_vector.h"
#include "utils/sparse_array_locator.h"
#include "utils/sparse_bitmap.h"
#include "utils/vector.h"
#include "utils/options_menu.h"
#include "stats/stats_vector.h"
#include "align/align_ond.h"
#include "archive/locator_builder.h"
#include "archive/archive.h"
#include "archive/archive_text.h"
#include "data_structures/cdna_text.h"
#include "filtering/region_profile.h"
#include "io/input_parser.h"
#include "matches/match_align_dto.h"
#include "matches/match_alignment.h"
#include "neighborhood_search/nsearch_hamming.h"
#include "neighborhood_search/nsearch_levenshtein.h"

/*
 * Generic parameters
 */
typedef struct {
  char *name_input_file;
  char *name_output_file;
  char *option;
  uint64_t number;
  uint64_t param1;
  uint64_t param2;
  uint64_t param3;
  uint64_t param4;
  uint64_t param5;
  uint64_t param6;
  uint64_t param7;
  uint64_t num_threads;
  uint64_t max_memory;
  char* tmp_folder;
} gem_map_filter_args;

gem_map_filter_args parameters = {
    .name_input_file=NULL,
    .name_output_file=NULL,
    .option="",
    .number=0,
    .num_threads=1,
    .max_memory=0,
    .tmp_folder=NULL
};

/*
 * Functions
 */
void constructor_write_fm() {
  /*
   * Write
   */
  fm_t* file = fm_open_file("test.fm",FM_WRITE);
  fm_write_uint64(file,0);
  fm_write_uint64(file,1);
  fm_write_uint64(file,UINT64_MAX);
  fm_close(file);

//  /*
//   * Read
//   */
//  file = fm_open_file("test.sbm",FM_READ);
}
void constructor_stats_vector() {
  stats_vector_t* const stats_raw = stats_vector_raw_new(5,10);
  stats_vector_inc(stats_raw,0);
  stats_vector_inc(stats_raw,1);
  stats_vector_inc(stats_raw,2);
  stats_vector_inc(stats_raw,3);
  stats_vector_inc(stats_raw,4);

  stats_vector_inc(stats_raw,5);
  stats_vector_inc(stats_raw,6);
  stats_vector_inc(stats_raw,7);
  stats_vector_inc(stats_raw,8);
  stats_vector_inc(stats_raw,9);

  stats_vector_inc(stats_raw,15);

  stats_vector_inc(stats_raw,100);
  stats_vector_display(stderr,stats_raw,false,false,NULL);
  fprintf(stderr,"\n");

  uint64_t values[] = {0,100,1000,10000};
  stats_vector_t* const stats_range = stats_vector_customed_range_new(values,3,10000);
  stats_vector_inc(stats_range,0);
  stats_vector_inc(stats_range,1);

  stats_vector_inc(stats_range,100);
  stats_vector_inc(stats_range,101);
  stats_vector_inc(stats_range,500);

  stats_vector_inc(stats_range,1050);
  stats_vector_inc(stats_range,10000);
  stats_vector_inc(stats_range,10001);

  stats_vector_inc(stats_range,20001);
  stats_vector_display(stderr,stats_range,false,false,NULL);
  fprintf(stderr,"\n");

  stats_vector_t* const stats_step = stats_vector_step_range_new(1000,10,1000);
  stats_vector_inc(stats_step,0);
  stats_vector_inc(stats_step,9);

  stats_vector_inc(stats_step,10);
  stats_vector_inc(stats_step,11);

  stats_vector_inc(stats_step,500);

  stats_vector_inc(stats_step,1000);
  stats_vector_inc(stats_step,1001);

  stats_vector_inc(stats_step,1010);
  stats_vector_display(stderr,stats_step,false,false,NULL);
  fprintf(stderr,"\n");
}

typedef struct {
  uint64_t count;
  uint64_t div_17;
} test_t;

void constructor_svector_load() {
  mm_slab_t* const slab = mm_slab_new_(BUFFER_SIZE_64M,BUFFER_SIZE_512M,MM_UNLIMITED_MEM,"");
  svector_t* const svector = svector_new(slab,test_t);
  svector_iterator_t iterator;

  // Writing
  svector_iterator_new(&iterator,svector,SVECTOR_WRITE_ITERATOR,0);
  uint64_t i;
  for (i=0;i<parameters.number;++i) {
    svector_iterator_get_element(&iterator,test_t)->count = i;
    svector_iterator_get_element(&iterator,test_t)->div_17 = i/17;
    // Next!
    svector_write_iterator_next(&iterator);
  }

  // Reading
  svector_iterator_new(&iterator,svector,SVECTOR_READ_ITERATOR,0);
  i=0;
  while (!svector_read_iterator_eoi(&iterator)) {
    test_t* const test = svector_iterator_get_element(&iterator,test_t);
    gem_cond_fatal_error_msg(test->count!=i,"Invalid count at %"PRIu64" :: count=%"PRIu64"\n",i,test->count);
    gem_cond_fatal_error_msg(test->div_17!=i/17,"Invalid DIV at %"PRIu64" :: div_17=%"PRIu64"\n",i,test->div_17);
    // Next!
    ++i;
    svector_read_iterator_next(&iterator);
  }
  gem_cond_fatal_error_msg(i!=parameters.number,"Wrong number of iterations (done %"PRIu64", should be %"PRIu64")",i,parameters.number);
}
//const uint64_t cdna_text_block_mask_left[256] =
//{
//  [0 ... 255] = 0xFFFFFFFFFFFFFFFFull,
//  [ 0] = 0x7FFFFFFFFFFFFFFFull,
//  [ 1] = 0x0FFFFFFFFFFFFFFFull,
//  [ 2] = 0x01FFFFFFFFFFFFFFull,
//  [ 3] = 0x003FFFFFFFFFFFFFull,
//  [ 4] = 0x0007FFFFFFFFFFFFull,
//  [ 5] = 0x0000FFFFFFFFFFFFull,
//  [ 6] = 0x00001FFFFFFFFFFFull,
//  [ 7] = 0x000003FFFFFFFFFFull,
//  [ 8] = 0x0000007FFFFFFFFFull,
//  [ 9] = 0x0000000FFFFFFFFFull,
//  [10] = 0x00000001FFFFFFFFull,
//  [11] = 0x000000003FFFFFFFull,
//  [12] = 0x0000000007FFFFFFull,
//  [13] = 0x0000000000FFFFFFull,
//  [14] = 0x00000000001FFFFFull,
//  [15] = 0x000000000003FFFFull,
//  [16] = 0x0000000000007FFFull,
//  [17] = 0x0000000000000FFFull,
//  [18] = 0x00000000000001FFull,
//  [19] = 0x000000000000003Full,
//  [20] = 0x0000000000000007ull,
//};
//const uint64_t cdna_text_block_mask_right[256] =
//{
//  [0 ... 255] = 0xFFFFFFFFFFFFFFFFull,
//  [ 0] = 0x8000000000000000ull,
//  [ 1] = 0xF000000000000000ull,
//  [ 2] = 0xFE00000000000000ull,
//  [ 3] = 0xFFC0000000000000ull,
//  [ 4] = 0xFFF8000000000000ull,
//  [ 5] = 0xFFFF000000000000ull,
//  [ 6] = 0xFFFFE00000000000ull,
//  [ 7] = 0xFFFFFC0000000000ull,
//  [ 8] = 0xFFFFFF8000000000ull,
//  [ 9] = 0xFFFFFFF000000000ull,
//  [10] = 0xFFFFFFFE00000000ull,
//  [11] = 0xFFFFFFFFC0000000ull,
//  [12] = 0xFFFFFFFFF8000000ull,
//  [13] = 0xFFFFFFFFFF000000ull,
//  [14] = 0xFFFFFFFFFFE00000ull,
//  [15] = 0xFFFFFFFFFFFC0000ull,
//  [16] = 0xFFFFFFFFFFFF8000ull,
//  [17] = 0xFFFFFFFFFFFFF000ull,
//  [18] = 0xFFFFFFFFFFFFFE00ull,
//  [19] = 0xFFFFFFFFFFFFFFC0ull,
//  [20] = 0xFFFFFFFFFFFFFFF8ull,
//};
//void constructor_show_cdna_text_block_mask() {
//  uint64_t i;
//  fprintf(stderr,"RIGHT\n");
//  for (i=0;i<21;++i) {
//    fprintf_uint64_binary(stderr,cdna_text_block_mask_left[i]);
//    fprintf(stderr,"\n");
//  }
//  fprintf(stderr,"LEFT\n");
//  for (i=0;i<21;++i) {
//    fprintf_uint64_binary(stderr,cdna_text_block_mask_right[i]);
//    fprintf(stderr,"\n");
//  }
//}
void constructor_sparse_bitmap_test() {
  mm_slab_t* const slab = mm_slab_new(BUFFER_SIZE_16M);

  /*
   * Create & add some elements
   */
  sparse_bitmap_builder_t* const sparse_bitmap_builder = sparse_bitmap_builder_new(slab);
  sparse_bitmap_builder_add_bitmap(sparse_bitmap_builder,1ull);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_add_bitmap(sparse_bitmap_builder,3ull);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); // 4
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); // 5
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); // 6
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); // 7
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); // 8
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); // 9
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); // 10

  sparse_bitmap_builder_add_bitmap(sparse_bitmap_builder,11ull); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);

  sparse_bitmap_builder_add_bitmap(sparse_bitmap_builder,21ull); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);

  sparse_bitmap_builder_add_bitmap(sparse_bitmap_builder,31ull); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);

  sparse_bitmap_builder_add_bitmap(sparse_bitmap_builder,41ull); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);

  sparse_bitmap_builder_add_bitmap(sparse_bitmap_builder,51ull); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);

  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder); sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_skip_bitmap(sparse_bitmap_builder);
  sparse_bitmap_builder_add_bitmap(sparse_bitmap_builder,64ull);

  sparse_bitmap_builder_add_bitmap(sparse_bitmap_builder,65ull);

  /*
   * Write
   */
  fm_t* file = fm_open_file("test.sbm",FM_WRITE);
  sparse_bitmap_builder_write(file,sparse_bitmap_builder);
  fm_close(file);

  /*
   * Read
   */
  file = fm_open_file("test.sbm",FM_READ);
  sparse_bitmap_t* const sparse_bitmap = sparse_bitmap_read(file);

  /*
   * Show some stats
   */
  sparse_bitmap_print(stderr,sparse_bitmap,true);

  /*
   * Show a random element
   */
//  fprintf(stderr,"Pos %"PRIu64" Contained %"PRIu64", bitmap %"PRIu64"\n",
//      63ull,sparse_bitmap_is_contained(sparse_bitmap,63ull),sparse_bitmap_get_bitmap(sparse_bitmap,63ull));
//  fprintf(stderr,"Pos %"PRIu64" Contained %"PRIu64", bitmap %"PRIu64"\n",
//      2ull,sparse_bitmap_is_contained(sparse_bitmap,2ull),sparse_bitmap_get_bitmap(sparse_bitmap,2ull));
}
void constructor_cdna_test() {
  mm_slab_t* const slab = mm_slab_new(BUFFER_SIZE_16M);
  cdna_text_t* const cdna_text = cdna_text_new(slab);

  // 0
  cdna_text_add_char(cdna_text,dna_encode('A')); /*  1 */
  cdna_text_add_char(cdna_text,dna_encode('C')); /*  2 */
  cdna_text_add_char(cdna_text,dna_encode('G')); /*  3 */
  cdna_text_add_char(cdna_text,dna_encode('T')); /*  4 */
  cdna_text_add_char(cdna_text,dna_encode('N')); /*  5 */
  cdna_text_add_char(cdna_text,dna_encode('|')); /*  6 */

  // 6
  cdna_text_add_char(cdna_text,dna_encode('A')); /*  1 */
  cdna_text_add_char(cdna_text,dna_encode('C')); /*  2 */
  cdna_text_add_char(cdna_text,dna_encode('G')); /*  3 */
  cdna_text_add_char(cdna_text,dna_encode('T')); /*  4 */
  cdna_text_add_char(cdna_text,dna_encode('N')); /*  5 */
  cdna_text_add_char(cdna_text,dna_encode('|')); /*  6 */

  // 12
  cdna_text_add_char(cdna_text,dna_encode('A')); /*  1 */
  cdna_text_add_char(cdna_text,dna_encode('C')); /*  2 */
  cdna_text_add_char(cdna_text,dna_encode('G')); /*  3 */
  cdna_text_add_char(cdna_text,dna_encode('T')); /*  4 */
  cdna_text_add_char(cdna_text,dna_encode('N')); /*  5 */
  cdna_text_add_char(cdna_text,dna_encode('|')); /*  6 */

  // 18
  cdna_text_add_char(cdna_text,dna_encode('A')); /*  1 */
  cdna_text_add_char(cdna_text,dna_encode('C')); /*  2 */
  cdna_text_add_char(cdna_text,dna_encode('G')); /*  3 */
  cdna_text_add_char(cdna_text,dna_encode('T')); /*  4 */
  cdna_text_add_char(cdna_text,dna_encode('N')); /*  5 */
  cdna_text_add_char(cdna_text,dna_encode('|')); /*  6 */

  // 24
  cdna_text_add_char(cdna_text,dna_encode('A')); /*  1 */
  cdna_text_add_char(cdna_text,dna_encode('C')); /*  2 */
  cdna_text_add_char(cdna_text,dna_encode('G')); /*  3 */
  cdna_text_add_char(cdna_text,dna_encode('T')); /*  4 */
  cdna_text_add_char(cdna_text,dna_encode('N')); /*  5 */
  cdna_text_add_char(cdna_text,dna_encode('|')); /*  6 */

  // 30
  cdna_text_add_char(cdna_text,dna_encode('A')); /*  1 */
  cdna_text_add_char(cdna_text,dna_encode('C')); /*  2 */
  cdna_text_add_char(cdna_text,dna_encode('G')); /*  3 */
  cdna_text_add_char(cdna_text,dna_encode('T')); /*  4 */
  cdna_text_add_char(cdna_text,dna_encode('N')); /*  5 */
  cdna_text_add_char(cdna_text,dna_encode('|')); /*  6 */

  // Close
  cdna_text_close(cdna_text);

  // Iterator
  cdna_text_iterator_t iterator;
  cdna_text_iterator_init(&iterator,cdna_text,0);
  while (!cdna_text_iterator_eoi(&iterator)) {
    fprintf(stderr,"%c",dna_decode(cdna_text_iterator_get_char_encoded(&iterator)));
    cdna_text_iterator_next_char(&iterator);
  }
}
void constructor_locator_test() {
  mm_slab_t* const slab = mm_slab_new(BUFFER_SIZE_16M);

  /*
   * Build some example
   */
  locator_builder_t* const locator_builder = locator_builder_new(slab);

//  locator_builder_add_sequence(locator_builder,"AAA",3);
//  locator_builder_add_interval(locator_builder,0,0,100,100,locator_interval_regular);
//  locator_builder_add_interval(locator_builder,0,0,200,200,locator_interval_regular);
//
//  locator_builder_add_sequence(locator_builder,"BBB",3);
//  locator_builder_add_interval(locator_builder,0,0,100,100,locator_interval_regular);
////  locator_builder_skip_text(locator_builder,100);
//  locator_builder_add_interval(locator_builder,0,0,200,200,locator_interval_regular);

  /*
   * Print some stats
   */
  locator_builder_print(stderr,locator_builder,true);

  /*
   * Write
   */
  fm_t* file = fm_open_file("test.loc",FM_WRITE);
  locator_builder_write(file,locator_builder);
  fm_close(file);

  /*
   * Read
   */
  file = fm_open_file("test.loc",FM_READ);
//  locator_t* const locator = locator_read(file);
//  locator_print(stderr,locator,true);

}
void constructor_packed_integer_array(const uint64_t int_length) {
  packed_integer_array_t* const array = packed_integer_array_new(64,int_length);
  fprintf(stderr,">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
                 " Printing %"PRIu64" bits-length\n",int_length);
  uint64_t i, integer;
  for (i=1,integer=1;i<int_length;++i,integer<<=1) {
    fprintf(stderr,"%"PRIu64"\n",integer);
    packed_integer_array_store(array,i-1,integer);
  }
  packed_integer_array_print(stderr,array,true);
  packed_integer_array_delete(array);
}
void constructor_packed_integer_array_test() {
  constructor_packed_integer_array(12);
  constructor_packed_integer_array(13);
  constructor_packed_integer_array(18);
  constructor_packed_integer_array(28);
  constructor_packed_integer_array(37);
  constructor_packed_integer_array(49);
  constructor_packed_integer_array(63);
  constructor_packed_integer_array(64);
}
void constructor_sparse_array_locator_test() {
  sparse_array_locator_t* sal[2];

  sal[0] = sparse_array_locator_new(0,64);
  sal[1] = sparse_array_locator_new(64,100);

  sparse_array_locator_mark(sal[0],2);
  sparse_array_locator_mark(sal[0],3);
  sparse_array_locator_mark(sal[0],47);
  sparse_array_locator_mark(sal[0],48);
  sparse_array_locator_mark(sal[0],63);

  sparse_array_locator_mark(sal[1],64);
  sparse_array_locator_mark(sal[1],66);
  sparse_array_locator_mark(sal[1],80);

  fm_t* file = fm_open_file("sal_test.sal",FM_WRITE);
  sparse_array_locator_merge__write(file,sal,2);
  fm_close(file);
  sparse_array_locator_delete(sal[0]);
  sparse_array_locator_delete(sal[1]);

  file = fm_open_file("sal_test.sal",FM_READ);
  sparse_array_locator_t* const locator = sparse_array_locator_read(file);

  sparse_array_locator_print(stderr,locator,true);
  if (sparse_array_locator_is_marked(locator,2) &&
      sparse_array_locator_is_marked(locator,3) &&
      sparse_array_locator_is_marked(locator,47) &&
      sparse_array_locator_is_marked(locator,48) &&
      sparse_array_locator_is_marked(locator,63) &&
      sparse_array_locator_is_marked(locator,64) &&
      sparse_array_locator_is_marked(locator,66) &&
      sparse_array_locator_is_marked(locator,80)) {
    fprintf(stderr,"OK - passed (%"PRIu64" == 7)",sparse_array_locator_get_erank(locator,80));
  }
  sparse_array_locator_delete(locator);
}
void constructor_cdna_text__reverse() {
  mm_slab_t* const slab = mm_slab_new(BUFFER_SIZE_16M);
  cdna_text_t* const text = cdna_text_new(slab);

  cdna_text_add_char(text,1); // 1
  cdna_text_add_char(text,1); //
  cdna_text_add_char(text,2); //
  cdna_text_add_char(text,2); //
  cdna_text_add_char(text,3); //
  cdna_text_add_char(text,3); //
  cdna_text_add_char(text,4); //
  cdna_text_add_char(text,4); //
  cdna_text_add_char(text,5); //
  cdna_text_add_char(text,5); //
  cdna_text_add_char(text,1); //
  cdna_text_add_char(text,1); //
  cdna_text_add_char(text,1); //
  cdna_text_add_char(text,1); //
  cdna_text_add_char(text,2); //
  cdna_text_add_char(text,2); //
  cdna_text_add_char(text,2); //
  cdna_text_add_char(text,2); //
  cdna_text_add_char(text,3); //
  cdna_text_add_char(text,3); //
  cdna_text_add_char(text,3); // 21

  cdna_text_add_char(text,3); // 22
  cdna_text_add_char(text,4); // 23
  cdna_text_add_char(text,4); // 24
  cdna_text_add_char(text,4); // 25
  cdna_text_add_char(text,4); // 26

//  fprintf(stderr,"Forward\n");
  cdna_text_iterator_t it;
//  cdna_text_iterator_init(&it,text,0);
//  while (!cdna_text_iterator_eoi(&it)) {
//    const uint8_t enc = cdna_text_iterator_get_char_encoded(&it);
//    fprintf(stderr,"%c",dna_decode(enc));
//    cdna_text_iterator_next_char(&it);
//  }
//  fprintf(stderr,"\n");

  int64_t i;
  for (i=text->text_length-1;i>=0;--i) {
    fprintf(stderr,"[%03"PRIu64"]Reverse  ",i);
    cdna_text_reverse_iterator_init(&it,text,i);
    while (!cdna_text_reverse_iterator_eoi(&it)) {
      const uint8_t enc = cdna_text_reverse_iterator_get_char_encoded(&it);
      fprintf(stderr,"%c",dna_decode(enc));
      cdna_text_reverse_iterator_next_char(&it);
    }
    fprintf(stderr,"\n");
  }

}
void constructor_fast_mapper_setup() {

}
void constructor_priority_queue() {
  pqueue_t* const queue = pqueue_new(15);

  pqueue_push(queue,NULL,3);
  pqueue_push(queue,NULL,14);
  pqueue_clear(queue);

  pqueue_push(queue,NULL,3);
  pqueue_push(queue,NULL,14);
  pqueue_push(queue,NULL,1);
  pqueue_push(queue,NULL,2);
  pqueue_push(queue,NULL,4);
  pqueue_push(queue,NULL,11);
  pqueue_push(queue,NULL,15);
  pqueue_push(queue,NULL,16);
  pqueue_push(queue,NULL,5);
  pqueue_push(queue,NULL,5);
  pqueue_push(queue,NULL,6);
  pqueue_push(queue,NULL,7);
  pqueue_push(queue,NULL,13);
  pqueue_push(queue,NULL,0);
  pqueue_push(queue,NULL,12);
  pqueue_push(queue,NULL,9);
  pqueue_push(queue,NULL,10);
  pqueue_push(queue,NULL,8);

  while (!pqueue_is_empty(queue)) {
    printf("Pqueue\t%"PRIu64"\n",pqueue_top_priority(queue));
    pqueue_pop(queue,void);
  }

  pqueue_delete(queue);
}
void constructor_itoa() {
  char buffer[200], check[200];
  uint64_t i;
  for (i=0;i<UINT64_MAX;++i) {
    // GEM itoa
    const uint64_t num_digits = integer_to_ascii(buffer,i);
    // Check
    sprintf(check,"%"PRIu64"",i);
    // CMP
    if (strncmp(check,buffer,num_digits)!=0) {
      fprintf(stderr,"Error at %"PRIu64" (%s,%s)\n",i,check,buffer);
    }
  }
  fprintf(stderr,"All good!!\n");
}
void constructor_swg() {
//  uint64_t match_position=100, cigar_length=0;
//  int64_t effective_length=0;
//  int32_t alignment_score=0;
//  vector_t* const cigar_vector = vector_new(100,cigar_element_t);
//  swg_penalties_t swg_penalties;
//  mm_stack_t* mm_stack = mm_stack_new(mm_slab_new(BUFFER_SIZE_8M));
//  swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_A] = +1;
//  swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_C] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_G] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_T] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_A][ENC_DNA_CHAR_N] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_A] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_C] = +1;
//  swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_G] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_T] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_C][ENC_DNA_CHAR_N] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_A] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_C] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_G] = +1;
//  swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_T] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_G][ENC_DNA_CHAR_N] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_A] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_C] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_G] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_T] = +1;
//  swg_penalties.matching_score[ENC_DNA_CHAR_T][ENC_DNA_CHAR_N] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_N][ENC_DNA_CHAR_A] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_N][ENC_DNA_CHAR_C] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_N][ENC_DNA_CHAR_G] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_N][ENC_DNA_CHAR_T] = -4;
//  swg_penalties.matching_score[ENC_DNA_CHAR_N][ENC_DNA_CHAR_N] = -4;
//  swg_penalties.generic_match_score = 1;
//  swg_penalties.generic_mismatch_score = -4;
//  swg_penalties.gap_open_score = -6;
//  swg_penalties.gap_extension_score = -1;
//  // Regular
//  uint8_t* const key = (uint8_t*)"ACGTACGT";
//  uint8_t* const text = (uint8_t*)"ACGTACGT";
//  const uint64_t key_length = strlen((char*)key);
//  const uint64_t text_length = strlen((char*)text);
//  uint8_t* const key_enc = malloc(key_length*sizeof(char));
//  uint8_t* const text_enc = malloc(text_length*sizeof(char));
//  uint64_t i;
//  for (i=0;i<key_length;++i) key_enc[i] = dna_encode(key[i]);
//  for (i=0;i<text_length;++i) text_enc[i] = dna_encode(text[i]);
//  // SWG
//  swg_align_match_base(key_enc,key_length,&swg_penalties,&match_position,
//      text_enc,text_length,cigar_vector,&cigar_length,
//      &effective_length,&alignment_score,mm_stack);
//  printf("\n");
//  // SWG SIMD
//  swg_query_profile_t swg_query_profile;
//  swg_init_query_profile(&swg_query_profile,&swg_penalties,key_length,mm_stack);
//  swg_align_match_int16_simd128(key_enc,key_length,&swg_query_profile,&swg_penalties,
//      &match_position,text_enc,text_length,true,true,cigar_vector,&cigar_length,
//      &effective_length,&alignment_score,mm_stack);
//  swg_align_match_int16_simd128(key_enc,key_length,&swg_query_profile,&swg_penalties,
//      &match_position,text_enc,text_length,true,true,cigar_vector,&cigar_length,
//      &effective_length,&alignment_score,mm_stack);
//  // Free
//  vector_delete(cigar_vector);
}
/*
 *
 */
void constructor_nsearch_region_permutations_n(
    region_profile_t* const region_profile,const uint64_t current_region,
    const uint64_t offset,const uint64_t left_length,
    uint8_t* const key,const uint64_t key_length,
    const uint64_t max_error,mm_stack_t* const mm_stack) {
  uint64_t i;
  if (left_length==0) return;
  if (current_region+1 == region_profile->num_filtering_regions) {
    // Last
    region_profile->filtering_region[current_region].begin = offset;
    region_profile->filtering_region[current_region].end = key_length;
    region_profile->filtering_region[current_region].min = 0;
    region_profile->filtering_region[current_region].max = max_error;
    // Print partition
    fprintf(stderr,"(");
    for (i=0;i<region_profile->num_filtering_regions;++i) {
      fprintf(stderr,"%s%02lu",(i==0) ? "" : ",",
          region_profile->filtering_region[i].end-region_profile->filtering_region[i].begin);
    }
    fprintf(stderr,")\t");
    // Search
    nsearch_hamming_preconditioned(NULL,region_profile,key,key_length,max_error,NULL,mm_stack);
  } else {
    for (i=1;i<left_length;++i) {
      region_profile->filtering_region[current_region].begin = offset;
      region_profile->filtering_region[current_region].end = offset+i;
      region_profile->filtering_region[current_region].min = 0;
      region_profile->filtering_region[current_region].max = max_error;
      // Next
      constructor_nsearch_region_permutations_n(region_profile,
          current_region+1,offset+i,left_length-i,key,key_length,max_error,mm_stack);
    }
  }
}
void constructor_nsearch_region_permutations() {
  // MM-Stack
  mm_slab_t* const slab = mm_slab_new(BUFFER_SIZE_16M);
  mm_stack_t* const mm_stack = mm_stack_new(slab);
  // Search Parameters
  const uint64_t max_error = parameters.number;
  const char* const key = parameters.name_input_file; // "ACGGTGACAAACGTCACGGTGGCACATGCAACCAAAAGCTG";//"AAAAAAA";
  const uint64_t key_length = strlen(key);
  // Encode key
  uint8_t* const enc_key = malloc(key_length+1);
  uint64_t i;
  for (i=0;i<key_length;++i) enc_key[i] = dna_encode(key[i]);
  enc_key[key_length] = '\0';
  // Configure Region Profile
  const uint64_t num_filtering_regions = 4;
  region_profile_t region_profile;
  region_profile.num_filtering_regions = num_filtering_regions;
  region_profile.filtering_region = mm_stack_calloc(mm_stack,num_filtering_regions,region_search_t,true);
  // Generate all possible partitions
  constructor_nsearch_region_permutations_n(&region_profile,0,0,key_length,enc_key,key_length,max_error,mm_stack);
}
/*
 *
 */
void constructor_ns_hamming_brute() {
  // MM-Stack
  mm_slab_t* const slab = mm_slab_new(BUFFER_SIZE_16M);
  mm_stack_t* const mm_stack = mm_stack_new(slab);
  // Search Parameters
  const uint64_t max_error = parameters.number;
  const char* const key = parameters.name_input_file; // "ACGGTGACAAACGTCACGGTGGCACATGCAACCAAAAGCTG";//"AAAAAAA";
  const uint64_t key_length = strlen(key);
  // Encode key
  uint8_t* const enc_key = malloc(key_length+1);
  uint64_t i;
  for (i=0;i<key_length;++i) enc_key[i] = dna_encode(key[i]);
  enc_key[key_length] = '\0';
  // Search
  nsearch_hamming_brute_force(NULL,enc_key,key_length,max_error,NULL,mm_stack);
}
void constructor_ns_hamming() {
  // MM-Stack
  mm_slab_t* const slab = mm_slab_new(BUFFER_SIZE_16M);
  mm_stack_t* const mm_stack = mm_stack_new(slab);
  // Search Parameters
  const uint64_t max_error = parameters.number;
  const char* const key = parameters.name_input_file; // "ACGGTGACAAACGTCACGGTGGCACATGCAACCAAAAGCTG";//"AAAAAAA";
  const uint64_t key_length = strlen(key);
  // Encode key
  uint8_t* const enc_key = malloc(key_length+1);
  uint64_t i;
  for (i=0;i<key_length;++i) enc_key[i] = dna_encode(key[i]);
  enc_key[key_length] = '\0';
  // Search
  nsearch_hamming(NULL,enc_key,key_length,max_error,NULL,mm_stack);
}
void constructor_ns_hamming_2regions() {
  // MM-Stack
  mm_slab_t* const slab = mm_slab_new(BUFFER_SIZE_16M);
  mm_stack_t* const mm_stack = mm_stack_new(slab);
  // Search Parameters
  const uint64_t max_error = parameters.number;
  const char* const key = parameters.name_input_file; // "ACGGTGACAAACGTCACGGTGGCACATGCAACCAAAAGCTG";//"AAAAAAA";
  const uint64_t key_length = strlen(key);
  // Encode key
  uint8_t* const enc_key = malloc(key_length+1);
  uint64_t i;
  for (i=0;i<key_length;++i) enc_key[i] = dna_encode(key[i]);
  enc_key[key_length] = '\0';
  // Configure Region Profile
  region_profile_t region_profile;
  region_profile.filtering_region = mm_stack_calloc(mm_stack,2,region_search_t,true);
  region_profile.filtering_region[0].begin = 0;
  region_profile.filtering_region[0].end = 4;
  region_profile.filtering_region[0].min = 1;
  region_profile.filtering_region[0].max = max_error;
  region_profile.filtering_region[1].begin = 4;
  region_profile.filtering_region[1].end = key_length;
  region_profile.filtering_region[1].min = 0;
  region_profile.filtering_region[1].max = max_error-1;
  region_profile.num_filtering_regions = 2;
  // Search
  nsearch_hamming_preconditioned(NULL,&region_profile,enc_key,key_length,max_error,NULL,mm_stack);
}
void constructor_ns_hamming_permutations() {
  // MM-Stack
  mm_slab_t* const slab = mm_slab_new(BUFFER_SIZE_16M);
  mm_stack_t* const mm_stack = mm_stack_new(slab);
  // Search Parameters
  const uint64_t max_error = parameters.number;
  const char* const key = parameters.name_input_file; // "ACGGTGACAAACGTCACGGTGGCACATGCAACCAAAAGCTG";//"AAAAAAA";
  const uint64_t key_length = strlen(key);
  // Encode key
  uint8_t* const enc_key = malloc(key_length+1);
  uint64_t i;
  for (i=0;i<key_length;++i) enc_key[i] = dna_encode(key[i]);
  enc_key[key_length] = '\0';
  // Configure Region Profile
  const uint64_t num_filtering_regions = parameters.param1;
  region_profile_t region_profile;
  region_profile.num_filtering_regions = num_filtering_regions;
  region_profile.filtering_region = mm_stack_calloc(mm_stack,num_filtering_regions,region_search_t,true);
  // Generate all possible partitions
  constructor_nsearch_region_permutations_n(&region_profile,0,0,key_length,enc_key,key_length,max_error,mm_stack);
}
void constructor_ns_edit_brute() {
  // MM-Stack
  mm_slab_t* const slab = mm_slab_new(BUFFER_SIZE_16M);
  mm_stack_t* const mm_stack = mm_stack_new(slab);
  // Search Parameters
  const uint64_t max_error = parameters.number;
  const char* const key = parameters.name_input_file;
  const uint64_t key_length = strlen(key);
  // Encode key
  uint8_t* const enc_key = malloc(key_length+1);
  uint64_t i;
  for (i=0;i<key_length;++i) enc_key[i] = dna_encode(key[i]);
  enc_key[key_length] = '\0';
  // Search
  nsearch_levenshtein_brute_force(NULL,enc_key,key_length,max_error,NULL,mm_stack);
}
void constructor_ns_edit_partition() {
  // MM-Stack
  mm_slab_t* const slab = mm_slab_new(BUFFER_SIZE_16M);
  mm_stack_t* const mm_stack = mm_stack_new(slab);
  // Search Parameters
  const uint64_t max_error = parameters.number;
  const char* const key = parameters.name_input_file;
  const uint64_t key_length = strlen(key);
  // Encode key
  uint8_t* const enc_key = malloc(key_length+1);
  uint64_t i;
  for (i=0;i<key_length;++i) enc_key[i] = dna_encode(key[i]);
  enc_key[key_length] = '\0';
  // Search
  nsearch_levenshtein(NULL,enc_key,key_length,max_error,NULL,mm_stack);
}
/*
 * LCS
 */
void constructor_lsc() {
  // MM-Stack
  mm_slab_t* const slab = mm_slab_new(BUFFER_SIZE_16M);
  mm_stack_t* const mm_stack = mm_stack_new(slab);
  // Text & key
  char* text = "ACAAGTA";
  char* key =  "ACGT";
//  char* text = "TTCAGATAGCCCTTAAAAGGAGTTTCATCATCTATACGGGAGGTTATCTAACAAAAATGCTATTTCTGTTTGAAAGAAACTATGTAAAAATTACAATTA";
//  char* key =  "AGATAGCCCTTCAAAGGAGTTTCATCATCTTTACGGGAGGTTATCTAACAAAAATGCTATTTCTGTTTGAAAGAAACTATGTAAAAATTACAATTAACGA";
  const uint64_t key_length = strlen(key);
  const uint64_t text_length = strlen(text);
  // Encode key
  uint64_t i;
  uint8_t* const enc_key = malloc(key_length+1);
  for (i=0;i<key_length;++i) enc_key[i] = dna_encode(key[i]);
  enc_key[key_length] = '\0';
  // Encode text
  uint8_t* const enc_text = malloc(text_length+1);
  for (i=0;i<text_length;++i) enc_text[i] = dna_encode(text[i]);
  enc_text[text_length] = '\0';

//  // Compute LCS
//  uint64_t lcs_distance, match_end_column;
//  align_ond_compute_lcs_distance(enc_key,key_length,enc_text,text_length,&lcs_distance,&match_end_column,0,mm_stack);
//  printf("  => LCS %lu (col=%lu)\n",lcs_distance,match_end_column);

  // Prepare input
  match_align_input_t align_input;
  align_input.key = enc_key;
  align_input.key_length = key_length;
  align_input.text = enc_text;
  align_input.text_length = text_length;
  const uint64_t max_distance = key_length+text_length;
  // Align O(ND)
  vector_t* cigar_vector = vector_new(100,cigar_element_t);
  match_alignment_t match_alignment;
  match_alignment.match_position = 0;
  align_ond_match(&align_input,max_distance,&match_alignment,cigar_vector,mm_stack);
  // Display
  match_alignment_print_pretty(stderr,&match_alignment,
      cigar_vector,enc_key,key_length,
      enc_text,text_length,mm_stack);
//      enc_text+match_alignment.match_position,text_length,mm_stack);
}
/*
 * Generic Menu
 */
void usage() {
  fprintf(stderr, "USE: ./gem-tools-examples -i input -o output \n"
                  "      Options::\n"
                  "        --input|i      <File>\n"
                  "        --output|o     <File>\n"
                  "        --select|s     <Number>\n"
                  "        --number|n     <Number>\n"
                  "        --tmp-folder   <Path>\n"
                  "        --help|h\n");
}
void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "select", required_argument, 0, 's' },
    { "number", required_argument, 0, 'n' },
    // Parameters
    { "param1", required_argument, 0, '1' },
    { "param2", required_argument, 0, '2' },
    { "param3", required_argument, 0, '3' },
    { "param4", required_argument, 0, '4' },
    { "param5", required_argument, 0, '5' },
    { "param6", required_argument, 0, '6' },
    { "param7", required_argument, 0, '7' },
    // Temporary folder
    { "tmp-folder", required_argument, 0, 100},
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  while (1) {
    c=getopt_long(argc,argv,"i:o:s:n:1:2:3:4:5:6:7:T:h",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    case 'i':
      parameters.name_input_file = optarg;
      break;
    case 'o':
      parameters.name_output_file = optarg;
      break;
    case 's':
      parameters.option = optarg;
      break;
    case 'n':
     parameters.number = atol(optarg);
     break;
    /* Numbered Params */
    case '1':
     parameters.param1 = atol(optarg);
     break;
    case '2':
     parameters.param2 = atol(optarg);
     break;
    case '3':
     parameters.param3 = atol(optarg);
     break;
    case '4':
     parameters.param4 = atol(optarg);
     break;
    case '5':
     parameters.param5 = atol(optarg);
     break;
    case '6':
     parameters.param6 = atol(optarg);
     break;
    case '7':
     parameters.param7 = atol(optarg);
     break;
    case '8':
      mm_set_tmp_folder(optarg);
      break;
    case 'T':
      parameters.num_threads = atol(optarg);
      break;
    case 'M': // --max-memory
      gem_cond_fatal_error_msg(input_text_parse_size(optarg,&(parameters.max_memory)),"Wrong Size '%s'",optarg);
      break;
    case 100: // --tmp-folder
      parameters.tmp_folder = optarg;
      break;
    case 'h':
      usage();
      exit(1);
    case '?': default:
      fprintf(stderr, "Option not recognized \n"); exit(1);
    }
  }
  /* System */
  if (parameters.tmp_folder!=NULL) mm_set_tmp_folder(parameters.tmp_folder);
}
int main(int argc,char** argv) {
  // GT error handler
  gem_handle_error_signals();

  // Parsing command-line options
  parse_arguments(argc,argv);

  // Load!
  // constructor_write_fm();
  // constructor_stats_vector();
  // constructor_svector_load();
  // constructor_show_cdna_text_block_mask();
  // constructor_sparse_bitmap_test();
  // constructor_cdna_test();
  // constructor_cdna_bitwise_test();
  // constructor_sparse_bitmap_test();
  // constructor_locator_test();

  // constructor_packed_integer_array_test();
  // constructor_packed_integer_arrays_test();
  // constructor_sparse_array_locator_test();
  // constructor_cdna_text__reverse();
  // constructor_packed_integer_arrays_test_bis();
  // constructor_fast_mapper_setup();
  // constructor_priority_queue();
  //  constructor_itoa();
  // constructor_swg();

  //  if (gem_strcaseeq(parameters.option,"hamming-brute")) {
  //    constructor_ns_hamming_brute();
  //  }
  //  if (gem_strcaseeq(parameters.option,"hamming-partition")) {
  //    constructor_ns_hamming();
  //  }
  //  if (gem_strcaseeq(parameters.option,"hamming-regions")) {
  //    constructor_ns_hamming_2regions();
  //  }
  //  if (gem_strcaseeq(parameters.option,"hamming-permutations")) {
  //    constructor_ns_hamming_permutations();
  //  }
  //  if (gem_strcaseeq(parameters.option,"edit-brute")) {
  //    constructor_ns_edit_brute();
  //  }
  //  if (gem_strcaseeq(parameters.option,"edit-partition")) {
  //    constructor_ns_edit_partition();
  //  }

  //  constructor_lsc();
  //  return 0;

  archive_t* archive = archive_read("./test.gem",false);
  fprintf(stderr,"Position=%lu\tProyection=%lu\n",
      0ul,archive_text_get_unitary_projection(archive->text,0));

  return 0;
}

