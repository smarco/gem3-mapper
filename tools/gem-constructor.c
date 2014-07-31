/*
 * PROJECT: GEMMapper
 * FILE: gem-constructor.c
 * DATE:5/12/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Loads/test things
 */

#include "gem_core.h"

/*
 * Generic parameters
 */
typedef struct {
  char *name_input_file;
  char *name_output_file;
  uint64_t option;
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
    .option=0,
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
    gem_cond_fatal_error_msg(test->count!=i,"Invalid count at %lu :: count=%lu\n",i,test->count);
    gem_cond_fatal_error_msg(test->div_17!=i/17,"Invalid DIV at %lu :: div_17=%lu\n",i,test->div_17);
    // Next!
    ++i;
    svector_read_iterator_next(&iterator);
  }
  gem_cond_fatal_error_msg(i!=parameters.number,"Wrong number of iterations (done %lu, should be %lu)",i,parameters.number);
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
  fprintf(stderr,"Pos %lu Contained %lu, bitmap %lu\n",
      63ul,(uint64_t)sparse_bitmap_is_contained(sparse_bitmap,63),sparse_bitmap_get_bitmap(sparse_bitmap,63));
  fprintf(stderr,"Pos %lu Contained %lu, bitmap %lu\n",
      2ul, (uint64_t)sparse_bitmap_is_contained(sparse_bitmap,2),sparse_bitmap_get_bitmap(sparse_bitmap,2));
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

  locator_builder_add_sequence(locator_builder,"AAA",3);
  locator_builder_add_interval(locator_builder,0,0,100,100,locator_interval_regular);
  locator_builder_add_interval(locator_builder,0,0,200,200,locator_interval_regular);

  locator_builder_add_sequence(locator_builder,"BBB",3);
  locator_builder_add_interval(locator_builder,0,0,100,100,locator_interval_regular);
//  locator_builder_skip_text(locator_builder,100);
  locator_builder_add_interval(locator_builder,0,0,200,200,locator_interval_regular);

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
  locator_t* const locator = locator_read(file);
  locator_print(stderr,locator,true);

}
GEM_INLINE void constructor_packed_integer_array(const uint64_t int_length) {
  packed_integer_array_t* const array = packed_integer_array_new(0,64,int_length);
  fprintf(stderr,">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
                 " Printing %lu bits-length\n",int_length);
  uint64_t i, integer;
  for (i=1,integer=1;i<int_length;++i,integer<<=1) {
    fprintf(stderr,"%lu\n",integer);
    packed_integer_array_store(array,i-1,integer);
  }
  packed_integer_array_print(stderr,array,true);
  packed_integer_array_delete(array);
}
GEM_INLINE void constructor_packed_integer_array_test() {
  constructor_packed_integer_array(12);
  constructor_packed_integer_array(13);
  constructor_packed_integer_array(18);
  constructor_packed_integer_array(28);
  constructor_packed_integer_array(37);
  constructor_packed_integer_array(49);
  constructor_packed_integer_array(63);
  constructor_packed_integer_array(64);
}
GEM_INLINE void constructor_sparse_array_locator_test() {
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
    fprintf(stderr,"OK - passed (%lu == 7)",sparse_array_locator_get_erank(locator,80));
  }
  sparse_array_locator_delete(locator);
}
GEM_INLINE void constructor_packed_integer_arrays_test() {
  packed_integer_array_t* pia[4];
  pia[0] = packed_integer_array_new(0,10,36);
  pia[1] = packed_integer_array_new(10,20,36);
  pia[2] = packed_integer_array_new(20,30,36);
  pia[3] = packed_integer_array_new(30,40,36);

  const uint64_t integer = ((((1ull<<35) | (1ull<<34)) | (1ull<1)) | 1ull);
  uint64_t i = 0;
  for (;i<10;++i) packed_integer_array_store(pia[0],i,integer);
  for (;i<20;++i) packed_integer_array_store(pia[1],i,integer);
  for (;i<30;++i) packed_integer_array_store(pia[2],i,integer);
  for (;i<40;++i) packed_integer_array_store(pia[3],i,integer);

  fm_t* file = fm_open_file("test.pia",FM_WRITE);
  packed_integer_array_sliced_write(file,pia,4);
  fm_close(file);
  packed_integer_array_delete(pia[0]);
  packed_integer_array_delete(pia[1]);
  packed_integer_array_delete(pia[2]);
  packed_integer_array_delete(pia[3]);

  file = fm_open_file("test.pia",FM_READ);
  packed_integer_array_t* const result = packed_integer_array_read(file);
  for (i=0;i<result->num_elements;++i) {
    fprintf(stderr,"%lu:%lu\n",i,packed_integer_array_load(result,i));
  }
}
GEM_INLINE void constructor_packed_integer_arrays_test_bis() {
  packed_integer_array_t* pia[4];
  pia[0] = packed_integer_array_new(0,2,6);
  pia[1] = packed_integer_array_new(2,4,6);
  pia[2] = packed_integer_array_new(4,6,6);
  pia[3] = packed_integer_array_new(6,8,6);

  uint64_t i = 0;
  for (;i<2;++i) packed_integer_array_store(pia[0],i,i);
  for (;i<4;++i) packed_integer_array_store(pia[1],i,i);
  for (;i<6;++i) packed_integer_array_store(pia[2],i,i);
  for (;i<8;++i) packed_integer_array_store(pia[3],i,i);

  fm_t* file = fm_open_file("test.pia",FM_WRITE);
  packed_integer_array_sliced_write(file,pia,4);
  fm_close(file);
  packed_integer_array_delete(pia[0]);
  packed_integer_array_delete(pia[1]);
  packed_integer_array_delete(pia[2]);
  packed_integer_array_delete(pia[3]);

  file = fm_open_file("test.pia",FM_READ);
  packed_integer_array_t* const result = packed_integer_array_read(file);
  for (i=0;i<result->num_elements;++i) {
    fprintf(stderr,"%lu:%lu\n",i,packed_integer_array_load(result,i));
  }
}
GEM_INLINE void constructor_cdna_text__reverse() {
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
    fprintf(stderr,"[%03lu]Reverse  ",i);
    cdna_text_reverse_iterator_init(&it,text,i);
    while (!cdna_text_reverse_iterator_eoi(&it)) {
      const uint8_t enc = cdna_text_reverse_iterator_get_char_encoded(&it);
      fprintf(stderr,"%c",dna_decode(enc));
      cdna_text_reverse_iterator_next_char(&it);
    }
    fprintf(stderr,"\n");
  }

}

GEM_INLINE void constructor_fast_mapper_setup() {

//  // GEM Runtime setup
//  gem_runtime_init(parameters.num_threads,0,NULL,NULL);
//
//  // Load GEM-archive
//  archive_t* const gem_archive = archive_read(parameters.name_input_file,false,true);




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
      parameters.option = atol(optarg);
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

  constructor_fast_mapper_setup();

  return 0;
}

