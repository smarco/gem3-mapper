/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2026 by Jacob Schreiber  <jacob.schreiber@umassmed.edu>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Jacob Schreiber <jacob.schreiber@umassmed.edu> (unit tests, written with Claude Code)
 * DESCRIPTION:
 *   Area D: the ordered output file shared by the mapper threads.
 *
 *   Threads take block ids in increasing order (as the mapper threads take
 *   input blocks), write each block's records through their own buffered
 *   output file with random pauses, and return it; records of a block may
 *   overflow the output buffer and continue in extension buffers. The file
 *   must hold a header written with ofprintf followed by every block in id
 *   order, whatever the thread count, number of buffers (fewer than threads
 *   included), buffer size and compression (plain, gzip, bzip2). The
 *   expected bytes are generated independently from the block ids.
 *
 *   Internal signatures called (a speed iteration may change them; update
 *   this driver with them):
 *     output_file_t* output_file_new(char* file_name,uint64_t max_output_buffers,uint64_t buffer_size)
 *     output_file_t* output_gzip_stream_new(FILE*,uint64_t,uint64_t), output_bzip_stream_new(FILE*,uint64_t,uint64_t)
 *     void output_file_close(output_file_t*)
 *     int ofprintf(output_file_t*,const char*,...)
 *     buffered_output_file_t* buffered_output_file_new(output_file_t*)
 *     void buffered_output_file_request_buffer(buffered_output_file_t*,uint32_t block_id)
 *     void buffered_output_file_reserve(buffered_output_file_t*,uint64_t num_chars)
 *     void buffered_output_file_dump_buffer(buffered_output_file_t*), buffered_output_file_close(...)
 *     bofprintf_uint64/bofprintf_char/bofprintf_string(buffered_output_file_t*,...)
 *   Threads start through gt_run_threads (gem3_test.h), which registers gem
 *   thread ids and profile slots as gem-mapper does.
 */

#include "gem3_test.h"
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>
#include <bzlib.h>
#include "io/output_file.h"
#include "io/buffered_output_file.h"

#define NUM_BLOCKS 160
#define MAX_RECORD 300

typedef struct {
  buffered_output_file_t* bof;
  pthread_mutex_t* mutex;
  uint64_t* next_block;
  uint64_t seed;
} worker_t;

/* Record j of block b: deterministic text of 1..MAX_RECORD-1 chars + '\n' */
static size_t record(const uint64_t b,const uint64_t j,char* const out) {
  gt_rng_t rng = gt_rng_new(b*1000003+j);
  const uint64_t len = gt_rng_range(&rng,1,MAX_RECORD-20);
  size_t n = (size_t)sprintf(out,"%" PRIu64 ":%" PRIu64 ":",b,j);
  uint64_t i;
  for (i=0;i<len;++i) out[n++] = "ACGT"[gt_rng_below(&rng,4)];
  out[n++] = '\n';
  return n;
}
static uint64_t records_in_block(const uint64_t b) {
  gt_rng_t rng = gt_rng_new(b*7+3);
  return gt_rng_below(&rng,25); // includes empty blocks
}

static void* worker(void* arg) {
  worker_t* const w = arg;
  gt_rng_t rng = gt_rng_new(w->seed);
  char rec[MAX_RECORD+64];
  while (true) {
    pthread_mutex_lock(w->mutex);
    const uint64_t b = (*w->next_block)++;
    pthread_mutex_unlock(w->mutex);
    if (b >= NUM_BLOCKS) break;
    buffered_output_file_request_buffer(w->bof,(uint32_t)b);
    const uint64_t n = records_in_block(b);
    uint64_t j;
    for (j=0;j<n;++j) {
      const size_t len = record(b,j,rec);
      // Write a record in pieces, as the SAM writer does field by field:
      // the block id as a number, then the rest in random-sized strings
      const size_t digits = (size_t)snprintf(NULL,0,"%" PRIu64,b);
      buffered_output_file_reserve(w->bof,20);
      bofprintf_uint64(w->bof,b);
      size_t off = digits;
      while (off < len) {
        const size_t r = (size_t)gt_rng_range(&rng,1,40);
        const size_t piece = MIN(len-off,r);
        buffered_output_file_reserve(w->bof,piece);
        bofprintf_string(w->bof,(int)piece,rec+off);
        off += piece;
      }
      if (gt_rng_below(&rng,4)==0) {
        struct timespec ts = { 0, (long)gt_rng_below(&rng,60000) };
        nanosleep(&ts,NULL);
      }
    }
    buffered_output_file_dump_buffer(w->bof);
  }
  return NULL;
}

static char* expected(size_t* const length) {
  char* const out = malloc((size_t)NUM_BLOCKS*25*(MAX_RECORD+64)+64);
  size_t n = (size_t)sprintf(out,"@HD\tVN:test\n");
  uint64_t b, j;
  for (b=0;b<NUM_BLOCKS;++b) {
    const uint64_t m = records_in_block(b);
    for (j=0;j<m;++j) n += record(b,j,out+n);
  }
  *length = n;
  return out;
}

typedef enum { plain_file, gzip_stream, bzip_stream } kind_t;

static char* read_back(const char* const path,const kind_t kind,size_t* const length) {
  const size_t cap = (size_t)NUM_BLOCKS*25*(MAX_RECORD+64)+64;
  char* const buf = malloc(cap);
  size_t n = 0;
  if (kind==gzip_stream) {
    gzFile g = gzopen(path,"rb");
    if (g==NULL) { free(buf); return NULL; }
    int r;
    while ((r=gzread(g,buf+n,(unsigned)(cap-n))) > 0) n += (size_t)r;
    gzclose(g);
  } else if (kind==bzip_stream) {
    FILE* const f = fopen(path,"rb");
    int err;
    BZFILE* const bz = BZ2_bzReadOpen(&err,f,0,0,NULL,0);
    while (err==BZ_OK) {
      const int r = BZ2_bzRead(&err,bz,buf+n,(int)(cap-n));
      if (r > 0) n += (size_t)r;
    }
    BZ2_bzReadClose(&err,bz);
    fclose(f);
  } else {
    FILE* const f = fopen(path,"rb");
    if (f==NULL) { free(buf); return NULL; }
    n = fread(buf,1,cap,f);
    fclose(f);
  }
  *length = n;
  return buf;
}

static void run_case(
    const int num_threads,const uint64_t num_buffers,
    const uint64_t buffer_size,const kind_t kind) {
  char path[] = "/tmp/gem3-test-d-output-XXXXXX";
  const int fd = mkstemp(path);
  GT_REQUIRE(fd >= 0);
  close(fd);
  output_file_t* of;
  if (kind==plain_file) {
    of = output_file_new(path,num_buffers,buffer_size);
  } else {
    FILE* const stream = fopen(path,"wb");
    GT_REQUIRE(stream!=NULL);
    of = (kind==gzip_stream) ?
        output_gzip_stream_new(stream,num_buffers,buffer_size) :
        output_bzip_stream_new(stream,num_buffers,buffer_size);
  }
  ofprintf(of,"@HD\tVN:%s\n","test");
  pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  uint64_t next_block = 0;
  worker_t workers[16];
  int i;
  for (i=0;i<num_threads;++i) {
    workers[i].bof = buffered_output_file_new(of);
    workers[i].mutex = &mutex;
    workers[i].next_block = &next_block;
    workers[i].seed = (uint64_t)(i+1)*0x9E3779B97F4A7C15ull ^ buffer_size ^ num_buffers;
  }
  // Gem thread ids 1..n and their profile slots, as the mapper threads have
  GT_REQUIRE(gt_run_threads((uint64_t)num_threads,worker,workers,sizeof(worker_t))==0);
  for (i=0;i<num_threads;++i) buffered_output_file_close(workers[i].bof);
  output_file_close(of);
  size_t want_len, got_len;
  char* const want = expected(&want_len);
  char* const got = read_back(path,kind,&got_len);
  unlink(path);
  GT_REQUIRE_MSG(got!=NULL,"could not read %s back",path);
  size_t first = 0;
  while (first < want_len && first < got_len && want[first]==got[first]) ++first;
  GT_CHECK_MSG(got_len==want_len && first==want_len,
      "threads=%d buffers=%" PRIu64 " size=%" PRIu64 " kind=%d: %zu bytes, expected %zu, first difference at %zu",
      num_threads,num_buffers,buffer_size,(int)kind,got_len,want_len,first);
  free(want);
  free(got);
}

GT_TEST(blocks_in_order_one_thread) {
  run_case(1,1,4096,plain_file);
  run_case(1,4,BUFFER_SIZE_4M,plain_file);
}
GT_TEST(blocks_in_order_threads_and_buffers) {
  run_case(2,1,1024,plain_file);
  run_case(4,2,1024,plain_file);
  run_case(4,3,512,plain_file);
  run_case(8,2,2048,plain_file);
  run_case(8,16,BUFFER_SIZE_64K,plain_file);
  run_case(16,5,1024,plain_file);
}
GT_TEST(blocks_in_order_gzip) {
  run_case(1,2,1024,gzip_stream);
  run_case(6,3,1024,gzip_stream);
}
GT_TEST(blocks_in_order_bzip2) {
  run_case(1,2,1024,bzip_stream);
  run_case(6,3,1024,bzip_stream);
}

GT_CASES(
  GT_CASE(blocks_in_order_one_thread),
  GT_CASE(blocks_in_order_threads_and_buffers),
  GT_CASE(blocks_in_order_gzip),
  GT_CASE(blocks_in_order_bzip2)
)
