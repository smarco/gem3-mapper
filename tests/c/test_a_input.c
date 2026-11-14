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
 *   Read input: FASTA/FASTQ records parsed through the sliced and buffered
 *   input files (plain, gzip, bzip2, streams and pipes; block and buffer
 *   sizes down to a few bytes, so records and lines straddle buffers),
 *   compared with a reference parser in this driver; command-line value
 *   parsers (input_text.c); the quality model; file and memory managers
 *   used to read inputs and load the index (fm.c, mm.c).
 *
 *   Reference parsing rules (FASTQ/FASTA as gem-mapper documents them):
 *   the tag runs from after '@' or '>' to the first space or tab and loses
 *   a trailing /1, /2 or /3 (end 1, 2, 2); a following CASAVA 1.8 field
 *   (three colons, "1:N:0:ACGT") gives the end from its first character;
 *   any other field ending in /1, /2 or /3 is appended to the tag after an
 *   underscore. Read bases A, C, G, T, N in either case are kept, anything
 *   else becomes 'N'. A trailing '\r' is dropped from every line.
 *
 *   Internal signatures called (a speed iteration may change them; update
 *   this driver with them):
 *     input_file_sliced_open(char* name,uint64_t num_blocks,uint64_t block_size)
 *     input_file_sliced_popen(char* command,uint64_t num_blocks,uint64_t block_size)
 *     input_stream_sliced_open / input_gzip_stream_sliced_open / input_bzip_stream_sliced_open(FILE*,uint64_t,uint64_t)
 *     input_file_sliced_close(input_file_sliced_t*), input_file_sliced_get_file_name(input_file_sliced_t*)
 *     buffered_input_file_new(input_file_sliced_t*,uint64_t prefered_read_size), buffered_input_file_close()
 *     buffered_input_file_reload(buffered_input_file_t*,uint64_t forced_read_lines), buffered_input_file_eob()
 *     buffered_input_file_get_file_name()
 *     input_fasta_parse_sequence(buffered_input_file_t*,sequence_t*,bool check_input_buffer)
 *     sequence_init(sequence_t*,bisulfite_read_t,mm_allocator_t*), sequence_get_* / sequence_has_*
 *     sequence_qualities_model_process(sequence_qualities_format_t,char* qualities,uint64_t length,uint8_t* mask)
 *     input_text_parse_integer/_double(const char**,...), input_text_parse_size(char*,uint64_t*),
 *     input_text_parse_csv_arguments(char*,uint64_t,...), input_text_parse_csv_array(char*,char**,uint64_t),
 *     input_text_parse_extended_uint64/_double(char*,...), input_text_parse_extended_bool(char*),
 *     input_text_parse_count_colons_in_field(const char*)
 *     fm_open_file(char*,fm_mode), fm_read_mem(fm_t*,void*,uint64_t), fm_eof(fm_t*), fm_close(fm_t*),
 *     fm_get_file_name(fm_t*), fmprintf(fm_t*,const char*,...), fm_bulk_read_file(char*,void*,uint64_t,uint64_t), gem_access(char*,fm_mode)
 *     mm_bulk_mload_file(char*), mm_bulk_mmalloc(uint64_t,bool), mm_bulk_malloc(uint64_t,bool),
 *     mm_bulk_mmalloc_temp(uint64_t), mm_bulk_free(mm_t*), mm_read_uint64(mm_t*), mm_read_mem(mm_t*,uint64_t),
 *     mm_skip_align_4KB(mm_t*), mm_get_mem(mm_t*), mm_get_tmp_folder(), mm_set_tmp_folder(char*),
 *     mm_get_mem_available_cached/_free/_total(), mm_malloc_nothrow, mm_realloc_nothrow, mm_free
 *     packed_integer_array_new/_store/_write, packed_integer_array_read_mem, packed_integer_array_load
 *     text_trace_compose_padded_text(text_trace_t*,uint64_t trim_left,uint64_t trim_right,mm_allocator_t*), text_trace_destroy()
 *   and reads these fields: sequence_t (tag, casava_tag.buffer, extra_tag.buffer)
 *   and text_trace_t (text, text_length, text_padded, text_padded_length,
 *   text_padded_left, text_padded_right; text_trace_padding sets them up).
 */

#include "gem3_test.h"
#include <unistd.h>
#include <zlib.h>
#include <bzlib.h>
#include "io/input_file_sliced.h"
#include "io/buffered_input_file.h"
#include "io/input_fasta.h"
#include "io/input_text.h"
#include "text/sequence.h"
#include "text/sequence_qualities_model.h"
#include "text/text_trace.h"
#include "text/dna_text.h"
#include "utils/packed_integer_array.h"

#define MAX_RECORDS 400
#define MAX_FIELD 512

/*
 * Temporary directory
 */
static char tmp_dir[512];
static void tmp_dir_create(void) {
  const char* tmp = getenv("TMPDIR");
  snprintf(tmp_dir,sizeof(tmp_dir),"%s/gem3-test-a-XXXXXX",(tmp && *tmp) ? tmp : "/tmp");
  GT_REQUIRE(mkdtemp(tmp_dir) != NULL);
}
static void tmp_path(char* const out,const size_t size,const char* const name) {
  snprintf(out,size,"%s/%s",tmp_dir,name);
}

/*
 * Records: what a file holds and what the parser must return
 */
typedef struct {
  char tag[MAX_FIELD];         // expected tag
  char read[MAX_FIELD];        // expected (normalized) read
  char qualities[MAX_FIELD];   // expected qualities ("" for FASTA)
  bool has_qualities;
  sequence_end_t end_info;
  char casava[MAX_FIELD];
  char extra[MAX_FIELD];
} record_t;
typedef struct {
  record_t records[MAX_RECORDS];
  uint64_t num_records;
  char* text;                  // the file content
  uint64_t text_length;
  uint64_t num_lines;          // lines in the file
} records_t;

static void text_append(records_t* const rs,const char* const s,const uint64_t n) {
  memcpy(rs->text+rs->text_length,s,n);
  rs->text_length += n;
  rs->text[rs->text_length] = '\0';
}
static void text_append_str(records_t* const rs,const char* const s) {
  text_append(rs,s,strlen(s));
}
static void text_eol(records_t* const rs,const bool crlf) {
  text_append_str(rs,crlf ? "\r\n" : "\n");
  ++rs->num_lines;
}
static char random_read_char(gt_rng_t* const rng,const bool last) {
  const uint64_t r = gt_rng_below(rng,100);
  if (r < 80) return "ACGT"[gt_rng_below(rng,4)];
  if (r < 88) return "acgt"[gt_rng_below(rng,4)];
  if (r < 92) return "Nn"[gt_rng_below(rng,2)];
  // The parser does not check a read's last character (known bug
  // known-bugs/A-fastq-last-character-unchecked.md), so only bases end a read
  if (last) return "ACGT"[gt_rng_below(rng,4)];
  if (r < 97) return "RYKMSWBDHVrykmswbdhvX"[gt_rng_below(rng,21)];
  return ".-*0"[gt_rng_below(rng,4)];
}
static char normalized(const char c) {
  return strchr("ACGTNacgtn",c) ? c : 'N';
}
/*
 * Generate @n records (FASTQ, FASTA or mixed) with random names, end
 * suffixes, CASAVA and extra fields, read lengths and characters
 */
static void records_generate(
    records_t* const rs,gt_rng_t* const rng,const uint64_t n,
    const int format /* 0 FASTQ, 1 FASTA, 2 mixed */,const bool crlf,const bool final_eol,
    const uint64_t max_length) {
  memset(rs,0,sizeof(records_t));
  rs->text = calloc(n*(4*MAX_FIELD+16)+16,1);
  uint64_t i, j;
  for (i=0;i<n;++i) {
    record_t* const r = rs->records + i;
    const bool fastq = (format==0) || (format==2 && gt_rng_chance(rng,1,2));
    r->has_qualities = fastq;
    // Header
    char name[64];
    snprintf(name,sizeof(name),"read%" PRIu64 "_%" PRIu64,i,gt_rng_below(rng,100000));
    char header[MAX_FIELD];
    strcpy(r->tag,name);
    r->end_info = single_end;
    const uint64_t style = gt_rng_below(rng,7);
    switch (style) {
      case 0: snprintf(header,sizeof(header),"%s",name); break;
      case 1: snprintf(header,sizeof(header),"%s/1",name); r->end_info = paired_end1; break;
      case 2: snprintf(header,sizeof(header),"%s/2",name); r->end_info = paired_end2; break;
      case 3: { // CASAVA 1.8
        const char end = "12"[gt_rng_below(rng,2)];
        snprintf(r->casava,sizeof(r->casava),"%c:N:0:ACGTAC",end);
        snprintf(header,sizeof(header),"%s %s",name,r->casava);
        r->end_info = (end=='1') ? paired_end1 : paired_end2;
        break;
      }
      case 4: // free text after a space or a tab
        snprintf(r->extra,sizeof(r->extra),"comment%" PRIu64,gt_rng_below(rng,1000));
        snprintf(header,sizeof(header),"%s%c%s",name,gt_rng_chance(rng,1,2) ? ' ' : '\t',r->extra);
        break;
      case 5: { // extra field carrying the end (SRA style): appended to the tag
        const char end = "123"[gt_rng_below(rng,3)];
        snprintf(header,sizeof(header),"%s HWI-ST382_0049:1:1:%" PRIu64 ":%" PRIu64 "/%c",name,gt_rng_below(rng,10000),gt_rng_below(rng,10000),end);
        const char* const field = strchr(header,' ')+1;
        snprintf(r->tag,sizeof(r->tag),"%s_%.*s",name,(int)(strlen(field)-2),field);
        r->end_info = (end=='1') ? paired_end1 : paired_end2;
        break;
      }
      default: // /3 is end 2
        snprintf(header,sizeof(header),"%s/3",name); r->end_info = paired_end2; break;
    }
    text_append_str(rs,fastq ? "@" : ">");
    text_append_str(rs,header);
    text_eol(rs,crlf);
    // Read
    const uint64_t length = gt_rng_range(rng,1,max_length);
    for (j=0;j<length;++j) {
      const char c = random_read_char(rng,j==length-1);
      text_append(rs,&c,1);
      r->read[j] = normalized(c);
    }
    r->read[length] = '\0';
    const bool last_line = (i==n-1);
    if (!fastq) {
      if (!last_line || final_eol) text_eol(rs,crlf);
      continue;
    }
    text_eol(rs,crlf);
    // Separator (sometimes repeating the name) and qualities
    text_append_str(rs,gt_rng_chance(rng,1,3) ? "+" : "+");
    if (gt_rng_chance(rng,1,3)) text_append_str(rs,name);
    text_eol(rs,crlf);
    for (j=0;j<length;++j) {
      const char q = (char)gt_rng_range(rng,33,126);
      r->qualities[j] = q;
    }
    r->qualities[length] = '\0';
    text_append(rs,r->qualities,length);
    if (!last_line || final_eol) text_eol(rs,crlf);
  }
  rs->num_records = n;
}
static void records_free(records_t* const rs) {
  free(rs->text);
}

/*
 * Writers for each compression
 */
static void write_plain(const char* const path,const records_t* const rs) {
  FILE* const f = fopen(path,"wb");
  GT_REQUIRE(f != NULL);
  fwrite(rs->text,1,rs->text_length,f);
  fclose(f);
}
static void write_gz(const char* const path,const records_t* const rs,const uint64_t members) {
  // members > 1: concatenated gzip members (as from `cat a.gz b.gz`)
  FILE* const f = fopen(path,"wb");
  GT_REQUIRE(f != NULL);
  fclose(f);
  uint64_t m, offset = 0;
  for (m=0;m<members;++m) {
    const uint64_t end = (m==members-1) ? rs->text_length : (rs->text_length*(m+1))/members;
    gzFile gz = gzopen(path,"ab6");
    GT_REQUIRE(gz != NULL);
    if (end > offset) gzwrite(gz,rs->text+offset,(unsigned)(end-offset));
    gzclose(gz);
    offset = end;
  }
}
static void write_bz2(const char* const path,const records_t* const rs) {
  FILE* const f = fopen(path,"wb");
  GT_REQUIRE(f != NULL);
  int bzerror;
  BZFILE* const bz = BZ2_bzWriteOpen(&bzerror,f,9,0,0);
  GT_REQUIRE(bzerror == BZ_OK);
  BZ2_bzWrite(&bzerror,bz,rs->text,(int)rs->text_length);
  BZ2_bzWriteClose(&bzerror,bz,0,NULL,NULL);
  fclose(f);
}

/*
 * Parse a whole input and compare with the records
 */
typedef enum { open_file, open_stream, open_gzip_stream, open_bzip_stream, open_pipe } open_mode_t;
static uint64_t parse_and_compare(
    const records_t* const rs,const char* const path,const open_mode_t mode,
    const uint64_t num_blocks,const uint64_t block_size,const uint64_t buffer_size,const char* const label) {
  mm_slab_t* const slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM);
  mm_allocator_t* const mm_allocator = mm_allocator_new(slab);
  input_file_sliced_t* sliced = NULL;
  FILE* stream = NULL;
  char command[700];
  switch (mode) {
    case open_file: sliced = input_file_sliced_open((char*)path,num_blocks,block_size); break;
    case open_stream: stream = fopen(path,"rb"); sliced = input_stream_sliced_open(stream,num_blocks,block_size); break;
    case open_gzip_stream: stream = fopen(path,"rb"); sliced = input_gzip_stream_sliced_open(stream,num_blocks,block_size); break;
    case open_bzip_stream: stream = fopen(path,"rb"); sliced = input_bzip_stream_sliced_open(stream,num_blocks,block_size); break;
    case open_pipe:
      snprintf(command,sizeof(command),"cat '%s'",path);
      sliced = input_file_sliced_popen(command,num_blocks,block_size);
      break;
  }
  buffered_input_file_t* const input = buffered_input_file_new(sliced,buffer_size);
  // File name as reported in error messages
  const char* const name = buffered_input_file_get_file_name(input);
  if (mode==open_file) {
    GT_CHECK_EQ_STR(name,path);
  } else if (mode==open_pipe) {
    GT_CHECK_EQ_STR(name,command);
  } else {
    GT_CHECK_EQ_STR(name,"stream");
  }
  sequence_t sequence;
  sequence_init(&sequence,bisulfite_disabled,mm_allocator);
  // The line counter (buffered_input_file_get_current_line_num) only feeds
  // error messages and starts uninitialized
  // (known-bugs/A-input-line-counter-uninitialized.md): not checked
  uint64_t i = 0, errors = 0;
  while (true) {
    if (buffered_input_file_eob(input)) {
      if (buffered_input_file_reload(input,0) == INPUT_STATUS_EOF) break;
    }
    const int status = input_fasta_parse_sequence(input,&sequence,false);
    if (status != INPUT_STATUS_OK) { ++errors; break; }
    if (i >= rs->num_records) { ++errors; break; }
    const record_t* const r = rs->records + i;
    // Tags, reads and qualities are counted strings (not NUL-terminated)
    const bool ok =
        string_get_length(&sequence.tag)==strlen(r->tag) &&
        strncmp(sequence_get_tag(&sequence),r->tag,strlen(r->tag))==0 &&
        sequence_get_length(&sequence)==strlen(r->read) &&
        strncmp(sequence_get_read(&sequence),r->read,strlen(r->read))==0 &&
        sequence_has_qualities(&sequence)==r->has_qualities &&
        (!r->has_qualities || strncmp(sequence_get_qualities(&sequence),r->qualities,strlen(r->qualities))==0) &&
        sequence_get_end_info(&sequence)==r->end_info &&
        sequence_has_casava_tag(&sequence)==(r->casava[0]!='\0') &&
        (r->casava[0]=='\0' || strncmp(sequence.casava_tag.buffer,r->casava,strlen(r->casava))==0) &&
        sequence_has_extra_tag(&sequence)==(r->extra[0]!='\0') &&
        (r->extra[0]=='\0' || strncmp(sequence.extra_tag.buffer,r->extra,strlen(r->extra))==0);
    if (!ok && errors++ < 3) {
      gt_note("%s: record %" PRIu64 " parsed as tag=%.*s end=%d length=%" PRIu64 "; expected tag=%s end=%d length=%zu read=%s",
          label,i,(int)string_get_length(&sequence.tag),sequence_get_tag(&sequence),(int)sequence_get_end_info(&sequence),
          sequence_get_length(&sequence),r->tag,(int)r->end_info,strlen(r->read),r->read);
    }
    ++i;
  }
  if (i != rs->num_records) {
    gt_note("%s: %" PRIu64 " records parsed, %" PRIu64 " written",label,i,rs->num_records);
    ++errors;
  }
  sequence_destroy(&sequence);
  buffered_input_file_close(input);
  input_file_sliced_close(sliced); // closes the stream too
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
  return errors;
}

/*
 * Cases
 */
// (num_blocks, block_size, buffer_size). A buffered read that stops at a
// block's end after a multiple of 8 lines drops the rest of the line that
// block ends in (known bug known-bugs/A-input-model-small-blocks-drop-lines.md);
// small blocks are therefore read either with a buffer larger than the
// input (one buffered read assembles lines across hundreds of blocks) or
// with forced line counts (buffered_input_forced_lines). The last model is
// gem-mapper's default shape.
static const uint64_t input_models[][3] = {
  {300,1000,1ull<<30}, {80,4096,1ull<<30}, {2,1024*1024,4*1024*1024},
};
#define NUM_INPUT_MODELS (sizeof(input_models)/sizeof(input_models[0]))

GT_TEST(fastq_plain_all_buffer_sizes) {
  tmp_dir_create();
  gt_rng_t rng = gt_rng_new(3);
  static records_t rs;
  char path[600];
  uint64_t format, crlf, m, errors = 0;
  for (format=0;format<3;++format) {
    for (crlf=0;crlf<2;++crlf) {
      records_generate(&rs,&rng,150,(int)format,crlf==1,format!=1,300);
      tmp_path(path,sizeof(path),"reads.fq");
      write_plain(path,&rs);
      for (m=0;m<NUM_INPUT_MODELS;++m) {
        char label[128];
        snprintf(label,sizeof(label),"format %" PRIu64 " crlf %" PRIu64 " model %" PRIu64,format,crlf,m);
        errors += parse_and_compare(&rs,path,open_file,input_models[m][0],input_models[m][1],input_models[m][2],label);
      }
      unlink(path);
      records_free(&rs);
    }
  }
  GT_CHECK_EQ_UINT(errors,0);
  rmdir(tmp_dir);
}
GT_TEST(fastq_compressed_streams_and_pipes) {
  tmp_dir_create();
  gt_rng_t rng = gt_rng_new(5);
  static records_t rs;
  records_generate(&rs,&rng,300,2,false,true,250);
  char plain[600], gz[600], gz3[600], bz[600];
  tmp_path(plain,sizeof(plain),"r.fq");
  tmp_path(gz,sizeof(gz),"r.fq.gz");
  tmp_path(gz3,sizeof(gz3),"r3.fq.gz");
  tmp_path(bz,sizeof(bz),"r.fq.bz2");
  write_plain(plain,&rs);
  write_gz(gz,&rs,1);
  write_gz(gz3,&rs,3);
  write_bz2(bz,&rs);
  uint64_t m, errors = 0;
  for (m=0;m<NUM_INPUT_MODELS;++m) {
    const uint64_t nb = input_models[m][0], bs = input_models[m][1], buf = input_models[m][2];
    // Compression detected from the file's magic bytes
    errors += parse_and_compare(&rs,gz,open_file,nb,bs,buf,"gzip file");
    errors += parse_and_compare(&rs,gz3,open_file,nb,bs,buf,"multi-member gzip file");
    errors += parse_and_compare(&rs,bz,open_file,nb,bs,buf,"bzip2 file");
    // Streams (stdin with and without -z / -j)
    errors += parse_and_compare(&rs,plain,open_stream,nb,bs,buf,"plain stream");
    errors += parse_and_compare(&rs,gz,open_gzip_stream,nb,bs,buf,"gzip stream");
    errors += parse_and_compare(&rs,bz,open_bzip_stream,nb,bs,buf,"bzip2 stream");
    // Pipe ("-i 'command |'")
    errors += parse_and_compare(&rs,plain,open_pipe,nb,bs,buf,"pipe");
  }
  GT_CHECK_EQ_UINT(errors,0);
  unlink(plain); unlink(gz); unlink(gz3); unlink(bz);
  records_free(&rs);
  rmdir(tmp_dir);
}
// A reload with forced lines takes exactly that many lines (how the second
// file of a split pair is kept in step with the first); lines straddle
// 64-byte blocks
GT_TEST(buffered_input_forced_lines) {
  tmp_dir_create();
  gt_rng_t rng = gt_rng_new(7);
  static records_t rs;
  records_generate(&rs,&rng,100,0,false,true,120);
  char path[600];
  tmp_path(path,sizeof(path),"r.fq");
  write_plain(path,&rs);
  mm_slab_t* const slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM);
  mm_allocator_t* const mm_allocator = mm_allocator_new(slab);
  input_file_sliced_t* const sliced = input_file_sliced_open(path,64,64);
  buffered_input_file_t* const input = buffered_input_file_new(sliced,1);
  sequence_t sequence;
  sequence_init(&sequence,bisulfite_disabled,mm_allocator);
  // A forced read holds every block it touches until the next reload, so a
  // chunk must fit in the 64 blocks (input that does not hangs: known bug
  // known-bugs/input-model-few-blocks-hang.md)
  static const uint64_t chunks[] = {8, 4, 40, 12};
  uint64_t c, parsed = 0, errors = 0;
  for (c=0;parsed<rs.num_records;++c) {
    const uint64_t chunk = chunks[c % (sizeof(chunks)/sizeof(chunks[0]))];
    GT_REQUIRE(buffered_input_file_reload(input,chunk) == INPUT_STATUS_OK);
    const uint64_t expected_lines = MIN(chunk,4*rs.num_records-4*parsed);
    GT_CHECK_EQ_UINT(buffered_input_file_get_num_lines(input),expected_lines);
    uint64_t k;
    for (k=0;k<expected_lines/4;++k) {
      GT_CHECK_EQ_INT(input_fasta_parse_sequence(input,&sequence,false),INPUT_STATUS_OK);
      const record_t* const r = rs.records + parsed;
      if (!(string_get_length(&sequence.tag)==strlen(r->tag) &&
            strncmp(sequence_get_tag(&sequence),r->tag,strlen(r->tag))==0 &&
            sequence_get_length(&sequence)==strlen(r->read) &&
            strncmp(sequence_get_read(&sequence),r->read,strlen(r->read))==0 &&
            strncmp(sequence_get_qualities(&sequence),r->qualities,strlen(r->qualities))==0)) ++errors;
      ++parsed;
    }
    GT_CHECK(buffered_input_file_eob(input));
  }
  GT_CHECK_EQ_UINT(errors,0);
  GT_CHECK_EQ_UINT(parsed,rs.num_records);
  GT_CHECK_EQ_INT(buffered_input_file_reload(input,0),INPUT_STATUS_EOF);
  sequence_destroy(&sequence);
  buffered_input_file_close(input);
  input_file_sliced_close(sliced);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
  unlink(path);
  records_free(&rs);
  rmdir(tmp_dir);
}
// Hand-written headers
GT_TEST(fastq_header_forms) {
  tmp_dir_create();
  static const char* const text =
      "@plain\nACGT\n+\nIIII\n"
      "@end1/1\nACGT\n+\nIIII\n"
      "@end2/2\nACGT\n+\nIIII\n"
      "@x/4\nACGT\n+\nIIII\n"                                   // not an end suffix
      "@EAS139:136:FC706VJ:2:5:1000:12850 2:Y:18:ATCACG\nACGT\n+\nIIII\n" // CASAVA 1.8
      "@SRR384920.1 HWI-ST382_0049:1:1:1217:1879/1\nACGT\n+\nIIII\n"     // end in an extra field
      "@t\tsome comment\nACGT\n+\nIIII\n"
      ">fasta/2 desc\nAC.T\n";
  static const struct { const char* tag; int end; const char* casava; const char* extra; const char* read; } expected[] = {
    { "plain", single_end, "", "", "ACGT" },
    { "end1", paired_end1, "", "", "ACGT" },
    { "end2", paired_end2, "", "", "ACGT" },
    { "x/4", single_end, "", "", "ACGT" },
    { "EAS139:136:FC706VJ:2:5:1000:12850", paired_end2, "2:Y:18:ATCACG", "", "ACGT" },
    { "SRR384920.1_HWI-ST382_0049:1:1:1217:1879", paired_end1, "", "", "ACGT" },
    { "t", single_end, "", "comment", "ACGT" },          // the last field is kept
    { "fasta", paired_end2, "", "desc", "ACNT" },
  };
  char path[600];
  tmp_path(path,sizeof(path),"h.fq");
  FILE* const f = fopen(path,"wb");
  fputs(text,f);
  fclose(f);
  mm_slab_t* const slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM);
  mm_allocator_t* const mm_allocator = mm_allocator_new(slab);
  input_file_sliced_t* const sliced = input_file_sliced_open(path,2,4096);
  buffered_input_file_t* const input = buffered_input_file_new(sliced,4096);
  sequence_t sequence;
  sequence_init(&sequence,bisulfite_disabled,mm_allocator);
  uint64_t i;
  for (i=0;i<sizeof(expected)/sizeof(expected[0]);++i) {
    GT_REQUIRE(input_fasta_parse_sequence(input,&sequence,true) == INPUT_STATUS_OK);
    GT_CHECK_MSG(string_get_length(&sequence.tag)==strlen(expected[i].tag) &&
        strncmp(sequence_get_tag(&sequence),expected[i].tag,strlen(expected[i].tag))==0,
        "tag %.*s expected %s",(int)string_get_length(&sequence.tag),sequence_get_tag(&sequence),expected[i].tag);
    GT_CHECK_EQ_INT(sequence_get_end_info(&sequence),expected[i].end);
    GT_CHECK_EQ_INT(sequence_has_casava_tag(&sequence),expected[i].casava[0]!='\0');
    if (expected[i].casava[0]) GT_CHECK(strncmp(sequence.casava_tag.buffer,expected[i].casava,strlen(expected[i].casava))==0);
    GT_CHECK_EQ_INT(sequence_has_extra_tag(&sequence),expected[i].extra[0]!='\0');
    if (expected[i].extra[0]) GT_CHECK(strncmp(sequence.extra_tag.buffer,expected[i].extra,strlen(expected[i].extra))==0);
    GT_CHECK_EQ_UINT(sequence_get_length(&sequence),strlen(expected[i].read));
    GT_CHECK(strncmp(sequence_get_read(&sequence),expected[i].read,strlen(expected[i].read))==0);
  }
  GT_CHECK_EQ_INT(input_fasta_parse_sequence(input,&sequence,true),INPUT_STATUS_EOF);
  sequence_destroy(&sequence);
  buffered_input_file_close(input);
  input_file_sliced_close(sliced);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
  unlink(path);
  rmdir(tmp_dir);
}

/*
 * Command-line value parsers (input_text.c)
 */
GT_TEST(input_text_parsers) {
  // Integers: decimal with sign, hexadecimal; stop at the first non-digit
  static const struct { const char* text; int ret; int64_t value; uint64_t consumed; } ints[] = {
    { "0", 0, 0, 1 }, { "42", 0, 42, 2 }, { "+7", 0, 7, 2 }, { "-15x", 0, -15, 3 },
    { "0x1F", 0, 31, 4 }, { "0Xff", 0, 255, 4 }, { "123456789012", 0, 123456789012ll, 12 },
    { "abc", -1, 0, 0 }, { "-", -1, 0, 1 }, { "0xg", -1, 0, 2 },
  };
  uint64_t i;
  for (i=0;i<sizeof(ints)/sizeof(ints[0]);++i) {
    const char* p = ints[i].text;
    int64_t value = 0;
    const int ret = input_text_parse_integer(&p,&value);
    GT_CHECK_MSG(ret==ints[i].ret,"parse_integer(%s) returned %d",ints[i].text,ret);
    if (ret==0) {
      GT_CHECK_MSG(value==ints[i].value,"parse_integer(%s) = %" PRId64,ints[i].text,value);
      GT_CHECK_MSG((uint64_t)(p-ints[i].text)==ints[i].consumed,"parse_integer(%s) consumed %d",ints[i].text,(int)(p-ints[i].text));
    }
  }
  // Doubles: [+-]?[0-9]*(.[0-9]+)?([eE][+-]?[0-9]+)?
  static const struct { const char* text; int ret; double value; } dbls[] = {
    { "0.5", 0, 0.5 }, { "12", 0, 12.0 }, { "-3.25", 0, -3.25 }, { ".125", 0, 0.125 },
    { "2e3", 0, 2000.0 }, { "1.5E-2", 0, 0.015 }, { "+4.0e+1", 0, 40.0 },
    { "1.", -1, 0 }, { "3e", -1, 0 },
  };
  for (i=0;i<sizeof(dbls)/sizeof(dbls[0]);++i) {
    const char* p = dbls[i].text;
    double value = 0;
    const int ret = input_text_parse_double(&p,&value);
    GT_CHECK_MSG(ret==dbls[i].ret,"parse_double(%s) returned %d",dbls[i].text,ret);
    if (ret==0) GT_CHECK_MSG(fabs(value-dbls[i].value) <= 1e-12*fabs(dbls[i].value),"parse_double(%s) = %g",dbls[i].text,value);
  }
  // Sizes: B, K, M, G (1024-based, optional trailing B)
  static const struct { const char* text; int ret; uint64_t value; } sizes[] = {
    { "100", 0, 100 }, { "100B", 0, 100 }, { "4K", 0, 4096 }, { "4k", 0, 4096 }, { "64KB", 0, 65536 },
    { "32M", 0, 32ull<<20 }, { "2mb", 0, 2ull<<20 }, { "1G", 0, 1ull<<30 }, { "3gB", 0, 3ull<<30 },
    { "4X", 1, 0 }, { "4KX", 1, 0 }, { "10BB", 1, 0 },
  };
  for (i=0;i<sizeof(sizes)/sizeof(sizes[0]);++i) {
    char buffer[32];
    strcpy(buffer,sizes[i].text);
    uint64_t value = 0;
    const int ret = input_text_parse_size(buffer,&value);
    GT_CHECK_MSG(ret==sizes[i].ret,"parse_size(%s) returned %d",sizes[i].text,ret);
    if (ret==0) GT_CHECK_MSG(value==sizes[i].value,"parse_size(%s) = %" PRIu64,sizes[i].text,value);
  }
  // CSV arguments: at most the number asked for
  char csv[64];
  strcpy(csv,"32M,c,4M,extra");
  char *a = NULL, *b = NULL, *c = NULL;
  GT_CHECK_EQ_INT(input_text_parse_csv_arguments(csv,3,&a,&b,&c),3);
  GT_CHECK_EQ_STR(a,"32M"); GT_CHECK_EQ_STR(b,"c"); GT_CHECK_EQ_STR(c,"4M");
  strcpy(csv,"fixed,20");
  a = b = c = NULL;
  GT_CHECK_EQ_INT(input_text_parse_csv_arguments(csv,3,&a,&b,&c),2);
  GT_CHECK_EQ_STR(a,"fixed"); GT_CHECK_EQ_STR(b,"20"); GT_CHECK(c==NULL);
  char* array[4];
  strcpy(csv,"a,bb,ccc,d,e");
  GT_CHECK_EQ_INT(input_text_parse_csv_array(csv,array,4),4);
  GT_CHECK_EQ_STR(array[0],"a"); GT_CHECK_EQ_STR(array[3],"d");
  // Extended values
  static const struct { const char* text; uint64_t value; } ext[] = {
    { "all", UINT64_MAX }, { "INF", UINT64_MAX }, { "unlimited", UINT64_MAX }, { "none", 0 },
    { "zero", 0 }, { "null", 0 }, { "17", 17 }, { "0x10", 16 },
  };
  for (i=0;i<sizeof(ext)/sizeof(ext[0]);++i) {
    char buffer[32];
    strcpy(buffer,ext[i].text);
    uint64_t value = 1;
    GT_CHECK_EQ_INT(input_text_parse_extended_uint64(buffer,&value),0);
    GT_CHECK_MSG(value==ext[i].value,"parse_extended_uint64(%s) = %" PRIu64,ext[i].text,value);
  }
  char buffer[32];
  double d = 0;
  strcpy(buffer,"all"); GT_CHECK_EQ_INT(input_text_parse_extended_double(buffer,&d),0); GT_CHECK(d==(double)INT64_MAX);
  strcpy(buffer,"-inf"); GT_CHECK_EQ_INT(input_text_parse_extended_double(buffer,&d),0); GT_CHECK(d==(double)INT64_MIN);
  strcpy(buffer,"none"); GT_CHECK_EQ_INT(input_text_parse_extended_double(buffer,&d),0); GT_CHECK(d==0.0);
  strcpy(buffer,"0.08"); GT_CHECK_EQ_INT(input_text_parse_extended_double(buffer,&d),0); GT_CHECK(fabs(d-0.08)<1e-15);
  GT_CHECK(input_text_parse_extended_bool(NULL));
  static const char* const yes[] = {"true","yes","allow","enabled"};
  static const char* const no[] = {"false","no","disabled","0","TRUE"};
  for (i=0;i<4;++i) { strcpy(buffer,yes[i]); GT_CHECK_MSG(input_text_parse_extended_bool(buffer),"bool(%s)",yes[i]); }
  for (i=0;i<5;++i) { strcpy(buffer,no[i]); GT_CHECK_MSG(!input_text_parse_extended_bool(buffer),"bool(%s)",no[i]); }
  // Colons in the first field (up to a space, tab or end of line)
  GT_CHECK_EQ_UINT(input_text_parse_count_colons_in_field("1:N:0:ACGT rest\n"),3);
  GT_CHECK_EQ_UINT(input_text_parse_count_colons_in_field("a:b\tc:d:e\n"),1);
  GT_CHECK_EQ_UINT(input_text_parse_count_colons_in_field("none\n"),0);
}

/*
 * Quality model: phred values from the ASCII offset
 */
GT_TEST(qualities_model) {
  char q33[] = "!+5?I~";
  char q64[] = "@JT^h~";
  uint8_t mask[8];
  sequence_qualities_model_process(sequence_qualities_offset_33,q33,6,mask);
  static const uint8_t expected33[] = {0,10,20,30,40,93};
  GT_CHECK_EQ_MEM(mask,expected33,6);
  sequence_qualities_model_process(sequence_qualities_offset_64,q64,6,mask);
  static const uint8_t expected64[] = {0,10,20,30,40,62};
  GT_CHECK_EQ_MEM(mask,expected64,6);
  // Only key_length values are written
  memset(mask,0xAA,sizeof(mask));
  sequence_qualities_model_process(sequence_qualities_offset_33,q33,3,mask);
  GT_CHECK(mask[2]==20 && mask[3]==0xAA);
}

/*
 * Padded candidate text (N on both sides)
 */
GT_TEST(text_trace_padding) {
  mm_slab_t* const slab = mm_slab_new_(BUFFER_SIZE_8M,BUFFER_SIZE_8M,MM_UNLIMITED_MEM);
  mm_allocator_t* const mm_allocator = mm_allocator_new(slab);
  uint8_t text[] = {0,1,2,3,3,2};
  text_trace_t trace;
  memset(&trace,0,sizeof(trace));
  trace.text = text;
  trace.text_length = 6;
  trace.text_padded = text;
  trace.text_padded_length = 6;
  text_trace_compose_padded_text(&trace,0,0,mm_allocator);
  GT_CHECK_EQ_UINT(trace.text_padded_length,6);
  GT_CHECK_EQ_MEM(trace.text_padded,text,6);
  GT_CHECK(trace.text_padded_left==0 && trace.text_padded_right==0);
  text_trace_compose_padded_text(&trace,2,3,mm_allocator);
  static const uint8_t expected[] = {4,4,0,1,2,3,3,2,4,4,4};
  GT_CHECK_EQ_UINT(trace.text_padded_length,11);
  GT_CHECK_EQ_MEM(trace.text_padded,expected,11);
  GT_CHECK(trace.text_padded_left==2 && trace.text_padded_right==3);
  text_trace_destroy(&trace,mm_allocator); // frees the padded copy only
  GT_CHECK(text[0]==0 && text[5]==2);
  mm_allocator_delete(mm_allocator);
  mm_slab_delete(slab);
}

/*
 * File manager and bulk memory (the index loader reads the .gem file with
 * mm_bulk_mload_file and walks it with mm_read_*)
 */
GT_TEST(file_and_memory_managers) {
  tmp_dir_create();
  char path[600], gz[600], bz[600];
  tmp_path(path,sizeof(path),"data.bin");
  tmp_path(gz,sizeof(gz),"data.gz");
  tmp_path(bz,sizeof(bz),"data.bz2");
  // 3 words, 4 KB alignment padding, 2 more words, then bytes
  const uint64_t size = 3*8 + (4096-24) + 2*8 + 100;
  uint8_t* const data = calloc(size,1);
  uint64_t words[5] = {5012, 7, UINT64_MAX-3, 0x0123456789ABCDEFull, 42};
  memcpy(data,words,24);
  memcpy(data+4096,words+3,16);
  uint64_t i;
  for (i=0;i<100;++i) data[4096+16+i] = (uint8_t)(i*7);
  FILE* f = fopen(path,"wb"); fwrite(data,1,size,f); fclose(f);
  // mm_bulk_mload_file: the whole file in memory; reads advance a cursor
  mm_t* const mm = mm_bulk_mload_file(path);
  GT_CHECK_EQ_UINT(mm_read_uint64(mm),5012);
  GT_CHECK_EQ_UINT(mm_read_uint64(mm),7);
  GT_CHECK_EQ_UINT(mm_read_uint64(mm),UINT64_MAX-3);
  mm_skip_align_4KB(mm);
  GT_CHECK_EQ_UINT(mm_read_uint64(mm),0x0123456789ABCDEFull);
  GT_CHECK_EQ_UINT(mm_read_uint64(mm),42);
  const uint8_t* const bytes = mm_read_mem(mm,100);
  GT_CHECK_EQ_MEM(bytes,data+4096+16,100);
  GT_CHECK((uint8_t*)mm_get_mem(mm) == bytes+100);
  mm_bulk_free(mm);
  // fm_bulk_read_file: from an offset
  uint8_t* const copy = malloc(size);
  fm_bulk_read_file(path,copy,4096,16+100);
  GT_CHECK_EQ_MEM(copy,data+4096,116);
  fm_bulk_read_file(path,copy,0,0);
  GT_CHECK_EQ_MEM(copy,data,size);
  // fm_open_file reads plain, gzip and bzip2 files alike
  FILE* gzf = fopen(gz,"wb"); fclose(gzf);
  gzFile g = gzopen(gz,"wb"); gzwrite(g,data,(unsigned)size); gzclose(g);
  FILE* bf = fopen(bz,"wb");
  int bzerror;
  BZFILE* bzw = BZ2_bzWriteOpen(&bzerror,bf,9,0,0);
  BZ2_bzWrite(&bzerror,bzw,data,(int)size);
  BZ2_bzWriteClose(&bzerror,bzw,0,NULL,NULL);
  fclose(bf);
  const char* const paths[] = {path,gz,bz};
  for (i=0;i<3;++i) {
    fm_t* const fm = fm_open_file((char*)paths[i],FM_READ);
    GT_CHECK_EQ_STR(fm_get_file_name(fm),paths[i]);
    memset(copy,0,size);
    uint64_t got = 0, chunk;
    while (!fm_eof(fm) && got < size + 10) {
      chunk = fm_read_mem(fm,copy+got,MIN(1000,size+10-got));
      got += chunk;
      if (chunk==0) break;
    }
    GT_CHECK_MSG(got==size,"%s: read %" PRIu64 " of %" PRIu64 " bytes",paths[i],got,size);
    GT_CHECK_EQ_MEM(copy,data,size);
    GT_CHECK(fm_eof(fm));
    fm_close(fm);
  }
  // fm_printf through a file manager (plain output)
  char out[600];
  tmp_path(out,sizeof(out),"out.txt");
  fm_t* const fo = fm_open_file(out,FM_WRITE);
  fmprintf(fo,"%s:%d:%c\n","chr1",12345,'+');
  fm_close(fo);
  f = fopen(out,"rb");
  char line[64] = {0};
  GT_CHECK(fgets(line,sizeof(line),f) != NULL);
  fclose(f);
  GT_CHECK_EQ_STR(line,"chr1:12345:+\n");
  // gem_access
  GT_CHECK(gem_access(path,FM_READ));
  GT_CHECK(!gem_access("/nonexistent/gem3/file",FM_READ));
  unlink(out); unlink(path); unlink(gz); unlink(bz);
  free(copy);
  free(data);
  // Bulk memory: writable (mmap'ed, or from malloc in GEM_DEBUG builds);
  // mm_bulk_malloc with init_mem is zero-filled
  mm_t* const anon = mm_bulk_mmalloc(1<<20,false);
  uint8_t* const am = mm_get_base_mem(anon);
  am[0] = 3; am[12345] = 9; am[(1<<20)-1] = 7;
  GT_CHECK(am[0]==3 && am[12345]==9 && am[(1<<20)-1]==7);
  mm_bulk_free(anon);
  mm_t* const heap = mm_bulk_malloc(1<<16,true);
  uint8_t* const hm = mm_get_base_mem(heap);
  GT_CHECK(hm[0]==0 && hm[(1<<16)-1]==0);
  mm_bulk_free(heap);
  // Temporary (disk-backed) memory goes to the temporary folder
  char* const saved = mm_get_tmp_folder();
  char folder[600];
  snprintf(folder,sizeof(folder),"%s/",tmp_dir);
  mm_set_tmp_folder(folder);
  GT_CHECK_EQ_STR(mm_get_tmp_folder(),folder);
  FILE* const log_stream = gem_log_get_stream();
  FILE* const devnull = fopen("/dev/null","w");
  if (devnull) gem_log_set_stream(devnull);
  mm_t* const temp = mm_bulk_mmalloc_temp(1<<20);
  gem_log_set_stream(log_stream);
  if (devnull) fclose(devnull);
  uint8_t* const tm = mm_get_base_mem(temp);
  tm[0] = 1; tm[(1<<20)-1] = 2;
  GT_CHECK(tm[0]==1 && tm[(1<<20)-1]==2);
  mm_bulk_free(temp);
  mm_set_tmp_folder(saved);
  // Small allocations
  uint64_t* const v = mm_malloc_nothrow(16,sizeof(uint64_t),true,0);
  GT_REQUIRE(v != NULL);
  for (i=0;i<16;++i) GT_CHECK(v[i]==0);
  v[15] = 77;
  uint64_t* const w = mm_realloc_nothrow(v,64*sizeof(uint64_t));
  GT_REQUIRE(w != NULL);
  GT_CHECK(w[15]==77);
  mm_free(w);
  // Memory figures agree with /proc/meminfo (free + cached)
  const int64_t total = mm_get_mem_available_total();
  const int64_t free_mem = mm_get_mem_available_free();
  const int64_t cached = mm_get_mem_available_cached();
  GT_CHECK(free_mem > 0 && cached >= 0);
  GT_CHECK(total > 0 && total <= mm_get_mem_total());
  rmdir(tmp_dir);
}

/*
 * Packed integer arrays (sampled SA and sampled RL storage)
 */
GT_TEST(packed_integer_array_roundtrip) {
  tmp_dir_create();
  char path[600];
  tmp_path(path,sizeof(path),"packed.bin");
  gt_rng_t rng = gt_rng_new(11);
  static const uint64_t bits[] = {1, 5, 19, 31, 32, 33, 47, 63, 64};
  uint64_t b, errors = 0;
  for (b=0;b<sizeof(bits)/sizeof(bits[0]);++b) {
    const uint64_t n = gt_rng_range(&rng,1,3000);
    uint64_t* const values = malloc(n*sizeof(uint64_t));
    packed_integer_array_t* const array = packed_integer_array_new(n,bits[b]);
    uint64_t i;
    for (i=0;i<n;++i) {
      values[i] = (bits[b]==64) ? gt_rng_u64(&rng) : gt_rng_below(&rng,1ull<<bits[b]);
      packed_integer_array_store(array,i,values[i]);
    }
    fm_t* const fm = fm_open_file(path,FM_WRITE);
    packed_integer_array_write(fm,array);
    fm_close(fm);
    packed_integer_array_delete(array);
    mm_t* const mm = mm_bulk_mload_file(path);
    packed_integer_array_t* const loaded = packed_integer_array_read_mem(mm);
    GT_CHECK(packed_integer_array_get_size(loaded) > 0);
    for (i=0;i<n;++i) {
      packed_integer_array_prefetch(loaded,i);
      if (packed_integer_array_load(loaded,i) != values[i]) ++errors;
    }
    packed_integer_array_delete(loaded);
    mm_bulk_free(mm);
    free(values);
  }
  GT_CHECK_EQ_UINT(errors,0);
  unlink(path);
  rmdir(tmp_dir);
}

GT_CASES(
  GT_CASE(fastq_plain_all_buffer_sizes),
  GT_CASE(fastq_compressed_streams_and_pipes),
  GT_CASE(buffered_input_forced_lines),
  GT_CASE(fastq_header_forms),
  GT_CASE(input_text_parsers),
  GT_CASE(qualities_model),
  GT_CASE(text_trace_padding),
  GT_CASE(file_and_memory_managers),
  GT_CASE(packed_integer_array_roundtrip)
)
