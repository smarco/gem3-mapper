/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
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
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION:
 *   Output module provides SAM output generation
 */

#include "io/output_sam.h"
#include "io/output_map.h"
#include "matches/matches.h"
#include "matches/matches_cigar.h"
#include "matches/paired_map.h"
#include "system/gruntime.h"

/*
 * Constants
 */
#define OUTPUT_SAM_FORMAT_VERSION "1.4"

// SAM FLAGS
#define SAM_FLAG_MULTIPLE_SEGMENTS 0x1
#define SAM_FLAG_PROPERLY_ALIGNED 0x2
#define SAM_FLAG_UNMAPPED 0x4
#define SAM_FLAG_NEXT_UNMAPPED 0x8
#define SAM_FLAG_REVERSE_COMPLEMENT 0x10
#define SAM_FLAG_NEXT_REVERSE_COMPLEMENT 0x20
#define SAM_FLAG_FIRST_SEGMENT 0x40
#define SAM_FLAG_LAST_SEGMENT 0x80
#define SAM_FLAG_SECONDARY_ALIGNMENT 0x100
#define SAM_FLAG_NOT_PASSING_QC 0x200
#define SAM_FLAG_PCR_OR_OPTICAL_DUPLICATE 0x400
#define SAM_FLAG_SUPPLEMENTARY_ALIGNMENT 0x800

/*
 * Setup
 */
void output_sam_parameters_set_defaults(output_sam_parameters_t* const sam_parameters) {
  /* Header & RG  */
  sam_parameters->read_group_header = NULL;
  sam_parameters->read_group_id = NULL;
  /* Read & Qualities */
  sam_parameters->omit_secondary_read__qualities = true;
  /* CIGAR */
  sam_parameters->print_mismatches = false;
  /* XA */
  sam_parameters->compact_xa = true;
  /* Bisulfite */
  sam_parameters->bisulfite_output = false;
  /* GEM compatibility */
  sam_parameters->print_gem_fields = false;
  /* Benchmark mode */
  sam_parameters->benchmark_mode = false;
}
/*
 * SAM Headers
 *
 * @HD  VN:1.0  SO:unsorted
 * @HD  VN:1.3
 *
 * @SQ  SN:chr10   LN:135534747  AS:hg19_ncbi37  SP:human
 * @SQ  SN:chr10   LN:135534747
 * @SQ  SN:chr11   LN:135006516
 *
 * @RG  ID:NOID   PG:GEM   SM:NOSM
 * @RG  ID:0      PG:GEM   PL:ILLUMINA  SM:0
 *
 */
void output_sam_print_header(
		output_file_t* const output_file,
		archive_t* const archive,
		output_sam_parameters_t* const sam_parameters,
		int argc,
		char** argv,
		char* const gem_version) {
	/*
	 * TODO
	 *   Other options like BWA/Bowtie2
	 */
	// Print all @HD lines (Header line)
	//   @HD  VN:1.0  SO:unsorted
	ofprintf(output_file,"@HD\tVN:"OUTPUT_SAM_FORMAT_VERSION"\tSO:unsorted\n");
	// Print all @SQ lines (Reference sequence dictionary)
	//   @SQ  SN:chr10 LN:135534747 AS:hg19_ncbi37  SP:human
	locator_t* locator = archive->locator;
	const uint64_t num_intervals = locator->num_intervals;
	const locator_interval_t* const intervals = locator->intervals;
	uint64_t i = 0;
	while (i<num_intervals) {
		const int64_t tag_id = intervals[i].tag_id;
		const strand_t strand = intervals[i].strand;
		const bs_strand_t bs_strand = intervals[i].bs_strand;
		// Skip Reverse or G2A contigs
		if (strand==Reverse || bs_strand==bs_strand_G2A) {
			++i; continue;
		}
		// Calculate the length of the sequence (compose by several intervals)
		while (i+1<num_intervals && intervals[i+1].tag_id==tag_id &&
				intervals[i+1].strand==strand && intervals[i+1].bs_strand==bs_strand) ++i; // Go to the last
		const uint64_t total_length = intervals[i].sequence_offset+intervals[i].sequence_length;
		const char* const tag = locator_interval_get_tag(locator,intervals+i);
		uint64_t tag_len = gem_strlen(tag);
		// Print SQ
		ofprintf(output_file,"@SQ\tSN:%"PRIs"\tLN:%"PRIu64"\n",tag_len,tag,total_length);
		// Next
		++i;
	}
	// Print @RG line (Read group)
	if (sam_parameters->read_group_header) ofprintf(output_file,"%s\n",sam_parameters->read_group_header);
	if(sam_parameters->benchmark_mode) {
		ofprintf(output_file,"@CO\tgemBS benchmark mode\n");
	} else {

		// Print all @PG lines (Program)
		//   @PG  ID:GEM PN:gem-2-sam  VN:3.0.0 CL:gem-mapper -I /home/u/hsapiens_v37 -i /h/uu.fastq -e 0.08
		ofprintf(output_file,"@PG\tID:GEM\tPN:gem-mapper\tVN:%s",gem_version);
		size_t t_len = 256;
		char *q = mm_malloc(t_len);
		for (i=0;i<argc;++i) {
			// Escape tabs etc. in command line to avoid confusing SAM parsers
			char* p = argv[i];
			size_t l = strlen(p) * 2 + 1;
			if(l>t_len) {
				t_len = l;
				q = mm_realloc(q , t_len);
			}
			char* r = q;
			while(*p) {
				switch(*p) {
				case '\t':
					*r++ = '\\';
					*r++ = 't';
					break;
				case '\n':
					*r++ = '\\';
					*r++ = 'n';
					break;
				case '\r':
					*r++ = '\\';
					*r++ = 'r';
					break;
				default:
					*r++ = *p;
					break;
				}
				p++;
			}
			*r=0;
			ofprintf(output_file,"%s%s",(i==0)?"\tCL:":" ",q);
		}
		mm_free(q);
		ofprintf(output_file,"\n");
		// Print all @CO lines (Comments)
		//   @CO  TM:Fri, 30 Nov 2012 14:14:13 CET        WD:/home/user/CMD      HN:cn38   UN:user
		// @CO Print Current Date
		time_t current_time=time(0);
		struct tm local_time;
		localtime_r(&current_time,&local_time);
		// @CO Print Current year
		ofprintf(output_file,"@CO\tTM:%4d/%d/%d %02d:%02d:%02d CET",
				1900+local_time.tm_year,local_time.tm_mon+1,local_time.tm_mday,
			  local_time.tm_hour,local_time.tm_min,local_time.tm_sec);
	  char* cwd = system_get_cwd();
	  ofprintf(output_file,"\tWD:%s",cwd);
	  free(cwd);
	  // @CO Print Host
	  char* hostname = system_get_hostname();
	  ofprintf(output_file,"\tHN:%s",hostname);
	  free(hostname);
	  // @CO Print User
	  const char* const user_name = system_get_user_name();
	  ofprintf(output_file,"\tUN:%s\n",user_name);
  }
}
/*
 * SAM CIGAR
 */
void output_sam_print_cigar_element(
    buffered_output_file_t* const buffered_output_file,
    const cigar_element_t* const cigar_element) {
  switch (cigar_element->type) {
    case cigar_match:
      bofprintf_uint64(buffered_output_file,(uint32_t)cigar_element->length);
      bofprintf_char(buffered_output_file,'M');
      break;
    case cigar_mismatch:
      bofprintf_char(buffered_output_file,'X');
      break;
    case cigar_ins:
      bofprintf_int64(buffered_output_file,cigar_element->length);
      bofprintf_char(buffered_output_file,'D');
      break;
    case cigar_del:
      if (cigar_element->attributes == cigar_attr_trim) {
        bofprintf_int64(buffered_output_file,cigar_element->length);
        bofprintf_char(buffered_output_file,'S');
      } else {
        bofprintf_int64(buffered_output_file,cigar_element->length);
        bofprintf_char(buffered_output_file,'I');
      }
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
void output_sam_print_reduced_cigar_element(
    buffered_output_file_t* const buffered_output_file,
    const cigar_element_t* const cigar_element,
    uint64_t* const matching) {
  if (cigar_element->type==cigar_match) {
    *matching += (uint32_t)cigar_element->length;
  } else if (cigar_element->type==cigar_mismatch) {
    ++(*matching);
  } else {
    if (*matching > 0) {
      bofprintf_uint64(buffered_output_file,*matching);
      bofprintf_char(buffered_output_file,'M');
      *matching = 0;
    }
    output_sam_print_cigar_element(buffered_output_file,cigar_element);
  }
}
void output_sam_print_md_cigar_element(
    buffered_output_file_t* const buffered_output_file,
    const uint8_t* const text,
    const uint64_t text_length,
    uint64_t* const text_pos,
    const strand_t match_strand,
    const cigar_element_t* const cigar_element,
    uint64_t* const matching_length) {
  // Handle Matches/Mismatches
  switch (cigar_element->type) {
    case cigar_match:
      *matching_length += cigar_element->length;
      *text_pos += cigar_element->length;
      break;
    case cigar_ins:
      if (*matching_length > 0) {
        bofprintf_uint64(buffered_output_file,*matching_length);
        *matching_length = 0;
      }
      // Print Indel
      bofprintf_char(buffered_output_file,'^');
      const uint64_t length = cigar_element->length;
      int64_t i; // Print inserted text (text from the reference)
      if (match_strand==Forward) {
        for (i=0;i<length;++i) {
          bofprintf_char(buffered_output_file,dna_decode(text[*text_pos+i]));
        }
      } else {
        for (i=1;i<=length;++i) {
          bofprintf_char(buffered_output_file,dna_complement(dna_decode(text[text_length-(*text_pos+i)])));
        }
      }
      *text_pos += length;
      break;
    case cigar_mismatch:
      bofprintf_uint64(buffered_output_file,*matching_length);
      *matching_length = 0;
      if (match_strand==Forward) {
        bofprintf_char(buffered_output_file,dna_decode(cigar_element->mismatch));
      } else {
        bofprintf_char(buffered_output_file,dna_complement(dna_decode(cigar_element->mismatch)));
      }
      ++(*text_pos);
      break;
    case cigar_del:
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
void output_sam_parse_read_group_header(
    char* const read_group_buffer,
    output_sam_parameters_t* const sam_parameters) {
  // Check TAG RG
  gem_cond_fatal_error_msg(strncmp(read_group_buffer,"@RG",3),"Option '-r|--sam-read-group-header' must start with '@RG'");
  // Translate tab and new line escape sequences
  char *p = read_group_buffer;
  char *q = mm_malloc(strlen(read_group_buffer) + 1);
  char *r = q;
  while (*p) {
    if (*p == '\\') {
      p++;
      switch (*p) {
      case 't':
        *q++ = '\t';
        break;
      case 'n':
        *q++ = '\n';
        break;
      default:
        *q++ = *p;
        break;
      }
    } else {
      *q++ = *p;
    }
    p++;
  }
  // Strip of trailing new line if present
  if (q != r && *(q - 1) == '\n') {
    q--;
  }
  *q = 0;
  p = strstr(r,"\tID:");
  gem_cond_fatal_error_msg(p == NULL,"Option '-r|--sam-read-group-header' must contain ID field");
  p += 4;
  q = p;
  while (*q && *q != '\n' && *q != '\t') {
    q++;
  }
  gem_cond_fatal_error_msg(q - p > 255,"Option '-r|--sam-read-group-header' has ID longer than 255 characters");
  sam_parameters->read_group_header = r;
  sam_parameters->read_group_id = mm_malloc(sizeof(string_t));
  string_init(sam_parameters->read_group_id, q - p + 1,NULL);
  string_set_buffer(sam_parameters->read_group_id, p, q - p);
}
/*
 * Optional fields
 */
//  AM  i  The smallest template-independent mapping quality of segments in the rest // FIXME
//  AS  i  Alignment score generated by aligner
void output_sam_print_opt_field_tag_AS(
    buffered_output_file_t* const buffered_output_file,
    const match_alignment_model_t match_alignment_model,
    const match_trace_t* const match_trace) {
  // Switch alignment model
  buffered_output_file_reserve(buffered_output_file,6+INT_MAX_LENGTH);
  switch (match_alignment_model) {
    case match_alignment_model_none: break;
    case match_alignment_model_hamming:
    case match_alignment_model_levenshtein:
      bofprintf_string_literal(buffered_output_file,"\tAS:i:");
      bofprintf_uint64(buffered_output_file,match_trace->edit_distance);
      break;
    case match_alignment_model_gap_affine:
      bofprintf_string_literal(buffered_output_file,"\tAS:i:");
      bofprintf_int64(buffered_output_file,match_trace->swg_score);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
//  BC  Z  Barcode sequence
//  BQ  Z  Offset to base alignment quality (BAQ), of the same length as the read sequence.
//         At the i-th read base, BAQi = Qi - (BQi - 64) where Qi is the i-th base quality. // FIXME
//  CC  Z  Reference name of the next hit; "=" for the same chromosome
//  CM  i  Edit distance between the color sequence and the color reference (see also NM)
//  CP  i  Leftmost coordinate of the next hit
//  CQ  Z  Color read quality on the original strand of the read. Same encoding as QUAL; same length as CS.
//  CS  Z  Color read sequence on the original strand of the read. The primer base must be included.
//  cs  Z  Casava TAG (if any)
void output_sam_print_opt_field_tag_cs() {}
//  E2  Z  The 2nd most likely base calls. Same encoding and same length as QUAL.
//  FI  i  The index of segment in the template.
//  FS  Z  Segment suffix.
//  FZ  B,S Flow signal intensities on the original strand of the read, stored as (uint16 t) round(value * 100.0).
//  LB  Z  Library. Value to be consistent with the header RG-LB tag if @RG is present.
void output_sam_print_opt_field_tag_LB() {}
//  H0  i  Number of perfect hits
//  H1  i  Number of 1-difference hits (see also NM)
//  H2  i  Number of 2-difference hits
//  HI  i  Query hit index, indicating the alignment record is the i-th one stored in SAM
//  IH  i  Number of stored alignments in SAM that contains the query in the current record
//  MD  Z  String for mismatching positions (complementary to CIGAR string). Regex : [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
void output_sam_print_opt_field_tag_MD(
    buffered_output_file_t* const buffered_output_file,
    const matches_t* const matches,
    const match_trace_t* const match_trace) {
  // Match CIGAR
  const uint64_t cigar_length = match_trace->match_alignment.cigar_length;
  const cigar_element_t* const cigar_array =
      vector_get_mem(matches->cigar_vector,cigar_element_t) + match_trace->match_alignment.cigar_offset;
  // Reserve (upper bound)
  buffered_output_file_reserve(buffered_output_file,6+cigar_length*(INT_MAX_LENGTH+1));
  // Print MD
  bofprintf_string_literal(buffered_output_file,"\tMD:Z:");
  // Print MD CIGAR
  const uint8_t* const text = match_trace->text;
  const uint64_t effective_length = match_trace->match_alignment.effective_length;
  int64_t i;
  uint64_t matching_length = 0, text_pos = 0;
  if (match_trace->strand==Forward) {
    for (i=0;i<cigar_length;++i) {
      output_sam_print_md_cigar_element(buffered_output_file,
          text,effective_length,&text_pos,Forward,cigar_array+i,&matching_length);
    }
  } else {
    for (i=cigar_length-1;i>=0;--i) {
      output_sam_print_md_cigar_element(buffered_output_file,
          text,effective_length,&text_pos,Reverse,cigar_array+i,&matching_length);
    }
  }
  if (matching_length > 0) bofprintf_uint64(buffered_output_file,matching_length);
}
//  md  Z  GEM CIGAR String
void output_sam_print_opt_field_tag_md() {}
//  MQ  i  Mapping quality of the mate/next segment
void output_sam_print_opt_field_tag_MQ() {}
//  NH  i  Number of reported alignments that contains the query in the current record
void output_sam_print_opt_field_tag_NH() {}
//  NM  i  Edit distance to the reference, including ambiguous bases but excluding clipping
void output_sam_print_opt_field_tag_NM(
    buffered_output_file_t* const buffered_output_file,
    const matches_t* const matches,
    const match_trace_t* const match_trace) {
  // Reserve
  buffered_output_file_reserve(buffered_output_file,6+INT_MAX_LENGTH);
  // Calculate edit-distance
  const uint64_t edit_distance =
      matches_cigar_compute_edit_distance__excluding_clipping(
          matches->cigar_vector,match_trace->match_alignment.cigar_offset,
          match_trace->match_alignment.cigar_length);
  bofprintf_string_literal(buffered_output_file,"\tNM:i:");
  bofprintf_uint64(buffered_output_file,edit_distance);
}
//  OQ  Z  Original base quality (usually before recalibration). Same encoding as QUAL.
//  OP  i  Original mapping position (usually before realignment)
//  OC  Z  Original CIGAR (usually before realignment)
//  PG  Z  Program. Value matches the header PG-ID tag if @PG is present.
//  PQ  i  Phred likelihood of the template, conditional on both the mapping being correct
void output_sam_print_opt_field_tag_PQ() {}
//  PU  Z  Platform unit. Value to be consistent with the header RG-PU tag if @RG is present.
//  Q2  Z  Phred quality of the mate/next segment. Same encoding as QUAL.
//  R2  Z  Sequence of the mate/next segment in the template.
//  RG  Z  Read group. Value matches the header RG-ID tag if @RG is present in the header.
void output_sam_print_opt_field_tag_RG(
    buffered_output_file_t* const buffered_output_file,
    string_t* const read_group) {
  const uint64_t rg_length = string_get_length(read_group);
  // Reserve
  buffered_output_file_reserve(buffered_output_file,6+rg_length);
  // Print Field
  bofprintf_string_literal(buffered_output_file,"\tRG:Z:");
  bofprintf_string(buffered_output_file,rg_length,string_get_buffer(read_group));
}
//  SA  Z  Supplementary alignment information for chimeric alignments
void output_sam_print_opt_field_tag_SA() {
// TODO
  /*
   * H.Sapiens.1M.Illumina.l100.low.000433219  0 chr1  162259007 0 64M2I13M21S * 0 0
   *   CGTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTTTGCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCC
   *   XF@BMCCLECD>K@ZMGADO=HFFNLAHFIOHLH?GILGOODFCGIIFFIMFOJEQCEGJHIGNFLFHHJIKGEHHKHIJHIIHHIIIGIHHHIHHIIII
   *   NM:i:3  MD:Z:1T75 AS:i:67 XS:i:65 SA:Z:chr1,38528339,+,64S36M,0,0;
   *
   * H.Sapiens.1M.Illumina.l100.low.000433219  2048  chr1  38528339  0 64H36M  * 0 0
   *   TGCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCC
   *   FLFHHJIKGEHHKHIJHIIHHIIIGIHHHIHHIIII
   *   NM:i:0  MD:Z:36 AS:i:36 XS:i:36 SA:Z:chr1,162259007,+,64M2I13M21S,0,3;
   */
}
//  SM  i  Template-independent mapping quality
//  TC  i  The number of segments in the template.
//  U2  Z  Phred probability of the 2nd call being wrong conditional on the best being wrong. The same encoding as QUAL.
//  UQ  i  Phred likelihood of the segment, conditional on the mapping being correct
void output_sam_print_opt_field_tag_UQ() {}
//  TQ  i  Custom tag - equivalent to MAPQ score for a template
void output_sam_print_opt_field_tag_TQ() {}
//  TP  i  Custom tag - as TQ but comparing to all possible pairings of mappings (not taking account of orientation, contig location or interval size)
void output_sam_print_opt_field_tag_TP() {}
/*
 *   X?  ?  Reserved fields for end users (together with Y? and Z?)
 */
//  X0  i  Number of best hits
//  X1  i  Number of suboptimal hits found by BWA
//  XN  i  Number of ambiguous bases in the reference
//  XM  i  Number of mismatches in the alignment
//  XO  i  Number of gap opens
//  XG  i  Number of gap extentions
/*
 *  XB  A  Conversion type for bisulfite reads
 *  i.e.
 *    XB:A:U => non_converted
 *    XB:A:C => C2T strand
 *    XB:A:G => G2A strand
 *    XB:A:M => Mixed (invalid)
 */
void output_sam_print_opt_field_tag_XB(
    buffered_output_file_t* buffered_output_file,
    const match_trace_t* const match_trace) {
	if (match_trace->bs_strand == bs_strand_C2T) {
	  bofprintf_string_literal(buffered_output_file,"\tXB:A:C");
	} else if (match_trace->bs_strand == bs_strand_G2A) {
	  bofprintf_string_literal(buffered_output_file,"\tXB:A:G");
	}
}
/*
 *  XT  A  Type: Unique/Repeat/N/Mate-sw
 *  i.e.
 *    XT:A:U => Unique alignment
 *    XT:A:R => Repeat
 *    XT:A:N => Not mapped
 *    XT:A:M => Mate-sw (Read is fixed due to paired end rescue)
 */
void output_sam_print_opt_field_tag_XT() {}
/* Similar to XT, but for a template */
void output_sam_print_opt_field_tag_XP() {}
//  XS  i  Suboptimal alignment score
void output_sam_print_opt_field_tag_XS(
    buffered_output_file_t* const buffered_output_file,
    const match_alignment_model_t match_alignment_model,
    const match_trace_t* const match_trace) {
  // Reserve
  buffered_output_file_reserve(buffered_output_file,6+INT_MAX_LENGTH);
  // Switch alignment model
  switch (match_alignment_model) {
    case match_alignment_model_none: break;
    case match_alignment_model_hamming:
    case match_alignment_model_levenshtein:
      bofprintf_string_literal(buffered_output_file,"\tXS:i:");
      bofprintf_uint64(buffered_output_file,match_trace->edit_distance);
      break;
    case match_alignment_model_gap_affine:
      bofprintf_string_literal(buffered_output_file,"\tXS:i:");
      bofprintf_uint64(buffered_output_file,match_trace->swg_score);
      break;
    default:
      GEM_INVALID_CASE();
      break;
  }
}
//  XF  ?  Support from forward/reverse alignment
//  XE  i  Number of supporting seeds

/*
 * SAM Flag
 *
 * 0x1   (Bit 0)  => Template having multiple segments in sequencing
 *                   (The read was part of a pair during sequencing) [read paired]
 * 0x2   (Bit 1)  => Each segment properly aligned according to the aligner
 *                   (The read is mapped in a pair) [read mapped in proper pair]
 * 0x4   (Bit 2)  => Segment unmapped. The query sequence is unmapped [read unmapped]
 * 0x8   (Bit 3)  => Next segment in the template unmapped. The mate is unmapped [mate unmapped]
 * 0x10  (Bit 4)  => SEQ being reverse complemented. Strand of query (0=forward 1=reverse) [read reverse strand]
 * 0x20  (Bit 5)  => SEQ of the next segment in the template being reversed [mate reverse strand]
 * 0x40  (Bit 6)  => The first segment in the template [first in pair]
 * 0x80  (Bit 7)  => The last segment in the template [second in pair]
 * 0x100 (Bit 8)  => Secondary alignment [not primary alignment]
 * 0x200 (Bit 9)  => Not passing quality controls [read fails platform/vendor quality checks]
 * 0x400 (Bit 10) => PCR or optical duplicate [read is PCR or optical duplicate]
 * 0x800 (Bit 11) => Supplementary alignment [second or subsequent segment in a chimeric alignment]
 *
 * - RULE1:: Bit 0x4 is the only reliable place to tell whether the segment is unmapped. If 0x4 is set,
 *     no assumptions can be made about RNAME, POS, CIGAR, MAPQ, bits 0x2, 0x10 and 0x100
 *     and the bit 0x20 of the next segment in the template.
 * - RULE2:: If 0x40 and 0x80 are both set, the segment is part of a linear template, but it is neither
 *     the first nor the last segment. If both 0x40 and 0x80 are unset, the index of the segment
 *     in the template is unknown. This may happen for a non-linear template or the index is
 *     lost in data processing.
 * - RULE3:: Bit 0x100 marks the alignment not to be used in certain analyses when the tools in use
 *     are aware of this bit.
 * - RULE4:: If 0x1 is unset, no assumptions can be made about 0x2, 0x8, 0x20, 0x40 and 0x80.
*/
uint16_t output_sam_calculate_flag_se(
    const bool read_mapped,
    const strand_t read_strand,
    const bool secondary_alignment,
    const bool supplementary_alignment,
    const bool not_passing_QC,
    const bool PCR_duplicate) {
  uint16_t sam_flag = 0; // (**RULE4)
  /* 0x4  */
  if (!read_mapped) { // (**RULE1)
    sam_flag |= SAM_FLAG_UNMAPPED;
  } else {
    /* 0x10 */  if (read_strand==Reverse) sam_flag |= SAM_FLAG_REVERSE_COMPLEMENT;
    /* 0x100 */ if (secondary_alignment) sam_flag |= SAM_FLAG_SECONDARY_ALIGNMENT;
    /* 0x800 */ if (supplementary_alignment) sam_flag |= SAM_FLAG_SUPPLEMENTARY_ALIGNMENT;
  }
  /* 0x200 */
  if (not_passing_QC) sam_flag |= SAM_FLAG_NOT_PASSING_QC;
  /* 0x400 */
  if (PCR_duplicate) sam_flag |= SAM_FLAG_PCR_OR_OPTICAL_DUPLICATE;
  return sam_flag;
}
uint16_t output_sam_calculate_flag_pe(
    const bool read_paired,
    const bool read_mapped,
    const bool mate_mapped,
    const strand_t read_strand,
    const strand_t mate_strand,
    const bool first_in_pair,
    const bool last_in_pair,
    const bool secondary_alignment,
    const bool supplementary_alignment,
    const bool not_passing_QC,
    const bool PCR_duplicate) {
  /* 0x1 */
  uint16_t sam_flag = SAM_FLAG_MULTIPLE_SEGMENTS;
  /* 0x4  */
  if (!read_mapped) { // (**RULE1)
    sam_flag |= SAM_FLAG_UNMAPPED;
    /* 0x8  */
    if (!mate_mapped) {
      sam_flag |= SAM_FLAG_NEXT_UNMAPPED;
    } else {
      /* 0x20 */if (mate_strand==Reverse) sam_flag |= SAM_FLAG_NEXT_REVERSE_COMPLEMENT;
    }
  } else {
    /* 0x10  */ if (read_strand==Reverse) sam_flag |= SAM_FLAG_REVERSE_COMPLEMENT;
    /* 0x100 */ if (secondary_alignment) sam_flag |= SAM_FLAG_SECONDARY_ALIGNMENT;
    /* 0x800 */ if (supplementary_alignment) sam_flag |= SAM_FLAG_SUPPLEMENTARY_ALIGNMENT;
    /* 0x8 */
    if (!mate_mapped) {
      sam_flag |= SAM_FLAG_NEXT_UNMAPPED;
    } else {
      /* 0x2 Each segment properly aligned can only take place if both ends are mapped */
      if (read_paired) sam_flag |= SAM_FLAG_PROPERLY_ALIGNED;
      /* 0x20 */
      if (mate_strand==Reverse) sam_flag |= SAM_FLAG_NEXT_REVERSE_COMPLEMENT;
    }
  }
  /* 0x40 */  if (first_in_pair)  sam_flag |= SAM_FLAG_FIRST_SEGMENT;
  /* 0x80 */  if (last_in_pair)   sam_flag |= SAM_FLAG_LAST_SEGMENT;
  /* 0x200 */ if (not_passing_QC) sam_flag |= SAM_FLAG_NOT_PASSING_QC;
  /* 0x400 */ if (PCR_duplicate)  sam_flag |= SAM_FLAG_PCR_OR_OPTICAL_DUPLICATE;
  return sam_flag;
}
uint16_t output_sam_calculate_flag_se_match(
    const match_trace_t* const match,
    const bool secondary_alignment,
    const bool supplementary_alignment,
    const bool not_passing_QC,
    const bool PCR_duplicate) {
  return (gem_expect_true(match!=NULL)) ?
      output_sam_calculate_flag_se(true,match->strand,secondary_alignment,
          supplementary_alignment,not_passing_QC,PCR_duplicate): // Mapped
      output_sam_calculate_flag_se(false,Forward,secondary_alignment,
          supplementary_alignment,not_passing_QC,PCR_duplicate); // Unmapped
}
uint16_t output_sam_calculate_flag_pe_map(
    const match_trace_t* const match,
    const match_trace_t* const mate_match,
    const bool is_map_first_in_pair,
    const bool secondary_alignment,
    const bool supplementary_alignment,
    const bool not_passing_QC,
    const bool PCR_duplicate,
    const bool paired) {
  return output_sam_calculate_flag_pe(
      match!=NULL && mate_match!=NULL && paired,         /* read_paired   */
      match!=NULL,                                       /* read_mapped   */
      mate_match!=NULL,                                  /* mate_strand   */
      (match!=NULL) ? match->strand : Forward,           /* read_strand   */
      (mate_match!=NULL) ? mate_match->strand : Forward, /* mate_strand   */
      is_map_first_in_pair,                              /* first_in_pair */
      !is_map_first_in_pair,                             /* last_in_pair  */
      secondary_alignment,supplementary_alignment,not_passing_QC,PCR_duplicate);
}
/*
 * SAM CORE fields
 *   (QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL). No EOL is printed
 */
void output_sam_print_match_cigar(
    buffered_output_file_t* const buffered_output_file,
    const matches_t* const matches,
    const match_trace_t* const match,
    const output_sam_parameters_t* const output_sam_parameters) {
  // Get CIGAR buffer/length
  const uint64_t cigar_length = match->match_alignment.cigar_length;
  const cigar_element_t* const cigar_array =
      vector_get_mem(matches->cigar_vector,cigar_element_t) + match->match_alignment.cigar_offset;
  // Reserve (upper-bound)
  buffered_output_file_reserve(buffered_output_file,cigar_length*(INT_MAX_LENGTH+1));
  // Generate CIGAR
  int64_t i;
  if (output_sam_parameters->print_mismatches) {
    // Traverse all CIGAR elements
    if (match->strand==Forward) {
      for (i=0;i<cigar_length;++i) {
        output_sam_print_cigar_element(buffered_output_file,cigar_array+i);
      }
    } else {
      for (i=cigar_length-1;i>=0;--i) {
        output_sam_print_cigar_element(buffered_output_file,cigar_array+i);
      }
    }
  } else {
    // Traverse all CIGAR elements
    uint64_t matching_length=0;
    if (match->strand==Forward) {
      for (i=0;i<cigar_length;++i) {
        output_sam_print_reduced_cigar_element(buffered_output_file,cigar_array+i,&matching_length);
      }
    } else {
      for (i=cigar_length-1;i>=0;--i) {
        output_sam_print_reduced_cigar_element(buffered_output_file,cigar_array+i,&matching_length);
      }
    }
    if (matching_length > 0) {
      bofprintf_uint64(buffered_output_file,matching_length);
      bofprintf_char(buffered_output_file,'M');
    }
  }
}
/*
 * X4/X5  Z  GEM-Compatibility Flags
 *   X4  Column 4 of GEM output (MAP) => MAP Counters
 *   X5  Column 5 of GEM output (MAP) => MAPs
 */
void output_sam_print_opt_field_tag_X4(
    buffered_output_file_t* const buffered_output_file,
    const matches_t* const matches) {
  // Print TAG
  buffered_output_file_reserve(buffered_output_file,6); // Reserve
  bofprintf_string_literal(buffered_output_file,"\tX4:Z:");
  // Print COUNTERS
  output_map_print_counters(buffered_output_file,matches->counters,matches->max_complete_stratum,false);
}
void output_sam_print_opt_field_tag_X4_pe(
    buffered_output_file_t* const buffered_output_file,
    paired_matches_t* const paired_matches) {
  // Print TAG
  buffered_output_file_reserve(buffered_output_file,6); // Reserve
  bofprintf_string_literal(buffered_output_file,"\tX4:Z:");
  // Print COUNTERS
  output_map_print_counters(
      buffered_output_file,paired_matches->counters,
      paired_matches_get_max_complete_stratum(paired_matches),false);
}
void output_sam_print_opt_field_tag_X5(
    buffered_output_file_t* const buffered_output_file,
    const matches_t* const matches) {
  // Print TAG
  buffered_output_file_reserve(buffered_output_file,6); // Reserve
  bofprintf_string_literal(buffered_output_file,"\tX5:Z:");
  // Print MAPs
  if (gem_expect_false(matches_get_num_match_traces(matches)==0)) {
    buffered_output_file_reserve(buffered_output_file,1);
    bofprintf_char(buffered_output_file,'-');
  } else {
    const uint64_t num_match_traces = matches_get_num_match_traces(matches);
    match_trace_t** const match_traces = matches_get_match_traces(matches);
    uint64_t i;
    for (i=0;i<num_match_traces;++i) {
      // Separator
      buffered_output_file_reserve(buffered_output_file,1);
      if (i > 0) bofprintf_char(buffered_output_file,',');
      // Print Match
      output_map_print_match(buffered_output_file,matches,match_traces[i],true,map_format_v2);
    }
  }
}
void output_sam_print_opt_field_tag_X5_pe(
    buffered_output_file_t* const buffered_output_file,
    paired_matches_t* const paired_matches,
    const sequence_end_t sequence_end) {
  // Print TAG
  buffered_output_file_reserve(buffered_output_file,6); // Reserve
  bofprintf_string_literal(buffered_output_file,"\tX5:Z:");
  // Print PE-MAPs
  if (paired_matches_get_num_maps(paired_matches)==0) {
    buffered_output_file_reserve(buffered_output_file,1);
    bofprintf_char(buffered_output_file,'-');
  } else {
    if (sequence_end == paired_end1) {
      matches_t* const matches_end1 = paired_matches->matches_end1;
      const uint64_t num_maps = paired_matches_get_num_maps(paired_matches);
      paired_map_t** const paired_maps = paired_matches_get_maps(paired_matches);
      uint64_t i;
      for (i=0;i<num_maps;++i) {
        // Separator
        buffered_output_file_reserve(buffered_output_file,1);
        if (i > 0) bofprintf_char(buffered_output_file,',');
        // Print Match
        match_trace_t* const match_end1 = paired_maps[i]->match_trace_end1;
        output_map_print_match(buffered_output_file,matches_end1,match_end1,true,map_format_v2);
      }
    } else {
      matches_t* const matches_end2 = paired_matches->matches_end2;
      const uint64_t num_maps = paired_matches_get_num_maps(paired_matches);
      paired_map_t** const paired_maps = paired_matches_get_maps(paired_matches);
      uint64_t i;
      for (i=0;i<num_maps;++i) {
        // Separator
        buffered_output_file_reserve(buffered_output_file,1);
        if (i > 0) bofprintf_char(buffered_output_file,',');
        // Print Match
        match_trace_t* const match_end2 = paired_maps[i]->match_trace_end2;
        output_map_print_match(buffered_output_file,matches_end2,match_end2,true,map_format_v2);
      }
    }
  }
}
/*
 * XA  Z  Alternative hits; format: (chr,pos,CIGAR,NM;)*
 *   XA maps (chr12,+91022,101M,0) =>
 *     XA:Z:chr10,+100306093,100M,0;chr2,+173184544,100M,0;chr2,+175137710,100M,0;
 */
void output_sam_print_opt_field_tag_XA_match(
    buffered_output_file_t* const buffered_output_file,
    const matches_t* const matches,
    match_trace_t* const subdominant_match,
    const output_sam_parameters_t* const output_sam_parameters) {
  uint64_t seq_length = gem_strlen(subdominant_match->sequence_name);
  // Print SeqName
  buffered_output_file_reserve(buffered_output_file,seq_length+INT_MAX_LENGTH+10);
  bofprintf_string(buffered_output_file,seq_length,subdominant_match->sequence_name);
  bofprintf_char(buffered_output_file,',');
  // Print Strand
  bofprintf_char(buffered_output_file,(subdominant_match->strand==Forward)?'+':'-');
  if (subdominant_match->bs_strand==bs_strand_C2T) {
    bofprintf_char(buffered_output_file,'C');
  } else if (subdominant_match->bs_strand==bs_strand_G2A) {
    bofprintf_char(buffered_output_file,'G');
  }
  // Print Position
  bofprintf_uint64(buffered_output_file,subdominant_match->text_position+1);
  // Print CIGAR
  bofprintf_char(buffered_output_file,',');
  output_sam_print_match_cigar(buffered_output_file,matches,subdominant_match,output_sam_parameters);
  // Print Distance
  buffered_output_file_reserve(buffered_output_file,INT_MAX_LENGTH+10);
  bofprintf_char(buffered_output_file,',');
  bofprintf_uint64(buffered_output_file,subdominant_match->edit_distance);
  bofprintf_char(buffered_output_file,';');
}
void output_sam_print_opt_field_tag_XA_se(
    buffered_output_file_t* const buffered_output_file,
    const matches_t* const matches,
    const output_sam_parameters_t* const output_sam_parameters) {
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  if (num_matches > 1) {
    match_trace_t** match_traces = matches_get_match_traces(matches);
    ++match_traces; // Skip primary
    // Reserve
    buffered_output_file_reserve(buffered_output_file,6);
    bofprintf_string_literal(buffered_output_file,"\tXA:Z:");
    uint64_t i;
    for (i=1;i<num_matches;++i,++match_traces) {
      output_sam_print_opt_field_tag_XA_match(buffered_output_file,
          matches,*match_traces,output_sam_parameters);
    }
  }
}
void output_sam_print_opt_field_tag_XA_pe_end1(
    buffered_output_file_t* const buffered_output_file,
    paired_matches_t* const paired_matches,
    const output_sam_parameters_t* const output_sam_parameters) {
  const uint64_t num_matches = paired_matches_get_num_maps(paired_matches);
  if (num_matches > 1) {
    matches_t* const matches = paired_matches->matches_end1;
    paired_map_t** const paired_map = paired_matches_get_maps(paired_matches);
    // Reserve
    buffered_output_file_reserve(buffered_output_file,6);
    bofprintf_string_literal(buffered_output_file,"\tXA:Z:");
    uint64_t i;
    for (i=1;i<num_matches;++i) { // Skip primary
      match_trace_t* const match_end1 = paired_map[i]->match_trace_end1;
      output_sam_print_opt_field_tag_XA_match(buffered_output_file,matches,match_end1,output_sam_parameters);
    }
  }
}
void output_sam_print_opt_field_tag_XA_pe_end2(
    buffered_output_file_t* const buffered_output_file,
    paired_matches_t* const paired_matches,
    const output_sam_parameters_t* const output_sam_parameters) {
  const uint64_t num_matches = paired_matches_get_num_maps(paired_matches);
  if (num_matches > 1) {
    matches_t* const matches = paired_matches->matches_end2;
    paired_map_t** const paired_map = paired_matches_get_maps(paired_matches);
    // Reserve
    buffered_output_file_reserve(buffered_output_file,6);
    bofprintf_string_literal(buffered_output_file,"\tXA:Z:");
    uint64_t i;
    for (i=1;i<num_matches;++i) { // Skip primary
      match_trace_t* const match_end2 = paired_map[i]->match_trace_end2;
      output_sam_print_opt_field_tag_XA_match(buffered_output_file,matches,match_end2,output_sam_parameters);
    }
  }
}
void output_sam_print_qname(
    buffered_output_file_t* const buffered_output_file,
    string_t* const read) {
  const char* const tag_buffer = string_get_buffer(read);
  const uint64_t tag_length = string_get_length(read);
  // Calculate the effective length
  int effective_tag_length;
  for (effective_tag_length=0;effective_tag_length<tag_length;++effective_tag_length) {
    if (tag_buffer[effective_tag_length]==SPACE) break;
  }
  // Print the plain tag (no read pair info, nor extra tag, ...nothing)
  buffered_output_file_reserve(buffered_output_file,effective_tag_length);
  bofprintf_string(buffered_output_file,effective_tag_length,string_get_buffer(read));
}
void output_sam_print_seq__qualities(
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read,
    const match_trace_t* const match,
    const bool secondary_alignment,
    const output_sam_parameters_t* const output_sam_parameters) {
  if (secondary_alignment && output_sam_parameters->omit_secondary_read__qualities) {
    buffered_output_file_reserve(buffered_output_file,4);
    bofprintf_string_literal(buffered_output_file,"\t*\t*");
  } else {
    string_t* const read = &seq_read->read;
    string_t* const qualities = &seq_read->qualities;
    const uint64_t read_length = string_get_length(read);
    if (!string_is_null(read) && !string_is_null(qualities)) {
      const uint64_t qualities_length = string_get_length(qualities);
      buffered_output_file_reserve(buffered_output_file,read_length+qualities_length+20);
      if (match==NULL || match->strand==Forward) {
        bofprintf_char(buffered_output_file,'\t');
        bofprintf_string(buffered_output_file,read_length,string_get_buffer(read));
        bofprintf_char(buffered_output_file,'\t');
        bofprintf_string(buffered_output_file,qualities_length,string_get_buffer(qualities));
      } else {
        int64_t i;
        bofprintf_char(buffered_output_file,'\t');
        const char* const read_buffer = string_get_buffer(read);
        for (i=read_length-1;i>=0;--i) {
          bofprintf_char(buffered_output_file,dna_complement(read_buffer[i]));
        }
        bofprintf_char(buffered_output_file,'\t');
        const char* const qualities_buffer = string_get_buffer(qualities);
        for (i=qualities_length-1;i>=0;--i) {
          bofprintf_char(buffered_output_file,qualities_buffer[i]);
        }
      }
    } else {
      GEM_INTERNAL_CHECK(!string_is_null(read),"Output SAM. Sequence read is null");
      buffered_output_file_reserve(buffered_output_file,read_length+20);
      if (match==NULL || match->strand==Forward) {
        bofprintf_char(buffered_output_file,'\t');
        bofprintf_string(buffered_output_file,read_length,string_get_buffer(read));
      } else {
        int64_t i;
        bofprintf_char(buffered_output_file,'\t');
        const char* const read_buffer = string_get_buffer(read);
        for (i=read_length-1;i>=0;--i) {
          bofprintf_char(buffered_output_file,dna_complement(read_buffer[i]));
        }
      }
      bofprintf_string_literal(buffered_output_file,"\t*");
    }
  }
}
void output_sam_print_core_fields_pe(
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read,
    const matches_t* const matches,
    match_trace_t* const match,
    match_trace_t* const mate,
    const uint8_t template_mapq_score,
    const int64_t template_length,
    const bool is_map_first_in_pair,
    const bool paired,
    const bool secondary_alignment,
    const bool supplementary_alignment,
    const bool not_passing_QC,
    const bool PCR_duplicate,
    const output_sam_parameters_t* const output_sam_parameters) {
  // (1) Print QNAME
  output_sam_print_qname(buffered_output_file,&seq_read->tag);
  // (2) Print FLAG
  const uint16_t flag = output_sam_calculate_flag_pe_map(
      match,mate,is_map_first_in_pair,secondary_alignment,
      supplementary_alignment,not_passing_QC,PCR_duplicate,paired);
  buffered_output_file_reserve(buffered_output_file,INT_MAX_LENGTH+1);
  bofprintf_char(buffered_output_file,'\t');
  bofprintf_uint64(buffered_output_file,flag);
  // (3) Print RNAME
  // (4) Print POS
  // (5) Print MAPQ
  // (6) Print CIGAR
  if (match!=NULL) {
    uint64_t match_seq_name_length = gem_strlen(match->sequence_name);
    buffered_output_file_reserve(buffered_output_file,match_seq_name_length+2*INT_MAX_LENGTH+10);
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_string(buffered_output_file,match_seq_name_length,match->sequence_name); // (3) RNAME
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_uint64(buffered_output_file,match->text_position+1); // (4) POS
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_int64(buffered_output_file,template_mapq_score); // (5) MAPQ
    bofprintf_char(buffered_output_file,'\t');
    output_sam_print_match_cigar(buffered_output_file,matches,match,output_sam_parameters); // (6) CIGAR
  } else if(mate!=NULL) {
    uint64_t mate_seq_name_length = gem_strlen(mate->sequence_name);
    buffered_output_file_reserve(buffered_output_file,mate_seq_name_length+INT_MAX_LENGTH+10);
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_string(buffered_output_file,mate_seq_name_length,mate->sequence_name); // (3) RNAME
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_uint64(buffered_output_file,mate->text_position+1); // (4) POS
    bofprintf_string_literal(buffered_output_file,"\t0\t*");
  } else {
    buffered_output_file_reserve(buffered_output_file,10);
    bofprintf_string_literal(buffered_output_file,"\t*\t0\t0\t*");
  }
  // (7) Print RNEXT
  // (8) Print PNEXT
  // (9) Print TLEN
  if (mate!=NULL) {
    if (match==NULL || !gem_streq(match->sequence_name,mate->sequence_name)) {
      uint64_t mate_seq_name_length = gem_strlen(mate->sequence_name);
      buffered_output_file_reserve(buffered_output_file,mate_seq_name_length+2*INT_MAX_LENGTH+10);
      bofprintf_char(buffered_output_file,'\t');
      bofprintf_string(buffered_output_file,mate_seq_name_length,mate->sequence_name);
      bofprintf_char(buffered_output_file,'\t');
      bofprintf_uint64(buffered_output_file,mate->text_position+1);
    } else {
      buffered_output_file_reserve(buffered_output_file,2*INT_MAX_LENGTH+10);
      bofprintf_string_literal(buffered_output_file,"\t=\t");
      bofprintf_uint64(buffered_output_file,mate->text_position+1);
    }
    bofprintf_char(buffered_output_file,'\t');
    const bool leftmost_segment = (match==NULL || match->text_position <= mate->text_position);
    bofprintf_int64(buffered_output_file,leftmost_segment ? template_length : -template_length);
  } else if(!secondary_alignment && match!=NULL) {
    buffered_output_file_reserve(buffered_output_file,2*INT_MAX_LENGTH);
    bofprintf_string_literal(buffered_output_file,"\t=\t");
    bofprintf_uint64(buffered_output_file,match->text_position+1);
    bofprintf_string_literal(buffered_output_file,"\t0");
  } else {
    buffered_output_file_reserve(buffered_output_file,10);
    bofprintf_string_literal(buffered_output_file,"\t*\t0\t0");
  }
  // (10) Print SEQ
  // (11) Print QUAL
  output_sam_print_seq__qualities(buffered_output_file,seq_read,match,secondary_alignment,output_sam_parameters);
}
void output_sam_print_core_fields_se(
    buffered_output_file_t* const buffered_output_file,
    sequence_t* const seq_read,
    const matches_t* const matches,
    match_trace_t* const match,
    const bool secondary_alignment,
    const bool supplementary_alignment,
    const bool not_passing_QC,
    const bool PCR_duplicate,
    const output_sam_parameters_t* const output_sam_parameters) {
  // (1) Print QNAME
  output_sam_print_qname(buffered_output_file,&seq_read->tag);
  // (2) Print FLAG
  const uint16_t flag = output_sam_calculate_flag_se_match(
      match,secondary_alignment,supplementary_alignment,not_passing_QC,PCR_duplicate);
  buffered_output_file_reserve(buffered_output_file,INT_MAX_LENGTH+1);
  bofprintf_char(buffered_output_file,'\t');
  bofprintf_uint64(buffered_output_file,flag);
  // Is mapped?
  if (match!=NULL) {
    // (3) Print RNAME
    // (4) Print POS
    // (5) Print MAPQ
    uint64_t sequence_name_length = gem_strlen(match->sequence_name);
    buffered_output_file_reserve(buffered_output_file,sequence_name_length+2*INT_MAX_LENGTH+10);
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_string(buffered_output_file,sequence_name_length,match->sequence_name);
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_uint64(buffered_output_file,match->text_position+1);
    bofprintf_char(buffered_output_file,'\t');
    bofprintf_int64(buffered_output_file,match->mapq_score);
    bofprintf_char(buffered_output_file,'\t');
    // (6) Print CIGAR
    output_sam_print_match_cigar(buffered_output_file,matches,match,output_sam_parameters);
  } else {
    // (3) Print RNAME
    // (4) Print POS
    // (5) Print MAPQ
    // (6) Print CIGAR
    buffered_output_file_reserve(buffered_output_file,20);
    bofprintf_string_literal(buffered_output_file,"\t*\t0\t0\t*");
  }
  //  (7) Print RNEXT
  //  (8) Print PNEXT
  //  (9) Print TLEN
  buffered_output_file_reserve(buffered_output_file,20);
  bofprintf_string_literal(buffered_output_file,"\t*\t0\t0");
  // (10) Print SEQ
  // (11) Print QUAL
  output_sam_print_seq__qualities(buffered_output_file,seq_read,match,secondary_alignment,output_sam_parameters);
}
void output_sam_print_optional_fields_se(
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search,
    const matches_t* const matches,
    const match_trace_t* const match_trace,
    const uint64_t match_number,
    const bool print_tag_XS,
    const output_sam_parameters_t* const sam_parameters) {
  /*
   * RG:Z:ReadGroup NM:i:2 MD:Z:98 AS:i:87 XS:i:47
   */
  const uint64_t num_matches = matches_get_num_match_traces(matches);
	// RG
	if (sam_parameters->read_group_id) {
	  output_sam_print_opt_field_tag_RG(buffered_output_file,sam_parameters->read_group_id);
	}
  // NM
  if (match_trace) output_sam_print_opt_field_tag_NM(buffered_output_file,matches,match_trace);
  /*
   * MD. We don't output MD flags for bisulfite matches because
   *   (a) they are not correct w.r.t. the stated reference
   *   (b) if we calculated them correctly they would be large as they show all mismatches
   */
  if (match_trace) {
    if (sam_parameters->bisulfite_output) {
	    output_sam_print_opt_field_tag_XB(buffered_output_file,match_trace); // XB
    } else {
	    output_sam_print_opt_field_tag_MD(buffered_output_file,matches,match_trace); // MD
    }
  }
  // AS
  const match_alignment_model_t match_alignment_model = archive_search->search_parameters.match_alignment_model;
  if (match_trace) output_sam_print_opt_field_tag_AS(buffered_output_file,match_alignment_model,match_trace);
  // XS
  if (print_tag_XS && match_number+1 < num_matches) {
    const match_trace_t* const next_match_trace = matches_get_match_traces(matches)[match_number+1];
    output_sam_print_opt_field_tag_XS(buffered_output_file,match_alignment_model,next_match_trace);
  }
}
void output_sam_print_optional_fields_pe(
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search,
    const paired_matches_t* const paired_matches,
    const sequence_end_t sequence_end,
    const match_trace_t* const match_trace,
    const uint64_t match_number,
    const bool print_tag_XS,
    const output_sam_parameters_t* const sam_parameters) {
  /*
   * RG:Z:ReadGroup NM:i:2 MD:Z:98 AS:i:87 XS:i:47
   */
  const uint64_t num_matches = paired_matches_get_num_maps(paired_matches);
  matches_t* const matches = (sequence_end==paired_end1) ? paired_matches->matches_end1 : paired_matches->matches_end2;
  // RG
  if (sam_parameters->read_group_id) {
    output_sam_print_opt_field_tag_RG(buffered_output_file,sam_parameters->read_group_id);
  }
  // NM
  if (match_trace) output_sam_print_opt_field_tag_NM(buffered_output_file,matches,match_trace);
  /*
   * MD We don't output MD flags for bisulfite matches because
   *   (a) they are not correct w.r.t. the stated reference
   *   (b) if we calculated them correctly they would be large as they would show all mismatches
   */
  if (match_trace) {
    if (sam_parameters->bisulfite_output) {
	    output_sam_print_opt_field_tag_XB(buffered_output_file,match_trace); // XB
    } else {
	    output_sam_print_opt_field_tag_MD(buffered_output_file,matches,match_trace); // MD
    }
  }
  // AS
  const match_alignment_model_t match_alignment_model = archive_search->search_parameters.match_alignment_model;
  if (match_trace) output_sam_print_opt_field_tag_AS(buffered_output_file,match_alignment_model,match_trace);
  // XS
  if (print_tag_XS && match_number+1 < num_matches) {
    const match_trace_t* const next_match_trace = matches_get_match_traces(matches)[match_number+1];
    output_sam_print_opt_field_tag_XS(buffered_output_file,match_alignment_model,next_match_trace);
  }
}
/*
 * SAM output SE
 */
void output_sam_single_end_matches(
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search,
    matches_t* const matches,
    const output_sam_parameters_t* const output_sam_parameters) {
  PROF_START_TIMER(GP_OUTPUT_SAM_SE);
  const bool supplementary_alignment = false; // TODO
  const bool not_passing_QC = false; // TODO
  const bool PCR_duplicate = false; // TODO
  // Print MATCHES
  const uint64_t vector_match_trace_used = matches_get_num_match_traces(matches);
  if (gem_expect_false(vector_match_trace_used==0)) {
    // Print Unmapped
    output_sam_print_core_fields_se(
        buffered_output_file,archive_search->sequence,matches,NULL,false,
        supplementary_alignment,not_passing_QC,PCR_duplicate,output_sam_parameters);
    // Print Optional Fields
    output_sam_print_optional_fields_se(buffered_output_file,
        archive_search,matches,NULL,0,false,output_sam_parameters);
    // GEM compatibility
    if (output_sam_parameters->print_gem_fields) {
      output_sam_print_opt_field_tag_X4(buffered_output_file,matches);
      output_sam_print_opt_field_tag_X5(buffered_output_file,matches);
    }
    // EOL
    buffered_output_file_reserve(buffered_output_file,1);
    bofprintf_char(buffered_output_file,'\n');
  } else {
    // Traverse all matches (Position-matches)
    const uint64_t num_match_traces = matches_get_num_match_traces(matches);
    match_trace_t** const match_traces = matches_get_match_traces(matches);
    uint64_t i;
    for (i=0;i<num_match_traces;++i) {
      // Print Core Fields
      const bool secondary_alignment = (i>0);
      output_sam_print_core_fields_se(buffered_output_file,archive_search->sequence,matches,match_traces[i],
          secondary_alignment,supplementary_alignment,not_passing_QC,PCR_duplicate,output_sam_parameters);
      // Print Optional Fields
      output_sam_print_optional_fields_se(buffered_output_file,archive_search,
          matches,match_traces[i],i,false,output_sam_parameters);
      // GEM compatibility
      if (!secondary_alignment) {
        if (output_sam_parameters->print_gem_fields) {
          output_sam_print_opt_field_tag_X4(buffered_output_file,matches);
          output_sam_print_opt_field_tag_X5(buffered_output_file,matches);
        } else if (output_sam_parameters->compact_xa) {
          // Print XA (sub-dominant alignments)
          output_sam_print_opt_field_tag_XA_se(buffered_output_file,matches,output_sam_parameters);
          // EOL
          buffered_output_file_reserve(buffered_output_file,1);
          bofprintf_char(buffered_output_file,'\n');
          // Exit
          break;
        }
      }
      // EOL
      buffered_output_file_reserve(buffered_output_file,1);
      bofprintf_char(buffered_output_file,'\n');
    }
  }
  PROF_STOP_TIMER(GP_OUTPUT_SAM_SE);
}
void output_sam_paired_end_matches_unpaired(
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search,
    const sequence_end_t sequence_end,
    matches_t* const matches,
    match_trace_t* const primary_match_mate,
    const output_sam_parameters_t* const output_sam_parameters) {
  const bool supplementary_alignment = false;   // TODO
  const bool not_passing_QC = false;   // TODO
  const bool PCR_duplicate = false;  // TODO
  const bool is_map_first_in_pair = sequence_end==paired_end1;
  const uint64_t num_matches = matches_get_num_match_traces(matches);
  if (num_matches == 0) { // End Unmapped
     // Print Core Fields
     output_sam_print_core_fields_pe(buffered_output_file,
        archive_search->sequence,matches,NULL,primary_match_mate,0,0,is_map_first_in_pair,
        false,false,supplementary_alignment,not_passing_QC,PCR_duplicate,output_sam_parameters);
     // Print Optional Fields
     output_sam_print_optional_fields_se(buffered_output_file,archive_search,matches,NULL,0,false,output_sam_parameters);
     // GEM compatibility
     if (output_sam_parameters->print_gem_fields) {
        output_sam_print_opt_field_tag_X4(buffered_output_file,matches);
        output_sam_print_opt_field_tag_X5(buffered_output_file,matches);
     }
     // EOL
     buffered_output_file_reserve(buffered_output_file,1);
     bofprintf_char(buffered_output_file,'\n');
  } else { // End Mapped
    const uint64_t num_match_traces = matches_get_num_match_traces(matches);
    match_trace_t** const match_traces = matches_get_match_traces(matches);
    uint64_t i;
    for (i=0;i<num_match_traces;++i) {
        // Print Core Fields
        const bool secondary_alignment = (i>0);
        output_sam_print_core_fields_pe(buffered_output_file,
           archive_search->sequence,matches,match_traces[i],i==0?primary_match_mate:NULL,
           match_traces[i]->mapq_score,0,is_map_first_in_pair,false,secondary_alignment,
           supplementary_alignment,not_passing_QC,PCR_duplicate,output_sam_parameters);
        // Print Optional Fields
        output_sam_print_optional_fields_se(buffered_output_file,archive_search,
            matches,match_traces[i],i,false,output_sam_parameters);
        // GEM compatibility
        if (!secondary_alignment) {
           if (output_sam_parameters->print_gem_fields) {
              output_sam_print_opt_field_tag_X4(buffered_output_file,matches);
              output_sam_print_opt_field_tag_X5(buffered_output_file,matches);
           } else if (output_sam_parameters->compact_xa) {
              output_sam_print_opt_field_tag_XA_se(buffered_output_file,matches,output_sam_parameters); // Print XA
              // EOL
              buffered_output_file_reserve(buffered_output_file,1);
              bofprintf_char(buffered_output_file,'\n');
              // Exit
              break;
           }
        }
        // EOL
        buffered_output_file_reserve(buffered_output_file,1);
        bofprintf_char(buffered_output_file,'\n');
     }
  }
}
void output_sam_paired_end_matches_paired(
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,
    const output_sam_parameters_t* const output_sam_parameters) {
  // Parameters
  const bool supplementary_alignment = false; // TODO
  const bool not_passing_QC = false; // TODO
  const bool PCR_duplicate = false; // TODO
  matches_t* const matches_end1 = paired_matches->matches_end1;
  matches_t* const matches_end2 = paired_matches->matches_end2;
  // Traverse all matches (Position-matches)
  const uint64_t num_paired_matches = paired_matches_get_num_maps(paired_matches);
  paired_map_t** const paired_maps = paired_matches_get_maps(paired_matches);
  uint64_t match_number;
  for (match_number=0;match_number<num_paired_matches;++match_number) {
    const bool secondary_alignment = (match_number>0);
    paired_map_t* const paired_map = paired_maps[match_number];
    match_trace_t* const match_end1 = paired_map->match_trace_end1;
    match_trace_t* const match_end2 = paired_map->match_trace_end2;
    /*
     * End/1
     */
    // Print Core Fields
    output_sam_print_core_fields_pe(
        buffered_output_file,archive_search_end1->sequence,
        matches_end1,match_end1,match_end2,paired_map->mapq_score,
        paired_map->template_length,true,true,secondary_alignment,
        supplementary_alignment,not_passing_QC,PCR_duplicate,output_sam_parameters);
    // Print Optional Fields
    output_sam_print_optional_fields_pe(buffered_output_file,archive_search_end1,paired_matches,
        paired_end1,match_end1,match_number,!output_sam_parameters->compact_xa,output_sam_parameters);
    // GEM compatibility
    if (!secondary_alignment) {
      if (output_sam_parameters->print_gem_fields) {
        output_sam_print_opt_field_tag_X4_pe(buffered_output_file,paired_matches);
        output_sam_print_opt_field_tag_X5_pe(buffered_output_file,paired_matches,paired_end1);
      } else if (output_sam_parameters->compact_xa) {
        output_sam_print_opt_field_tag_XA_pe_end1(buffered_output_file,paired_matches,output_sam_parameters); // Print XA
      }
    }
    // EOL
    buffered_output_file_reserve(buffered_output_file,1);
    bofprintf_char(buffered_output_file,'\n');
    /*
     * End/2
     */
    // Print Core Fields
    output_sam_print_core_fields_pe(
        buffered_output_file,archive_search_end2->sequence,
        matches_end2,match_end2,match_end1,paired_map->mapq_score,
        paired_map->template_length,false,true,secondary_alignment,
        supplementary_alignment,not_passing_QC,PCR_duplicate,output_sam_parameters);
    // Print Optional Fields
    output_sam_print_optional_fields_pe(buffered_output_file,archive_search_end2,paired_matches,
        paired_end2,match_end2,match_number,!output_sam_parameters->compact_xa,output_sam_parameters);
    // GEM compatibility
    if (!secondary_alignment) {
      if (output_sam_parameters->print_gem_fields) {
        output_sam_print_opt_field_tag_X4_pe(buffered_output_file,paired_matches);
        output_sam_print_opt_field_tag_X5_pe(buffered_output_file,paired_matches,paired_end2);
      } else if (output_sam_parameters->compact_xa) {
        output_sam_print_opt_field_tag_XA_pe_end2(buffered_output_file,paired_matches,output_sam_parameters); // Print XA
        // EOL
        buffered_output_file_reserve(buffered_output_file,1);
        bofprintf_char(buffered_output_file,'\n');
        // Exit
        break;
      }
    }
    // EOL
    buffered_output_file_reserve(buffered_output_file,1);
    bofprintf_char(buffered_output_file,'\n');
  }
}
void output_sam_paired_end_matches(
    buffered_output_file_t* const buffered_output_file,
    archive_search_t* const archive_search_end1,
    archive_search_t* const archive_search_end2,
    paired_matches_t* const paired_matches,
    const output_sam_parameters_t* const output_sam_parameters) {
  PROF_START_TIMER(GP_OUTPUT_SAM_PE);
  if (gem_expect_false(!paired_matches_is_mapped(paired_matches))) {
    // Parameters
    matches_t* const matches_end1 = paired_matches->matches_end1;
    matches_t* const matches_end2 = paired_matches->matches_end2;
    const uint64_t num_matches_end1 = matches_get_num_match_traces(matches_end1);
    const uint64_t num_matches_end2 = matches_get_num_match_traces(matches_end2);
    match_trace_t* const primary_match_end1 = num_matches_end1 ? matches_get_primary_match(matches_end1) : NULL;
    match_trace_t* const primary_match_end2 = num_matches_end2 ? matches_get_primary_match(matches_end2) : NULL;
    // End1
    output_sam_paired_end_matches_unpaired(buffered_output_file,archive_search_end1,paired_end1,
        matches_end1,primary_match_end2,output_sam_parameters);
    // End2
    output_sam_paired_end_matches_unpaired(buffered_output_file,archive_search_end2,paired_end2,
        matches_end2,primary_match_end1,output_sam_parameters);
  } else {
    output_sam_paired_end_matches_paired(buffered_output_file,
        archive_search_end1,archive_search_end2,paired_matches,output_sam_parameters);
  }
  PROF_STOP_TIMER(GP_OUTPUT_SAM_PE);
}
