/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Simon Heath  <simon.heath@gmail.com>
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
 * AUTHOR(S): Simon Heath  <simon.heath@gmail.com>
 */

#include "stats/report_stats.h"
#include "archive/archive.h"
#include "archive/locator.h"
#include "archive/search/archive_search.h"
#include "archive/search/archive_search_se_parameters.h"
#include "matches/matches.h"
#include "matches/matches_cigar.h"
#include "stats/report_stats_mstats.h"
#include "text/dna_text.h"
#include "text/text_trace.h"
#include "utils/string_buffer.h"
#include "utils/vector.h"

int btab_fwd[256]={ ['A'] = 1, ['C'] = 2, ['G'] = 3, ['T'] = 4 };
int btab_rev[256]={ ['A'] = 4, ['C'] = 3, ['G'] = 2, ['T'] = 1 };

void update_counts(
    sequence_t* const seq_read,
    mapping_stats_t* const mstats,
    const int end) {
	 string_t* const read = &seq_read->read;
	 const uint64_t len = string_get_length(read);
	 uint64_t* count;
	 ihash_element_t* ih = ihash_get_ihash_element(mstats->read_length_dist[end],len);
	 if(ih == NULL) {
			count=mm_alloc(uint64_t);
			*count=1;
			ihash_insert_element(mstats->read_length_dist[end],len,count);
	 } else {
			count = ih->element;
			(*count)++;
	 }
	 char *p = string_get_buffer(read);
	 uint64_t *q = mstats->overall_counts.counts[end];
	 while(*p) q[btab_fwd[(int)*p++]]++;
}
void update_distance_counts(
    matches_t const * matches,
    match_trace_t* const match_trace,
    mapping_stats_t* mstats, int end) {
	 vector_t* const cigar_vector = matches->cigar_vector;
	 const uint64_t cigar_buffer_offset = match_trace->match_alignment.cigar_offset;
	 const uint64_t cigar_length = match_trace->match_alignment.cigar_length;
	 const cigar_element_t* const cigar_buffer = vector_get_elm(cigar_vector,cigar_buffer_offset,cigar_element_t);
	 uint64_t i, distance = 0;
	 for (i=0;i<cigar_length;++i) {
			if(cigar_buffer[i].type == cigar_mismatch) distance++;
	 }
	 uint64_t* count;
	 ihash_element_t* ih = ihash_get_ihash_element(mstats->distance_dist[end],distance);
	 if(ih == NULL) {
			count=mm_alloc(uint64_t);
			*count=1;
			ihash_insert_element(mstats->distance_dist[end],distance,count);
	 } else {
			count = ih->element;
			(*count)++;
	 }
}

int get_text(const archive_search_t *archive, const match_trace_t *match, text_trace_t *text, bs_strand_t bs) {
    
    // Get locator interval
    locator_interval_t* const locator_interval = 
        locator_inverse_map(archive->archive->locator, (const uint8_t *)match->sequence_name, Forward, bs, match->text_position);
    if(locator_interval==NULL || locator_interval->type==locator_interval_uncalled) return -1;
    
    // Adjust sequence boundaries
    const uint64_t index_begin_position = locator_interval->begin_position +
        match->text_position-locator_interval->sequence_offset;
    uint64_t index_end_position = index_begin_position + match->text_length;
    if (index_end_position > locator_interval->end_position) {
      index_end_position = locator_interval->end_position;
    }
    uint64_t text_length = index_end_position-index_begin_position;
    if(text_length!=match->text_length) return -1;
    
    // Retrieve the sequence
    archive_text_retrieve(
        archive->archive->text,index_begin_position,
        text_length,match->strand==Reverse,false,
        text,archive->mm_allocator);
    
    return 0;
}

const uint8_t encode_prev[49]= {
    3,3,3,3,3,3,3, 2,0,2,1,2,2,2, 3,3,3,3,3,3,3, 3,3,3,3,3,3,3, 4,4,4,4,4,4,4, 4,4,4,4,4,4,4, 4,4,4,4,4,4,4
};

const uint8_t encode_curr[49]= {
    0x31,0x32,0x33,0x34,0x30,0x30,0x30,
    0x31,0x32,0x33,0x34,0x30,0x30,0x30,
    0x11,0x22,0x03,0x24,0x20,0x20,0x20,
    0x31,0x32,0x33,0x34,0x30,0x30,0x30,
    0x41,0x42,0x43,0x44,0x40,0x40,0x40,
    0x41,0x42,0x43,0x44,0x40,0x40,0x40,
    0x41,0x42,0x43,0x44,0x40,0x40,0x40,
};

const uint8_t encode_cts_c2t[25] = {7,7,7,5,0, 8,8,8,6,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
const uint8_t encode_cts_g2a[25] = {7,8,0,0,0, 7,8,0,0,0, 7,8,0,0,0, 5,6,0,0,0, 0,0,0,0,0};

void update_conversion_counts(
    const archive_search_t *archive,
    const matches_t* const matches,
    const match_trace_t* const match,
    const mapping_stats_t* mstats,
    const int read_idx,
    const int end) {

     /*
      * bs strand
      *
      * 0 => None
      * 1 => C2T
      * 2 => G2A
      * 3 => mixed
      */
      
     strand_t strand = match->strand;
     bs_strand_t bs = match->bs_strand;
	 if(bs==bs_strand_C2T || bs==bs_strand_G2A) {	
        sequence_t * seq_read = archive->sequence;
        const uint64_t cigar_length = match->match_alignment.cigar_length;
        const cigar_element_t* const cigar_array=vector_get_mem(matches->cigar_vector,cigar_element_t) + match->match_alignment.cigar_offset;
   		const search_parameters_t *params = &archive->search_parameters;
        //const int (*cnv_idx)[25] = bs==1?&cnv_idx_c2t:&cnv_idx_g2a;
        
        int clip = params->conversion_clip_start;
        string_t *read = &seq_read->read;
		const char *p = string_get_buffer(read);
		uint64_t *ct = mapping_stats_base_counts(mstats, bs - 1, read_idx)->counts[end];
		int initial_skip = cigar_array[0].type==cigar_del?cigar_array[0].length : 0;
		int final_skip = cigar_array[cigar_length-1].type==cigar_del?cigar_array[cigar_length-1].length : 0;
		if(clip < initial_skip) clip = initial_skip;

		int l = read->length - final_skip;
		/* We don't check that l is not negative (it shouldn't be!), but since clip is >=0, if it is negative it
		 * will be detected in the next line.
		 */
		if(l<=clip) return;
		
		// Get text from other bisulfite strand
		text_trace_t c2t_text, g2a_text;
		if(get_text(archive, match, &c2t_text, bs_strand_C2T)) return;
		if(get_text(archive, match, &g2a_text, bs_strand_G2A)) {
            text_trace_destroy(&c2t_text, archive->mm_allocator);
            return;
		}
		const uint8_t *t1, *t2, (*encode_cts)[];
		if(strand == Forward) {
		    t1 = c2t_text.text;
			t2 = g2a_text.text;	
			encode_cts = bs==bs_strand_C2T?&encode_cts_c2t:&encode_cts_g2a;	
		} else {
            t2 = c2t_text.text;
            t1 = g2a_text.text;
            encode_cts = bs==bs_strand_C2T?&encode_cts_g2a:&encode_cts_c2t;	
		}
		uint64_t text_len = match->text_length;
		
		uint8_t prev;
		int i=0,j=0,el;
		const cigar_element_t *elem;
		const uint8_t min_base_qual = params->conversion_min_base_qual + 33;
        const char *q = seq_read->has_qualities?string_get_buffer(&seq_read->qualities):NULL;  
        for(int k=0;k<cigar_length;k++) {
            elem=cigar_array+k;
            el = elem->length;
            switch(elem->type) {
                case cigar_match:
                case cigar_mismatch:
                    prev=4;
                    while(el-->0 && j<text_len) {
                        if(i >= clip) {
                            int rf = t1[j];
                            if(rf==3) rf = t2[j];
                            uint8_t base=p[i];
                            uint8_t qual = q==NULL?min_base_qual:q[i];
                            uint8_t sq = dna_encode_table[(int)(qual>=min_base_qual?base:'N')];
                            int ix = rf * 7 + sq; 
                            uint8_t curr = encode_curr[ix];
                            // Base counts
                            ct[(int)(curr&0xf)]++;
                            //Conversion counts
                            int cts_ix = (*encode_cts)[prev * 5 + (curr>>4)];
                            // fprintf(stderr,"OOOK\t%d\tbs %d\tstrand %s\trf %d\tsq %d\tprev %d\tcurr %d\tcts_ix %d\n", i, (int)bs, strand==Forward?"F":"R", rf, sq, prev, curr, cts_ix);
                            if(cts_ix) ct[cts_ix]++;
                            prev=encode_prev[ix];
                        }
                        i++;
                        j++;
                    }
                    break;
                case cigar_ins:
                    j+=el;
                    break;
                case cigar_del:
                    i+=el;                        
                    break;
                default:
                    GEM_INVALID_CASE();
                    break;
            }
        }		
	 }
}

/*
 * Returns index + 1 to control_sequence_t in control_seq_vector or 0 if not found
 */
int get_read_control_index(match_trace_t* const match, const vector_t *v, bool bisulfite_mode) {
	 char* seq = match->sequence_name;
	 int ix = -1;
	 int nc = vector_get_used(v);
	 for(int i = 0; i < nc; i++) {
		control_sequence_t *cs = vector_get_elm(v, i, control_sequence_t);
		if((bisulfite_mode || cs->sequence_type == SequenceControl) && !strcmp(seq, string_get_buffer(&(cs->sequence_name)))) {
		    ix = i;
			break;
		}
	 }
	 return ix + 1;
}

void collect_se_mapping_stats(
    const archive_search_t* archive_search,
    const matches_t* matches,
    mapping_stats_t* mstats) {
	 update_counts(archive_search->sequence,mstats,0);
	 bool bisulfite_mode = archive_search->archive->type == archive_dna_bisulfite;
	 int read_idx = -1;
	 uint8_t min_mapq = archive_search->search_parameters.conversion_min_mapq;
	 const uint64_t num_match_traces = matches_get_num_match_traces(matches);
	 if (gem_expect_false(num_match_traces==0)) { // Unmapped
			mstats->unmapped[0]++;
	 } else {
			// We just look at primary alignments
			match_trace_t* match = matches_get_primary_match(matches);
			update_distance_counts(matches, match, mstats, 0);
			bs_strand_t bs = match->bs_strand;
			read_idx = get_read_control_index(match, archive_search->search_parameters.control_sequences, bisulfite_mode);
			if(match->mapq_score>=min_mapq) {
				 update_conversion_counts(archive_search, matches, match, mstats, read_idx, 0);
			}
			mstats->hist_mapq[(int)match->mapq_score]++;
			if(bs == bs_strand_C2T) mstats->BSreads[0][0]++;
			else if(bs == bs_strand_G2A) mstats->BSreads[0][1]++;
			mapping_stats_reads_inc(mstats, 0, read_idx);
	 }
}
void collect_pe_mapping_stats(
    const archive_search_t* archive_search1,
    const archive_search_t* archive_search2,
 	  paired_matches_t* paired_matches,
 	  mapping_stats_t* mstats) {
     bool bisulfite_mode = archive_search1->archive->type == archive_dna_bisulfite;
	 update_counts(archive_search1->sequence,mstats,0);
	 update_counts(archive_search2->sequence,mstats,1);
	 matches_t* const matches_end1 = paired_matches->matches_end1;
	 matches_t* const matches_end2 = paired_matches->matches_end2;
	 vector_t const *control_sequences = archive_search1->search_parameters.control_sequences;
	 bs_strand_t bs1,bs2;
	 bs1 = bs2 = bs_strand_none;
	 int read_idx1 = -1, read_idx2 = -1;
	 uint8_t min_mapq = archive_search1->search_parameters.conversion_min_mapq;
	 if (gem_expect_false(!paired_matches_is_mapped(paired_matches))) { // Non paired
			const uint64_t vector_match_trace_used_end1 = matches_get_num_match_traces(matches_end1);
			const uint64_t vector_match_trace_used_end2 = matches_get_num_match_traces(matches_end2);
			if(vector_match_trace_used_end1 && vector_match_trace_used_end2) {
 				 match_trace_t* prim_match_end1 = matches_get_primary_match(matches_end1);
				 match_trace_t* prim_match_end2 = matches_get_primary_match(matches_end2);
				 update_distance_counts(matches_end1, prim_match_end1, mstats, 0);
				 update_distance_counts(matches_end2, prim_match_end2, mstats, 1);
				 bs1 = prim_match_end1 -> bs_strand;
				 bs2 = prim_match_end2 -> bs_strand;
				 read_idx1 = get_read_control_index(prim_match_end1,control_sequences, bisulfite_mode);
				 read_idx2 = get_read_control_index(prim_match_end2,control_sequences, bisulfite_mode);
			} else if(vector_match_trace_used_end1) {
				 mstats->unmapped[1]++;
				 match_trace_t* prim_match_end1 = matches_get_primary_match(matches_end1);
				 update_distance_counts(matches_end1, prim_match_end1, mstats, 0);
				 bs1 = prim_match_end1->bs_strand;
				 read_idx1 = get_read_control_index(prim_match_end1,control_sequences, bisulfite_mode);
			} else if(vector_match_trace_used_end2) {
				 mstats->unmapped[0]++;
				 match_trace_t* prim_match_end2 = matches_get_primary_match(matches_end2);
				 update_distance_counts(matches_end2, prim_match_end2, mstats, 1);
				 bs2 = prim_match_end2->bs_strand;
				 read_idx2 = get_read_control_index(prim_match_end2,control_sequences, bisulfite_mode);
			} else {
				 mstats->unmapped[0]++;
				 mstats->unmapped[1]++;
			}
			mstats->hist_mapq[0]++;
	 } else {
			// We just look at primary alignments
			const paired_map_t * paired_map = paired_matches_get_primary_map(paired_matches);
			if(paired_map->pair_relation == pair_relation_concordant) { // Only collect template length stats for concordant pairs
				 mstats->correct_pairs++;
				 int64_t tlen=paired_map->template_length;
 				 uint64_t* count;
				 ihash_element_t* ih = ihash_get_ihash_element(mstats->insert_size_dist,tlen);
				 if(ih == NULL) {
						count=mm_alloc(uint64_t);
						*count=1;
						ihash_insert_element(mstats->insert_size_dist,tlen,count);
				 } else {
						count = ih->element;
						(*count)++;
				 }
			}
			match_trace_t* const match_end1 = paired_map->match_trace_end1;
			match_trace_t* const match_end2 = paired_map->match_trace_end2;
			bs1 = match_end1->bs_strand;
			bs2 = match_end2->bs_strand;
			update_distance_counts(paired_matches->matches_end1, match_end1, mstats, 0);
			update_distance_counts(paired_matches->matches_end2, match_end2, mstats, 1);

			read_idx1 = read_idx2 = get_read_control_index(match_end1,control_sequences, bisulfite_mode);
			// Get read type from index
			mstats->hist_mapq[(int)paired_map->mapq_score]++;
			if(paired_map->mapq_score>=min_mapq && paired_map->pair_relation == pair_relation_concordant) {
				update_conversion_counts(archive_search1, matches_end1, match_end1, mstats, read_idx1, 0);
				update_conversion_counts(archive_search2, matches_end2, match_end2, mstats, read_idx2, 1);
			}
	 }
	 if(bs1 == bs_strand_C2T) mstats->BSreads[0][0]++;
	 else if(bs1 == bs_strand_G2A) mstats->BSreads[0][1]++;
	 if(bs2 == bs_strand_C2T) mstats->BSreads[1][0]++;
	 else if(bs2 == bs_strand_G2A) mstats->BSreads[1][1]++;
	 if(read_idx1>=0) mapping_stats_reads_inc(mstats, 0, read_idx1);
	 if(read_idx2>=0) mapping_stats_reads_inc(mstats, 1, read_idx2);
}
char *indent_str="\t\t\t\t\t\t\t\t";
#define MAX_JSON_ARRAY_LINE 8
void output_json_uint_element(
    FILE *fp,
    char *key,
    uint64_t value,
    int indent,
    bool last) {
	 fprintf(fp,"%.*s\"%s\": %"PRIu64"%s",indent,indent_str,key,value,last?"\n":",\n");
}
void output_json_uint_array(
    FILE *fp,
    char *key,
    uint64_t const *values,
    int n,
    int indent,
    bool last) {
	 if(n<=MAX_JSON_ARRAY_LINE) {
			fprintf(fp,"%.*s\"%s\": ",indent,indent_str,key);
			int i;
			for(i=0;i<n;i++) fprintf(fp,"%s%" PRIu64,i?", ":"[",values[i]);
			fputs(last?"]\n":"],\n",fp);
	 } else {
			fprintf(fp,"%.*s\"%s\": [\n",indent,indent_str,key);
			int i;
			for(i=0;i<n;i+=MAX_JSON_ARRAY_LINE) {
				 fprintf(fp,"%.*s",indent+1,indent_str);
				 int j;
				 for(j=0;j<MAX_JSON_ARRAY_LINE;j++) {
						if(i+j<n-1) fprintf(fp,"%s%" PRIu64",",j?" ":"",values[i+j]);
						else {
							 fprintf(fp,"%s%" PRIu64,j?" ":"",values[i+j]);
							 break;
						}
				 }
				 fputc('\n',fp);
			}
			fprintf(fp,"%.*s%s",indent,indent_str,last?"]\n":"],\n");
	 }
}

const char *ct_desc[N_BASE_COUNTS] = {"N", "A", "C", "G", "T", "non_conv_C", "conv_C", "non_conv_CG", "conv_CG"};

void _output_base_counts_pe(FILE *fp, base_counts_t const *ct, int indent, int n) {
    bool first=true;
    for(int k=0;k<n;k++) {
        if(ct->counts[0][k]+ct->counts[1][k]>0) {
            fprintf(fp,"%s%.*s\"%s\": [%" PRIu64", %" PRIu64"]",first?"":",\n",indent,indent_str,ct_desc[k],ct->counts[0][k],ct->counts[1][k]);
            first=false;
        }   
	}
    fprintf(fp,"\n");
}

int get_n_counts(bs_strand_t bs) {
    if(bs==bs_strand_C2T || bs==bs_strand_G2A) return N_BASE_COUNTS;
    return N_REDUCED_BASE_COUNTS;
}

void output_base_counts_pe(FILE *fp, base_counts_t const *ct, int indent, bs_strand_t bs) {
    int n = get_n_counts(bs);
    // const char* (*cp)[N_BASE_COUNTS]=get_bs_strand_desc(bs, &n);
    _output_base_counts_pe(fp, ct, indent, n);
}

void _output_base_counts_se(FILE *fp, base_counts_t const *ct, int indent, int n) {
    bool first=true;
    for(int k=0;k<n;k++) {
        if(ct->counts[0][k]>0) {
            fprintf(fp,"%s%.*s\"%s\": [%" PRIu64"]",first?"":",\n",indent,indent_str,ct_desc[k],ct->counts[0][k]);
            first=false;
        }
	}
    fprintf(fp,"\n");
}

void output_base_counts_se(FILE *fp, base_counts_t const *ct, int indent, bs_strand_t bs) {
    int n=get_n_counts(bs);
    _output_base_counts_se(fp, ct, indent, n);
}

char * read_type[4] = {"SequencingControl", "UnderConversionControl", "OverConversionControl", "ConversionControl"};
void output_read_counts_pe(FILE *fp, mapper_parameters_t const * parameters, mapping_stats_t const * mstats, int i, int indent) {
    uint64_t n1, n2;
    if(i<0) {
        n1 = mstats->unmapped[0];
        n2 = mstats->unmapped[1];
    } else {
        n1 = *mapping_stats_reads_get(mstats, 0, i);
        n2 = *mapping_stats_reads_get(mstats, 1, i);
    }
	if(i<=0 || (n1 + n2) > 0) {
	    char *s = NULL;
		control_sequence_t *cs = NULL;
		if(i<0) {
		    s = "Unmapped";
		} else if(!i) {
		    s = "General";
		} else {
		    cs = vector_get_elm(parameters->search_parameters.control_sequences, i - 1, control_sequence_t);
			s = string_get_buffer(&cs->sequence_name);
		}
		fprintf(fp,"%.*s\"%s\": [%" PRIu64", %" PRIu64"]%s\n",indent,indent_str,s,n1,n2,i<0?"":",");
	}
}

void output_read_counts_se(FILE *fp, mapper_parameters_t const * parameters, mapping_stats_t const * mstats, int i, int indent) {
    uint64_t n;
    if(i<0) {
        n = mstats->unmapped[0];
    } else {
        n = *mapping_stats_reads_get(mstats, 0, i);
    }

	if(i<=0 || n > 0) {
	    char *s = NULL;
		control_sequence_t *cs = NULL;
		if(i<0) {
		    s = "Unmapped";
		} else if(!i) {
		    s = "General";
		} else {
		    cs = vector_get_elm(parameters->search_parameters.control_sequences, i - 1, control_sequence_t);
			s = string_get_buffer(&cs->sequence_name);
		}
		fprintf(fp,"%.*s\"%s\": [%" PRIu64"]%s\n",indent,indent_str,s,n,i<0?"":",");
	}
}

void output_mapping_stats(
    const mapper_parameters_t* parameters,
    const mapping_stats_t* mstats) {
	 char *output_file = parameters->io.report_file_name;
	 const bool bisulfite_index = (parameters->archive->type == archive_dna_bisulfite);
	 FILE *fp = fopen(output_file,"w");
	 if(!fp) return;
	 fputs("{\n",fp);
	 int indent=1;
	 int paired = 1;
	 char *mapper_type_s;
	 switch (parameters->mapper_type) {
		case mapper_se:
		  mapper_type_s="Single";
			paired = 0;
			break;
		case mapper_pe:
		  mapper_type_s="Paired";
			break;
		default:
		  mapper_type_s="Unknown";
			break;
	 }
	 fprintf(fp,"%.*s\"MapperType\": \"%s\",\n",indent,indent_str,mapper_type_s);
	 if(parameters->io.sam_parameters.read_group_header) {
			char *p = parameters->io.sam_parameters.read_group_header;
			size_t l = strlen(p) * 2 + 1;
			char *q = mm_malloc(l);
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
			fprintf(fp,"%.*s\"ReadGroup\": \"%s\",\n",indent,indent_str,q);
			mm_free(q);
	 }
	 int nc = mstats->n_control_seq;
	 vector_t const *control_sequences=parameters->search_parameters.control_sequences;
	 fprintf(fp,"%.*s\"ControlSequences\": {\n",indent++,indent_str);
	 bool first = true;
	 int jmax = bisulfite_index ? 4 : 1;
	 for(int j=0;j<jmax;j++) {
        for(int i=0;i<nc;i++) {
            control_sequence_t *cs = vector_get_elm(control_sequences, i, control_sequence_t);
            if(cs->sequence_type==j) {
                if(!first) fprintf(fp,",\n");
                else first=false;
                fprintf(fp,"%.*s\"%s\": {\n",indent++,indent_str, string_get_buffer(&cs->sequence_name));
                fprintf(fp,"%.*s\"type\": \"%s\"",indent,indent_str, read_type[cs->sequence_type]);
                if(cs->alt_name) {
                }
                fprintf(fp,"\n%.*s}",--indent,indent_str);
            }
        }
     }
	 fprintf(fp,"\n%.*s},\n",--indent,indent_str);
	 uint64_t tot_BSreads=0;
	 char* conv_type[]={"C2T","G2A"};
	 int i,j,k;
	 fprintf(fp,"%.*s\"Reads\": {\n",indent++,indent_str);

	 if(paired) {
		    // General reads
            output_read_counts_pe(fp, parameters, mstats, 0, indent);
            int jmax = bisulfite_index ? 4 : 1;
            for(j=0;j<jmax;j++) {
                for(i=0;i<nc;i++) {
                    control_sequence_t *cs = vector_get_elm(control_sequences, i, control_sequence_t);
                    if(cs->sequence_type==j) {
                        output_read_counts_pe(fp, parameters, mstats, i+1, indent);
                    }
                }
            }
            // Unmapped reads
			output_read_counts_pe(fp, parameters, mstats, -1, indent);
			fprintf(fp,"%.*s},\n",--indent,indent_str);

			for(i=0;i<2;i++) tot_BSreads+=mstats->BSreads[0][i]+mstats->BSreads[1][i];
			if(tot_BSreads) {
				 fprintf(fp,"%.*s\"NumReadsBS\": {\n",indent++,indent_str);
				 for(i=0;i<2;i++) {
						fprintf(fp,"%.*s\"%s\": [%"PRIu64", %"PRIu64"]%s",indent,indent_str,conv_type[i],mstats->BSreads[0][i],mstats->BSreads[1][i],i?"\n":",\n");
				 }
				 fprintf(fp,"%.*s},\n",--indent,indent_str);
			}
			output_json_uint_element(fp,"CorrectPairs",mstats->correct_pairs,indent,false);
			fprintf(fp,"%.*s\"BaseCounts\": {\n",indent++,indent_str);

			fprintf(fp,"%.*s\"Overall\": {\n",indent,indent_str);
			output_base_counts_pe(fp, &mstats->overall_counts, indent+1, bs_strand_mixed);
			fprintf(fp,"%.*s}",indent,indent_str);

			for(i=0;i<=nc;i++) {
				 base_counts_t const *c2t = mstats->base_counts[0] + i;
				 base_counts_t const *g2a = mstats->base_counts[1] + i;
				 uint64_t tot=0;
				 for(k=0;k<5;k++) {
					tot+=c2t->counts[0][k]+c2t->counts[1][k]+g2a->counts[0][k]+g2a->counts[1][k];
				 }

				 if(!i || tot>0) {
				    char const *s = NULL;
					if(i==0) {
					    s="General";
				    } else {
						s=string_get_buffer(&(vector_get_elm(control_sequences, i-1, control_sequence_t)->sequence_name));
					}
					fprintf(fp,",\n%.*s\"%s_C2T\": {\n",indent,indent_str,s);
					output_base_counts_pe(fp, c2t, indent+1, bs_strand_C2T);
					fprintf(fp,"%.*s}",indent,indent_str);
					fprintf(fp,",\n%.*s\"%s_G2A\": {\n",indent,indent_str,s);
					output_base_counts_pe(fp, g2a, indent+1, bs_strand_G2A);
					fprintf(fp,"%.*s}",indent,indent_str);
				 }
			}
			fprintf(fp,"\n%.*s},\n",--indent,indent_str);
			for(i=255;i>0;i--) if(mstats->hist_mapq[i]) break;
			output_json_uint_array(fp,"HistMapq",mstats->hist_mapq,i+1,indent,false);
			fprintf(fp,"%.*s\"HistReadLen\": [\n",indent++,indent_str);
			for(i=0;i<2;i++) {
				 fprintf(fp,"%.*s{\n",indent++,indent_str);
				 ihash_sort_by_key(mstats->read_length_dist[i]);
				 ihash_element_t* ih;
				 for(ih=mstats->read_length_dist[i]->head;ih;ih=ih->hh.next) {
						fprintf(fp,"%.*s\"%"PRId64"\": %"PRIu64"%s",indent,indent_str,ih->key,*((uint64_t *)ih->element),ih->hh.next?",\n":"\n");
				 }
				 fprintf(fp,"%.*s}%s",--indent,indent_str,i?"\n":",\n");
			}
			fprintf(fp,"%.*s],\n",--indent,indent_str);
			fprintf(fp,"%.*s\"HistMismatch\": [\n",indent++,indent_str);
			for(i=0;i<2;i++) {
				 fprintf(fp,"%.*s{\n",indent++,indent_str);
				 ihash_sort_by_key(mstats->distance_dist[i]);
				 ihash_element_t* ih;
				 for(ih=mstats->distance_dist[i]->head;ih;ih=ih->hh.next) {
						fprintf(fp,"%.*s\"%"PRId64"\": %"PRIu64"%s",indent,indent_str,ih->key,*((uint64_t *)ih->element),ih->hh.next?",\n":"\n");
				 }
				 fprintf(fp,"%.*s}%s",--indent,indent_str,i?"\n":",\n");
			}
			fprintf(fp,"%.*s],\n",--indent,indent_str);
			ihash_sort_by_key(mstats->insert_size_dist);
			fprintf(fp,"%.*s\"HistTemplateLen\": {\n",indent++,indent_str);
			ihash_element_t* ih;
			for(ih=mstats->insert_size_dist->head;ih;ih=ih->hh.next) {
				fprintf(fp,"%.*s\"%"PRId64"\": %"PRIu64"%s",indent,indent_str,ih->key,*((uint64_t *)ih->element),ih->hh.next?",\n":"\n");
			}
			fprintf(fp,"%.*s}\n",--indent,indent_str);
	 } else {
		    output_read_counts_se(fp, parameters, mstats, 0, indent);
			int jmax = bisulfite_index ? 4 : 1;
            for(j=0;j<jmax;j++) {
                for(i=0;i<nc;i++) {
                    control_sequence_t *cs = vector_get_elm(parameters->search_parameters.control_sequences, i, control_sequence_t);
                    if(cs->sequence_type==j) {
                        output_read_counts_se(fp, parameters, mstats, i+1, indent);
                    }
                }
            }
            output_read_counts_se(fp, parameters, mstats, -1, indent);
			fprintf(fp,"%.*s},\n",--indent,indent_str);
			for(i=0;i<2;i++) tot_BSreads+=mstats->BSreads[0][i];
			if(tot_BSreads) {
				 fprintf(fp,"%.*s\"NumReadsBS\": {\n",indent++,indent_str);
				 for(i=0;i<2;i++) {
						fprintf(fp,"%.*s\"%s\": [%"PRIu64"]%s",indent,indent_str,conv_type[i],mstats->BSreads[0][i],i?"\n":",\n");
				 }
				 fprintf(fp,"%.*s},\n",--indent,indent_str);
			}
			fprintf(fp,"%.*s\"BaseCounts\": {\n",indent++,indent_str);

			fprintf(fp,"%.*s\"Overall\": {\n",indent,indent_str);
			output_base_counts_se(fp, &mstats->overall_counts, indent+1, bs_strand_mixed);
			fprintf(fp,"%.*s}",indent,indent_str);

			for(i=0;i<=nc;i++) {
				 base_counts_t const *c2t = mstats->base_counts[0] + i;
				 base_counts_t const *g2a = mstats->base_counts[1] + i;
				 uint64_t tot=0;
				 for(k=0;k<5;k++) {
					tot+=c2t->counts[0][k]+g2a->counts[0][k];
				 }
				 if(!i || tot>0) {
					char const *s = NULL;
				    if(i==0) {
						s="General";
					} else {
					    s=string_get_buffer(&(vector_get_elm(control_sequences, i-1, control_sequence_t)->sequence_name));
				    }
					fprintf(fp,",\n%.*s\"%sC2T\": {\n",indent,indent_str,s);
					output_base_counts_se(fp, c2t, indent+1, bs_strand_C2T);
					fprintf(fp,"%.*s}",indent,indent_str);
					fprintf(fp,",\n%.*s\"%sG2A\": {\n",indent,indent_str,s);
					output_base_counts_se(fp, g2a, indent+1, bs_strand_G2A);
					fprintf(fp,"%.*s}",indent,indent_str);
				 }
			}
			fprintf(fp,"%.*s},\n",--indent,indent_str);
			for(i=255;i>0;i--) if(mstats->hist_mapq[i]) break;
			output_json_uint_array(fp,"HistMapq",mstats->hist_mapq,i+1,indent,false);
			fprintf(fp,"%.*s\"HistReadLen\": [\n",indent++,indent_str);
			fprintf(fp,"%.*s{\n",indent++,indent_str);
			ihash_sort_by_key(mstats->read_length_dist[0]);
			ihash_element_t* ih;
			for(ih=mstats->read_length_dist[0]->head;ih;ih=ih->hh.next) {
				 fprintf(fp,"%.*s\"%"PRId64"\": %"PRIu64"%s",indent,indent_str,ih->key,*((uint64_t *)ih->element),ih->hh.next?",\n":"\n");
			}
			fprintf(fp,"%.*s}\n",--indent,indent_str);
			fprintf(fp,"%.*s],\n",--indent,indent_str);
			fprintf(fp,"%.*s\"HistMismatch\": [\n",indent++,indent_str);
			fprintf(fp,"%.*s{\n",indent++,indent_str);
			ihash_sort_by_key(mstats->distance_dist[0]);
			for(ih=mstats->distance_dist[0]->head;ih;ih=ih->hh.next) {
				 fprintf(fp,"%.*s\"%"PRId64"\": %"PRIu64"%s",indent,indent_str,ih->key,*((uint64_t *)ih->element),ih->hh.next?",\n":"\n");
			}
			fprintf(fp,"%.*s}\n",--indent,indent_str);
			fprintf(fp,"%.*s]\n",--indent,indent_str);
	 }
	 fputs("}\n",fp);
	 fclose(fp);
}
