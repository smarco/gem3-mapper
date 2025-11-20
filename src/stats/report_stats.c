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
#include "archive/search/archive_search_se_parameters.h"
#include "matches/matches.h"
#include "matches/matches_cigar.h"
#include "stats/report_stats_mstats.h"
#include "text/dna_text.h"
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

const int cnv_idx_c2t[25]={0,0,0,0,0, 1,1,1,1,1, 0,6,6,5,6, 3,3,3,3,3, 0,8,8,7,8};
const int cnv_idx_g2a[25]={0,0,2,0,4, 0,8,2,6,4, 0,7,2,5,4, 0,8,2,6,4, 0,8,2,6,4};

void update_conversion_counts(
    const archive_search_t *archive,
    const matches_t* const matches,
    const match_trace_t* const match,
    const mapping_stats_t* mstats,
    const bs_strand_t bs,
	const strand_t strand,
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
	 if(bs==1 || bs==2) {
        sequence_t * seq_read = archive->sequence;
        const uint64_t cigar_length = match->match_alignment.cigar_length;
        const cigar_element_t* const cigar_array=vector_get_mem(matches->cigar_vector,cigar_element_t) + match->match_alignment.cigar_offset;
   		const search_parameters_t *params = &archive->search_parameters;
        const int (*cnv_idx)[25] = bs==1?&cnv_idx_c2t:&cnv_idx_g2a;

        int clip = params->conversion_clip_start;
        string_t *read = &seq_read->read;
		const char *p = string_get_buffer(read);
		uint64_t *ct = mapping_stats_base_counts(mstats, bs - 1, read_idx)->counts[end];
		int initial_skip = cigar_array[0].type==cigar_del?cigar_array[0].length : 0;
		int final_skip = cigar_array[cigar_length-1].type==cigar_del?cigar_array[cigar_length-1].length : 0;
		if(clip < initial_skip) clip = initial_skip;

		/* We don't check that l is not negative, but since clip is >=0, if it is negative it
		 * will be detected in the next line.
		 */
		int l = read->length -= final_skip;
		if(l<=clip) return;
		int prev=0;
		if(strand==Forward) {
		    if(seq_read->has_qualities) {
		        const uint8_t min_base_qual = params->conversion_min_base_qual + 33;
				const char *q = string_get_buffer(&seq_read->qualities);
				for(int i=clip; i < l; i++) {
				    const int b=q[i]>=min_base_qual?btab_fwd[(int)p[i]]:0;
					const int k=(*cnv_idx)[prev*5+b];
					ct[k]++;
					prev=b;
				}
			} else {
		        for(int i=clip; i < l; i++) {
					const int b = btab_fwd[(int)p[i]];
					const int k=(*cnv_idx)[prev*5+b];
					ct[k]++;
					prev=b;
				}
			}
		} else {
		    if(seq_read->has_qualities) {
				const uint8_t min_base_qual = params->conversion_min_base_qual + 33;
				const char *q = string_get_buffer(&seq_read->qualities);
				for(int i=l-1; i >=clip; i--) {
                    const int b=q[i]>=min_base_qual?btab_rev[(int)p[i]]:0;
                    const int k=(*cnv_idx)[prev*5+b];
                    ct[k]++;
                    prev=b;
				}
			} else {
				for(int i=l-1; i >=clip; i--) {
                    const int b=btab_rev[(int)p[i]];
                    const int k=(*cnv_idx)[prev*5+b];
                    ct[k]++;
                    prev=b;
				}
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
	 bs_strand_t bs = bs_strand_none;
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
			bs = match->bs_strand;
			read_idx = get_read_control_index(match, archive_search->search_parameters.control_sequences, bisulfite_mode);
			if(match->mapq_score>=min_mapq) {
				 update_conversion_counts(archive_search, matches, match, mstats, bs, match->strand, read_idx, 0);
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
			update_distance_counts(paired_matches->matches_end1, match_end1, mstats, 0);
			update_distance_counts(paired_matches->matches_end2, match_end2, mstats, 1);
			bs1 = match_end1 -> bs_strand;
			bs2 = match_end2 -> bs_strand;

			read_idx1 = read_idx2 = get_read_control_index(match_end1,control_sequences, bisulfite_mode);
			// Get read type from index
			mstats->hist_mapq[(int)paired_map->mapq_score]++;
			if(paired_map->mapq_score>=min_mapq && paired_map->pair_relation == pair_relation_concordant) {
				 update_conversion_counts(archive_search1, matches_end1, paired_map->match_trace_end1, mstats, bs1, match_end1->strand, read_idx1, 0);
				 update_conversion_counts(archive_search2, matches_end2, paired_map->match_trace_end2, mstats, bs2, match_end2->strand, read_idx2, 1);
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

const char *ct_desc_c2t[N_BASE_COUNTS] = {"N", "A", "C", "G", "T", "CG", "CH", "TG", "TH"};
const char *ct_desc_g2a[N_BASE_COUNTS] = {"N", "A", "C", "G", "T", "CG", "DG", "CA", "DA"};

void _output_base_counts_pe(FILE *fp, base_counts_t const *ct, int indent, int n, const char* (*ct_desc)[N_BASE_COUNTS]) {
    for(int k=0;k<n;k++) {
        if(ct->counts[0][k]+ct->counts[1][k]>0)
            fprintf(fp,"%.*s\"%s\": [%" PRIu64", %" PRIu64"]%s",indent,indent_str,(*ct_desc)[k],ct->counts[0][k],ct->counts[1][k],k==n-1?"\n":",\n");
	}
}

const char* (* get_bs_strand_desc(bs_strand_t bs, int *n))[N_BASE_COUNTS]  {
    const char* (*cp)[N_BASE_COUNTS];
    *n=N_BASE_COUNTS;
    if(bs==bs_strand_C2T) cp=&ct_desc_c2t;
    else if(bs==bs_strand_G2A) cp=&ct_desc_g2a;
    else {
        cp=&ct_desc_c2t;
        *n=N_REDUCED_BASE_COUNTS;
    }
    return cp;
}

void output_base_counts_pe(FILE *fp, base_counts_t const *ct, int indent, bs_strand_t bs) {
    int n;
    const char* (*cp)[N_BASE_COUNTS]=get_bs_strand_desc(bs, &n);
    _output_base_counts_pe(fp, ct, indent, n, cp);
}

void _output_base_counts_se(FILE *fp, base_counts_t const *ct, int indent, int n, const char* (*ct_desc)[N_BASE_COUNTS]) {
    for(int k=0;k<n;k++) {
        if(ct->counts[0][k]>0)
            fprintf(fp,"%.*s\"%s\": [%" PRIu64"]%s",indent,indent_str,(*ct_desc)[k],ct->counts[0][k],k==n-1?"\n":",\n");
	}
}

void output_base_counts_se(FILE *fp, base_counts_t const *ct, int indent, bs_strand_t bs) {
    int n;
    const char* (*cp)[N_BASE_COUNTS]=get_bs_strand_desc(bs, &n);
    _output_base_counts_se(fp, ct, indent, n, cp);
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
