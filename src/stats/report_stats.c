/*
 * PROJECT: GEMMapper
 * FILE: report_stats.c
 * DATE: 06/06/2012
 */

#include "stats/report_stats.h"
#include "matches/matches.h"

/*
 * Initialize mapping stats
 */
void init_mapping_stats(mapping_stats_t* mstats) {
	 int i,read;
	 for(read=0;read<2;read++) {
			mstats->unmapped[read]=0;
			mstats->read_length_dist[read]=ihash_new();
			int j;
			for(j=0;j<7;j++) for(i=0;i<5;i++) mstats->base_counts[j][read][i]=0;
			for(j=0;j<4;j++) mstats->reads[read][j]=0;
			for(j=0;j<2;j++) mstats->BSreads[read][j]=0;
	 }
	 for(i=0;i<256;i++) mstats->hist_mapq[i]=0;
	 mstats->correct_pairs=0;
	 mstats->insert_size_dist=ihash_new();
}

void merge_ihash(ihash_t* ihash1, ihash_t* ihash2) {
	 ihash_element_t* ih;
	 for(ih=ihash2->head;ih;ih=ih->hh.next) {
			uint64_t* count;
			ihash_element_t* ih1 = ihash_get_ihash_element(ihash1,ih->key);
			if(ih1 == NULL) {
				 count = mm_alloc(uint64_t);
				 *count = *((uint64_t*)ih->element);
				 ihash_insert_element(ihash1,ih->key,count);
			} else {
				 count = ih1->element;
				 *count += *((uint64_t*)ih->element);
			}
	 }
}

void merge_mapping_stats(mapping_stats_t* global_mstats, mapping_stats_t* mstats, uint64_t num_threads) {
	 init_mapping_stats(global_mstats);
	 uint64_t i;
	 for(i=0;i<num_threads;i++) {
			global_mstats->correct_pairs+=mstats[i].correct_pairs;
			int rd,j;
			for(rd=0;rd<2;rd++) {
				 for(j=0;j<4;j++) global_mstats->reads[rd][j]+=mstats[i].reads[rd][j];
				 for(j=0;j<2;j++) global_mstats->BSreads[rd][j]+=mstats[i].BSreads[rd][j];
				 global_mstats->unmapped[rd]+=mstats[i].unmapped[rd];
				 int k;
				 for(k=0;k<7;k++) for(j=0;j<5;j++) global_mstats->base_counts[k][rd][j]+=mstats[i].base_counts[k][rd][j];
			}
			for(j=0;j<256;j++) global_mstats->hist_mapq[j]+=mstats[i].hist_mapq[j];
			for(j=0;j<2;j++) merge_ihash(global_mstats->read_length_dist[j],mstats[i].read_length_dist[j]);
			merge_ihash(global_mstats->insert_size_dist,mstats[i].insert_size_dist);
	 }
}

int btab[256]={ ['A'] = 1, ['C'] = 2, ['G'] = 3, ['T'] = 4 };

void update_counts(sequence_t* const seq_read, mapping_stats_t* mstats,int end) {
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
	 while(*p) mstats->base_counts[0][end][btab[(int)*p++]]++;
}

void update_conversion_counts(sequence_t* const seq_read, mapping_stats_t* mstats,int end,bs_strand_t bs,int read_type) {
	 int cnv_idx[4][4] = {{0,0,0,0},{1,0,3,5},{2,0,4,6},{0,0,0,0}};
	 
	 int idx=cnv_idx[bs][read_type];
	 if(idx) {
			string_t* const read = &seq_read->read;
			char *p = string_get_buffer(read);
			while(*p) mstats->base_counts[idx][end][btab[(int)*p++]]++;
	 }
}

int get_read_type(match_trace_t* match) {
	 char *control_seqs[]={SEQUENCING_CONTROL, UNDERCONVERSION_CONTROL, OVERCONVERSION_CONTROL,0};
	 char* seq = match->sequence_name;
	 char* p;
	 int i = 0;
	 while((p = control_seqs[i++])) if(strstr(seq, p)) break;
	 return i%4;
}

void collect_SE_mapping_stats(archive_search_t* const archive_search, matches_t* const matches, mapping_stats_t* mstats) {
	 update_counts(&archive_search->sequence,mstats,0);
	 bs_strand_t bs = bs_strand_none;
	 int read_type = -1;

	 const uint64_t vector_match_trace_used = vector_get_used(matches->position_matches);
	 if (gem_expect_false(vector_match_trace_used==0)) { // Unmapped
			mstats->unmapped[0]++;
	 } else {
			// We just look at primary alignments
			match_trace_t* match = vector_get_mem(matches->position_matches,match_trace_t);
			bs = match->bs_strand;
			read_type = get_read_type(match);
			if(match->mapq_score>0) {
				 update_conversion_counts(&archive_search->sequence, mstats, 0, bs, read_type);
			}
			mstats->hist_mapq[(int)match->mapq_score]++;
			if(bs == bs_strand_C2T) mstats->BSreads[0][0]++;
			else if(bs == bs_strand_G2A) mstats->BSreads[0][1]++;
			if(read_type>=0) mstats->reads[0][read_type]++;
	 }
}

void collect_PE_mapping_stats(archive_search_t* const archive_search1, archive_search_t* const archive_search2,
 	 paired_matches_t* const paired_matches, mapping_stats_t* mstats) {

	 update_counts(&archive_search1->sequence,mstats,0);
	 update_counts(&archive_search2->sequence,mstats,1);
	 matches_t* const matches_end1 = paired_matches->matches_end1;
	 matches_t* const matches_end2 = paired_matches->matches_end2;	 
	 bs_strand_t bs1,bs2;
	 bs1 = bs2 = bs_strand_none;
	 int read_type1 = -1, read_type2 = -1;
	 
	 if (gem_expect_false(!paired_matches_is_mapped(paired_matches))) { // Non paired
			const uint64_t vector_match_trace_used_end1 = vector_get_used(matches_end1->position_matches);
			const uint64_t vector_match_trace_used_end2 = vector_get_used(matches_end2->position_matches);
			if(vector_match_trace_used_end1 && vector_match_trace_used_end2) {
 				 match_trace_t* prim_match_end1 = vector_get_mem(matches_end1->position_matches,match_trace_t);
				 match_trace_t* prim_match_end2 = vector_get_mem(matches_end2->position_matches,match_trace_t);
				 bs1 = prim_match_end1 -> bs_strand;
				 bs2 = prim_match_end2 -> bs_strand;
				 read_type1 = get_read_type(prim_match_end1);
				 read_type2 = get_read_type(prim_match_end2);
			} else if(vector_match_trace_used_end1) {
				 mstats->unmapped[1]++;
				 match_trace_t* prim_match_end1 = vector_get_mem(matches_end1->position_matches,match_trace_t);
				 bs1 = prim_match_end1->bs_strand;
				 read_type1 = get_read_type(prim_match_end1);
			} else if(vector_match_trace_used_end2) {
				 mstats->unmapped[0]++;
				 match_trace_t* prim_match_end2 = vector_get_mem(matches_end2->position_matches,match_trace_t);
				 bs2 = prim_match_end2->bs_strand;
				 read_type2 = get_read_type(prim_match_end2);
			} else {
				 mstats->unmapped[0]++;
				 mstats->unmapped[1]++;
			}
			mstats->hist_mapq[0]++;
	 } else {
			// We just look at primary alignments
			paired_map_t* paired_map = vector_get_mem(paired_matches->paired_maps,paired_map_t);
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
			match_trace_t* const match_end1 = matches_get_match_trace(paired_matches->matches_end1,paired_map->match_end1_offset);
			match_trace_t* const match_end2 = matches_get_match_trace(paired_matches->matches_end2,paired_map->match_end2_offset);
			bs1 = match_end1 -> bs_strand;
			bs2 = match_end2 -> bs_strand;
			read_type1 = read_type2 = get_read_type(match_end1);
			mstats->hist_mapq[(int)paired_map->mapq_score]++;
			if(paired_map->mapq_score>0) {
				 update_conversion_counts(&archive_search1->sequence, mstats, 0, bs1, read_type1);
				 update_conversion_counts(&archive_search2->sequence, mstats, 1, bs2, read_type2);
			}
	 }
	 if(bs1 == bs_strand_C2T) mstats->BSreads[0][0]++;
	 else if(bs1 == bs_strand_G2A) mstats->BSreads[0][1]++;
	 if(bs2 == bs_strand_C2T) mstats->BSreads[1][0]++;
	 else if(bs2 == bs_strand_G2A) mstats->BSreads[1][1]++;
	 if(read_type1>=0) mstats->reads[0][read_type1]++;
	 if(read_type2>=0) mstats->reads[1][read_type2]++;
}

char *indent_str="\t\t\t\t\t\t\t\t";
#define MAX_JSON_ARRAY_LINE 8

void output_json_uint_element(FILE *fp,char *key,uint64_t value,int indent,bool last) {
	 fprintf(fp,"%.*s\"%s\": %"PRIu64"%s",indent,indent_str,key,value,last?"\n":",\n");
}

void output_json_uint_array(FILE *fp,char *key,uint64_t *values,int n,int indent,bool last) {
	 if(n<=MAX_JSON_ARRAY_LINE) {
			fprintf(fp,"%.*s\"%s\": ",indent,indent_str,key);
			int i;
			for(i=0;i<n;i++) fprintf(fp,"%s%"PRIu64,i?", ":"[",values[i]);
			fputs(last?"]\n":"],\n",fp);
	 } else {
			fprintf(fp,"%.*s\"%s\": [\n",indent,indent_str,key);
			int i;
			for(i=0;i<n;i+=MAX_JSON_ARRAY_LINE) {
				 fprintf(fp,"%.*s",indent+1,indent_str);
				 int j;
				 for(j=0;j<MAX_JSON_ARRAY_LINE;j++) {
						if(i+j<n-1) fprintf(fp,"%s%"PRIu64",",j?" ":"",values[i+j]);
						else {
							 fprintf(fp,"%s%"PRIu64,j?" ":"",values[i+j]);
							 break;
						}
				 }
				 fputc('\n',fp);
			}
			fprintf(fp,"%.*s%s",indent,indent_str,last?"]\n":"],\n");
	 }
}

void output_mapping_stats(mapper_parameters_t *parameters, mapping_stats_t* mstats) {
	 char *output_file = parameters->io.report_file_name;
	 FILE *fp = fopen(output_file,"w");
	 if(!fp) return;
	 fputs("{\n",fp);
	 int indent=1;
	 int paired = 1;
	 char *mapper_type;
	 switch (parameters->mapper_type) {
		case mapper_se:
			mapper_type="Single";
			paired = 0;
			break;
		case mapper_pe:
			mapper_type="Paired";
			break;
		default:
			mapper_type="Unknown";
			break;
	 }
	 fprintf(fp,"%.*s\"MapperType\": \"%s\",\n",indent,indent_str,mapper_type);
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
	 uint64_t tot_BSreads=0;
	 char* read_type[]={"General","SequencingControl","UnderConversionControl","OverConversionControl"};
	 char* conv_type[]={"C2T","G2A"};
	 int i,j,k;

	 fprintf(fp,"%.*s\"Reads\": {\n",indent++,indent_str);
	 if(paired) {
			for(i=0;i<4;i++) {
				 if(!i || mstats->reads[0][i] || mstats->reads[1][i])
					 fprintf(fp,"%.*s\"%s\": [%"PRIu64", %"PRIu64"],\n",indent,indent_str,read_type[i],mstats->reads[0][i],mstats->reads[1][i]);
			}
			fprintf(fp,"%.*s\"Unmapped\": [%"PRIu64", %"PRIu64"]\n",indent,indent_str,mstats->unmapped[0],mstats->unmapped[1]);
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
			char* bc_type[] = {"Overall","GeneralC2T","GeneralG2A","UnderConversionControlC2T","UnderConversionControlG2A","OverConversionControlC2T","OverConversionControlG2A"};
			char* base = "NACGT";
			uint64_t tot[6];
			for(i=1;i<7;i++) {
				 tot[i-1]=0;
				 for(j=0;j<2;j++) for(k=0;k<5;k++) tot[i-1]+=mstats->base_counts[i][j][k];
			}
			for(i=5;i>=0;i--) if(tot[i]) break;
			j=i+1;
			for(i=0;i<=j;i++) {
				 bool last = (i == j);
				 if(!i || tot[i-1]) {
						fprintf(fp,"%.*s\"%s\": {\n",indent++,indent_str,bc_type[i]);
						for(k=0;k<5;k++) {
							 int k1 = (k+1)%5;
							 fprintf(fp,"%.*s\"%c\": [%"PRIu64", %"PRIu64"]%s",indent,indent_str,base[k1],mstats->base_counts[i][0][k1],mstats->base_counts[i][1][k1],k==4?"\n":",\n");
						}
						fprintf(fp,"%.*s}%s",--indent,indent_str,last?"\n":",\n");
				 }
			}
			fprintf(fp,"%.*s},\n",--indent,indent_str);
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
			ihash_sort_by_key(mstats->insert_size_dist);
			fprintf(fp,"%.*s\"HistTemplateLen\": {\n",indent++,indent_str);
			ihash_element_t* ih;
			for(ih=mstats->insert_size_dist->head;ih;ih=ih->hh.next) {
				fprintf(fp,"%.*s\"%"PRId64"\": %"PRIu64"%s",indent,indent_str,ih->key,*((uint64_t *)ih->element),ih->hh.next?",\n":"\n");
			}
			fprintf(fp,"%.*s}\n",--indent,indent_str);
	 } else {
			for(i=0;i<4;i++) {
				 if(!i || mstats->reads[0][i])
					 fprintf(fp,"%.*s\"%s\": [%"PRIu64"],\n",indent,indent_str,read_type[i],mstats->reads[0][i]);
			}
			fprintf(fp,"%.*s\"Unmapped\": [%"PRIu64"]\n",indent,indent_str,mstats->unmapped[0]);
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
			char* bc_type[] = {"Overall","GeneralC2T","GeneralG2A","UnderConversionControlC2T","UnderConversionControlG2A","OverConversionControlC2T","OverConversionControlG2A"};
			char* base = "NACGT";
			uint64_t tot[6];
			for(i=1;i<7;i++) {
				 tot[i-1]=0;
				 for(k=0;k<5;k++) tot[i-1]+=mstats->base_counts[i][0][k];
			}
			for(i=5;i>=0;i--) if(tot[i]) break;
			j=i+1;
			for(i=0;i<=j;i++) {
				 bool last = (i == j);
				 if(!i || tot[i-1]) {
						fprintf(fp,"%.*s\"%s\": {\n",indent++,indent_str,bc_type[i]);
						for(k=0;k<5;k++) {
							 int k1 = (k+1)%5;
							 fprintf(fp,"%.*s\"%c\": [%"PRIu64"]%s",indent,indent_str,base[k1],mstats->base_counts[i][0][k1],k==4?"\n":",\n");
						}
						fprintf(fp,"%.*s}%s",--indent,indent_str,last?"\n":",\n");
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
	 }
	 fputs("}\n",fp);
	 fclose(fp);
}
