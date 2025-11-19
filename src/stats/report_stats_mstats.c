#include "stats/report_stats_mstats.h"

void clear_base_counts(base_counts_t *p) {
    for(int i=0; i < 2; i++) {
        for(int j = 0; j < N_BASE_COUNTS; j++) {
            p->counts[i][j]=0;
        }
    }
}

void add_base_counts(base_counts_t *p, base_counts_t const *p1) {
    for(int i=0; i < 2; i++) {
        for(int j = 0; j < N_BASE_COUNTS; j++) {
            p->counts[i][j]+=p1->counts[i][j];
        }
    }
}

/*
 * Initialize mapping stats
 */
void init_mapping_stats(
    mapping_stats_t* const mstats) {
	 int i,j,read;
	 int nc = mstats->n_control_seq;
	 memset(mstats->reads, 0, sizeof(uint64_t) * 2 * (nc + 1));
	 
	 clear_base_counts(&mstats->overall_counts);
	 for(i=0;i<2;i++) {
		for(j=0;j<=nc;j++) {
            clear_base_counts(mstats->base_counts[i] + j);
		}
	 }
	 for(read=0;read<2;read++) {
			mstats->unmapped[read]=0;
			mstats->read_length_dist[read]=ihash_new(NULL);
			mstats->distance_dist[read]=ihash_new(NULL);
			int j;
			for(j=0;j<2;j++) mstats->BSreads[read][j]=0;
	 }
	 for(i=0;i<256;i++) mstats->hist_mapq[i]=0;
	 mstats->correct_pairs=0;
	 mstats->insert_size_dist=ihash_new(NULL);
}

void merge_ihash(
    ihash_t* ihash1,
    ihash_t* ihash2) {
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

void merge_mapping_stats(
    mapping_stats_t* const global_mstats,
    mapping_stats_t* const mstats,
    const uint64_t num_threads) {
	 init_mapping_stats(global_mstats);
	 uint64_t i;
	 int nc = global_mstats->n_control_seq;
	 for(i=0;i<num_threads;i++) {
			global_mstats->correct_pairs+=mstats[i].correct_pairs;
			for(int j=0;j<2*(nc+1);j++) {
                global_mstats->reads[j] += mstats[i].reads[j]; 
			} 
			
			int rd,j;
			for(rd=0;rd<2;rd++) {
				 for(j=0;j<2;j++) global_mstats->BSreads[rd][j]+=mstats[i].BSreads[rd][j];
				 global_mstats->unmapped[rd]+=mstats[i].unmapped[rd];
				 add_base_counts(&global_mstats->overall_counts, &mstats[i].overall_counts);
				 int k;
				 for(k=0;k<=nc;k++) {
				 	add_base_counts(global_mstats->base_counts[0] + k, mstats[i].base_counts[0] + k);
					add_base_counts(global_mstats->base_counts[1] + k, mstats[i].base_counts[1] + k);
				 } 
			}
			for(j=0;j<256;j++) global_mstats->hist_mapq[j]+=mstats[i].hist_mapq[j];
			for(j=0;j<2;j++) {
				 merge_ihash(global_mstats->read_length_dist[j],mstats[i].read_length_dist[j]);
				 merge_ihash(global_mstats->distance_dist[j],mstats[i].distance_dist[j]);
			}
			merge_ihash(global_mstats->insert_size_dist,mstats[i].insert_size_dist);
	 }
}

void setup_mapping_stats(mapping_stats_t *ms, int n_control_seq) {
    ms->n_control_seq = n_control_seq;
    ms->reads = mm_calloc(2 * (n_control_seq + 1), uint64_t, false);
    ms->base_counts[0] = mm_calloc(n_control_seq + 1, base_counts_t, false);
    ms->base_counts[1] = mm_calloc(n_control_seq + 1, base_counts_t, false);
}

mapping_stats_t *new_mapping_stats(int n_control_seq) {
    mapping_stats_t *ms = mm_alloc(mapping_stats_t);
    setup_mapping_stats(ms, n_control_seq);
    return ms;
}

