/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
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
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 *            Simon Heath <simon.heath@gmail.com>
 * DESCRIPTION:
 *   Mapper module encapsulates and provides accessors to all
 *   the parameters used by the mapper
 */

#include "text/restriction_text.h"
#include "text/dna_text.h"

static uint8_t iupac[256] = {
  ['A'] = 1, ['C'] = 2, ['G'] = 4, ['T'] = 8,
  ['R'] = 5, ['Y'] = 10, ['S'] = 6, ['W'] = 9, ['K'] = 12, ['M'] = 3,
  ['B'] = 14, ['D'] = 13, ['H'] = 11, ['V'] = 7, ['N'] = 15,
  ['a'] = 1, ['c'] = 2, ['g'] = 4, ['t'] = 8,
  ['r'] = 5, ['y'] = 10, ['s'] = 6, ['w'] = 9, ['k'] = 12, ['m'] = 3,
  ['b'] = 14, ['d'] = 13, ['h'] = 11, ['v'] = 7, ['n'] = 15,
};

static uint8_t iupac_compl[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

static char *restriction_enzyme_defs[][2] = {
		{"HindIII", "A-AGCTT"},
		{"AluI", "AG-CT"},
		{"DpnI", "-GATC"},
		{"DpnII", "-GATC"},
		{"MboI", "-GATC"},
		{"Sau3AI", "-GATC"},
		{"MseI", "T-TAA"},
		{"HaeIII", "GG-CC"},
		{"NcoI", "C-CATGG"},
		{"TaqI", "T-CGA"},
		{"BfaI", "C-TAG"},
		{"HpyCH4V", "TG-CA"},
		{"MluCI", "-AATT"},
		{"ApeKI", "G-CWGC"},
		{"CviQI", "G-TAC"},
		{"CviAII", "C-ATG"},
		{"AseI", "AT-TAAT"},
		{NULL, NULL}
};

void restriction_make_search_patterns(restriction_site_t * const rest, const char * const pattern, const uint64_t search_pattern_len) {
//	fprintf(stderr,"Junction: %.*s\n", (int)search_pattern_len, pattern);
	uint64_t x = (uint64_t)1;
	for(uint64_t i = 0; i < search_pattern_len; i++) {
		uint8_t enc = iupac[(int)pattern[search_pattern_len - 1 - i]];
		if(enc & 1) { // A
			rest->search_pattern[bs_strand_none][ENC_DNA_CHAR_A] |= x;
			rest->search_pattern[bs_strand_C2T][ENC_DNA_CHAR_A] |= x;
			rest->search_pattern[bs_strand_G2A][ENC_DNA_CHAR_A] |= x;
		}
		if(enc & 2) { // C - for bs_strand_C2T, both C and T in the read should match against a C in the pattern
			rest->search_pattern[bs_strand_none][ENC_DNA_CHAR_C] |= x;
			rest->search_pattern[bs_strand_C2T][ENC_DNA_CHAR_C] |= x;
			rest->search_pattern[bs_strand_C2T][ENC_DNA_CHAR_T] |= x;
			rest->search_pattern[bs_strand_G2A][ENC_DNA_CHAR_C] |= x;
		}
		if(enc & 4) { // G - for bs_strand_G2A, both G and A in the read should match against a G in the pattern
			rest->search_pattern[bs_strand_none][ENC_DNA_CHAR_G] |= x;
			rest->search_pattern[bs_strand_C2T][ENC_DNA_CHAR_G] |= x;
			rest->search_pattern[bs_strand_G2A][ENC_DNA_CHAR_G] |= x;
			rest->search_pattern[bs_strand_G2A][ENC_DNA_CHAR_A] |= x;
		}
		if(enc & 8) { // T
			rest->search_pattern[bs_strand_none][ENC_DNA_CHAR_T] |= x;
			rest->search_pattern[bs_strand_C2T][ENC_DNA_CHAR_T] |= x;
			rest->search_pattern[bs_strand_G2A][ENC_DNA_CHAR_T] |= x;
		}
		x <<= 1;
	}
 	uint64_t mask = ((uint64_t)1 << search_pattern_len) - (uint64_t)1;
	for(int i = 0; i < DNA_EXT_RANGE; i++) {
		rest->search_pattern[bs_strand_none][i] ^= mask;
		rest->search_pattern[bs_strand_C2T][i] ^= mask;
		rest->search_pattern[bs_strand_G2A][i] ^= mask;
	}
	rest->search_pattern_mask = mask;
}

restriction_site_t *restriction_new(char * const rest_char) {
  restriction_site_t *rest = NULL;
  if(rest_char != NULL) {
	char *tp = rest_char;
    rest = mm_alloc(restriction_site_t);

    string_init(&rest->restriction_site, 0, NULL);
    rest->restriction_enzyme = NULL;

    // First check if the supplied string matches a known enzyme
    // If set tp to point to the stored restriction site for that enzyme
    int i = 0;
    while(restriction_enzyme_defs[i][0] != NULL) {
    	if(!strcasecmp(rest_char, restriction_enzyme_defs[i][0])) {
    		rest->restriction_enzyme = restriction_enzyme_defs[i][0];
    		tp = restriction_enzyme_defs[i][1];
    	}
    	i++;
    }
    char c;
    int index = -1;
    i = 0;
    while((c = tp[i++])) {
      if(iupac[(int)c]) string_append_char(&rest->restriction_site, c);
      else if(c == '.' || c == ':' || c == '_' || c == '-' || c == '|' || isspace((int)c)) {
        if(index == -1) index = i - 1;
        else {
          fprintf(stderr,"Multiple cut sites (%s) passed to restriction_new()\n", rest_char);
          break;
        }
      } else {
        fprintf(stderr,"Could not parse restriction site '%s' passed to restriction_new()\n", rest_char);
        break;
      }
    }
    if(c) {
      index = -1;
    } else if(index < 0) {
      fprintf(stderr,"No cut site found in string '%s' passed to restriction_new()\n", rest_char);    } else {
      int sz = string_get_length(&rest->restriction_site);
      if(sz == 0) {
        fprintf(stderr,"No sequence found in string '%s' passed to restriction_new()\n", rest_char);
        index = -1;
      } else {
        char const *p = string_get_buffer(&rest->restriction_site);
        for(int i = 0; i <= (sz >> 1); i++) {
          if(iupac[(int)p[i]] != iupac_compl[(int)iupac[(int)p[sz - 1 - i]]]) {
            fprintf(stderr,"Sequence '%s' passed to restriction_new is not palindromic()\n", rest_char);
            index = -1;
            break;
          }
        }
      }
    }
    if(index < 0) {
      restriction_site_delete(rest);
      rest = NULL;
    } else {
      rest->cut_site_index = index;
//      if(rest->restriction_enzyme != NULL) {
//      	fprintf(stderr,"Added restriction site: '%.*s' (%s), cut index = %"PRIu64"\n", PRIs_content(&rest->restriction_site), rest->restriction_enzyme, rest->cut_site_index);
//
//      } else {
//      	fprintf(stderr,"Added restriction site: '%.*s', cut index = %"PRIu64"\n", PRIs_content(&rest->restriction_site), rest->cut_site_index);
//      }
    }
    if(rest) {
    	rest->interaction_site = false;
    	// Generate search patterns
    	for(uint64_t i = 0; i < DNA_EXT_RANGE; i++) {
    		rest->reference_pattern[i] = (uint64_t)0;
    		for(uint64_t k = 0; k < 3; k++) {
    			rest->search_pattern[k][i] = (uint64_t)0;
    		}
    	}
    	const uint64_t ref_pattern_len = string_get_length(&rest->restriction_site);
    	const uint64_t search_pattern_len = (ref_pattern_len - rest->cut_site_index) << 1;
    	gem_cond_error_msg(search_pattern_len > 64,"Restriction site sequence too long (limit 64 bases)");
    	rest->search_pattern_len = search_pattern_len;
    	rest->reference_pattern_len = ref_pattern_len;
    	const char * const p = string_get_buffer(&rest->restriction_site);
    	// First we construct the pattern for the reference genome
     	uint64_t x = (uint64_t)1;
     	for(uint64_t i = 0; i < ref_pattern_len; i++) {
     		uint8_t enc = iupac[(int)p[ref_pattern_len - 1 - i]];
     		if(enc & 1) rest->reference_pattern[ENC_DNA_CHAR_A] |= x; // A
     		if(enc & 2) rest->reference_pattern[ENC_DNA_CHAR_C] |= x; // C
     		if(enc & 4) rest->reference_pattern[ENC_DNA_CHAR_G] |= x; // G
     		if(enc & 8) rest->reference_pattern[ENC_DNA_CHAR_T] |= x; // T
     		x <<= 1;
     	}
     	uint64_t mask = ((uint64_t)1 << ref_pattern_len) - (uint64_t)1;
     	for(int i = 0; i < DNA_EXT_RANGE; i++)
     		rest->reference_pattern[i] ^= mask;
     	rest->reference_pattern_mask = mask;

     	// And now for patterns for searching the reads for ligation junctions
     	// In the reads we should see the reverse strand sequence from the end of the site until the cut site fused to
     	// the forward strand sequence after the cut site.
     	// For example, with HindIII, the sequence is A|AGCTT, so in the reads we should see AAGCTAGCTT

     	char *const junction = mm_malloc(search_pattern_len);
     	uint64_t k = 0;
     	for(uint64_t i = 0; i < ref_pattern_len - rest->cut_site_index; i++) junction[k++] = p[i];
     	rest->search_pattern_split_offset = k;
     	for(uint64_t i = rest->cut_site_index; i < ref_pattern_len; i++) junction[k++] = p[i];
     	restriction_make_search_patterns(rest, junction, search_pattern_len);
     	mm_free(junction);
   }
  } else fprintf(stderr,"NULL restriction site description passed to restriction_new()\n");
  return rest;
}

/*
 * If we have multiple restriction sites in 3C / HiC experiments, we can have hybrid ligation junctions
 * We will create these as additional restriction_site_t structures with the interaction_site flag set
 */

void restriction_create_interaction_junctions(vector_t * const restriction_sites) {
	const uint64_t num_sites = vector_get_used(restriction_sites);
  restriction_site_t ** sites = vector_get_mem(restriction_sites, restriction_site_t *);
	for(uint64_t i = 1; i < num_sites; i++) {
		for(uint64_t j = 0; j < i; j++) {
			uint64_t len = sites[i]->reference_pattern_len - sites[i]->cut_site_index +
					sites[j]->reference_pattern_len - sites[j]->cut_site_index;
			char *const junction = mm_malloc(len);
			const char * const p1 = string_get_buffer(&sites[i]->restriction_site);
			const char * const p2 = string_get_buffer(&sites[j]->restriction_site);
			// First do junction site[i] : site[j]
			restriction_site_t *rest = mm_calloc(1, restriction_site_t, true);
			rest->interaction_site = true;
			rest->search_pattern_len = len;
			uint64_t l = 0;
			for(uint64_t k = 0; k < sites[i]->reference_pattern_len - sites[i]->cut_site_index; k++) junction[l++] = p1[k];
     	rest->search_pattern_split_offset = l;
			for(uint64_t k = sites[j]->cut_site_index; k < sites[j]->reference_pattern_len; k++) junction[l++] = p2[k];
     	restriction_make_search_patterns(rest, junction, len);
     	vector_insert(restriction_sites, rest, restriction_site_t *);
			// And then junction site[j] : site[i]
			rest = mm_calloc(1, restriction_site_t, true);
			rest->interaction_site = true;
			rest->search_pattern_len = len;
			l = 0;
			for(uint64_t k = 0; k < sites[j]->reference_pattern_len - sites[j]->cut_site_index; k++) junction[l++] = p2[k];
     	rest->search_pattern_split_offset = l;
			for(uint64_t k = sites[i]->cut_site_index; k < sites[i]->reference_pattern_len; k++) junction[l++] = p1[k];
     	restriction_make_search_patterns(rest, junction, len);
     	vector_insert(restriction_sites, rest, restriction_site_t *);
		}
	}
}

void restriction_site_delete(restriction_site_t * const rest) {
  if(rest != NULL) {
    string_destroy(&rest->restriction_site);
    mm_free(rest);
  }
}

/*
 * Find matches to ligation junctions within a read (allowing for bisulfite conviersion if necessary)
 */
void find_restriction_site_matches(
		pattern_t const * const pattern,
		vector_t const * const restriction_sites,
		vector_t * const restriction_hits,
		bisulfite_conversion_t const bisulfite_conversion) {

	vector_clear(restriction_hits);

//	uint64_t number_hits = 0;
	const uint64_t num_junctions = vector_get_used(restriction_sites);
	if(num_junctions == 0) return;
	restriction_site_t **junctions = vector_get_mem(restriction_sites, restriction_site_t *);
	const uint64_t len = pattern->key_length;
	uint64_t num_hits = 0;
	for(uint64_t j = 0; j < num_junctions; j++) {
		uint8_t const *key = pattern->key_non_bs;
		uint64_t const * const search_pattern = junctions[j]->search_pattern[bisulfite_conversion];
  	uint64_t D = junctions[j]->search_pattern_mask;
  	uint64_t const off = junctions[j]->search_pattern_split_offset;
		for(uint64_t i = 0; i < len; i++) {
			D = (D >> 1) | search_pattern[(int)*key++];
			if(!(D & 1)) { // Full length match found
				uint64_t hit = i - off;
				if(gem_expect_true(num_hits == 0)) {
					vector_insert(restriction_hits, hit, uint64_t);
					num_hits++;
				} else {
					uint64_t k;
					bool insert = true;
					for(k = 0; k < num_hits; k++) {

						uint64_t hit1 = *(vector_get_elm(restriction_hits, k, uint64_t));
						if(hit < hit1) {
							vector_reserve_additional(restriction_hits, 1);
							uint64_t * const hits = vector_get_mem(restriction_hits, uint64_t);
							memmove(hits + k + 1, hits + k,	sizeof(uint64_t) * (num_hits - k));
							vector_set_used(restriction_hits, num_hits + 1);
							break;
						} else if(hit == hit1) {
							insert = false;
							break;
						}
					}
					if(insert) {
						if(k == num_hits) {
							vector_insert(restriction_hits, hit, uint64_t);
						} else {
							vector_set_elm(restriction_hits, k, uint64_t, hit);
						}
						num_hits++;
					}
				}
			}
		}
	}
}

/*
 * Hands out next available sequence (contig) to be processed.
 * Returns NULL when nothing left to process
 */
restriction_site_sequence_locator_t * sequence_locator_sharer(
		restriction_site_locator_builder_t* const builder) {

	restriction_site_sequence_locator_t *sequence_locator = NULL;
	pthread_mutex_lock(&builder->mutex);
	if(builder->next_tag_idx < builder->num_tags)
		sequence_locator = builder->restriction_site_locator->restriction_site_sequence_locator + (builder->next_tag_idx++);
	pthread_mutex_unlock(&builder->mutex);
	return sequence_locator;
}

void *restriction_site_write(restriction_site_locator_builder_t* const builder) {
	const archive_t * const archive = builder->archive;
	const locator_t * const locator = archive->locator;
	const restriction_site_locator_t * const restriction_locator = builder->restriction_site_locator;
	FILE *site_output = builder->output_handle;
	const uint64_t num_tags = builder->num_tags;
	for(uint64_t tag = 0; tag < num_tags; tag++) {
		restriction_site_sequence_locator_t * const sequence_locator = restriction_locator->restriction_site_sequence_locator + tag;

		while(!sequence_locator->completed) {
			nanosleep((const struct timespec[]){{0, 100000000L}}, NULL); // Sleep for 0.1s if sequence not ready
		}
		const locator_interval_t * interval = sequence_locator->start_interval;
		fprintf(site_output,"%s",locator_interval_get_tag(locator, interval));
		const uint64_t num_blocks = sequence_locator->num_blocks;
		uint64_t prev_idx = 0;
		const uint32_t * const blk_index = sequence_locator->restriction_site_block_index;
		const uint8_t * const off = vector_get_mem(sequence_locator->restriction_site_offsets, uint8_t);
		const uint64_t start_pos = interval->begin_position;
		for(uint64_t blk = 0; blk < num_blocks; blk++) {
			if(blk_index[blk] > prev_idx) {
				for(uint64_t idx = prev_idx; idx < blk_index[blk]; idx++) {
					const uint64_t pos = start_pos + ((blk << 8) | off[idx]);
					while(pos > interval->end_position) interval++;
					fprintf(site_output," %" PRIu64, pos - interval->begin_position + interval->sequence_offset + 1);
				}
				prev_idx = blk_index[blk];
			}
		}
		fputc('\n', site_output);
	}
	return NULL;
}

void *restriction_site_init_sequence(restriction_site_locator_builder_t* const builder) {
	const archive_t * const archive = builder->archive;
	const vector_t * const restriction_sites = builder->restriction_sites;
  restriction_site_t ** sites = vector_get_mem(restriction_sites, restriction_site_t *);

  const bool bisulfite = archive->type==archive_dna_bisulfite;
	const uint8_t * text = bisulfite ?
			archive->text->bisulfite_enc_text :
			dna_text_get_text(archive->text->enc_text);

	uint64_t * D = NULL;
	uint64_t * offset = NULL;
	uint64_t ** search_patterns = NULL;

	// We only consider the individual (non-interaction) sites when searching the reference
	uint64_t num;
	for(num = 0; num < vector_get_used(restriction_sites); num++) {
		if(sites[num]->interaction_site) break;
	}
	const uint64_t num_restriction_sites = num;

	if(num_restriction_sites > 1) {
		D = mm_calloc(num_restriction_sites, uint64_t, false);
		offset = mm_calloc(num_restriction_sites, uint64_t, false);
		search_patterns = mm_calloc(num_restriction_sites, uint64_t *, false);
	}
	while(true) {
		restriction_site_sequence_locator_t * const sequence_locator = sequence_locator_sharer(builder);
		if(sequence_locator == NULL) break;
		const uint64_t text_start = sequence_locator->start_interval->begin_position;
		const uint64_t text_end = sequence_locator->end_interval->end_position;
		const uint64_t text_len = text_end - text_start + 1;
		const uint64_t num_blocks = (text_len >> 8) + (text_len & 255 ? 1 : 0);
		sequence_locator->num_blocks = num_blocks;
		uint32_t * const restriction_idx =
				sequence_locator->restriction_site_block_index = mm_calloc(num_blocks, uint32_t, false);
		uint64_t est_num_hits = 0;
		for(uint64_t i = 0; i < num_restriction_sites; i++) {
			est_num_hits += (text_len >> 2 * sites[i]->reference_pattern_len) + 1;
		}
		vector_t * restriction_off =
				sequence_locator->restriction_site_offsets = vector_new(est_num_hits, uint32_t);
		const uint8_t *ctext = text + text_start;
		if(num_restriction_sites == 1) {
			uint64_t D = sites[0]->reference_pattern_mask;
			const uint64_t offset = sites[0]->reference_pattern_len - 1;
			uint64_t const * const search_pattern = sites[0]->reference_pattern;
			for(uint64_t i = 0; i < offset; i++) D = (D >> 1) | search_pattern[(int)*ctext++];
			uint64_t len = text_len - offset;
			uint64_t pos = 0;
			uint32_t * const restriction_idx = sequence_locator->restriction_site_block_index;
			for(uint64_t blk = 0; blk <= num_blocks; blk++) {
				for(int i = 0; i < 256 && pos < len; i++, pos++) {
					D = (D >> 1) | search_pattern[(int)*ctext++];
					if(!(D & 1)) {
						vector_insert(restriction_off, (uint8_t)i, uint8_t);
					}
				}
				restriction_idx[blk] =  (uint32_t)vector_get_used(restriction_off);
			}
		} else {
			// We want these ordered so the largest offset site comes first
			for(uint64_t i = 0; i < num_restriction_sites; i++) {
				uint64_t off = sites[i]->reference_pattern_len - 1;
				uint64_t k;
				for(k = 0; k < i; k++) {
					if(off > offset[k]) {
						for(uint64_t k1 = i; k1 > k; k1--) {
							offset[k1] = offset[k1 - 1];
							search_patterns[k1] = search_patterns[k1 - 1];
							D[k1] = D[k1 - 1];
						}
						break;
					}
				}
				offset[k] = off;
				D[k] = sites[i]->reference_pattern_mask;
				search_patterns[k] = sites[i]->reference_pattern;
			}
			uint32_t last_block = 0;
			uint64_t prev_site = UINT64_MAX;
			restriction_idx[0]= 0;
			for(uint64_t pos = 0; pos < text_len; pos++) {
				uint8_t enc = *ctext++;
				for(uint64_t k = 0; k < num_restriction_sites; k++) {
					D[k] = (D[k] >> 1) | search_patterns[k][(int)enc];
					if(!(D[k] & 1)) {
						// We don't have to worry that the pos - offset[k] is in the previous interval as a match
						// can not cross an interval boundary
						uint64_t site = pos - offset[k];
						if(site != prev_site) {
							uint32_t blk = site >> 8;
							const uint32_t idx = (uint32_t)vector_get_used(restriction_off);
							for(uint32_t j = last_block + 1; j < blk; j++) restriction_idx[j] = idx;
							last_block = blk;
							restriction_idx[blk] = idx + 1;
							vector_insert(restriction_off, (uint8_t)(site & 0xff), uint8_t);
							prev_site = site;
						}
					}
				}
			}
		}
		sequence_locator->completed = true;
	}
	if(num_restriction_sites > 1) {
		mm_free(search_patterns);
		mm_free(offset);
		mm_free(D);
	}

	return NULL;
}

/*
 * Build list of restriction sites found in the reference
 */
restriction_site_locator_builder_t * const restriction_text_init_locator(
		archive_t* const archive,
		vector_t* const restriction_sites,
		restriction_site_locator_t * const restriction_site_locator,
		const char * const output_sites_name,
		const bool verbose_user,
		const uint64_t num_threads) {

	const uint64_t num_restriction_sites = vector_get_used(restriction_sites);
	if(!num_restriction_sites) return NULL;
	gem_timer_t elapsed_time;
  TIMER_START(&elapsed_time);

  restriction_site_locator_builder_t * const builder = mm_alloc(restriction_site_locator_builder_t) ;
  builder->restriction_site_locator = restriction_site_locator;
  builder->archive = archive;
	builder->restriction_sites = restriction_sites;
	builder->next_tag_idx = 0;
	builder->output_handle = NULL;

  if(output_sites_name != NULL) {
  	builder->output_handle = fopen(output_sites_name, "w");
    gem_cond_log(verbose_user && builder->output_handle != NULL,"[Writing sites out to '%s']", output_sites_name);
  }
  MUTEX_INIT(builder->mutex);
  restriction_site_t ** sites = vector_get_mem(restriction_sites, restriction_site_t *);
  gem_cond_log(verbose_user,"[Finding restriction sites for:]");
  for(uint64_t i = 0; i < num_restriction_sites; i++) {
  	if(sites[i]->interaction_site) break;
  	if(sites[i]->restriction_enzyme != NULL) gem_cond_log(verbose_user,"[...'%s']", sites[i]->restriction_enzyme);
  	else gem_cond_log(verbose_user,"[...'%.*s']", string_get_length(&sites[i]->restriction_site), string_get_buffer(&sites[i]->restriction_site));
  }
	const locator_t * const locator = archive->locator;
	const locator_interval_t * interval = locator->intervals;
	const uint64_t num_intervals = locator->num_intervals;
	const uint64_t num_tags = locator->num_tags;
	restriction_site_locator->forward_text_length = archive->text->forward_text_length;
	restriction_site_locator->bisulfite_index = archive->type == archive_dna_bisulfite;
	builder->num_tags = num_tags;
	restriction_site_locator->restriction_site_sequence_locator = mm_calloc(num_tags, restriction_site_sequence_locator_t, true);
	restriction_site_sequence_locator_t * const sequence_locator = restriction_site_locator->restriction_site_sequence_locator;
	uint64_t ctag = UINT64_MAX;
	const bs_strand_t bs_strand = archive->type==archive_dna_bisulfite ? bs_strand_C2T : bs_strand_none;
	uint64_t idx;
	for(idx = 0; idx < num_intervals; idx++) {
		if(interval[idx].strand != Forward || interval[idx].bs_strand != bs_strand) break;
		if(interval[idx].tag_id != ctag) {
			if(ctag != UINT64_MAX) sequence_locator[ctag].end_interval = interval + idx - 1;
			ctag = interval[idx].tag_id;
			sequence_locator[ctag].start_interval = interval + idx;
			sequence_locator[ctag].completed = false;
			sequence_locator[ctag].restriction_site_block_index = NULL;
			sequence_locator[ctag].restriction_site_offsets = NULL;
		}
	}
	if(ctag != UINT64_MAX)sequence_locator[ctag].end_interval = interval + idx - 1;
	if(builder->output_handle)
		gem_cond_fatal_error(
				pthread_create(&builder->restriction_site_writer, 0, (pthread_handler_t)restriction_site_write, builder), SYS_THREAD_CREATE);

	builder->threads = mm_calloc(num_threads, pthread_t, false);
	for(uint64_t i = 0; i < num_threads; i++) {
	  gem_cond_fatal_error(
	  		pthread_create(builder->threads + i, 0, (pthread_handler_t)restriction_site_init_sequence, builder), SYS_THREAD_CREATE);
	}
	for(uint64_t i = 0; i < num_threads; i++) {
		gem_cond_fatal_error(pthread_join(builder->threads[i], 0),SYS_THREAD_JOIN);
	}
	// if(builder.output_handle != NULL) {
//		gem_cond_fatal_error(pthread_join(restriction_site_writer, 0),SYS_THREAD_JOIN);
//		fclose(builder.output_handle);
//	}
	mm_free(builder->threads);
	MUTEX_DESTROY(builder->mutex);
  TIMER_STOP(&elapsed_time);
	const double calc_time = TIMER_GET_TOTAL_S(&elapsed_time);
	uint64_t num_sites = 0;
	for(uint64_t i = 0; i < num_tags; i++) {
		if(sequence_locator[i].restriction_site_offsets != NULL)
			num_sites += vector_get_used(sequence_locator[i].restriction_site_offsets);
	}
  gem_cond_log(verbose_user,"[Completed in %.2fs; Sites found = %"PRIu64"]", calc_time, num_sites);
  return builder;
}

/*
 * Find the restriction fragment interval(s) covered by a match_trace
 */
void restriction_match_trace_locate(
		match_trace_t * const match_trace,
		restriction_site_locator_t const * const restriction_site_locator) {
	if(restriction_site_locator->restriction_site_sequence_locator != NULL) {
		uint64_t start_pos = match_trace->match_alignment.match_position;

		const uint64_t forward_text_len = restriction_site_locator->forward_text_length;
		const uint64_t length =  match_trace->match_alignment.effective_length;
		// Normalize to forward strand
		if(start_pos > forward_text_len) {
			start_pos = 2 * forward_text_len - start_pos - length - 1;
		}
		// And to C2T strand if bisulfite index
		if(restriction_site_locator->bisulfite_index) {
			if(start_pos > forward_text_len >> 1) start_pos -= forward_text_len >> 1;
		}
		// Get data for sequence
		const restriction_site_sequence_locator_t * const sequence_locator = restriction_site_locator->restriction_site_sequence_locator + match_trace->tag_id;

		const uint64_t text_pos1 = start_pos - sequence_locator->start_interval->begin_position;
		const uint64_t text_pos2 = text_pos1 + length - 1;
		const uint64_t blk1 = text_pos1 >> 8;
		const uint64_t blk2 = text_pos2 >> 8;

		const uint32_t * const blocks = sequence_locator->restriction_site_block_index;
		const uint8_t * const off = vector_get_mem(sequence_locator->restriction_site_offsets, uint8_t);
//		fprintf(stderr,"%s\t%lu\n", match_trace->sequence_name, match_trace->text_position + 1);
//		fprintf(stderr,"text_pos: %lu - %lu (%lu, %lu)\n", text_pos1, text_pos2, blk1, blk2);
//		fprintf(stderr,"idx: %u %u %u\n", blocks[blk1-1], blocks[blk1], blocks[blk2]);
//		for(uint32_t i  = blocks[blk1-1]; i < blocks[blk1]; i++) {
//			fprintf(stderr,"  A -> %lu\n",(blk1 << 8) | off[i]);
//		}
//		for(uint32_t i  = blocks[blk2-1]; i < blocks[blk2]; i++) {
//			fprintf(stderr,"  B -> %lu\n",(blk1 << 8) | off[i]);
//		}
		uint32_t idx1 = blk1 > 0 ? blocks[blk1 - 1] : 0;
		uint32_t first_index = idx1;
		uint32_t second_index = idx1;
		uint32_t idx = idx1;
		uint64_t blk;
		for(blk = blk1; blk <= blk2; blk++) {
			const uint32_t idx2 = blocks[blk];
			uint64_t blk_start = blk << 8;
			for(; idx < idx2; idx++) {
				uint64_t pos = blk_start | (uint64_t)off[idx];
				if(pos <= text_pos1) first_index = second_index = idx + 1;
				else break;
			}
			if(idx < idx2) break;
		}
		for(;blk <= blk2; blk++) {
			const uint32_t idx2 = blocks[blk];
			uint64_t blk_start = blk << 8;
			for(; idx < idx2; idx++) {
				uint64_t pos = blk_start | (uint64_t)off[idx];
				if(pos >= text_pos1 && pos < text_pos2) second_index = idx + 1;
				else break;
			}
			if(idx < idx2) break;
		}
		match_trace->first_restriction_site = first_index;
		match_trace->last_restriction_site = second_index;
	} else {
		match_trace->first_restriction_site = match_trace->last_restriction_site = 0;
	}
}

