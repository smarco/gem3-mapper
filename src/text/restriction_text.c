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
		{NULL, NULL}
};

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
    	char const *p = string_get_buffer(&rest->restriction_site);
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

     	// And now for patterns for searching the reads
     	// In the reads we should see the reverse strand sequence from the end of the site until the cut site fused to
     	// the forward strand sequence after the cut site.
     	// For example, with HindIII, the sequence is A|AGCTT, so in the reads we should see AAGCTAGCTT

    	x = (uint64_t)1;
    	const uint64_t half_len = search_pattern_len >> 1;
    	const uint64_t cut_index = rest->cut_site_index;
    	for(uint64_t i = 0; i < search_pattern_len; i++) {
    		uint8_t enc = i < half_len ?
    				iupac[(int)p[ref_pattern_len - 1 - i]] :
    				iupac[(int)p[ref_pattern_len - 1 - i - cut_index + half_len]];
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
     	mask = ((uint64_t)1 << search_pattern_len) - (uint64_t)1;
    	for(int i = 0; i < DNA_EXT_RANGE; i++) {
    		rest->search_pattern[bs_strand_none][i] ^= mask;
    		rest->search_pattern[bs_strand_C2T][i] ^= mask;
    		rest->search_pattern[bs_strand_G2A][i] ^= mask;
    	}
    	rest->search_pattern_mask = mask;
    }
  } else fprintf(stderr,"NULL restriction site description passed to restriction_new()\n");
  return rest;
}

void restriction_site_delete(restriction_site_t * const rest) {
  if(rest != NULL) {
    string_destroy(&rest->restriction_site);
    mm_free(rest);
  }
}

void find_restriction_site_matches(
		pattern_t const * const pattern,
		vector_t const * const restriction_sites,
		vector_t * const restriction_hits,
		bisulfite_conversion_t const bisulfite_conversion) {

	vector_clear(restriction_hits);
	uint64_t number_hits = 0;
	const uint64_t num_restriction_sites = vector_get_used(restriction_sites);
	if(num_restriction_sites == 0) return;
	restriction_site_t **sites = vector_get_mem(restriction_sites, restriction_site_t *);
	const uint64_t len = pattern->key_length;
	uint8_t const *key = pattern->key_non_bs;
	if(num_restriction_sites == 1) {
		uint64_t const * const search_pattern = sites[0]->search_pattern[bisulfite_conversion];
  	uint64_t D = sites[0]->search_pattern_mask;
  	uint64_t const site_len = sites[0]->search_pattern_len;
		for(uint64_t i = 0; i < len; i++) {
			D = (D >> 1) | search_pattern[(int)*key++];
			if(!(D & 1)) { // Full length match found
				vector_reserve_additional(restriction_hits, 1);
				uint64_t *hit = vector_get_free_elm(restriction_hits, uint64_t);
				*hit = i + 1 - (site_len >> 1);
				vector_set_used(restriction_hits, ++number_hits);
			}
		}
	}
}

void restriction_text_init_locator(
		const archive_t* const archive,
		const vector_t* const restriction_sites,
		restriction_site_locator_t * const restriction_site_locator,
		const char * const output_sites_name,
		const bool verbose_user) {

	const uint64_t num_restriction_sites = vector_get_used(restriction_sites);
	if(!num_restriction_sites) return;
	restriction_site_t ** sites = vector_get_mem(restriction_sites, restriction_site_t *);
	const bool bisulfite = archive->type==archive_dna_bisulfite;
	const uint8_t * text = bisulfite ?
			archive->text->bisulfite_enc_text :
			dna_text_get_text(archive->text->enc_text);
	const uint64_t len = bisulfite ?
			archive->text->forward_text_length >> 1 :
			archive->text->forward_text_length;
	gem_timer_t elapsed_time;

  TIMER_START(&elapsed_time);
  // Load archive
  gem_cond_log(verbose_user,"[Finding restriction sites for:]");
  for(uint64_t i = 0; i < num_restriction_sites; i++) {
  	if(sites[i]->restriction_enzyme != NULL) gem_cond_log(verbose_user,"[...'%s']", sites[i]->restriction_enzyme);
  	else gem_cond_log(verbose_user,"[...'%.*s']", string_get_length(&sites[i]->restriction_site), string_get_buffer(&sites[i]->restriction_site));
  }
  FILE *site_output = NULL;
  if(output_sites_name != NULL) {
  	site_output = fopen(output_sites_name, "w");
    gem_cond_log(verbose_user,"[Writing sites out to '%s']", output_sites_name);
  }
	const locator_t * const locator = archive->locator;
	const uint64_t num_tags = locator->num_tags;
	const locator_interval_t * interval = locator->intervals;
	restriction_site_locator->forward_text_length = archive->text->forward_text_length;
	restriction_site_locator->bisulfite_index = bisulfite;
	restriction_site_locator->sequence_offsets = mm_calloc(num_tags, uint64_t, false);
	vector_t * const restriction_off = vector_new(8192, uint8_t);
	uint64_t tag_id = UINT64_MAX;
	uint64_t begin = interval->begin_position;
	uint64_t end = interval->end_position;
	uint64_t seq_offset = interval->sequence_offset;
	restriction_site_locator->restriction_site_offsets = restriction_off;
	if(num_restriction_sites == 1) {
		uint64_t D = sites[0]->reference_pattern_mask;
		const uint64_t offset = sites[0]->reference_pattern_len - 1;
		uint64_t const * const search_pattern = sites[0]->reference_pattern;
		for(uint64_t i = 0; i < offset; i++) D = (D >> 1) | search_pattern[(int)*text++];
		uint64_t len1 = len - offset;
		uint64_t pos = 0;
		const uint64_t number_blocks = len1 >> 8;
		vector_t * const restriction_idx = vector_new(number_blocks + 1, uint32_t);
		restriction_site_locator->restriction_site_block_index = restriction_idx;
		for(uint64_t blk = 0; blk <= number_blocks; blk++) {
			const uint64_t lim = blk < number_blocks ? 256 : len1 & 255;
			for(int i = 0; i < lim; i++, pos++) {
				if(pos >= end) {
					interval++;
					if(interval->tag_id != tag_id) {
						if(site_output != NULL) {
						if(tag_id != UINT64_MAX) fprintf(site_output," %" PRIu64 "\n", end - begin + seq_offset);
							fprintf(site_output,"%s",locator_interval_get_tag(locator, interval));
						}
						tag_id = interval->tag_id;
						restriction_site_locator->sequence_offsets[tag_id] = vector_get_used(restriction_off);
					}
					end = interval->end_position;
					begin = interval->begin_position;
					seq_offset = interval->sequence_offset;
				}
				D = (D >> 1) | search_pattern[(int)*text++];
				if(!(D & 1)) {
					vector_insert(restriction_off, (uint8_t)i, uint8_t);
					if(site_output != NULL) fprintf(site_output," %" PRIu64, pos - begin + seq_offset + 1);
				}
			}
			if(lim > 0)	vector_insert(restriction_idx, (uint32_t)(vector_get_used(restriction_off)), uint32_t);
		}
	} else {
		uint64_t * const D = mm_calloc(num_restriction_sites, uint64_t, false);
		uint64_t * const offset = mm_calloc(num_restriction_sites, uint64_t, false);
		uint64_t ** search_pattern = mm_calloc(num_restriction_sites, uint64_t *, false);
		// We want these ordered so the largest offset site comes first
		for(uint64_t i = 0; i < num_restriction_sites; i++) {
			uint64_t off = sites[i]->reference_pattern_len - 1;
			uint64_t k;
			for(k = 0; k < i; k++) {
				if(off > offset[k]) {
					for(uint64_t k1 = i; k1 > k; k1--) {
						offset[k1] = offset[k1 - 1];
						search_pattern[k1] = search_pattern[k1 - 1];
						D[k1] = D[k1 - 1];
					}
					break;
				}
			}
			offset[k] = off;
			D[k] = sites[i]->reference_pattern_mask;
			search_pattern[k] = sites[i]->reference_pattern;
		}
		const uint64_t number_blocks = len >> 8;
		vector_t * const restriction_idx = vector_new(number_blocks + 1, uint32_t);
		restriction_site_locator->restriction_site_block_index = restriction_idx;
		uint32_t *block_idx = vector_get_mem(restriction_idx, uint32_t);
		uint32_t last_block = 0;
		uint64_t prev_site = UINT64_MAX;
		block_idx[0] = 0;
		for(uint64_t pos = 0; pos < len; pos++) {
			if(pos > end) {
				interval++;
				if(interval->tag_id != tag_id) {
					if(site_output != NULL) {
						if(tag_id != UINT64_MAX) fprintf(site_output," %" PRIu64 "\n", end - begin + seq_offset);
						fprintf(site_output,"%s",locator_interval_get_tag(locator, interval));
					}
					tag_id = interval->tag_id;
					restriction_site_locator->sequence_offsets[tag_id] = vector_get_used(restriction_off);
				}
				begin = interval->begin_position;
				end = interval->end_position;
				seq_offset = interval->sequence_offset;
			}
			uint8_t enc = *text++;
			for(uint64_t k = 0; k < num_restriction_sites; k++) {
				D[k] = (D[k] >> 1) | search_pattern[k][(int)enc];
				if(!(D[k] & 1)) {
					// We don't have to worry that the pos - offset[k] is in the previous interval as a match
					// can not cross an interval boundary
					uint64_t site = pos - offset[k];
					if(site != prev_site) {
						uint32_t blk = site >> 8;
						uint32_t idx = (uint32_t)vector_get_used(restriction_idx);
						for(uint32_t j = last_block + 1; j < blk; j++) block_idx[j] = idx;
						last_block = blk;
						block_idx[blk] = idx + 1;
						vector_insert(restriction_off, (uint8_t)(site & 0xff), uint8_t);
						prev_site = site;
						if(site_output != NULL) fprintf(site_output," %" PRIu64, site - begin + seq_offset + 1);
					}
				}
			}
		}
		vector_set_used(restriction_idx, last_block + 1);
		mm_free(search_pattern);
		mm_free(offset);
		mm_free(D);
	}
	if(site_output != NULL) {
		fprintf(site_output," %" PRIu64 "\n", len - begin + seq_offset);
		fclose(site_output);
	}
  TIMER_STOP(&elapsed_time);
	const double calc_time = TIMER_GET_TOTAL_S(&elapsed_time);
  gem_cond_log(verbose_user,"[Completed in %.2fs; Sites found = %"PRIu64"]", calc_time, vector_get_used(restriction_off));
}

void restriction_match_trace_locate(
		match_trace_t * const match_trace,
		restriction_site_locator_t const * const restriction_site_locator) {
	if(restriction_site_locator->restriction_site_block_index != NULL) {
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
		// Get first index for chromosome/contig
		const uint64_t index_offset = restriction_site_locator->sequence_offsets[match_trace->tag_id];
		const uint64_t text_pos1 = start_pos;
		const uint64_t text_pos2 = text_pos1 + length - 1;
		const uint64_t blk1 = text_pos1 >> 8;
		const uint64_t blk2 = text_pos2 >> 8;

		const uint32_t * const blocks = vector_get_mem(restriction_site_locator->restriction_site_block_index, uint32_t);
		const uint8_t * const off = vector_get_mem(restriction_site_locator->restriction_site_offsets, uint8_t);
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
		match_trace->first_restriction_site = first_index - index_offset;
		match_trace->last_restriction_site = second_index - index_offset;
	} else {
		match_trace->first_restriction_site = match_trace->last_restriction_site = 0;
	}
}

