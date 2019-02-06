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

static char iupac[256] = {
  ['A'] = 1, ['C'] = 2, ['G'] = 4, ['T'] = 8,
  ['R'] = 5, ['Y'] = 10, ['S'] = 6, ['W'] = 9, ['K'] = 12, ['M'] = 3,
  ['B'] = 14, ['D'] = 13, ['H'] = 11, ['V'] = 7, ['N'] = 15,
  ['a'] = 1, ['c'] = 2, ['g'] = 4, ['t'] = 8,
  ['r'] = 5, ['y'] = 10, ['s'] = 6, ['w'] = 9, ['k'] = 12, ['m'] = 3,
  ['b'] = 14, ['d'] = 13, ['h'] = 11, ['v'] = 7, ['n'] = 15,
};

static char iupac_compl[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

restriction_t *restriction_new(char * const rest_char) {
  restriction_t *rest = NULL;
  if(rest_char != NULL) {
    rest = mm_alloc(restriction_t);
    string_init(&rest->restriction_site, 0, NULL);
    char * const tp = rest_char;
    char c;
    int index = -1;
    int i = 0;
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
      restriction_delete(rest);
      rest = NULL;
    } else {
      rest->cut_site_index = index;
      fprintf(stderr,"Added restriction site: '%.*s', cut index = %d\n", PRIs_content(&rest->restriction_site), rest->cut_site_index);
    }
  } else fprintf(stderr,"NULL restriction site description passed to restriction_new()\n");
  return rest;
}

void restriction_delete(restriction_t * const rest) {
  if(rest != NULL) {
    string_destroy(&rest->restriction_site);
    mm_free(rest);
  }
}
