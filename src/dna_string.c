/*
 * PROJECT: GEMMapper
 * FILE: dna_string.c
 * DATE: 20/08/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */

#include "dna_string.h"

const bool dna_table[256] =
{
    [0 ... 255] = false,
    ['A'] = true, ['C'] = true, ['G'] = true, ['T'] = true, ['N'] = true,
    ['a'] = true, ['c'] = true, ['g'] = true, ['t'] = true, ['n'] = true
};
const bool dna_encoded_table[DNA_EXT_RANGE] =
{
  [ENC_DNA_CHAR_A] = true,
  [ENC_DNA_CHAR_C] = true,
  [ENC_DNA_CHAR_G] = true,
  [ENC_DNA_CHAR_T] = true,
  [ENC_DNA_CHAR_N] = false,
  [ENC_DNA_CHAR_SEP] = false,
  [ENC_DNA_CHAR_JUMP] = false
};
const bool extended_dna_table[256] =
{
    [0 ... 255] = false,
    ['A'] = true, ['C'] = true, ['G'] = true, ['T'] = true, ['N'] = true,
    ['a'] = true, ['c'] = true, ['g'] = true, ['t'] = true, ['n'] = true,
    [DNA_CHAR_SEP] = true,
    [DNA_CHAR_JUMP] = true
};
const bool extended_dna_encoded_table[DNA_EXT_RANGE] =
{
  [ENC_DNA_CHAR_A] = true,
  [ENC_DNA_CHAR_C] = true,
  [ENC_DNA_CHAR_G] = true,
  [ENC_DNA_CHAR_T] = true,
  [ENC_DNA_CHAR_N] = true,
  [ENC_DNA_CHAR_SEP] = true,
  [ENC_DNA_CHAR_JUMP] = true
};
const bool unmasked_dna_table[256] =
{
    [0 ... 255] = false,
    ['A'] = true, ['C'] = true, ['G'] = true, ['T'] = true,
    ['N'] = true,
};
const bool iupac_code_table[256] =
{
  [0 ... 255] = false,
  /* Upper case */
  ['A'] = true, ['C'] = true, ['G'] = true, ['T'] = true,
  ['N'] = true, ['R'] = true, ['N'] = true, ['D'] = true,
  ['E'] = true, ['H'] = true, ['I'] = true, ['L'] = true,
  ['K'] = true, ['M'] = true, ['F'] = true, ['P'] = true,
  ['S'] = true, ['W'] = true, ['Y'] = true, ['V'] = true, ['Q'] = true,
  ['B'] = true, ['Z'] = true, ['X'] = true, ['U'] = true, ['R'] = true,
  /* Lower case */
  ['a'] = true, ['c'] = true, ['g'] = true, ['t'] = true,
  ['n'] = true, ['r'] = true, ['n'] = true, ['d'] = true,
  ['e'] = true, ['h'] = true, ['i'] = true, ['l'] = true,
  ['k'] = true, ['m'] = true, ['f'] = true, ['p'] = true,
  ['s'] = true, ['w'] = true, ['y'] = true, ['v'] = true, ['q'] = true,
  ['b'] = true, ['z'] = true, ['x'] = true, ['u'] = true, ['r'] = true
};
const char dna_normalized_table[256] =
{
    [0 ... 255] = 'N',
    ['A'] = 'A', ['C'] = 'C', ['G'] = 'G', ['T'] = 'T',
    ['a'] = 'A', ['c'] = 'C', ['g'] = 'G', ['t'] = 'T',
};
const char dna_strictly_normalized_table[256] =
{
    [0 ... 255] = 'N',
    ['A'] = 'A', ['C'] = 'C', ['G'] = 'G', ['T'] = 'T',
};
const char dna_complement_table[256] =
{
  [0 ... 255] = '~',
  ['A'] = 'T', ['C'] = 'G', ['G'] = 'C',  ['T'] = 'A', ['N'] = 'N',
  ['a'] = 'T', ['c'] = 'G', ['g'] = 'C',  ['t'] = 'A', ['n'] = 'N',
  [DNA_CHAR_SEP] = DNA_CHAR_SEP,
  [DNA_CHAR_JUMP] = DNA_CHAR_JUMP
};
const uint8_t dna_encoded_complement_table[DNA_EXT_RANGE] =
{
  [ENC_DNA_CHAR_A] = ENC_DNA_CHAR_T,
  [ENC_DNA_CHAR_C] = ENC_DNA_CHAR_G,
  [ENC_DNA_CHAR_G] = ENC_DNA_CHAR_C,
  [ENC_DNA_CHAR_T] = ENC_DNA_CHAR_A,
  [ENC_DNA_CHAR_N] = ENC_DNA_CHAR_N,
  [ENC_DNA_CHAR_SEP] = ENC_DNA_CHAR_SEP,
  [ENC_DNA_CHAR_JUMP] = ENC_DNA_CHAR_JUMP
};
const uint8_t dna_encode_table[256] =
{
  [0 ... 255] = 4,
  ['A'] = 0, ['C'] = 1, ['G'] = 2,  ['T'] = 3, ['N'] = 4,
  ['a'] = 0, ['c'] = 1, ['g'] = 2,  ['t'] = 3, ['n'] = 4,
  [DNA_CHAR_SEP] = ENC_DNA_CHAR_SEP,
  [DNA_CHAR_JUMP] = ENC_DNA_CHAR_JUMP
};
const char dna_decode_table[DNA_EXT_RANGE] =
{

  [ENC_DNA_CHAR_A] = DNA_CHAR_A,
  [ENC_DNA_CHAR_C] = DNA_CHAR_C,
  [ENC_DNA_CHAR_G] = DNA_CHAR_G,
  [ENC_DNA_CHAR_T] = DNA_CHAR_T,
  [ENC_DNA_CHAR_N] = DNA_CHAR_N,
  [ENC_DNA_CHAR_SEP] = DNA_CHAR_SEP,
  [ENC_DNA_CHAR_JUMP] = DNA_CHAR_JUMP
};
const uint8_t dna_encoded_colorspace_table[DNA_EXT_RANGE][DNA_EXT_RANGE] = {
   /* A */  {ENC_DNA_CHAR_A,   ENC_DNA_CHAR_C,   ENC_DNA_CHAR_G,   ENC_DNA_CHAR_T,   ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_N},
   /* C */  {ENC_DNA_CHAR_C,   ENC_DNA_CHAR_A,   ENC_DNA_CHAR_T,   ENC_DNA_CHAR_G,   ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_N},
   /* G */  {ENC_DNA_CHAR_G,   ENC_DNA_CHAR_T,   ENC_DNA_CHAR_A,   ENC_DNA_CHAR_C,   ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_N},
   /* T */  {ENC_DNA_CHAR_T,   ENC_DNA_CHAR_G,   ENC_DNA_CHAR_C,   ENC_DNA_CHAR_A,   ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_N},
  /* SEP */ {ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_SEP},
   /* N */  {ENC_DNA_CHAR_N,   ENC_DNA_CHAR_N,   ENC_DNA_CHAR_N,   ENC_DNA_CHAR_N,   ENC_DNA_CHAR_SEP, ENC_DNA_CHAR_N}
};
