/*
 * PROJECT: GEM-Tools library
 * FILE: gt_json.c
 * DATE: 01/06/2012
 * DESCRIPTION: Helper functions for json output
 */

#include "gt_json.h"

GT_INLINE JsonNode* gt_json_int_named_tuple(const uint64_t num_elements,...) {
  uint64_t i=0;
  va_list v_args;
  va_start(v_args,num_elements);
  JsonNode* a = json_mkobject();
  char* key = NULL;
  uint64_t value = 0;
  for(i=0;i<num_elements;i++) {
    key = va_arg(v_args,char*);
    value = va_arg(v_args,uint64_t);
    json_append_member(a,key,json_mknumber(value));
  }
  va_end(v_args);
  return a;
}
GT_INLINE JsonNode* gt_json_int_array(const uint64_t start,const uint64_t len,uint64_t* const data) {
  JsonNode* a = json_mkarray();
  uint64_t i= 0;
  for (i=0;i<len;++i) {
    json_append_element(a,json_mknumber((double)data[start + i]));
  }
  return a;
}
GT_INLINE JsonNode* gt_json_int_hash(gt_shash* const data) {
  JsonNode* a = json_mkobject();
  GT_SHASH_BEGIN_ITERATE(data,key,e,uint64_t) {
    json_append_member(a,key,json_mknumber((double)(*e)));
  } GT_SHASH_END_ITERATE
  return a;
}
