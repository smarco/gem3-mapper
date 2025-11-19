#include "text/control_sequence.h"

int control_sequence_init(control_sequence_t *cs, control_sequence_type ctype, char const *s, char const **p) {
    if(s != NULL && *s != 0) {
        cs->alt_name = NULL;
        cs->sequence_type = ctype;
        uint64_t l = strlen(s);
        string_init(&cs->sequence_name, l, NULL);
        char c;
        
        // Copy sequence name
        while((c = *s)) {
            if(c == ':' || c == ',') {
                s++;
                break;
            }
            string_append_char(&cs->sequence_name, c);
            s++;
        }
        string_append_eos(&cs->sequence_name);
        // If there is additional text, store as the alt_name
        if(*s != 0 && c==':') {
            cs->alt_name=mm_alloc(string_t);
            string_init(cs->alt_name, l - cs->sequence_name.length, NULL);
            while((c = *s++) && c!=',') {
                string_append_char(cs->alt_name, c);
            }
            string_append_eos(cs->alt_name);
            if(c==',') s++;
        }
        if(p) {
            *p = s;
        }
        return 0;
    }
    return 1;
}

char const *add_control_sequence(vector_t *v, control_sequence_type ct, char const *s) {
    vector_reserve_additional(v, 1);
    control_sequence_t *cs = vector_get_free_elm(v, control_sequence_t);
    if(control_sequence_init(cs, ct, s, &s)!=0) return NULL;
    vector_inc_used(v);
    return s;
}

void get_control_sequences(vector_t *v, control_sequence_type ct, char *emsg, char *s) {
    char const *p = s;
    while(*p) {
        p = add_control_sequence(v, ct, p);
        if(p == NULL) gem_fatal_error_msg("Error setting --%s option", emsg);
    }
}