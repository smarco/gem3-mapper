#include "text/control_sequence.h"

control_sequence_t *control_sequence_new(control_sequence_type ctype, char const *s, char const **p) {
    control_sequence_t *cs = NULL;
    if(s != NULL && *s != 0) {
        cs = mm_alloc(control_sequence_t);
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
            string_init(cs->alt_name, l - cs->sequence_name.length, NULL);
            while((c = *s++)) {
                string_append_char(cs->alt_name, c);
            }
            string_append_eos(cs->alt_name);
        }
        if(p) {
            *p = s;
        }
    }
    return cs;
}

void get_control_sequences(vector_t *v, control_sequence_type ct, char *emsg, char *s) {
    char const *p = s;
    while(*p) {
        control_sequence_t *cs = control_sequence_new(ct, p, &p);
        if(cs != NULL) {
            vector_insert(v, cs, control_sequence_t *);
        } else {
            gem_fatal_error_msg("Error setting --%s option", emsg);
        }
    }
}