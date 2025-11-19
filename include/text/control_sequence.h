#include <utils/essentials.h>

typedef enum {
	SequenceControl=0, UnderConversion=1, OverConversion=2, Conversion=3
} control_sequence_type;

typedef struct {
	string_t sequence_name;
	string_t *alt_name;
	control_sequence_type sequence_type;
} control_sequence_t;

char const *add_control_sequence(vector_t *v, control_sequence_type ct, char const *s);
int control_sequence_init(control_sequence_t *,control_sequence_type, char const *, char const **);
void get_control_sequences(vector_t *v, control_sequence_type ct, char *emsg, char *s);