#define _BW_ALPHA_SIZE 256

typedef unsigned char uchar;
typedef unsigned char uint8;
typedef int int32;


// ---- struct containing the (uncompressed) bwt 
typedef struct {
  uchar *bwt;
  int size;
  int eof_pos;
} bwt_data;


// prototypes of bwt procedures defined in bwtlcp.a
void _bw_sa2bwt(uchar *t, int32 n, int32 *sa, bwt_data *b);

int32 _bw_bwt2ranknext(bwt_data *b, int32* occ, int32 *rank_next);
int32 _bw_sa2ranknext(uchar *t,int32 n,int32 *sa,int32 *occ,int32 *rank_next);
void _bw_ranknext2t(int32 *rank_next, int32 r0, bwt_data *b, uchar *t);
void _bw_ranknext2sa(int32 *rank_next, int32 r0, int32 *sa);

int32 _bw_bwt2rankprev(bwt_data *b, int32* occ, int32 *rank_prev);
int32 _bw_sa2rankprev(uchar *t,int32 n,int32 *sa,int32 *occ,int32 *rank_prev);
void _bw_rankprev2t(int32 *rank_prev, int32 rn1, bwt_data *b, uchar *t);
void _bw_rankprev2sa(int32 *rank_prev, int32 n, int32 rn1, int32 *sa);
