/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   global defintion for the ds suffix-sort algorithm 
   Giovanni Manzini 
   2-apr 2001
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/* ----------- use assertion if DEBUG!=0 ------------- */
#ifndef DEBUG
#define DEBUG 1   /* set DEBUG to 0 to remove assertions and extra checks */
#endif
#if !DEBUG
#define NDEBUG 1  /* do not compile assertions */
#endif
#include <stdint.h>
#include <assert.h>

/* ---------- types and costants ----------- */
typedef int32_t        Int32;
typedef uint32_t UInt32;
typedef uint16_t UInt16;
typedef char          Char;
typedef unsigned char UChar;
typedef unsigned char Bool;
typedef uint64_t UInt64;
#define True   ((Bool)1)
#define False  ((Bool)0)
#define Cmp_overshoot 16
#define Max_thresh 30


#ifndef min
#define min(a, b) ((a)<=(b) ? (a) : (b))
#endif

#ifndef max
#define max(a, b) ((a)>=(b) ? (a) : (b))
#endif


#define MIN(a, b) ((a)<=(b) ? (a) : (b))
#define MAX(a, b) ((a)>=(b) ? (a) : (b))


// constant and macro for marking groups
#define SETMASK (1 << 30)
#define CLEARMASK (~(SETMASK))
#define IS_SORTED_BUCKET(sb) (ftab[sb] & SETMASK)
#define BUCKET_FIRST(sb) (ftab[sb]&CLEARMASK)
#define BUCKET_LAST(sb) ((ftab[sb+1]&CLEARMASK)-1)
#define BUCKET_SIZE(sb) ((ftab[sb+1]&CLEARMASK)-(ftab[sb]&CLEARMASK))

int scmp3(unsigned char *p, unsigned char *q, int *l, int maxl);
void pretty_putchar(int c);







