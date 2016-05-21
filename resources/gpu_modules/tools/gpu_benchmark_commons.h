/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2013-2016 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#include <time.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#define MASK_ZEROS					0x00000000
#define	UINT32_LENGTH				32 			  // 32 bits
#define SEED_CHAR_LENGTH	  2 		    //  2 bits
#define SEED_FIELD_SIZE			8			    //  8 bits
#define	UINT64_LENGTH				64 			  // 64 bits

/* Encoded DNA Nucleotides */
#define ENC_DNA_CHAR_A      0
#define ENC_DNA_CHAR_C      1
#define ENC_DNA_CHAR_G      2
#define ENC_DNA_CHAR_T      3

#define GPU_DIV_CEIL(NUMERATOR,DENOMINATOR) (((NUMERATOR)+((DENOMINATOR)-1))/(DENOMINATOR))
#define MIN(_a, _b) (((_a) < (_b)) ? (_a) : (_b))

inline double sample_time()
{
	struct timespec tv;

	#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
	clock_serv_t cclock;
	mach_timespec_t mts;
	host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
	clock_get_time(cclock, &mts);
	mach_port_deallocate(mach_task_self(), cclock);
	tv.tv_sec = mts.tv_sec;
	tv.tv_nsec = mts.tv_nsec;

	#else
	clock_gettime(CLOCK_REALTIME, &tv);
	#endif

	return((tv.tv_sec+tv.tv_nsec/1000000000.0));
}

inline char BinASCIItoChar(uint32_t base)
{
	switch(base)
	{
    	case ENC_DNA_CHAR_A:
    		return('A');
    	case ENC_DNA_CHAR_C:
    		return('C');
    	case ENC_DNA_CHAR_G:
    		return('G');
    	case ENC_DNA_CHAR_T:
    		return('T');
    	default :
    		return('X');
	}
}

uint64_t charToBinASCII(unsigned char base)
{
	switch(base)
	{
    	case 'A':
    	case 'a':
    	    return(ENC_DNA_CHAR_A);
    	case 'C':
    	case 'c':
    	    return(ENC_DNA_CHAR_C);
    	case 'G':
    	case 'g':
    	    return(ENC_DNA_CHAR_G);
    	case 'T':
    	case 't':
    	    return(ENC_DNA_CHAR_T);
    	default :
    	    return(ENC_DNA_CHAR_A);
	}
}
