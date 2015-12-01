

#ifndef GPU_SAMPLE_H_
#define GPU_SAMPLE_H_

#ifndef _POSIX_C_SOURCE
	#define _POSIX_C_SOURCE	 199309L
#endif

#include "gpu_commons.h"
#include <time.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

double gpu_sample_time();

#endif /* GPU_SAMPLE_H_ */
