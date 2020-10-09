#include "../include/cu_timer.h"

double get_ms(
		struct timespec start,
		struct timespec stop){

	double ms=0;
	ms = (double)(stop.tv_sec - start.tv_sec)*1000*1000 + (double)(stop.tv_nsec - start.tv_nsec)/1000;
	return ms;
}
