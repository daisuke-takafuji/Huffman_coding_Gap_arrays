#ifndef CU_TIMER_
#define CU_TIMER_

#include<time.h>

double get_ms(struct timespec start, struct timespec stop);
static double sum_of_time = 0;

#define TIMER_START(id)									\
	cudaEvent_t start##id, stop##id;						\
	cudaEventCreate(&start##id);							\
	cudaEventCreate(&stop##id);								\
	cudaEventRecord(start##id);

#define TIMER_STOP(id)										\
	cudaEventRecord(stop##id);								\
	cudaEventSynchronize(stop##id);							\
	float time##id;											\
	cudaEventElapsedTime(&time##id, start##id, stop##id);	\
	cudaEventDestroy(start##id);							\
	cudaEventDestroy(stop##id);								\
	sum_of_time += time##id;								\
	printf("%s,%f, ", #id,time##id);

#define CPU_TIMER_START(id)										\
	struct timespec start##id, stop##id;					\
	clock_gettime(CLOCK_REALTIME, &start##id);				\


#define CPU_TIMER_STOP(id)									\
	clock_gettime(CLOCK_REALTIME, &stop##id);				\
	sum_of_time += get_ms(start##id, stop##id);					\
	printf("%s(cpu),%lf, ", #id, get_ms(start##id, stop##id));

#define CPU_TIMER_START2(id)										\
	struct timespec start##id, stop##id;					\
	clock_gettime(CLOCK_REALTIME, &start##id);				\


#define CPU_TIMER_STOP2(id)									\
	clock_gettime(CLOCK_REALTIME, &stop##id);				\
	double time##id = get_ms(start##id, stop##id);			\
	printf("%lf, ", get_ms(start##id, stop##id));		\
	printf("%s(cpu), ", #id);

#endif // TIMER_
