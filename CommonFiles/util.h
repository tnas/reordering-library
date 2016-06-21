#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include "heads.h"

#define iseven(n) ((n)%(2)==(0)?(1):(0))
#define isdivisor(d, n) ((n)%(d)==(0)?(1):(0))

typedef struct 
{
	int initial_prefix_sum;
	int curr_prefix_sum;
	int curr_total_sum;
	int last_prefix_sum;
	int last_total_sum;
} status_prefix_sum;

extern  int COMPARE_int_ASC (const void * a, const void * b);
extern 	int get_random_integer(int max);
extern 	int pow_uint(int base, const int exp);
extern 	double get_time();
extern 	void write_output_after(const MAT *A);
extern 	void write_output_before(const MAT *A);
extern void  prefix_sum(const int* counts, int** sums, const int max_level);