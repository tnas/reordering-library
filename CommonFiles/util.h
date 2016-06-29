#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>

#define min(a,b) (((a)>(b))?(b):(a))
#define max(a,b) (((a)>(b))?(a):(b))

extern  int COMPARE_int_ASC (const void * a, const void * b);
extern 	int get_random_integer(int max);
extern 	int pow_uint(int base, const int exp);
extern 	double get_time();
