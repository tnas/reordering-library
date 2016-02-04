#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include "heads.h"

extern 	int get_random_integer(int max);
extern 	int pow_uint(int base, const int exp);
extern 	double get_time();
extern 	void write_output_after(const MAT *A);
extern 	void write_output_before(const MAT *A);