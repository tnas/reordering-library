/*
 * Copyright 2016 Thiago Nascimento nascimenthiago@gmail.com
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#include "util.h"

int COMPARE_int_ASC (const void * a, const void * b)
{ 
	return ((int*)a) - ((int*)b);
}


int get_random_integer(int max)
{
	srand((unsigned int) time(NULL));
	return rand() % (max + 1);
}


int pow_uint(int base, const int exp)
{
	int count, val_pow;
	if (exp == 0) return 1;
	val_pow = base;
	for (count = 1; count < exp; ++count)
		val_pow *= base;
	return val_pow;
}

double get_time ()
{
	struct timeval tv; gettimeofday(&tv, 0);
	return (double)(tv.tv_sec * 100.0 + tv.tv_usec / 10000.0);
}









