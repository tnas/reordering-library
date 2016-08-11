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
#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>

#define min(a,b) (((a)>(b))?(b):(a))
#define max(a,b) (((a)>(b))?(a):(b))

typedef struct 
{
	double average_value;
	double standard_deviation;
} statistic;


int    COMPARE_int_ASC		      (const void * a, const void * b);
int    get_random_integer	      (int max);
int    pow_uint			      (int base, const int exp);
double get_time			      ();
void   normalize_results	      (const double* values, const int length, statistic* norm_values);
void   normalize_int_results          (const long int* values, const int length, statistic* norm_values);
void   normalize_cutbound_results     (const double* results, const int length, statistic* norm_values);
void   normalize_cutbound_int_results (const long int* values, const int length, statistic* norm_values);

#endif