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


/**
 * This function calculates average value and standard deviation for an 
 * array of length elements. The minimum and maximum values are
 * excluded from the calculus.
 * 
 */
void normalize_cutbound_results(const double* values, const int length, statistic* norm_values)
{
	int exec_min_val, exec_max_val, exec;
	double sum_val, max_val, min_val, sum_pow_val;
	
	max_val = sum_val = sum_pow_val= 0;
	min_val = INT_MAX;
	
	if (length == 1) 
	{
		norm_values->average_value = values[0];
		norm_values->standard_deviation = 0;
	}
	
	// Finding out max/min value
	for (exec = 0; exec < length; ++exec)
	{
		if (values[exec] > max_val)
		{
			exec_max_val = exec;
			max_val = values[exec];
		}
		
		if (values[exec] < min_val)
		{
			exec_min_val = exec;
			min_val = values[exec];
		}
	}
	
	// Discarding max and min results value
	for (exec = 0; exec < length; ++exec)
	{
		if (exec != exec_min_val && exec != exec_max_val)
			sum_val += values[exec];
	}
	
	norm_values->average_value = sum_val / (length - 2);
	
	for (exec = 0; exec < length; ++exec)
	{
		if (exec != exec_min_val && exec != exec_max_val)
			sum_pow_val += pow(values[exec] - norm_values->average_value, 2);
	}
	
	norm_values->standard_deviation = sqrt(sum_pow_val / (length - 3));
}


/**
 * This function calculates average value and standard deviation for an 
 * array of length elements.
 * 
 * @since 29-07-2016
 */
void normalize_results(const double* values, const int length, statistic* norm_values)
{
	int exec;
	double sum_val, sum_pow_val;

	if (length == 1)
	{
		norm_values->average_value = values[0];
		norm_values->standard_deviation = 0;
	}
	
	sum_val = sum_pow_val = 0;
	
	for (exec = 0; exec < length; ++exec)
		sum_val += values[exec];
	
	norm_values->average_value = sum_val / length;
	
	for (exec = 0; exec < length; ++exec)
		sum_pow_val += pow(values[exec] - norm_values->average_value, 2);
	
	norm_values->standard_deviation = sqrt(sum_pow_val / (length - 1));
}
