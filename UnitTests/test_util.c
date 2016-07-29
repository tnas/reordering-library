/*
 * Copyright 2016 Thiago Nascimento <nascimenthiago@gmail.com>
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

#include "test_util.h"

void test_normalize_results_1()
{
	statistic expected_results;
	int length = 10;
	double values[10] = {
		2, 3, 3, 4, 5, 6, 7, 8, 9, 10
	};
	double avg, stdev;
	
	normalize_results(values, length, &expected_results);
	
	avg   = round(expected_results.average_value * 100)/100;
	stdev = round(expected_results.standard_deviation*100)/100;
	
	assert(avg == 5.7);
	assert(stdev == 2.75);
}

void test_normalize_results_2()
{
	statistic expected_results;
	int length = 10;
	double values[10] = {
		5, 5, 7, 8, 10, 3, 7, 8, 9, 7 
	};
	double stdev;
	
	normalize_results(values, length, &expected_results);
	
	stdev = round(expected_results.standard_deviation*100 + .5)/100;
	
	assert(stdev == 2.08);
}


void test_normalize_cutbounds_results_1()
{
	statistic expected_results;
	int length = 10;
	double values[10] = {
		2, 3, 3, 4, 5, 6, 7, 8, 9, 10
	};
	double stdev;
	
	normalize_cutbound_results(values, length, &expected_results);
	
	stdev = round(expected_results.standard_deviation*100 + .5)/100;
	
	assert(stdev == 2.27);
}

void run_all_util_tests()
{
	   test_normalize_results_1();
	   test_normalize_results_2();
	   test_normalize_cutbounds_results_1();
}