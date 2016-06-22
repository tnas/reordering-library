/*
 * Copyright 2016 Thiago Nascimento <thiago.rodrigues@tre-es.jus.br>
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

#include "test_suite_prefixsum.h"

void test_prefix_sum()
{
	int* sums;
	int max_level, i;
	
	max_level = 15;
	int counts[15] =       { 9, 8, 3, 2, 7, 1, 6, 4, 5, 5, 8, 3, 4, 9, 1 };
	int correct_sums[15] = { 9, 17, 20, 22, 29, 30, 36, 40, 45, 50, 58, 61, 65, 74, 75 };
	
	prefix_sum(counts, &sums, max_level);
	
	for (i = 0; i < max_level; ++i)
		assert(sums[i] == correct_sums[i]);
}


void test_prefix_sum_parallel_serial(int* counts, int max_level)
{
	int i;
	int* p_sums;
	int* s_sums;
	
	// Parallel prefix sum
	prefix_sum(counts, &p_sums, max_level);
	
	// Serial prefix sum
	s_sums = calloc(max_level, sizeof(int));
	s_sums[0] = counts[0];
	for (i = 1; i < max_level; ++i)
		s_sums[i] = s_sums[i-1] + counts[i];
	
	// Comparing parallel and serial prefix sum
	for (i = 0; i < max_level; ++i) 
	{
		if (s_sums[i] != p_sums[i]) 
		{
			printf("Serial prefix sum not equals Parallel prefix sum: position %d: s = %d and p = %d\n",
				i, s_sums[i], p_sums[i]);
			exit(1);
		}
	}
		
	free(s_sums);
	free(p_sums);
}
