#include "rcm_parallel_test.h"

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
	for (i = 0; i < max_level; ++i)
		s_sums[i+1] = s_sums[i] + counts[i+1];
	
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