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
#include <unistd.h>
#include <ctype.h>
#include "./UnitTests/test_suite.h"
#include "./UnitTests/test_suite_reordering.h"
#include "./UnitTests/test_suite_prefixsum.h"
#include "./UnitTests/test_suite_matrix.h"
#include "./UnitTests/test_graph_parallel.h"
#include "./UnitTests/test_util.h"
#include "./UnitTests/test_linked_list.h"

#define DYNAMIC_OFF 0
#define OVERWRITE_VARIABLE 1

/* *****************************************************
 * Program Parameters:
 * 
 * -m <path of matrix .mtx file>
 * -a <algorithm>:
 * 	0: Serial RCM
 * 	1: Serial Sloan
 * 	2: HSL RCM
 * 	3: HSL Spectral
 * 	4: HSL Sloan
 * 	5: Unordered RCM
 * 	6: Leveled RCM
 * 	7: Bucket RCM
 * 	8: Relaxed Order Sloan
 * 	9: Boost RCM
 *     10: Boost Sloan
 *     11: Logical Bag Sloan
 *     12: Shrinked RCM
 * -p <number of threads>
 * -b <percent of chunk for Unordered RCM>
 * -t <test suite>:
 * 	0: all tests
 * 	1: reordering algorithms 
 * 	2: prefix sum 
 * 	3: matrix operations
 * 	4: graph parallel
 * 	5: util functions
 * 	6: linked list operations
 * 	7: TEMA journal
 * 	8: Dissertation Biggest Matrices
 *      9: Dissertation Smallest Matrices
 *     10: SBPO-2017 conference test
 * *****************************************************
 */
int main (int argc, char* argv[]) {
  
	int opt, num_threads, algorithm, test_suite;
	float bfs_chunk_size;
	char* matrix_name;
	test_scope scope;
	test defs;
	
	scope = TEST_CASE;
	
	while ((opt = getopt(argc, argv, "m:p:b:t:a:")) != -1)
	{
		switch (opt)
		{
			case 'm' :
				matrix_name = optarg;
				scope = TEST_CASE;
				break;
			case 'a' :
				algorithm = atoi(optarg);
				break;		
				
			case 'p':
				num_threads = atoi(optarg);
				break;
				
			case 'b':
				bfs_chunk_size = atof(optarg);
				break;
				
			case 't' :
				test_suite = atoi(optarg);
				scope = TEST_SUITE;
				break;
		}
	}
	
	// Disabling dynamic adjustment of the number of threads
	omp_set_dynamic(DYNAMIC_OFF);
	
	// Hinting idle threads to slep rather than spin
	setenv("OMP_WAIT_POLICY", "PASSIVE", OVERWRITE_VARIABLE);
	
	// Prevent threads migrating between cores
	setenv("OMP_PROC_BIND", "TRUE", OVERWRITE_VARIABLE);
	
	if (scope == TEST_SUITE)
	{
		switch (test_suite) 
		{
			case ALL : 
				run_reordering_tests();
				test_prefix_sum();
				run_all_test_matrix();
				run_all_test_GRAPH_parallel();
				run_all_util_tests();
				break;
				
			case REORDERING :
				run_reordering_tests();
				break;
			
			case PREFIX_SUM:
				test_prefix_sum();
				break;
				
			case MATRIX :
				run_all_test_matrix();
				break;
			
			case GRAPH_PARALLEL :
				run_all_test_GRAPH_parallel();
				run_all_comparison_GRAPH_parallel();
				break;
			
			case UTIL_FUNCIONS :
				run_all_util_tests();
				break;
				
			case LINKED_LIST :
				run_all_linked_list_tests();
				break;
			
			case TEMA :
				run_tema_journal_tests();
				break;
			
			case DISSERTATION_LARGEST :
				run_dissertation_largest_matrices();
				break;
				
			case DISSERTATION_SMALLEST :
				run_dissertation_smallest_matrices();
				break;
				
			case SBPO2017 :
				run_sbpo2017_tests();
				break;
		}
	}
	else if (scope == TEST_CASE)
	{
		defs.path_matrix_file = matrix_name;
		defs.algorithm        = algorithm;
		defs.percent_chunk    = bfs_chunk_size;
		defs.num_threads      = num_threads;
		
		test_reorder_algorithm(defs);
	}
	else
	{
		printf("*** [Error] Test not identified ***\n");
		exit(1);	
	}
	
	return 0;
}