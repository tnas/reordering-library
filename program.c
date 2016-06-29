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
// #include "./UnitTests/test_suite_matrix.h"

#define DYNAMIC_OFF 0
#define OVERWRITE_VARIABLE 1

/*
 * Program Parameters:
 * -m <path of matrix .mtx file>
 * -t <algorithm>:
 * 	0: Serial RCM
 * 	1: Serial Sloan
 * 	2: HSL RCM
 * 	3: HSL Spectral
 * 	4: HSL Sloan
 * 	5: Unordered RCM
 * 	6: Leveled RCM
 * 	7: Bucket RCM
 * 	8: Parallel Sloan
 * -p <number of threads>
 * -b <percent of chunk for Unordered RCM>
 */
int main (int argc, char* argv[]){
  
	int opt, num_threads, algorithm, test_suite;
	float bfs_chunk_size;
	char* matrix_name;
	EXECUTION exec_type;
	test defs;
	
	while ((opt = getopt(argc, argv, "m:p:b:t:a:")) != -1)
	{
		switch (opt)
		{
			case 'm' :
				matrix_name = optarg;
				exec_type = ONE_INSTANCE;
				break;
			case 't' :
				algorithm = atoi(optarg);
				break;		
				
			case 'p':
				num_threads = atoi(optarg);
				break;
				
			case 'b':
				bfs_chunk_size = atof(optarg);
				break;
				
			case 'a' :
				test_suite = atoi(optarg);
				exec_type = TEST_SUITE;
				break;
		}
	}
	
	// Disabling dynamic adjustment of the number of threads
	omp_set_dynamic(DYNAMIC_OFF);
	
	// Hinting idle threads to spin rather than sleep
	setenv("OMP_WAIT_POLICY", "PASSIVE", OVERWRITE_VARIABLE);
	
	// Prevent threads migrating between cores
	setenv("OMP_PROC_BIND", "TRUE", OVERWRITE_VARIABLE);
	
	if (exec_type == TEST_SUITE)
	{
		switch (test_suite) 
		{
			case REORDERING :
				run_all_tests();
				break;
				
			case MATRIX :
				break;
		}
	}
	else 
	{
		defs.path_matrix_file = matrix_name;
		defs.algorithm        = algorithm;
		defs.percent_chunk    = bfs_chunk_size;
		defs.num_threads      = num_threads;
		
		test_reorder_algorithm(defs);
	}
	
	return 0;
}