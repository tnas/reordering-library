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
#include <assert.h>
#include <limits.h>
#include <float.h>
#include "../CommonFiles/matrix_parallel.h"
#include "../CommonFiles/graph_hsl.h"
#include "../Reorderings/hsl.h"
#include "../Reorderings/rcm.h"
#include "../Reorderings/sloan.h"
#include "../Reorderings/spectral.h"
#include "../Reorderings/rcm_parallel.h"
#include "../Reorderings/rcm_hsl.h"
#include "../Reorderings/sloan_parallel.h"



#define TEST_EXEC_TIMES 1

typedef enum {
	serial_rcm,
	serial_sloan,
	hsl_rcm,
	hsl_spectral,
	hsl_sloan,
	unordered_rcm,
	leveled_rcm,
	bucket_rcm,
	parallel_sloan
} reorder_algorithm;


typedef struct {
	double time_peripheral;
	double time_reordering;
	double time_permutation;
	long int original_band;
	long int reorder_band;
	long int wavefront;
	const char* path_matrix_file;
	int root;
	int start_node;
	int end_node;
	char* algorithm_name;
	int num_threads;
	reorder_algorithm algorithm;
	float percent_chunk;
} test;


typedef enum { START, END } PERIPHERAL_NODES;


test test_reorder_algorithm  (test defs);
void run_all_reordering_tests();