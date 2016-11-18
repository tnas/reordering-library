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

#include "test_suite_matrix.h"

void test_parallel_max_wavefront()
{
	MAT* matrix;
	int num_matrices, size_set_threads, mat, th;
	long int calculated_wavefront;
	
	int nthreads[] = { 1, 2, 4, 6, 8 };
	
	char* matrices[] = {
		"../Matrices/hsl.mtx",
	};
	
	long int expected_wavefront[] = {
		5
	};
	
	num_matrices     = sizeof(matrices)/sizeof(matrices[0]);
	size_set_threads = sizeof(nthreads)/sizeof(nthreads[0]);
	
	for (mat = 0; mat < num_matrices; ++mat)
	{
		MATRIX_read_from_path(matrices[mat], &matrix);
		
		for (th = 0; th < size_set_threads; ++th)
		{
			omp_set_num_threads(nthreads[th]);
			
			calculated_wavefront = MATRIX_PARALLEL_max_wavefront(matrix);
			
			printf("Calculated wavefront for matrix %s: %ld\n", 
			       matrices[mat], calculated_wavefront);fflush(stdout);
			       
			assert(calculated_wavefront == expected_wavefront[mat]);
			
			printf("Test of wavefront of matrix %s and %d threads ----- OK\n", 
			       matrices[mat], nthreads[th]);fflush(stdout);
		}
		
		MATRIX_clean(matrix);
	}
}


void test_parallel_rms_wavefront()
{
	MAT* matrix;
	int num_matrices, size_set_threads, mat, th;
	long int calculated_wavefront;
	
	int nthreads[] = { 1, 2, 4, 6, 8 };
	
	char* matrices[] = {
		"../Big-Matrices/inline_1.mtx",
		"../Big-Matrices/audikw_1.mtx",
		"../Big-Matrices/dielFilterV3real.mtx",
		"../Big-Matrices/G3_circuit.mtx"
	};
	
	long int expected_wavefront[] = {
		130603, 477149, 598949, 85516
	};
	
	num_matrices     = sizeof(matrices)/sizeof(matrices[0]);
	size_set_threads = sizeof(nthreads)/sizeof(nthreads[0]);
	
	for (mat = 0; mat < num_matrices; ++mat)
	{
		MATRIX_read_from_path(matrices[mat], &matrix);
		
		for (th = 0; th < size_set_threads; ++th)
		{
			omp_set_num_threads(nthreads[th]);
			
			calculated_wavefront = MATRIX_PARALLEL_rms_wavefront(matrix);
			
			printf("Calculated wavefront for matrix %s: %ld\n", 
			       matrices[mat], calculated_wavefront);fflush(stdout);
			assert(calculated_wavefront == expected_wavefront[mat]);
			printf("Test of wavefront of matrix %s and %d threads ----- OK\n", 
			       matrices[mat], nthreads[th]);fflush(stdout);
		}
		
		MATRIX_clean(matrix);
	}
}


void run_all_test_matrix()
{
	test_parallel_max_wavefront();
	test_parallel_rms_wavefront();
}
