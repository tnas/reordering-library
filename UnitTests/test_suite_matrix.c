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

void load_matrix(const char* path_matrix_file, MAT** matrix)
{
	FILE* matrix_file;
	
	if ((matrix_file = fopen(path_matrix_file, "r")) == NULL) 
		exit(1);
	
	matrix = (MAT**) malloc (sizeof(MAT*));
	MATRIX_readCSR (*matrix, matrix_file);
	fclose(matrix_file);
}



void test_parallel_wavefront_inline_1()
{
	MAT* matrix;
	char* matrix_name = "../Big-Matrices/inline_1.mtx";
	long int expected_wavefront = 130603;
	long int calculated_wavefront;
	
	load_matrix(matrix_name, &matrix);
	calculated_wavefront = MATRIX_PARALLEL_wavefront(matrix);
	
	assert(calculated_wavefront == expected_wavefront);
}


void test_parallel_wavefront_audikw_1()
{
}

void test_parallel_wavefront_dielFilterV3real()
{
}

void test_parallel_wavefront_G3_circuit()
{
}

