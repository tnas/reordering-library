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

#include "matrix.h"
#include "util_parallel.h"

#ifndef __MATRIX_PARALLEL_H__
#define __MATRIX_PARALLEL_H__

#define MATRIX_NUM_THREADS 8

long int MATRIX_PARALLEL_bandwidth   	(MAT* A);
long int MATRIX_PARALLEL_max_wavefront  (MAT* A);
long int MATRIX_PARALLEL_rms_wavefront  (MAT* A);
void     MATRIX_PARALLEL_permutation    (MAT* A, int* p);
void 	 MATRIX_PARALLEL_wavefront	(const MAT* A, int** wf_per_row);

#endif /* __MATRIX_PARALLEL_H__ */
