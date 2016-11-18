/*
 * Copyright 2016 Brenno Lugon brennolugon@gmail.com
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

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct
{
	double*     AA;
	double*      D;
	int*        JA;
	int*        IA;
	int     m,n,nz;
} MAT;

typedef struct 
{
	double arr1;
	int    arr2;
	int    arr3;
} ARRAY;

/*----------------------------------------------------------------------------
 * MATRIX HEADER FUNCTIONS PROTOTYPE
 *--------------------------------------------------------------------------*/
extern void     MATRIX_readCSR           (MAT* A, FILE* f);
extern void 	MATRIX_read_from_path	 (const char* path_matrix_file, MAT** matrix);
extern double   MATRIX_aij               (MAT* A, int i, int j);
extern void     MATRIX_printCSR          (MAT* A);
extern void     MATRIX_printFULL         (MAT* A);
extern long int MATRIX_envelope          (MAT* A);
extern long int MATRIX_bandwidth         (MAT* A);
extern void     MATRIX_clean             (MAT* A);
extern void     MATRIX_matvec            (MAT* A, double* x, double* b);
extern void     MATRIX_forward           (MAT* L, double* b, double* y);
extern void     MATRIX_backward          (MAT* U, double* y, double* x);
extern void     MATRIX_permutation       (MAT* A, int* p);
extern void     MATRIX_writeCSR          (MAT* A, double* f, int* s, int nP, int bandwidth);
extern void 	MATRIX_readCSR_SymmUpper (MAT* A, FILE* f);
extern long int MATRIX_wavefront	 (MAT* A);
extern void     MATRIX_write_gnuplot     (const MAT *A, const char* file_name);
long int 	MATRIX_avg_nnz_per_row 	 (MAT* A);
extern int      COMPARE_array            (const void * a, const void * b);
extern int      COMPARE_eig              (const void * a, const void * b);


#endif