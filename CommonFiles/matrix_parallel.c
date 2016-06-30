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

#include "matrix_parallel.h"


/*---------------------------------------------------------------------------
 * Compute matrix bandwidth
 *--------------------------------------------------------------------------*/
long int MATRIX_PARALLEL_bandwidth (MAT* A)
{
	int size = A->n;
	unsigned long int bandwidth = 0;
	
	#pragma omp parallel
	{
		int row, band, max_band;
		
		max_band = 0;
		
		#pragma omp for schedule(static)
		for (row = 0; row < size; ++row)
		{
			band = row - A->JA[A->IA[row]];
			if (band > max_band) max_band = band;
		}
		
		if (max_band > bandwidth)
		{
			#pragma omp critical
			bandwidth = max_band;
		}
	}
	
	return bandwidth;
}


/*---------------------------------------------------------------------------
 * Compute matrix maximum wavefront
 *--------------------------------------------------------------------------*/
long int MATRIX_PARALLEL_wavefront (MAT* A)
{
	int size = A->n;
	int nnz  = A->nz;
	unsigned long int wavefront = 0;
	
	#pragma omp parallel
	{
		int row, slice_row, ja, wfront, max_wfront;
		
		max_wfront = 0; 
		
		#pragma omp for schedule(static)
		for (row = 0; row < size; ++row)
		{
			wfront = 0;
			
			// Computing row-th wavefront 
			for (slice_row = row; slice_row < size-1; ++slice_row) 
			{
				for (ja = A->IA[slice_row]; ja < A->IA[slice_row+1]; ++ja)
				{
					if (A->JA[ja] <= row) 
					{
						++wfront;
						ja = A->IA[slice_row+1];
					}
				}
			}
			
			for (ja = A->IA[slice_row]; ja < nnz; ++ja)
			{
				if (A->JA[ja] <= row) 
				{
					++wfront;
					ja = nnz;
				}
			}
				
// 			printf("Wavefront of row %d: %d\n", row, wfront);fflush(stdout);
			
			if (wfront > max_wfront) max_wfront = wfront;
		}
		
// 		printf("Max wavefront from thread %d: %d\n", omp_get_thread_num(), max_wfront);fflush(stdout);
		
		if (max_wfront > wavefront)
		{
			#pragma omp critical
			wavefront = max_wfront;
		}
	}
	
	return wavefront;
}


/*----------------------------------------------------------------------------
 * Perform operation P*A*P' in CSR format
 *--------------------------------------------------------------------------*/
void MATRIX_PARALLEL_permutation (MAT* A, int* p)
{
	int n, m, nz, n_nz, num_threads, chunk_size, chunk_size_nz;
	MAT*  B;
	ARRAY* a;
	int* q;
	int* counts;
	
	n  = A->n;
	m  = A->m;
	nz = A->nz;  
	B  = malloc(sizeof(MAT));
	
	#pragma omp parallel 
	{
		int i;
		
		#pragma omp sections
		{
			#pragma omp section
			{
				B->n  = n;
				B->m  = m;
				B->nz = nz;
				n_nz  = nz - n;
				num_threads = omp_get_num_threads();
				chunk_size     = num_threads > n ? 1 : n/num_threads;
				chunk_size_nz  = num_threads > n_nz ? 1 : n_nz/num_threads;
			}
			
			#pragma omp section
			B->AA = malloc(nz * sizeof(double));
			
			#pragma omp section
			B->JA = malloc(nz * sizeof(int));
			
			#pragma omp section
			B->D = malloc(n * sizeof(double));
			
			#pragma omp section
			a = malloc (nz*sizeof(ARRAY));
			
			#pragma omp section
			q = malloc (n *sizeof(int));
			
			#pragma omp section
			counts = calloc(n + 1, sizeof(int));
		}

		#pragma omp for schedule(static, chunk_size)
		for (i = 0; i < n; ++i)
		{
			counts[i + 1] = (A->IA[p[i]+1]) - (A->IA[p[i]]);
			q[p[i]] = i;
		}
	}
	
	prefix_sum(counts, &B->IA, n + 1);
	
	#pragma omp parallel
	{
		int i, j, k, index;
		
		#pragma omp for schedule(static, chunk_size)
		for (i = 0; i < n; ++i)
		{
			k = B->IA[i];
			index = A->IA[p[i]];
			
			for (j = B->IA[i]; j < B->IA[i+1]; ++j)
			{
				a[k].arr1 = A->AA[index];
				a[k].arr2 = q[A->JA[index]];
				a[k].arr3 = i;
				++k;
				++index;
			}
		}
		
		#pragma omp single
		qsort(a, nz, sizeof(ARRAY), COMPARE_array);
		
		#pragma omp for schedule(static, chunk_size) nowait
		for (i = 0; i < n; ++i)
		{
			A->AA[i]   = a[i].arr1;
			A->JA[i]   = a[i].arr2;
			A->IA[i+1] = B->IA[i+1];
		}
		
		#pragma omp for schedule(static, chunk_size_nz)
		for (i = n; i < nz; ++i)
		{
			A->AA[i] = a[i].arr1;
			A->JA[i] = a[i].arr2;
		}
		
		#pragma omp sections
		{
			#pragma omp section
			free(a);
			
			#pragma omp section
			free(q);
			
			#pragma omp section
			free(counts);
			
			#pragma omp section
			free(B->AA);
			
			#pragma omp section
			free(B->JA);
			
			#pragma omp section
			free(B->IA);
			
			#pragma omp section
			free(B->D);
		}
		
		#pragma omp single
		free(B);
	}
}