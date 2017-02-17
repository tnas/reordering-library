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
 * Compute the wavefront of each row of matrix A
 * @since 01-07-2016
 *--------------------------------------------------------------------------*/
void inline MATRIX_PARALLEL_wavefront(const MAT* A, int** wf_per_row)
{
	int size, nnz;
	int* ext_IA;
	
	size = A->n;
	nnz  = A->nz;
	ext_IA = calloc(size+1, sizeof(int));
	
	#pragma omp parallel
	{
		int row, ja, cand_col, virtual_row;
		int* computed_colu_per_row;
	
		#pragma omp sections
		{
			#pragma omp section
			*wf_per_row = calloc(size, sizeof(int));
			
			#pragma omp section
			memcpy(ext_IA, A->IA, size*sizeof(int));
			
			#pragma omp section
			ext_IA[size] = nnz;
		}
		
		computed_colu_per_row = calloc(size, sizeof(int));
		
		#pragma omp for schedule(static)
		for (row = 0; row < size; ++row)
		{
			memset(computed_colu_per_row, 0, size * sizeof(int));
			
			for (ja = ext_IA[row]; ja < ext_IA[row+1]; ++ja)
			{
				cand_col = A->JA[ja];
				
				for (virtual_row = 0; virtual_row <= row; ++virtual_row)
				{
					if (cand_col <= virtual_row && computed_colu_per_row[virtual_row] == 0) 
					{
						#pragma omp atomic
						++((*wf_per_row)[virtual_row]);
						
						computed_colu_per_row[virtual_row] = 1;
					}
				}
			}
		}
		
		free(computed_colu_per_row);
		
		#pragma omp single nowait
		free(ext_IA);
	}
}



/*---------------------------------------------------------------------------
 * Compute matrix maximum wavefront
 *--------------------------------------------------------------------------*/
long int MATRIX_PARALLEL_max_wavefront(MAT* A)
{
	int size, row, max_wavefront;
	int* wf_per_row;
	
	size = A->n;
	max_wavefront = 0;
	
	MATRIX_PARALLEL_wavefront(A, &wf_per_row);
	
	#pragma omp parallel for schedule(static) private(row)
	for (row = 0; row < size; ++row)
	{
		if (wf_per_row[row] > max_wavefront)
		{
			#pragma omp critical
			{
				if (wf_per_row[row] > max_wavefront)
					max_wavefront = wf_per_row[row];
			}
		}
		
	}
	
	free(wf_per_row);
	
	return max_wavefront;
	
}



/*---------------------------------------------------------------------------
 * Compute matrix root mean-square (RMS) wavefront
 *--------------------------------------------------------------------------*/
long int MATRIX_PARALLEL_rms_wavefront(MAT* A)
{
	int size, row, rms_wavefront;
	int* wf_per_row;
	
	size = A->n;
	rms_wavefront = 0;
	
	MATRIX_PARALLEL_wavefront(A, &wf_per_row);
	
	#pragma omp parallel for schedule(static) private(row) reduction(+:rms_wavefront)
	for (row = 0; row < size; ++row)
	{
		wf_per_row[row] = pow(wf_per_row[row], 2);
		rms_wavefront += wf_per_row[row];
	}
	
	#pragma omp single nowait
	free(wf_per_row);
	
	return rms_wavefront / size;
}



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
		
		#pragma omp critical
		{
			if (max_band > bandwidth)
				bandwidth = max_band;
		}
	}
	
	return bandwidth;
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