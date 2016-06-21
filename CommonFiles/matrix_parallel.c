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
#include "protos_parallel.h"
#include "matrix_parallel.h"

/*----------------------------------------------------------------------------
 * Perform the operation P*A*P' in CSR format
 *--------------------------------------------------------------------------*/
void MATRIX_PARALLEL_permutation (MAT* A, int* p)
{
	int n, m, nz, num_threads, chunk_size, offset;
	MAT*  B;
	ARRAY* a;
	int* q;
	
	n   = A->n;
	m   = A->m;
	nz  = A->nz;  
	B = malloc(sizeof(MAT));
	
	#pragma omp parallel 
	{
		int i, j, k;
		
		#pragma omp sections
		{
			#pragma omp section
			{
				B->n  = n;
				B->m  = m;
				B->nz = nz;
				num_threads = omp_get_num_threads();
				chunk_size = n/num_threads;
				offset = 0;
			}
			
			#pragma omp section
			B->AA = malloc(nz * sizeof(double));
			
			#pragma omp section
			B->JA = malloc(nz * sizeof(int));
			
			#pragma omp section
			{
				B->IA = malloc((n+1) * sizeof(int));
				B->IA[0] = 0;
			}
			
			#pragma omp section
			B->D = malloc(n * sizeof(double));
			
			#pragma omp section
			a = malloc (nz*sizeof(ARRAY));
			
			#pragma omp section
			q = malloc (n *sizeof(int));
		}
		
		#pragma omp critical
		{
			k = offset;
			offset += chunk_size;
		}
		
		#pragma omp for schedule(static, chunk_size)
		for (i = 0; i < n; ++i)
		{
			for (j = A->IA[p[i]]; j <= A->IA[p[i]+1] - 1; ++j)
			{
				B->AA[k] = A->AA[j];
				B->JA[k] = A->JA[j];
				k  = k + 1;
			}
			
			B->IA[i+1] = k;   
			
			q[p[i]] = i;
		}
		
		//TODO: to parallelize after here!
		
		k = 0;
	
		for (i = 0; i < n; ++i)
		{
			for (j = B->IA[i]; j <= B->IA[i+1] - 1; ++j)
			{
				a[k].arr1 = B->AA[j];
				a[k].arr2 = q[B->JA[j]];
				a[k].arr3 = i;
				k = k + 1;
			}
			
			A->IA[i+1] = k;    
		}
		
		qsort(a,nz,sizeof(ARRAY),COMPARE_array);
		
		for (i = 0; i < nz; ++i)
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
			MATRIX_clean(B);
		}
	}
}