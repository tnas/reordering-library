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

#include "util.h"

int COMPARE_int_ASC (const void * a, const void * b)
{ 
	return ((int*)a) - ((int*)b);
}


int get_random_integer(int max)
{
	srand((unsigned int) time(NULL));
	return rand() % (max + 1);
}


int pow_uint(int base, const int exp)
{
	int count, val_pow;
	if (exp == 0) return 1;
	val_pow = base;
	for (count = 1; count < exp; ++count)
		val_pow *= base;
	return val_pow;
}

double get_time ()
{
	struct timeval tv; gettimeofday(&tv, 0);
	return (double)(tv.tv_sec * 100.0 + tv.tv_usec / 10000.0);
}



/*---------------------------------------------*/
/*---WRITING OUTPUT----------------------------*/
/*---------------------------------------------*/
void write_output_after(const MAT *A)
{
	int i,j;
	
	FILE *f = fopen("output_after.mtx", "w");
	fprintf(f,"%d %d %d\n", A->n, A->m, A->nz);
	for (i = 0; i < A->n; ++i)
	{
		for (j = A->IA[i]; j <= A->IA[i+1]-1; j++)
			fprintf (f,"%d %d %le\n",A->JA[j],i,A->AA[j]); 
	}
	fclose(f);
	FILE *g = fopen("plot_after.plt", "w");
	fprintf(g,"#!/bin/sh\n");
	fprintf(g,"set xrange [0:%d]\n", A->n);
	fprintf(g,"set yrange [0:%d] reverse\n", A->n);
	fprintf(g,"set size square\n");
	fprintf(g,"set terminal png\n");
	fprintf(g,"set output \"output_after.png\"\n");
	fprintf(g,"plot \"output_after.mtx\" using 1:2 notitle with dots linecolor rgb \"#000000\"\n");
	fclose(g);
}



void write_output_before(const MAT *A)
{
	int i,j;
	
	FILE *f = fopen("output_before.mtx", "w");
	fprintf(f,"%d %d %d\n", A->n, A->m, A->nz);
	for (i = 0; i < A->n; ++i)
	{
		for (j = A->IA[i]; j <= A->IA[i+1]-1; j++)
			fprintf (f,"%d %d %le\n",A->JA[j],i,A->AA[j]); 
	}
	fclose(f);
	FILE *g = fopen("plot_before.plt", "w");
	fprintf(g,"#!/bin/sh\n");
	fprintf(g,"set xrange [0:%d]\n", A->n);
	fprintf(g,"set yrange [0:%d] reverse\n", A->n);
	fprintf(g,"set size square\n");
	fprintf(g,"set terminal png\n");
	fprintf(g,"set output \"output_before.png\"\n");
	fprintf(g,"plot \"output_before.mtx\" using 1:2 notitle with dots linecolor rgb \"#000000\"\n");
	fclose(g);
}



void prefix_sum(const int* counts, int** sums, const int max_level)
{
	int level, local_level, chunk_size, index_processors, id_proc, offset_level, 
		coef_target_proc, target_proc, count_thread;
		
	const int num_threads = omp_get_max_threads();
	status_prefix_sum status_ps[num_threads];
	
	if (num_threads > max_level)
	{
		chunk_size = 1;
		omp_set_num_threads(max_level);
	} 
	else
	{
		chunk_size = isdivisor(num_threads, max_level) ? max_level/num_threads : max_level/num_threads + 1;	
		
		for (count_thread = 1, offset_level = 0; 
		     count_thread <= num_threads && offset_level < max_level; 
		    ++count_thread, offset_level += chunk_size);
			
		omp_set_num_threads(count_thread-1);
	}
	
	offset_level = -chunk_size;
	count_thread = 0;
	
	#pragma omp parallel private (level, local_level, id_proc, coef_target_proc, target_proc)
	{
		#pragma omp single nowait
		index_processors = ceil(log2(omp_get_num_threads()));
		
		#pragma omp single
		*sums = calloc(max_level, sizeof(int));
		
		#pragma omp critical
		{
			#pragma omp flush (offset_level, count_thread)
			
			offset_level += chunk_size;
			level = offset_level;
			id_proc = count_thread++;
		}

		
		(*sums)[level] = counts[level];

		
		for (local_level = level+1; local_level < level+chunk_size && local_level < max_level; ++local_level)
			(*sums)[local_level] = (*sums)[local_level-1] + counts[local_level];
		
		
		status_ps[id_proc].initial_prefix_sum 
			= status_ps[id_proc].curr_prefix_sum 
			= status_ps[id_proc].curr_total_sum 
			= status_ps[id_proc].last_prefix_sum
			= status_ps[id_proc].last_total_sum
			= (*sums)[local_level-1];
		
		#pragma omp barrier
			
		for (coef_target_proc = 0; coef_target_proc < index_processors; ++coef_target_proc)
		{
			target_proc = id_proc ^ pow_uint(2, coef_target_proc);
			
			if (target_proc < omp_get_num_threads() && target_proc != id_proc)
			{
				if (id_proc < target_proc) {
					status_ps[target_proc].last_prefix_sum += status_ps[id_proc].curr_total_sum;
					status_ps[target_proc].last_total_sum += status_ps[id_proc].curr_total_sum;
				}
				else  
				{
					status_ps[target_proc].last_total_sum += status_ps[id_proc].curr_total_sum;
				}
			}
			
			#pragma omp barrier
			
			status_ps[id_proc].curr_prefix_sum = status_ps[id_proc].last_prefix_sum;
			status_ps[id_proc].curr_total_sum  = status_ps[id_proc].last_total_sum;
			
			#pragma omp barrier
		}
		
		status_ps[id_proc].last_prefix_sum -= status_ps[id_proc].initial_prefix_sum;
		
		#pragma omp barrier
		
		level = id_proc * chunk_size;
		
		for (local_level = level; local_level < level+chunk_size && local_level < max_level; ++local_level)
			(*sums)[local_level] += status_ps[id_proc].last_prefix_sum;
	}
}

