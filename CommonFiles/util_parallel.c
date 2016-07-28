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

#include "util_parallel.h"


int count_nodes_by_level(const int* levels, const int n_nodes, int** counts)
{
	int node, count_thread, max_level, level, id_thread;
	int** local_count;
	int* local_max;
	
	max_level = 0;
	const int num_threads = omp_get_max_threads();
	
	#pragma omp parallel private(node, count_thread, level, id_thread) 
	{
		id_thread = omp_get_thread_num();
		
		#pragma omp single nowait
		local_max = calloc(num_threads, sizeof(int*));
		
		#pragma omp single
		local_count = calloc(num_threads, sizeof(int*));
		
		#pragma omp barrier
		
		#pragma omp for nowait
		for (count_thread = 0; count_thread < num_threads; ++count_thread)
			local_count[count_thread] = calloc(n_nodes, sizeof(int));
		
		#pragma omp for nowait
		for (node = 0; node < n_nodes; ++node)
		{
			++local_count[id_thread][levels[node]];
			local_max[id_thread] = max(local_max[id_thread], levels[node]); 
		}
		
// 		#pragma omp for reduction(max:max_level)
// 		for (count_thread = 0; count_thread < num_threads; ++count_thread)
// 			max_level = max(max_level, local_max[count_thread]);

		#pragma omp critical
		max_level = max(max_level, local_max[id_thread]);
		
		#pragma omp barrier
		
		#pragma omp single	
		max_level += 2;
		
		#pragma omp flush(max_level)
		
		#pragma omp single
		{
			*counts = calloc(max_level, sizeof(int));
			(*counts)[0] = 0;
		}
		
		#pragma omp flush(counts, local_count)
		
		#pragma omp for 
		for (level = 0; level < max_level-1; ++level) 
		{
			for (count_thread = 0; count_thread < num_threads; ++count_thread) 
				(*counts)[level+1] += local_count[count_thread][level];
		}
		
		#pragma omp for			
		for (count_thread = 0; count_thread < num_threads; ++count_thread) 
			free(local_count[count_thread]);
		
		#pragma omp single nowait
		free(local_count);
		
		#pragma omp single nowait
		free(local_max);
	}
	
	return max_level;
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
