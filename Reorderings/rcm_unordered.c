/*----------------------------------------------------------------------------
 * UNORDERED RCM REORDERING SOLVER
 *--------------------------------------------------------------------------*/
#include "../CommonFiles/protos_parallel.h"


int count_nodes_by_level(const int* levels, const int n_nodes, int** counts)
{
	int node, count_thread, max_level, level;
	int** local_count;
	int* local_max;
	
	max_level = 0;
		
	#pragma omp parallel private(node, count_thread, level)
	{
		#pragma omp single nowait
		local_max = calloc(NUM_THREADS, sizeof(int*));
		
		#pragma omp single
		local_count = calloc(NUM_THREADS, sizeof(int*));
		
		#pragma omp barrier
		
		#pragma omp for nowait
		for (count_thread = 0; count_thread < NUM_THREADS; ++count_thread)
			local_count[count_thread] = calloc(n_nodes, sizeof(int));
		
		#pragma omp for nowait
		for (node = 0; node < n_nodes; ++node)
		{
			++local_count[omp_get_thread_num()][levels[node]];
			local_max[omp_get_thread_num()] = max(local_max[omp_get_thread_num()], levels[node]); 
		}
		
		#pragma omp for reduction(max:max_level)
		for (count_thread = 0; count_thread < NUM_THREADS; ++count_thread)
			max_level = max(max_level, local_max[count_thread]);
		
		#pragma omp single	
		++max_level;
		
		#pragma omp single
		*counts = calloc(max_level, sizeof(int));
		
		#pragma omp flush(max_level, counts, local_count)
		
		#pragma omp for 
		for (level = 0; level < max_level; ++level) 
		{
			for (count_thread = 0; count_thread < NUM_THREADS; ++count_thread) 
				(*counts)[level] += local_count[count_thread][level];
		}
		
		#pragma omp for nowait			
		for (count_thread = 0; count_thread < NUM_THREADS; ++count_thread) 
			free(local_count[count_thread]);
		
		#pragma omp single nowait
		free(local_count);
		
		#pragma omp single nowait
		free(local_max);
	}
	
	return max_level;
}

int pow_uint(int base, const int exp)
{
	int count;
	if (exp == 0) return 1;
	for (count = 0; count < exp; ++count)
		base *= base;
	return base;
}

void prefix_sum(const int* counts, int** sums, int max_level)
{
	int level, local_level, chunk_size, index_processors, id_proc, offset_level, 
		count_thread, coef_target_proc, target_proc;
	status_prefix_sum status_ps[NUM_THREADS];
/*	
	*sums = calloc(max_level + 1, sizeof(int));
	(*sums)[0] = 0;*/
	
	max_level = 15;
	*sums = calloc(15, sizeof(int));
	
	int tcounts[15] = { 9, 8, 3, 2, 7, 1, 6, 4, 5, 5, 8, 3, 4, 9, 1}; 
	
	chunk_size = iseven(max_level) ? max_level/NUM_THREADS : max_level/NUM_THREADS + 1;
	offset_level = -chunk_size;
	
	#pragma omp parallel private (level, local_level, id_proc, coef_target_proc, target_proc)
	{
		#pragma omp single
		index_processors = (int) log2f(NUM_THREADS);
		
		#pragma omp critical
		{
			#pragma omp flush (offset_level)
			offset_level += chunk_size;
			level = offset_level;
		}
	
		#pragma omp for
		for (count_thread = 0; count_thread < NUM_THREADS; ++count_thread)
		{	
			(*sums)[level] = tcounts[level];
			id_proc = level / NUM_THREADS;
			
			for (local_level = level+1; local_level < level+chunk_size && local_level < max_level; ++local_level)
				(*sums)[local_level] = (*sums)[local_level-1] + tcounts[local_level];
			
			status_ps[id_proc].initial_prefix_sum = status_ps[id_proc].total_sum = (*sums)[local_level-1];
			
			for (coef_target_proc = 0; coef_target_proc < index_processors; ++coef_target_proc)
			{
				target_proc = id_proc ^ pow_uint(2, coef_target_proc);
				if (target_proc > NUM_THREADS) continue;
				
				status_ps[target_proc].total_sum += status_ps[id_proc].total_sum;
				
// 				if (target_proc < id_proc)
// 					status_ps
			}
		}
	}
}

/*----------------------------------------------------------------------------
 * Unordered RCM reordering from the LEVEL STRUCTURE 
 *--------------------------------------------------------------------------*/
void REORDERING_RCM_parallel(MAT* A, int** perm)
{ 
	int n_nodes, root, count_nodes, max_level;
	int* levels;
	int* counts;
	int* sums;
	
	n_nodes = A->n;
	levels = calloc (n_nodes, sizeof(int));
	
// 	int* g = GRAPH_LS_peripheral (A, &root, &e);
	
	root = 1;
	levels = GRAPH_parallel_fixedpoint_bfs(A, root, levels);
	
	max_level = count_nodes_by_level(levels, n_nodes, &counts);
	
	prefix_sum(counts, &sums, max_level);
	
	/* Reverse order */
// 	for (count_nodes = 0; count_nodes < n_nodes; ++count_nodes) 
// 		(*perm)[n_nodes-1-count_nodes] = tperm[count_nodes]; 
	
	printf("Sums vector: ");
	for (count_nodes = 0; count_nodes < 15; ++count_nodes) 
		printf("%d ", sums[count_nodes]);
	printf("\n");
// 	
// 	printf("******************************************\n");
// 	
// 	for (count_nodes = 0; count_nodes < n_nodes; ++count_nodes) 
// 		if (levels[count_nodes] == 2147483647) printf("Node: %d is in level %d\n", count_nodes, levels[count_nodes]);
// 	printf("******************************************\n");
	
	free(levels);
	free(counts);
	free(sums);
	
		
// 	free(g);
}



