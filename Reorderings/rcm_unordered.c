/*----------------------------------------------------------------------------
 * UNORDERED RCM REORDERING SOLVER
 *--------------------------------------------------------------------------*/
#include "../CommonFiles/protos_parallel.h"


void count_nodes_by_level(const int* levels, const int n_nodes, int** counts);

/*----------------------------------------------------------------------------
 * Unordered RCM reordering from the LEVEL STRUCTURE 
 *--------------------------------------------------------------------------*/
void REORDERING_RCM_parallel(MAT* A, int** perm)
{ 
	int n_nodes, root, count_nodes;
	int* levels;
	int* counts;
	
	n_nodes = A->n;
	levels = calloc (n_nodes, sizeof(int));
	
// 	int* g = GRAPH_LS_peripheral (A, &root, &e);
	
	root = 1;
	levels = GRAPH_parallel_fixedpoint_bfs(A, root, levels);
	
	count_nodes_by_level(levels, n_nodes, &counts);
	
	/* Reverse order */
// 	for (count_nodes = 0; count_nodes < n_nodes; ++count_nodes) 
// 		(*perm)[n_nodes-1-count_nodes] = tperm[count_nodes]; 
	
	
	for (count_nodes = 0; count_nodes < n_nodes; ++count_nodes) 
		printf("Node: %d is in level %d\n", count_nodes, levels[count_nodes]);
	
	printf("******************************************\n");
	
	for (count_nodes = 0; count_nodes < n_nodes; ++count_nodes) 
		if (levels[count_nodes] == 2147483647) printf("Node: %d is in level %d\n", count_nodes, levels[count_nodes]);
	printf("******************************************\n");
	
	free(levels);
	free(counts);
	
		
// 	free(g);
}


void count_nodes_by_level(const int* levels, const int n_nodes, int** counts)
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
}
