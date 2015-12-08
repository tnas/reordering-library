/*----------------------------------------------------------------------------
 * UNORDERED RCM REORDERING SOLVER
 *--------------------------------------------------------------------------*/
// #include "../CommonFiles/protos_parallel.h"
#include "../UnitTests/rcm_parallel_test.h"


int count_nodes_by_level(const int* levels, const int n_nodes, int** counts)
{
	int node, count_thread, max_level, level;
	int** local_count;
	int* local_max;
	
	max_level = 0;
	omp_set_num_threads(NUM_THREADS);
	
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
		
		#pragma omp flush(max_level)
		
		#pragma omp single
		*counts = calloc(max_level, sizeof(int));
		
		#pragma omp flush(counts, local_count)
		
		#pragma omp for 
		for (level = 0; level < max_level; ++level) 
		{
			for (count_thread = 0; count_thread < NUM_THREADS; ++count_thread) 
				(*counts)[level] += local_count[count_thread][level];
		}
		
		#pragma omp for			
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
	for (count = 1; count < exp; ++count)
		base *= base;
	return base;
}

void prefix_sum(const int* counts, int** sums, const int max_level)
{
	int level, local_level, chunk_size, index_processors, id_proc, offset_level, 
		coef_target_proc, target_proc, count_thread;
	status_prefix_sum status_ps[NUM_THREADS];
	
	if (NUM_THREADS > max_level)
	{
		chunk_size = 1;
		omp_set_num_threads(max_level);
	} 
	else
	{
		chunk_size = isdivisor(NUM_THREADS, max_level) ? max_level/NUM_THREADS : max_level/NUM_THREADS + 1;	
		
		for (count_thread = 1, offset_level = 0; 
		     count_thread <= NUM_THREADS && offset_level < max_level; 
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
			
			if (target_proc < NUM_THREADS && target_proc != id_proc)
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


void inline init_mem(mem_write_next_level* mem_write, const int* sums, const int level)
{
	int pos;
	
	mem_write->size = sums[level+2]-sums[level+1];
	mem_write->mem = calloc(mem_write->size, sizeof(int));
	mem_write->free_position = 0;
	
	for (pos = 0; pos < mem_write->size; ++pos)
		mem_write->mem[pos] = -1;
}

inline int is_in_mem(const mem_write_next_level mem_write, const int val)
{
	int pos;
	for (pos = 0; pos < mem_write.size; ++pos)
		if (mem_write.mem[pos] == val) return 1;
	return 0;
}


void place(MAT* graph, const int source_node, const int* sums, const int max_dist, int** perm, const int* levels)
{
	int level, node, degree, count;
	int* read_offset;
	int* write_offset;
	GRAPH* children;
	mem_write_next_level mem_write;
	
// 	#pragma omp parallel sections num_threads(3)
// 	{
// 		#pragma omp section
// 		{
// 			read_offset  = calloc(max_dist, sizeof(int));
// 			bcopy(sums, read_offset, max_dist * sizeof(int));
// 		}
// 		
// 		#pragma omp section
// 		{
// 			write_offset = calloc(max_dist, sizeof(int));
// 			bcopy(sums, write_offset, max_dist * sizeof(int));
// 			write_offset[0] = 1;
// 		}
// 		
// 		#pragma omp section
// 		(*perm)[0] =  source_node;
// 	}
	
	read_offset  = calloc(max_dist, sizeof(int));
	bcopy(sums, read_offset, max_dist * sizeof(int));
	
	write_offset = calloc(max_dist, sizeof(int));
	bcopy(sums, write_offset, max_dist * sizeof(int));
	write_offset[0] = 1;
	
	(*perm)[0] =  source_node;
	
	omp_set_num_threads(max_dist - 2);
	
	#pragma omp parallel private (level, node, children, degree, count, mem_write)
	{
		level = omp_get_thread_num();
		
		init_mem(&mem_write, sums, level);
		
		while (read_offset[level] != sums[level+1]) // There are nodes to read
		{
			#pragma omp flush (write_offset)
			
			while (read_offset[level] == write_offset[level]) { } // Spin
			
			node = (*perm)[read_offset[level]];
			++read_offset[level];
			
			// Edges of node with dist == level + 1
			degree = GRAPH_degree_per_level(graph, node, levels, level+1);
			children = GRAPH_adjacent_per_level(graph, node, levels, level+1);
			
			// Sorting children by degree
			qsort(children, degree, sizeof(GRAPH), COMPARE_degr_ASC);
			
			for (count = 0; count < degree; ++count)
			{
				if (!is_in_mem(mem_write, children[count].label))
				{
					
					(*perm)[write_offset[level+1]] = mem_write.mem[mem_write.free_position] = children[count].label;
					++write_offset[level+1];
					++mem_write.free_position;
				}
			}
			
			free(children);
		}
		
		free(mem_write.mem);
	}
	
	free(read_offset);
	free(write_offset);
}


/*----------------------------------------------------------------------------
 * Unordered RCM reordering from the LEVEL STRUCTURE 
 *--------------------------------------------------------------------------*/
void Unordered_RCM(MAT* A, int** perm)
{ 
	int n_nodes, root, max_level, index, count_nodes;
	int* levels;
	int* counts;
	int* sums;
	int* tcounts;
	
	n_nodes = A->n;
	levels = calloc(n_nodes, sizeof(int));
	(*perm) = calloc(n_nodes, sizeof(int));
	
// 	int* g = GRAPH_LS_peripheral (A, &root, &e);
	
	root = 1;
	levels = GRAPH_parallel_fixedpoint_bfs(A, root, levels);
	
	max_level = count_nodes_by_level(levels, n_nodes, &counts);
	++max_level;
	
	tcounts = calloc(max_level, sizeof(int));
	tcounts[0] = 0;
	
	#pragma omp parallel for private (index)
	for (index = 1; index < max_level; ++index)
		tcounts[index] = counts[index-1];
	
	printf("Vetor de tcounts:\n"); fflush(stdout);
	for (count_nodes = 0; count_nodes < max_level; ++count_nodes) 
		printf("%d ", tcounts[count_nodes]); fflush(stdout);
	printf("\n\n");fflush(stdout);
	
	prefix_sum(tcounts, &sums, max_level);
	
	printf("Vetor de sums:\n"); fflush(stdout);
	for (count_nodes = 0; count_nodes < max_level; ++count_nodes) 
		printf("%d ", sums[count_nodes]); fflush(stdout);
	printf("\n\n");fflush(stdout);
	
	test_prefix_sum_parallel_serial(tcounts, max_level);
	
	place(A, root, sums, max_level, perm, levels);
	
	printf("Vetor de permutação:\n"); fflush(stdout);
	for (count_nodes = 0; count_nodes < n_nodes; ++count_nodes) 
		printf("%d ", (*perm)[count_nodes]); fflush(stdout);
	printf("\n\n"); fflush(stdout);
	
	/* Reverse order */
// 	for (count_nodes = 0; count_nodes < n_nodes; ++count_nodes) 
// 		(*perm)[n_nodes-1-count_nodes] = tperm[count_nodes]; 
	
	
	free(levels);
	free(counts);
	free(tcounts);
	free(sums);
	
		
// 	free(g);
}



