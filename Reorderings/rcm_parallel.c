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
	
	#pragma omp parallel private(node, count_thread, level) num_threads(NUM_THREADS)
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
// 		++max_level;
		max_level += 2;
		
		#pragma omp flush(max_level)
		
		#pragma omp single
		{
			*counts = calloc(max_level, sizeof(int));
			(*counts)[0] = 0;
		}
		
		#pragma omp flush(counts, local_count)
		
		#pragma omp for 
		for (level = 0; level < max_level; ++level) 
		{
			for (count_thread = 0; count_thread < NUM_THREADS; ++count_thread) 
				(*counts)[level+1] += local_count[count_thread][level];
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
			
// 		#pragma omp single
// 		{
// 			printf("First generation of sums: ");fflush(stdout);
// 			for (temp = 0; temp < max_level; ++temp)
// 				printf("%d ", (*sums)[temp]);fflush(stdout);
// 			printf("\n");fflush(stdout);
// 		
// 		}
		
// 		printf("index_processors: %d\n", index_processors);fflush(stdout);	
		
		for (coef_target_proc = 0; coef_target_proc < index_processors; ++coef_target_proc)
		{
			target_proc = id_proc ^ pow_uint(2, coef_target_proc);
			
// 			if (id_proc == 0)
// 			printf("id_proc %d sending to target_proc %d because coef_target_proc: %d and pow_uint is %d\n", id_proc, target_proc, coef_target_proc, pow_uint(2, coef_target_proc)); fflush(stdout);
// 			
// 			if (target_proc < omp_get_num_threads() && target_proc != id_proc)
// 			{
// 				
// 				if (id_proc < target_proc) 
// 				{
// 					#pragma omp critical
// 					{
// 						status_ps[target_proc].last_prefix_sum += status_ps[id_proc].curr_total_sum;
// 						status_ps[target_proc].last_total_sum += status_ps[id_proc].curr_total_sum;
// // 						if (target_proc == 8)
// // 						printf("proc: %d updating proc %d => cpsum: %d, ctsum: %d, lpsum: %d, ltsum: %d\n",
// // 							id_proc, target_proc, status_ps[target_proc].curr_prefix_sum, status_ps[target_proc].curr_total_sum, status_ps[target_proc].last_prefix_sum, status_ps[target_proc].last_total_sum);fflush(stdout);
// 							
// 					}
// 				}
// 				else  
// 				{
// 					#pragma omp critical
// 					{
// 						status_ps[target_proc].last_total_sum += status_ps[id_proc].curr_total_sum;
// // 						if (target_proc == 8)
// // 						printf("proc: %d updated by proc %d => cpsum: %d, ctsum: %d, lpsum: %d, ltsum: %d\n",
// // 							target_proc, id_proc, status_ps[target_proc].curr_prefix_sum, status_ps[target_proc].curr_total_sum, status_ps[target_proc].last_prefix_sum, status_ps[target_proc].last_total_sum);fflush(stdout);
// 							
// 					}
// 				}
// 			}
			
			
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





// void place(MAT* graph, const int source_node, const int* sums, const int max_dist, int** perm, const int* levels)
// {
// 	int level, node, degree, count;
// 	int* read_offset;
// 	int* write_offset;
// 	GRAPH* children;
// 	mem_write_next_level mem_write;
// 	
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
// 
// // 	printf("max_dist: %d\n", max_dist); fflush(stdout);
// 	#pragma omp parallel private (level, node, children, degree, count, mem_write) num_threads(NUM_THREADS)
// 	{
// 		level = omp_get_thread_num();
// 		
// 		while (level < max_dist - 1)
// 		{
// // 			printf("Thread %d running\n", level); fflush(stdout);
// 			
// 			init_mem(&mem_write, sums, level);
// 			
// 			while (read_offset[level] != sums[level+1]) // There are nodes to read
// 			{
// 				#pragma omp flush (write_offset)
// 				
// 				while (read_offset[level] == write_offset[level]) { } // Spin
// 				
// 				
// 				node = (*perm)[read_offset[level]];
// 				++read_offset[level];
// 				
// // 				if (level == 6) printf("Thread %d get node %d\n", level, node); fflush(stdout);
// 				
// 				// Edges of node with dist == level + 1
// 				degree   = GRAPH_degree_per_level(graph, node, levels, level+1);
// 				children = GRAPH_adjacent_per_level(graph, node, levels, level+1);
// 				
// // 				if (level == 6) printf("Thread %d get get %d degree %d\n", level, node, degree); fflush(stdout);
// 				
// 				// Sorting children by degree
// 				qsort(children, degree, sizeof(GRAPH), COMPARE_degr_ASC);
// 				
// 				for (count = 0; count < degree; ++count)
// 				{
// 					if (!is_in_mem(mem_write, children[count].label))
// 					{
// // 						if (level == 6) printf("Thread %d write_offset: %d\n", level, write_offset[level+1]); fflush(stdout);
// 						(*perm)[write_offset[level+1]] = mem_write.mem[mem_write.free_position] = children[count].label;
// // 						if (level == 6) printf("Thread %d writing node %d\n", level, children[count].label); fflush(stdout);
// 						++write_offset[level+1];
// 						++mem_write.free_position;
// // 						if (level == 6) printf("Thread %d finish writing\n", level); fflush(stdout);
// 					}
// 				}
// 				
// 				free(children);
// 			}
// 			
// // 			printf("Thread %d finish\n", level); fflush(stdout);
// 			
// 			free(mem_write.mem);
// 			level += NUM_THREADS;
// 		}
// 		
// // 		printf("Thread %d waiting to die.\n", level-NUM_THREADS < 0 ? level : level-NUM_THREADS); fflush(stdout);
// 	}
// 	
// 	free(read_offset);
// 	free(write_offset);
// }

// void inline init_mem(mem_write_next_level* mem_write, const int* sums, const int level)
// {
// 	int pos;
// 	
// 	mem_write->size = sums[level+2]-sums[level+1];
// 	mem_write->mem = calloc(mem_write->size, sizeof(int));
// 	mem_write->free_position = 0;
// 	
// 	for (pos = 0; pos < mem_write->size; ++pos)
// 		mem_write->mem[pos] = -1;
// }
// 
// 
// 
// inline int is_in_mem(const mem_write_next_level mem_write, const int val)
// {
// 	int pos;
// 	for (pos = 0; pos < mem_write.size; ++pos)
// 		if (mem_write.mem[pos] == val) return 1;
// 	return 0;
// }

// Tentativa usando omp_set e omp_unset lock - ultimo ok
// void place(MAT* graph, const int source_node, const int* sums, const int max_dist, int** perm, const int* levels)
// {
// 	int level, node, degree, count;
// 	int* read_offset;
// 	int* write_offset;
// 	GRAPH* children;
// 	mem_write_next_level mem_write;
// 	omp_lock_t locks[NUM_THREADS];
// 	
// 	#pragma omp parallel sections num_threads(4)
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
// 		
// 		#pragma omp section 
// 		for (count = 0; count < NUM_THREADS; ++count)
// 			omp_init_lock(&locks[count]);
// 	}
// 	
// 	
// 
// // 	printf("max_dist: %d\n", max_dist); fflush(stdout);
// 	#pragma omp parallel private (level, node, children, degree, count, mem_write) num_threads(NUM_THREADS)
// 	{
// 		level = omp_get_thread_num();
// 		
// 		omp_set_lock(&locks[level]);
// 		#pragma omp barrier
// 		
// 		while (level < max_dist - 1)
// 		{
// // 			printf("Thread %d running\n", level); fflush(stdout);
// 			
// 			init_mem(&mem_write, sums, level);
// 			
// 			while (read_offset[level] != sums[level+1]) // There are nodes to read
// 			{
// 				#pragma omp flush (write_offset)
// 				
// 				// Spin
// 				while (read_offset[level] == write_offset[level]) 
// 				{
// // 					printf("Thread %d blocked.\n", level); fflush(stdout);
// 					omp_set_lock(&locks[(level-1)%NUM_THREADS]);
// 				} 
// 				
// 				if (level > 0) omp_unset_lock(&locks[(level-1)%NUM_THREADS]);
// 				
// 				node = (*perm)[read_offset[level]];
// 				++read_offset[level];
// 				
// // 				printf("Thread %d get node %d\n", level, node); fflush(stdout);
// 				
// 				// Edges of node with dist == level + 1
// 				degree   = GRAPH_degree_per_level(graph, node, levels, level+1);
// 				children = GRAPH_adjacent_per_level(graph, node, levels, level+1);
// 				
// // 				printf("Thread %d get %d degree %d\n", level, node, degree); fflush(stdout);
// 				
// 				// Sorting children by degree
// 				qsort(children, degree, sizeof(GRAPH), COMPARE_degr_ASC);
// 				
// 				for (count = 0; count < degree; ++count)
// 				{
// // 					if (!is_in_mem(mem_write, children[count].label))
// // 					{
// // 						printf("Thread %d write_offset: %d\n", level, write_offset[level+1]); fflush(stdout);
// // 						printf("Thread %d getting lock\n", level); fflush(stdout);
// 						omp_test_lock(&locks[level%NUM_THREADS]);
// // 						printf("Thread %d got lock\n", level); fflush(stdout);
// 						(*perm)[write_offset[level+1]] = mem_write.mem[mem_write.free_position] = children[count].label;
// // 						printf("Thread %d writing node %d\n", level, children[count].label); fflush(stdout);
// 						++write_offset[level+1];
// 						++mem_write.free_position;
// 						omp_unset_lock(&locks[level%NUM_THREADS]);
// // 						printf("Thread %d finish writing\n", level); fflush(stdout);
// // 					}
// 				}
// 				
// 				free(children);
// 			}
// 			
// // 			printf("Thread %d finish\n", level); fflush(stdout);
// 			
// 			free(mem_write.mem);
// 			level += NUM_THREADS;
// 		}
// 		
// // 		printf("Thread %d waiting to die.\n", level-NUM_THREADS < 0 ? level : level-NUM_THREADS); fflush(stdout);
// 	}
// 	
// 	#pragma omp parallel num_threads(3)
// 	{
// 		#pragma omp single nowait
// 		free(read_offset);
// 		
// 		#pragma omp single nowait
// 		free(write_offset);
// 		
// 		#pragma omp single
// 		for (count = 0; count < NUM_THREADS; ++count)
// 			omp_destroy_lock(&locks[count]);
// 	}
// }



void place(MAT* graph, const int source_node, const int* sums, const int max_dist, int** perm, const int* levels)
{
	int level, node, degree, count;
	int* read_offset;
	int* write_offset;
	int colors[graph->n];
	GRAPH* children;
// 	mem_write_next_level mem_write;
	omp_lock_t locks[NUM_THREADS];
	
	#pragma omp parallel sections num_threads(NUM_THREADS)
	{
		#pragma omp section 
		{
			read_offset  = calloc(max_dist, sizeof(int));
			bcopy(sums, read_offset, max_dist * sizeof(int));
		}
		
		#pragma omp section 
		{
			write_offset = calloc(max_dist, sizeof(int));
			bcopy(sums, write_offset, max_dist * sizeof(int));
			write_offset[0] = 1;
		}
		
		#pragma omp section
		(*perm)[0] =  source_node;
		
		#pragma omp section 
		for (count = 0; count < NUM_THREADS; ++count)
			omp_init_lock(&locks[count]);
	}
	
	
	#pragma omp parallel for num_threads(NUM_THREADS)
	for (count = 0; count < graph->n; ++count)
		colors[count] = UNREACHED;

// 	printf("max_dist: %d\n", max_dist); fflush(stdout);
	#pragma omp parallel private (level, node, children, degree, count) num_threads(NUM_THREADS)
	{
		level = omp_get_thread_num();
		
		omp_set_lock(&locks[level]);
		#pragma omp barrier
		
		while (level < max_dist - 1)
		{
// 			printf("Thread %d running\n", level); fflush(stdout);
			
// 			init_mem(&mem_write, sums, level);
			
			while (read_offset[level] != sums[level+1]) // There are nodes to read
			{
				#pragma omp flush (write_offset)
				
				// Spin
				while (read_offset[level] == write_offset[level]) 
				{
// 					printf("Thread %d blocked.\n", level); fflush(stdout);
					omp_set_lock(&locks[(level-1)%NUM_THREADS]);
				} 
				
				if (level > 0) omp_unset_lock(&locks[(level-1)%NUM_THREADS]);
				
				node = (*perm)[read_offset[level]];
				++read_offset[level];
				
// 				printf("Thread %d get node %d\n", level, node); fflush(stdout);
				
				// Edges of node with dist == level + 1
				degree   = GRAPH_degree_per_level(graph, node, levels, level+1, colors);
				children = GRAPH_adjacent_per_level(graph, node, levels, level+1, colors);
				
// 				printf("Thread %d get %d degree %d\n", level, node, degree); fflush(stdout);
				
				// Sorting children by degree
				qsort(children, degree, sizeof(GRAPH), COMPARE_degr_ASC);
				
				for (count = 0; count < degree; ++count)
				{
// 					if (!is_in_mem(mem_write, children[count].label))
// 					{
// 						printf("Thread %d write_offset: %d\n", level, write_offset[level+1]); fflush(stdout);
// 						printf("Thread %d getting lock\n", level); fflush(stdout);
						omp_test_lock(&locks[level%NUM_THREADS]);
// 						printf("Thread %d got lock\n", level); fflush(stdout);
// 						(*perm)[write_offset[level+1]] = mem_write.mem[mem_write.free_position] = children[count].label;
						(*perm)[write_offset[level+1]] = children[count].label;
						colors[children[count].label] = LABELED;
// 						printf("Thread %d writing node %d\n", level, children[count].label); fflush(stdout);
						++write_offset[level+1];
// 						++mem_write.free_position;
						omp_unset_lock(&locks[level%NUM_THREADS]);
// 						printf("Thread %d finish writing\n", level); fflush(stdout);
// 					}
				}
				
				free(children);
			}
			
// 			printf("Thread %d finish\n", level); fflush(stdout);
			
// 			free(mem_write.mem);
			level += NUM_THREADS;
		}
		
// 		printf("Thread %d waiting to die.\n", level-NUM_THREADS < 0 ? level : level-NUM_THREADS); fflush(stdout);
	}
	
	#pragma omp parallel sections num_threads(NUM_THREADS)
	{
		#pragma omp section
		free(read_offset);
		
		#pragma omp section
		free(write_offset);
		
		#pragma omp section
		for (count = 0; count < NUM_THREADS; ++count)
			omp_destroy_lock(&locks[count]);
	}
}





/*
void *alloc_perm_vector(void *params)
{
	int node, degree, count;
	mem_write_next_level mem_write;
	GRAPH* children;
	thread_args *args = params;
	
	while (args->level < args->max_dist - 1)
	{
		init_mem(&mem_write, args->sums, args->level);
		
		while (args->read_offset[args->level] != args->sums[args->level+1]) // There are nodes to read
		{

			while (args->read_offset[args->level] == args->write_offset[args->level]) 
			{
				
			}
			
			node = (*(args->perm))[args->read_offset[args->level]];
			++(args->read_offset[args->level]);
			
			// Edges of node with dist == level + 1
			degree   = GRAPH_degree_per_level(args->graph, node, args->levels, args->level+1);
			children = GRAPH_adjacent_per_level(args->graph, node, args->levels, args->level+1);
			
			// Sorting children by degree
			qsort(children, degree, sizeof(GRAPH), COMPARE_degr_ASC);
			
			for (count = 0; count < degree; ++count)
			{
				if (!is_in_mem(mem_write, children[count].label))
				{
					(*(args->perm))[args->write_offset[args->level+1]] = mem_write.mem[mem_write.free_position] = children[count].label;
					++(args->write_offset[args->level+1]);
					++mem_write.free_position;
				}
			}
			
			free(children);
		}
		
		free(mem_write.mem);
		args->level += NUM_THREADS;
	}
	
	return NULL;
}*/



// void place(MAT* graph, const int source_node, const int* sums, const int max_dist, int** perm, const int* levels)
// {
// 	int count, rc;
// 	int* read_offset;
// 	int* write_offset;
// 	pthread_t placers[NUM_THREADS];
// 	thread_args targs;
// 	
// 	#pragma omp parallel sections num_threads(4)
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
// 		
// 		#pragma omp section
// 		{
// 			targs.graph        = graph;
// 			targs.levels       = levels;
// 			targs.max_dist     = max_dist;
// 			targs.perm         = perm;
// 			targs.read_offset  = read_offset;
// 			targs.write_offset = write_offset;
// 			targs.sums         = sums;
// 		}
// 	}
// 
// 	
// 	for (count = 0; count < NUM_THREADS; ++count)
// 	{
// 		targs.level = count;
// 		rc = pthread_create(&placers[count], NULL, alloc_perm_vector, (void *)&targs);
// 		
// 		if (rc)
// 		{
// 			printf("ERROR; return code from pthread_create() is %d\n", rc);
// 			exit(-1);
// 		}
// 	}
// 	
// 	for (count = 0; count < NUM_THREADS; ++count)
// 		pthread_join(placers[count], NULL);
// 	
// 	pthread_exit(NULL);
// 	
// 	#pragma omp parallel num_threads(2)
// 	{
// 		#pragma omp single nowait
// 		free(read_offset);
// 		
// 		#pragma omp single
// 		free(write_offset);
// 	}
// }
// 



/*----------------------------------------------------------------------------
 * Unordered RCM reordering from the LEVEL STRUCTURE 
 *--------------------------------------------------------------------------*/
void Unordered_RCM(MAT* A, int** perm, int root)
{ 
	int n_nodes, max_level, count_nodes;
	int* levels;
	int* counts;
	int* sums;
	int* tperm;
	int* graph_ls;
// 	double time;
	
	n_nodes = A->n;
	
	#pragma omp parallel sections num_threads(NUM_THREADS)
	{
		#pragma omp section
		(*perm) = calloc(n_nodes, sizeof(int));
		
		#pragma omp section
		tperm = calloc(n_nodes, sizeof(int));
		
		#pragma omp section
		levels = calloc(n_nodes, sizeof(int));
	}
	
// 	time = get_time(); 
// 	graph_ls = GRAPH_LS_peripheral_PARALLEL(A, &root, &e);
// 	time = (get_time() - time)/100.0;
// 	printf("Parallel peripheral - Elapsed time: %.6f sec\n\n", time);
	
// 	printf("Iniciando GRAPH_parallel_fixedpoint_bfs\n"); fflush(stdout);
// 	time = get_time(); 
	GRAPH_parallel_fixedpoint_bfs(A, root, &levels);
// 	time = (get_time() - time)/100.0;
// 	printf("Parallel BFS - Elapsed time: %.6f sec\n\n", time);
	
// 	time = get_time();
	max_level = count_nodes_by_level(levels, n_nodes, &counts);
// 	time = (get_time() - time)/100.0;
// 	printf("Parallel Count nodes by level - Elapsed time: %.6f sec\n\n", time);
// 	++max_level;
	
// 	printf("Iniciando prefix_sum\n"); fflush(stdout);
// 	time = get_time();
	prefix_sum(counts, &sums, max_level);
// 	time = (get_time() - time)/100.0;
// 	printf("Parallel Prefix sum - Elapsed time: %.6f sec\n\n", time);
	
// 	printf("Vetor de sums:\n"); fflush(stdout);
// 	for (count_nodes = 0; count_nodes < max_level; ++count_nodes) 
// 		printf("%d ", sums[count_nodes]); fflush(stdout);
// 	printf("\n\n");fflush(stdout);
	
// 	printf("Iniciando place\n"); fflush(stdout);
// 	time = get_time();
	place(A, root, sums, max_level, &tperm, levels);
// 	time = (get_time() - time)/100.0;
// 	printf("Parallel Place - Elapsed time: %.6f sec\n\n", time);
	
	#pragma omp parallel num_threads(NUM_THREADS)
	{
		/* Reverse order */
		#pragma omp for private (count_nodes)
		for (count_nodes = 0; count_nodes < n_nodes; ++count_nodes) 
			(*perm)[n_nodes-1-count_nodes] = tperm[count_nodes]; 
		
		#pragma omp single nowait
		free(levels);
		
		#pragma omp single nowait
		free(counts);
		
		#pragma omp single nowait
		free(sums);
		
		#pragma omp single nowait
		free(tperm);
		
// 		#pragma omp single nowait
// 		free(graph_ls);
	}
	
}



