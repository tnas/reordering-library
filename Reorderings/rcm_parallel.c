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
	const int num_threads = omp_get_max_threads();
	
	#pragma omp parallel private(node, count_thread, level) 
	{
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
			++local_count[omp_get_thread_num()][levels[node]];
			local_max[omp_get_thread_num()] = max(local_max[omp_get_thread_num()], levels[node]); 
		}
		
		#pragma omp for reduction(max:max_level)
		for (count_thread = 0; count_thread < num_threads; ++count_thread)
			max_level = max(max_level, local_max[count_thread]);
		
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



void place(MAT* graph, const int source_node, const int* sums, const int max_dist, int** perm, const int* levels)
{
	int level, node, degree, count;
	int* read_offset;
	int* write_offset;
	int colors[graph->n];
	GRAPH* children;
	const int num_threads = omp_get_max_threads();
	omp_lock_t locks[num_threads];
	
	#pragma omp parallel sections 
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
		for (count = 0; count < num_threads; ++count)
			omp_init_lock(&locks[count]);
	}
	
	
	#pragma omp parallel for num_threads(num_threads)
	for (count = 0; count < graph->n; ++count)
		colors[count] = UNREACHED;

	#pragma omp parallel private (level, node, children, degree, count) 
	{
		level = omp_get_thread_num();
		
		omp_set_lock(&locks[level]);
		#pragma omp barrier
		
		while (level < max_dist - 1)
		{
			
			while (read_offset[level] != sums[level+1]) // There are nodes to read
			{
				#pragma omp flush (write_offset)
				
				// Spin
				while (read_offset[level] == write_offset[level]) 
				{
					omp_set_lock(&locks[(level-1)%num_threads]);
				} 
				
				if (level > 0) omp_unset_lock(&locks[(level-1)%num_threads]);
				
				node = (*perm)[read_offset[level]];
				++read_offset[level];
				
				// Edges of node with dist == level + 1
				degree   = GRAPH_degree_per_level(graph, node, levels, level+1, colors);
				children = GRAPH_adjacent_per_level(graph, node, levels, level+1, colors);
				
				// Sorting children by degree
				qsort(children, degree, sizeof(GRAPH), COMPARE_degr_ASC);
				
				for (count = 0; count < degree; ++count)
				{
					omp_test_lock(&locks[level%num_threads]);
					
					(*perm)[write_offset[level+1]] = children[count].label;
					colors[children[count].label] = LABELED;
					++write_offset[level+1];
					
					omp_unset_lock(&locks[level%num_threads]);
				}
				
				free(children);
			}
			
			level += num_threads;
		}
	}
	
	#pragma omp parallel sections 
	{
		#pragma omp section
		free(read_offset);
		
		#pragma omp section
		free(write_offset);
		
		#pragma omp section
		for (count = 0; count < num_threads; ++count)
			omp_destroy_lock(&locks[count]);
	}
}



/*----------------------------------------------------------------------------
 * Unordered RCM reordering from the LEVEL STRUCTURE 
 *--------------------------------------------------------------------------*/
void Unordered_RCM(MAT* A, int** perm, int root, const float percent_chunk)
{ 
	int n_nodes, max_level, count_nodes;
	int* levels;
	int* counts;
	int* sums;
	int* tperm;
	
	n_nodes = A->n;
	
	#pragma omp parallel sections
	{
		#pragma omp section
		(*perm) = calloc(n_nodes, sizeof(int));
		
		#pragma omp section
		tperm = calloc(n_nodes, sizeof(int));
		
		#pragma omp section
		levels = calloc(n_nodes, sizeof(int));
	}
	
	GRAPH_parallel_fixedpoint_bfs(A, root, &levels, percent_chunk);
	
	max_level = count_nodes_by_level(levels, n_nodes, &counts);
	
	prefix_sum(counts, &sums, max_level);
	
	place(A, root, sums, max_level, &tperm, levels);
	
	#pragma omp parallel 
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
	}
}



void Leveled_RCM(MAT* mat, int** perm, int root) 
{
	int graph_size, perm_size, node, perm_offset, size_children, size_offset, node_offset;
	GRAPH* graph;
	int* counts;
	int* psum;
	int* parent_index;
	LIST* children;
	LIST* ch_pointer;
	
	graph_size = mat->n;
	perm_size  = perm_offset = 0;
	
	#pragma omp parallel
	{
		#pragma omp single nowait
		*perm  = calloc(graph_size, sizeof(int));
		
		#pragma omp single 
		graph  = calloc(graph_size, sizeof(GRAPH));
		
		#pragma omp for private(node)
		for (node = 0; node < graph_size; ++node)
		{
			graph[node].distance = INFINITY_LEVEL;
			graph[node].parent   = -1;
			graph[node].chnum    = 0;
			graph[node].label    = node;
			graph[node].degree   = GRAPH_degree(mat, node);
		}
		
		#pragma omp single nowait
		graph[root].distance = 0;
		
		#pragma omp single nowait
		graph[root].parent = ORPHAN_NODE;
		
		#pragma omp single nowait
		graph[root].chnum = 0;
		
		#pragma omp single nowait
		(*perm)[0] = root;
		
		#pragma omp single nowait
		perm_size++;
	}
	
	while (perm_size < graph_size)
	{
		#pragma omp parallel
		{
			int n_par, n_ch, degree;
			int* neighbors;
			
			// ******************************
			// Step 1: Expansion
			// ******************************
			
			#pragma omp single nowait
			children = NULL;
			
// 			#pragma omp single nowait
// 			printf("[offset, size] = [%d, %d]\n", perm_offset, perm_size);fflush(stdout);
			
			#pragma omp for schedule(static) ordered
			for (n_par = perm_offset; n_par < perm_size; ++n_par)
			{
				degree    = graph[(*perm)[n_par]].degree;
				neighbors = GRAPH_adjacent(mat, (*perm)[n_par]);
				
// 				int i;
// 				printf("Processing node %d with degree %d and children: ", (*perm)[n_par], degree);fflush(stdout);
// 				for (i = 0; i< degree; i++) printf("%d ", neighbors[i]); fflush(stdout);
// 				printf("\n");fflush(stdout);
				
				for (n_ch = 0; n_ch < degree; ++n_ch)
				{
					if (graph[neighbors[n_ch]].parent == ORPHAN_NODE) 
					{
						#pragma omp critical
						graph[neighbors[n_ch]].parent = (*perm)[n_par];
					}
					
					
					if (graph[neighbors[n_ch]].distance > 
						graph[(*perm)[n_par]].distance)
					{
						if (graph[neighbors[n_ch]].distance > 
							graph[(*perm)[n_par]].distance + 1)
						{
							#pragma omp ordered
							{
								graph[neighbors[n_ch]].distance =
									graph[(*perm)[n_par]].distance + 1;
								children = LIST_insert_IF_NOT_EXIST(children, neighbors[n_ch]);
// 								printf("Adding child %d\n", neighbors[n_ch]);fflush(stdout);
							}
						}
						
						if (graph[(*perm)[n_par]].distance < graph[graph[neighbors[n_ch]].parent].distance)
						{
							#pragma omp critical
							graph[neighbors[n_ch]].parent = (*perm)[n_par];
						}
					}
				}
				
				free(neighbors);
			}

		}
		
		// ******************************
		// Step 2: Reduction
		// ******************************
		
		size_children = children->size;
		ch_pointer    = children;
		
// 		printf("size children: %d\n", size_children);fflush(stdout);
		
		while (ch_pointer != NULL)
		{
			graph[graph[ch_pointer->data].parent].chnum++;
			ch_pointer = ch_pointer->next;
		}
		
		size_offset = perm_size - perm_offset;
		counts = calloc(size_offset+1, sizeof(int));
		counts[0] = perm_size;
		node_offset = 0;
		
		for (node = 0; node < size_offset; ++node)
		{
			if (graph[(*perm)[node+perm_offset]].chnum > 0)
			{
				++node_offset;
				counts[node_offset] = graph[(*perm)[node+perm_offset]].chnum;
				graph[(*perm)[node+perm_offset]].index = node_offset - 1;
// 				printf("parent %d with index %d in psum\n", graph[(*perm)[node+perm_offset]].label, node_offset-1);fflush(stdout);
			}
// 			else
// 			{
// 				printf("parent %d with chnum = %d\n", graph[(*perm)[node+perm_offset]].label, graph[(*perm)[node+perm_offset]].chnum);fflush(stdout);
// 			}
		}
		
// 		printf("Counts vector: ");
// 		int j;
// 		for (j = 0; j < size_offset+1; ++j)
// 			printf("%d ", counts[j]);
// 		printf("\n");fflush(stdout);
		
		// *********************
		// Step 3: Prefix sum
		// *********************
		prefix_sum(counts, &psum, size_offset+1);
		
// 		printf("Sums vector: ");
// 			int k;
// 			for (k = 0; k < size_offset+1; ++k)
// 				printf("%d ", psum[k]);
// 			printf("\n");fflush(stdout);
		
		#pragma omp parallel
		{
			int child, ch, num_children, index_children;
			GRAPH* p_children;
			
			#pragma omp single
			{
				parent_index  = calloc(size_offset+1, sizeof(int));
				bcopy(psum, parent_index, (size_offset+1) * sizeof(int));
			}
			
			// *********************
			// Step 4: Placement
			// *********************
			while (children != NULL)
			{
				child = -1;
				
				#pragma omp critical
				{
					if (children != NULL)
					{
						child    = LIST_first(children);
						children = LIST_remove(children, child);
// 						printf("thread %d get child %d\n", omp_get_thread_num(), child);fflush(stdout);
					}
				}
				
				if (child != -1)
				{
					#pragma omp critical
					{
// 						printf("thread %d set child %d of parent %d at position %d\n", omp_get_thread_num(), child, graph[child].parent, parent_index[graph[graph[child].parent].index]);fflush(stdout);
						(*perm)[parent_index[graph[graph[child].parent].index]++] = child;
					}
					
					if (parent_index[graph[graph[child].parent].index] == psum[graph[graph[child].parent].index + 1])
					{
						num_children   = graph[graph[child].parent].chnum;
						index_children = psum[graph[graph[child].parent].index];
						
						if (num_children > 1)
						{
							// Sorting children by degree
// 							printf("Thread %d Sorting from position %d to %d\n", omp_get_thread_num(), index_children, psum[graph[graph[child].parent].index + 1]);fflush(stdout);
							
							p_children = calloc(num_children, sizeof(GRAPH));
							
							for (ch = 0; ch < num_children; ++ch)
							{
								p_children[ch].label  = (*perm)[index_children + ch];
								p_children[ch].degree = graph[(*perm)[index_children + ch]].degree;
							}
							
							qsort(&((*perm)[index_children]), num_children, sizeof(int), COMPARE_degr_ASC);
							
							for (ch = 0; ch < num_children; ++ch)
								(*perm)[index_children + ch] = p_children[ch].label;
							
							free(p_children);
						}
					}
				}
			}
			
			#pragma omp barrier
			
			#pragma omp single
			perm_size += size_children;
			
			#pragma omp single
			perm_offset = perm_size - size_children;
			
			#pragma omp single nowait
			free(counts);

			#pragma omp single nowait
			free(psum);
			
			#pragma omp single nowait
			free(parent_index);
			
// 			#pragma omp single
// 			printf("------------------------------------\n");fflush(stdout);
		}
		
// 		printf("offset, permsize: [%d, %d]\n", perm_offset, perm_size);fflush(stdout);
// 		printf("Permutation vector: ");fflush(stdout);
// 		int i;
// 		for (i = 0; i < perm_size; ++i)
// 			printf("%d ", (*perm)[i]);
// 		printf("\n");fflush(stdout);
	}
	
	free(graph);
}


void Leveled_RCM_v1(MAT* mat, int** perm, int root) 
{
	int graph_size, perm_size, perm_offset, chunk_size, total_size_children, num_threads;
	GRAPH* graph;
	
	graph_size = mat->n;
	perm_size  = perm_offset = total_size_children = 0;
	
	#pragma omp parallel
	{
		int node, n_par, n_ch, degree, child, size_children_per_parent, 
			size_children, th_id, th_ini;
		int* neighbors;
		LIST* children;
		GRAPH* children_per_parent;
		
		#pragma omp single nowait
		num_threads = omp_get_num_threads();
		
		#pragma omp single nowait
		*perm  = calloc(graph_size, sizeof(int));
		
		#pragma omp single 
		graph  = calloc(graph_size, sizeof(GRAPH));
		
		#pragma omp for
		for (node = 0; node < graph_size; ++node)
		{
			graph[node].label    = node;
			graph[node].distance = INFINITY_LEVEL;
			graph[node].parent   = ORPHAN_NODE;
			graph[node].status   = UNREACHED;
			graph[node].degree   = GRAPH_degree(mat, node);
		}
		
		#pragma omp single nowait
		graph[root].distance = 0;
		
		#pragma omp single nowait
		graph[root].status = LABELED;
		
		#pragma omp single nowait
		(*perm)[0] = root;
		
		#pragma omp single nowait
		perm_size++;
		
		#pragma omp barrier
	
		while (perm_size < graph_size)
		{
			children      = NULL;
			size_children = 0;
			
			#pragma omp single 
			chunk_size = ceil((float) (perm_size - perm_offset) / num_threads);
			
			#pragma omp for schedule(static, chunk_size)
			for (n_par = perm_offset; n_par < perm_size; ++n_par)
			{
// 				printf("Thread %d Processing n_par: %d\n", omp_get_thread_num(), n_par);fflush(stdout);
				degree    = graph[(*perm)[n_par]].degree;
				neighbors = GRAPH_adjacent(mat, (*perm)[n_par]);
				
				children_per_parent      = calloc(degree, sizeof(GRAPH));
				size_children_per_parent = 0;
				
// 				int i;
// 				printf("Thread %d Processing node %d with degree %d and children: ", omp_get_thread_num(), (*perm)[n_par], degree);fflush(stdout);
// 				for (i = 0; i< degree; i++) printf("%d ", neighbors[i]); fflush(stdout);
// 				printf("\n");fflush(stdout);
				
				for (n_ch = 0; n_ch < degree; ++n_ch)
				{
					if (graph[neighbors[n_ch]].parent == ORPHAN_NODE)
						graph[neighbors[n_ch]].parent = (*perm)[n_par];
					
					if (graph[neighbors[n_ch]].distance > graph[(*perm)[n_par]].distance)
					{
						if (graph[neighbors[n_ch]].distance > graph[(*perm)[n_par]].distance + 1)
						{
							#pragma omp critical
							{
								graph[neighbors[n_ch]].distance = graph[(*perm)[n_par]].distance + 1;
									
								if (graph[neighbors[n_ch]].status == UNREACHED)
								{
// 									children = LIST_insert_IF_NOT_EXIST(children, neighbors[n_ch]);
									children_per_parent[size_children_per_parent++] = graph[neighbors[n_ch]]; 
									graph[neighbors[n_ch]].status = LABELED;
								}
							}
						}
						
						if (graph[(*perm)[n_par]].distance < graph[graph[neighbors[n_ch]].parent].distance)
						{
							#pragma omp critical
							graph[neighbors[n_ch]].parent = (*perm)[n_par];
						}
					}
				}
				
				if (size_children_per_parent > 0)
				{
					if (size_children_per_parent > 1) {
// 						printf("Thread %d Sorting position %d to %d\n", omp_get_thread_num(), 0, size_children_per_parent);
// 						printf("Children vector Before: ");fflush(stdout);
// 						int i;
// 						for (i = 0; i < size_children_per_parent; ++i)
// 							printf("%d ", children_per_parent[i].label);
// 						printf("\n");fflush(stdout);
						qsort(children_per_parent, size_children_per_parent, sizeof(GRAPH), COMPARE_degr_ASC);
// 						printf("Children vector After: ");fflush(stdout);
// 						for (i = 0; i < size_children_per_parent; ++i)
// 							printf("%d ", children_per_parent[i].label);
// 						printf("\n");fflush(stdout);
					}
					
					for (node = 0; node < size_children_per_parent; ++node)
						children = LIST_insert_IF_NOT_EXIST(children, children_per_parent[node].label);
				}
				
				free(neighbors);
				free(children_per_parent);
			}
			
			if (children != NULL) size_children = children->size;
			
			#pragma omp for schedule(static, 1)
			for (th_id = 0; th_id < num_threads; ++th_id)
			{
// 				printf("Thread %d Processing th_id: %d\n", omp_get_thread_num(), th_id);fflush(stdout);
				
				#pragma omp critical
				{
// 					printf("Thread %d setting perm vector\n", omp_get_thread_num());fflush(stdout);
					th_ini     = perm_size;
					perm_size += size_children;
					total_size_children += size_children;
				}
				
				while (children != NULL)
				{
					child    = LIST_first(children);
					children = LIST_remove(children, child);
// 					printf("Thread %d Setting child %d(%d) at position %d\n", omp_get_thread_num(), child, GRAPH_degree(mat, child), th_ini);fflush(stdout);
					(*perm)[th_ini++] = child;
				}
				
			}
			
			#pragma omp single
			perm_offset = perm_size - total_size_children;
			
// 			#pragma omp single
// 			printf("-----------------------------\n");fflush(stdout);
		}
	}
		
	free(graph);
/*		
	printf("Permutation vector: ");fflush(stdout);
	int i;
	for (i = 0; i < perm_size; ++i)
		printf("%d ", (*perm)[i]);
	printf("\n");fflush(stdout);*/
}


void Leveled_RCM_v2(MAT* mat, int** perm, int root) 
{
	int graph_size, perm_size, perm_offset, chunk_size, total_size_children, num_threads,
		pos_mem_size_children;
	GRAPH* graph;
	int* mem_th_size_children;
	offspring* children_memory;
	
	graph_size = mat->n;
	perm_size  = perm_offset = total_size_children = 0;
	
	#pragma omp parallel
	{
		int node, n_par, n_ch, degree, size_children, size_ch_per_parent, th_id, th_size_children, 
			th_degree, th_mem_children;
		int* neighbors;
		GRAPH* children;
		
		#pragma omp single nowait
		num_threads = omp_get_num_threads();
		
		#pragma omp single nowait
		*perm  = calloc(graph_size, sizeof(int));
		
		#pragma omp single 
		graph  = calloc(graph_size, sizeof(GRAPH));
		
		#pragma omp single 
		mem_th_size_children = calloc(num_threads, sizeof(int));
		
		#pragma omp single 
		children_memory = calloc(num_threads, sizeof(offspring));
		
		#pragma omp for
		for (node = 0; node < graph_size; ++node)
		{
			graph[node].label    = node;
			graph[node].distance = INFINITY_LEVEL;
			graph[node].parent   = ORPHAN_NODE;
			graph[node].status   = UNREACHED;
			graph[node].degree   = GRAPH_degree(mat, node);
		}
		
		#pragma omp single nowait
		graph[root].distance = 0;
		
		#pragma omp single nowait
		graph[root].status = LABELED;
		
		#pragma omp single nowait
		(*perm)[0] = root;
		
		#pragma omp single nowait
		perm_size++;
		
		#pragma omp barrier
	
		while (perm_size < graph_size)
		{
			children         = NULL;
			size_children    = 0;
			th_degree        = 0;
			th_size_children = 0;
			th_mem_children  = 0;
			
			#pragma omp single nowait
			chunk_size = ceil((float) (perm_size - perm_offset) / num_threads);
			
// 			#pragma omp single 
// 			printf("[offset = %d, size = %d, chunksize = %d]\n", perm_offset, perm_size, chunk_size);fflush(stdout);
			
// 			#pragma omp single 
// 			printf("------------------------------------------------------- processing\n");fflush(stdout);

			#pragma omp single 
			pos_mem_size_children = 0;
			
			#pragma omp for schedule(static, chunk_size)
			for (n_par = perm_offset; n_par < perm_size; ++n_par)
			{
				#pragma omp critical
				{
					th_mem_children = pos_mem_size_children++;
					children_memory[th_mem_children].num = 0;
				}
				
				if (th_degree == 0)
				{
					for (node = n_par; node < n_par + chunk_size; ++node)
					{
// 						printf("Adding th_degree from node %d with degree %d\n", (*perm)[node], graph[(*perm)[node]].degree);
						th_degree += graph[(*perm)[node]].degree;
					}
					
// 					printf("th_degree: %d\n", th_degree);fflush(stdout);
					children = calloc(th_degree, sizeof(GRAPH));
				}
				
				size_ch_per_parent = 0;
// 				printf("Thread %d Processing n_par: %d from [%d, %d] chunk_size %d\n", omp_get_thread_num(), n_par, perm_offset, perm_size, chunk_size);fflush(stdout);
				degree    = graph[(*perm)[n_par]].degree;
				neighbors = GRAPH_adjacent(mat, (*perm)[n_par]);
				
// 				printf("Thread %d Processing node %d with degree %d and children: ", omp_get_thread_num(), (*perm)[n_par], graph[(*perm)[n_par]].degree);fflush(stdout);
// 				for (i = 0; i< graph[(*perm)[n_par]].degree; i++) printf("%d(%d) ", neighbors[i], graph[neighbors[i]].degree); fflush(stdout);
// 				printf("\n");fflush(stdout);
				
				// Processing children from a parent
				for (n_ch = 0; n_ch < degree; ++n_ch)
				{
					if (graph[neighbors[n_ch]].parent == ORPHAN_NODE)
						graph[neighbors[n_ch]].parent = (*perm)[n_par];
					
					if (graph[neighbors[n_ch]].distance > graph[(*perm)[n_par]].distance)
					{
						if (graph[neighbors[n_ch]].distance > graph[(*perm)[n_par]].distance + 1)
						{
							#pragma omp critical
							{
								graph[neighbors[n_ch]].distance = graph[(*perm)[n_par]].distance + 1;
									
								if (graph[neighbors[n_ch]].status == UNREACHED)
								{
									children[size_children++] = graph[neighbors[n_ch]]; 
									graph[neighbors[n_ch]].status = LABELED;
									size_ch_per_parent++;
// 									printf("Adding child %d at position %d\n", graph[neighbors[n_ch]].label, size_children-1);fflush(stdout);
// 									printf("Child %d added at position %d\n", children[size_children-1].label, size_children-1);
								}
							}
						}
						
						if (graph[(*perm)[n_par]].distance < graph[graph[neighbors[n_ch]].parent].distance)
						{
							#pragma omp critical
							graph[neighbors[n_ch]].parent = (*perm)[n_par];
						}
					}
				}
				
				if (size_ch_per_parent > 1) 
				{
// 					int i;
// 					printf("Before sorting\n");fflush(stdout);
// 					for (i = 0; i< size_ch_per_parent; i++) printf("%d(%d) ", children[i].label, children[i].degree); fflush(stdout);
// 					printf("\n");fflush(stdout);
					
					qsort(&(children[th_size_children]), size_ch_per_parent, sizeof(GRAPH), COMPARE_degr_ASC);
					
// 					printf("After sorting\n");fflush(stdout);
// 					for (i = 0; i< size_ch_per_parent; i++) printf("%d ", children[i].label); fflush(stdout);
// 					printf("\n");fflush(stdout);
				}
				
				th_size_children += size_ch_per_parent;
				
// 				printf("Children: ");fflush(stdout);
// 				for (i = 0; i< size_children; i++) printf("%d ", children[i].label); fflush(stdout);
// 				printf("\n");fflush(stdout);
				
				free(neighbors);
			}
			
			children_memory[th_mem_children].children = children;
			children_memory[th_mem_children].num = size_children;
			
			#pragma omp single nowait
			total_size_children = 0;
			
// 			#pragma omp single 
// 			printf("------------------------------------------------------- alloc\n");fflush(stdout);
			
// 			int th_ini, th_fim, pos;
			
			
			#pragma omp single 
			pos_mem_size_children = 0;
			
// 			#pragma omp single 
// 			printf("Before LOADING ................\n");fflush(stdout);
			
			#pragma omp for schedule(static, 1) 
			for (th_id = 0; th_id < num_threads; ++th_id)
			{
// 				#pragma omp critical
// 				{
// 				printf("Thread %d loading mem\n", omp_get_thread_num());fflush(stdout);
// 				}

// 				th_ini = perm_size;
// 				
// 				for (pos = 0; pos < th_id; ++pos)
// 					th_ini += mem_th_size_children[pos];
// 				
// 				for (node = 0; node < size_children; ++node)
// 					(*perm)[th_ini++] = children[node].label;


				
				#pragma omp critical
				{
					th_mem_children = pos_mem_size_children++;
					th_size_children = perm_size;
					perm_size += size_children;
					total_size_children += size_children;
				}
				
				for (node = 0; node < size_children; ++node)
					(*perm)[th_size_children++] = children[node].label;
				


				
				
				
				
// 				#pragma omp critical
// 				{
// // 					printf("Thread %d setting perm vector\n", omp_get_thread_num());fflush(stdout);
// 					th_size_children = perm_size;
// 					perm_size += size_children;
// 					total_size_children += size_children;
// 				}
// 				
// 				for (node = 0; node < size_children; ++node)
// 					(*perm)[th_size_children++] = children[node].label;
			}
			
					
// 			printf("Permutation vector: ");fflush(stdout);
// 			for (i = 0; i < perm_size; ++i)
// 				printf("%d ", (*perm)[i]);
// 			printf("\n");fflush(stdout);
			
// 			#pragma omp critical
// 			{
// 				perm_size += size_children;
// 				total_size_children += size_children;
// 			}
			
// 			#pragma omp single 
// 			printf("AFTER LOADING ................\n");fflush(stdout);
			
			#pragma omp single nowait
			free(children);
			
			#pragma omp single
			perm_offset = perm_size - total_size_children;
			
// 			#pragma omp single 
// 			printf("------------------------------------------------------- finish\n");fflush(stdout);
		}
	}
		
	free(graph);
		
// 	printf("Permutation vector of size %d: ", perm_size);fflush(stdout);
// 	int i;
// 	for (i = 0; i < perm_size; ++i)
// 		printf("%d ", (*perm)[i]);
// 	printf("\n");fflush(stdout);
}



