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
			
			#pragma omp for schedule(static) ordered
			for (n_par = perm_offset; n_par < perm_size; ++n_par)
			{
				degree    = graph[(*perm)[n_par]].degree;
				neighbors = GRAPH_adjacent(mat, (*perm)[n_par]);
				
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
			}
		}
		
		// *********************
		// Step 3: Prefix sum
		// *********************
		prefix_sum(counts, &psum, size_offset+1);
		
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
					}
				}
				
				if (child != -1)
				{
					#pragma omp critical
					{
						(*perm)[parent_index[graph[graph[child].parent].index]++] = child;
					}
					
					if (parent_index[graph[graph[child].parent].index] == psum[graph[graph[child].parent].index + 1])
					{
						num_children   = graph[graph[child].parent].chnum;
						index_children = psum[graph[graph[child].parent].index];
						
						if (num_children > 1)
						{
							// Sorting children by degree
							
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
		}
	}
	
	free(graph);
}



void Bucket_RCM(MAT* mat, int** perm, int root) 
{ 

	int graph_size, perm_size, perm_offset, chunk_size, total_size_children, num_threads,
		children_offset_gen, generation_size, gen_pos_reference;
	GRAPH* graph;
	genealogy* generation;
	
	graph_size = mat->n;
	perm_size  = perm_offset = total_size_children = 0;
	
	#pragma omp parallel
	{
		int node, n_par, n_ch, degree, pos_child_gen, pos_parent_gen, gen_pos, num_children;
		int* neighbors;
		
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
			pos_child_gen = 0;
			
			#pragma omp single
			generation_size = perm_size - perm_offset;
			
			#pragma omp single nowait
			chunk_size = ceil((float) generation_size / num_threads);
			
			#pragma omp single 
			generation = calloc(generation_size, sizeof(genealogy));
			
			#pragma omp for schedule(static, chunk_size) ordered
			for (n_par = perm_offset; n_par < perm_size; ++n_par)
			{
				degree    = graph[(*perm)[n_par]].degree;
				neighbors = GRAPH_adjacent(mat, (*perm)[n_par]);
				
				pos_parent_gen = n_par - perm_offset;
				generation[pos_parent_gen].num_children = 0;
				generation[pos_parent_gen].children = calloc(degree, sizeof(GRAPH));
				
				// Processing children from a parent
				for (n_ch = 0; n_ch < degree; ++n_ch)
				{
					if (graph[neighbors[n_ch]].parent == ORPHAN_NODE)
						graph[neighbors[n_ch]].parent = (*perm)[n_par];
					
					if (graph[neighbors[n_ch]].distance > graph[(*perm)[n_par]].distance)
					{
						if (graph[neighbors[n_ch]].distance > graph[(*perm)[n_par]].distance + 1)
						{
							#pragma omp ordered
							{
								graph[neighbors[n_ch]].distance = graph[(*perm)[n_par]].distance + 1;
									
								if (graph[neighbors[n_ch]].status == UNREACHED)
								{
									graph[neighbors[n_ch]].status = LABELED;
									generation[pos_parent_gen].children[generation[pos_parent_gen].num_children++] =  graph[neighbors[n_ch]];
								}
							}
						}
					}
				}
				
				free(neighbors);
			}
			
			#pragma omp single nowait
			total_size_children = 0;
			
			#pragma omp single nowait
			children_offset_gen = 0;
			
			#pragma omp single
			gen_pos_reference = 0;
			
			while (gen_pos_reference < generation_size)
			{
				gen_pos = generation_size;
				
				#pragma omp critical
				{
					if (gen_pos_reference < generation_size)
					{
						gen_pos       = gen_pos_reference++;
						num_children  = generation[gen_pos].num_children;
						pos_child_gen = gen_pos + perm_size + children_offset_gen;
						children_offset_gen += num_children - 1; // excluding itself node
						total_size_children += num_children;
					}
				}
				
				if (gen_pos < generation_size)
				{
					qsort(generation[gen_pos].children, num_children, sizeof(GRAPH), COMPARE_degr_ASC);
				
					for (node = 0; node < num_children; ++node)
						(*perm)[pos_child_gen + node] = generation[gen_pos].children[node].label;
				}
			}
				
			#pragma omp barrier
			
			#pragma omp for
			for (node = 0; node < generation_size; ++node)
				free(generation[node].children);
			
			#pragma omp single nowait
			free(generation);
			
			#pragma omp single
			perm_size += total_size_children;
			
			#pragma omp single
			perm_offset = perm_size - total_size_children;

		}
	}
		
	free(graph);
}
