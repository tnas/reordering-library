/*
 * Copyright 2016 Thiago Nascimento <nascimenthiago@gmail.com>
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

#include "graph_parallel.h"


/**
 * It executes a Breadth-First-Search using the fixed point concept. However, 
 * it uses a optimized data structure METAGRAPH, in order to increase the
 * performance of operations like take degree and neighboors of nodes.
 */
void GRAPH_parallel_fixedpoint_BFS(const METAGRAPH* mgraph, int root, int** levels, const float percent_chunk)
{
	int node, n_nodes, has_unreached_nodes, size_workset;
	LIST* work_set;
	omp_lock_t lock;
	
	n_nodes = mgraph->size;
  
	#pragma omp parallel for private (node)
	for (node = 0; node < n_nodes; ++node) 
		(*levels)[node]   = INFINITY_LEVEL;
	(*levels)[root] = 0;
	
	work_set = NULL;
	work_set = LIST_insert_IF_NOT_EXIST(work_set, root);
	has_unreached_nodes = n_nodes - 1;
	size_workset = 1;
	
	omp_init_lock(&lock);
	
	#pragma omp parallel 
	{
		int active_node, adj_node, node_degree, count_nodes, level, count_chunk;
		int* neighboors;
		LIST* cache_work_set;
		LIST* active_chunk_ws;
		double size_chunk;
		
		cache_work_set = active_chunk_ws = NULL;
		
		while (work_set != NULL || has_unreached_nodes > 0) 
		{

			if (work_set != NULL) 
			{
				count_chunk = 0;
				
				omp_set_lock(&lock);
				
				size_chunk = percent_chunk * size_workset;
				
				while (work_set != NULL && count_chunk < size_chunk)
				{
					active_node     = LIST_first(work_set);
					work_set        = LIST_remove(work_set, active_node);
					--size_workset;
					active_chunk_ws = LIST_insert_IF_NOT_EXIST(active_chunk_ws, active_node);
					++count_chunk;
				}
				
				omp_unset_lock(&lock);
			}
			
			while (active_chunk_ws != NULL)
			{
				active_node      = LIST_first(active_chunk_ws);
				active_chunk_ws  = LIST_remove(active_chunk_ws, active_node);
				
				neighboors  = GRAPH_adjacent(mgraph->mat, active_node);
				node_degree = mgraph->graph[active_node].degree;
				
				for (count_nodes = 0; count_nodes < node_degree; ++count_nodes)
				{
					adj_node = neighboors[count_nodes];
					level    = (*levels)[active_node] + 1;
					
					if (level < (*levels)[adj_node])
					{
						if ((*levels)[adj_node] == INFINITY_LEVEL) 
						{
							#pragma omp atomic
							--has_unreached_nodes;
						}
						
						#pragma omp critical
						(*levels)[adj_node] = level;
						
						cache_work_set = LIST_insert_IF_NOT_EXIST(cache_work_set, adj_node);
					}
				}
				
				free(neighboors);
			}
			
			if (cache_work_set != NULL)
			{
				omp_set_lock(&lock);
				
				while (cache_work_set != NULL)
				{
					active_node    = LIST_first(cache_work_set);
					cache_work_set = LIST_remove(cache_work_set, active_node);
					work_set       = LIST_insert_IF_NOT_EXIST(work_set, active_node);
					++size_workset;
				}
				
				omp_unset_lock(&lock);
			}
		}
	}
	
	omp_destroy_lock(&lock);
}




void GRAPH_parallel_fixedpoint_static_BFS(const METAGRAPH* mgraph, int root, int** levels, const float percent_chunk)
{
	int node, n_nodes, ws_size, has_unreached_nodes, tail, head;
	int* work_set;
	omp_lock_t lock;
	
	n_nodes = mgraph->size;
	ws_size = (omp_get_max_threads() + 1) * n_nodes; // oversizing estimate
	
	#pragma omp parallel 
	{
		int active_node, adj_node, node_degree, count_nodes, level, 
		    count_chunk, size_chunk;
		int* neighboors;
		int* cache_work_set;
		int* active_chunk_ws;
		int th_head, th_tail, active_head, active_tail, cache_head, cache_tail;
		
		#pragma omp for private(node)
		for (node = 0; node < n_nodes; ++node) 
			(*levels)[node] = INFINITY_LEVEL;
		
		#pragma omp sections 
		{
			#pragma omp section
			(*levels)[root] = 0;
			
			#pragma omp section
			{
				work_set = calloc(ws_size, sizeof(int)); 
				tail = head = 0;
				QUEUE_enque(&work_set, ws_size, &tail, root);
			}
			
			#pragma omp section
			has_unreached_nodes = n_nodes - 1;
			
			#pragma omp section
			omp_init_lock(&lock);
		}
		
		active_chunk_ws = calloc(n_nodes, sizeof(int));
		active_head = active_tail = 0;
		
		cache_work_set = calloc(n_nodes, sizeof(int));
		cache_head = cache_tail = 0;
		
		while ((!QUEUE_empty(work_set, head, tail)) || has_unreached_nodes > 0) 
		{
			if (!QUEUE_empty(work_set, head, tail)) 
			{
				#pragma omp critical
				{
					th_tail = tail;
					th_head = head;	
					size_chunk = ceil(percent_chunk * (th_tail - th_head));
					head += size_chunk;
				}
					
				for (count_chunk = 0; count_chunk < size_chunk; ++count_chunk)
				{
					active_node = QUEUE_deque(&work_set, ws_size, &th_head);
					QUEUE_enque(&active_chunk_ws, n_nodes, &active_tail, active_node);
				}
			}
			
			// Fixed Point Iteration
			while (!QUEUE_empty(active_chunk_ws, active_head, active_tail))
			{
				active_node = QUEUE_deque(&active_chunk_ws, n_nodes, &active_head);
				
				neighboors  = GRAPH_adjacent(mgraph->mat, active_node);
				node_degree = mgraph->graph[active_node].degree;
				level       = (*levels)[active_node] + 1;
				
				for (count_nodes = 0; count_nodes < node_degree; ++count_nodes)
				{
					adj_node = neighboors[count_nodes];
					
					if (level < (*levels)[adj_node])
					{
						if ((*levels)[adj_node] == INFINITY_LEVEL) 
						{
							#pragma omp atomic
							--has_unreached_nodes;
						}
						
						#pragma omp critical
						if (level < (*levels)[adj_node]) (*levels)[adj_node] = level;
						
						QUEUE_enque(&cache_work_set, n_nodes, &cache_tail, adj_node);
					}
				}
				
				free(neighboors);
			}
			
			if (!QUEUE_empty(cache_work_set, cache_head, cache_tail))
			{
				omp_set_lock(&lock);
				
				while (!QUEUE_empty(cache_work_set, cache_head, cache_tail))
				{
					active_node = QUEUE_deque(&cache_work_set, n_nodes, &cache_head);
					QUEUE_enque(&work_set, ws_size, &tail, active_node);
				}
				
				omp_unset_lock(&lock);
			}
		}
		
		free(active_chunk_ws);
		free(cache_work_set);
	}
	
	omp_destroy_lock(&lock);
	free(work_set);
}



void GRAPH_parallel_fixedpoint_sloan_BFS(METAGRAPH* mgraph, int root, const float percent_chunk)
{
	int node, n_nodes, ws_size, has_unreached_nodes, tail, head;
	int* work_set;
	int* levels;
	omp_lock_t lock;
	n_nodes = mgraph->size;
	ws_size = (omp_get_max_threads() + 1) * n_nodes; // oversizing estimate
	
	#pragma omp parallel 
	{
		int active_node, adj_node, node_degree, count_nodes, level, 
		    count_chunk, size_chunk;
		int* neighboors;
		int* cache_work_set;
		int* active_chunk_ws;
		int th_head, th_tail, active_head, active_tail, cache_head, cache_tail;
		
		#pragma omp sections
		{
			#pragma omp section
			levels = calloc(n_nodes, sizeof(int));
			
			#pragma omp section
			mgraph->sloan_priority = calloc(n_nodes, sizeof(int));
			
			#pragma omp section
			mgraph->min_sloan_priority = INFINITY_LEVEL;
			
			#pragma omp section
			mgraph->max_sloan_priority = -INFINITY_LEVEL;
		}
		
		#pragma omp for private(node)
		for (node = 0; node < n_nodes; ++node) 
			levels[node] = INFINITY_LEVEL;
		
		#pragma omp sections 
		{
			#pragma omp section
			levels[root] = 0;
			
			#pragma omp section
			{
				work_set = calloc(ws_size, sizeof(int)); 
				tail = head = 0;
				QUEUE_enque(&work_set, ws_size, &tail, root);
			}
			
			#pragma omp section
			has_unreached_nodes = n_nodes - 1;
			
			#pragma omp section
			omp_init_lock(&lock);
		}
		
		active_chunk_ws = calloc(n_nodes, sizeof(int));
		active_head = active_tail = 0;
		
		cache_work_set = calloc(n_nodes, sizeof(int));
		cache_head = cache_tail = 0;
		
		while ((!QUEUE_empty(work_set, head, tail)) || has_unreached_nodes > 0) 
		{
			if (!QUEUE_empty(work_set, head, tail)) 
			{
				#pragma omp critical
				{
					th_tail = tail;
					th_head = head;	
					size_chunk = ceil(percent_chunk * (th_tail - th_head));
					head += size_chunk;
				}
					
				for (count_chunk = 0; count_chunk < size_chunk; ++count_chunk)
				{
					active_node = QUEUE_deque(&work_set, ws_size, &th_head);
					QUEUE_enque(&active_chunk_ws, n_nodes, &active_tail, active_node);
				}
			}
			
			// Fixed Point Iteration
			while (!QUEUE_empty(active_chunk_ws, active_head, active_tail))
			{
				active_node = QUEUE_deque(&active_chunk_ws, n_nodes, &active_head);
				
				neighboors  = GRAPH_adjacent(mgraph->mat, active_node);
				node_degree = mgraph->graph[active_node].degree;
				
				for (count_nodes = 0; count_nodes < node_degree; ++count_nodes)
				{
					adj_node = neighboors[count_nodes];
					level    = levels[active_node] + 1;
					
					if (level < levels[adj_node])
					{
						if (levels[adj_node] == INFINITY_LEVEL) 
						{
							#pragma omp atomic
							--has_unreached_nodes;
						}
						
						#pragma omp critical
						if (level < levels[adj_node])
						{
							levels[adj_node] = level;
							
							mgraph->sloan_priority[adj_node] = 
								SLOAN_W1*levels[adj_node] - SLOAN_W2*(mgraph->graph[adj_node].degree + 1);
								
							if (mgraph->sloan_priority[adj_node] > mgraph->max_sloan_priority)
								mgraph->max_sloan_priority = mgraph->sloan_priority[adj_node];
							
							if (mgraph->sloan_priority[adj_node] < mgraph->min_sloan_priority)
								mgraph->min_sloan_priority = mgraph->sloan_priority[adj_node];
						}
						
						QUEUE_enque(&cache_work_set, n_nodes, &cache_tail, adj_node);
					}
				}
				
				free(neighboors);
			}
			
			if (!QUEUE_empty(cache_work_set, cache_head, cache_tail))
			{
				omp_set_lock(&lock);
				
				while (!QUEUE_empty(cache_work_set, cache_head, cache_tail))
				{
					active_node = QUEUE_deque(&cache_work_set, n_nodes, &cache_head);
					QUEUE_enque(&work_set, ws_size, &tail, active_node);
				}
				
				omp_unset_lock(&lock);
			}
		}
		
		free(active_chunk_ws);
		free(cache_work_set);
	}
	
	omp_destroy_lock(&lock);
	free(levels);
	free(work_set);
}




/**
 * Build a METAGRAPH from an adjacency matrix using multithreads. 
 * 
 * @since 12-07-2016
 */
inline METAGRAPH* GRAPH_parallel_build_METAGRAPH(MAT* mat)
{
	METAGRAPH* meta_graph;
	int size_graph, min_degree;
	
	size_graph = mat->n;
	min_degree = INFINITY_LEVEL;
	
	meta_graph                 = malloc(sizeof(METAGRAPH));
	meta_graph->size           = size_graph;
	meta_graph->graph          = malloc(size_graph * sizeof(GRAPH));
	meta_graph->mat            = mat;
	meta_graph->edges          = mat->nz;
	meta_graph->lock_node      = malloc(size_graph * sizeof(omp_lock_t));
	meta_graph->sloan_priority = NULL;
	
	#pragma omp parallel 
	{
		int node, local_node_min_degree, local_min_degree;
	
		local_node_min_degree = 0;
		local_min_degree      = INFINITY_LEVEL;
		
		#pragma omp for schedule(static) 
		for (node = 0; node < size_graph; ++node)
		{
			meta_graph->graph[node].distance   = INFINITY_LEVEL;
			meta_graph->graph[node].parent     = NON_VERTEX;
			meta_graph->graph[node].status     = UNREACHED;
			meta_graph->graph[node].chnum      = 0;
			meta_graph->graph[node].label      = node;
			meta_graph->graph[node].degree     = GRAPH_degree(mat, node);
			omp_init_lock(&meta_graph->lock_node[node]); 
			
			if (meta_graph->graph[node].degree < local_min_degree)
			{
				local_node_min_degree = node;
				local_min_degree      = meta_graph->graph[node].degree;
			}
		}
		
		if (local_min_degree < min_degree)
		{
			#pragma omp critical
			{
				if (local_min_degree < min_degree)
				{
					min_degree = local_min_degree;
					meta_graph->vertex_min_degree = local_node_min_degree;
				}
			}
		}
	}
	
	return meta_graph;
}


/**
 * Build a BFS data structure from a METAGRAPH parameter. For
 * this, the function execute a Breadht-First-Search throught
 * the GRAPH_parallel_fixedpoint_BFS procedure. A root node
 * must to be provided.
 * 
 */
inline BFS* GRAPH_parallel_build_BFS(const METAGRAPH* mgraph, int root)
{
	int n_nodes, max_level;
	BFS* bfs;
	int* levels;
	int* counts;
	
	n_nodes = mgraph->size;
	bfs     = malloc(sizeof(BFS));
	levels  = calloc(n_nodes, sizeof(int));
	
	GRAPH_parallel_fixedpoint_static_BFS(mgraph, root, &levels, BFS_PERCENT_CHUNK);
	
	max_level   = count_nodes_by_level(levels, n_nodes, &counts);
	// decrementing two levels added to max_level by count_nodes_by_level
	bfs->height = max_level - 2;
	bfs->width  = 0;
	
	#pragma omp parallel
	{
		GRAPH graph_node;
		int level, count_level, node;
		
		#pragma omp sections
		{
			#pragma omp section
			bfs->vertices_at_level  = calloc((bfs->height + 1), sizeof(GRAPH*));
			
			#pragma omp section
			bfs->num_nodes_at_level = calloc((bfs->height + 1), sizeof(int));
		}
		
		// Allocating memory for BFS structure
		#pragma omp for schedule(static)
		for (count_level = 1; count_level < max_level; ++count_level)
		{
			level = count_level - 1;
			
			bfs->vertices_at_level[level]  = calloc(counts[count_level], sizeof(GRAPH));
			bfs->num_nodes_at_level[level] = counts[count_level];
			
			#pragma omp critical
			if (counts[count_level] > bfs->width) 
				bfs->width = counts[count_level];
		}
		
		// Filling BFS structure with respective datas
		#pragma omp for schedule(static) 
		for (node = 0; node < n_nodes; ++node)
		{
			level       = levels[node];
			count_level = level + 1;
			
			#pragma omp critical
			{
				graph_node.label  = node;
				graph_node.degree = mgraph->graph[node].degree; 
				bfs->vertices_at_level[level][counts[count_level]-1] = graph_node;
				counts[count_level]--;
			}
		}

		#pragma omp sections
		{
			#pragma omp section
			free(levels);
			
			#pragma omp section
			free(counts);
		}
	}
	
	return bfs;
}


/**
 * Shrinking strategy in which a single node of each degree is chosen.
 * 
 * @since 22-07-2016
 */
inline GRAPH* GRAPH_shrinking_strategy_vertex_by_degree(GRAPH* nodes, int* length)
{
	int node, shrink_length, last_degree;
	GRAPH* shrinked_nodes;
	
	if (*length == 1)
	{
		shrinked_nodes = calloc(*length, sizeof(GRAPH));
		memcpy(shrinked_nodes, nodes, (*length) * sizeof(GRAPH));
		
		return shrinked_nodes;
	}
	
	qsort(nodes, *length, sizeof(GRAPH), COMPARE_degr_ASC);
	shrinked_nodes = calloc(*length, sizeof(GRAPH));
	
	shrink_length  = 0;
	last_degree = -1;
	
	for (node = 0; node < *length; ++node)
	{
		if (nodes[node].degree != last_degree)
		{
			last_degree = nodes[node].degree;
			shrinked_nodes[shrink_length++] = nodes[node];
		}
	}
	
	*length = shrink_length;
	
	return shrinked_nodes;
}


/**
 * Shrinking strategy in which the candidate set is sorted by node degree
 * and the first five nodes are chosen. Each selected vertex must not be adjacent 
 * to any previously chosen nodes.
 * 
 * @since 22-07-2016
 */
inline GRAPH* GRAPH_shrinking_strategy_five_non_adjacent(GRAPH* nodes, int* length)
{
	int node, num_chosen, length_shrink_nodes, chosen, neigh, chose_node;
	GRAPH* shrinked_nodes;
	
	length_shrink_nodes = 5; // According to Kumfert
	
	if (*length < length_shrink_nodes) 
	{
		shrinked_nodes = calloc(*length, sizeof(GRAPH));
		memcpy(shrinked_nodes, nodes, (*length) * sizeof(GRAPH));
		
		return shrinked_nodes;
	}
	
	num_chosen     = 0;
	shrinked_nodes = calloc(length_shrink_nodes, sizeof(GRAPH));
	
	qsort(nodes, *length, sizeof(GRAPH), COMPARE_degr_ASC);
	
	for (node = 0; node < *length; ++node)
	{
		if (num_chosen == length_shrink_nodes) break;
		
		chose_node = 1;
		
		for (chosen = 0; chosen < num_chosen && chose_node == 0; ++chosen)
		{
			for (neigh = 0; neigh < shrinked_nodes[chosen].degree; ++neigh)
			{
				if (shrinked_nodes[chosen].neighboors[neigh] == nodes[node].label) 
				{
					chose_node = 0;
					break;
				}
			}
		}
		
		if (chose_node) shrinked_nodes[num_chosen++] = nodes[node];
	}
	
	*length = num_chosen;
	
	return shrinked_nodes;
}


/**
 * Shrinking strategy in which nodes are sorted by degree and half
 * of them are chosen.
 * 
 * @since 14-07-2016
 */
inline GRAPH* GRAPH_shrinking_strategy_half_sorted(GRAPH* nodes, int* length)
{
	GRAPH* shrinked_nodes;
	
	if (*length == 1)
	{
		shrinked_nodes = calloc(*length, sizeof(GRAPH));
		memcpy(shrinked_nodes, nodes, (*length) * sizeof(GRAPH));
		
		return shrinked_nodes;
	}
		
	qsort(nodes, *length, sizeof(GRAPH), COMPARE_degr_ASC);
	*length /= 2;
	shrinked_nodes = calloc(*length, sizeof(GRAPH));
	memcpy(shrinked_nodes, nodes, (*length) * sizeof(GRAPH));
	
	return shrinked_nodes;
}

static GRAPH* (*strategy[3])(GRAPH* nodes, int* length) = {
	GRAPH_shrinking_strategy_half_sorted,
	GRAPH_shrinking_strategy_vertex_by_degree,
	GRAPH_shrinking_strategy_five_non_adjacent
};

/**
 * It calculates the pseudo diameter of a graph represented by the
 * METAGRAPH data structure. This algorithm was proposed by 
 * Gary karl Kumfert in his Phd Thesis - An Object-Oriented Algorithmic
 * Laboratory for Ordering Sparse Matrices (2000).
 * 
 */
graph_diameter* GRAPH_parallel_pseudodiameter(const METAGRAPH* meta_graph, Shrinking_Strategy shrink_type)
{
	// Create two breadth first search engines
	BFS* forwardBFS;
	BFS* reverseBFS;
	graph_diameter* diameter;
	int local_diameter, size_cand_set, min_width, cand, candidate, swap;
	GRAPH* candidate_set;
	
	diameter = malloc(sizeof(graph_diameter));
	forwardBFS = reverseBFS = NULL;
	
	// Initialize start and end vertices of pseudo-diameter
	diameter->start = meta_graph->vertex_min_degree;
	diameter->end   = NON_VERTEX;
	
	do 
	{
		if (forwardBFS != NULL) GRAPH_parallel_destroy_BFS(forwardBFS);
		
		// do BFS starting at start node 'start'
		forwardBFS = GRAPH_parallel_build_BFS(meta_graph, diameter->start);
		
		// get candidate set of end nodes
		local_diameter = forwardBFS->height;
		candidate_set  = forwardBFS->vertices_at_level[local_diameter];
		size_cand_set  = forwardBFS->num_nodes_at_level[local_diameter];
		
		// shrink candidate set to a manageable number
		if (shrink_type == FIVE_NON_ADJACENT)
		{
			for (cand = 0; cand < size_cand_set; ++cand)
				candidate_set[cand].neighboors = GRAPH_adjacent(meta_graph->mat, candidate_set[cand].label);
		}
		
		candidate_set = strategy[shrink_type](candidate_set, &size_cand_set);
		
		min_width = INT_MAX;
		
		for (cand = 0; cand < size_cand_set; ++cand)
		{
			candidate = candidate_set[cand].label;
			
			// do BFS from each candidate
			if (reverseBFS != NULL) GRAPH_parallel_destroy_BFS(reverseBFS);
			reverseBFS = GRAPH_parallel_build_BFS(meta_graph, candidate);
			
			if (reverseBFS->width < min_width)
			{
				if (reverseBFS->height > local_diameter)
				{
					// reverseBFS is better than the forwardBFS
					// reset algorithm with candidate as new start
					diameter->start = candidate;
					diameter->end   = NON_VERTEX;
					cand = size_cand_set; // break
				}
				else
				{
					// reverseBFS is narrower than any others
					// make this new end node
					min_width = reverseBFS->width;
					diameter->end = candidate;
				}
			}
		}
		
		if (shrink_type == FIVE_NON_ADJACENT)
		{
			for (cand = 0; cand < size_cand_set; ++cand)
				free(candidate_set[cand].neighboors);
		}
		
		free(candidate_set);
		
	} while (diameter->end == NON_VERTEX);
	
	// swap start & end if the reverseBFS is narrower than forwardBFS
	if (forwardBFS->width > reverseBFS->width)
	{
		swap = diameter->start;
		diameter->start = diameter->end;
		diameter->end   = swap;
	}
	
	GRAPH_parallel_destroy_BFS(forwardBFS);
	GRAPH_parallel_destroy_BFS(reverseBFS);
	
	return diameter;
}


/**
 * Clean the memory used by METAGRAPH structure.
 * 
 * @since 19-07-2016
 */
void inline GRAPH_parallel_destroy_METAGRAPH(METAGRAPH* mgraph)
{
	int node, size_graph;
	
	size_graph = mgraph->size;
	
	#pragma omp parallel for schedule(static) private(node)
	for (node = 0; node < size_graph; ++node) 
		omp_destroy_lock(&mgraph->lock_node[node]);
	
	if (mgraph->sloan_priority != NULL)
		free(mgraph->sloan_priority);
	
	MATRIX_clean(mgraph->mat);
	free(mgraph->graph);
	free(mgraph->lock_node);
	free(mgraph);
}


/**
 * Clean the memory used by BFS structure.
 * 
 * @since 19-07-2016
 */
void inline GRAPH_parallel_destroy_BFS(BFS* bfs)
{
	int level, max_level;
	
	max_level = bfs->height + 1;
	
	#pragma omp parallel for schedule(static) private(level)
	for (level = 0; level < max_level; ++level)
		free(bfs->vertices_at_level[level]);
	
	free(bfs->vertices_at_level);
	free(bfs->num_nodes_at_level);
	
	free(bfs);
}