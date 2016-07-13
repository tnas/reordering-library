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


void GRAPH_parallel_fixedpoint_bfs(MAT* adjacency, int root, int** levels, const float percent_chunk)
{
	int node, n_nodes, has_unreached_nodes, size_workset;
	LIST* work_set;
	omp_lock_t lock;
	
	n_nodes = adjacency->n;
  
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
// 				while (work_set != NULL && count_chunk < BFS_WORK_CHUNK)
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
				
				neighboors  = GRAPH_adjacent(adjacency, active_node);
				node_degree = GRAPH_degree(adjacency, active_node);
				
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


void GRAPH_parallel_fixedpoint_BFS(METAGRAPH mgraph, int root, int** levels, const float percent_chunk)
{
	int node, n_nodes, has_unreached_nodes, size_workset;
	LIST* work_set;
	omp_lock_t lock;
	
	n_nodes = mgraph.size;
  
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
				
				neighboors  = mgraph.graph[active_node].neighboors;
				node_degree = mgraph.graph[active_node].degree;
				
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


/**
 * Build a METAGRAPH from an adjacency matrix using multithreads. 
 * 
 * @since 12-07-2016
 */
METAGRAPH* GRAPH_parallel_build(MAT* mat)
{
	METAGRAPH* meta_graph;
	int size_graph, min_degree;
	
	size_graph = mat->n;
	min_degree = INFINITY_LEVEL;
	
	meta_graph        = malloc(sizeof(METAGRAPH));
	meta_graph->size  = size_graph;
	meta_graph->graph = malloc(size_graph * sizeof(GRAPH));
	
	#pragma omp parallel 
	{
		int node, local_node_min_degree, local_min_degree;
	
		local_node_min_degree = 0;
		local_min_degree      = INFINITY_LEVEL;
		
		#pragma omp for schedule(static) 
		for (node = 0; node < size_graph; ++node)
		{
			meta_graph->graph[node].distance   = INFINITY_LEVEL;
			meta_graph->graph[node].parent     = -1;
			meta_graph->graph[node].chnum      = 0;
			meta_graph->graph[node].label      = node;
			meta_graph->graph[node].degree     = GRAPH_degree(mat, node);
			meta_graph->graph[node].neighboors = GRAPH_adjacent(mat, node);
			
			if (meta_graph->graph[node].degree < local_min_degree)
				local_node_min_degree = node;
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


BFS* GRAPH_parallel_execute_BFS(METAGRAPH mgraph, int root)
{
	int n_nodes, max_level, level, node;
	BFS* bfs;
	int* levels;
	int* counts;
	
	GRAPH_parallel_fixedpoint_BFS(mgraph, root, &levels, BFS_PERCENT_CHUNK);
	
	max_level   = count_nodes_by_level(levels, n_nodes, &counts);
	bfs->height = --max_level;
	
	#pragma omp parallel
	{
		#pragma omp single
		{
			bfs = malloc(sizeof(BFS));
			bfs->vertices_at_level = malloc(max_level * sizeof(int*));
		}
		
		#pragma omp for schedule(static) private(level)
		for (level = 0; level < max_level; ++level)
			bfs->vertices_at_level[level] = calloc(counts[level], sizeof(int));
		
		#pragma omp for schedule(static) private(node, level) 
		for (node = 0; node < n_nodes; ++node)
		{
			level = levels[node];
			
			#pragma omp critical
			{
				bfs->vertices_at_level[level][counts[level]] = node;
				--counts[level];
			}
		}
		
		#pragma omp single nowait
		
		
		#pragma omp single nowait
		free(levels);
		
		#pragma omp single nowait
		free(counts);
	}
	
	return bfs;
}


void GRAPH_parallel_shrinking_strategy_half_sorting()
{
}

void GRAPH_parallel_shrinking_strategy_vertex_by_degree()
{
}

void GRAPH_parallel_shrinking_strategy_five_non_adjacent()
{
}

graph_diameter* GRAPH_parallel_pseudodiameter(const METAGRAPH meta_graph)
{
	BFS* forwardBFS;
	BFS* reverseBFS;
	graph_diameter* diameter;
	int local_diameter;
	int* candidate_set;
	
	diameter = malloc(sizeof(graph_diameter));
	
	// Create two breadth first search engines
	
	// Initialize start and end vertices of pseudo-diameter
	diameter->start = meta_graph.vertex_min_degree;
	diameter->end   = NON_VERTEX;
	
	do 
	{
		// do BFS starting at start node 'start'
		forwardBFS = GRAPH_parallel_execute_BFS(meta_graph, diameter->start);
		
		// get candidate set of end nodes
		local_diameter = forwardBFS->height;
		candidate_set  = forwardBFS->vertices_at_level[local_diameter];
		
		// shrink candidate set to a manageable number
		
		
	} while (diameter->end == NON_VERTEX);
	
	return diameter;
}

