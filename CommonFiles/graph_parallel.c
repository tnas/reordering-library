/*----------------------------------------------------------------------------
 * PARALLEL GRAPH FUNCTIONS
 *--------------------------------------------------------------------------*/
#include "protos_parallel.h"

// /*
// void GRAPH_parallel_fixedpoint_bfs(MAT* adjacency, int root, int** levels)
// {
// 	int node, n_nodes, has_unreached_nodes;
// 	LIST* work_set;
// 	omp_lock_t lock;
// 	
// 	n_nodes = adjacency->n;
//   
// 	omp_set_num_threads(NUM_THREADS);
// 	
// 	#pragma omp parallel for private (node)
// 	for (node = 0; node < n_nodes; ++node) 
// 		(*levels)[node]   = INFINITY_LEVEL;
// 	(*levels)[root] = 0;
// 	
// 	work_set = NULL;
// 	work_set = LIST_insert_IF_NOT_EXIST(work_set, root);
// 	has_unreached_nodes = n_nodes - 1;
// 	
// 	omp_init_lock(&lock);
// 	
// 	#pragma omp parallel 
// 	{
// 		int active_node, adj_node, node_degree, count_nodes, level, count_chunk;
// 		int* neighboors;
// 		LIST* cache_work_set;
// 		LIST* active_chunk_ws;
// 		
// 		cache_work_set = active_chunk_ws = NULL;
// 		
// 		while (work_set != NULL || has_unreached_nodes > 0) 
// 		{
// 
// 			if (work_set != NULL) 
// 			{
// 				count_chunk = 0;
// 				
// 				omp_set_lock(&lock);
// 				
// 				while (work_set != NULL && count_chunk < BFS_WORK_CHUNK)
// 				{
// 					active_node     = LIST_first(work_set);
// 					work_set        = LIST_remove(work_set, active_node);
// 					active_chunk_ws = LIST_insert_IF_NOT_EXIST(active_chunk_ws, active_node);
// 					++count_chunk;
// 				}
// 				
// 				omp_unset_lock(&lock);
// 			}
// 			
// 			while (active_chunk_ws != NULL)
// 			{
// 				active_node      = LIST_first(active_chunk_ws);
// 				active_chunk_ws  = LIST_remove(active_chunk_ws, active_node);
// 				
// 				neighboors  = GRAPH_adjacent(adjacency, active_node);
// 				node_degree = GRAPH_degree(adjacency, active_node);
// 				
// 				for (count_nodes = 0; count_nodes < node_degree; ++count_nodes)
// 				{
// 					adj_node = neighboors[count_nodes];
// 					level    = (*levels)[active_node] + 1;
// 					
// 					if (level < (*levels)[adj_node])
// 					{
// 						if ((*levels)[adj_node] == INFINITY_LEVEL) 
// 						{
// 							#pragma omp atomic
// 							--has_unreached_nodes;
// 						}
// 						
// 						#pragma omp critical
// 						(*levels)[adj_node] = level;
// 						
// 						cache_work_set = LIST_insert_IF_NOT_EXIST(cache_work_set, adj_node);
// 					}
// 					
// 				}
// 				
// 				
// 				free(neighboors);
// 			}
// 			
// 			
// 			if (cache_work_set != NULL)
// 			{
// 				omp_set_lock(&lock);
// 				
// 				while (cache_work_set != NULL)
// 				{
// 					active_node    = LIST_first(cache_work_set);
// 					cache_work_set = LIST_remove(cache_work_set, active_node);
// 					work_set       = LIST_insert_IF_NOT_EXIST(work_set, active_node);
// 				}
// 				
// 				omp_unset_lock(&lock);
// 			}
// 		}
// 	}
// 	
// 	omp_destroy_lock(&lock);
// 	
// 	//   if (work_set != NULL) 
// 	//   {
// 	// 	printf ("Error: Work set has nodes not processed. Exiting.. [GRAPH_parallel_fixedpoint]\n");
// 	// 	exit(0);
// 	//   }
// }*/


void GRAPH_parallel_fixedpoint_bfs(MAT* adjacency, int root, int** levels, const float percent_chunk)
{
	int node, n_nodes, has_unreached_nodes, size_workset;
	double size_chunk;
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
	
	//   if (work_set != NULL) 
	//   {
	// 	printf ("Error: Work set has nodes not processed. Exiting.. [GRAPH_parallel_fixedpoint]\n");
	// 	exit(0);
	//   }
}




/*----------------------------------------------------------------------------
 * Find the depth of an LEVEL STRUCTURE
 * (i.e. the total number of levels - same as find the maximum of an array)
 *--------------------------------------------------------------------------*/
int GRAPH_LS_depth_PARALLEL(int* LS, int n)
{
	if (LS == NULL || n <= 0)
	{
		printf ("error: Array is empty. Exiting.. [GRAPH_LS_depth]\n");
		exit(0);
	}

	int level, depth = 0;
	
	#pragma omp parallel for private(level) reduction(max: depth)
	for (level = 0; level < n; ++level)
		if (LS[level] > depth) depth = LS[level];

	return depth + 1;
}


/*----------------------------------------------------------------------------
 * Find the width of an LEVEL STRUCTURE
 * (i.e. maximum number of nodes which belong to a single level)
 *--------------------------------------------------------------------------*/
int GRAPH_LS_width_PARALLEL(int* LS, int n)
{
	int i, width, size_L;
	int* c;
	
	if (LS == NULL || n <= 0)
	{
		printf ("error: Array is empty. Exiting.. [GRAPH_LS_width]\n");
		exit(0);
	}
	
	width  = 0;	
	size_L = GRAPH_LS_depth_PARALLEL(LS,n);
	c      = calloc (size_L,sizeof(int));
	
	#pragma omp parallel
	{
		#pragma omp for private(i)
		for (i = 0; i < n; ++i)
			c[LS[i]]++;
		
		#pragma omp for private(i) reduction(max: width)
		for (i = 0; i < size_L; ++i)
			if (width < c[i])
				width = c[i];
			
		#pragma omp single
		free(c);
	}
	
	return width;
}



/*----------------------------------------------------------------------------
 * Return the Last level shrinked from a LEVEL STRUCTURE
 *--------------------------------------------------------------------------*/
LIST* GRAPH_LS_last_level_PARALLEL(MAT* A, int* LS, int n)
{
	int i, k1, k2, last;
	GRAPH *G;
	LIST  *L;
	
	#pragma omp parallel
	{
		#pragma omp sections
		{
			#pragma omp section
			last = GRAPH_LS_depth(LS, n) - 1;
			
			#pragma omp section
			G = (GRAPH*) malloc (n*sizeof(GRAPH));
			
			#pragma omp section
			L = NULL;
		}
		
		#pragma omp for private(i) 
		for (i = 0; i < n; ++i)
		{
			G[i].label    = i;
			G[i].distance = LS[i];
			G[i].degree   = GRAPH_degree (A,i);			
		}
	}
	
	qsort (G, n, sizeof(GRAPH), COMPARE_dist_degr_DES);

	k1 = k2 = 0;
	
	for (i = 0; i < n; ++i)
	{
		k1 = G[i].degree; 
		if (G[i].distance != last) 
			break;
		if (k1 != k2)
		{ 
			L  = LIST_insert_IF_NOT_EXIST (L, G[i].label);
			k2 = k1; 
		}
	}
	
	free(G);
	
	return L;
}



/*----------------------------------------------------------------------------
 * Find the pseudo-peripheral nodes in the graph (i.e. matrix) A 
 *--------------------------------------------------------------------------*/
int* GRAPH_LS_peripheral_PARALLEL (MAT* A, int *node_s, int* node_e)
{
	int n, dgr, vertex, min_vertex, e, x, min_dgr, width, new_width;
	int*  LS1;
	int*  LS2;
	LIST* level_list;
	
	n = width = A->n;
	
	#pragma omp parallel 
	{
		#pragma omp single nowait
		LS1 =  (int*) calloc (n, sizeof (int));
		
		#pragma omp single nowait
		LS2 =  (int*) calloc (n, sizeof (int));
		
		#pragma omp single nowait
		level_list = NULL;
		
		#pragma omp single
		{
			min_dgr    = n;
			min_vertex = INFINITY_LEVEL;
		}
		
		/* Choose a vertex min_vertex with minimum degree */
		#pragma omp for private(vertex, dgr) reduction(min: min_dgr, min_vertex)
		for (vertex = 0; vertex < n; ++vertex)
		{
			dgr = GRAPH_degree(A, vertex);
			
			if (min_dgr > dgr)
			{
				min_dgr    = dgr;
				min_vertex = vertex;
			}
		}
	}
	
	/* Construct the Level Strucure of s */
	GRAPH_bfs (A, min_vertex, LS1);
// 	GRAPH_parallel_fixedpoint_bfs(A, min_vertex, &LS1);
	level_list = GRAPH_LS_last_level_PARALLEL(A, LS1, n);
	
	while (level_list != NULL)
	{
		x          = LIST_first(level_list);
		level_list = LIST_remove(level_list, x);
		
		GRAPH_bfs(A, x, LS2);
// 		GRAPH_parallel_fixedpoint_bfs(A, x, &LS2);
		new_width = GRAPH_LS_width_PARALLEL(LS2, n);

		if (new_width < width)
		{
			if (GRAPH_LS_depth_PARALLEL(LS2, n) > GRAPH_LS_depth_PARALLEL(LS1, n))
			{
				min_vertex = x; 
				GRAPH_bfs(A, min_vertex, LS1);
// 				GRAPH_parallel_fixedpoint_bfs(A, min_vertex, &LS1);
				level_list = GRAPH_LS_last_level_PARALLEL(A, LS1, n);
				width = n;
			}
			else
			{
				e = x; 
				width = new_width;
			}
			
		}
	}
	
	free(LS2);
	
	*node_s = min_vertex;
	*node_e = e;
	
	return LS1;	
}
