/*----------------------------------------------------------------------------
 * PARALLEL GRAPH FUNCTIONS
 *--------------------------------------------------------------------------*/
#include "protos_parallel.h"


// int inline has_work(const int* status_th)
// {
// 	int working_ok, count;
// 	for (count = 0, working_ok = 0; count < NUM_THREADS; working_ok |= status_th[count++]);
// 	return working_ok;
// }


// void inline update_work_set(const OPERATION op, LIST** ws, int* node)
// {
// 	#pragma omp critical (workset)
// 	{
// 		
// 		switch (op)
// 		{
// 			case READ :
// 				if ((*ws) != NULL)
// 				{
// 					*node = LIST_first(*ws);
// 					(*ws) = LIST_remove(*ws, *node);
// 				}
// 				break;
// 			case WRITE :
// 				*ws = LIST_insert_IF_NOT_EXIST(*ws, *node);
// // 				*ws = LIST_insert(*ws, *node);
// 				break;
// 		}
// 	}
// }



/*----------------------------------------------------------------------------
 * Perform an unordered Breadth-First Search for the RCM reordering
 * 
 * param adjacency: adjacency matrix of graph
 * param root: start node
 * param levels: vector of levels of nodes
 --------------------------------------------------------------------------*/
// void GRAPH_parallel_fixedpoint_bfs(MAT* adjacency, int root, int** levels)
// {
//   int node, n_nodes, adj_node, node_degree, count_nodes, level, count;
//   int* neighboors;
//   int status_threads[NUM_THREADS];
//   
//   n_nodes = adjacency->n;
//   *levels = calloc(n_nodes, sizeof(int));
//   
//   for (node = 0; node < n_nodes; ++node) (*levels)[node] = INFINITY;
//   (*levels)[root] = 0;
//   
//   LIST* work_set = NULL;
//   work_set = LIST_insert_IF_NOT_EXIST(work_set, root);
//   
//   omp_set_num_threads(NUM_THREADS);
//   for (count = 0; count < NUM_THREADS; status_threads[count++] = THREAD_ON);
//   
//   #pragma omp parallel private(node, neighboors, node_degree, count_nodes, adj_node, level)	
//   {
// 	while (has_work(status_threads)) 
// 	{
// 		node = -1;
// 		update_work_set(READ, &work_set, &node);
// 		
// 		if (node != -1)
// 		{
// 			status_threads[omp_get_thread_num()] = 1;
// 			
// 			neighboors  = GRAPH_adjacent(adjacency, node);
// 			node_degree = GRAPH_degree(adjacency, node);
// 		
// 			for (count_nodes = 0; count_nodes < node_degree; ++count_nodes)
// 			{
// 				adj_node = neighboors[count_nodes];
// 				
// 				#pragma omp critical
// 				{
// 					#pragma omp flush (levels)
// 					
// 					level = (*levels)[node] + 1;
// 					if (level < (*levels)[adj_node])
// 					{
// 						(*levels)[adj_node] = level;
// 						update_work_set(WRITE, &work_set, &adj_node);
// 					}
// 				}
// 			}
// 			
// 			free(neighboors);
// 		}
// 		else 
// 		{
// 			status_threads[omp_get_thread_num()] = THREAD_OFF;
// 		}
// 	}
//   }
//   
// //   if (work_set != NULL) 
// //   {
// // 	printf ("Error: Work set has nodes not processed. Exiting.. [GRAPH_parallel_fixedpoint]\n");
// // 	exit(0);
// //   }
// }


// void GRAPH_parallel_fixedpoint_bfs(MAT* adjacency, int root, int** levels)
// {
// 	int node, n_nodes, adj_node, node_degree, count_nodes, level, has_work, tid;
// 	int* neighboors;
// 	int status_threads[NUM_THREADS];
// 	LIST* work_set;
// 	omp_lock_t lock;
//   
// 	n_nodes = adjacency->n;
// 	*levels = calloc(n_nodes, sizeof(int));
//   
// 	omp_set_num_threads(NUM_THREADS);
// 	omp_init_lock(&lock);
// 	
// 	#pragma omp parallel for private (node)
// 	for (node = 0; node < n_nodes; ++node) 
// 		(*levels)[node] = INFINITY;
// 	(*levels)[root] = 0;
//   
// 	#pragma omp parallel private(node, neighboors, node_degree, count_nodes, adj_node, level, tid)	
// 	{
// 		tid = omp_get_thread_num();
// 		status_threads[tid] = THREAD_ON;
// 		
// 		#pragma omp single nowait
// 		{
// 			work_set = NULL;
// 			work_set = LIST_insert_IF_NOT_EXIST(work_set, root);
// 		}
// 		
// 		#pragma omp single
// 		has_work = NUM_THREADS;
// 		
// 		while (has_work) 
// 		{
// 			node = UNDEF_NODE;
// 			
// 			omp_set_lock(&lock);
// 			
// 			if (work_set != NULL)
// 			{
// 				node     = LIST_first(work_set);
// 				work_set = LIST_remove(work_set, node);
// 			}
// 			
// 			omp_unset_lock(&lock);
// 			
// 			if (node != UNDEF_NODE)
// 			{
// 				
// 				if (status_threads[tid] != THREAD_ON)
// 				{
// 					status_threads[tid] = THREAD_ON;
// 					
// 					#pragma omp atomic
// 					++has_work;
// 				}	
// 				
// 				neighboors  = GRAPH_adjacent(adjacency, node);
// 				node_degree = GRAPH_degree(adjacency, node);
// 				
// 				omp_set_lock(&lock);
// 				
// 				for (count_nodes = 0; count_nodes < node_degree; ++count_nodes)
// 				{
// 					adj_node = neighboors[count_nodes];
// 					level    = (*levels)[node] + 1;
// 					
// 					if (level < (*levels)[adj_node])
// 					{
// 						(*levels)[adj_node] = level;
// 						work_set = LIST_insert_IF_NOT_EXIST(work_set, adj_node);
// 					}
// 				}
// 				
// 				omp_unset_lock(&lock);
// 				
// 				free(neighboors);
// 			}
// 			else 
// 			{
// 				if (status_threads[tid] != THREAD_OFF)
// 				{
// 					status_threads[tid] = THREAD_OFF;
// 					
// 					#pragma omp atomic
// 					--has_work;
// 				}	
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
// }




// VERSÃO TOP
// void GRAPH_parallel_fixedpoint_bfs(MAT* adjacency, int root, int** levels)
// {
// 	int node, n_nodes, adj_node, node_degree, count_nodes, level, has_unreached_nodes;
// 	int* neighboors;
// 	LIST* work_set;
// 	omp_lock_t lock;
//   
// 	n_nodes = adjacency->n;
// 	*levels = calloc(n_nodes, sizeof(int));
//   
// 	omp_set_num_threads(NUM_THREADS);
// 	omp_init_lock(&lock);
// 	
// 	#pragma omp parallel for private (node)
// 	for (node = 0; node < n_nodes; ++node) 
// 		(*levels)[node] = INFINITY_LEVEL;
// 	(*levels)[root] = 0;
// 	
// 	work_set = NULL;
// 	work_set = LIST_insert_IF_NOT_EXIST(work_set, root);
// 	has_unreached_nodes = n_nodes - 1;
//   
// 	#pragma omp parallel private(node, neighboors, node_degree, count_nodes)
// 	{
// 		while (work_set != NULL || has_unreached_nodes > 0) 
// 		{
// 			node = UNDEF_NODE;
// 			
// 			omp_set_lock(&lock);
// 			
// 			if (work_set != NULL)
// 			{
// 				node     = LIST_first(work_set);
// 				work_set = LIST_remove(work_set, node);
// 			}
// 			
// 			omp_unset_lock(&lock);
// 			
// 			if (node != UNDEF_NODE)
// 			{
// 				
// 				neighboors  = GRAPH_adjacent(adjacency, node);
// 				node_degree = GRAPH_degree(adjacency, node);
// 				
// 				omp_set_lock(&lock);
// 				
// 				for (count_nodes = 0; count_nodes < node_degree; ++count_nodes)
// 				{
// 					adj_node = neighboors[count_nodes];
// 					level    = (*levels)[node] + 1;
// 					
// 					if (level < (*levels)[adj_node])
// 					{
// 						if ((*levels)[adj_node] == INFINITY_LEVEL)
// 							--has_unreached_nodes;
// 						
// 						(*levels)[adj_node] = level;
// 						work_set = LIST_insert_IF_NOT_EXIST(work_set, adj_node);
// 					}
// 				}
// 				
// 				omp_unset_lock(&lock);
// 				
// 				free(neighboors);
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
// }


// Sequencial mais rápido de todos
// void GRAPH_parallel_fixedpoint_bfs(MAT* adjacency, int root, int** levels)
// {
// 	int node, n_nodes, adj_node, node_degree, count_nodes, level;
// 	int* neighboors;
// 	LIST* work_set;
//   
// 	n_nodes = adjacency->n;
// 	*levels = calloc(n_nodes, sizeof(int));
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
// 	
// 	while (work_set != NULL) 
// 	{
// 		node     = LIST_first(work_set);
// 		work_set = LIST_remove(work_set, node);
// 		
// 		neighboors  = GRAPH_adjacent(adjacency, node);
// 		node_degree = GRAPH_degree(adjacency, node);
// 		
// 		for (count_nodes = 0; count_nodes < node_degree; ++count_nodes)
// 		{
// 			adj_node = neighboors[count_nodes];
// 			level    = (*levels)[node] + 1;
// 			
// 			if (level < (*levels)[adj_node])
// 			{
// 				
// 				(*levels)[adj_node] = level;
// 				work_set = LIST_insert_IF_NOT_EXIST(work_set, adj_node);
// 			}
// 			
// 		}
// 		
// 		
// 		free(neighboors);
// 	}
// 	
// 	
// 	//   if (work_set != NULL) 
// 	//   {
// 	// 	printf ("Error: Work set has nodes not processed. Exiting.. [GRAPH_parallel_fixedpoint]\n");
// 	// 	exit(0);
// 	//   }
// }



void GRAPH_parallel_fixedpoint_bfs(MAT* adjacency, int root, int** levels)
{
	int node, n_nodes, has_unreached_nodes, work_chunking;
	LIST* work_set;
	omp_lock_t lock;
	
	
	work_chunking = 128; // big matrices
// 	work_chunking = 64;  // little matrices
	n_nodes = adjacency->n;
  
	omp_set_num_threads(NUM_THREADS);
	
	#pragma omp parallel for private (node)
	for (node = 0; node < n_nodes; ++node) 
		(*levels)[node]   = INFINITY_LEVEL;
	(*levels)[root] = 0;
	
	work_set = NULL;
	work_set = LIST_insert_IF_NOT_EXIST(work_set, root);
	has_unreached_nodes = n_nodes - 1;
	
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
				
				while (work_set != NULL && count_chunk < work_chunking)
				{
					active_node     = LIST_first(work_set);
					work_set        = LIST_remove(work_set, active_node);
					active_chunk_ws = LIST_insert_IF_NOT_EXIST(active_chunk_ws, active_node);
					++count_chunk;
				}
				
				omp_unset_lock(&lock);
			}
			
			
// 			#pragma omp critical
// 			{
// 				while (cache_work_set != NULL)
// 				{
// 					active_node    = LIST_first(cache_work_set);
// 					cache_work_set = LIST_remove(cache_work_set, active_node);
// 					work_set       = LIST_insert_IF_NOT_EXIST(work_set, active_node);
// 				}
// 				
// 				count_chunk = 0;
// 				
// 				while (work_set != NULL && count_chunk < work_chunking)
// 				{
// 					active_node     = LIST_first(work_set);
// 					work_set        = LIST_remove(work_set, active_node);
// 					active_chunk_ws = LIST_insert_IF_NOT_EXIST(active_chunk_ws, active_node);
// 					++count_chunk;
// 				}
// 			}
			
			
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
 * Return the Last level shrinked from a LEVEL STRUCTURE
 *--------------------------------------------------------------------------*/
LIST* GRAPH_LS_last_level_PARALLEL(MAT* A, int* LS, int n)
{
	int i, k1, k2, last;
	GRAPH *G;
	LIST  *L;
	
	omp_set_num_threads(NUM_THREADS);
	
	#pragma omp parallel sections
	{
		#pragma omp section
		last = GRAPH_LS_depth(LS, n) - 1;
		
		#pragma omp section
		G = (GRAPH*) malloc (n*sizeof(GRAPH));
		
		#pragma omp section
		L = NULL;
	}
	
	#pragma omp parallel for private(i) 
	for (i = 0; i < n; ++i)
	{
		G[i].label    = i;
		G[i].distance = LS[i];
		G[i].degree   = GRAPH_degree (A,i);			
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
	int vertex, min_vertex, e, x, dgr, min_dgr, width;
	int n = width = A->n;
	int*  LS1;
	LIST* L;
	omp_lock_t lock_width, lock_L;
	
	omp_set_num_threads(NUM_THREADS);
	
	#pragma omp parallel sections
	{
		#pragma omp section
		LS1 =  (int*) calloc (n,sizeof (int));
		
		#pragma omp section
		{
			omp_init_lock(&lock_L);
			omp_init_lock(&lock_width);
		}
		
		#pragma omp section
		{
			L = NULL;
			min_dgr = n;
			min_vertex = INFINITY_LEVEL;
		}
	}
	
	/* Choose a vertex min_vertex with minimum degree */
	#pragma omp parallel for private(vertex, dgr) reduction(min: min_dgr)
	for (vertex = 0; vertex < n; ++vertex)
	{
		dgr = GRAPH_degree (A, vertex);
	
		if (min_dgr > dgr)
		{
			min_dgr    = dgr;
			min_vertex = vertex;
		}
	}
	
	/* Construct the Level Strucure of s */
	GRAPH_bfs (A, min_vertex, LS1);
// 	L = GRAPH_LS_last_level (A,LS1,n);
	L = GRAPH_LS_last_level_PARALLEL(A, LS1, n);
	
	
	
	
	#pragma omp parallel private(x)
	{
		int new_width;
		LIST* new_L;
		int*  LS2;
		
		new_width = width;
		new_L     = NULL;
		LS2 =  (int*) calloc (n,sizeof (int));
		
		while (L != NULL)
		{
			printf("Thread %d working\n", omp_get_thread_num()); fflush(stdout);
			
			omp_set_lock(&lock_L);
			
			if (L != NULL)
			{
				x   = LIST_first(L);
				L   = LIST_remove(L, x);
			}
			
			omp_unset_lock(&lock_L);
			
			GRAPH_bfs(A, x, LS2);
			
			new_width = GRAPH_LS_width(LS2, n);

			if (new_width < width)
			{
				if (GRAPH_LS_depth(LS2, n) > GRAPH_LS_depth(LS1, n))
				{
					#pragma omp critical
					min_vertex = x; 
					
					GRAPH_bfs(A, min_vertex, LS1);
					
					new_L = GRAPH_LS_last_level(A, LS1, n);
					
					omp_set_lock(&lock_L);
					L = new_L;
					omp_unset_lock(&lock_L);

					omp_set_lock(&lock_width);
					width = n;
					omp_unset_lock(&lock_width);
				}
				else
				{
					#pragma omp critical
					e = x; 
					
					omp_set_lock(&lock_width);
					width = new_width;
					omp_unset_lock(&lock_width);
				}
				
			}
			
// 			if (GRAPH_LS_depth (LS2,n) > GRAPH_LS_depth (LS1,n) && GRAPH_LS_width (LS2,n) < width)
// 			{
// 				min_vertex     = x; 
// 				GRAPH_bfs (A,min_vertex,LS1);
// 				L     = GRAPH_LS_last_level (A,LS1,n);
// 				width = n;
// 			}
// 			else if (GRAPH_LS_width (LS2,n) < width)
// 			{
// 				e     = x; 
// 				width = GRAPH_LS_width (LS2,n);
// 			}
		}
		
		free(LS2);
	}
	
	omp_destroy_lock(&lock_L);
	omp_destroy_lock(&lock_width);
	
	*node_s = min_vertex;
	*node_e = e;
	
	return LS1;	
}
