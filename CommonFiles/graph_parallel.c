/*----------------------------------------------------------------------------
 * PARALLEL GRAPH FUNCTIONS
 *--------------------------------------------------------------------------*/
#include "protos.h"
#include <time.h>
#include <omp.h>

/*----------------------------------------------------------------------------
 * Perform an unordered Breadth-First Search for the RCM reordering
 * 
 * param adjacency: adjacency matrix of graph
 * param root: start node
 * param levels: vector of levels of nodes
 *--------------------------------------------------------------------------*/
int* GRAPH_fixedpoint_bfs(MAT* adjacency, int root, int* levels)
{
  int node, n_nodes, adj_node, node_degree, count_nodes, count_visited_nodes, level;
  int* neighboors;
  
  n_nodes = adjacency->n;
  levels = calloc(n_nodes, sizeof(int));
  
  for (node = 0; node < n_nodes; ++node) levels[node] = INFINITY;
  levels[root] = 0;
  
//   printf("Number of nodes: %d, some node: %d\n", n_nodes, rand() % n_nodes);
//   srand(time(NULL));
  
  LIST* work_set = NULL;
  work_set = LIST_insert_IF_NOT_EXIST(work_set, root);
  
  while (count_visited_nodes != n_nodes) 
  {
	  node = LIST_first(work_set);
// 	  TODO: work_set = LIST_remove_first(work_set);
	  work_set = LIST_remove(work_set, node);
	  neighboors = GRAPH_adjacent(adjacency, node);
	  node_degree = GRAPH_degree(adjacency, node);
    
	  for (count_nodes = 0; count_nodes < node_degree; ++count_nodes)
	  {
		adj_node = neighboors[count_nodes];
		level = levels[node] + 1;
		if (level < levels[adj_node])
		{
			levels[adj_node] = level;
			work_set = LIST_insert_IF_NOT_EXIST(work_set, adj_node);
		}
		
	  }
	  
	  ++count_visited_nodes;
  }
  
  return levels;
}



int* GRAPH_parallel_fixedpoint_bfs(MAT* adjacency, int root, int* levels)
{
  int node, n_nodes, adj_node, node_degree, count_nodes, count_visited_nodes, level;
  int* neighboors;
  
  n_nodes = adjacency->n;
  levels = calloc(n_nodes, sizeof(int));
  
  for (node = 0; node < n_nodes; ++node) levels[node] = INFINITY;
  levels[root] = 0;
  
  LIST* work_set = NULL;
  work_set = LIST_insert_IF_NOT_EXIST(work_set, root);
  
  omp_set_num_threads(4);
  
  #pragma omp parallel private(node, neighboors, node_degree, count_nodes, adj_node, level)	
  {
	while (count_visited_nodes < n_nodes) 
	{
		node = -1;
		
		#pragma omp critical
		{
			if (work_set != NULL)
			{
				node = LIST_first(work_set);
				work_set = LIST_remove(work_set, node);
				++count_visited_nodes;
// 				printf("Thread: %d got node %d -> count_visited_nodes: ** %d **\n", omp_get_thread_num(), node, count_visited_nodes);
			}
		}
		
		if (node != -1) 
		{
// 			printf("Thread: %d processing node %d\n", omp_get_thread_num(), node);
			neighboors = GRAPH_adjacent(adjacency, node);
			node_degree = GRAPH_degree(adjacency, node);
		
			for (count_nodes = 0; count_nodes < node_degree; ++count_nodes)
			{
				adj_node = neighboors[count_nodes];
				
				#pragma omp critical
				{
					level = levels[node] + 1;
					if (level < levels[adj_node])
					{
// 						printf("Thread: %d setting level %d for node %d\n", omp_get_thread_num(), level, adj_node);
						levels[adj_node] = level;
						work_set = LIST_insert_IF_NOT_EXIST(work_set, adj_node);
					}
				}
			}
		}
	}
  }
  
  return levels;
}
