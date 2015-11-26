/*----------------------------------------------------------------------------
 * PARALLEL GRAPH FUNCTIONS
 *--------------------------------------------------------------------------*/
#include "protos.h"
#include <time.h>
#include <omp.h>

#define NUM_THREADS 8
#define THREAD_ON 1
#define THREAD_OFF 0

typedef enum OPERATION {READ, WRITE} OPERATION;



int inline has_work(const int* status_th)
{
	int working_ok, count;
	for (count = 0, working_ok = 0; count < NUM_THREADS; working_ok |= status_th[count++]);
	return working_ok;
}


void inline update_work_set(const OPERATION op, LIST** ws, int* node)
{
	#pragma omp critical (workset)
	{
		#pragma omp flush (ws)
		
		switch (op)
		{
			case READ :
				if ((*ws) != NULL)
				{
					*node = LIST_first(*ws);
					(*ws) = LIST_remove(*ws, *node);
				}
				break;
			case WRITE :
				*ws = LIST_insert_IF_NOT_EXIST(*ws, *node);
				break;
		}
	}
}



/*----------------------------------------------------------------------------
 * Perform an unordered Breadth-First Search for the RCM reordering
 * 
 * param adjacency: adjacency matrix of graph
 * param root: start node
 * param levels: vector of levels of nodes
 --------------------------------------------------------------------------*/
int* GRAPH_parallel_fixedpoint_bfs(MAT* adjacency, int root, int* levels)
{
  int node, n_nodes, adj_node, node_degree, count_nodes, level, count;
  int* neighboors;
  int status_threads[NUM_THREADS];
  
  n_nodes = adjacency->n;
  levels = calloc(n_nodes, sizeof(int));
  
  for (node = 0; node < n_nodes; ++node) levels[node] = INFINITY;
  levels[root] = 0;
  
  LIST* work_set = NULL;
  work_set = LIST_insert_IF_NOT_EXIST(work_set, root);
  
  omp_set_num_threads(NUM_THREADS);
  for (count = 0; count < NUM_THREADS; status_threads[count++] = THREAD_ON);
  
  #pragma omp parallel private(node, neighboors, node_degree, count_nodes, adj_node, level)	
  {
	while (has_work(status_threads)) 
	{
		node = -1;
		update_work_set(READ, &work_set, &node);
		
		if (node != -1)
		{
			status_threads[omp_get_thread_num()] = 1;
			
			neighboors  = GRAPH_adjacent(adjacency, node);
			node_degree = GRAPH_degree(adjacency, node);
		
			for (count_nodes = 0; count_nodes < node_degree; ++count_nodes)
			{
				adj_node = neighboors[count_nodes];
				
				#pragma omp critical
				{
					#pragma omp flush (levels)
					
					level = levels[node] + 1;
					if (level < levels[adj_node])
					{
						levels[adj_node] = level;
						update_work_set(WRITE, &work_set, &adj_node);
					}
				}
			}
		}
		else 
		{
			status_threads[omp_get_thread_num()] = THREAD_OFF;
		}
	}
  }
  
  if (work_set != NULL) 
  {
	printf ("Error: Work set has nodes not processed. Exiting.. [GRAPH_parallel_fixedpoint]\n");
	exit(0);
  }
  
  return levels;
}

