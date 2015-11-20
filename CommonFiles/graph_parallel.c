/*----------------------------------------------------------------------------
 * PARALLEL GRAPH FUNCTIONS
 *--------------------------------------------------------------------------*/
#include "protos.h"
#include <time.h>

/*----------------------------------------------------------------------------
 * Perform an unordered Breadth-First Search for the RCM reordering
 * 
 * param adjacency: adjacency matrix of graph
 * param root: start node
 * param levels: vector of levels of nodes
 *--------------------------------------------------------------------------*/
int* GRAPH_unordered_bfs(MAT* adjacency, int root, int* levels)
{
  int node, n_nodes, count_visited_nodes;
  
  n_nodes = adjacency->n;
  levels = calloc(n_nodes, sizeof(int));
  
  for (node = 0; node < n_nodes; ++node) levels[node] = INFINITY;
  levels[root] = 0;
  
  srand(time(NULL));
  printf("Number of nodes: %d, some node: %d\n", n_nodes, rand() % n_nodes);
  
  LIST* work_set = NULL;
  LIST_insert_IF_NOT_EXIST(work_set, root);
  
  while (count_visited_nodes != n_nodes) 
  {
    
  }
  
  return 0;
}