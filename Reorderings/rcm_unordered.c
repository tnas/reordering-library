/*----------------------------------------------------------------------------
 * UNORDERED RCM REORDERING SOLVER
 *--------------------------------------------------------------------------*/
#include "../CommonFiles/protos.h"

/*----------------------------------------------------------------------------
 * Unordered RCM reordering from the LEVEL STRUCTURE 
 *--------------------------------------------------------------------------*/
void REORDERING_RCM_parallel(MAT* A, int** perm)
{ 
	int n_nodes, root, e, count_nodes;
	int* tperm;
	
	n_nodes = A->n;
	tperm = calloc (n_nodes, sizeof(int));
	
// 	int* g = GRAPH_LS_peripheral (A, &root, &e);
	
	root = 4;
	tperm = GRAPH_parallel_fixedpoint_bfs(A, root, tperm);
	
	/* Reverse order */
	for (count_nodes = 0; count_nodes < n_nodes; ++count_nodes) 
		(*perm)[n_nodes-1-count_nodes] = tperm[count_nodes]; 
	
	
// 	for (count_nodes = 0; count_nodes < n_nodes; ++count_nodes) 
// 		printf("Node: %d is in level %d\n", count_nodes, tperm[count_nodes]);
	
	printf("******************************************\n");
	
	for (count_nodes = 0; count_nodes < n_nodes; ++count_nodes) 
		if (tperm[count_nodes] == 2147483647) printf("Node: %d is in level %d\n", count_nodes, tperm[count_nodes]);
	
// 	free(g);
}