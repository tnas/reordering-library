/*----------------------------------------------------------------------------
 * UNORDERED RCM REORDERING SOLVER
 *--------------------------------------------------------------------------*/
#include "../CommonFiles/protos.h"

/*----------------------------------------------------------------------------
 * Unordered RCM reordering from the LEVEL STRUCTURE 
 *--------------------------------------------------------------------------*/
void REORDERING_RCM_parallel(MAT* A, int** Fp)
{ 
  int s;
  int* q = calloc (A->n,sizeof(int));
  GRAPH_unordered_bfs(A, s, q);
}
