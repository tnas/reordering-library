/*----------------------------------------------------------------------------
 * UNORDERED RCM REORDERING SOLVER
 *--------------------------------------------------------------------------*/
#include "../CommonFiles/protos.h"

/*----------------------------------------------------------------------------
 * Unordered RCM reordering from the LEVEL STRUCTURE 
 *--------------------------------------------------------------------------*/
void UNORDERED_RCM (MAT* A, int** Fp)
{ 
  int s;
  int* q = calloc (A->n,sizeof(int));
  graph_unordered_bfs(A, s, q);
}
