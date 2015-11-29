#include <omp.h>
#include "protos.h"

#define NUM_THREADS 4
#define THREAD_ON 1
#define THREAD_OFF 0

#define iseven(n) ((n)%(2)==(0)?(1):(0))

typedef enum OPERATION {READ, WRITE} OPERATION;


extern int*	GRAPH_parallel_fixedpoint_bfs(MAT* adjacency, int root, int* levels);