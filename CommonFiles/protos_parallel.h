#include <omp.h>
#include "protos.h"

#define NUM_THREADS 4
#define THREAD_ON 1
#define THREAD_OFF 0

#define iseven(n) ((n)%(2)==(0)?(1):(0))
#define isdivisor(d, n) ((n)%(d)==(0)?(1):(0))

typedef enum OPERATION {READ, WRITE} OPERATION;

typedef struct 
{
	int initial_prefix_sum;
	int curr_prefix_sum;
	int curr_total_sum;
	int last_prefix_sum;
	int last_total_sum;
} status_prefix_sum;


extern int*	GRAPH_parallel_fixedpoint_bfs(MAT* adjacency, int root, int* levels);