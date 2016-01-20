#include <omp.h>
#include "protos.h"
#include "util.h"

#define THREAD_ON 1
#define THREAD_OFF 0
#define UNDEF_NODE -1
#define UNDEF_THREAD -1
#define INFINITY_LEVEL 2147483647
#define BFS_WORK_CHUNK 1024


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


extern void	GRAPH_parallel_fixedpoint_bfs(MAT* adjacency, int root, int** levels, const float percent_chunk);
extern int* 	GRAPH_LS_peripheral_PARALLEL (MAT* A, int *node_s, int* node_e);
extern int 	GRAPH_LS_depth_PARALLEL(int* LS, int n);
extern int 	GRAPH_LS_width_PARALLEL(int* LS, int n);
extern LIST* 	GRAPH_LS_last_level_PARALLEL (MAT* A, int* LS, int n);


extern void	Unordered_RCM(MAT* A, int** perm, int root, const float percent_chunk);
extern void 	prefix_sum(const int* counts, int** sums, const int max_level);