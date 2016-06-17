#include <omp.h>
#include "protos.h"
#include "util.h"

#define THREAD_ON 1
#define THREAD_OFF 0
#define UNDEF_THREAD -1

#define INFINITY_LEVEL 2147483647
#define MIN_PRIORITY -9999999

#define BFS_WORK_CHUNK 1024
#define BFS_PERCENT_CHUNK 0.5

#define UNDEF_NODE -1
#define ORPHAN_NODE -1

#define SLOAN_PRIORITY_FACTOR 10
#define SLOAN_CURR_PRIOR 0
#define SLOAN_NEW_PRIOR 1

#define iseven(n) ((n)%(2)==(0)?(1):(0))
#define isdivisor(d, n) ((n)%(d)==(0)?(1):(0))

typedef enum REORDERING_ALGORITHM { SERIAL_RCM, UNORDERED_RCM, LEVELED_RCM } 
	REORDERING_ALGORITHM;

typedef enum SLOAN_STATE { INACTIVE, PREACTIVE, ACTIVE, NUMBERED }
	SLOAN_STATE;

typedef enum UPDATE { OFF, ON } UPDATE;
	
typedef struct 
{
	int initial_prefix_sum;
	int curr_prefix_sum;
	int curr_total_sum;
	int last_prefix_sum;
	int last_total_sum;
} status_prefix_sum;


typedef struct
{
	int num_children;
	GRAPH* children;
} genealogy;


extern void	GRAPH_parallel_fixedpoint_bfs(MAT* adjacency, int root, int** levels, const float percent_chunk);
extern int* 	GRAPH_LS_peripheral_PARALLEL (MAT* A, int *node_s, int* node_e);
extern int 	GRAPH_LS_depth_PARALLEL(int* LS, int n);
extern int 	GRAPH_LS_width_PARALLEL(int* LS, int n);
extern LIST* 	GRAPH_LS_last_level_PARALLEL (MAT* A, int* LS, int n);


extern void	Unordered_RCM(MAT* A, int** perm, int root, const float percent_chunk);
extern void	Leveled_RCM(MAT* mat, int** perm, int root);
extern void	Bucket_RCM(MAT* mat, int** perm, int root);
extern void 	Parallel_Sloan (MAT* adjacency, int** Fp, int start_node, int end_node);
extern void 	prefix_sum(const int* counts, int** sums, const int max_level);
