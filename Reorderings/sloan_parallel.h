#include "../CommonFiles/graph_parallel.h"
#include "reordering.h"

#define SLOAN_PRIORITY_FACTOR 10
#define SLOAN_CURR_PRIOR 0
#define SLOAN_NEW_PRIOR 1

typedef enum { 
	INACTIVE, PREACTIVE, ACTIVE, NUMBERED 
	
} SLOAN_STATE;

void Parallel_Sloan (MAT* adjacency, int** permutation, int start_node, int end_node);