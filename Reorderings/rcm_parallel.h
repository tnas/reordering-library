#include "../CommonFiles/graph_parallel.h"
#include "../CommonFiles/util_parallel.h"
#include "reordering.h"


extern void	Unordered_RCM(MAT* A, int** perm, int root, const float percent_chunk);
extern void	Leveled_RCM(MAT* mat, int** perm, int root);
extern void	Bucket_RCM(MAT* mat, int** perm, int root);