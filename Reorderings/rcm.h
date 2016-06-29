#include "hsl.h"
#include "../CommonFiles/graph.h"



extern void     REORDERING_RCM_opt          (MAT* A, int** p, int s);
extern void     REORDERING_RCM              (MAT* A, int** p);
extern long int REORDERING_HSL_RCM          (MAT* A);