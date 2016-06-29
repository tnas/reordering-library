#include "hsl.h"
#include "reordering.h"
#include "../CommonFiles/graph.h"


extern void     REORDERING_SLOAN            (MAT* A, int** Fp, int node_s, int node_e);
extern long int	REORDERING_SLOAN_HSL 	    (MAT* A);