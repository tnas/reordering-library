#include "matrix.h"
#include "linked_list.h"


#ifndef __GRAPH_H__
#define __GRAPH_H__

typedef enum LABEL { UNREACHED, LABELED } LABEL;

typedef struct 
{
	int label;
	int degree;
	int distance;
	int parent;
	int status;
	int chnum;
	int index;
} GRAPH;


/*----------------------------------------------------------------------------
 * GRAPH FUNCTIONS PROTOTYPE IN CSR FORMAT
 *--------------------------------------------------------------------------*/
extern int      GRAPH_degree             (MAT* A, int x);
extern int*     GRAPH_adjacent           (MAT* A, int x);
extern int 	GRAPH_degree_per_level   (MAT* A, int x, const int* levels, const int adjacency_level, const int* colors);
extern GRAPH* 	GRAPH_adjacent_per_level (MAT* A, int x, const int* levels, const int adjacency_level, const int* colors);
extern void     GRAPH_bfs                (MAT* A, int x, int* dist);
extern int*     GRAPH_bfs_RCM            (MAT* A, int x, int* dist);
extern int      GRAPH_LS_depth           (int* LS, int n);
extern int      GRAPH_LS_width           (int* LS, int n);
extern LIST*    GRAPH_LS_last_level      (MAT* A, int* LS, int n);
extern int*     GRAPH_LS_peripheral      (MAT* A, int *node_s, int* node_e);

extern int      COMPARE_degr_ASC         (const void * a, const void * b);
extern int      COMPARE_dist_degr_DES    (const void * a, const void * b);
extern int      COMPARE_dist_degr_ASC    (const void * a, const void * b);

#endif