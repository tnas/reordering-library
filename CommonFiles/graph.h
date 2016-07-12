/*
 * Copyright 2016 Brenno Lugon brennolugon@gmail.com
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#include "matrix.h"
#include "linked_list.h"


#ifndef __GRAPH_H__
#define __GRAPH_H__

typedef enum LABEL { UNREACHED, LABELED } LABEL;

typedef struct 
{
	int label;
	int degree;
	int* neighboors;
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