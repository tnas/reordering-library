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


/* ***************************************************************
 * **************** A simple GRAPH static queue ******************
 * ***************************************************************/
/*
inline static void GRAPH_enque(GRAPH** queue, const int size, int* tail_index, const GRAPH value)
{
	(*queue)[*tail_index] = value;
	*tail_index = (*tail_index + 1) % size;
}

// Dequeue operation (removes a node from the head)
inline static GRAPH GRAPH_deque(GRAPH** queue, int size, int* head_index)
{
	GRAPH value = (*queue)[*head_index];
	*head_index = (*head_index + 1) % size;
	return value;
}*/


/*----------------------------------------------------------------------------
 * GRAPH FUNCTIONS PROTOTYPE IN CSR FORMAT
 *--------------------------------------------------------------------------*/
int      GRAPH_degree             (MAT* A, int x);
int*     GRAPH_adjacent           (MAT* A, int x);
int 	 GRAPH_degree_per_level   (MAT* A, int x, const int* levels, const int adjacency_level, const int* colors);
GRAPH*   GRAPH_adjacent_per_level (MAT* A, int x, const int* levels, const int adjacency_level, const int* colors);
void     GRAPH_bfs                (MAT* A, int x, int* dist);
int*     GRAPH_bfs_RCM            (MAT* A, int x, int* dist);
int      GRAPH_LS_depth           (int* LS, int n);
int      GRAPH_LS_width           (int* LS, int n);
LIST*    GRAPH_LS_last_level      (MAT* A, int* LS, int n);
int*     GRAPH_LS_peripheral      (MAT* A, int *node_s, int* node_e);
int* 	 GRAPH_neighboors 	  (MAT* A, int node, int degree);

int      COMPARE_degr_ASC         (const void * a, const void * b);
int      COMPARE_dist_degr_DES    (const void * a, const void * b);
int      COMPARE_dist_degr_ASC    (const void * a, const void * b);

#endif