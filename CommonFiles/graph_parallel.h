/*
 * Copyright 2016 Thiago Nascimento <nascimenthiago@gmail.com>
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
#include <omp.h>
#include "graph.h"

#ifndef __GRAPH_PARALLEL_H__
#define __GRAPH_PARALLEL_H__

#define INFINITY_LEVEL 2147483647
#define NON_VERTEX -1


typedef struct
{
	int num_children;
	GRAPH* children;
} genealogy;

typedef struct
{
	int start;
	int end;
} graph_diameter;


typedef struct
{
	GRAPH* graph;
	int vertex_min_degree;
	int size;
} METAGRAPH; 

void       GRAPH_parallel_fixedpoint_bfs (MAT* adjacency, int root, int** levels, const float percent_chunk);
METAGRAPH* GRAPH_parallel_build		 (MAT* mat);

#endif