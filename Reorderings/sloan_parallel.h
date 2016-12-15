/*
 * Copyright 2016 Thiago Nascimento nascimenthiago@gmail.com
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

#include "../CommonFiles/graph_parallel.h"
#include "reordering.h"

#define SLOAN_PRIORITY_FACTOR 10
#define SLOAN_CURR_PRIOR 0
#define SLOAN_NEW_PRIOR 1

typedef enum 
{ 
	INACTIVE, PREACTIVE, ACTIVE, NUMBERED 
	
} SLOAN_STATE;

typedef struct 
{
	int head;
	int tail;
	int* elements;
	ARRAY_LIST* list_elements;
	omp_lock_t lock;
} bag;


typedef struct 
{
	int label;
	int status;
	int distance;
	int* priorities;
} SLOAN_GRAPH;



inline static void GRAPH_enque(SLOAN_GRAPH** queue, const int size, int* tail_index, const SLOAN_GRAPH value)
{
	(*queue)[*tail_index] = value;
	*tail_index = (*tail_index + 1) % size;
}

// Dequeue operation (removes a node from the head)
inline static SLOAN_GRAPH GRAPH_deque(SLOAN_GRAPH** queue, int size, int* head_index)
{
	SLOAN_GRAPH value = (*queue)[*head_index];
	*head_index = (*head_index + 1) % size;
	return value;
}

void Parallel_Logical_Bag_Sloan (METAGRAPH* mgraph, int** permutation, int start_node, int end_node);
void Parallel_Bag_Sloan	        (METAGRAPH* mgraph, int** permutation, int start_node, int end_node);
void Parallel_Sloan	        (METAGRAPH* mgraph, int** permutation, int start_node, int end_node);