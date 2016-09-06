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

#include "rcm_octave.h"

#include <config.h>
#include "ov.h"
#include "defun-dld.h"
#include "error.h"
#include "gripes.h"
#include "utils.h"
#include "oct-locbuf.h"
#include "ov-re-mat.h"
#include "ov-re-sparse.h"
#include "ov-cx-sparse.h"
#include "oct-sparse.h"


/**
 * Dependencies: libsuitesparse-dev
 */

typedef int octave_idx_type;

// A node struct for the Cuthill-McKee algorithm
typedef struct 
{
	// the node's id (matrix row index)
	octave_idx_type id;
	// the node's degree
	octave_idx_type deg;
	// minimal distance to the root of the spanning tree
	octave_idx_type dist;
} CMK_Node;


// A simple queue. 
// Queues Q have a fixed maximum size N (rows, cols of the matrix) and are
// stored in an array. qh and qt point to queue head and tail.
inline static void Q_enq(CMK_Node *Q, octave_idx_type N, octave_idx_type& qt, const CMK_Node& o)
{
	Q[qt] = o;
	qt = (qt + 1) % (N + 1);
}

// Dequeue operation (removes a node from the head)
inline static CMK_Node Q_deq(CMK_Node *Q, octave_idx_type N, octave_idx_type& qh)
{
	CMK_Node r = Q[qh];
	qh = (qh + 1) % (N + 1);
	return r;
}

// Predicate (queue empty)
#define Q_empty(Q, N, qh, qt)  ((qh) == (qt))


// A simple, array-based binary heap (used as a priority queue for nodes)

// the left descendant of entry i
#define LEFT(i)		(((i) << 1) + 1)		// = (2*(i) + 1)
//the right descendant of entry i
#define RIGHT(i)	(((i) << 1) + 2)		// = (2*(i) + 2)
//the parent of entry i
#define PARENT(i)	(((i) - 1) >> 1)		// = floor(((i)-1)/2)

// Builds a min-heap (the root contains the smallest element). A is an array
// with the graph's nodes, i is a starting position, size is the length of A.
 
static void H_heapify_min (CMK_Node *A, octave_idx_type i, octave_idx_type size)
{
   octave_idx_type j = i;
   
   for (;;)
   {
       octave_idx_type l = LEFT(j);
       octave_idx_type r = RIGHT(j);
 
       octave_idx_type smallest;
       if (l < size && A[l].deg < A[j].deg)
         smallest = l;
       else
         smallest = j;
 
       if (r < size && A[r].deg < A[smallest].deg)
         smallest = r;
 
       if (smallest != j)
         {
           CMK_Node tmp = A[j];
           A[j] = A[smallest];
           A[smallest] = tmp;
           j = smallest;
         }
       else 
         break;
     }
}
 
 // Heap operation insert. Running time is O(log(n))
 
static void H_insert (CMK_Node *H, octave_idx_type& h, const CMK_Node& o)
{
   octave_idx_type i = h++;
 
   H[i] = o;
 
   if (i == 0) 
     return;
   do
     {
       octave_idx_type p = PARENT(i);
       if (H[i].deg < H[p].deg)
         {
           CMK_Node tmp = H[i];
           H[i] = H[p];
           H[p] = tmp;
 
           i = p;
         }
       else 
         break;
     }
   while (i > 0);
 }
 
 // Heap operation remove-min. Removes the smalles element in O(1) and
 // reorganizes the heap optionally in O(log(n))
 
inline static CMK_Node H_remove_min (CMK_Node *H, octave_idx_type& h, int reorg/*=1*/)
{
   CMK_Node r = H[0];
   H[0] = H[--h];
   if (reorg) 
     H_heapify_min(H, 0, h);
   return r;
}
 
// Predicate (heap empty)
#define H_empty(H, h)   ((h) == 0)


// Calculates the node's degrees. 
static octave_idx_type calc_degrees (METAGRAPH* mgraph)
{
	octave_idx_type max_deg = 0;
	
	for (octave_idx_type i = 0; i < mgraph->size; i++)
	{
		if (mgraph->graph[i].degree > max_deg)
			max_deg = mgraph->graph[i].degree;
	}
   
	return max_deg;
}


double Octave_RCM(METAGRAPH* mgraph, int** perm, int root) 
{
	// size of heaps
	octave_idx_type s = 0;
	
	// head- and tail-indices for the queue
	octave_idx_type qt = 0, qh = 0;
	
	CMK_Node v, w;
	
	// dimension of the matrix
	octave_idx_type N = mgraph->size;
	
	// the permutation vector
// 	NDArray P (dim_vector (1, N));
	
	octave_idx_type max_deg = calc_degrees(mgraph);
	
	// a heap for the a node's neighbors. The number of neighbors is
	// limited by the maximum degree max_deg:
	OCTAVE_LOCAL_BUFFER (CMK_Node, S, max_deg);
	
	// a queue for the BFS. The array is always one element larger than
	// the number of entries that are stored.
	OCTAVE_LOCAL_BUFFER (CMK_Node, Q, N+1);
	
	// a counter (for building the permutation)
	octave_idx_type c = -1;
	
	// upper bound for the bandwidth (=quality of solution)
	// initialize the bandwidth of the graph with 0. B contains the
	// the maximum of the theoretical lower limits of the subgraphs
	// bandwidths.
	octave_idx_type B = 0;
	
	// mark all nodes as unvisited; with the exception of the nodes
	// that have degree==0 and build a CC of the graph.
	
	boolNDArray btmp (dim_vector (1, N), false);
	bool *visit = btmp.fortran_vec ();
	
}

