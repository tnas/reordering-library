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

// #ifdef HAVE_CONFIG_H
#include <config.h>
// #endif

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
 * 		 liboctave-dev
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



// Transpose of the structure of a square sparse matrix
static void transpose (octave_idx_type N, const octave_idx_type *ridx,
            const octave_idx_type *cidx, octave_idx_type *ridx2,
            octave_idx_type *cidx2)
{
	octave_idx_type nz = cidx[N];
	
	OCTAVE_LOCAL_BUFFER (octave_idx_type, w, N + 1);
	for (octave_idx_type i = 0; i < N; i++)
		w[i] = 0;
	
	for (octave_idx_type i = 0; i < nz; i++)
		w[ridx[i]]++;
	
	nz = 0;
	
	for (octave_idx_type i = 0; i < N; i++)
	{
		OCTAVE_QUIT;
		cidx2[i] = nz;
		nz += w[i];
		w[i] = cidx2[i];
	}
	
	cidx2[N] = nz;
	w[N] = nz;
	
	for (octave_idx_type j = 0; j < N; j++)
		for (octave_idx_type k = cidx[j]; k < cidx[j + 1]; k++)
		{
			OCTAVE_QUIT;
			octave_idx_type q = w[ridx[k]]++;
			ridx2[q] = j;
		}
}


/**
 * It does not work.
 */
double Octave_RCM(METAGRAPH* mgraph, int** perm, int root) 
{
	octave_idx_type *cidx2 = mgraph->mat->JA;
	octave_idx_type *ridx2 = mgraph->mat->IA;
	
	// size of heaps
	octave_idx_type s = 0;
	
	// head- and tail-indices for the queue
	octave_idx_type qt = 0, qh = 0;
	
	CMK_Node v, w;
	
	// dimension of the matrix
	octave_idx_type N = mgraph->size;
	
	OCTAVE_LOCAL_BUFFER (octave_idx_type, cidx, N + 1);
	OCTAVE_LOCAL_BUFFER (octave_idx_type, ridx, cidx2[N]);
	transpose (N, ridx2, cidx2, ridx, cidx);
	
	// the permutation vector
	(*perm) = (int*) calloc(N, sizeof(int));
	
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
	bool *visit = btmp.fortran_vec();
	
	do
	{
		// locate an unvisited starting node of the graph
		octave_idx_type i;
		for (i = 0; i < N; i++)
			if (! visit[i]) break;
		
		// locate a probably better starting node
		v.id = root;
		
		// mark the node as visited and enqueue it (a starting node
		// for the BFS). Since the node will be a root of a spanning
		// tree, its dist is 0.
		v.deg = mgraph->graph[v.id].degree;
		v.dist = 0;
		visit[v.id] = true;
		Q_enq (Q, N, qt, v);
		
		// lower bound for the bandwidth of a subgraph
		// keep a "level" in the spanning tree (= min. distance to the
		// root) for determining the bandwidth of the computed
		// permutation P
		octave_idx_type Bsub = 0;
		// min. dist. to the root is 0
		octave_idx_type level = 0;
		// the root is the first/only node on level 0
		octave_idx_type level_N = 1;
		
		 while (! Q_empty (Q, N, qh, qt))
		{
			v = Q_deq (Q, N, qh);
			i = v.id;
		
			c++;
		
			// for computing the inverse permutation P where
			// A(inv(P),inv(P)) or P'*A*P is banded
			//         P(i) = c;
		
			// for computing permutation P where
			// A(P(i),P(j)) or P*A*P' is banded
			(*perm)[c] = i;
		
			// put all unvisited neighbors j of node i on the heap
			s = 0;
			octave_idx_type j1 = cidx[i];
			octave_idx_type j2 = cidx2[i];
		
			OCTAVE_QUIT;
			while (j1 < cidx[i+1] || j2 < cidx2[i+1])
			{
				OCTAVE_QUIT;
				if (j1 == cidx[i+1])
				{
					octave_idx_type r2 = ridx2[j2++];
					
					if (! visit[r2])
					{
						// the distance of node j is dist(i)+1
						w.id = r2;
						w.deg = mgraph->graph[r2].degree;
						w.dist = v.dist+1;
						H_insert (S, s, w);
						visit[r2] = true;
					}
				}
				else if (j2 == cidx2[i+1])
				{
					octave_idx_type r1 = ridx[j1++];
					
					if (! visit[r1])
					{
						w.id = r1;
						w.deg = mgraph->graph[r1].degree;
						w.dist = v.dist+1;
						H_insert (S, s, w);
						visit[r1] = true;
					}
				}
				else
				{
					octave_idx_type r1 = ridx[j1];
					octave_idx_type r2 = ridx2[j2];
					
					if (r1 <= r2)
					{
						if (! visit[r1])
						{
							w.id = r1;
							w.deg = mgraph->graph[r1].degree;
							w.dist = v.dist+1;
							H_insert (S, s, w);
							visit[r1] = true;
						}
						
						j1++;
						
						if (r1 == r2) j2++;
					}
					else
					{
						if (! visit[r2])
						{
							w.id = r2;
							w.deg = mgraph->graph[r2].degree;
							w.dist = v.dist+1;
							H_insert (S, s, w);
							visit[r2] = true;
						}
						
						j2++;
					}
				}
			}
			
			// add the neighbors to the queue (sorted by node degree)
			while (! H_empty (S, s))
			{
				OCTAVE_QUIT;
			
				// locate a neighbor of i with minimal degree in O(log(N))
				v = H_remove_min (S, s, 1);
			
				// entered the BFS a new level?
				if (v.dist > level)
				{
					// adjustment of bandwith:
					// "[...] the minimum bandwidth that
					// can be obtained [...] is the
					//  maximum number of nodes per level"
					if (Bsub < level_N) Bsub = level_N;
			
					level = v.dist;
					// v is the first node on the new level
					level_N = 1;
				}
				else
				{
					// there is no new level but another node on
					// this level:
					level_N++;
				}
			
				// enqueue v in O(1)
				Q_enq (Q, N, qt, v);
			}
		
			// synchronize the bandwidth with level_N once again:
			if (Bsub < level_N) Bsub = level_N;
		}
		
		// finish of BFS. If there are still unvisited nodes in the graph
		// then it is split into CCs. The computed bandwidth is the maximum
		// of all subgraphs. Update:
		if (Bsub > B) B = Bsub;
	}
	// are there any nodes left?
	while (c+1 < N);
 
	// compute the reverse-ordering
	s = N / 2 - 1;
	for (octave_idx_type i = 0, j = N - 1; i <= s; i++, j--)
		std::swap ((*perm)[i], (*perm)[j]);
	
	std::cout << "Permutation array: ";
	
	for (int i = 0; i < N; ++i)
		printf("%d ", (*perm)[i]);fflush(stdout);
	std::cout << "\n";
	
	return 0;
}

