#include <stdlib.h>
#include <stdio.h>

#ifndef __LINKED_LIST_H__
#define __LINKED_LIST_H__

// A simple static queue. 

// Queues Q have a fixed maximum size N (rows, cols of the matrix) and are
// stored in an array. qh and qt point to queue head and tail.
inline static void QUEUE_enque(int** queue, const int size, int* tail_index, const int value)
{
	(*queue)[*tail_index] = value;
	*tail_index = (*tail_index + 1) % (size + 1);
}

// Dequeue operation (removes a node from the head)
inline static int QUEUE_deque(int** queue, int size, int* head_index)
{
	int value = (*queue)[*head_index];
	*head_index = (*head_index + 1) % (size + 1);
	return value;
}

// Predicate (queue empty)
#define QUEUE_empty(queue, head, tail) ((head) == (tail))


typedef struct node
{
	int          data;
	int	     size;
	int	     value;
	struct node* next;
} LIST;

extern LIST*    LIST_insert_IF_NOT_EXIST (LIST* L, int x);
extern LIST*    LIST_remove              (LIST* L, int x);
extern LIST*    LIST_remove_first        (LIST* L);
extern void     LIST_print               (LIST* L);
extern int      LIST_first               (LIST* L);
extern void     LIST_destroy             (LIST* L);
extern LIST*    LIST_add_IF_NOT_EXIST	 (LIST* list, int node, int val);
extern int 	LIST_contains 		 (LIST* list, int element);

#endif