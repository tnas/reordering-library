#include <stdlib.h>
#include <stdio.h>

#ifndef __LINKED_LIST_H__
#define __LINKED_LIST_H__

/* ***************************************************************
 * **************** A simple static queue ************************
 * ***************************************************************/

// Queues Q have a fixed maximum size N (rows, cols of the matrix) and are
// stored in an array. qh and qt point to queue head and tail.
inline static void QUEUE_enque(int** queue, const int size, int* tail_index, const int value)
{
	(*queue)[*tail_index] = value;
	*tail_index = (*tail_index + 1) % size;
}

// Dequeue operation (removes a node from the head)
inline static int QUEUE_deque(int** queue, int size, int* head_index)
{
	int value = (*queue)[*head_index];
	*head_index = (*head_index + 1) % size;
	return value;
}

// Predicate (queue empty)
#define QUEUE_empty(queue, head, tail) ((head) == (tail))


/* ***************************************************************
 * ***************** A Dynamic FIFO queue ************************
 * ***************************************************************/

typedef struct node
{
	int          data;
	int	     size; // it does not make sense
	int	     value;
	struct node* next;
} LIST;


LIST*    LIST_insert_IF_NOT_EXIST (LIST* L, int x);
LIST*    LIST_remove              (LIST* L, int x);
LIST*    LIST_remove_first        (LIST* L);
void     LIST_print               (LIST* L);
int      LIST_first               (LIST* L);
void     LIST_destroy             (LIST* L);
LIST*    LIST_add_IF_NOT_EXIST	  (LIST* list, int node, int val);
int   	 LIST_contains 		  (LIST* list, int element);

/* ***************************************************************
 * ***************** A Dynamic Array List ************************
 * ***************************************************************/

#define NON_ELEMENT -2147483647

typedef struct
{
	LIST* first;
	LIST* last;
	int   size;
} ARRAY_LIST;

void 	ARRAY_LIST_init		  (ARRAY_LIST** array_list);
void 	ARRAY_LIST_insert	  (ARRAY_LIST** array_list, int data);
int 	ARRAY_LIST_remove_first   (ARRAY_LIST** array_list);
void 	ARRAY_LIST_destroy	  (ARRAY_LIST** array_list);

#endif