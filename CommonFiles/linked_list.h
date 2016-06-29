#include <stdlib.h>
#include <stdio.h>

#ifndef __LINKED_LIST_H__
#define __LINKED_LIST_H__

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