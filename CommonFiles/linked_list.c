/*----------------------------------------------------------------------------
 * LINKED LIST FUNCTIONS
 *--------------------------------------------------------------------------*/
#include "protos.h"

/*----------------------------------------------------------------------------
 * Insert a element (if not already exist) in the LIST structure 
 *--------------------------------------------------------------------------*/
LIST* LIST_insert_IF_NOT_EXIST (LIST* L, int x)
{
	/* creating new list element */
	if (L == NULL)
	{
		LIST *N = (LIST*) malloc (sizeof(LIST));
		N->data = x;
		N->next = NULL;
		N->size = 1;
		N->value = 0;
		return N;		
	}
		
	LIST *P;
	
	/* if already exist, return */
	for (P = L; P->next != NULL; P = P->next)
	{
		if (P->data == x)
			return L;
	}
	if (P->data == x) 
		return L;

	LIST *N = (LIST*) malloc (sizeof(LIST));
	N->data = x;
	N->next = NULL;
	N->value = 0;
	
	P->next = N;
	L->size++;
	
	return L;
}


/*----------------------------------------------------------------------------
 * Add the value val to the element cell of the LIST list 
 *--------------------------------------------------------------------------*/
LIST* LIST_add_IF_NOT_EXIST(LIST* list, int node, int val)
{
	LIST* cell;
	
	// if already exist, update and return
	for (cell = list; cell != NULL; cell = cell->next)
	{
		if (cell->data == node)
		{
			cell->value += val;
			return list;
		}
	}
	
	// Creating a new element when list is NULL
	cell = (LIST*) malloc (sizeof(LIST));
	cell->data = node;
	cell->next = NULL;
	cell->value = val;
	cell->next = NULL;
	
	return list;
}


/*----------------------------------------------------------------------------
 * Remove the element x from the LIST structure
 *--------------------------------------------------------------------------*/
LIST* LIST_remove (LIST* L, int x)
{
	if (L == NULL)
	{
		printf("warning: Empty LIST. Returning.. [LIST_remove]\n");
		return L;
	}
	LIST* Q = NULL;
	LIST* P = L;
	while (P != NULL && P->data != x)
	{
		Q = P;
		P = P->next;
	}

	if (P == NULL)
	{
		printf ("warning: Element %d does not exist in this LIST. Returning.. [LIST_remove]\n",x);
		return L;		
	}
	
	L->size--;
	
	if (Q == NULL)
		L = P->next;
	else
		Q->next = P->next;
	
	free(P);
	
	return L;
}


/*----------------------------------------------------------------------------
 * Print all elements in the LIST structure
 *--------------------------------------------------------------------------*/
void LIST_print (LIST* L)
{
	if (L == NULL)
	{
		printf("warning: Empty LIST. Returning.. [LIST_print]\n");
		return;
	}
	LIST *current = L;
	while (current != NULL)
	{
		printf ("%d ", current->data);
		current = current->next;
	}
	printf("\n");
}

/*----------------------------------------------------------------------------
 * Return the first element in the LIST structure
 *--------------------------------------------------------------------------*/
int LIST_first (LIST* L)
{
	return L->data;
}

/*----------------------------------------------------------------------------
 * Destroy the LIST structure from memory
 *--------------------------------------------------------------------------*/
void LIST_destroy (LIST* L)
{
	if (L == NULL)
	{
		printf("warning: Empty LIST. Returning.. [LIST_destroy]\n");
		return;
	}
	LIST *P;
	while (L != NULL)
	{
		P = L->next;
		free(L);
		L->size--;
		L = P;
	}
}
