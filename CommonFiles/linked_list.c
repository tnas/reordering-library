/*----------------------------------------------------------------------------
 * LINKED LIST FUNCTIONS
 *--------------------------------------------------------------------------*/
#include "linked_list.h"

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
	LIST* new_cell;
	
	cell = list;
	
	if (cell == NULL)
	{
		new_cell = (LIST*) malloc (sizeof(LIST));
		new_cell->data = node;
		new_cell->next = NULL;
		new_cell->value = val;
		new_cell->next = NULL;
		
		return new_cell;
	}
	
	for (; cell->next != NULL; cell = cell->next)
	{
		if (cell->data == node)
		{
			cell->value += val;
			return list;
		}
	}
	
	if (cell->data == node)
	{
		cell->value += val;
		return list;
	}
	
	// Creating a new element when list is NULL
	new_cell = (LIST*) malloc (sizeof(LIST));
	new_cell->data = node;
	new_cell->next = NULL;
	new_cell->value = val;
	new_cell->next = NULL;
	
	cell->next = new_cell;
	
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
 * Returns 1 if the LIST contains the specified element
 *--------------------------------------------------------------------------*/
int LIST_contains (LIST* list, int element)
{
	LIST* cell = list;
	
	if (list == NULL) return 0;
	
	while (cell != NULL && cell->data != element) cell = cell->next;

	if (cell == NULL) return 0;		
	else return 1;
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
		return;
	}
	LIST *P;
	while (L != NULL)
	{
		P = L->next;
		L->size--;
		free(L);
		L = P;
	}
}


/*----------------------------------------------------------------------------
 * Initialize an ARRAY_LIST structure. 
 * 
 * @since 03-11-2016
 *--------------------------------------------------------------------------*/
void ARRAY_LIST_init(ARRAY_LIST** array_list)
{
	*array_list = malloc(sizeof(ARRAY_LIST));
	(*array_list)->first = 
	(*array_list)->last  = NULL;
	(*array_list)->size  = 0;
}


/*----------------------------------------------------------------------------
 * Insert an element in the ARRAY_LIST structure. 
 * Duplicate elements are permitted.
 * 
 * @since 03-11-2016
 *--------------------------------------------------------------------------*/
void ARRAY_LIST_insert(ARRAY_LIST** array_list, int data)
{
	LIST* new_node;
	
	if ((*array_list)->first == NULL)
	{
		(*array_list)->first = malloc(sizeof(LIST));
		(*array_list)->first->data = data;
		(*array_list)->first->next = NULL;
		(*array_list)->last  = (*array_list)->first;
		(*array_list)->size = 1;
	}
	else 
	{
		new_node = (*array_list)->last;
		new_node->next = malloc(sizeof(LIST));
		new_node->next->data = data;
		new_node->next->next = NULL;
		(*array_list)->last = new_node->next;
		(*array_list)->size++;
	}
}


/*----------------------------------------------------------------------------
 * Insert an element in the ARRAY_LIST structure in descending order . 
 * Duplicate elements are no permitted.
 * 
 * @since 04-11-2016
 *--------------------------------------------------------------------------*/
void ARRAY_LIST_add_desc_order(ARRAY_LIST** array_list, int data)
{
	LIST* new_node;
	LIST* searcher;
	
	if ((*array_list)->first == NULL)
	{
		ARRAY_LIST_insert(array_list, data);
	}
	else 
	{
		if (data == (*array_list)->first->data)
		{
			return;
		}
		else if (data > (*array_list)->first->data)
		{
			new_node = malloc(sizeof(LIST));
			new_node->data = data;
			new_node->next = (*array_list)->first;
			(*array_list)->first = new_node;
			(*array_list)->size++;
			
			return;
		}
		else {
			for (searcher = (*array_list)->first; searcher->next != NULL; searcher = searcher->next)
			{
				if (data == searcher->next->data)
				{
					return;
				}
				else if (data > searcher->next->data)
				{
					// Inserting in the middle of array list
					new_node = malloc(sizeof(LIST));
					new_node->data = data;
					new_node->next = searcher->next;
					searcher->next = new_node;
					(*array_list)->size++;
					
					return;
				}
			}
			
			// Inserting as new last position
			new_node = malloc(sizeof(LIST));
			new_node->data = data;
			new_node->next = NULL;
			searcher->next = new_node;
			(*array_list)->last = new_node;
			(*array_list)->size++;
		}
	}
}


/*----------------------------------------------------------------------------
 * Remove the first element from ARRAY_LIST structure. 
 * If there is no elements (ARRAY_LIST is null) a NON_ELEMENT is
 * returned.
 * 
 * @since 03-11-2016
 *--------------------------------------------------------------------------*/
int ARRAY_LIST_remove_first(ARRAY_LIST** array_list)
{
	int data;
	LIST* garbage;
	
	if ((*array_list)->first == NULL)
	{
		data = -NON_ELEMENT;
	}
	else 
	{
		data = (*array_list)->first->data;
		garbage = (*array_list)->first;
		(*array_list)->first = garbage->next;
		garbage->next = NULL;
		(*array_list)->size--;
		free(garbage);
		
		if ((*array_list)->size == 0) 
			(*array_list)->last = NULL;
	}
	
	return data;
}


/*----------------------------------------------------------------------------
 * Destroy an ARRAY_LIST structure. 
 * 
 * @since 04-11-2016
 *--------------------------------------------------------------------------*/
void ARRAY_LIST_destroy(ARRAY_LIST** array_list)
{
	LIST_destroy((*array_list)->first);
	free(*array_list);
}
