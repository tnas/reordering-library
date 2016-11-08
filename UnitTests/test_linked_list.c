/*
 * Copyright 2016 Thiago Nascimento - nascimenthiago@gmail.com
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

#include "test_linked_list.h"

void test_array_list_init()
{
	ARRAY_LIST* list;
	
	ARRAY_LIST_init(&list);
	
	assert(list->first == NULL);
	assert(list->last == NULL);
	assert(list->size == 0);
	
	ARRAY_LIST_destroy(&list);
	
	printf("[test_array_list_init] --- OK.\n");
	fflush(stdout);
}


void test_array_list_insert_remove_one()
{
	ARRAY_LIST* list;
	int data = 10;

	ARRAY_LIST_init(&list);
	
	// Inserting data
	ARRAY_LIST_insert(&list, data);
	assert(list->size == 1);
	assert(list->first->data == data);
	assert(list->first == list->last);
	assert(list->last->next == NULL);
	
	// Removing data
	ARRAY_LIST_remove_first(&list);
	assert(list->size == 0);
	assert(list->last == NULL);
	assert(list->last == list->first);
	
	ARRAY_LIST_destroy(&list);
	
	printf("[test_array_list_insert_remove_one] --- OK.\n");
	fflush(stdout);
}



void test_array_list_insert_remove_multiples()
{
	ARRAY_LIST* list;
	int removed;
	int data1 = 10;
	int data2 = 20;
	int data3 = 10;

	ARRAY_LIST_init(&list);
	
	// Inserting data
	ARRAY_LIST_insert(&list, data1);
	ARRAY_LIST_insert(&list, data2);
	ARRAY_LIST_insert(&list, data3);
	assert(list->size == 3);
	assert(list->first->data == data1);
	assert(list->last->next == NULL);
	
	// Removing data1
	removed = ARRAY_LIST_remove_first(&list);
	assert(removed == data1);
	assert(list->size == 2);
	assert(list->first->data == data2);
	assert(list->last->data == data3);
	assert(list->last->next == NULL);
	
	// Removing data2
	removed = ARRAY_LIST_remove_first(&list);
	assert(removed == data2);
	assert(list->size == 1);
	assert(list->first->data == data3);
	assert(list->last->data == data3);
	assert(list->last->next == NULL);
	
	// Removing data3
	removed = ARRAY_LIST_remove_first(&list);
	assert(removed == data3);
	assert(list->size == 0);
	assert(list->last == NULL);
	assert(list->last == list->first);
	
	ARRAY_LIST_destroy(&list);
	
	printf("[test_array_list_insert_remove_multiples] --- OK.\n");
	fflush(stdout);
}


void test_array_list_add_desc_order()
{
	ARRAY_LIST* list;
	int removed;
	int data1 = 10;
	int data2 = 20;
	int data3 = 30;

	ARRAY_LIST_init(&list);
	
	// Inserting data
	ARRAY_LIST_add_desc_order(&list, data1);
	ARRAY_LIST_add_desc_order(&list, data2);
	ARRAY_LIST_add_desc_order(&list, data3);
	assert(list->size == 3);
	assert(list->first->data == data3);
	assert(list->last->data == data1);
	assert(list->last->next == NULL);
	
	// Removing data3
	removed = ARRAY_LIST_remove_first(&list);
	assert(removed == data3);
	assert(list->size == 2);
	assert(list->first->data == data2);
	assert(list->last->data == data1);
	assert(list->last->next == NULL);
	
	// Removing data2
	removed = ARRAY_LIST_remove_first(&list);
	assert(removed == data2);
	assert(list->size == 1);
	assert(list->first->data == data1);
	assert(list->last->data == data1);
	assert(list->last->next == NULL);
	
	// Removing data1
	removed = ARRAY_LIST_remove_first(&list);
	assert(removed == data1);
	assert(list->size == 0);
	assert(list->last == NULL);
	assert(list->last == list->first);
	
	ARRAY_LIST_destroy(&list);
	
	printf("[test_array_list_add_desc_order] --- OK.\n");
	fflush(stdout);
}



void run_all_linked_list_tests()
{
	test_array_list_init();
	test_array_list_insert_remove_one();
	test_array_list_insert_remove_multiples();
	test_array_list_add_desc_order();
}
