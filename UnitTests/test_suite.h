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

/* inclusion guard */
#ifndef __TEST_SUITE_H__
#define __TEST_SUITE_H__

#include "../CommonFiles/matrix.h"

typedef enum {
	ALL,
	REORDERING, 
	PREFIX_SUM,
	MATRIX,
	GRAPH_PARALLEL,
	UTIL_FUNCIONS,
	LINKED_LIST,
	TEMA
} test_suite;


typedef enum { 
	TEST_SUITE, 
	TEST_CASE 
} test_scope;


#endif /* __TEST_SUITE_H__ */