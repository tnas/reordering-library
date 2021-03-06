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

#ifndef SLOAN_BOOST_H
#define SLOAN_BOOST_H

#ifdef __cplusplus
extern "C" {
#endif

#include "../CommonFiles/graph_parallel.h"
	
double Boost_Sloan(const METAGRAPH* mgraph, int** permutation, int start_node, int end_node);

#ifdef __cplusplus
}
#endif

#endif // SLOAN_BOOST_H
