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

#define BFS_WORK_CHUNK 1024
#define BFS_PERCENT_CHUNK 0.5

#define THREAD_ON 1
#define THREAD_OFF 0
#define UNDEF_THREAD -1

#define INFINITY_LEVEL 2147483647
#define MIN_PRIORITY -9999999

#define UNDEF_NODE -1
#define ORPHAN_NODE -1

#define SLOAN_W1 1
#define SLOAN_W2 2

#ifndef __REORDERING_H__
#define __REORDERING_H__

typedef enum { OFF, ON } UPDATE;

#endif