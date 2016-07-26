/*
 * Copyright 2016 Thiago Nascimento <nascimenthiago@gmail.com>
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
#ifndef __SLOAN_HSL_H__
#define __SLOAN_HSL_H__

#include "hsl.h"
#include "reordering.h"
#include "../CommonFiles/matrix.h"

long int REORDERING_SLOAN_HSL                (MAT* A);
long int REORDERING_SLOAN_pseudodiameter_HSL (MAT* A, int start_node, int end_node);

#endif /* __SLOAN_HSL_H__ */
