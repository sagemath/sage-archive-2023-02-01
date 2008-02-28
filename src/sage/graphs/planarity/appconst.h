/********************************************************************
Copyright 2005 John M. Boyer

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
 ********************************************************************/

#ifndef APPCONST_H
#define APPCONST_H

//#define PROFILE
#ifdef PROFILE
#include "platformTime.h"
#endif

/* If the DEBUG macro is not defined, then low-level functions are replaced by faster macros */

#ifndef DEBUG
#define SPEED_MACROS
#endif

/* Return status values */

#define OK              0
#define NOTOK           -2
#define NONPLANAR       -3
#define NONOUTERPLANAR	-4

/* Array indices are used as pointers, and this means bad pointer */

#define NIL		-1
#define NIL_CHAR	0xFF

/* Defines fopen strings for reading and writing text files on PC and UNIX */

#ifdef WINDOWS
#define READTEXT        "rt"
#define WRITETEXT       "wt"
#else
#define READTEXT        "r"
#define WRITETEXT       "w"
#endif

/* This macro controls whether DFS puts all child edges at the beginning of
    the adjacency list (link[0]) and all forward edges at the end of the
    adjacency list (link[1]).  This allows several subsequent methods to
    find the desired edges more easily and to terminate without doing
    as much unnecessary traversal of adjacency lists */

#define ORDER_EDGES

#endif
