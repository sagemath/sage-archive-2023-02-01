#ifndef APPCONST_H
#define APPCONST_H

/*
Planarity-Related Graph Algorithms Project
Copyright (c) 1997-2010, John M. Boyer
All rights reserved. Includes a reference implementation of the following:

* John M. Boyer. "Simplified O(n) Algorithms for Planar Graph Embedding,
  Kuratowski Subgraph Isolation, and Related Problems". Ph.D. Dissertation,
  University of Victoria, 2001.

* John M. Boyer and Wendy J. Myrvold. "On the Cutting Edge: Simplified O(n)
  Planarity by Edge Addition". Journal of Graph Algorithms and Applications,
  Vol. 8, No. 3, pp. 241-273, 2004.

* John M. Boyer. "A New Method for Efficiently Generating Planar Graph
  Visibility Representations". In P. Eades and P. Healy, editors,
  Proceedings of the 13th International Conference on Graph Drawing 2005,
  Lecture Notes Comput. Sci., Volume 3843, pp. 508-511, Springer-Verlag, 2006.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of the Planarity-Related Graph Algorithms Project nor the names
  of its contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// When PROFILE is defined, prints out run-time stats on a number of subordinate
// routines in the embedder

//#define PROFILE
#ifdef PROFILE
#include "platformTime.h"
#endif

/* Define DEBUG to get additional debugging. The default is to define it when MSC does */

#ifdef _DEBUG
#define DEBUG
#endif

/* Some low-level functions are replaced by faster macros, except when debugging */

#define SPEED_MACROS
#ifdef DEBUG
#undef SPEED_MACROS
#endif

/* Return status values; OK/NOTOK behave like Boolean true/false,
   not like program exit codes. */

#define OK              1
#define NOTOK           0

#ifdef DEBUG
#undef NOTOK
extern int debugNOTOK();
#include <stdio.h>
#define NOTOK           (printf("NOTOK on Line %d of %s\n", __LINE__, __FILE__), debugNOTOK())
#endif

#define NONEMBEDDABLE   -3

#ifndef TRUE
#define TRUE            1
#endif

#ifndef FALSE
#define FALSE           0
#endif

#ifndef NULL
#define NULL			0L
#endif

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

#endif
