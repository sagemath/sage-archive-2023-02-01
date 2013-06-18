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

#ifndef STACK_H
#define STACK_H

#ifdef __cplusplus
extern "C" {
#endif

// includes mem functions like memcpy
#include <string.h>

typedef struct
{
        int *S;
        int size, capacity;
} stack;

typedef stack * stackP;

stackP sp_New(int);
void sp_Free(stackP *);

int  sp_Copy(stackP, stackP);

int  sp_CopyContent(stackP stackDst, stackP stackSrc);
stackP sp_Duplicate(stackP theStack);

#define sp_GetCapacity(theStack) (theStack->capacity)

#ifndef SPEED_MACROS

int  sp_ClearStack(stackP);
int  sp_GetCurrentSize(stackP theStack);
int  sp_SetCurrentSize(stackP theStack, int top);

int  sp_IsEmpty(stackP);
int  sp_NonEmpty(stackP);

#define sp_Push(theStack, a) { if (sp__Push(theStack, (a)) != OK) return NOTOK; }
#define sp_Push2(theStack, a, b) { if (sp__Push2(theStack, (a), (b)) != OK) return NOTOK; }

int  sp__Push(stackP, int);
int  sp__Push2(stackP, int, int);

#define sp_Pop(theStack, a) { if (sp__Pop(theStack, &(a)) != OK) return NOTOK; }
#define sp_Pop2(theStack, a, b) { if (sp__Pop2(theStack, &(a), &(b)) != OK) return NOTOK; }

int  sp__Pop(stackP, int *);
int  sp__Pop2(stackP, int *, int *);

int  sp_Top(stackP);
int  sp_Get(stackP, int);
int  sp_Set(stackP, int, int);

#else

#define sp_ClearStack(theStack) theStack->size=0
#define sp_GetCurrentSize(theStack) (theStack->size)
#define sp_SetCurrentSize(theStack, Size) ((Size) > theStack->capacity ? NOTOK : (theStack->size = (Size), OK))

#define sp_IsEmpty(theStack) !theStack->size
#define sp_NonEmpty(theStack) theStack->size

#define sp_Push(theStack, a) theStack->S[theStack->size++] = a
#define sp_Push2(theStack, a, b) {sp_Push(theStack, a); sp_Push(theStack, b);}

#define sp_Pop(theStack, a) a=theStack->S[--theStack->size]
#define sp_Pop2(theStack, a, b) {sp_Pop(theStack, b);sp_Pop(theStack, a);}

#define sp_Top(theStack) (theStack->size ? theStack->S[theStack->size-1] : NIL)
#define sp_Get(theStack, pos) (theStack->S[pos])
#define sp_Set(theStack, pos, val) (theStack->S[pos] = val)

#endif

#ifdef __cplusplus
}
#endif

#endif
