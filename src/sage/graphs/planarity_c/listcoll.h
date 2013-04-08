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

#ifndef _LISTCOLL_H
#define _LISTCOLL_H

#ifdef __cplusplus
extern "C" {
#endif

/* This include is needed for memset and memcpy */
#include <string.h>

typedef struct
{
        int prev, next;
} lcnode;

typedef struct
{
        int N;
        lcnode *List;
} listCollectionRec;

typedef listCollectionRec * listCollectionP;

listCollectionP LCNew(int N);
void LCFree(listCollectionP *pListColl);

void LCInsertAfter(listCollectionP listColl, int theAnchor, int theNewNode);
void LCInsertBefore(listCollectionP listColl, int theAnchor, int theNewNode);

#ifndef SPEED_MACROS

void LCReset(listCollectionP listColl);
void LCCopy(listCollectionP dst, listCollectionP src);

int  LCGetNext(listCollectionP listColl, int theList, int theNode);
int  LCGetPrev(listCollectionP listColl, int theList, int theNode);

int  LCPrepend(listCollectionP listColl, int theList, int theNode);
int  LCAppend(listCollectionP listColl, int theList, int theNode);
int  LCDelete(listCollectionP listColl, int theList, int theNode);

#else

/* void LCReset(listCollectionP listColl); */

#define LCReset(listColl) memset(listColl->List, NIL_CHAR, listColl->N*sizeof(lcnode))

/* void LCCopy(listCollectionP dst, listCollectionP src) */

#define LCCopy(dst, src) memcpy(dst->List, src->List, src->N*sizeof(lcnode))

/* int  LCGetNext(listCollectionP listColl, int theList, int theNode);
	Return theNode's successor, unless it is theList head pointer */

#define LCGetNext(listColl, theList, theNode) listColl->List[theNode].next==theList ? NIL : listColl->List[theNode].next

/* int  LCGetPrev(listCollectionP listColl, int theList, int theNode);
	Return theNode's predecessor unless theNode is theList head.
	To start going backwards, use NIL for theNode, which returns theList head's predecessor
	Usage: Obtain last node, loop while NIL not returned, process node then get predecessor.
		After theList head processed, get predecessor returns NIL because we started with
		theList head's predecessor. */

#define LCGetPrev(listColl, theList, theNode) \
        (theNode==NIL \
	 ? listColl->List[theList].prev \
         : theNode==theList ? NIL : listColl->List[theNode].prev)

/* int  LCPrepend(listCollectionP listColl, int theList, int theNode);
    After an append, theNode is last, which in a circular list is the direct predecessor
	of the list head node, so we just back up one. For singletons, this has no effect.*/

#define LCPrepend(listColl, theList, theNode) listColl->List[LCAppend(listColl, theList, theNode)].prev

/* int  LCAppend(listCollectionP listColl, int theList, int theNode);
	If theList is empty, then theNode becomes its only member and is returned.
	Otherwise, theNode is placed before theList head, which is returned. */

#define LCAppend(listColl, theList, theNode) \
        (theList==NIL \
         ? (listColl->List[theNode].prev = listColl->List[theNode].next = theNode) \
         : (listColl->List[theNode].next = theList, \
            listColl->List[theNode].prev = listColl->List[theList].prev, \
            listColl->List[listColl->List[theNode].prev].next = theNode, \
            listColl->List[theList].prev = theNode, \
	    theList))

/* int  LCDelete(listCollectionP listColl, int theList, int theNode);
	If theList contains only one node, then NIL it out and return NIL meaning empty list
	Otherwise, join the predecessor and successor, then
	return either the list head or its successor if the deleted node is the list head
	(in that case, the caller makes the successor become the new list head).*/


#define LCDelete(listColl, theList, theNode) \
        listColl->List[theList].next == theList \
	? (listColl->List[theList].prev = listColl->List[theList].next = NIL) \
        : (listColl->List[listColl->List[theNode].prev].next = listColl->List[theNode].next, \
           listColl->List[listColl->List[theNode].next].prev = listColl->List[theNode].prev, \
	   (theList==theNode ? listColl->List[theNode].next : theList))

#endif

#ifdef __cplusplus
}
#endif

#endif
