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

#define _LISTCOLL_C

#include "appconst.h"
#include "listcoll.h"
#include <stdlib.h>

/*****************************************************************************
 The data structure defined by this module manages a set of N objects
 arranged as a collection of circular lists, each containing distinct
 elements from the set.

 On construction, LCNew() creates an array of N nodes, each containing a
 prev and next pointer.  The identity of the node is given by its array index.
 Each node's prev and next pointers are set to NIL, indicating that the node
 is not currently part of a list.  LCReset() can be called to reset all
 pointers to NIL.

 The function LCFree() deallocates the collection of lists and clears the
 pointer variable used to pass the collection.

 An empty list is indicated by NIL.  To begin a list with node I, call
 LCPrepend() or LCAppend() with the NIL list and with I as the node.  The prev
 and next pointers in node I are set to I and I is returned as the head of
 the list.

 Future calls to LCPrepend() add a node J as the new first element of the list,
 so the list given as input is pointed to by J's next, and J is returned as
 the head of the list.

 Future calls to LCAppend() add a node J as the new last element, so the prev
 pointer of the list given as input will indicate node J, and the input list
 is returned as the head of the list.

 The function LCDelete() removes a node I from a list L.  If node I is in the
 list alone, then its pointers are set to NIL, and NIL is returned as the list.
 If node I is not alone in the list, but it is the head of the list (in other
 words, I is equal to L), then L's sucessor is returned as the new head of the
 list. Whether or not I equals L, node I is deleted by joining its predecessor
 and successor nodes.

 LCCopy() copies the contents of one collection to another if both are of
 equal size.

 LCGetNext() is used for forward iteration through a list in the collection.
 The expected iteration pattern is to process the node one has then call
 LCGetNext() to get the next node, so if the result of LCGetNext() would be the
 head of the list, then NIL is returned instead.  This simplifies most
 coding operations involving LCGetNext().

 LCGetPrev() is used for backward iteration through a list in the collection.
 The expected iteration pattern is that the last list element will be obtained
 by an initial call to LCGetPrev() with theNode equal to NIL.  This call
 should appear outside of the iteration loop.  The iteration loop then
 proceeds while the current node is not NIL.  The loop body processes the
 current node, then LCGetPrev() is called with theNode equal to the current
 node.  LCGetPrev() returns NIL if theNode is equal to theList.  Otherwise,
 the predecessor of theNode is returned.

 *****************************************************************************/

/*****************************************************************************
 LCNew()
 *****************************************************************************/

listCollectionP LCNew(int N)
{
listCollectionP theListColl = NULL;

     if (N <= 0) return theListColl;

     theListColl = (listCollectionP) malloc(sizeof(listCollectionRec));
     if (theListColl != NULL)
     {
         theListColl->List = (lcnode *) malloc(N*sizeof(lcnode));
         if (theListColl->List == NULL)
         {
             free(theListColl);
             theListColl = NULL;
         }
         else
         {
             theListColl->N = N;
			 LCReset(theListColl);
         }
     }
     return theListColl;
}

/*****************************************************************************
 LCFree()
 *****************************************************************************/

void LCFree(listCollectionP *pListColl)
{
     if (pListColl==NULL || *pListColl==NULL) return;

     if ((*pListColl)->List != NULL)
         free((*pListColl)->List);

     free(*pListColl);
     *pListColl = NULL;
}


#ifndef SPEED_MACROS

/*****************************************************************************
 LCReset()
 *****************************************************************************/

void LCReset(listCollectionP listColl)
{
int  I;

     for (I=0; I < listColl->N; I++)
          listColl->List[I].prev = listColl->List[I].next = NIL;
}

/*****************************************************************************
 LCCopy()
 *****************************************************************************/

void LCCopy(listCollectionP dst, listCollectionP src)
{
int  I;

     if (dst==NULL || src==NULL || dst->N != src->N) return;

     for (I=0; I<dst->N; I++)
          dst->List[I] = src->List[I];

}

/*****************************************************************************
 LCGetNext()
 *****************************************************************************/

int  LCGetNext(listCollectionP listColl, int theList, int theNode)
{
int  next;

     if (listColl==NULL || theList==NIL || theNode==NIL) return NIL;
     next = listColl->List[theNode].next;
     return next==theList ? NIL : next;
}

/*****************************************************************************
 LCGetPrev()
 *****************************************************************************/

int  LCGetPrev(listCollectionP listColl, int theList, int theNode)
{
     if (listColl==NULL || theList==NIL) return NIL;
     if (theNode == NIL) return listColl->List[theList].prev;
     if (theNode == theList) return NIL;
     return listColl->List[theNode].prev;
}

/*****************************************************************************
 LCPrepend()
 *****************************************************************************/

int  LCPrepend(listCollectionP listColl, int theList, int theNode)
{
int  newList = LCAppend(listColl, theList, theNode);

     /* If the append worked, then theNode is last, which in a circular
        list is the direct predecessor of the list head node, so we
        just back up one. For singletons, the result is unchanged. */

     if (newList != NOTOK)
         newList = listColl->List[newList].prev;

     return newList;
}

/*****************************************************************************
 LCAppend()
 *****************************************************************************/

int  LCAppend(listCollectionP listColl, int theList, int theNode)
{
     /* If the given list is empty, then the given node becomes the
        singleton list output */

     if (theList == NIL)
     {
         listColl->List[theNode].prev = listColl->List[theNode].next = theNode;
         theList = theNode;
     }

     /* Otherwise, make theNode the predecessor of head node of theList,
        which is where the last node goes in a circular list. */

     else
     {
     int pred = listColl->List[theList].prev;

         listColl->List[theList].prev = theNode;
         listColl->List[theNode].next = theList;
         listColl->List[theNode].prev = pred;
         listColl->List[pred].next = theNode;
     }

     /* Return the list (only really important if it was NIL) */

     return theList;
}

/*****************************************************************************
 LCDelete()
 *****************************************************************************/

int  LCDelete(listCollectionP listColl, int theList, int theNode)
{
     /* If the list is a singleton, then NIL its pointers and
        return NIL for theList*/

     if (listColl->List[theList].next == theList)
     {
         listColl->List[theList].prev = listColl->List[theList].next = NIL;
         theList = NIL;
     }

     /* Join predecessor and successor, dropping theNode from the list.
        If theNode is the head of the list, then return the successor as
        the new head node. */

     else
     {
     int pred=listColl->List[theNode].prev,
         succ=listColl->List[theNode].next;

         listColl->List[pred].next = succ;
         listColl->List[succ].prev = pred;

         listColl->List[theNode].prev = listColl->List[theNode].next = NIL;

         if (theList == theNode)
             theList = succ;
     }

     return theList;
}

#endif // SPEED_MACROS
