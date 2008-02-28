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

#define GRAPHNONPLANAR_C

#include "graph.h"

/* Imported functions */

extern void _ClearIsolatorContext(graphP theGraph);
extern void _FillVisitedFlags(graphP, int);
extern void _FillVisitedFlagsInBicomp(graphP theGraph, int BicompRoot, int FillValue);
extern void _SetVertexTypeInBicomp(graphP theGraph, int BicompRoot, int theType);

extern int  _GetNextVertexOnExternalFace(graphP theEmbedding, int curVertex, int *pPrevLink);
extern int  _GetPertinentChildBicomp(graphP theEmbedding, int W);
extern void _WalkDown(graphP theEmbedding, int I, int RootVertex);
extern void _OrientVerticesInEmbedding(graphP theEmbedding);
extern void _OrientVerticesInBicomp(graphP theEmbedding, int BicompRoot, int PreserveSigns);

/* Private functions (exported to system) */

int  _ChooseTypeOfNonplanarityMinor(graphP theEmbedding, int I, int R);
int  _InitializeNonplanarityContext(graphP theEmbedding, int I, int R);

int  _FindNonplanarityBicompRoot(graphP theEmbedding);
void _FindActiveVertices(graphP theEmbedding, int R, int *pX, int *pY);
int  _FindPertinentVertex(graphP theEmbedding);

void _PopAndUnmarkVerticesAndEdges(graphP theEmbedding, int Z);

int  _MarkHighestXYPath(graphP theEmbedding);
int  _MarkZtoRPath(graphP theEmbedding);
int  _FindExtActivityBelowXYPath(graphP theEmbedding);

/****************************************************************************
 _ChooseTypeOfNonplanarityMinor()
 ****************************************************************************/

int  _ChooseTypeOfNonplanarityMinor(graphP theEmbedding, int I, int R)
{
int  N, X, Y, W, Px, Py, Z, DFSChild, RootId;

/* Create the initial non-planarity minor state in the isolator context */

     if (_InitializeNonplanarityContext(theEmbedding, I, R) != OK)
         return NOTOK;

     N = theEmbedding->N;
     R = theEmbedding->IC.r;
     X = theEmbedding->IC.x;
     Y = theEmbedding->IC.y;
     W = theEmbedding->IC.w;

/* If the root copy is not a root copy of the current vertex I,
        then the Walkdown terminated because it couldn't find
        a viable path along a child bicomp, which is Minor A. */

     if (theEmbedding->V[R - N].DFSParent != I)
     {
         theEmbedding->IC.minorType |= FLAGS_MINOR_A;
         return OK;
     }

/* If W has an externally active pertinent child bicomp, then
     we've found Minor B */

     if (theEmbedding->V[W].pertinentBicompList != NIL)
     {
         RootId = LCGetPrev(theEmbedding->BicompLists,
                            theEmbedding->V[W].pertinentBicompList, NIL);
         DFSChild = RootId;
         if (theEmbedding->V[DFSChild].Lowpoint < I)
         {
             theEmbedding->IC.minorType |= FLAGS_MINOR_B;
             return OK;
         }
     }

/* Find the highest obstructing X-Y path */

     if (_MarkHighestXYPath(theEmbedding) != OK)
         return NOTOK;

     Px = theEmbedding->IC.px;
     Py = theEmbedding->IC.py;

/* If either point of attachment is 'high' (P_x closer to R than X
     or P_y closer to R than Y along external face), then we've
     matched Minor C. */

     if (theEmbedding->G[Px].type == VERTEX_HIGH_RXW ||
         theEmbedding->G[Py].type == VERTEX_HIGH_RYW)
     {
            theEmbedding->IC.minorType |= FLAGS_MINOR_C;
            return OK;
     }

/* For Minor D, we search for a path from an internal
     vertex Z along the X-Y path up to the root R of the bicomp. */

     if (_MarkZtoRPath(theEmbedding) == NOTOK)
         return NOTOK;

     if (theEmbedding->IC.z != NIL)
     {
         theEmbedding->IC.minorType |= FLAGS_MINOR_D;
         return OK;
     }

/* For Minor E, we search for an externally active vertex Z
     below the points of attachment of the X-Y path */

     Z = _FindExtActivityBelowXYPath(theEmbedding);
     if (Z != NIL)
     {
         theEmbedding->IC.z = Z;
         theEmbedding->IC.minorType |= FLAGS_MINOR_E;
         return OK;
     }

     return NOTOK;
}

/****************************************************************************
 _InitializeNonplanarityContext()
 ****************************************************************************/

int  _InitializeNonplanarityContext(graphP theEmbedding, int I, int R)
{
int  e, X, Y, W, Z, ZPrevLink, ZType;
int  singleBicompMode =  (R == NIL) ? 0 : 1;

/* For the embedding or in a given bicomp, orient the vertices,
    and clear the visited members of all vertex and edge records. */

     if (!singleBicompMode)
     {
         _OrientVerticesInEmbedding(theEmbedding);
         _FillVisitedFlags(theEmbedding, 0);
     }
     else
     {
         _OrientVerticesInBicomp(theEmbedding, R, 1);
         _FillVisitedFlagsInBicomp(theEmbedding, R, 0);
     }

/* Blank out the isolator context, then assign the input graph reference
     and the current vertext I into the context. */

     _ClearIsolatorContext(theEmbedding);
     theEmbedding->IC.v = I;

/* Now we find a root copy R of the current vertex on which the Walkdown failed
   (i.e. there is an unembedded back edge between an ancestor of the current
   vertex and descendant of the current vertex in the subtree rooted by
   the DFS child C=R-N. */

     R = _FindNonplanarityBicompRoot(theEmbedding);
     if (R == NIL) return NOTOK;

/* Now we find the active vertices along both external face paths extending
     from R. If either is not a stopping vertex, then we call Walkdown to
     reconstruct the stack to the root of the descendant bicomp that blocked
     the Walkdown. Otherwise, R is the desired root. Either way, the stopping
     vertices x and y are calculated. */

     _FindActiveVertices(theEmbedding, R, &X, &Y);

     if (theEmbedding->V[X].pertinentBicompList != NIL ||
         theEmbedding->V[Y].pertinentBicompList != NIL)
     {
          /* We are dealing with minor A, which has a main bicomp rooted by
                a descendant of R.  So we put the bicomp rooted by R back
                the way we found it. */

          if (singleBicompMode)
              _OrientVerticesInBicomp(theEmbedding, R, 1);

          /* Now we do the Walkdown to find the descendant bicomp. */

          _WalkDown(theEmbedding, I, R);
          if (sp_IsEmpty(theEmbedding->theStack))
              return NOTOK;

          sp_Pop2(theEmbedding->theStack, R, e);

          /* Now we give the proper orientation to the descendant bicomp,
                but we make it reversible (last param=1) to make processing
                easier for the K3,3 search, which simply reinitializes
                everything when it finds minor A. */

          if (singleBicompMode)
              _OrientVerticesInBicomp(theEmbedding, R, 1);

          /* Now we find the stopping vertices along the external face
                paths of the descendant bicomp. */

          _FindActiveVertices(theEmbedding, R, &X, &Y);
     }

/* We store the info we have gathered so far in the isolator context */

     theEmbedding->IC.r = R;
     theEmbedding->IC.x = X;
     theEmbedding->IC.y = Y;

/* Now, we obtain the pertinent vertex W on the lower external face
    path between X and Y (that path that does not include R). */

     theEmbedding->IC.w = W = _FindPertinentVertex(theEmbedding);

/* In the current bicomp, we clear the type flags */

/* Now we can classify the vertices along the external face of the bicomp
    rooted at R as 'high RXW', 'low RXW', 'high RXY', 'low RXY' */

     _SetVertexTypeInBicomp(theEmbedding, R, TYPE_UNKNOWN);

     ZPrevLink = 1;
     Z = _GetNextVertexOnExternalFace(theEmbedding, R, &ZPrevLink);
     ZType = VERTEX_HIGH_RXW;
     while (Z != W)
     {
         if (Z == X) ZType = VERTEX_LOW_RXW;
         theEmbedding->G[Z].type = ZType;
         Z = _GetNextVertexOnExternalFace(theEmbedding, Z, &ZPrevLink);
     }

     ZPrevLink = 0;
     Z = _GetNextVertexOnExternalFace(theEmbedding, R, &ZPrevLink);
     ZType = VERTEX_HIGH_RYW;
     while (Z != W)
     {
         if (Z == Y) ZType = VERTEX_LOW_RYW;
         theEmbedding->G[Z].type = ZType;
         Z = _GetNextVertexOnExternalFace(theEmbedding, Z, &ZPrevLink);
     }

/* All work is done, so return success */

     return OK;
}

/****************************************************************************
 _FindNonplanarityBicompRoot()

 This procedure finds the root copy R of the current vertex on which the
 Walkdown failed (whether it failed while traversing the bicomp rooted by
 R or some descendant bicomp is determined later).

 We iterate the forward cycle edges of the vertex I looking for a forward
 edge (I, W) that was not embedded.  Once it is found, we figure out which
 bicomp rooted by a root copy of I contains W or contains a DFS ancestor of W.

 This turns out to be an easy test.  The desired bicomp is rooted by the DFS
 tree edge (I, C) with the largest value of C that does not exceed W.  C is
 a DFS ancestor of Z.

 Return: The desired root copy, or NIL on error.
 ****************************************************************************/

int  _FindNonplanarityBicompRoot(graphP theEmbedding)
{
int  R, tempChild, fwdArc, W=NIL, C=NIL, I=theEmbedding->IC.v;

/* Obtain the forward arc of an unembedded back edge from I to one of its
    descendants (edges are removed from the forward arc list as they are
    embedded, so the list will be empty if all edges were embedded). */

    if ((fwdArc = theEmbedding->V[I].fwdArcList) == NIL)
        return NIL;

    W = theEmbedding->G[fwdArc].v;

/* Find the greatest DFS child C of I that is less than W.  This will
    give us the ancestor of W that is a child of I.  Since the
    ancestors of I have not been processed by the planarity algorithm,
    the separatedDFSChildList of I contains all the children of I. */

    tempChild = theEmbedding->V[I].separatedDFSChildList;

    while (tempChild != NIL)
    {
        if (tempChild > C && tempChild < W)
            C = tempChild;

        tempChild = LCGetNext(theEmbedding->DFSChildLists,
                              theEmbedding->V[I].separatedDFSChildList, tempChild);
    }

    if (C == NIL) return NIL;

/* The root vertex of a bicomp rooted by edge (I, C) is located at
        position C+N in our data structures */

     R = C + theEmbedding->N;
     return R;
}

/****************************************************************************
 _FindStoppingVertices()

 Descends from the root of a bicomp R in both the link[0] and link[1]
 directions, returning the first active vertex appearing in either direction.
 ****************************************************************************/

void _FindActiveVertices(graphP theEmbedding, int R, int *pX, int *pY)
{
int  XPrevLink=1, YPrevLink=0, I=theEmbedding->IC.v;

     *pX = _GetNextVertexOnExternalFace(theEmbedding, R, &XPrevLink);
     *pY = _GetNextVertexOnExternalFace(theEmbedding, R, &YPrevLink);

     while (_VertexActiveStatus(theEmbedding, *pX, I) == VAS_INACTIVE)
        *pX = _GetNextVertexOnExternalFace(theEmbedding, *pX, &XPrevLink);

     while (_VertexActiveStatus(theEmbedding, *pY, I) == VAS_INACTIVE)
        *pY = _GetNextVertexOnExternalFace(theEmbedding, *pY, &YPrevLink);
}

/****************************************************************************
 _FindPertinentVertex()

 Get the first vertex after x. Since x was obtained using a prevlink of 1 on r,
 we use the same prevlink so we don't go back to r (works because all vertices
 have the same orientation).
 Then, we proceed around the lower path until we find a vertex W that either
 has pertinent child bicomps or is directly adjacent to the current vertex I.
 ****************************************************************************/

int  _FindPertinentVertex(graphP theEmbedding)
{
int  W=theEmbedding->IC.x, WPrevLink=1, I=theEmbedding->IC.v;

     W = _GetNextVertexOnExternalFace(theEmbedding, W, &WPrevLink);

     while (W != theEmbedding->IC.y)
     {
         if (PERTINENT(theEmbedding, W, I))
             return W;

         W = _GetNextVertexOnExternalFace(theEmbedding, W, &WPrevLink);
     }

     return NIL;
}

/****************************************************************************
 _PopAndUnmarkVerticesAndEdges()

 Pop all vertex/edge pairs from the top of the stack up to a terminating
 vertex Z and mark as unvisited.  If Z is NIL, then all vertex/edge pairs
 are popped and marked as unvisited.
 ****************************************************************************/

void _PopAndUnmarkVerticesAndEdges(graphP theEmbedding, int Z)
{
int  V, e;

     while (sp_NonEmpty(theEmbedding->theStack))
     {
            sp_Pop(theEmbedding->theStack, e);

            /* If we popped a vertex other than the termination vertex Z, then
                we also pop the edge we pushed, and we clear the visited flags
                for the vertex and the edge's two edge records. */

            if (e < 2*theEmbedding->N && e != Z)
            {
                V = e;
                sp_Pop(theEmbedding->theStack, e);
                theEmbedding->G[V].visited = 0;
                theEmbedding->G[e].visited = 0;
                theEmbedding->G[gp_GetTwinArc(theEmbedding, e)].visited = 0;
            }

            /* If we popped an edge or the terminating vertex Z, then put it
                back and break */

            else
            {
                sp_Push(theEmbedding->theStack, e);
                break;
            }
     }
}

/****************************************************************************
 _MarkHighestXYPath()

 An X-Y path in the bicomp rooted by R is a path attached to the external
 face at points Px and Py that separates W from R such that a back edge (R, W)
 cannot be embedded within the bicomp. Recall that R is a root copy of I, so
 (R, W) is the representative of (I, W).  Also, note that W is pertinent if
 either W *or* one of its descendants in a separate bicomp has, in the input
 graph, a back edge to I.

 If no X-Y path separating W from R is found, then NOTOK is returned because
 this procedure is only called once an X-Y path is guaranteed (by our proof
 of correctness) to exist.

 The desired output is to set the 'visited' flags of the X-Y path with
 highest points of attachment to the external face (i.e. the points of
 attachment that are closest to R along the external face).  This includes
 marking both the vertices and edges along the X-Y path.

 As a function of initialization, the vertices along the external face
 (other than R and W) have been classified as 'high RXW', 'low RXW', 'high RXY',
 or 'low RXY'. Once the vertices have been categorized, we proceed with trying
 to set the visitation flags in the way described above.  To do this, we first
 remove all edges incident to R except the two edges that join R to the external
 face. The result is that R and its two remaining edges are a 'corner' in the
 external face but also in a single proper face whose boundary includes the
 X-Y path with the highest attachment points. Thus, we simply need to walk
 this proper face to find the desired X-Y path. Note, however, that the
 resulting face boundary may have attached cut vertices.  Any such separable
 component contains a vertex neighbor of R, but the edge to R has been
 temporarily removed.  The algorithm removes loop of vertices and edges along
 the proper face so that only a path is identified.

 To walk the proper face containing R, we begin with its link[0] successor,
 then take the link[1] corner at every subsequent turn.  For each vertex,
 we mark as visited the vertex as well as the edge used to enter the vertex
 (except for the edge used to enter the RXW vertex).  We also push the visited
 vertices and edges onto a stack.

 As we walk the proper face, we keep track of the last vertex P_x we visited of
 type RXW (high or low).  Each time we encounter a vertex of type RXW (high or
 low), we pop the stack and unmark all of the edges and vertices visited because
 they were part of a path parallel to the external face that does not obstruct
 W from reaching R within the bicomp.  If we encounter vertex W, then there is
 no obstructing X-Y path since we removed only edges incident to R, so we pop
 the stack unmarking everything then return NOTOK as stated above.  If we
 encounter a vertex Z previously visited, then we pop the stack, unmarking the
 vertices and edges popped, until we find the prior occurence of Z on the stack.

 Otherwise, the first time we encounter a vertex P_y of type 'RYW', we stop
 because the obstructing X-Y path has been marked visited and its points of
 connection P_x and P_y have been found.

 Once the X-Y path is identified, we restore the edges incident to R.

 The stack space used is no greater than 3N.  The first N accounts for removing
 the edges incident to R.  The other 2N accounts for the fact that each
 iteration of the main loop visits a vertex, pushing its index and the
 location of an edge record.  If a vertex is encountered that is already
 on the stack, then it is not pushed again (and in fact part of the stack
 is removed).
 ****************************************************************************/

int  _MarkHighestXYPath(graphP theEmbedding)
{
int J, Z, e;
int R, X, Y, W;

/* Initialization */

     R = theEmbedding->IC.r;
     X = theEmbedding->IC.x;
     Y = theEmbedding->IC.y;
     W = theEmbedding->IC.w;
     theEmbedding->IC.px = theEmbedding->IC.py = NIL;

     sp_ClearStack(theEmbedding->theStack);

/* Remove the internal edges incident to vertex R */

     J = theEmbedding->G[R].link[0];
     J = theEmbedding->G[J].link[0];
     while (J != theEmbedding->G[R].link[1])
     {
          sp_Push(theEmbedding->theStack, J);
          gp_HideEdge(theEmbedding, J);
          J = theEmbedding->G[J].link[0];
     }

/* Walk the proper face containing R to find and mark the highest
        X-Y path. Note that if W is encountered, then there is no
        intervening X-Y path, so we would return NOTOK. */

     Z = R;
     J = theEmbedding->G[R].link[1];
     while (theEmbedding->G[Z].type != VERTEX_HIGH_RYW &&
            theEmbedding->G[Z].type != VERTEX_LOW_RYW)
     {
          /* Advance J and Z along the proper face containing R */

          J = theEmbedding->G[J].link[1];
          if (J < 2*theEmbedding->N)
              J = theEmbedding->G[J].link[1];
          Z = theEmbedding->G[J].v;
          J = gp_GetTwinArc(theEmbedding, J);

          /* If Z is already visited, then pop everything since the last time
                we visited Z because its all part of a separable component. */

          if (theEmbedding->G[Z].visited)
          {
              _PopAndUnmarkVerticesAndEdges(theEmbedding, Z);
          }

          /* If we have not visited this vertex before... */

          else
          {
              /* If we found another vertex along the RXW path, then blow off
                 all the vertices we visited so far because they're not part of
                 the obstructing path */

              if (theEmbedding->G[Z].type == VERTEX_HIGH_RXW ||
                  theEmbedding->G[Z].type == VERTEX_LOW_RXW)
              {
                  theEmbedding->IC.px = Z;
                  _PopAndUnmarkVerticesAndEdges(theEmbedding, NIL);
              }

              /* Push the current vertex onto the stack of vertices visited
                 since the last RXW vertex was encountered */

              sp_Push(theEmbedding->theStack, J);
              sp_Push(theEmbedding->theStack, Z);

              /* Mark the vertex Z as visited as well as its edge of entry
                 (except the entry edge for P_x).*/

              theEmbedding->G[Z].visited = 1;
              if (Z != theEmbedding->IC.px)
              {
                  theEmbedding->G[J].visited = 1;
                  theEmbedding->G[gp_GetTwinArc(theEmbedding, J)].visited = 1;
              }

              /* If we found an RYW vertex, then we have successfully finished
                 identifying the highest X-Y path, so we record the point of
                 attachment and break the loop. */

              if (theEmbedding->G[Z].type == VERTEX_HIGH_RYW ||
                  theEmbedding->G[Z].type == VERTEX_LOW_RYW)
              {
                 theEmbedding->IC.py = Z;
                 break;
              }
          }
     }

/* Restore the internal edges incident to R that were previously removed,
        ignoring any leftover vertices that might be on the stack. */

     while (sp_NonEmpty(theEmbedding->theStack))
     {
          sp_Pop(theEmbedding->theStack, e);
          if (e < 2*theEmbedding->N)
               sp_Pop(theEmbedding->theStack, e);
          else gp_RestoreEdge(theEmbedding, e);
     }

/* Return the result */

     return theEmbedding->IC.py==NIL ? NOTOK : OK;
}

/****************************************************************************
 _MarkZtoRPath()

 This function assumes that _MarkHighestXYPath() has already been called,
 which marked as visited the vertices and edges along the X-Y path.

 We begin at the point of attachment P_x and traverse its link[1] edge records
 until we find one marked visited, which leads to the first internal vertex
 along the X-Y path.  We begin with this vertex (and its edge of entry), and
 we run until we find P_y.  For each internal vertex Z and its edge of entry
 ZPrevArc, we take the link[1] successor edge record of ZPrevArc (skipping Z
 if it intervenes).  This is called ZNextArc.  If ZNextArc is marked visited
 then it is along the X-Y path, so we use it to exit Z and go to the next
 vertex on the X-Y path.

 If ZNextArc is not visited, then when _MarkHighestXYPath() ran, it exited
 Z from ZNextArc, then eventually reentered Z.  In other words, Z became a
 cut vertex when we removed the internal edges incident to R. Thus, ZNextArc
 indicates the first edge in an internal path to R.

 When we find an unvisited ZNextArc, we stop running the X-Y path and instead
 begin marking the Z to R path.  We move to successive vertices using a
 twin arc then a link[1] successor edge record, only this time we have not
 removed the internal edges incident to R, so this technique does eventually
 lead us all the way to R.

 If we do not find an unvisited ZNextArc for any vertex Z on the X-Y path and
 inside the bicomp, then there is no Z to R path, so we return.
 ****************************************************************************/

int  _MarkZtoRPath(graphP theEmbedding)
{
int ZPrevArc, ZNextArc, Z, R, Px, Py;

/* Initialize */

    R = theEmbedding->IC.r;
    Px = theEmbedding->IC.px;
    Py = theEmbedding->IC.py;
    theEmbedding->IC.z = NIL;

/* Begin at Px and search its adjacency list for the edge leading to
   the first internal vertex of the X-Y path. */

    Z = Px;
    ZNextArc = theEmbedding->G[Z].link[1];
    while (ZNextArc != theEmbedding->G[Z].link[0])
    {
       if (theEmbedding->G[ZNextArc].visited) break;

       ZNextArc = theEmbedding->G[ZNextArc].link[1];
    }

    if (!theEmbedding->G[ZNextArc].visited)
        return NOTOK;

/* For each internal vertex Z, determine whether it has a path to root. */

    while (theEmbedding->G[ZNextArc].visited)
    {
        ZPrevArc = gp_GetTwinArc(theEmbedding, ZNextArc);
        ZNextArc = theEmbedding->G[ZPrevArc].link[1];
        if (ZNextArc < 2*theEmbedding->N)
            ZNextArc = theEmbedding->G[ZNextArc].link[1];
    }

    ZPrevArc = gp_GetTwinArc(theEmbedding, ZNextArc);
    Z = theEmbedding->G[ZPrevArc].v;

/* If there is no Z to R path, return */

    if (Z == Py) return OK;

/* Otherwise, store Z in the isolation context */

    theEmbedding->IC.z = Z;

/* Walk the proper face starting with (Z, ZNextArc) until we reach R, marking
        the vertices and edges encountered along the way, then Return OK. */

    while (Z != R)
    {
        /* If we ever encounter a non-internal vertex (other than the root R),
                then corruption has occured, so we return NOTOK */

        if (theEmbedding->G[Z].type != TYPE_UNKNOWN)
            return NOTOK;

        /* Go to the next vertex indicated by ZNextArc */

        Z = theEmbedding->G[ZNextArc].v;

        /* Mark the next vertex and the edge leading to it as visited. */

        theEmbedding->G[ZNextArc].visited = 1;
        theEmbedding->G[ZPrevArc].visited = 1;
        theEmbedding->G[Z].visited = 1;

        /* Go to the next edge in the proper face */

        ZNextArc = theEmbedding->G[ZPrevArc].link[1];
        if (ZNextArc < 2*theEmbedding->N)
            ZNextArc = theEmbedding->G[ZNextArc].link[1];

        ZPrevArc = gp_GetTwinArc(theEmbedding, ZNextArc);
    }

/* Found Z to R path, so indicate as much to caller */

    return OK;
}

/****************************************************************************
 _FindExtActivityBelowXYPath()

 Get an externally active vertex along the lower external face path between
 the points of attachment P_x and P_y of a 'low' X-Y Path.
 NOTE: By the time this function is called, Px and Py have already been found
        to be at or below X and Y.
 ****************************************************************************/

int  _FindExtActivityBelowXYPath(graphP theEmbedding)
{
int  Z=theEmbedding->IC.px, ZPrevLink=1,
     Py=theEmbedding->IC.py, I=theEmbedding->IC.v;

     Z = _GetNextVertexOnExternalFace(theEmbedding, Z, &ZPrevLink);

     while (Z != Py)
     {
         if (_VertexActiveStatus(theEmbedding, Z, I) == VAS_EXTERNAL)
             return Z;

         Z = _GetNextVertexOnExternalFace(theEmbedding, Z, &ZPrevLink);
     }

     return NIL;
}
