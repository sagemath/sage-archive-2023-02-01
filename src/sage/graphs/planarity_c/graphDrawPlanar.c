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

#include "graphDrawPlanar.h"
#include "graphDrawPlanar.private.h"

extern int DRAWPLANAR_ID;

#include "graph.h"

#if 0
#include <malloc.h>
#else
#if	!defined(_XOPEN_SOURCE) && !defined(_ISOC99_SOURCE)
#define	_XOPEN_SOURCE	600
#endif
#include <stdlib.h>		/* ISO C99 also defines malloc() to be there. */
#endif
#include <string.h>
#include <stdio.h>

extern void _FillVisitedFlags(graphP theGraph, int FillValue);

/* Private functions exported to system */

void _CollectDrawingData(DrawPlanarContext *context, int RootVertex, int W, int WPrevLink);
int  _BreakTie(DrawPlanarContext *context, int BicompRoot, int W, int WPrevLink);

int  _ComputeVisibilityRepresentation(DrawPlanarContext *context);
int  _CheckVisibilityRepresentationIntegrity(DrawPlanarContext *context);

/* Private functions */
int _ComputeVertexPositions(DrawPlanarContext *context);
int _ComputeVertexPositionsInComponent(DrawPlanarContext *context, int root, int *pIndex);
int _ComputeEdgePositions(DrawPlanarContext *context);
int _ComputeVertexRanges(DrawPlanarContext *context);
int _ComputeEdgeRanges(DrawPlanarContext *context);

/********************************************************************
 _ComputeVisibilityRepresentation()

  Compute vertex positions
  Compute edge positions
  Assign horizontal ranges of vertices
  Assign vertical ranges of edges

 ********************************************************************/

int _ComputeVisibilityRepresentation(DrawPlanarContext *context)
{
    if (sp_NonEmpty(context->theGraph->edgeHoles))
        return NOTOK;

    if (_ComputeVertexPositions(context) != OK)
        return NOTOK;

    if (_ComputeEdgePositions(context) != OK)
        return NOTOK;

    if (_ComputeVertexRanges(context) != OK)
        return NOTOK;

    if (_ComputeEdgeRanges(context) != OK)
        return NOTOK;

    return OK;
}

/********************************************************************
 _ComputeVertexPositions()

  Computes the vertex positions in the graph.  This method accounts
  for disconnected graphs by finding the DFS tree roots and then,
  for each, invoking _ComputeVertexPositionsInComponent().
  The index variable for the positioning is maintained by this method
  so that the vertices in separate components still get distinct
  vertex positions.
 ********************************************************************/

int _ComputeVertexPositions(DrawPlanarContext *context)
{
graphP theEmbedding = context->theGraph;
int I, index;

    for (I = index = 0; I < theEmbedding->N; I++)
    {
        // For each DFS tree root in the embedding, we
        // compute the vertex positions
        if (theEmbedding->V[I].DFSParent == NIL)
        {
            if (_ComputeVertexPositionsInComponent(context, I, &index) != OK)
                return NOTOK;
        }
    }

    return OK;
}

/********************************************************************
 _ComputeVertexPositionsInComponent()

  The vertical positions of the vertices are computed based in part
  on the information compiled during the planar embedding.

  Each vertex is marked as being between its parent and some ancestor
  or beyond the parent relative to the ancestor.  The localized,
  intuitive notion is that the vertex is either below the parent
  or above the parent, but the bicomp containing the vertex, its
  parent and the ancestor may be turned upside-down as the result
  of a global sequence of operations, resulting in a between or beyond
  generalization.

  As the core planarity algorithm constructs successively larger
  bicomps out of smaller ones, the bicomp root and its DFS child
  are marked as 'tied' in vertex position using markers along the
  external face.  The marking of the DFS child may be indirect.
  Since the child may not be on the external face, its descendant
  that is next along the external face is marked instead.

  Later (possibly in the same step or possibly many vertices later),
  the Walkdown proceeds around the bicomp and returns to each merge
  point, and the tie is broken based on the direction of approach.

  As the Walkdown passes a vertex to its successor, the external
  face is short-circuited to remove the vertex from it.  Immediately
  before this occurs, the new drawing method resolves the tie.  Since
  the vertex is going to the internal face, its vertex position should
  be 'between' its successor and the current vertex being processed
  by the Walkdown.

  If the vertex is a child of its external face successor, then it
  is simply marked as being 'between' that successor and the current
  vertex being processed by the planarity method.  But if the vertex
  is the parent of its external face successor, then the successor
  is placed 'beyond' the vertex.  Recall that the successor is either
  the DFS child of the vertex or a descendant of that DFS child that
  was specially marked because it, not the DFS child, was on the
  external face.

  This explains the information that has been collected by the
  planarity embedder, which will now be turned into a vertex ordering
  system.  The idea is to proceed with a pre-order traversal of
  the DFS tree, determining the relative orders of the ancestors of
  a vertex by the time we get to a vertex.  This will allow us to
  convert between/beyond into above/below based on the known relative
  order of the parent and some given ancestor of the vertex.  A vertex
  would then be added immediately above or below its parent in the
  total ordering, and then the algorithm proceeds to the descendants.

  Consider a depth-first pre-order visitation of vertices.  If the
  full order of all vertices visited so far is dynamically maintained,
  then it is easy to decide whether a vertex goes above or below
  its parent based on the between/beyond indicator and the relative
  positions in the order of the parent and given ancestor of the
  vertex.  If the ancestor is above the parent, then 'between' means
  put the vertex immediately above its parent and 'beyond' means put
  the vertex immediately below its parent in the order.  And if the
  ancestor is below the parent, then the meaning of between and
  beyond are simply reversed.

  Once a vertex is known to be above or below its parent, the drawing
  flag is changed from between/beyond to above/below, and processing
  proceeds to the next vertex in pre-order depth first search.

  The difficulty lies in keeping an up-to-date topological ordering
  that can be queried in constant time to find the relative positions
  of two vertices.  By itself, this is an instance of "online" or
  dynamic topological sorting and has been proven not to be achievable
  in linear total time.  But this is a special case of the problem and
  is therefore solvable through the collection and maintenance of some
  additional information.

  Recall that the ancestor V of a vertex is recorded when the setting
  for between/beyond is made for a vertex. However, the Walkdown is
  invoked on the bicomp rooted by edge (V', C), so the child C of V
  that roots the subtree containing the vertex being marked is known.

  Note that when a DFS child is placed above its parent, the entire
  DFS subtree of vertices is placed above the parent.  Hence, to
  determine whether the parent P of a vertex W is above the ancestor
  V, where W is marked either between or beyond P and V, we need
  only determine the relationship between V and C, which has already
  been directly determined due to previous steps of the algorithm
  (because V and C are ancestors of P and W).  If C is above/below V
  then so is P.

  As mentioned above, once the position of P is known relative to V,
  it is a simple matter to decide whether to put W above or below P
  based on the between/beyond indicator stored in W during embedding.
 ********************************************************************/

int _ComputeVertexPositionsInComponent(DrawPlanarContext *context, int root, int *pIndex)
{
graphP theEmbedding = context->theGraph;
listCollectionP theOrder = LCNew(theEmbedding->N);
int W, P, C, V, J;

    if (theOrder == NULL)
        return NOTOK;

    // Determine the vertex order using a depth first search with
    // pre-order visitation.

    sp_ClearStack(theEmbedding->theStack);
    sp_Push(theEmbedding->theStack, root);
    while (!sp_IsEmpty(theEmbedding->theStack))
    {
        sp_Pop(theEmbedding->theStack, W);

        P = theEmbedding->V[W].DFSParent;
        V = context->V[W].ancestor;
        C = context->V[W].ancestorChild;

        // For the special case that we just popped the DFS tree root,
        // we simply add the root to its own position.
        if (P == NIL)
        {
            // Put the DFS root in the list by itself
            LCAppend(theOrder, NIL, W);
            // The children of the DFS root have the root as their
            // ancestorChild and 'beyond' as the drawingFlag, so this
            // causes the root's children to be placed below the root
            context->V[W].drawingFlag = DRAWINGFLAG_BELOW;
        }

        // Determine vertex W position relative to P
        else
        {
            // An unresolved tie is an error
            if (context->V[W].drawingFlag == DRAWINGFLAG_TIE)
                return NOTOK;

            // If C below V, then P below V, so interpret W between
            // P and V as W above P, and interpret W beyond P relative
            // to V as W below P.
            if (context->V[C].drawingFlag == DRAWINGFLAG_BELOW)
            {
                if (context->V[W].drawingFlag == DRAWINGFLAG_BETWEEN)
                    context->V[W].drawingFlag = DRAWINGFLAG_ABOVE;
                else
                    context->V[W].drawingFlag = DRAWINGFLAG_BELOW;
            }

            // If C above V, then P above V, so interpret W between
            // P and V as W below P, and interpret W beyond P relative
            // to V as W above P.
            else
            {
                if (context->V[W].drawingFlag == DRAWINGFLAG_BETWEEN)
                    context->V[W].drawingFlag = DRAWINGFLAG_BELOW;
                else
                    context->V[W].drawingFlag = DRAWINGFLAG_ABOVE;
            }

            if (context->V[W].drawingFlag == DRAWINGFLAG_BELOW)
                LCInsertAfter(theOrder, P, W);
            else
                LCInsertBefore(theOrder, P, W);
        }

        // Push DFS children
        J = gp_GetFirstArc(theEmbedding, W);
        while (gp_IsArc(theEmbedding, J))
        {
            if (theEmbedding->G[J].type == EDGE_DFSCHILD)
                sp_Push(theEmbedding->theStack, theEmbedding->G[J].v);

            J = gp_GetNextArc(theEmbedding, J);
        }
    }

    // Use the order to assign vertical positions
    V = root;
    while (V != NIL)
    {
        context->G[V].pos = *pIndex;
        (*pIndex)++;
        V = LCGetNext(theOrder, root, V);
    }

    // Clean up and return

    LCFree(&theOrder);
    return OK;
}


#ifdef LOGGING
/********************************************************************
 _LogEdgeList()
 Used to show the progressive calculation of the edge position list.
 ********************************************************************/
void _LogEdgeList(graphP theEmbedding, listCollectionP edgeList, int edgeListHead)
{
    int e = edgeListHead, J, JTwin;

    gp_Log("EdgeList: [ ");

    while (e != NIL)
    {
        J = theEmbedding->edgeOffset + 2*e;
        JTwin = gp_GetTwinArc(theEmbedding, J);

        gp_Log(gp_MakeLogStr2("(%d, %d) ",
        		theEmbedding->G[theEmbedding->G[J].v].v,
        		theEmbedding->G[theEmbedding->G[JTwin].v].v));

        e = LCGetNext(edgeList, edgeListHead, e);
    }

    gp_LogLine("]");
}
#endif

/********************************************************************
 _ComputeEdgePositions()

  Performs a vertical sweep of the combinatorial planar embedding,
  developing the edge order in the horizontal sweep line as it
  advances through the vertices according to their assigned
  vertical positions.

  For expedience, the 'visited' flag for each vertex shall be used
  instead to indicate the location in the edge order list of the
  generator edge for the vertex, i.e. the first edge added to the
  vertex from a higher vertex (with lower position number).
  All edges added from this vertex to the neighbors below it are
  added immediately after the generator edge for the vertex.
 ********************************************************************/

int _ComputeEdgePositions(DrawPlanarContext *context)
{
graphP theEmbedding = context->theGraph;
int *vertexOrder = NULL;
listCollectionP edgeList = NULL;
int edgeListHead, edgeListInsertPoint;
int I, J, Jcur, e, v, vpos;
int eIndex, JTwin;

	gp_LogLine("\ngraphDrawPlanar.c/_ComputeEdgePositions() start");

    // Sort the vertices by vertical position (in linear time)

    if ((vertexOrder = (int *) malloc(theEmbedding->N * sizeof(int))) == NULL)
        return NOTOK;

    for (I = 0; I < theEmbedding->N; I++)
        vertexOrder[context->G[I].pos] = I;

    // Allocate the edge list of size M.
    //    This is an array of (prev, next) pointers.
    //    An edge at position X corresponds to the edge
    //    at position X in the graph structure, which is
    //    represented by a pair of adjacent graph nodes
    //    starting at index 2N + 2X.

    if (theEmbedding->M > 0 && (edgeList = LCNew(theEmbedding->M)) == NULL)
    {
        free(vertexOrder);
        return NOTOK;
    }

    edgeListHead = NIL;

    // Each vertex starts out with a NIL generator edge.

    for (I=0; I < theEmbedding->N; I++)
        theEmbedding->G[I].visited = NIL;

    // Perform the vertical sweep of the combinatorial embedding, using
    // the vertex ordering to guide the sweep.
    // For each vertex, each edge leading to a vertex with a higher number in
    // the vertex order is recorded as the "generator edge", or the edge of
    // first discovery of that higher numbered vertex, unless the vertex already has
    // a recorded generator edge
    for (vpos=0; vpos < theEmbedding->N; vpos++)
    {
        // Get the vertex associated with the position
        v = vertexOrder[vpos];
        gp_LogLine(gp_MakeLogStr3("Processing vertex %d with DFI=%d at position=%d",
    				 theEmbedding->G[v].v, v, vpos));

        // The DFS tree root of a connected component is always the least
        // number vertex in the vertex ordering.  We have to give it a
        // false generator edge so that it is still "visited" and then
        // all of its edges are generators for its neighbor vertices because
        // they all have greater numbers in the vertex order.
        if (theEmbedding->V[v].DFSParent == NIL)
        {
            // False generator edge, so the vertex is distinguishable from
            // a vertex with no generator edge when its neighbors are visited
            // This way, an edge from a neighbor won't get recorded as the
            // generator edge of the DFS tree root.
            theEmbedding->G[v].visited = 1;

            // Now we traverse the adjacency list of the DFS tree root and
            // record each edge as the generator edge of the neighbors
            J = gp_GetFirstArc(theEmbedding, v);
            while (gp_IsArc(theGraph, J))
            {
                e = (J - theEmbedding->edgeOffset) / 2;

                edgeListHead = LCAppend(edgeList, edgeListHead, e);
                gp_LogLine(gp_MakeLogStr2("Append generator edge (%d, %d) to edgeList",
                		theEmbedding->G[v].v, theEmbedding->G[theEmbedding->G[J].v].v));

                // Set the generator edge for the root's neighbor
                theEmbedding->G[theEmbedding->G[J].v].visited = J;

                // Go to the next node of the root's adj list
                J = gp_GetNextArc(theEmbedding, J);
            }
        }

        // Else, if we are not on a DFS tree root...
        else
        {
            // Get the generator edge of the vertex
            if ((JTwin = theEmbedding->G[v].visited) == NIL)
                return NOTOK;
            J = gp_GetTwinArc(theEmbedding, JTwin);

            // Traverse the edges of the vertex, starting
            // from the generator edge and going counterclockwise...

            e = (J - theEmbedding->edgeOffset) / 2;
            edgeListInsertPoint = e;

            Jcur = gp_GetNextArcCircular(theEmbedding, J);

            while (Jcur != J)
            {
                // If the neighboring vertex's position is greater
                // than the current vertex (meaning it is lower in the
                // diagram), then add that edge to the edge order.

                if (context->G[theEmbedding->G[Jcur].v].pos > vpos)
                {
                    e = (Jcur - theEmbedding->edgeOffset) / 2;
                    LCInsertAfter(edgeList, edgeListInsertPoint, e);

                    gp_LogLine(gp_MakeLogStr4("Insert (%d, %d) after (%d, %d)",
                    		theEmbedding->G[v].v,
                    		theEmbedding->G[theEmbedding->G[Jcur].v].v,
                    		theEmbedding->G[theEmbedding->G[gp_GetTwinArc(theEmbedding, J)].v].v,
                    		theEmbedding->G[theEmbedding->G[J].v].v));

                    edgeListInsertPoint = e;

                    // If the vertex does not yet have a generator edge, then set it.
                    if (theEmbedding->G[theEmbedding->G[Jcur].v].visited == NIL)
                    {
                        theEmbedding->G[theEmbedding->G[Jcur].v].visited = Jcur;
                        gp_LogLine(gp_MakeLogStr2("Generator edge (%d, %d)",
                        		theEmbedding->G[theEmbedding->G[gp_GetTwinArc(theEmbedding, J)].v].v,
                        		theEmbedding->G[theEmbedding->G[Jcur].v].v));
                    }
                }

                // Go to the next node in v's adjacency list
                Jcur = gp_GetNextArcCircular(theEmbedding, Jcur);
            }
        }

#ifdef LOGGING
        _LogEdgeList(theEmbedding, edgeList, edgeListHead);
#endif
    }

    // Now iterate through the edgeList and assign positions to the edges.
    eIndex = 0;
    e = edgeListHead;
    while (e != NIL)
    {
        J = theEmbedding->edgeOffset + 2*e;
        JTwin = gp_GetTwinArc(theEmbedding, J);

        context->G[J].pos = context->G[JTwin].pos = eIndex;

        eIndex++;

        e = LCGetNext(edgeList, edgeListHead, e);
    }

    // Clean up and return
    LCFree(&edgeList);
    free(vertexOrder);

	gp_LogLine("graphDrawPlanar.c/_ComputeEdgePositions() end\n");

    return OK;
}

/********************************************************************
 _ComputeVertexRanges()

   Assumes edge positions are known (see _ComputeEdgePositions()).
   A vertex spans horizontally the positions of the edges incident
   to it.
 ********************************************************************/

int _ComputeVertexRanges(DrawPlanarContext *context)
{
graphP theEmbedding = context->theGraph;
int I, J, min, max;

    for (I = 0; I < theEmbedding->N; I++)
    {
        min = theEmbedding->M + 1;
        max = -1;

        // Iterate the edges, except in the isolated vertex case we just
        // set the min and max to 1 since there no edges controlling where
        // it gets drawn.
        J = gp_GetFirstArc(theEmbedding, I);
        if (!gp_IsArc(theEmbedding, J))
        {
        	min = max = 0;
        }
        else
        {
            while (gp_IsArc(theEmbedding, J))
            {
                if (min > context->G[J].pos)
                    min = context->G[J].pos;

                if (max < context->G[J].pos)
                    max = context->G[J].pos;

                J = gp_GetNextArc(theEmbedding, J);
            }
        }

        context->G[I].start = min;
        context->G[I].end = max;
    }

    return OK;
}

/********************************************************************
 _ComputeEdgeRanges()

    Assumes vertex positions are known (see _ComputeVertexPositions()).
    An edges spans the vertical range of its endpoints.
 ********************************************************************/

int _ComputeEdgeRanges(DrawPlanarContext *context)
{
graphP theEmbedding = context->theGraph;
int e, J, JTwin, v1, v2, pos1, pos2;

    for (e = 0; e < theEmbedding->M; e++)
    {
        J = theEmbedding->edgeOffset + 2*e;
        JTwin = gp_GetTwinArc(theEmbedding, J);

        v1 = theEmbedding->G[J].v;
        v2 = theEmbedding->G[JTwin].v;

        pos1 = context->G[v1].pos;
        pos2 = context->G[v2].pos;

        if (pos1 < pos2)
        {
            context->G[J].start = pos1;
            context->G[J].end = pos2;
        }
        else
        {
            context->G[J].start = pos2;
            context->G[J].end = pos1;
        }

        context->G[JTwin].start = context->G[J].start;
        context->G[JTwin].end = context->G[J].end;
    }

    return OK;
}

/********************************************************************
 _GetNextExternalFaceVertex()
 Uses the extFace links to traverse to the next vertex on the external
 face given a current vertex and the link that points to its predecessor.
 ********************************************************************/
int _GetNextExternalFaceVertex(graphP theGraph, int curVertex, int *pPrevLink)
{
    int nextVertex, nextPrevLink;

    nextVertex = theGraph->extFace[curVertex].vertex[1 ^ *pPrevLink];

    // If the two links in the new vertex are not equal, then only one points
    // back to the current vertex, and it is the new prev link.
    if (theGraph->extFace[nextVertex].vertex[0] != theGraph->extFace[nextVertex].vertex[1])
    {
        nextPrevLink = theGraph->extFace[nextVertex].vertex[0]==curVertex ? 0 : 1;
    }
    else
    {
        int inverted = 0;

        // One of the two vertices is the root of the bicomp.  The non-root may
        // not have the same orientation as the root because a consistent orientation
        // is not imposed until post-processing of the planarity test.  But the
        // planarity test does special-case tracking of orientation inversions of
        // singleton bicomps in the inversionFlag of the non-root vertex.

        if (nextVertex < theGraph->N)
             inverted = theGraph->extFace[nextVertex].inversionFlag;
        else inverted = theGraph->extFace[curVertex].inversionFlag;

        // If the curVertex and nextVertex are in a singleton and have the same
        // orientation, then the prev link from the current vertex is the prev
        // link from the next vertex because you'd just be going around in a
        // small circle using the same links for prev every time.
        // If their orientations are inverted, then we just invert the prev link.

        nextPrevLink = inverted ? (1^*pPrevLink) : (*pPrevLink);
    }

    // Return the information
    *pPrevLink = nextPrevLink;
    return nextVertex;
}

/********************************************************************
 _CollectDrawingData()
 To be called by core planarity Walkdown immediately before merging
 bicomps and embedding a new back edge.

 Each bicomp is rooted by a DFS tree edge.  The parent vertex in
 that edge is the bicomp root, and the bicomp contains one DFS child
 of the vertex, which is on the child end of the 'root edge'.

 Here we decide whether the DFS child is to be embedded between or
 beyond its parent relative to vertex I, the one currently being
 processed (and the ancestor endpoint of a back edge being embedded,
 where the descendant endpoint is also an endpoint of the bicomp
 root being merged).
 ********************************************************************/

void _CollectDrawingData(DrawPlanarContext *context, int RootVertex, int W, int WPrevLink)
{
graphP theEmbedding = context->theGraph;
//int ancestorChild = RootVertex - theEmbedding->N;
//int ancestor = theEmbedding->V[ancestorChild].DFSParent;
int K, Parent, BicompRoot, DFSChild, direction, descendant;

    gp_LogLine("\ngraphDrawPlanar.c/_CollectDrawingData() start");
    gp_LogLine(gp_MakeLogStr3("_CollectDrawingData(RootVertex=%d, W=%d, W_in=%d)",
				 RootVertex, W, WPrevLink));

    /* Process all of the merge points to set their drawing flags. */

    for (K = 0; K < sp_GetCurrentSize(theEmbedding->theStack); K += 4)
    {
         /* Get the parent and child that are about to be merged from
            the 4-tuple in the merge stack */
         Parent = theEmbedding->theStack->S[K];
         BicompRoot = theEmbedding->theStack->S[K+2];
         DFSChild = BicompRoot - theEmbedding->N;

         /* We get the active descendant vertex in the child bicomp that
            will be adjacent to the parent along the external face.
            This vertex is guaranteed to be found in one step
            due to external face 'short-circuiting' that was done in
            step 'Parent' of the planarity algorithm.
            We pass theEmbedding->N for the second parameter because
            of this; we use this function to signify need of extFace
            links in the other implementation.*/

         direction = theEmbedding->theStack->S[K+3];
         descendant = _GetNextExternalFaceVertex(theEmbedding, BicompRoot, &direction);

         /* Now we set the tie flag in the DFS child, and mark the
            descendant and parent with non-NIL pointers to the child
            whose tie flag is to be resolved as soon as one of the
            two is connected to by an edge or child bicomp merge. */

         context->V[DFSChild].drawingFlag = DRAWINGFLAG_TIE;

         context->V[descendant].tie[direction] = DFSChild;

         direction = theEmbedding->theStack->S[K+1];
         context->V[Parent].tie[direction] = DFSChild;

         gp_LogLine(gp_MakeLogStr5("V[Parent=%d]=.tie[%d] = V[descendant=%d].tie[%d] = (child=%d)",
					 Parent, direction, descendant, theEmbedding->theStack->S[K+3], DFSChild));
    }

    gp_LogLine("graphDrawPlanar.c/_CollectDrawingData() end\n");
}

/********************************************************************
 _BreakTie()

 The given vertex W has just been arrived at by the core planarity
 algorithm.  Using WPrevLink, we seek its predecessor WPred on the
 external face and test whether the two are involved in a tie that
 can be resolved.

 Since the planarity algorithm has just passed by WPred, it is
 safe to conclude that WPred can go between W and the current vertex.

 Of course, if W was the parent to some DFS child whose subtree
 contains WPred, then the DFS child is marked 'between', placing
 the whole subtree including WPred between W and the current vertex.
 On the other hand, if WPred was the parent of some DFS child whose
 subtree contained W, then we achieve the same effect of putting WPred
 'between' W and the curent vertex by marking the DFS child 'beyond'.
 Since the DFS child and hence W are beyond W relative to the current
 vertex, WPred is also between W and the current vertex.

 Thus the certain positional relationship between W and WPred
 relative to a specific ancestor, the current vertex, is used to
 indirectly break the positional tie between MIN(W, WPred) and the
 DFS child of MIN(W, WPred) whose subtree contains MAX(W, WPred).

 The ancestorChild is the DFS child of the current vertex whose DFS
 subtree contains W and WPred, and it is recorded here in order to
 optimize the post-processing calculation of vertex positions.
 ********************************************************************/

int _BreakTie(DrawPlanarContext *context, int BicompRoot, int W, int WPrevLink)
{
graphP theEmbedding = context->theGraph;

    /* First we get the predecessor of W. */

int WPredNextLink = 1^WPrevLink,
    WPred = _GetNextExternalFaceVertex(theEmbedding, W, &WPredNextLink);

	gp_LogLine("\ngraphDrawPlanar.c/::_BreakTie() start");
    gp_LogLine(gp_MakeLogStr3("_BreakTie(BicompRoot=%d, W=%d, W_in=%d)",
				 BicompRoot, W, WPrevLink));

    /* Ties happen only within a bicomp (i.e. between two non-root vertices) */
    if (W > theEmbedding->N || WPred >= theEmbedding->N)
    {
    	gp_LogLine("graphDrawPlanar.c/_BreakTie() end\n");
        return OK;
    }

    /* The two vertices are either tied or not; having one tied and the other
        not is an error */

    if (context->V[W].tie[WPrevLink] != context->V[WPred].tie[WPredNextLink])
        return NOTOK;

    /* If there is a tie, it can now be resolved. */
    if (context->V[W].tie[WPrevLink] != NIL)
    {
        int DFSChild = context->V[W].tie[WPrevLink];

        /* Set the two ancestor variables that contextualize putting W 'between'
            or 'beyond' its parent relative to what. */

        context->V[DFSChild].ancestorChild = BicompRoot - theEmbedding->N;
        context->V[DFSChild].ancestor = theEmbedding->V[BicompRoot - theEmbedding->N].DFSParent;

        gp_LogLine(gp_MakeLogStr4("V[child=%d]=.ancestorChild = %d, V[child=%d]=.ancestor = %d",
					 DFSChild, context->V[DFSChild].ancestorChild, DFSChild, context->V[DFSChild].ancestor));

        /* If W is the ancestor of WPred, then the DFSChild subtree contains
            WPred, and so must go between W and some ancestor. */
        if (W < WPred)
        {
            context->V[DFSChild].drawingFlag = DRAWINGFLAG_BETWEEN;
            gp_LogLine(gp_MakeLogStr3("Child=%d is 'between' ancestorChild=%d and ancestor=%d",
    					 DFSChild, context->V[DFSChild].ancestorChild, context->V[DFSChild].ancestor));
        }

        /* If W is the descendant, so we achieve the effect of putting WPred
           between DFSChild and ancestor by putting the DFSChild 'beyond' WPred. */
        else
        {
            context->V[DFSChild].drawingFlag = DRAWINGFLAG_BEYOND;
            gp_LogLine(gp_MakeLogStr3("Child=%d is 'beyond' ancestorChild=%d relative to ancestor=%d",
    					 DFSChild, context->V[DFSChild].ancestorChild, context->V[DFSChild].ancestor));
        }

        /* The tie is resolved so clear the flags*/
        context->V[W].tie[WPrevLink] = context->V[WPred].tie[WPredNextLink] = NIL;
    }
    else
        return NOTOK;

	gp_LogLine("graphDrawPlanar.c/_BreakTie() end\n");
    return OK;
}

/********************************************************************
 _RenderToString()
 Draws the previously calculated visibility representation in a
 string of size (M+1)*2N + 1 characters, which should be deallocated
 with free().

 Returns NULL on failure, or the string containing the visibility
 representation otherwise.  The string can be printed using %s,
 ********************************************************************/

char *_RenderToString(graphP theEmbedding)
{
    DrawPlanarContext *context = NULL;
    gp_FindExtension(theEmbedding, DRAWPLANAR_ID, (void *) &context);

    if (context != NULL)
    {
        int N = theEmbedding->N;
        int M = theEmbedding->M;
        int I, J, e, Mid, Pos;
        char *visRep = (char *) malloc(sizeof(char) * ((M+1) * 2*N + 1));
        char numBuffer[32];

        if (visRep == NULL)
            return NULL;

        if (sp_NonEmpty(context->theGraph->edgeHoles))
        {
            free(visRep);
            return NULL;
        }

        // Clear the space
        for (I = 0; I < N; I++)
        {
            for (J=0; J < M; J++)
            {
                visRep[(2*I) * (M+1) + J] = ' ';
                visRep[(2*I+1) * (M+1) + J] = ' ';
            }

            visRep[(2*I) * (M+1) + M] = '\n';
            visRep[(2*I+1) * (M+1) + M] = '\n';
        }

        // Draw the vertices
        for (I = 0; I < N; I++)
        {
            Pos = context->G[I].pos;
            for (J=context->G[I].start; J<=context->G[I].end; J++)
                visRep[(2*Pos) * (M+1) + J] = '-';

            // Draw vertex label
            Mid = (context->G[I].start + context->G[I].end)/2;
            sprintf(numBuffer, "%d", I);
            if ((unsigned)(context->G[I].end - context->G[I].start + 1) >= strlen(numBuffer))
            {
                strncpy(visRep + (2*Pos) * (M+1) + Mid, numBuffer, strlen(numBuffer));
            }
            // If the vertex width is less than the label width, then fail gracefully
            else
            {
                if (strlen(numBuffer)==2)
                    visRep[(2*Pos) * (M+1) + Mid] = numBuffer[0];
                else
                    visRep[(2*Pos) * (M+1) + Mid] = '*';

                visRep[(2*Pos+1) * (M+1) + Mid] = numBuffer[strlen(numBuffer)-1];
            }
        }

        // Draw the edges
        for (e=0; e<M; e++)
        {
            J = 2*N + 2*e;
            Pos = context->G[J].pos;
            for (I=context->G[J].start; I<context->G[J].end; I++)
            {
                if (I > context->G[J].start)
                    visRep[(2*I) * (M+1) + Pos] = '|';
                visRep[(2*I+1) * (M+1) + Pos] = '|';
            }
        }

        // Null terminate string and return it
        visRep[(M+1) * 2*N] = '\0';
        return visRep;
    }

    return NULL;
}

/********************************************************************
 gp_DrawPlanar_RenderToFile()
 Creates a rendition of the planar graph visibility representation
 as a string, then dumps the string to the file.
 ********************************************************************/
int gp_DrawPlanar_RenderToFile(graphP theEmbedding, char *theFileName)
{
    if (sp_IsEmpty(theEmbedding->edgeHoles))
    {
        FILE *outfile;
        char *theRendition;

        if (strcmp(theFileName, "stdout") == 0)
             outfile = stdout;
        else if (strcmp(theFileName, "stderr") == 0)
             outfile = stderr;
        else outfile = fopen(theFileName, WRITETEXT);

        if (outfile == NULL)
            return NOTOK;

        theRendition = _RenderToString(theEmbedding);
        if (theRendition != NULL)
        {
            fprintf(outfile, "%s", theRendition);
            free(theRendition);
        }

        if (strcmp(theFileName, "stdout") == 0 || strcmp(theFileName, "stderr") == 0)
            fflush(outfile);

        else if (fclose(outfile) != 0)
            return NOTOK;

        return theRendition ? OK : NOTOK;
    }

    return NOTOK;
}

/********************************************************************
 _CheckVisibilityRepresentationIntegrity()
 ********************************************************************/

int _CheckVisibilityRepresentationIntegrity(DrawPlanarContext *context)
{
graphP theEmbedding = context->theGraph;
int I, e, J, JTwin, JPos, JIndex;

    if (sp_NonEmpty(context->theGraph->edgeHoles))
        return NOTOK;

    _FillVisitedFlags(theEmbedding, 0);

/* Test whether the vertex values make sense and
        whether the vertex positions are unique. */

    for (I = 0; I < theEmbedding->N; I++)
    {
    	if (theEmbedding->M > 0)
    	{
            if (context->G[I].pos < 0 ||
                context->G[I].pos >= theEmbedding->N ||
                context->G[I].start < 0 ||
                context->G[I].start > context->G[I].end ||
                context->G[I].end >= theEmbedding->M)
                return NOTOK;
    	}

        // Has the vertex position (context->G[I].pos) been used by a
        // vertex before vertex I?
        if (theEmbedding->G[context->G[I].pos].visited)
            return NOTOK;

        // Mark the vertex position as used by vertex I.
        // Note that this marking is made on some other vertex unrelated to I
        // We're just reusing the vertex visited array as cheap storage for a
        // detector of reusing vertex position integers.
        theEmbedding->G[context->G[I].pos].visited = 1;
    }

/* Test whether the edge values make sense and
        whether the edge positions are unique */

    for (e = 0; e < theEmbedding->M; e++)
    {
        /* Each edge has an index location J in the graph structure */
        J = theEmbedding->edgeOffset + 2*e;
        JTwin = gp_GetTwinArc(theEmbedding, J);

        if (context->G[J].pos != context->G[JTwin].pos ||
            context->G[J].start != context->G[JTwin].start ||
            context->G[J].end != context->G[JTwin].end ||
            context->G[J].pos < 0 ||
            context->G[J].pos >= theEmbedding->M ||
            context->G[J].start < 0 ||
            context->G[J].start > context->G[J].end ||
            context->G[J].end >= theEmbedding->N)
            return NOTOK;

        /* Get the recorded horizontal position of that edge,
            a number between 0 and M-1 */

        JPos = context->G[J].pos;

        /* Convert that to an index in the graph structure so we
            can use the visited flags in the graph's edges to
            tell us whether the positions are being reused. */

        JIndex = theEmbedding->edgeOffset + 2*JPos;
        JTwin = gp_GetTwinArc(theEmbedding, JIndex);

        if (theEmbedding->G[JIndex].visited || theEmbedding->G[JTwin].visited)
            return NOTOK;

        theEmbedding->G[JIndex].visited = theEmbedding->G[JTwin].visited = 1;
    }

/* Test whether any edge intersects any vertex position
    for a vertex that is not an endpoint of the edge. */

    for (e = 0; e < theEmbedding->M; e++)
    {
        J = theEmbedding->edgeOffset + 2*e;
        JTwin = gp_GetTwinArc(theEmbedding, J);

        for (I = 0; I < theEmbedding->N; I++)
        {
            /* If the vertex is an endpoint of the edge, then... */

            if (theEmbedding->G[J].v == I || theEmbedding->G[JTwin].v == I)
            {
                /* The vertical position of the vertex must be
                   at the top or bottom of the edge,  */
                if (context->G[J].start != context->G[I].pos &&
                    context->G[J].end != context->G[I].pos)
                    return NOTOK;

                /* The horizontal edge position must be in the range of the vertex */
                if (context->G[J].pos < context->G[I].start ||
                    context->G[J].pos > context->G[I].end)
                    return NOTOK;
            }

            /* If the vertex is not an endpoint of the edge... */

            else // if (theEmbedding->G[J].v != I && theEmbedding->G[JTwin].v != I)
            {
                /* If the vertical position of the vertex is in the
                    vertical range of the edge ... */

                if (context->G[J].start <= context->G[I].pos &&
                    context->G[J].end >= context->G[I].pos)
                {
                    /* And if the horizontal position of the edge is in the
                        horizontal range of the vertex, then return an error. */

                    if (context->G[I].start <= context->G[J].pos &&
                        context->G[I].end >= context->G[J].pos)
                        return NOTOK;
                }
            }
        }
    }


/* All tests passed */

    return OK;
}
