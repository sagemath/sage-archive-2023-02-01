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

#define GRAPH_C

#include <stdlib.h>

#include "graph.h"

/* Imported functions */

extern void _FillVisitedFlags(graphP, int);

extern int _IsolateKuratowskiSubgraph(graphP theGraph, int I, int R);
extern int _IsolateOuterplanarObstruction(graphP theGraph, int I, int R);

/* Private functions (some are exported to system only) */

void _CreateSortedSeparatedDFSChildLists(graphP theGraph);
int  _CreateFwdArcLists(graphP theGraph);
void _CreateDFSTreeEmbedding(graphP theGraph);

void _EmbedBackEdgeToDescendant(graphP theGraph, int RootSide, int RootVertex, int W, int WPrevLink);

int  _GetNextVertexOnExternalFace(graphP theGraph, int curVertex, int *pPrevLink);

void _InvertVertex(graphP theGraph, int V);
void _MergeVertex(graphP theGraph, int W, int WPrevLink, int R);
int  _MergeBicomps(graphP theGraph, int I, int RootVertex, int W, int WPrevLink);

void _WalkUp(graphP theGraph, int I, int J);
int  _WalkDown(graphP theGraph, int I, int RootVertex);

int  _HandleBlockedEmbedIteration(graphP theGraph, int I);
int  _EmbedPostprocess(graphP theGraph, int I, int edgeEmbeddingResult);

int  _OrientVerticesInEmbedding(graphP theGraph);
int  _OrientVerticesInBicomp(graphP theGraph, int BicompRoot, int PreserveSigns);
int  _JoinBicomps(graphP theGraph);

/********************************************************************
 _CreateSortedSeparatedDFSChildLists()
 We create a separatedDFSChildList in each vertex which contains the
 Lowpoint values of the vertex's DFS children sorted in non-descending order.
 To accomplish this in linear time for the whole graph, we must not
 sort the DFS children in each vertex, but rather bucket sort the
 Lowpoint values of all vertices, then traverse the buckets sequentially,
 adding each vertex to its parent's separatedDFSChildList.
 Note that this is a specialized bucket sort that achieves O(n)
 worst case rather than O(n) expected time due to the simplicity
 of the sorting problem.  Specifically, we know that the Lowpoint values
 are between 0 and N-1, so we create buckets for each value.
 Collisions occur only when two keys are equal, so there is no
 need to sort the buckets (hence O(n) worst case).
 ********************************************************************/

void  _CreateSortedSeparatedDFSChildLists(graphP theGraph)
{
int *buckets;
listCollectionP bin;
int I, J, N, DFSParent, theList;

     N = theGraph->N;
     buckets = theGraph->buckets;
     bin = theGraph->bin;

     /* Initialize the bin and all the buckets to be empty */

     LCReset(bin);
     for (I=0; I < N; I++)
          buckets[I] = NIL;

     /* For each vertex, add it to the bucket whose index is equal to
        the Lowpoint of the vertex. */

     for (I=0; I < N; I++)
     {
          J = theGraph->V[I].Lowpoint;
          buckets[J] = LCAppend(bin, buckets[J], I);
     }

     /* For each bucket, add each vertex in the bucket to the
        separatedDFSChildList of its DFSParent.  Since lower numbered buckets
        are processed before higher numbered buckets, vertices with lower
        Lowpoint values are added before those with higher Lowpoint values,
        so the separatedDFSChildList of each vertex is sorted by Lowpoint */

     for (I = 0; I < N; I++)
     {
          if ((J=buckets[I]) != NIL)
          {
              while (J != NIL)
              {
                  DFSParent = theGraph->V[J].DFSParent;

                  if (DFSParent != NIL && DFSParent != J)
                  {
                      theList = theGraph->V[DFSParent].separatedDFSChildList;
                      theList = LCAppend(theGraph->DFSChildLists, theList, J);
                      theGraph->V[DFSParent].separatedDFSChildList = theList;
                  }

                  J = LCGetNext(bin, buckets[I], J);
              }
          }
     }
}

/********************************************************************
 _CreateFwdArcLists()

 Puts the forward arcs (back edges from a vertex to its descendants)
 into a circular list indicated by the fwdArcList member, a task
 simplified by the fact that they have already been placed in
 succession at the end of the adjacency lists by gp_CreateDFSTree().

  Returns OK for success, NOTOK for internal code failure
 ********************************************************************/

int _CreateFwdArcLists(graphP theGraph)
{
int I, Jfirst, Jnext, Jlast;

    for (I=0; I < theGraph->N; I++)
    {
    	// The forward arcs are already in succession at the end of the adjacency list
    	// Skip this vertex if it has no edges

    	Jfirst = gp_GetLastArc(theGraph, I);
    	if (!gp_IsArc(theGraph, Jfirst))
    		continue;

        // If the vertex has any forward edges at all, then the last edge
    	// will be a forward edge.  So if we have any forward edges, ...

        if (theGraph->G[Jfirst].type == EDGE_FORWARD)
        {
            // Find the end of the forward edge list

            Jnext = Jfirst;
            while (theGraph->G[Jnext].type == EDGE_FORWARD)
                Jnext = gp_GetPrevArc(theGraph, Jnext);
            Jlast = gp_GetNextArc(theGraph, Jnext);

            // Remove the forward edges from the adjacency list of I
            gp_BindLastArc(theGraph, I, Jnext);

            // Make a circular forward edge list
            theGraph->V[I].fwdArcList = Jfirst;
            gp_SetNextArc(theGraph, Jfirst, Jlast);
            gp_SetPrevArc(theGraph, Jlast, Jfirst);
        }
    }

    return OK;
}

/********************************************************************
 ********************************************************************/

#ifdef DEBUG
int  TestIntegrity(graphP theGraph)
{
    int I, Jcur, result = 1;

        for (I=0; I < theGraph->N; I++)
        {
            Jcur = theGraph->V[I].fwdArcList;
            while (Jcur != NIL)
            {
                if (theGraph->G[Jcur].visited)
                {
                    printf("Found problem with fwdArcList of vertex %d.\n", I);
                    result = 0;
                    break;
                }

                theGraph->G[Jcur].visited = 1;

                Jcur = gp_GetNextArc(theGraph, Jcur);
                if (Jcur == theGraph->V[I].fwdArcList)
                    Jcur = NIL;
            }

            Jcur = theGraph->V[I].fwdArcList;
            while (Jcur != NIL)
            {
                if (!theGraph->G[Jcur].visited)
                    break;

                theGraph->G[Jcur].visited = 0;

                Jcur = gp_GetNextArc(theGraph, Jcur);
                if (Jcur == theGraph->V[I].fwdArcList)
                    Jcur = NIL;
            }
        }

        return result;
}
#endif

/********************************************************************
 _CreateDFSTreeEmbedding()

 Each vertex receives only its parent arc in the adjacency list, and
 the corresponding child arc is placed in the adjacency list of a root
 copy of the parent.  Each root copy of a vertex is uniquely associated
 with a child C, so it is simply stored at location C+N.

 The forward arcs are not lost because they are already in the
 fwdArcList of each vertex.  Each back arc can be reached as the
 twin arc of a forward arc, and the two are embedded together when
 the forward arc is processed.  Finally, the child arcs are initially
 placed in root copies of vertices, not the vertices themselves, but
 the child arcs are merged into the vertices as the embedder progresses.
 ********************************************************************/

void _CreateDFSTreeEmbedding(graphP theGraph)
{
int N, I, J, Jtwin, R;

    N = theGraph->N;

    // Embed all tree edges.  For each DFS tree child, we move
    // the child arc to a root copy of vertex I that is uniquely
    // associated with the DFS child, and we remove all edges
    // from the child except the parent arc

    for (I=0, R=N; I < N; I++, R++)
    {
        if (theGraph->V[I].DFSParent == NIL)
        {
        	gp_SetFirstArc(theGraph, I, gp_AdjacencyListEndMark(I));
        	gp_SetLastArc(theGraph, I, gp_AdjacencyListEndMark(I));
        }
        else
        {
            J = gp_GetFirstArc(theGraph, I);
            while (theGraph->G[J].type != EDGE_DFSPARENT)
                J = gp_GetNextArc(theGraph, J);

        	gp_SetFirstArc(theGraph, I, J);
        	gp_SetLastArc(theGraph, I, J);

        	gp_SetNextArc(theGraph, J, gp_AdjacencyListEndMark(I));
        	gp_SetPrevArc(theGraph, J, gp_AdjacencyListEndMark(I));

        	theGraph->G[J].v = R;

            Jtwin = gp_GetTwinArc(theGraph, J);

        	gp_SetFirstArc(theGraph, R, Jtwin);
        	gp_SetLastArc(theGraph, R, Jtwin);

        	gp_SetNextArc(theGraph, Jtwin, gp_AdjacencyListEndMark(R));
        	gp_SetPrevArc(theGraph, Jtwin, gp_AdjacencyListEndMark(R));

            theGraph->extFace[R].vertex[0] = theGraph->extFace[R].vertex[1] = I;
            theGraph->extFace[I].vertex[0] = theGraph->extFace[I].vertex[1] = R;
        }
    }
}

/********************************************************************
 _EmbedBackEdgeToDescendant()
 The Walkdown has found a descendant vertex W to which it can
 attach a back edge up to the root of the bicomp it is processing.
 The RootSide and WPrevLink indicate the parts of the external face
 that will be replaced at each endpoint of the back edge.
 ********************************************************************/

void _EmbedBackEdgeToDescendant(graphP theGraph, int RootSide, int RootVertex, int W, int WPrevLink)
{
int fwdArc, backArc, parentCopy;

    /* We get the two edge records of the back edge to embed.
        The Walkup recorded in W's adjacentTo the index of the forward arc
        from the root's parent copy to the descendant W. */

    fwdArc = theGraph->V[W].adjacentTo;
    backArc = gp_GetTwinArc(theGraph, fwdArc);

    /* The forward arc is removed from the fwdArcList of the root's parent copy. */

    parentCopy = theGraph->V[RootVertex - theGraph->N].DFSParent;

    gp_LogLine(gp_MakeLogStr5("graphEmbed.c/_EmbedBackEdgeToDescendant() V=%d, R=%d, R_out=%d, W=%d, W_in=%d",
    		parentCopy, RootVertex, RootSide, W, WPrevLink));

    if (theGraph->V[parentCopy].fwdArcList == fwdArc)
    {
    	theGraph->V[parentCopy].fwdArcList = gp_GetNextArc(theGraph, fwdArc);
        if (theGraph->V[parentCopy].fwdArcList == fwdArc)
            theGraph->V[parentCopy].fwdArcList = NIL;
    }

    gp_SetNextArc(theGraph, gp_GetPrevArc(theGraph, fwdArc), gp_GetNextArc(theGraph, fwdArc));
    gp_SetPrevArc(theGraph, gp_GetNextArc(theGraph, fwdArc), gp_GetPrevArc(theGraph, fwdArc));

    // The forward arc is added to the adjacency list of the RootVertex.
    // Note that we're guaranteed that the RootVertex adjacency list is non-empty,
    // so tests for NIL are not needed
    gp_SetAdjacentArc(theGraph, fwdArc, 1^RootSide, gp_AdjacencyListEndMark(RootVertex));
    gp_SetAdjacentArc(theGraph, fwdArc, RootSide, gp_GetArc(theGraph, RootVertex, RootSide));
    gp_SetAdjacentArc(theGraph, gp_GetArc(theGraph, RootVertex, RootSide), 1^RootSide, fwdArc);
    gp_SetArc(theGraph, RootVertex, RootSide, fwdArc);

    // The back arc is added to the adjacency list of W.
    // The adjacency list of W is also guaranteed non-empty
    gp_SetAdjacentArc(theGraph, backArc, 1^WPrevLink, gp_AdjacencyListEndMark(W));
    gp_SetAdjacentArc(theGraph, backArc, WPrevLink, gp_GetArc(theGraph, W, WPrevLink));
    gp_SetAdjacentArc(theGraph, gp_GetArc(theGraph, W, WPrevLink), 1^WPrevLink, backArc);
    gp_SetArc(theGraph, W, WPrevLink, backArc);

    theGraph->G[backArc].v = RootVertex;

    /* Link the two endpoint vertices together on the external face */

    theGraph->extFace[RootVertex].vertex[RootSide] = W;
    theGraph->extFace[W].vertex[WPrevLink] = RootVertex;
}

/********************************************************************
 _GetNextVertexOnExternalFace()
 Each vertex contains two 'link' index pointers that indicate the
 first and last adjacency list arc.  If the vertex is on the external face,
 then these two arcs are also on the external face.  We want to take one of
 those edges to get to the next vertex on the external face.
 On input *pPrevLink indicates which link we followed to arrive at
 curVertex.  On output *pPrevLink will be set to the link we follow to
 get into the next vertex.
 To get to the next vertex, we use the opposite link from the one used
 to get into curVertex.  This takes us to an edge node.  The twinArc
 of that edge node, carries us to an edge node in the next vertex.
 At least one of the two links in that edge node will lead to a vertex
 node in G, which is the next vertex.  Once we arrive at the next
 vertex, at least one of its links will lead back to the edge node, and
 that link becomes the output value of *pPrevLink.

 NOTE: This method intentionally ignores the extFace optimization
       links. It is invoked when the "real" external face must be
       traversed and hence when the constant time guarantee is not
       needed from the extFace short-circuit that connects the
       bicomp root to the first active vertices along each external
       face path emanating from the bicomp root.
 ********************************************************************/

int  _GetNextVertexOnExternalFace(graphP theGraph, int curVertex, int *pPrevLink)
{
     /* Exit curVertex from whichever link was not previously used to enter it */

     int arc = gp_GetArc(theGraph, curVertex, 1^(*pPrevLink));
     int nextVertex = theGraph->G[arc].v;

     /* This if stmt assigns the new prev link that tells us which edge
        record was used to enter nextVertex (so that we exit from the
        opposing edge record).

        However, if we are in a singleton bicomp, then both links in nextVertex
        lead back to curVertex.  We want the two arcs of a singleton bicomp to
        act like a cycle, so we just don't change the prev link in this case.

        But when nextVertex has more than one edge, we need to figure out
        whether the first edge or last edge (which are the two on the external
        face) was used to enter nextVertex so we can exit from the other one
        as traversal of the external face continues later. */

     if (gp_GetFirstArc(theGraph, nextVertex) != gp_GetLastArc(theGraph, nextVertex))
         *pPrevLink = gp_GetTwinArc(theGraph, arc) == gp_GetFirstArc(theGraph, nextVertex) ? 0 : 1;

     return nextVertex;
}

/********************************************************************
 _InvertVertex()
 This function flips the orientation of a single vertex such that
 instead of using link successors to go clockwise (or counterclockwise)
 around a vertex's adjacency list, link predecessors would be used.
 ********************************************************************/

void _InvertVertex(graphP theGraph, int V)
{
int J, temp;

	 gp_LogLine(gp_MakeLogStr1("graphEmbed.c/_InvertVertex() V=%d", V));

     // Swap the links in all the arcs of the adjacency list
     J = gp_GetFirstArc(theGraph, V);
     while (gp_IsArc(theGraph, J))
     {
    	 temp = gp_GetNextArc(theGraph, J);
    	 gp_SetNextArc(theGraph, J, gp_GetPrevArc(theGraph, J));
    	 gp_SetPrevArc(theGraph, J, temp);

         J = temp;
     }

     // Swap the first/last edge record indicators in the vertex
     temp = gp_GetFirstArc(theGraph, V);
     gp_SetFirstArc(theGraph, V, gp_GetLastArc(theGraph, V));
     gp_SetLastArc(theGraph, V, temp);

     // Swap the first/last external face indicators in the vertex
     temp = theGraph->extFace[V].vertex[0];
     theGraph->extFace[V].vertex[0] = theGraph->extFace[V].vertex[1];
     theGraph->extFace[V].vertex[1] = temp;
}

/********************************************************************
 _MergeVertex()
 The merge step joins the vertex W to the root R of a child bicompRoot,
 which is a root copy of W appearing in the region N to 2N-1.

 Actually, the first step of this is to redirect all of the edges leading
 into R so that they indicate W as the neighbor instead of R.
 For each edge node pointing to R, we set the 'v' field to W.  Once an
 edge is redirected from a root copy R to a parent copy W, the edge is
 never redirected again, so we associate the cost of the redirection
 as constant per edge, which maintains linear time performance.

 After this is done, a regular circular list union occurs. The only
 consideration is that WPrevLink is used to indicate the two edge
 records e_w and e_r that will become consecutive in the resulting
 adjacency list of W.  We set e_w to W's link [WPrevLink] and e_r to
 R's link [1^WPrevLink] so that e_w and e_r indicate W and R with
 opposing links, which become free to be cross-linked.  Finally,
 the edge record e_ext, set equal to R's link [WPrevLink], is the edge
 that, with e_r, held R to the external face.  Now, e_ext will be the
 new link [WPrevLink] edge record for W.  If e_w and e_r become part
 of a proper face, then e_ext and W's link [1^WPrevLink] are the two
 edges that attach W to the external face cycle of the containing bicomp.
 ********************************************************************/

void _MergeVertex(graphP theGraph, int W, int WPrevLink, int R)
{
int  J, JTwin;
int  e_w, e_r, e_ext;

	 gp_LogLine(gp_MakeLogStr4("graphEmbed.c/_MergeVertex() W=%d, W_in=%d, R=%d, R_out=%d",
			 W, WPrevLink, R, 1^WPrevLink));

     // All arcs leading into R from its neighbors must be changed
     // to say that they are leading into W
     J = gp_GetFirstArc(theGraph, R);
     while (gp_IsArc(theGraph, J))
     {
         JTwin = gp_GetTwinArc(theGraph, J);
         theGraph->G[JTwin].v = W;

    	 J = gp_GetNextArc(theGraph, J);
     }

     // Obtain the edge records involved in the list union
     e_w = gp_GetArc(theGraph, W, WPrevLink);
     e_r = gp_GetArc(theGraph, R, 1^WPrevLink);
     e_ext = gp_GetArc(theGraph, R, WPrevLink);

     // If W has any edges, then join the list with that of R
     if (gp_IsArc(theGraph, e_w))
     {
         // The WPrevLink arc of W is e_w, so the 1^WPrevLink arc in e_w leads back to W.
         // Now it must lead to e_r.  Likewise, e_r needs to lead back to e_w with the
         // opposing link, which is WPrevLink
         // Note that the adjacency lists of W and R are guaranteed non-empty, which is
         // why these linkages can be made without NIL tests.
         gp_SetAdjacentArc(theGraph, e_w, 1^WPrevLink, e_r);
         gp_SetAdjacentArc(theGraph, e_r, WPrevLink, e_w);

         // Cross-link W's WPrevLink arc and the 1^WPrevLink arc in e_ext
         gp_SetArc(theGraph, W, WPrevLink, e_ext);
         gp_SetAdjacentArc(theGraph, e_ext, 1^WPrevLink, gp_AdjacencyListEndMark(W));
     }
     // Otherwise, W just receives R's list.  This happens, for example, on a
     // DFS tree root vertex during JoinBicomps()
     else
     {
         // Cross-link W's 1^WPrevLink arc and the WPrevLink arc in e_r
         gp_SetArc(theGraph, W, 1^WPrevLink, e_r);
         gp_SetAdjacentArc(theGraph, e_r, WPrevLink, gp_AdjacencyListEndMark(W));

         // Cross-link W's WPrevLink arc and the 1^WPrevLink arc in e_ext
         gp_SetArc(theGraph, W, WPrevLink, e_ext);
         gp_SetAdjacentArc(theGraph, e_ext, 1^WPrevLink, gp_AdjacencyListEndMark(W));
     }

     // Erase the entries in R, which is a root copy that is no longer needed
     theGraph->functions.fpInitGraphNode(theGraph, R);
}

/********************************************************************
 _MergeBicomps()

 Merges all biconnected components at the cut vertices indicated by
 entries on the stack.

 theGraph contains the stack of bicomp roots and cut vertices to merge

 I, RootVertex, W and WPrevLink are not used in this routine, but are
          used by overload extensions

 Returns OK, but an extension function may return a value other than
         OK in order to cause Walkdown to terminate immediately.
********************************************************************/

int  _MergeBicomps(graphP theGraph, int I, int RootVertex, int W, int WPrevLink)
{
int  R, Rout, Z, ZPrevLink, J;
int  theList, RootID_DFSChild;
int  extFaceVertex;

     while (sp_NonEmpty(theGraph->theStack))
     {
         sp_Pop2(theGraph->theStack, R, Rout);
         sp_Pop2(theGraph->theStack, Z, ZPrevLink);

         /* The external faces of the bicomps containing R and Z will
            form two corners at Z.  One corner will become part of the
            internal face formed by adding the new back edge. The other
            corner will be the new external face corner at Z.
            We first want to update the links at Z to reflect this. */

         extFaceVertex = theGraph->extFace[R].vertex[1^Rout];
         theGraph->extFace[Z].vertex[ZPrevLink] = extFaceVertex;

         if (theGraph->extFace[extFaceVertex].vertex[0] == theGraph->extFace[extFaceVertex].vertex[1])
            theGraph->extFace[extFaceVertex].vertex[Rout ^ theGraph->extFace[extFaceVertex].inversionFlag] = Z;
         else
            theGraph->extFace[extFaceVertex].vertex[theGraph->extFace[extFaceVertex].vertex[0] == R ? 0 : 1] = Z;

         /* If the path used to enter Z is opposed to the path
            used to exit R, then we have to flip the bicomp
            rooted at R, which we signify by inverting R
            then setting the sign on its DFS child edge to
            indicate that its descendants must be flipped later */

         if (ZPrevLink == Rout)
         {
             Rout = 1^ZPrevLink;

             if (gp_GetFirstArc(theGraph, R) != gp_GetLastArc(theGraph, R))
                _InvertVertex(theGraph, R);

             J = gp_GetFirstArc(theGraph, R);
             while (gp_IsArc(theGraph, J))
             {
                 if (theGraph->G[J].type == EDGE_DFSCHILD)
                 {
                	 // A bicomp root edge cannot be inverted in the core planarity algorithm
                	 // but extensions may perform edge reductions on tree edges, resulting in
                	 // an inversion sign being promoted to the root edge.  So, now we reverse
                	 // the inversion flag on the root edge if the bicomp root must be
                	 // inverted before it is merged.
                	 if (GET_EDGEFLAG_INVERTED(theGraph, J))
                		 CLEAR_EDGEFLAG_INVERTED(theGraph, J);
                	 else
                		 SET_EDGEFLAG_INVERTED(theGraph, J);
                     break;
                 }

                 J = gp_GetNextArc(theGraph, J);
             }
         }

         // The endpoints of a bicomp's "root edge" are the bicomp root R and a
         // DFS child of the parent copy of the bicomp root R.
         // The GraphNode location of the root vertices is in the range N to 2N-1
         // at the offset indicated by the associated DFS child.  So, the location
         // of the root vertex R, less N, is the location of the DFS child and also
         // a convenient identifier for the bicomp root.
         RootID_DFSChild = R - theGraph->N;

         /* R is no longer pertinent to Z since we are about to
            merge R into Z, so we delete R from its pertinent
            bicomp list (Walkdown gets R from the head of the list). */

         theList = theGraph->V[Z].pertinentBicompList;
         theList = LCDelete(theGraph->BicompLists, theList, RootID_DFSChild);
         theGraph->V[Z].pertinentBicompList = theList;

         /* As a result of the merge, the DFS child of Z must be removed
            from Z's SeparatedDFSChildList because the child has just
            been joined directly to Z, rather than being separated by a
            root copy. */

         theList = theGraph->V[Z].separatedDFSChildList;
         theList = LCDelete(theGraph->DFSChildLists, theList, RootID_DFSChild);
         theGraph->V[Z].separatedDFSChildList = theList;

         /* Now we push R into Z, eliminating R */

         _MergeVertex(theGraph, Z, ZPrevLink, R);
     }

     return OK;
}

/********************************************************************
 _WalkUp()
 I is the vertex currently being embedded
 J is the forward arc to the descendant W on which the Walkup begins

 The Walkup establishes pertinence for step I.  It marks W as
 'adjacentTo' I so that the Walkdown will embed an edge to W when
 it is encountered.

 The Walkup also determines the pertinent child bicomps that should be
 set up as a result of the need to embed edge (I, W). It does this by
 recording the pertinent child biconnected components of all cut
 vertices between W and the child of I that is a descendant of W.
 Note that it stops the traversal if it finds a visited flag set to I,
 which indicates that a prior walkup call in step I has already done
 the work.

 Zig and Zag are so named because one goes around one side of a
 bicomp and the other goes around the other side, yet we have
 as yet no notion of orientation for the bicomp.
 The edge J from vertex I gestures to an adjacent descendant vertex W
 (possibly in some other bicomp).  Zig and Zag start out at W.
 They go around alternate sides of the bicomp until its root is found.
 Recall that the root vertex is just a copy in region N to 2N-1.
 We want to hop from the root copy to the parent copy of the vertex
 in order to record which bicomp we just came from and also to continue
 the walk-up to vertex I.
 If the parent copy actually is I, then the walk-up is done.
 ********************************************************************/

void _WalkUp(graphP theGraph, int I, int J)
{
int  Zig, Zag, ZigPrevLink, ZagPrevLink;
int  N, R, ParentCopy, nextVertex, W;
int  RootID_DFSChild, BicompList;

     W = theGraph->G[J].v;
     theGraph->V[W].adjacentTo = J;

     /* Shorthand for N, due to frequent use */

     N = theGraph->N;

     /* Start at the vertex W and walk around the both sides of the external face
        of a bicomp until we get back to vertex I. */

     Zig = Zag = W;
     ZigPrevLink = 1;
     ZagPrevLink = 0;

     while (Zig != I)
     {
        /* A previous walk-up may have been this way already */

        if (theGraph->G[Zig].visited == I) break;
        if (theGraph->G[Zag].visited == I) break;

        /* Mark the current vertices as visited during the embedding of vertex I. */

        theGraph->G[Zig].visited = I;
        theGraph->G[Zag].visited = I;

        /* Determine whether either Zig or Zag has landed on a bicomp root */

        if (Zig >= N) R = Zig;
        else if (Zag >= N) R = Zag;
        else R = NIL;

        // If we have a bicomp root, then we want to hop up to the parent copy and
        // record a pertinent child bicomp.
        // Prepends if the bicomp is internally active, appends if externally active.

        if (R != NIL)
        {
            // The endpoints of a bicomp's "root edge" are the bicomp root R and a
            // DFS child of the parent copy of the bicomp root R.
            // The GraphNode location of the root vertices is in the range N to 2N-1
            // at the offset indicated by the associated DFS child.  So, the location
            // of the root vertex R, less N, is the location of the DFS child and also
            // a convenient identifier for the bicomp root.
            RootID_DFSChild = R - N;

            // It is extra unnecessary work to record pertinent bicomps of I
            if ((ParentCopy = theGraph->V[RootID_DFSChild].DFSParent) != I)
            {
                 // Get the BicompList of the parent copy vertex.
                 BicompList = theGraph->V[ParentCopy].pertinentBicompList;

                 /* Put the new root vertex in the BicompList.  It is prepended if internally
                    active and appended if externally active so that all internally
                    active bicomps are processed before any externally active bicomps
                    by virtue of storage.

                    NOTE: The activity status of a bicomp is computed using the lowpoint of
                            the DFS child in the bicomp's root edge because we want to know
                            whether the DFS child or any of its descendants are joined by a
                            back edge to ancestors of I. If so, then the bicomp rooted
                            at RootVertex must contain an externally active vertex so the
                            bicomp must be kept on the external face. */

                 if (theGraph->V[RootID_DFSChild].Lowpoint < I)
                      BicompList = LCAppend(theGraph->BicompLists, BicompList, RootID_DFSChild);
                 else BicompList = LCPrepend(theGraph->BicompLists, BicompList, RootID_DFSChild);

                 /* The head node of the parent copy vertex's bicomp list may have changed, so
                    we assign the head of the modified list as the vertex's pertinent
                    bicomp list */

                 theGraph->V[ParentCopy].pertinentBicompList = BicompList;
            }

            Zig = Zag = ParentCopy;
            ZigPrevLink = 1;
            ZagPrevLink = 0;
        }

        /* If we did not encounter a bicomp root, then we continue traversing the
            external face in both directions. */

        else
        {
            nextVertex = theGraph->extFace[Zig].vertex[1^ZigPrevLink];
            ZigPrevLink = theGraph->extFace[nextVertex].vertex[0] == Zig ? 0 : 1;
            Zig = nextVertex;

            nextVertex = theGraph->extFace[Zag].vertex[1^ZagPrevLink];
            ZagPrevLink = theGraph->extFace[nextVertex].vertex[0] == Zag ? 0 : 1;
            Zag = nextVertex;
        }
     }
}


/********************************************************************
 _HandleBlockedDescendantBicomp()
 The core planarity/outerplanarity algorithm handles the blockage
 by pushing the root of the blocked bicomp onto the top of the stack
 because it is the central focus for obstruction minor A.
 Then NONEMBEDDABLE is returned so that the WalkDown can terminate,
 and the embedder can proceed to isolate the obstruction.
 Some algorithms may be able to clear the blockage, in which case
 a function overload would set Rout, W and WPrevLink, then return OK
 to indicate that the WalkDown may proceed.

 NOTE: When returning OK (blockage cleared), the overload implementation
       should NOT call this base implementation nor otherwise push R
       onto the stack because the core WalkDown implementation will push
       the appropriate stack entries based on R, Rout, W and WPrevLink
       Similarly, when returning NONEMBEDDABLE, it is typically not
       necessary to call this base implementation because pushing
       the bicomp root R is not usually necessary, i.e. the overload
       implementation usually does all embed post-processing before
       returning NONEMBEDDABLE.

 Returns OK to proceed with WalkDown at W,
         NONEMBEDDABLE to terminate WalkDown of Root Vertex
         NOTOK for internal error
 ********************************************************************/

int  _HandleBlockedDescendantBicomp(graphP theGraph, int I, int RootVertex, int R, int *pRout, int *pW, int *pWPrevLink)
{
    sp_Push2(theGraph->theStack, R, 0);
	return NONEMBEDDABLE;
}

/********************************************************************
 _HandleInactiveVertex()
 ********************************************************************/

int  _HandleInactiveVertex(graphP theGraph, int BicompRoot, int *pW, int *pWPrevLink)
{
     int X = theGraph->extFace[*pW].vertex[1^*pWPrevLink];
     *pWPrevLink = theGraph->extFace[X].vertex[0] == *pW ? 0 : 1;
     *pW = X;

     return OK;
}

/********************************************************************
 _GetPertinentChildBicomp()
 Returns the root of a pertinent child bicomp for the given vertex.
 Note: internally active roots are prepended by _Walkup()
 ********************************************************************/

#define _GetPertinentChildBicomp(theGraph, W) \
        (theGraph->V[W].pertinentBicompList==NIL \
         ? NIL \
         : theGraph->V[W].pertinentBicompList + theGraph->N)

/********************************************************************
 _WalkDown()
 Consider a circular shape with small circles and squares along its perimeter.
 The small circle at the top the root vertex of the bicomp.  The other small
 circles represent internally active vertices, and the squares represent
 externally active vertices.  The root vertex is a root copy of I, the
 vertex currently being processed.

 The Walkup previously marked all vertices adjacent to I by setting their
 adjacentTo flags.  Basically, we want to walk down both external face
 paths emanating from RootVertex, embedding edges between the RootVertex
 (a root copy of vertex I) and descendants of vertex I that have the
 adjacentTo flag set.

 During each walk down, it is sometimes necessary to hop from a vertex
 to one of its child biconnected components in order to reach the desired
 vertices.  In such cases, the biconnected components are merged together
 such that adding the back edge forms a new proper face in the biconnected
 component rooted at RootVertex (which, again, is a root copy of I).

 The outer loop performs both walks, unless the first walk got all the way
 around to RootVertex (only happens when bicomp contains no external activity,
 such as when processing the last vertex), or when non-planarity is
 discovered (in a pertinent child bicomp such that the stack is non-empty).

 For the inner loop, each iteration visits a vertex W.  If W is adjacentTo I,
 we call MergeBicomps to merge the biconnected components whose cut vertices
 have been collecting in theStack.  Then, we add the back edge (RootVertex, W)
 and clear the adjacentTo flag in W.

 Next, we check whether W has a pertinent child bicomp.  If so, then we figure
 out which path down from the root of the child bicomp leads to the next vertex
 to be visited, and we push onto the stack information on the cut vertex and
 the paths used to enter into it and exit from it.  Alternately, if W
 had no pertinent child bicomps, then we check to see if it is inactive.
 If so, we find the next vertex along the external face, then short-circuit
 its inactive predecessor (under certain conditions).  Finally, if W is not
 inactive, but it has no pertinent child bicomps, then we already know its
 adjacentTo flag is clear so both criteria for internal activity also fail.
 Therefore, W must be a stopping vertex.

 A stopping vertex X is an externally active vertex that has no pertinent
 child bicomps and no unembedded back edge to the current vertex I.
 The inner loop of Walkdown stops walking when it reaches a stopping vertex X
 because if it were to proceed beyond X and embed a back edge, then X would be
 surrounded by the bounding cycle of the bicomp.  This is clearly incorrect
 because X has a path leading from it to an ancestor of I (which is why it's
 externally active), and this path would have to cross the bounding cycle.

 After the loop, if the stack is non-empty, then the Walkdown halted because
 it could not proceed down a pertinent child biconnected component along either
 path from its root, which is easily shown to be evidence of a K_3,3, so
 we break the outer loop.  The caller performs further tests to determine
 whether Walkdown has embedded all back edges.  If the caller does not embed
 all back edges to descendants of the root vertex after walking both RootSide
 0 then 1 in all bicomps containing a root copy of I, then the caller can
 conclude that the input graph is non-planar.

  Returns OK if all possible edges were embedded, NONEMBEDDABLE if less
          than all possible edges were embedded, and NOTOK for an internal
          code failure
 ********************************************************************/

int  _WalkDown(graphP theGraph, int I, int RootVertex)
{
int  RetVal, W, WPrevLink, R, Rout, X, XPrevLink, Y, YPrevLink, RootSide, RootEdgeChild;

#ifdef DEBUG
     // Resolves typical watch expressions
     R = RootVertex;
#endif

     RootEdgeChild = RootVertex - theGraph->N;

     sp_ClearStack(theGraph->theStack);

     for (RootSide = 0; RootSide < 2; RootSide++)
     {
         W = theGraph->extFace[RootVertex].vertex[RootSide];

         // If the main bicomp rooted by RootVertex is a single tree edge,
         // (always the case for core planarity) then the external face links
         // of W will be equal
         if (theGraph->extFace[W].vertex[0] == theGraph->extFace[W].vertex[1])
         {
        	 // In this case, we treat the bicomp external face as if it were
        	 // a cycle of two edges and as if RootVertex and W had the same
        	 // orientation. Thus, the edge record leading back to RootVertex
        	 // would be indicated by link[1^RootSide] as this is the reverse of
        	 // link[RootSide], which was used to exit RootVertex and get to W
             WPrevLink = 1^RootSide;
             // We don't bother with the inversionFlag here because WalkDown is
             // never called on a singleton bicomp with an inverted orientation
             // Before the first Walkdown, the bicomp truly is a single edge
             // with proper orientation, and an extension algorithm does call
             // Walkdown again in post-processing, it wouldn't do so on this
             // bicomp because a singleton, whether inverted or not, would no
             // longer be pertinent (until a future vertex step).
             // Thus only the inner loop below accommodates the inversionFlag
             // when it walks down to a *pertinent* child biconnected component
             //WPrevLink = theGraph->extFace[W].inversionFlag ? RootSide : 1^RootSide;
         }
         // Otherwise, Walkdown has been called on a bicomp with two distinct
         // external face paths from RootVertex (a possibility in extension
         // algorithms), so both external face path links from W do not indicate
         // the RootVertex.
         else
         {
        	 WPrevLink = theGraph->extFace[W].vertex[0] == RootVertex ? 0 : 1;
        	 if (theGraph->extFace[W].vertex[WPrevLink] != RootVertex)
        		 return NOTOK;
         }

         while (W != RootVertex)
         {
             /* If the vertex W is the descendant endpoint of an unembedded
                back edge to I, then ... */

             if (theGraph->V[W].adjacentTo != NIL)
             {
                /* Merge bicomps at cut vertices on theStack and add the back edge,
                    creating a new proper face. */

                if (sp_NonEmpty(theGraph->theStack))
                {
                    if ((RetVal = theGraph->functions.fpMergeBicomps(theGraph, I, RootVertex, W, WPrevLink)) != OK)
                        return RetVal;
                }
                theGraph->functions.fpEmbedBackEdgeToDescendant(theGraph, RootSide, RootVertex, W, WPrevLink);

                /* Clear W's AdjacentTo flag so we don't add another edge to W if
                    this invocation of Walkdown visits W again later (and more
                    generally, so that no more back edges to W are added until
                    a future Walkup sets the flag to non-NIL again). */

                theGraph->V[W].adjacentTo = NIL;
             }

             /* If there is a pertinent child bicomp, then we need to push it onto the stack
                along with information about how we entered the cut vertex and how
                we exit the root copy to get to the next vertex. */

             if (theGraph->V[W].pertinentBicompList != NIL)
             {
                 sp_Push2(theGraph->theStack, W, WPrevLink);
                 R = _GetPertinentChildBicomp(theGraph, W);

                 /* Get next active vertices X and Y on ext. face paths emanating from R */

                 X = theGraph->extFace[R].vertex[0];
                 XPrevLink = theGraph->extFace[X].vertex[1]==R ? 1 : 0;
                 Y = theGraph->extFace[R].vertex[1];
                 YPrevLink = theGraph->extFace[Y].vertex[0]==R ? 0 : 1;

                 /* If this is a bicomp with only two ext. face vertices, then
                    it could be that the orientation of the non-root vertex
                    doesn't match the orientation of the root due to our relaxed
                    orientation method. */

                 if (X == Y && theGraph->extFace[X].inversionFlag)
                 {
                     XPrevLink = 0;
                     YPrevLink = 1;
                 }

                 /* Now we implement the Walkdown's simple path selection rules!
                    If either X or Y is internally active (pertinent but not
                    externally active), then we pick it first.  Otherwise,
                    we choose a pertinent vertex. If neither are pertinent,
                    then we let a handler decide.  The default handler for
                    core planarity/outerplanarity decides to stop the WalkDown
                    with the current blocked bicomp at the top of the stack. */

                 if (_VertexActiveStatus(theGraph, X, I) == VAS_INTERNAL)
                 {
                      W = X;
                      WPrevLink = XPrevLink;
                      Rout = 0;
                 }
                 else if (_VertexActiveStatus(theGraph, Y, I) == VAS_INTERNAL)
                 {
                      W = Y;
                      WPrevLink = YPrevLink;
                      Rout = 1;
                 }
                 else if (PERTINENT(theGraph, X))
                 {
                      W = X;
                      WPrevLink = XPrevLink;
                      Rout = 0;
                 }
                 else if (PERTINENT(theGraph, Y))
                 {
                	 W = Y;
                     WPrevLink = YPrevLink;
                     Rout = 1;
                 }
                 else
                 {
                	 // Both the X and Y sides of the bicomp are blocked.
                	 // Let the application decide whether it can unblock the bicomp.
                	 // The core planarity embedder simply pushes (R, 0) onto the top of
                	 // the stack and returns NONEMBEDDABLE, which causes a return here
                	 // and enables isolation of planarity/outerplanary obstruction minor A
                     if ((RetVal = theGraph->functions.fpHandleBlockedDescendantBicomp(theGraph, I, RootVertex, R, &Rout, &W, &WPrevLink)) != OK)
                         return RetVal;
                 }

                 sp_Push2(theGraph->theStack, R, Rout);
             }

             /* Skip inactive vertices, which will be short-circuited
                later by our fast external face linking method (once
                upon a time, we added false edges called short-circuit
                edges to eliminate inactive vertices, but the extFace
                links can do the same job and also give us the ability
                to more quickly test planarity without creating an embedding). */

             else if (_VertexActiveStatus(theGraph, W, I) == VAS_INACTIVE)
             {
                 if (theGraph->functions.fpHandleInactiveVertex(theGraph, RootVertex, &W, &WPrevLink) != OK)
                     return NOTOK;
             }

             /* At this point, we know that W is not inactive, but its adjacentTo flag
                is clear, and it has no pertinent child bicomps.  Therefore, it
                is an externally active stopping vertex. */

             else break;
         }

         /* We short-circuit the external face of the bicomp by hooking the root
            to the terminating externally active vertex so that inactive vertices
            are not visited in future iterations.  This setting obviates the need
            for those short-circuit edges mentioned above.

            NOTE: We skip the step if the stack is non-empty since in that case
                    we did not actually merge the bicomps necessary to put
                    W and RootVertex into the same bicomp. */

         theGraph->extFace[RootVertex].vertex[RootSide] = W;
         theGraph->extFace[W].vertex[WPrevLink] = RootVertex;

         /* If the bicomp is reduced to having only two external face vertices
             (the root and W), then we need to record whether the orientation
             of W is inverted relative to the root.  This is used later when a
             future Walkdown descends to and merges the bicomp containing W.
             Going from the root to W, we only get the correct WPrevLink if
             we know whether or not W is inverted.
             NOTE: We clear the flag because it may have been set in W if W
                 previously became part of a bicomp with only two ext. face
                 vertices, but then was flipped and merged into a larger bicomp
                 that is now again becoming a bicomp with only two ext. face vertices. */

         if (theGraph->extFace[W].vertex[0] == theGraph->extFace[W].vertex[1] &&
             WPrevLink == RootSide)
              theGraph->extFace[W].inversionFlag = 1;
         else theGraph->extFace[W].inversionFlag = 0;

         /* If we got back around to the root, then all edges
            are embedded, so we stop. */

         if (W == RootVertex)
             break;
     }

     return OK;
}


/********************************************************************
 gp_Embed()

  First, a DFS tree is created in the graph (if not already done).
  Then, the graph is sorted by DFI.

  Either a planar embedding is created in theGraph, or a Kuratowski
  subgraph is isolated.  Either way, theGraph remains sorted by DFI
  since that is the most common desired result.  The original vertex
  numbers are available in the 'v' members of the vertex graph nodes.
  Moreover, gp_SortVertices() can be invoked to put the vertices in
  the order of the input graph, at which point the 'v' members of the
  vertex graph nodes will contain the vertex DFIs.

 return OK if the embedding was successfully created or no subgraph
            homeomorphic to a topological obstruction was found.

        NOTOK on internal failure

        NONEMBEDDABLE if the embedding couldn't be created due to
                the existence of a subgraph homeomorphic to a
                topological obstruction.

  For core planarity, OK is returned when theGraph contains a planar
  embedding of the input graph, and NONEMBEDDABLE is returned when a
  subgraph homeomorphic to K5 or K3,3 has been isolated in theGraph.

  Extension modules can overload functions used by gp_Embed to achieve
  alternate algorithms.  In those cases, the return results are
  similar.  For example, a K3,3 search algorithm would return
  NONEMBEDDABLE if it finds the K3,3 obstruction, and OK if the graph
  is planar or only contains K5 homeomorphs.  Similarly, an
  outerplanarity module can return OK for an outerplanar embedding or
  NONEMBEDDABLE when a subgraph homeomorphic to K2,3 or K4 has been
  isolated.

  The algorithm extension for gp_Embed() is encoded in the embedFlags,
  and the details of the return value can be found in the extension
  module that defines the embedding flag.

 ********************************************************************/

int gp_Embed(graphP theGraph, int embedFlags)
{
int N, I, J, child;
int RetVal = OK;

    /* Basic parameter checks */

    if (theGraph==NULL)
    	return NOTOK;

    /* A little shorthand for the size of the graph */

    N = theGraph->N;

    /* Preprocessing */

    theGraph->embedFlags = embedFlags;

    if (gp_CreateDFSTree(theGraph) != OK)
        return NOTOK;

    if (!(theGraph->internalFlags & FLAGS_SORTEDBYDFI))
        if (gp_SortVertices(theGraph) != OK)
            return NOTOK;

    if (gp_LowpointAndLeastAncestor(theGraph) != OK)
    	return NOTOK;

    _CreateSortedSeparatedDFSChildLists(theGraph);

    if (theGraph->functions.fpCreateFwdArcLists(theGraph) != OK)
        return NOTOK;

    theGraph->functions.fpCreateDFSTreeEmbedding(theGraph);

    /* In reverse DFI order, process each vertex by embedding its
         the 'back edges' from the vertex to its DFS descendants. */

    for (I = 0; I < theGraph->edgeOffset; I++)
        theGraph->G[I].visited = N;

    for (I = theGraph->N-1; I >= 0; I--)
    {
          RetVal = OK;

          /* Do the Walkup for each cycle edge from I to a DFS descendant W. */

          J = theGraph->V[I].fwdArcList;
          while (J != NIL)
          {
        	  theGraph->functions.fpWalkUp(theGraph, I, J);

              J = gp_GetNextArc(theGraph, J);
              if (J == theGraph->V[I].fwdArcList)
                  J = NIL;
          }

          /* For each DFS child C of the current vertex with a pertinent
                child bicomp, do a Walkdown on each side of the bicomp rooted
                by tree edge (R, C), where R is a root copy of the current
                vertex stored at C+N and uniquely associated with the bicomp
                containing C. (NOTE: if C has no pertinent child bicomps, then
                there are no cycle edges from I to descendants of C). */

          child = theGraph->V[I].separatedDFSChildList;
          while (child != NIL)
          {
              if (theGraph->V[child].pertinentBicompList != NIL)
              {
                  // _Walkdown returns OK even if it couldn't embed all
                  // back edges from I to the subtree rooted by child
                  // It only returns NONEMBEDDABLE when it was blocked
            	  // on a descendant bicomp with stopping vertices along
            	  // both external face paths emanating from the bicomp root
            	  // Some extension algorithms are able to clear some such
            	  // blockages with a reduction, and those algorithms only
            	  // return NONEMBEDDABLE when unable to clear the blockage
                  if ((RetVal = theGraph->functions.fpWalkDown(theGraph, I, child + N)) != OK)
                  {
                      if (RetVal == NONEMBEDDABLE)
                    	  break;
                      else
                    	  return NOTOK;
                  }
              }
              child = LCGetNext(theGraph->DFSChildLists,
                                theGraph->V[I].separatedDFSChildList, child);
          }

          /* If the Walkdown sequence is completed but not all forward edges
             are embedded or an explicit NONEMBEDDABLE result was returned,
             then the graph is not planar/outerplanar.
             The handler below is invoked because some extension algorithms are
             able to clear the blockage to planarity/outerplanarity and continue
             the embedder iteration loop (they return OK below).
             The default implementation simply returns NONEMBEDDABLE, which stops
             the embedding process. */

          if (theGraph->V[I].fwdArcList != NIL || RetVal == NONEMBEDDABLE)
          {
              RetVal = theGraph->functions.fpHandleBlockedEmbedIteration(theGraph, I);
              if (RetVal != OK)
                  break;
          }
    }

    /* Postprocessing to orient the embedding and merge any remaining separated bicomps,
       or to isolate an obstruction to planarity/outerplanarity.  Some extension algorithms
       either do nothing if they have already isolated a subgraph of interest, or they may
       do so now based on information collected by their implementations of
       HandleBlockedDescendantBicomp or HandleBlockedEmbedIteration */

    return theGraph->functions.fpEmbedPostprocess(theGraph, I, RetVal);
}

/********************************************************************
 HandleBlockedEmbedIteration()

  At the end of each embedding iteration, this function is invoked
  if there are any unembedded cycle edges from the current vertex I
  to its DFS descendants. Specifically, the forward arc list of I is
  non-empty at the end of the edge addition processing for I.

  We return NONEMBEDDABLE to cause iteration to stop because the
  graph is non-planar if any edges could not be embedded.

  Extensions may overload this function and decide to proceed with or
  halt embedding iteration for application-specific reasons.
  For example, a search for K_{3,3} homeomorphs could reduce an
  isolated K5 homeomorph to something that can be ignored, and then
  return OK in order to continue the planarity algorithm in order to
  search for a K_{3,3} homeomorph elsewhere in the graph.  On the
  other hand, if such an algorithm found a K_{3,3} homeomorph,
  perhaps alone or perhaps entangled with the K5 homeomorph, it would
  return NONEMBEDDABLE since there is no need to continue with
  embedding iterations once the desired embedding obstruction is found.

  If this function returns OK, then embedding will proceed to the
  next iteration, or return OK if it finished the last iteration.

  If this function returns NONEMBEDDABLE, then the embedder will
  stop iteration and return NONEMBEDDABLE.  Note that the function
  _EmbedPostprocess() is still called in this case, allowing for
  further processing of the non-embeddable result, e.g. isolation
  of the desired embedding obstruction.

  This function can return NOTOK to signify an internal error.
 ********************************************************************/

int  _HandleBlockedEmbedIteration(graphP theGraph, int I)
{
     return NONEMBEDDABLE;
}

/********************************************************************
 _EmbedPostprocess()

 After the loop that embeds the cycle edges from each vertex to its
 DFS descendants, this method is invoked to postprocess the graph.
 If the graph is planar, then a consistent orientation is imposed
 on the vertices of the embedding, and any remaining separated
 biconnected components are joined together.
 If the graph is non-planar, then a subgraph homeomorphic to K5
 or K3,3 is isolated.
 Extensions may override this function to provide alternate
 behavior.

  @param theGraph - the graph ready for postprocessing
  @param I - the last vertex processed by the edge embedding loop
  @param edgeEmbeddingResult -
         OK if all edge embedding iterations returned OK
         NONEMBEDDABLE if an embedding iteration failed to embed
             all edges for a vertex

  @return NOTOK on internal failure
          NONEMBEDDABLE if a subgraph homeomorphic to a topological
              obstruction is isolated in the graph
          OK otherwise (for example if the graph contains a
             planar embedding or if a desired topological obstruction
             was not found)

 *****************************************************************/

int  _EmbedPostprocess(graphP theGraph, int I, int edgeEmbeddingResult)
{
int  RetVal = edgeEmbeddingResult;

    /* If an embedding was found, then post-process the embedding structure
        to eliminate root copies and give a consistent orientation to all vertices. */

    if (edgeEmbeddingResult == OK)
    {
    	if (_OrientVerticesInEmbedding(theGraph) != OK ||
    		_JoinBicomps(theGraph) != OK)
    		RetVal = NOTOK;
    }

    /* If the graph was found to be unembeddable, then we want to isolate an
        obstruction.  But, if a search flag was set, then we have already
        found a subgraph with the desired structure, so no further work is done. */

    else if (edgeEmbeddingResult == NONEMBEDDABLE)
    {
        if (theGraph->embedFlags == EMBEDFLAGS_PLANAR)
        {
            if (_IsolateKuratowskiSubgraph(theGraph, I, NIL) != OK)
                RetVal = NOTOK;
        }
        else if (theGraph->embedFlags == EMBEDFLAGS_OUTERPLANAR)
        {
            if (_IsolateOuterplanarObstruction(theGraph, I, NIL) != OK)
                RetVal = NOTOK;
        }
    }

    return RetVal;
}

/********************************************************************
 _OrientVerticesInEmbedding()

 Each vertex will then have an orientation, either clockwise or
 counterclockwise.  All vertices in each bicomp need to have the
 same orientation.
 This method clears the stack, and the stack is clear when it
 is finished.
 Returns OK on success, NOTOK on implementation failure.
 ********************************************************************/

int  _OrientVerticesInEmbedding(graphP theGraph)
{
int  R=0, edgeOffset = theGraph->edgeOffset;

     sp_ClearStack(theGraph->theStack);

/* Run the array of root copy vertices.  For each that is not defunct
        (i.e. has not been merged during embed), we orient the vertices
        in the bicomp for which it is the root vertex. */

     for (R = theGraph->N; R < edgeOffset; R++)
     {
          if (gp_IsArc(theGraph, gp_GetFirstArc(theGraph, R)))
          {
        	  if (_OrientVerticesInBicomp(theGraph, R, 0) != OK)
        		  return NOTOK;
          }
     }
     return OK;
}

/********************************************************************
 _OrientVerticesInBicomp()
  As a result of the work done so far, the edges around each vertex have
 been put in order, but the orientation may be counterclockwise or
 clockwise for different vertices within the same bicomp.
 We need to reverse the orientations of those vertices that are not
 oriented the same way as the root of the bicomp.

 During embedding, a bicomp with root edge (v', c) may need to be flipped.
 We do this by inverting the root copy v' and implicitly inverting the
 orientation of the vertices in the subtree rooted by c by assigning -1
 to the sign of the DFSCHILD edge record leading to c.

 We now use these signs to help propagate a consistent vertex orientation
 throughout all vertices that have been merged into the given bicomp.
 The bicomp root contains the orientation to be imposed on all parent
 copy vertices.  We perform a standard depth first search to visit each
 vertex.  A vertex must be inverted if the product of the edge signs
 along the tree edges between the bicomp root and the vertex is -1.

 Finally, the PreserveSigns flag, if set, performs the inversions
 but does not change any of the edge signs.  This allows a second
 invocation of this function to restore the state of the bicomp
 as it was before the first call.

 This method uses the stack but preserves whatever may have been
 on it.  In debug mode, it will return NOTOK if the stack overflows.
 This method pushes at most two integers per vertext in the bicomp.

 Returns OK on success, NOTOK on implementation failure.
 ********************************************************************/

int  _OrientVerticesInBicomp(graphP theGraph, int BicompRoot, int PreserveSigns)
{
int  V, J, invertedFlag;
int  stackBottom = sp_GetCurrentSize(theGraph->theStack);

     sp_Push2(theGraph->theStack, BicompRoot, 0);

     while (sp_GetCurrentSize(theGraph->theStack) > stackBottom)
     {
         /* Pop a vertex to orient */
         sp_Pop2(theGraph->theStack, V, invertedFlag);

         /* Invert the vertex if the inverted flag is set */
         if (invertedFlag)
             _InvertVertex(theGraph, V);

         /* Push the vertex's DFS children that are in the bicomp */
         J = gp_GetFirstArc(theGraph, V);
         while (gp_IsArc(theGraph, J))
         {
             if (theGraph->G[J].type == EDGE_DFSCHILD)
             {
                 sp_Push2(theGraph->theStack, theGraph->G[J].v,
                		  invertedFlag ^ GET_EDGEFLAG_INVERTED(theGraph, J));

                 if (!PreserveSigns)
                	 CLEAR_EDGEFLAG_INVERTED(theGraph, J);
             }

             J = gp_GetNextArc(theGraph, J);
         }
     }
     return OK;
}

/********************************************************************
 _JoinBicomps()
 The embedding algorithm works by only joining bicomps once the result
 forms a larger bicomp.  However, if the original graph was separable
 or disconnected, then the result of the embed function will be a
 graph that contains each bicomp as a distinct entity.  The root of
 each bicomp will be in the region N to 2N-1.  This function merges
 the bicomps into one connected graph.
 ********************************************************************/

int  _JoinBicomps(graphP theGraph)
{
int  R, N, edgeOffset=theGraph->edgeOffset;

     for (R=N=theGraph->N; R < edgeOffset; R++)
          if (gp_IsArc(theGraph, gp_GetFirstArc(theGraph, R)))
              _MergeVertex(theGraph, theGraph->V[R-N].DFSParent, 0, R);

     return OK;
}

/****************************************************************************
 _OrientExternalFacePath()

 The vertices along the path (v ... w) are assumed to be degree two vertices
 in an external face path connecting u and x.  This method imparts the
 orientation of u and x onto the vertices v ... w.
 The work done is on the order of the path length.
 Returns OK if the external face path was oriented, NOTOK on implementation
 error (i.e. if a condition arises providing the path is not on the
 external face).
 ****************************************************************************/

int  _OrientExternalFacePath(graphP theGraph, int u, int v, int w, int x)
{
int  e_u, e_v, e_ulink, e_vlink;

    // Get the edge record in u that indicates v; uses the twinarc method to
    // ensure the cost is dominated by the degree of v (which is 2), not u
    // (which can be any degree).
    e_u = gp_GetTwinArc(theGraph, gp_GetNeighborEdgeRecord(theGraph, v, u));

    do {
        // Get the external face link in vertex u that indicates the
        // edge e_u which connects to the next vertex v in the path
    	// As a sanity check, we determine whether e_u is an
    	// external face edge, because there would be an internal
    	// implementation error if not
    	if (gp_GetFirstArc(theGraph, u) == e_u)
    		e_ulink = 0;
    	else if (gp_GetLastArc(theGraph, u) == e_u)
    		e_ulink = 1;
    	else return NOTOK;

        v = theGraph->G[e_u].v;

        // Now get the external face link in vertex v that indicates the
        // edge e_v which connects back to the prior vertex u.
        e_v = gp_GetTwinArc(theGraph, e_u);

    	if (gp_GetFirstArc(theGraph, v) == e_v)
    		e_vlink = 0;
    	else if (gp_GetLastArc(theGraph, v) == e_v)
    		e_vlink = 1;
    	else return NOTOK;

        // The vertices u and v are inversely oriented if they
        // use the same link to indicate the edge [e_u, e_v].
        if (e_vlink == e_ulink)
        {
            _InvertVertex(theGraph, v);
            e_vlink = 1^e_vlink;
        }

        // This update of the extFace short-circuit is polite but unnecessary.
        // This orientation only occurs once we know we can isolate a K_{3,3},
        // at which point the extFace data structure is not used.
        theGraph->extFace[u].vertex[e_ulink] = v;
        theGraph->extFace[v].vertex[e_vlink] = u;

        u = v;
        e_u = gp_GetArc(theGraph, v, 1^e_vlink);
    } while (u != x);

    return OK;
}

/****************************************************************************
 _SetVisitedOnPath()
 This method sets the visited flags to 'visited' on the vertices and edges on
 the path (u, v, ..., w, x) in which all vertices except the endpoints u and x
 are degree 2.  This method avoids performing more than constant work at the
 path endpoints u and x, so the total work is on the order of the path length.

 Returns OK on success, NOTOK on internal failure
 ****************************************************************************/

int  _SetVisitedOnPath(graphP theGraph, int u, int v, int w, int x, int visited)
{
int  e, eTwin, pathLength=0;

     // We want to exit u from e, but we get eTwin first here in order to avoid
     // work, in case the degree of u is greater than 2.
     eTwin = gp_GetNeighborEdgeRecord(theGraph, v, u);
     if (!gp_IsArc(theGraph, eTwin))
    	 return NOTOK;
     e = gp_GetTwinArc(theGraph, eTwin);

     v = u;

     do {
    	 // Mark the vertex and the exiting edge
         theGraph->G[v].visited = visited;
         theGraph->G[e].visited = visited;
         theGraph->G[eTwin].visited = visited;

    	 // Get the next vertex
         v = theGraph->G[e].v;
         e = gp_GetNextArcCircular(theGraph, eTwin);
         eTwin = gp_GetTwinArc(theGraph, e);

         // A simple reality check on the preconditions of this method
         if (++pathLength > theGraph->N)
        	 return NOTOK;

     } while (v != x);

     // Mark the last vertex with 'visited'
     theGraph->G[x].visited = visited;

     return OK;
}
