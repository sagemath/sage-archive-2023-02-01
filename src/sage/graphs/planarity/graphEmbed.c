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

#define GRAPH_C

#include <stdlib.h>

#include "graph.h"

/* Imported functions */

extern void _InitGraphNode(graphP theGraph, int I);
extern void _FillVisitedFlags(graphP, int);

extern int _IsolateKuratowskiSubgraph(graphP theEmbedding, int I);

/* Private functions (some are exported to system only) */

void _CreateSortedSeparatedDFSChildLists(graphP theEmbedding);
void _CreateFwdArcLists(graphP theGraph);
void _CreateDFSTreeEmbedding(graphP theGraph);

void _EmbedBackEdgeToDescendant(graphP theEmbedding, int RootSide, int RootVertex, int W, int WPrevLink);

int  _GetNextVertexOnExternalFace(graphP theEmbedding, int curVertex, int *pPrevLink);

void _InvertVertex(graphP theEmbedding, int V);
void _MergeVertex(graphP theEmbedding, int W, int WPrevLink, int R);
void _MergeBicomps(graphP theEmbedding);

void _RecordPertinentChildBicomp(graphP theEmbedding, int I, int RootVertex);
#ifndef SPEED_MACROS
int  _GetPertinentChildBicomp(graphP theEmbedding, int W);
#endif

void _WalkUp(graphP theEmbedding, int I, int W);
void _WalkDown(graphP theEmbedding, int I, int RootVertex);

void _OrientVerticesInEmbedding(graphP theEmbedding);
void _OrientVerticesInBicomp(graphP theEmbedding, int BicompRoot, int PreserveSigns);
int  _JoinBicomps(graphP theEmbedding);

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

void  _CreateSortedSeparatedDFSChildLists(graphP theEmbedding)
{
int *buckets;
listCollectionP bin;
int I, J, N, DFSParent, theList;

     N = theEmbedding->N;
     buckets = theEmbedding->buckets;
     bin = theEmbedding->bin;

     /* Initialize the bin and all the buckets to be empty */

     LCReset(bin);
     for (I=0; I < N; I++)
          buckets[I] = NIL;

     /* For each vertex, add it to the bucket whose index is equal to
        the Lowpoint of the vertex. */

     for (I=0; I < N; I++)
     {
          J = theEmbedding->V[I].Lowpoint;
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
                  DFSParent = theEmbedding->V[J].DFSParent;

                  if (DFSParent != NIL && DFSParent != J)
                  {
                      theList = theEmbedding->V[DFSParent].separatedDFSChildList;
                      theList = LCAppend(theEmbedding->DFSChildLists, theList, J);
                      theEmbedding->V[DFSParent].separatedDFSChildList = theList;
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
 simplified by the fact that they have already been placed in link[1]
 succession.
 ********************************************************************/

void _CreateFwdArcLists(graphP theGraph)
{
int I, Jfirst, Jnext, Jlast;

#ifndef ORDER_EDGES
#error The edges must be ordered by the DFS; otherwise alternate code is needed here
#endif

    for (I=0; I < theGraph->N; I++)
    {
        Jfirst = theGraph->G[I].link[1];

        /* If the vertex has any forward edges, then ... */

        if (theGraph->G[Jfirst].type == EDGE_FORWARD)
        {
            /* Find the end of the forward edge list */

            Jnext = Jfirst;
            while (theGraph->G[Jnext].type == EDGE_FORWARD)
                Jnext = theGraph->G[Jnext].link[1];
            Jlast = theGraph->G[Jnext].link[0];

            /* Remove the forward edges from the adjacency list of I */

            theGraph->G[Jnext].link[0] = I;
            theGraph->G[I].link[1] = Jnext;

            /* Make a circular forward edge list */

            theGraph->V[I].fwdArcList = Jfirst;
            theGraph->G[Jfirst].link[0] = Jlast;
            theGraph->G[Jlast].link[1] = Jfirst;
        }
    }
}

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
    for (I=0, R=N; I < N; I++, R++)
    {
        if (theGraph->V[I].DFSParent == NIL)
        {
            theGraph->G[I].link[0] = theGraph->G[I].link[1] = I;
        }
        else
        {
             J = theGraph->G[I].link[0];
             while (theGraph->G[J].type != EDGE_DFSPARENT)
                 J = theGraph->G[J].link[0];

             theGraph->G[I].link[0] = theGraph->G[I].link[1] = J;
             theGraph->G[J].link[0] = theGraph->G[J].link[1] = I;
             theGraph->G[J].v = R;

             Jtwin = gp_GetTwinArc(theGraph, J);

             theGraph->G[R].link[0] = theGraph->G[R].link[1] = Jtwin;
             theGraph->G[Jtwin].link[0] = theGraph->G[Jtwin].link[1] = R;

             theGraph->extFace[R].link[0] = theGraph->extFace[R].link[1] = I;
             theGraph->extFace[I].link[0] = theGraph->extFace[I].link[1] = R;
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

    if (theGraph->V[parentCopy].fwdArcList == fwdArc)
    {
        if (theGraph->G[fwdArc].link[0] == fwdArc)
             theGraph->V[parentCopy].fwdArcList = NIL;
        else theGraph->V[parentCopy].fwdArcList = theGraph->G[fwdArc].link[0];
    }

    theGraph->G[theGraph->G[fwdArc].link[0]].link[1] = theGraph->G[fwdArc].link[1];
    theGraph->G[theGraph->G[fwdArc].link[1]].link[0] = theGraph->G[fwdArc].link[0];

    /* The forward arc is added to the adjacency list of the RootVertex. */

    theGraph->G[fwdArc].link[1^RootSide] = RootVertex;
    theGraph->G[fwdArc].link[RootSide] = theGraph->G[RootVertex].link[RootSide];
    theGraph->G[theGraph->G[RootVertex].link[RootSide]].link[1^RootSide] = fwdArc;
    theGraph->G[RootVertex].link[RootSide] = fwdArc;

    /* The back arc is added to the adjacency list of W. */

    theGraph->G[backArc].v = RootVertex;

    theGraph->G[backArc].link[1^WPrevLink] = W;
    theGraph->G[backArc].link[WPrevLink] = theGraph->G[W].link[WPrevLink];
    theGraph->G[theGraph->G[W].link[WPrevLink]].link[1^WPrevLink] = backArc;
    theGraph->G[W].link[WPrevLink] = backArc;

    /* Link the two endpoint vertices together on the external face */

    theGraph->extFace[RootVertex].link[RootSide] = W;
    theGraph->extFace[W].link[WPrevLink] = RootVertex;
}

/********************************************************************
 _VertexActiveStatus()
 Returns the active status of theVertex (VAS_INACTIVE, VAS_INTERNAL
 or VAS_EXTERNAL).
 This function is only called on vertices on the external face
 of a bicomp, and only the parent copy of a vertex.
 If the vertex has no (unembedded) edges to I or I's ancestors, and if it
 has no child bicomps with such connections, then it is classified inactive
 and can be short-circuited.  Note that the vertex may have DFS
 children in the same bicomp that are connected to I or I's ancestors,
 but when a DFS child joins its parent in the same bicomp, we remove
 that DFS child from the vertex's SeparatedDFSChildList because the child
 can represent itself when doing a walkdown of a bicomp.
 A vertex V should only be considered externally active if it has a
 least ancestor less than I or if it has a child bicomp containing
 such a least ancestor.  The latter test is done by comparing against
 the least Lowpoint from among all DFS children that have not yet merged
 into the parent bicomp of the vertex V.
 ********************************************************************/

#ifndef SPEED_MACROS

int  _VertexActiveStatus(graphP theEmbedding, int theVertex, int I)
{
int  leastLowpoint, DFSChild;

     if ((DFSChild=theEmbedding->V[theVertex].separatedDFSChildList) == NIL)
          leastLowpoint = theVertex;
     else leastLowpoint = theEmbedding->V[DFSChild].Lowpoint;

     if (leastLowpoint > theEmbedding->V[theVertex].leastAncestor)
         leastLowpoint = theEmbedding->V[theVertex].leastAncestor;

     if (leastLowpoint < I)
         return VAS_EXTERNAL;

     /* If we used (leastLowpoint == I) to return VAS_INTERNAL,
        a vertex would not switch from internally active to inactive
        until after step I, when the step variable is reduced.
        The notion of activity means 'pertinent to the future
        embedding of back edges'.  Once the back edges that
        cause a vertex to be internally active are embedded,
        we need the vertex to become inactive. */

     if (theEmbedding->V[theVertex].adjacentTo != NIL ||
         theEmbedding->V[theVertex].pertinentBicompList != NIL)
         return VAS_INTERNAL;

     return VAS_INACTIVE;
}

#endif

/********************************************************************
 _GetNextVertexOnExternalFace()
 Each vertex contains a link[0] and link[1] that link it into its
 list of edges.  If the vertex is on the external face, then the two
 edge nodes pointed to by link[0] and link[1] are also on the
 external face.  We want to take one of those edges to get to the
 next vertex on the external face.
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
 ********************************************************************/

int  _GetNextVertexOnExternalFace(graphP theEmbedding, int curVertex, int *pPrevLink)
{
int  arc, nextArc, nextVertex, newPrevLink;

     /* Exit curVertex from whichever link was not previously used to enter it */

     arc = theEmbedding->G[curVertex].link[1^(*pPrevLink)];

     nextArc = gp_GetTwinArc(theEmbedding, arc);

     nextVertex = theEmbedding->G[nextArc].link[newPrevLink=0];
     if (nextVertex >= 2*theEmbedding->N)
         nextVertex = theEmbedding->G[nextArc].link[newPrevLink=1];

     /* The setting above is how we exited an edge record to get to the
        next vertex.  The reverse pointer leads back from the vertex to
        the edge record. */

     newPrevLink = 1^newPrevLink;

     /* This if stmt assigns the new prev link that tells us which edge
        record was used to enter nextVertex (so that we exit from the
        opposing edge record).
        However, if we are in a singleton bicomp, then both links in nextVertex
        lead back to curVertex, so newPrevLink may get stop at the zero setting
        when it should become one.
        We want the two arcs of a singleton bicomp to act like a cycle, so the
        edge record given as the prev link for curVertex should be the same as
        the prev link for nextVertex.
        So, we only need to modify the prev link if the links in nextVertex
        are not equal. */

     if (theEmbedding->G[nextVertex].link[0] != theEmbedding->G[nextVertex].link[1])
         *pPrevLink = newPrevLink;

     return nextVertex;
}

/********************************************************************
 _InvertVertex()
 This function flips the orientation of a single vertex such that
 instead of using link[0] successors to go clockwise (or counterclockwise)
 around a vertex's adjacency list, link[1] successors would be used.
 The loop is constructed using do-while so we can swap the links
 in the vertex node as well as each arc node.
 ********************************************************************/

void _InvertVertex(graphP theEmbedding, int V)
{
int J, JTemp;

     J = V;
     do {
        JTemp = theEmbedding->G[J].link[0];
        theEmbedding->G[J].link[0] = theEmbedding->G[J].link[1];
        theEmbedding->G[J].link[1] = JTemp;

        J = theEmbedding->G[J].link[0];
     }  while (J >= 2*theEmbedding->N);

     JTemp = theEmbedding->extFace[V].link[0];
     theEmbedding->extFace[V].link[0] = theEmbedding->extFace[V].link[1];
     theEmbedding->extFace[V].link[1] = JTemp;
}

/********************************************************************
 _SetSignOfChildEdge()
 Finds the DFSCHILD edge of the vertex, and sets its sign.
 ********************************************************************/

void _SetSignOfChildEdge(graphP theEmbedding, int V, int sign)
{
int  J;

     J = theEmbedding->G[V].link[0];
     while (J >= 2*theEmbedding->N)
     {
         if (theEmbedding->G[J].type == EDGE_DFSCHILD)
         {
             theEmbedding->G[J].sign = sign;
             break;
         }

         J = theEmbedding->G[J].link[0];
     }
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
 adjacency list of W.  We set e_w to W's link[WPrevLink] and e_r to
 R's link[1^WPrevLink] so that e_w and e_r indicate W and R with
 opposing links, which become free to be cross-linked.  Finally,
 the edge record e_ext, set equal to R's link[WPrevLink], is the edge
 that, with e_r, held R to the external face.  Now, e_ext will be the
 new link[WPrevLink] edge record for W.  If e_w and e_r become part
 of a proper face, then e_ext and W's link[1^WPrevLink] are the two
 edges that hold W to the external face.
 ********************************************************************/

void _MergeVertex(graphP theEmbedding, int W, int WPrevLink, int R)
{
int  J, JTwin, N;
int  e_w, e_r, e_ext;

     N = theEmbedding->N;

     /* All arcs leading into R from its neighbors must be changed
        to say that they are leading into W */

     J = theEmbedding->G[R].link[0];
     while (J >= 2*N)
     {
         JTwin = gp_GetTwinArc(theEmbedding, J);
         theEmbedding->G[JTwin].v = W;

         J = theEmbedding->G[J].link[0];
     }

     /* Obtain the edge records involved in the circular list union */

     e_w = theEmbedding->G[W].link[WPrevLink];
     e_r = theEmbedding->G[R].link[1^WPrevLink];
     e_ext = theEmbedding->G[R].link[WPrevLink];

     /* WPrevLink leads away from W to e_w, so 1^WPrevLink in e_w leads back to W.
        Now it must lead to e_r.  Likewise, e_r needs to lead back to e_w
        with the opposing link, which is link[WPrevLink] */

     theEmbedding->G[e_w].link[1^WPrevLink] = e_r;
     theEmbedding->G[e_r].link[WPrevLink] = e_w;

     /* Now we cross-link W's link[WPrevLink] and link[1^WPrevLink] in the
        edge record e_ext */

     theEmbedding->G[W].link[WPrevLink] = e_ext;
     theEmbedding->G[e_ext].link[1^WPrevLink] = W;

     /* Erase the entries in R, which a root copy that is no longer needed. */

     _InitGraphNode(theEmbedding, R);
}

/********************************************************************
 _MergeBicomps()
 Merges all biconnected components at the cut vertices indicated by
 entries on the stack.
 ********************************************************************/

void _MergeBicomps(graphP theEmbedding)
{
int  R, Rout, Z, ZPrevLink;
int  theList, DFSChild, RootId;
int  extFaceVertex;

     while (sp_NonEmpty(theEmbedding->theStack))
     {
         sp_Pop2(theEmbedding->theStack, R, Rout);
         sp_Pop2(theEmbedding->theStack, Z, ZPrevLink);

         /* The external faces of the bicomps containing R and Z will
            form two corners at Z.  One corner will become part of the
            internal face formed by adding the new back edge. The other
            corner will be the new external face corner at Z.
            We first want to update the links at Z to reflect this. */

         extFaceVertex = theEmbedding->extFace[R].link[1^Rout];
         theEmbedding->extFace[Z].link[ZPrevLink] = extFaceVertex;

         if (theEmbedding->extFace[extFaceVertex].link[0] == theEmbedding->extFace[extFaceVertex].link[1])
            theEmbedding->extFace[extFaceVertex].link[Rout ^ theEmbedding->extFace[extFaceVertex].inversionFlag] = Z;
         else
            theEmbedding->extFace[extFaceVertex].link[theEmbedding->extFace[extFaceVertex].link[0] == R ? 0 : 1] = Z;

         /* If the path used to enter Z is opposed to the path
            used to exit R, then we have to flip the bicomp
            rooted at R, which we signify by inverting R
            then setting the sign on its DFS child edge to
            indicate that its descendants must be flipped later */

         if (ZPrevLink == Rout)
         {
             if (theEmbedding->G[R].link[0] != theEmbedding->G[R].link[1])
                _InvertVertex(theEmbedding, R);
             _SetSignOfChildEdge(theEmbedding, R, -1);
             Rout = 1^ZPrevLink;
         }

         /* R is no longer pertinent to Z since we are about to
            merge R into Z, so we delete R from its pertinent
            bicomp list (Walkdown gets R from the head of the list). */

         RootId = R - theEmbedding->N;
         theList = theEmbedding->V[Z].pertinentBicompList;
         theList = LCDelete(theEmbedding->BicompLists, theList, RootId);
         theEmbedding->V[Z].pertinentBicompList = theList;

         /* As a result of the merge, the DFS child of Z must be removed
            from Z's SeparatedDFSChildList because the child has just
            been joined directly to Z, rather than being separated by a
            root copy. */

         DFSChild = R - theEmbedding->N;
         theList = theEmbedding->V[Z].separatedDFSChildList;
         theList = LCDelete(theEmbedding->DFSChildLists, theList, DFSChild);
         theEmbedding->V[Z].separatedDFSChildList = theList;

         /* Now we push R into Z, eliminating R */

         _MergeVertex(theEmbedding, Z, ZPrevLink, R);
     }
}

/********************************************************************
 _RecordPertinentChildBicomp()
 Adds a child bicomp root into a vertex's pertinent bicomp list.
 ********************************************************************/

void _RecordPertinentChildBicomp(graphP theEmbedding, int I, int RootVertex)
{
int  ParentCopy, DFSChild, RootId, BicompList;

     /* The endpoints of the bicomp's root edge are the BicompRoot and a DFS Child
        of the parent copy of the bicomp root.
        NOTE: The location of the root vertices is in the range N to 2N-1 at an
                offset indicated by the DFS child with which it is associated */

     DFSChild = RootId = RootVertex - theEmbedding->N;
     ParentCopy = theEmbedding->V[DFSChild].DFSParent;

     /* Get the BicompList of the parent copy vertex. */

     BicompList = theEmbedding->V[ParentCopy].pertinentBicompList;

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

     if (theEmbedding->V[DFSChild].Lowpoint < I)
          BicompList = LCAppend(theEmbedding->BicompLists, BicompList, RootId);
     else BicompList = LCPrepend(theEmbedding->BicompLists, BicompList, RootId);

     /* The head node of the parent copy vertex's bicomp list may have changed, so
        we assign the head of the modified list as the vertex's pertinent
        bicomp list */

     theEmbedding->V[ParentCopy].pertinentBicompList = BicompList;
}

/********************************************************************
 _GetPertinentChildBicomp()
 This function returns the root of a pertinent child bicomp for the
 given vertex, with preference for an internally active child bicomp
 if one exists for the vertex.
 Note that always choosing the first member of the pertinentBicompList
 does this because, when we stored them, all of the internally active
 bicomps were stored before the externally active ones.
 ********************************************************************/

#ifndef SPEED_MACROS

int _GetPertinentChildBicomp(graphP theEmbedding, int W)
{
int  RootId;

    /* If the bicomp list is empty, then we just return NIL */

    if ((RootId=theEmbedding->V[W].pertinentBicompList) == NIL)
        return NIL;

    /* Return the RootVertex, which is computed by adding N because we
        subtracted N before storing it in the bicomp list */

    return RootId + theEmbedding->N;
}

#else

#define _GetPertinentChildBicomp(theEmbedding, W) \
        (theEmbedding->V[W].pertinentBicompList==NIL ? NIL : theEmbedding->V[W].pertinentBicompList + theEmbedding->N)

#endif

/********************************************************************
 _WalkUp()
 I is the vertex currently being embedded
 W is the vertex adjacent to I via the forward arc of back edge J=(I,W)

 The Walkup determines the pertinent child bicomps that should be
 set up as a result of the need to embed edge (I, W).

 It does this by recording the pertinent child biconnected components of
 all cut vertices between W and the child of I that is a descendant of W.
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

void _WalkUp(graphP theEmbedding, int I, int W)
{
int  Zig, Zag, ZigPrevLink, ZagPrevLink;
int  N, R, ParentCopy;
int  nextVertex;

     /* Shorthand for N, due to frequent use */

     N = theEmbedding->N;

     /* Start at the vertex W and walk around the both sides of the external face
        of a bicomp until we get back to vertex I. */

     Zig = Zag = W;
     ZigPrevLink = 1;
     ZagPrevLink = 0;

     while (Zig != I)
     {
        /* A previous walk-up may have been this way already */

        if (theEmbedding->G[Zig].visited == I) break;
        if (theEmbedding->G[Zag].visited == I) break;

        /* Mark the current vertices as visited during the embedding of vertex I. */

        theEmbedding->G[Zig].visited = I;
        theEmbedding->G[Zag].visited = I;

        /* Determine whether either Zig or Zag has landed on a bicomp root */

        if (Zig >= N) R = Zig;
        else if (Zag >= N) R = Zag;
        else R = NIL;

        /* If we have a bicomp root, then we want to hop up to the parent copy and
            record a pertinent child bicomp (except that we could but don't need
            to record pertinent child bicomps when the parent copy is I). */

        if (R != NIL)
        {
            ParentCopy = theEmbedding->V[R-N].DFSParent;
            if (ParentCopy != I)
                _RecordPertinentChildBicomp(theEmbedding, I, R);
            Zig = Zag = ParentCopy;
            ZigPrevLink = 1;
            ZagPrevLink = 0;
        }

        /* If we did not encounter a bicomp root, then we continue traversing the
            external face in both directions. */

        else
        {
            nextVertex = theEmbedding->extFace[Zig].link[1^ZigPrevLink];
            ZigPrevLink = theEmbedding->extFace[nextVertex].link[0] == Zig ? 0 : 1;
            Zig = nextVertex;

            nextVertex = theEmbedding->extFace[Zag].link[1^ZagPrevLink];
            ZagPrevLink = theEmbedding->extFace[nextVertex].link[0] == Zag ? 0 : 1;
            Zag = nextVertex;
        }
     }
}

/********************************************************************
 _WalkDown()
 Consider a circular shape with small circles and squares along its perimeter.
 The small circle at the top the root vertex of the bicomp.  The other small
 circles represent internally active vertices, and the squares represent
 externally active vertices.  The root vertex is a root copy of I, the
 vertex currently being processed.

 The Walkup previously marked all vertices adjacent to I by setting their
 adjacentTo flags.  Basically, we want to walkdown both the link[0] and
 then the link[1] sides of the bicomp rooted at RootVertex, embedding edges
 between it and descendants of I with the adjacentTo flag set.  It is sometimes
 necessary to hop to child biconnected components in order to reach the desired
 vertices and, in such cases, the biconnected components are merged together
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
 ********************************************************************/

void _WalkDown(graphP theEmbedding, int I, int RootVertex)
{
int  W, WPrevLink, R, Rout, X, XPrevLink, Y, YPrevLink, RootSide, RootEdgeChild;

     RootEdgeChild = RootVertex - theEmbedding->N;

     sp_ClearStack(theEmbedding->theStack);

     for (RootSide = 0; RootSide < 2; RootSide++)
     {
         WPrevLink = 1^RootSide;

         W = theEmbedding->extFace[RootVertex].link[RootSide];

         while (W != RootVertex)
         {
             /* If the vertex is adjacent to vertex I, then merge bicomps at cut
                vertices on theStack and add the back edge.  This creates a
                proper face. Then, clear W's AdjacentTo flag so we don't add an
                edge to W if we visit it again later. */

             if (theEmbedding->V[W].adjacentTo != NIL)
             {
                 _MergeBicomps(theEmbedding);
                 _EmbedBackEdgeToDescendant(theEmbedding, RootSide, RootVertex, W, WPrevLink);
                 theEmbedding->V[W].adjacentTo = NIL;
             }

             /* If there is a pertinent child bicomp, then we need to push it onto the stack
                along with information about how we entered the cut vertex and how
                we exit the root copy to get to the next vertex. */

             if (theEmbedding->V[W].pertinentBicompList != NIL)
             {
                 sp_Push2(theEmbedding->theStack, W, WPrevLink);
                 R = _GetPertinentChildBicomp(theEmbedding, W);

                 /* Get next active vertices X and Y on ext. face paths emanating from R */

                 X = theEmbedding->extFace[R].link[0];
                 XPrevLink = theEmbedding->extFace[X].link[1]==R ? 1 : 0;
                 Y = theEmbedding->extFace[R].link[1];
                 YPrevLink = theEmbedding->extFace[Y].link[0]==R ? 0 : 1;

                 /* If this is a bicomp with only two ext. face vertices, then
                    it could be that the orientation of the non-root vertex
                    doesn't match the orientation of the root due to our relaxed
                    orientation method. */

                 if (X == Y && theEmbedding->extFace[X].inversionFlag)
                 {
                     XPrevLink = 0;
                     YPrevLink = 1;
                 }

                 /* Now we implement the Walkdown's simple path selection rules!
                    If either X or Y is internally active (pertinent but not
                    externally active), then we pick it first.  Otherwise,
                    we choose a pertinent vertex. If neither are pertinent,
                    then we pick a vertex since the next iteration of the
                    loop will terminate on that vertex with a non-empty stack. */

                 if (_VertexActiveStatus(theEmbedding, X, I) == VAS_INTERNAL)
                      W = X;
                 else if (_VertexActiveStatus(theEmbedding, Y, I) == VAS_INTERNAL)
                      W = Y;
                 else if (PERTINENT(theEmbedding, X, I))
                      W = X;
                 else W = Y;

                 WPrevLink = W == X ? XPrevLink : YPrevLink;

                 Rout = W == X ? 0 : 1;
                 sp_Push2(theEmbedding->theStack, R, Rout);
             }

             /* Skip inactive vertices, which will be short-circuited
                later by our fast external face linking method (once
                upon a time, we added false edges called short-circuit
                edges to eliminate inactive vertices, but the extFace
                links can do the same job and also give us the ability
                to more quickly test planarity without creating an embedding). */

             else if (_VertexActiveStatus(theEmbedding, W, I) == VAS_INACTIVE)
             {
                 X = theEmbedding->extFace[W].link[1^WPrevLink];
                 WPrevLink = theEmbedding->extFace[X].link[0] == W ? 0 : 1;
                 W = X;
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

         if (sp_IsEmpty(theEmbedding->theStack))
         {
            theEmbedding->extFace[RootVertex].link[RootSide] = W;
            theEmbedding->extFace[W].link[WPrevLink] = RootVertex;

            /* If the bicomp is reduced to having only two external face vertices
                (the root and W), then we need to record whether the orientation
                of W is inverted relative to the root.  This is used later when a
                future Walkdown descends to and merges the bicomp containing W.
                Going from the root to W, we only get the correct WPrevLink if
                we know whether or not W is inverted.
                NOTE: Prior code based on short-circuit edges did not have this problem
                    because the root and W would be joined by two separate short-circuit
                    edges, so G[W].link[0] != G[W].link[1].
                NOTE: We clear the flag because it may have been set in W if W
                    previously became part of a bicomp with only two ext. face
                    vertices, but then was flipped and merged into a larger bicomp
                    that is now again becoming a bicomp with only two ext. face vertices. */

            if (theEmbedding->extFace[W].link[0] == theEmbedding->extFace[W].link[1] &&
                WPrevLink == RootSide)
                 theEmbedding->extFace[W].inversionFlag = 1;
            else theEmbedding->extFace[W].inversionFlag = 0;
         }

         /* If the stack is non-empty, then we had a non-planarity condition,
            so we stop.  If we got back around to the root, then all edges
            are embedded, so we stop. */

         if (sp_NonEmpty(theEmbedding->theStack) || W == RootVertex)
             break;
     }
}


/********************************************************************
 gp_Embed()

 return OK if the embedding was successfully created.

        NOTOK on internal failure

        NONPLANAR if the planar embedding couldn't be created due to
                crossing edges, and theGraph receives a Kuratowski
                subgraph (a.k.a. planar obstruction).
 ********************************************************************/

int gp_Embed(graphP theGraph, int embedFlags)
{
int N, I, J, W, child, RetVal;

    /* Basic parameter checks */

    if (theGraph==NULL) return NOTOK;
    if (embedFlags != 0 && embedFlags != EMBEDFLAGS_PLANAR) return NOTOK;

    /* A little shorthand for the size of the graph */

    N = theGraph->N;

    /* Preprocessing */

    theGraph->embedFlags = embedFlags;

    if (gp_CreateDFSTree(theGraph) != OK)
        return NOTOK;

    if (!(theGraph->internalFlags & FLAGS_SORTEDBYDFI))
        if (gp_SortVertices(theGraph) != OK)
            return NOTOK;

    gp_LowpointAndLeastAncestor(theGraph);

    _CreateSortedSeparatedDFSChildLists(theGraph);
    _CreateFwdArcLists(theGraph);

    _CreateDFSTreeEmbedding(theGraph);

    /* In reverse DFI order, process each vertex by embedding its
         the 'back edges' from the vertex to its DFS descendants. */

    _FillVisitedFlags(theGraph, N);

    for (I = theGraph->N-1; I >= 0; I--)
    {
          RetVal = OK;

          /* Do the Walkup for each cycle edge from I to a DFS descendant W. */

          J = theGraph->V[I].fwdArcList;
          while (J != NIL)
          {
              W = theGraph->G[J].v;
              theGraph->V[W].adjacentTo = J;
              _WalkUp(theGraph, I, W);

              J = theGraph->G[J].link[0];
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
                  _WalkDown(theGraph, I, child + N);

              child = LCGetNext(theGraph->DFSChildLists,
                                theGraph->V[I].separatedDFSChildList, child);
          }

          /* If all Walkdown calls succeed, but they don't embed all of the
                forward edges, then the graph is non-planar. */

          if (theGraph->V[I].fwdArcList != NIL)
          {
              RetVal = NONPLANAR;
              break;
          }
    }

    /* Post-process the embedding structure if an embedding was
        created to give a consistent orientation to all vertices. */

    if (RetVal == OK)
    {
        _OrientVerticesInEmbedding(theGraph);
        _JoinBicomps(theGraph);
    }

    else if (RetVal == NONPLANAR)
        _IsolateKuratowskiSubgraph(theGraph, I);

    /* Return indication of whether the graph is planar or non-planar. */

    return RetVal;
}

/********************************************************************
 _OrientVerticesInEmbedding()

 Each vertex will then have an orientation, either clockwise or
 counterclockwise.  All vertices in each bicomp need to have the
 same orientation.
 ********************************************************************/

void _OrientVerticesInEmbedding(graphP theEmbedding)
{
int  R, N;

     N = theEmbedding->N;

     sp_ClearStack(theEmbedding->theStack);

/* Run the array of root copy vertices.  For each that is not defunct
        (i.e. has not been merged during embed), we orient the vertices
        in the bicomp for which it is the root vertex. */

     for (R = N; R < 2*N; R++)
          if (theEmbedding->G[R].link[0] != NIL)
              _OrientVerticesInBicomp(theEmbedding, R, 0);
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
 ********************************************************************/

void _OrientVerticesInBicomp(graphP theEmbedding, int BicompRoot, int PreserveSigns)
{
int  V, J, curSign;

     sp_Push2(theEmbedding->theStack, BicompRoot, 1);

     while (sp_NonEmpty(theEmbedding->theStack))
     {
         /* Pop a vertex to orient */
         sp_Pop2(theEmbedding->theStack, V, curSign);

         /* Invert the vertex if the sign is -1 */
         if (curSign == -1)
             _InvertVertex(theEmbedding, V);

         /* Push the vertex's DFS children that are in the bicomp */
         J = theEmbedding->G[V].link[0];
         while (J >= 2*theEmbedding->N)
         {
             if (theEmbedding->G[J].type == EDGE_DFSCHILD)
             {
                 sp_Push2(theEmbedding->theStack, theEmbedding->G[J].v,
                                                  curSign * theEmbedding->G[J].sign);

                 if (!PreserveSigns)
                    theEmbedding->G[J].sign = 1;
             }

             J = theEmbedding->G[J].link[0];
         }
     }
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

int  _JoinBicomps(graphP theEmbedding)
{
int  R, N;

     for (R=N=theEmbedding->N; R < 2*N; R++)
          if (theEmbedding->G[R].link[0] != NIL)
              _MergeVertex(theEmbedding, theEmbedding->V[R-N].DFSParent, 0, R);

     return OK;
}
