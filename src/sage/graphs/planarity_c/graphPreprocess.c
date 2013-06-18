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

#define GRAPHPREPROCESS_C

#include "graph.h"

/********************************************************************
 gp_CreateDFSTree
 Assigns Depth First Index (DFI) to each vertex.  Also records parent
 of each vertex in the DFS tree, and marks DFS tree edges that go from
 parent to child.  Forward arc cycle edges are also distinguished from
 edges leading from a DFS tree descendant to an ancestor-- both DFS tree
 edges and back arcs.  The forward arcs are moved to the end of the
 adjacency list to make the set easier to find and process.
 ********************************************************************/

#include "platformTime.h"

int  gp_CreateDFSTree(graphP theGraph)
{
stackP theStack;
int N, DFI = 0, I, uparent, u, e, J;

#ifdef PROFILE
platform_time start, end;
platform_GetTime(start);
#endif

     if (theGraph==NULL) return NOTOK;
     if (theGraph->internalFlags & FLAGS_DFSNUMBERED) return OK;

     gp_LogLine("\ngraphPreprocess.c/gp_CreateDFSTree() start");

     N = theGraph->N;
     theStack  = theGraph->theStack;

/* There are 2M edge records (arcs) and for each we can push 2 integers,
        so a stack of 2 * arcCapacity integers suffices.
        This is already in theGraph structure, so we make sure it's empty,
        then clear all visited flags in prep for the Depth first search. */

     if (sp_GetCapacity(theStack) < 2*gp_GetArcCapacity(theGraph))
    	 return NOTOK;

     sp_ClearStack(theStack);

     for (I=0; I < N; I++)
          theGraph->G[I].visited = 0;

/* This outer loop causes the connected subgraphs of a disconnected
        graph to be numbered */

     for (I=0; I < N && DFI < N; I++)
     {
          if (theGraph->V[I].DFSParent != NIL)
              continue;

          sp_Push2(theStack, NIL, NIL);
          while (sp_NonEmpty(theStack))
          {
              sp_Pop2(theStack, uparent, e);
              u = uparent == NIL ? I : theGraph->G[e].v;

              if (!theGraph->G[u].visited)
              {
            	  gp_LogLine(gp_MakeLogStr3("V=%d, DFI=%d, Parent=%d", u, DFI, uparent));

            	  theGraph->G[u].visited = 1;
                  theGraph->G[u].v = DFI++;
                  theGraph->V[u].DFSParent = uparent;
                  if (e != NIL)
                  {
                      theGraph->G[e].type = EDGE_DFSCHILD;
                      theGraph->G[gp_GetTwinArc(theGraph, e)].type = EDGE_DFSPARENT;

                      // We want the child arcs to be at the beginning
                      // of the adjacency list.
                      gp_MoveArcToFirst(theGraph, uparent, e);
                  }

                  /* Push edges to all unvisited neighbors. These will be either
                        tree edges to children or forward arcs of back edges */

                  J = gp_GetFirstArc(theGraph, u);
                  while (gp_IsArc(theGraph, J))
                  {
                      if (!theGraph->G[theGraph->G[J].v].visited)
                          sp_Push2(theStack, u, J);
                      J = gp_GetNextArc(theGraph, J);
                  }
              }
              else
              {
                  // If the edge leads to a visited vertex, then it is
            	  // the forward arc of a back edge.
                  theGraph->G[e].type = EDGE_FORWARD;
                  theGraph->G[gp_GetTwinArc(theGraph, e)].type = EDGE_BACK;

                  // We want all of the forward edges to descendants to
                  // be at the end of the adjacency list.
                  // The tree edge to the parent and the back edges to ancestors
                  // are in the middle, between the child edges and forward edges.
                  gp_MoveArcToLast(theGraph, uparent, e);
              }
          }
     }

     gp_LogLine("graphPreprocess.c/gp_CreateDFSTree() end\n");

     theGraph->internalFlags |= FLAGS_DFSNUMBERED;

#ifdef PROFILE
platform_GetTime(end);
printf("DFS in %.3lf seconds.\n", platform_GetDuration(start,end));
#endif

     return OK;
}

/********************************************************************
 gp_SortVertices()
 Once depth first numbering has been applied to the graph, the v member
 of each vertex contains the DFI.  This routine can reorder the vertices
 in linear time so that they appear in ascending order by DFI.  Note
 that the field v is then used to store the original number of the
 vertex.
 Note that this function is not underscored (i.e. not private).  We
 export it because it can be called again at a later point to reorder
 the vertices back into the original numbering with the DFI values
 stored in the v fields (in other words this function is its own inverse).
 ********************************************************************/

int  gp_SortVertices(graphP theGraph)
{
     return theGraph->functions.fpSortVertices(theGraph);
}

int  _SortVertices(graphP theGraph)
{
int  I, N, M, e, J, srcPos, dstPos;
vertexRec tempV;
graphNode tempG;

#ifdef PROFILE
platform_time start, end;
platform_GetTime(start);
#endif

     if (theGraph == NULL) return NOTOK;
     if (!(theGraph->internalFlags&FLAGS_DFSNUMBERED))
         if (gp_CreateDFSTree(theGraph) != OK)
             return NOTOK;

/* Cache number of vertices and edges into local variables */

     N = theGraph->N;
     M = theGraph->M + sp_GetCurrentSize(theGraph->edgeHoles);

/* Change labels of edges from v to DFI(v)-- or vice versa
        Also, if any links go back to locations 0 to n-1, then they
        need to be changed because we are reordering the vertices */

     for (e=0, J=theGraph->edgeOffset; e < M; e++, J+=2)
     {
    	 if (theGraph->G[J].v != NIL)
    	 {
    		 theGraph->G[J].v = theGraph->G[theGraph->G[J].v].v;
    		 theGraph->G[J+1].v = theGraph->G[theGraph->G[J+1].v].v;
    	 }
     }

/* Convert DFSParent from v to DFI(v) or vice versa */

     for (I=0; I < N; I++)
          if (theGraph->V[I].DFSParent != NIL)
              theGraph->V[I].DFSParent = theGraph->G[theGraph->V[I].DFSParent].v;

/* Sort by 'v using constant time random access. Move each vertex to its
        destination 'v', and store its source location in 'v'. */

     /* First we clear the visitation flags.  We need these to help mark
        visited vertices because we change the 'v' field to be the source
        location, so we cannot use index==v as a test for whether the
        correct vertex is in location 'index'. */

     for (I=0; I < N; I++)
          theGraph->G[I].visited = 0;

     /* We visit each vertex location, skipping those marked as visited since
        we've already moved the correct vertex into that location. The
        inner loop swaps the vertex at location I into the correct position,
        G[I].v, marks that location as visited, then sets its 'v' field to
        be the location from whence we obtained the vertex record. */

     for (I=0; I < N; I++)
     {
          srcPos = I;
          while (!theGraph->G[I].visited)
          {
              dstPos = theGraph->G[I].v;

              tempG = theGraph->G[dstPos];
              tempV = theGraph->V[dstPos];
              theGraph->G[dstPos] = theGraph->G[I];
              theGraph->V[dstPos] = theGraph->V[I];
              theGraph->G[I] = tempG;
              theGraph->V[I] = tempV;

              theGraph->G[dstPos].visited = 1;
              theGraph->G[dstPos].v = srcPos;

              srcPos = dstPos;
          }
     }

/* Invert the bit that records the sort order of the graph */

     if (theGraph->internalFlags & FLAGS_SORTEDBYDFI)
          theGraph->internalFlags &= ~FLAGS_SORTEDBYDFI;
     else theGraph->internalFlags |= FLAGS_SORTEDBYDFI;

#ifdef PROFILE
platform_GetTime(end);
printf("SortVertices in %.3lf seconds.\n", platform_GetDuration(start,end));
#endif

     return OK;
}

/********************************************************************
 gp_LowpointAndLeastAncestor()
        leastAncestor: min(DFI of neighbors connected by backedge)
        Lowpoint: min(leastAncestor, Lowpoint of DFSChildren)

 This implementation requires that the vertices be sorted in DFI order
 (such that the edge records contain DFI values in their v fields).  An
 implementation could be made to run before sorting using the fact that
 the value G[G[e].v].v before sorting is equal to G[e].v after the sort.

 For computing Lowpoint, we must do a post-order traversal of the DFS tree.
 We push the root of the DFS tree, then we loop while the stack is not empty.
 We pop a vertex; if it is not marked, then we are on our way down the DFS
 tree, so we mark it and push it back on, followed by pushing its
 DFS children.  The next time we pop the node, all of its children
 will have been popped, marked+children pushed, and popped again.  On
 the second pop of the vertex, we can therefore compute the Lowpoint
 values based on the childrens' Lowpoints and the edges in the vertex's
 adjacency list.

 A stack of size N suffices because we push each vertex only once.
 ********************************************************************/

int  gp_LowpointAndLeastAncestor(graphP theGraph)
{
stackP theStack = theGraph->theStack;
int I, u, uneighbor, J, L, leastAncestor;
int totalVisited = 0;

#ifdef PROFILE
platform_time start, end;
platform_GetTime(start);
#endif

     sp_ClearStack(theStack);

     for (I=0; I < theGraph->N; I++)
          theGraph->G[I].visited = 0;

/* This outer loop causes the connected subgraphs of a disconnected
        graph to be processed */

     for (I=0; I < theGraph->N && totalVisited < theGraph->N; I++)
     {
          if (theGraph->G[I].visited)
              continue;

          sp_Push(theStack, I);
          while (sp_NonEmpty(theStack))
          {
              sp_Pop(theStack, u);
              if (!theGraph->G[u].visited)
              {
                  /* Mark u as visited, then push it back on the stack */
                  theGraph->G[u].visited = 1;
                  totalVisited++;
                  sp_Push(theStack, u);

                  /* Push DFS children */
                  J = gp_GetFirstArc(theGraph, u);
                  while (gp_IsArc(theGraph, J))
                  {
                      if (theGraph->G[J].type == EDGE_DFSCHILD)
                      {
                          sp_Push(theStack, theGraph->G[J].v);
                      }
                      else break;

                      J = gp_GetNextArc(theGraph, J);
                  }
              }
              else
              {
                  /* Start with high values because we are doing a min function */
                  L = leastAncestor = u;

                  /* Compute L and leastAncestor */
                  J = gp_GetFirstArc(theGraph, u);
                  while (gp_IsArc(theGraph, J))
                  {
                      uneighbor = theGraph->G[J].v;
                      if (theGraph->G[J].type == EDGE_DFSCHILD)
                      {
                          if (L > theGraph->V[uneighbor].Lowpoint)
                              L = theGraph->V[uneighbor].Lowpoint;
                      }
                      else if (theGraph->G[J].type == EDGE_BACK)
                      {
                          if (leastAncestor > uneighbor)
                              leastAncestor = uneighbor;
                      }
                      else if (theGraph->G[J].type == EDGE_FORWARD)
                          break;

                      J = gp_GetNextArc(theGraph, J);
                  }

                  /* Assign leastAncestor and Lowpoint to the vertex */
                  theGraph->V[u].leastAncestor = leastAncestor;
                  theGraph->V[u].Lowpoint = leastAncestor < L ? leastAncestor : L;
              }
         }
     }

#ifdef PROFILE
platform_GetTime(end);
printf("Lowpoint in %.3lf seconds.\n", platform_GetDuration(start,end));
#endif

     return OK;
}

