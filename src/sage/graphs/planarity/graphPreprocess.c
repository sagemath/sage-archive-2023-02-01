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

#define GRAPHPREPROCESS_C

#include "graph.h"

/********************************************************************
 gp_CreateDFSTree
 Assigns Depth First Index (DFI) to each vertex.  Also records parent
 of each vertex in the DFS tree, and marks DFS tree edges that go from
 parent to child.  Forwared arc cycle edges are also distinguished from
 edges leading from a DFS tree descendant to an ancestor-- both DFS
 edges and back arcs.  This last type is marked with EDGE_BACK, which the
 embedder ignores.  Because the embedder works in reverse DFI order,
 only forward arcs to children are needed because only the DFS descendants
 of a vertex have already been embedded.
 ********************************************************************/

int  gp_CreateDFSTree(graphP theGraph)
{
stackP theStack = theGraph->theStack;
int DFI = 0, I, uparent, u, e, J;

     if (theGraph==NULL) return NOTOK;
     if (theGraph->internalFlags & FLAGS_DFSNUMBERED) return OK;

/* There are 2M edge records and for each we can push 2 integers,
        and M can be up to EDGE_LIMIT * N, so a stack of 4 * EDGE_LIMIT * N
        integers suffices. This is already in theGraph structure, so we
        make sure it's empty, then clear all visited flags in prep for the
        Depth first search. */

     sp_ClearStack(theStack);

     for (I=0; I < theGraph->N; I++)
          theGraph->G[I].visited = 0;

/* This outer loop causes the connected subgraphs of a disconnected
        graph to be numbered */

     for (I=0; I < theGraph->N && DFI < theGraph->N; I++)
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
                  theGraph->G[u].visited = 1;
                  theGraph->G[u].v = DFI++;
                  theGraph->V[u].DFSParent = uparent;
                  if (e != NIL)
                  {
                      theGraph->G[e].type = EDGE_DFSCHILD;
#ifdef ORDER_EDGES
                     // We want the child edges to be at the beginning
                     // of the adjacency list (link[0] side).

                     // Delete the edge from the list
                     theGraph->G[theGraph->G[e].link[0]].link[1] = theGraph->G[e].link[1];
                     theGraph->G[theGraph->G[e].link[1]].link[0] = theGraph->G[e].link[0];

                     // Tell the edge where it belongs now
                     theGraph->G[e].link[0] = theGraph->G[uparent].link[0];
                     theGraph->G[e].link[1] = uparent;

                     // Tell the rest of the list where the edge belongs
                     theGraph->G[uparent].link[0] = e;
                     theGraph->G[theGraph->G[e].link[0]].link[1] = e;
#endif
                  }

                  /* Push all neighbors */
                  J = theGraph->G[u].link[0];
                  while (J >= theGraph->N)
                  {
                      sp_Push2(theStack, u, J);
                      J = theGraph->G[J].link[0];
                  }
              }
              else
              {
             /* If the edge leads to a visited vertex, then it is either a
                forward cycle edge to a DFS descendant that has already been
                numbered by the depth first search, or it is a 'back' edge
                from a node to a DFS ancestor.  We distinguish between
                back cycle edges and the edge from a node to its DFS parent
                (which is a DFS tree edge not a cycle edge) by testing
                whether the twinArc is marked with EDGE_DFSCHILD. */

                 if (theGraph->G[uparent].v < theGraph->G[u].v)
                 {
                      theGraph->G[e].type = EDGE_FORWARD;
#ifdef ORDER_EDGES
                     // We want all of the forward edges to descendants to
                     // be at the end of the adjacency list (link[1]).
                     // The tree edge to the parent and the back edges to
                     // ancestors are in the middle, between the child edges
                     // and forward edges.

                     // Delete the edge from the list
                     theGraph->G[theGraph->G[e].link[0]].link[1] = theGraph->G[e].link[1];
                     theGraph->G[theGraph->G[e].link[1]].link[0] = theGraph->G[e].link[0];

                     // Tell the edge where it belongs now
                     theGraph->G[e].link[0] = uparent;
                     theGraph->G[e].link[1] = theGraph->G[uparent].link[1];

                     // Tell the rest of the list where the edge belongs
                     theGraph->G[uparent].link[1] = e;
                     theGraph->G[theGraph->G[e].link[1]].link[0] = e;
#endif
                 }
                 else if (theGraph->G[gp_GetTwinArc(theGraph, e)].type == EDGE_DFSCHILD)
                     theGraph->G[e].type = EDGE_DFSPARENT;
                 else
                     theGraph->G[e].type = EDGE_BACK;
              }
          }
     }

     theGraph->internalFlags |= FLAGS_DFSNUMBERED;
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
int  I, N, M, e, J, srcPos, dstPos;
vertexRec tempV;
graphNode tempG;

     if (theGraph == NULL) return NOTOK;
     if (!(theGraph->internalFlags&FLAGS_DFSNUMBERED))
         if (gp_CreateDFSTree(theGraph) != OK)
             return NOTOK;

/* Cache number of vertices and edges into local variables */

     N = theGraph->N;
     M = theGraph->M;

/* Change labels of edges from v to DFI(v)-- or vice versa
        Also, if any links go back to locations 0 to n-1, then they
        need to be changed because we are reordering the vertices */

     for (e=0, J=2*N; e < M; e++, J+=2)
     {
          theGraph->G[J].v = theGraph->G[theGraph->G[J].v].v;
          if (theGraph->G[J].link[0] < N)
              theGraph->G[J].link[0] = theGraph->G[theGraph->G[J].link[0]].v;
          if (theGraph->G[J].link[1] < N)
              theGraph->G[J].link[1] = theGraph->G[theGraph->G[J].link[1]].v;

          theGraph->G[J+1].v = theGraph->G[theGraph->G[J+1].v].v;
          if (theGraph->G[J+1].link[0] < N)
              theGraph->G[J+1].link[0] = theGraph->G[theGraph->G[J+1].link[0]].v;
          if (theGraph->G[J+1].link[1] < N)
              theGraph->G[J+1].link[1] = theGraph->G[theGraph->G[J+1].link[1]].v;
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

void gp_LowpointAndLeastAncestor(graphP theGraph)
{
stackP theStack = theGraph->theStack;
int I, u, uneighbor, J, L, leastAncestor;

     sp_ClearStack(theStack);

     for (I=0; I < theGraph->N; I++)
          theGraph->G[I].visited = 0;

/* This outer loop causes the connected subgraphs of a disconnected
        graph to be processed */

     for (I=0; I < theGraph->N; I++)
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
                  sp_Push(theStack, u);

                  /* Push DFS children */
                  J = theGraph->G[u].link[0];
                  while (J >= theGraph->N)
                  {
                      if (theGraph->G[J].type == EDGE_DFSCHILD)
                          sp_Push(theStack, theGraph->G[J].v);
#ifdef ORDER_EDGES
                      else break;
#endif
                      J = theGraph->G[J].link[0];
                  }
              }
              else
              {
                  /* Start with high values because we are doing a min function */
                  L = leastAncestor = u;

                  /* Compute L and leastAncestor */
                  J = theGraph->G[u].link[0];
                  while (J >= theGraph->N)
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
#ifdef ORDER_EDGES
                      else if (theGraph->G[J].type == EDGE_FORWARD)
                          break;
#endif

                      J = theGraph->G[J].link[0];
                  }

                  /* Assign leastAncestor and Lowpoint to the vertex */
                  theGraph->V[u].leastAncestor = leastAncestor;
                  theGraph->V[u].Lowpoint = leastAncestor < L ? leastAncestor : L;
              }
         }
     }
}
