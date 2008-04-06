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

#define GRAPHTEST_C

#include "graph.h"
#include "stack.h"

/* Private function declarations */

int  _TestPath(graphP theGraph, int U, int V);
int  _TryPath(graphP theGraph, int J, int V);
void _MarkPath(graphP theGraph, int J);
int  _TestSubgraph(graphP theSubgraph, graphP theGraph);

/********************************************************************
 gp_CheckEmbeddingIntegrity()

 This function traverses all faces of a graph structure containing
 the planar embedding that results from gp_Embed().  The algorithm
 begins by placing all of the graph's arcs onto a stack and marking
 all of them as unvisited.  For each arc popped, if it is visited,
 it is immediately discarded and the next arc is popped.  Popping an
 unvisited arc J begins a face traversal.  We move to the true twin
 arc K of J, and obtain its link[0] successor arc L.  This amounts to
 always going clockwise or counterclockwise (depending on how the
 graph is drawn on the plane, or alternately whether one is above
 or below the plane).  This traversal continues until we make it
 back to J.  Each arc along the way is marked as visited.  Further,
 if L has been visited, then there is an error since an arc can only
 appear in one face (the twin arc appears in a separate face, which
 is traversed in the opposing direction).
 If this algorithm succeeds without double visiting any arcs, and it
 produces the correct face count according to Euler's formula, then
 the embedding has all vertices oriented the same way.
 NOTE:  Because a vertex is in its adj. list, if we go to a link[0]
        and it is a vertex, we simply take the vertex's link[0].
 NOTE:  In disconnected graphs, the face reader counts the external
        face of each connected component.  So, we adjust the face
        count by subtracting one for each component, then we add one
        to count the external face shared by all components.
 ********************************************************************/

int  gp_CheckEmbeddingIntegrity(graphP theEmbedding)
{
stackP theStack = theEmbedding->theStack;
int I, e, J, JTwin, K, L, NumFaces, connectedComponents;

     if (theEmbedding == NULL) return NOTOK;

/* The stack need only contain 2M entries, one for each edge record. With
        max M at 3N, this amounts to 6N integers of space.  The embedding
        structure already contains this stack, so we just make sure it
        starts out empty. */

     sp_ClearStack(theStack);

/* Push all arcs and set them to unvisited */

     for (e=0, J=2*theEmbedding->N; e < theEmbedding->M; e++, J+=2)
     {
          sp_Push(theStack, J);
          theEmbedding->G[J].visited = 0;
          JTwin = gp_GetTwinArc(theEmbedding, J);
          sp_Push(theStack, JTwin);
          theEmbedding->G[JTwin].visited = 0;
     }

/* Read faces until every arc is used */

     NumFaces = 0;
     while (sp_NonEmpty(theStack))
     {
            /* Get an arc; if it has already been used by a face, then
                don't use it to traverse a new face */
            sp_Pop(theStack, J);
            if (theEmbedding->G[J].visited) continue;

            L = NIL;
            JTwin = J;
            while (L != J)
            {
                K = gp_GetTwinArc(theEmbedding, JTwin);
                L = theEmbedding->G[K].link[0];
                if (L < 2*theEmbedding->N)
                    L = theEmbedding->G[L].link[0];
                if (theEmbedding->G[L].visited)
                    return NOTOK;
                theEmbedding->G[L].visited++;
                JTwin = L;
            }
            NumFaces++;
     }

/* Count the external face once rather than once per connected component;
    each connected component is detected by the fact that it has no
    DFS parent, except in the case of isolated vertices, no face was counted
    so we do not subtract one. */

     for (I=connectedComponents=0; I < theEmbedding->N; I++)
          if (theEmbedding->V[I].DFSParent == NIL)
          {
              if (gp_GetVertexDegree(theEmbedding, I) > 0)
                  NumFaces--;
              connectedComponents++;
          }

     NumFaces++;

/* Test number of faces using the extended Euler's formula.
     For connected components, Euler's formula is f=m-n+2, but
     for disconnected graphs it is extended to f=m-n+1+c where
     c is the number of connected components.*/

     return NumFaces == theEmbedding->M - theEmbedding->N + 1 + connectedComponents
            ? OK : NOTOK;
}

/********************************************************************
 gp_CheckKuratowskiSubgraphIntegrity()

 This function checks whether theGraph received as input contains
 either a K_5 or K_{3,3} homeomorph.

 To be a K_5 homeomorph, there must be exactly 5 vertices of degree 4,
 which are called 'image' vertices, and all other vertices must be
 either degree 2 or degree 0. Furthermore, each of the image vertices
 must be able to reach all of the other image vertices by a single edge
 or a path of degree two vertices.

 To be a K_{3,3} homeomorph, there must be exactly 6 vertices of degree 3,
 which are called 'image' vertices, and all other vertices must be either
 degree 2 or degree 0.  Furthermore, the image vertices must be connected
 by edges or paths of degree two vertices in the manner suggested by
 a K_{3,3}.  To test this, we select an image vertex U, and determine
 three image vertices X, Y and Z reachable from U by single edges or
 paths of degree 2 vertices.  Then, we check that the two remaining image
 vertices, V and W, can also reach X, Y and Z by single edges or paths of
 degree 2 vertices.

 It is not necessary to check that the paths between the image vertices
 are distinct since if the paths had a common vertex, then the common
 vertex would not be degree 2.
 ********************************************************************/

int  gp_CheckKuratowskiSubgraphIntegrity(graphP theGraph)
{
int  degrees[5], imageVerts[6];
int  I, degree, imageVertPos, temp, success;

/* Is theGraph even a subgraph of the original graph? */

//     if (_TestSubgraph(theGraph, origGraph) != OK)
//         return NOTOK;

/* Init the variables that count the number of verts of each degree
        and the location of each image vert. */

     for (I = 0; I < 5; I++)
          degrees[I] = 0;

     for (I = 0; I < 6; I++)
          imageVerts[I] = NIL;

     imageVertPos = 0;

/* Count the number of verts of each degree and find the locations of
        the image verts.  Quit if any vertex of degree higher than 4
        or equal to 1 is encountered. */

     for (I = 0; I < theGraph->N; I++)
     {
          degree = gp_GetVertexDegree(theGraph, I);
          if (degree == 1 || degree > 4) return NOTOK;
          degrees[degree]++;
          if (imageVertPos < 6 && degree > 2)
              imageVerts[imageVertPos++] = I;
     }

/* If there are five vertices of degree 4 (and no degree 3 vertices),
        then test for a K_5 */

     if (degrees[4] == 5 && degrees[3] == 0)
     {
         for (I = 0; I < theGraph->N; I++)
              theGraph->G[I].visited = 0;

         for (imageVertPos = 0; imageVertPos < 5; imageVertPos++)
              for (I = 0; I < 5; I++)
                   if (imageVertPos != I)
                       if (_TestPath(theGraph, imageVerts[imageVertPos],
                                               imageVerts[I]) != OK)
                           return NOTOK;

         for (I = 0; I < theGraph->N; I++)
              if (theGraph->G[I].visited)
                  degrees[2]--;

         /* If every degree 2 vertex is used in a path between image
                vertices, then there are no extra pieces of the graph
                in theGraph.  Specifically, the prior tests identify
                a K_5 and ensure that nothing else could exist in the
                graph except extra degree 2 vertices, which must be
                joined in a cycle so that all are degree 2. */

         return degrees[2] == 0 ? OK : NOTOK;
     }

/* Otherwise, if there are six vertices of degree 3 (and no degree 4 vertices),
        then test for a K_{3,3} */

     else if (degrees[3] == 6 && degrees[4] == 0)
     {
         /* Partition the six image vertices into two sets of 3
                (or report failure) */

         for (imageVertPos = 3; imageVertPos < 6; imageVertPos++)
         {
              I = success = 0;
              do {
                 if (_TestPath(theGraph, imageVerts[imageVertPos], imageVerts[0]) == OK)
                 {
                     success = 1;
                     break;
                 }

                 I++;
                 temp = imageVerts[I];
                 imageVerts[I] = imageVerts[imageVertPos];
                 imageVerts[imageVertPos] = temp;
              }  while (I < 3);

              if (!success) return NOTOK;
         }

         /* Now test the paths between each of the first three vertices and
                each of the last three vertices */

         for (I = 0; I < theGraph->N; I++)
              theGraph->G[I].visited = 0;

         for (imageVertPos=0; imageVertPos<3; imageVertPos++)
              for (I=3; I<6; I++)
                   if (_TestPath(theGraph, imageVerts[imageVertPos],
                                           imageVerts[I]) != OK)
                       return NOTOK;

         for (I = 0; I < theGraph->N; I++)
              if (theGraph->G[I].visited)
                  degrees[2]--;

         /* If every degree 2 vertex is used in a path between image
                vertices, then there are no extra pieces of the graph
                in theGraph.  Specifically, the prior tests identify
                a K_{3,3} and ensure that nothing else could exist in the
                graph except extra degree 2 vertices, which must be
                joined in a cycle so that all are degree 2. */

         return degrees[2] == 0 ? OK : NOTOK;
     }

/* We get here only if we failed to recognize a Kuratowski subgraph, so
        we return failure */

     return NOTOK;
}

/********************************************************************
 _TestPath()

 This function determines whether there exists a path of degree two
 vertices between two given vertices.  The function marks each
 degree two vertex as visited.  It returns NOTOK if it cannot find
 such a path.
 ********************************************************************/

int  _TestPath(graphP theGraph, int U, int V)
{
int  J;

     J = theGraph->G[U].link[0];

     while (J > theGraph->N)
     {
         if (_TryPath(theGraph, J, V) == OK)
         {
             _MarkPath(theGraph, J);
             return OK;
         }

         J = theGraph->G[J].link[0];
     }

     return NOTOK;
 }

/********************************************************************
 _TryPath()

 This function seeks a given path to a vertex V starting with a
 given edge out of a starting vertex U.  The path is allowed to
 contain zero or more degree two vertices, but we stop as soon as
 a vertex of degree higher than two is encountered.
 The function returns OK if that vertex is V, and NOTOK otherwise.
 ********************************************************************/

int  _TryPath(graphP theGraph, int J, int V)
{
int  Jin, nextVertex;

     nextVertex = theGraph->G[J].v;
     while (gp_GetVertexDegree(theGraph, nextVertex) == 2)
     {
         Jin = gp_GetTwinArc(theGraph, J);
         J = theGraph->G[nextVertex].link[0];
         if (J == Jin)
             J = theGraph->G[nextVertex].link[1];

         nextVertex = theGraph->G[J].v;
     }

     return nextVertex == V ? OK : NOTOK;
}

/********************************************************************
 _MarkPath()

 This function sets the visitation flag on all degree two vertices
 along a path to a vertex V that starts with a given edge out of
 a starting vertex U.
 ********************************************************************/

void _MarkPath(graphP theGraph, int J)
{
int  Jin, nextVertex;

     nextVertex = theGraph->G[J].v;
     while (gp_GetVertexDegree(theGraph, nextVertex) == 2)
     {
         theGraph->G[nextVertex].visited = 1;

         Jin = gp_GetTwinArc(theGraph, J);
         J = theGraph->G[nextVertex].link[0];
         if (J == Jin)
             J = theGraph->G[nextVertex].link[1];

         nextVertex = theGraph->G[J].v;
     }
}

/********************************************************************
 _TestSubgraph()
 Checks whether theSubgraph is in face a subgraph of theGraph.
 For each vertex v in graph G and subgraph H, we iterate the adjacency
 list of H(v) and, for each neighbor w, we mark G(w).  Then, we
 iterate the adjacency list of G(v) and unmark each neighbor.  Then,
 we iterate the adjacency list of H(v) again to ensure that every
 neighbor w was unmarked.  If there exists a marked neighbor, then
 H(v) contains an incident edge that is not incident to G(v).

 Returns OK if theSubgraph contains only edges from theGraph,
         NOTOK otherwise
 ********************************************************************/

int  _TestSubgraph(graphP theSubgraph, graphP theGraph)
{
int I, J;

/* We clear all visitation flags */

     for (I = 0; I < theGraph->N; I++)
          theGraph->G[I].visited = 0;

/* For each vertex... */

     for (I = 0; I < theSubgraph->N; I++)
     {
          /* For each neighbor w in the adjacency list of vertex I in the
                subgraph, set the visited flag in w in the graph */

          J = theSubgraph->G[I].link[0];
          while (J >= 2*theSubgraph->N)
          {
              theGraph->G[theSubgraph->G[J].v].visited = 1;
              J = theSubgraph->G[J].link[0];
          }

          /* For each neighbor w in the adjacency list of vertex I in the graph,
                clear the visited flag in w in the graph */

          J = theGraph->G[I].link[0];
          while (J >= 2*theGraph->N)
          {
              theGraph->G[theGraph->G[J].v].visited = 0;
              J = theGraph->G[J].link[0];
          }

          /* For each neighbor w in the adjacency list of vertex I in the
                subgraph, set the visited flag in w in the graph */

          J = theSubgraph->G[I].link[0];
          while (J >= 2*theSubgraph->N)
          {
              if (theGraph->G[theSubgraph->G[J].v].visited)
                  return NOTOK;
              J = theSubgraph->G[J].link[0];
          }
     }

     return OK;
}
