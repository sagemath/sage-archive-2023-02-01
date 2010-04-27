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

#define GRAPHTEST_C

#include "graph.h"
#include "stack.h"

/* Private function declarations */

int  _TestPath(graphP theGraph, int U, int V);
int  _TryPath(graphP theGraph, int J, int V);
void _MarkPath(graphP theGraph, int J);
int  _TestSubgraph(graphP theSubgraph, graphP theGraph);

int  _CheckEmbeddingIntegrity(graphP theGraph, graphP origGraph);
int  _CheckEmbeddingFacialIntegrity(graphP theGraph);
int  _CheckObstructionIntegrity(graphP theGraph, graphP origGraph);

int  _CheckKuratowskiSubgraphIntegrity(graphP theGraph);
int  _CheckOuterplanarObstructionIntegrity(graphP theGraph);

int _CheckAllVerticesOnExternalFace(graphP theGraph);
void _MarkExternalFaceVertices(graphP theGraph, int startVertex);

/********************************************************************
 gp_TestEmbedResultIntegrity()

  This function tests the integrity of the graph result returned
  from gp_Embed().

  The caller of gp_Embed() does not have to save the original graph
  because, for efficiency, gp_Embed() operates on the input graph.
  However, to test the integrity of the result relative to the input,
  a copy of the input graph is required.

  Modules that extend/alter the behavior of gp_Embed() beyond the
  core planarity embedder and planarity obstruction isolator should
  also provide overriding integrity test routines appropriate to the
  extension algorithm.

  The main method first calls gp_SortVertices on theGraph, if the
  origGraph is not in DFI order (the common case).  Therefore,
  extension integrity tests can count on a consistent numbering
  between theGraph and the origGraph, either DFI order or pre-DFS
  order if that is the state of the origGraph.

  After all tests, the main method ensures theGraph is restored to
  DFI order by invoking gp_SortVertices if needed, thus ensuring
  that theGraph has the documented post-condition of gp_Embed().

  For an embedResult of OK, fpCheckEmbeddingIntegrity is invoked.
  The core planarity implementation does a face walk of all faces
  of the embedding.  It ensures that all edges were used in the face
  walk and that the right number of faces exist for the number of
  vertices and edges. Also, we ensure that all adjacencies expressed
  in the original graph still exist in the result graph.

  For an embedResult of NONEMBEDDABLE, fpCheckObstructionIntegrity
  is invoked.  The core planarity algorithm checks that the result
  graph is homeomorphic to K5 or K3,3 and that it is in fact a
  subgraph of the input graph.  Other algorithms use overloads to
  make appropriate checks.

  Returns NOTOK on integrity check failure or embedResult of NOTOK
          OK for successful integrity check of OK embedResult
          NONEMBEDDABLE for successful integrity check of an
                        embedResult of NONEMBEDDABLE
 ********************************************************************/

int gp_TestEmbedResultIntegrity(graphP theGraph, graphP origGraph, int embedResult)
{
int RetVal = embedResult;

    if (theGraph == NULL || origGraph == NULL)
        return NOTOK;

    if (embedResult == OK)
    {
        RetVal = theGraph->functions.fpCheckEmbeddingIntegrity(theGraph, origGraph);
    }
    else if (embedResult == NONEMBEDDABLE)
    {
        RetVal = theGraph->functions.fpCheckObstructionIntegrity(theGraph, origGraph);
    }

    if (RetVal == OK)
    	RetVal = embedResult;

    return RetVal;
}

/********************************************************************
 _CheckEmbeddingIntegrity()

  The core planarity implementation does a face walk of all faces
  of the embedding.  It ensures that all edges were used in the face
  walk and that the right number of faces exist for the number of
  vertices and edges. Also, we ensure that all adjacencies expressed
  in the original graph still exist in the result graph, accounting
  for the fact that the result graph is sorted by DFI, but the input
  may or may not be sorted by DFI.

  returns OK if all integrity tests passed, NOTOK otherwise
 ********************************************************************/

int _CheckEmbeddingIntegrity(graphP theGraph, graphP origGraph)
{
    if (theGraph == NULL || origGraph == NULL)
        return NOTOK;

    if (_TestSubgraph(theGraph, origGraph) != TRUE)
        return NOTOK;

    if (_TestSubgraph(origGraph, theGraph) != TRUE)
        return NOTOK;

    if (_CheckEmbeddingFacialIntegrity(theGraph) != OK)
        return NOTOK;

    if (theGraph->embedFlags == EMBEDFLAGS_OUTERPLANAR)
    {
        if (_CheckAllVerticesOnExternalFace(theGraph) != OK)
            return NOTOK;
    }

    return OK;
}

/********************************************************************
 _CheckEmbeddingFacialIntegrity()

 This function traverses all faces of a graph structure containing
 the planar embedding that results from gp_Embed().  The algorithm
 begins by placing all of the graph's arcs onto a stack and marking
 all of them as unvisited.  For each arc popped, if it is visited,
 it is immediately discarded and the next arc is popped.  Popping an
 unvisited arc J begins a face traversal.  We move to the true twin
 arc K of J, and obtain its successor arc L.  This amounts to
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
 NOTE:  In disconnected graphs, the face reader counts the external
        face of each connected component.  So, we adjust the face
        count by subtracting one for each component, then we add one
        to count the external face shared by all components.
 ********************************************************************/

int  _CheckEmbeddingFacialIntegrity(graphP theGraph)
{
stackP theStack = theGraph->theStack;
int I, e, J, JTwin, K, L, NumFaces, connectedComponents;

     if (theGraph == NULL)
         return NOTOK;

/* The stack need only contain 2M entries, one for each edge record. With
        max M at 3N, this amounts to 6N integers of space.  The embedding
        structure already contains this stack, so we just make sure it
        starts out empty. */

     sp_ClearStack(theStack);

/* Push all arcs and set them to unvisited */

     for (e=0, J=theGraph->edgeOffset; e < theGraph->M + sp_GetCurrentSize(theGraph->edgeHoles); e++, J+=2)
     {
          if (theGraph->G[J].v == NIL)
              continue;

          sp_Push(theStack, J);
          theGraph->G[J].visited = 0;
          JTwin = gp_GetTwinArc(theGraph, J);
          sp_Push(theStack, JTwin);
          theGraph->G[JTwin].visited = 0;
     }

     // There are M edges, so we better have pushed 2M arcs just now
     // i.e. testing that the continue above skipped only edge holes
     if (sp_GetCurrentSize(theStack) != 2*theGraph->M)
    	 return NOTOK;


/* Read faces until every arc is used */

     NumFaces = 0;
     while (sp_NonEmpty(theStack))
     {
            /* Get an arc; if it has already been used by a face, then
                don't use it to traverse a new face */
            sp_Pop(theStack, J);
            if (theGraph->G[J].visited) continue;

            L = NIL;
            JTwin = J;
            while (L != J)
            {
                K = gp_GetTwinArc(theGraph, JTwin);
                L = gp_GetNextArcCircular(theGraph, K);
                if (theGraph->G[L].visited)
                    return NOTOK;
                theGraph->G[L].visited++;
                JTwin = L;
            }
            NumFaces++;
     }

/* Count the external face once rather than once per connected component;
    each connected component is detected by the fact that it has no
    DFS parent, except in the case of isolated vertices, no face was counted
    so we do not subtract one. */

     for (I=connectedComponents=0; I < theGraph->N; I++)
          if (theGraph->V[I].DFSParent == NIL)
          {
              if (gp_GetVertexDegree(theGraph, I) > 0)
                  NumFaces--;
              connectedComponents++;
          }

     NumFaces++;

/* Test number of faces using the extended Euler's formula.
     For connected components, Euler's formula is f=m-n+2, but
     for disconnected graphs it is extended to f=m-n+1+c where
     c is the number of connected components.*/

     return NumFaces == theGraph->M - theGraph->N + 1 + connectedComponents
            ? OK : NOTOK;
}

/********************************************************************
 _CheckAllVerticesOnExternalFace()

  Determines whether or not any vertices have been embedded within
  the bounding cycle of the external face.
  The input graph may be disconnected, so this routine walks the
  external face starting at each vertex with no DFSParent.

  return OK if all vertices visited on external face walks, NOTOK otherwise
 ********************************************************************/

int _CheckAllVerticesOnExternalFace(graphP theGraph)
{
    int I;

    // Mark all vertices unvisited
    for (I=0; I < theGraph->N; I++)
        theGraph->G[I].visited = 0;

    // For each connected component, walk its external face and
    // mark the vertices as visited
    for (I=0; I < theGraph->N; I++)
    {
         if (theGraph->V[I].DFSParent == NIL)
        	 _MarkExternalFaceVertices(theGraph, I);
    }

    // If any vertex is unvisited, then the embedding is not an outerplanar
    // embedding, so we return NOTOK
    for (I=0; I < theGraph->N; I++)
        if (!theGraph->G[I].visited)
        {
            return NOTOK;
        }

    // All vertices were found on external faces of the connected components
    // so the embedding is an outerplanar embedding and we return OK
    return OK;
}

/********************************************************************
 _MarkExternalFaceVertices()

  Walks the external face of the connected component containing the
  start vertex, and marks the visited flag of all vertices found.
  The start vertex is assumed to be on the external face.
  This method assumed the embedding integrity has already been
  verified to be correct.
  This method correctly handles components that have cut vertices,
  i.e. it does not assume that the outer face is a simple cycle;
  it only assumes that all vertices are reachable by walking a
  single face that starts with startVertex.
 ********************************************************************/

void _MarkExternalFaceVertices(graphP theGraph, int startVertex)
{
    int nextVertex = startVertex;
    int Jout = gp_GetFirstArc(theGraph, nextVertex);
    int Jin;

    // Handle the case of an isolated vertex
    if (!gp_IsArc(theGraph, Jout))
    {
    	theGraph->G[startVertex].visited = 1;
    	return;
    }

    // Process a non-trivial connected component
    do {
        theGraph->G[nextVertex].visited = 1;

        // The arc out of the vertex just visited points to the next vertex
        nextVertex = theGraph->G[Jout].v;

        // Arc used to enter the next vertex is needed so we can get the
        // next edge in rotation order.
        // Note: for bicomps, first and last arcs of all external face vertices
        //       indicate the edges that hold them to the external face
        //       But _JoinBicomps() has already occurred, so cut vertices
        //       will have external face edges other than the first and last arcs
        //       Hence we need this more sophisticated traversal method
        Jin = gp_GetTwinArc(theGraph, Jout);

        // Now we get the next arc in rotation order as the new arc out to the
        // vertex after nextVertex.  This sets us up for the next iteration.
        // Note: We cannot simply follow the chain of nextVertex first arcs
        //       as we started out doing at the top of this method.  This is
        //       because we are no longer dealing with bicomps only.
        //       Since _JoinBicomps() has already been invoked, there may now
        //       be cut vertices on the external face whose adjacency lists
        //       contain external face arcs in positions other than the first and
        //       and last arcs.  We will visit those vertices multiple times,
        //       which is OK (just that we have to explain why we're calculating
        //       jout in this way).
        Jout = gp_GetNextArcCircular(theGraph, Jin);

        // Now things get really interesting.  The DFS root (startVertex) may
        // itself be a cut vertex to which multiple bicomps have been joined.
        // So we cannot simply stop when the external face walk gets back to
        // startVertex.  We must actually get back to startVertex using its
        // last arc.  This ensures that we've looped down into all the DFS
        // subtrees rooted at startVertex and walked their external faces.

        // Since we started the whole external face walk with the first arc
        // of startVertex, we need to proceed until we reenter startVertex
        // using its last arc.

    } while (Jin != gp_GetLastArc(theGraph, startVertex));
}


/********************************************************************
 _CheckObstructionIntegrity()

  Returns OK if theGraph is a subgraph of origGraph and it contains
    an allowed homeomorph, and NOTOK otherwise.

  For core planarity, the allowed homeomorphs are K_5 or K_{3,3}

  Extension modules may overload this method to implement different
  tests.  For example, K_{3,3} search allows only K_{3,3} and
  outerplanarity allows only K_4 or K_{2,3}.
 ********************************************************************/

int _CheckObstructionIntegrity(graphP theGraph, graphP origGraph)
{
    if (theGraph == NULL || origGraph == NULL)
        return NOTOK;

    if (_TestSubgraph(theGraph, origGraph) != TRUE)
    {
        return NOTOK;
    }

    if (theGraph->embedFlags == EMBEDFLAGS_PLANAR)
        return _CheckKuratowskiSubgraphIntegrity(theGraph);

    else if (theGraph->embedFlags == EMBEDFLAGS_OUTERPLANAR)
        return _CheckOuterplanarObstructionIntegrity(theGraph);

    return NOTOK;
}

/********************************************************************
 _getImageVertices()

 Count the number of vertices of each degree and find the locations of
 the image vertices (also sometimes called the corners of the obstruction).
 An image vertex is a vertex of degree three or higher because degree
 2 vertices are generally internal to the paths between the image
 vertices.

 The notable exception is K_{2,3}, an obstruction to outerplanarity.
 This routine does not know the obstruction it is looking for, so the
 caller must decide whether there are any degree 2 vertices that should
 be added to imageVerts.

 Return NOTOK if any vertex of degree 1 or higher than the max is found
        NOTOK if more than the max number of image vertices is found.
        Return OK otherwise.
 ********************************************************************/

int  _getImageVertices(graphP theGraph, int *degrees, int maxDegree,
                       int *imageVerts, int maxNumImageVerts)
{
int I, imageVertPos, degree;

     for (I = 0; I <= maxDegree; I++)
          degrees[I] = 0;

     for (I = 0; I < maxNumImageVerts; I++)
          imageVerts[I] = NIL;

     imageVertPos = 0;

     for (I = 0; I < theGraph->N; I++)
     {
          degree = gp_GetVertexDegree(theGraph, I);
          if (degree == 1)
              return NOTOK;
          if (degree > maxDegree)
              return NOTOK;

          degrees[degree]++;

          if (imageVertPos < maxNumImageVerts && degree > 2)
              imageVerts[imageVertPos++] = I;
          else if (degree > 2)
              return NOTOK;
     }

     return OK;
}

/********************************************************************
 _TestForCompleteGraphObstruction()

 theGraph - the graph to test
 numVerts - the number of image vertices (corners) of the complete
            graph homeomorph being tested for (e.g. 5 for a K5)
 degrees -  array of counts of the number of vertices of each degree
            given by the array index.  Only valid up to numVerts-1
 imageVerts - the vertices of degree numVerts-1

 This routine tests whether theGraph is a K_{numVerts} homeomorph for
 numVerts >= 4.

 returns FALSE if numVerts < 4,
               if theGraph has other than numVerts image vertices
               if theGraph contains other than degree 2 vertices plus
                  the image vertices
               if any pair of image vertices lacks a connecting path
               if any degree two vertices are not in the connecting paths
         TRUE  otherwise
 ********************************************************************/

int _TestForCompleteGraphObstruction(graphP theGraph, int numVerts,
                                     int *degrees, int *imageVerts)
{
    int I, J;

    // We need to make sure we have numVerts vertices of degree numVerts-1
    // For example, if numVerts==5, then we're looking for a K5, so we
    // need to have degrees[4] == 5 (5 vertices of degree 4)
    if (degrees[numVerts-1] != numVerts)
        return FALSE;

    // All vertices need to be degree 0, degree 2 or degree numVerts-1
    if (degrees[0]+degrees[2]+degrees[numVerts-1] != theGraph->N)
        return FALSE;

    // We clear all the vertex visited flags
    for (I = 0; I < theGraph->N; I++)
        theGraph->G[I].visited = 0;

    // For each pair of image vertices, we test that there is a path
    // between the two vertices.  If so, the visited flags of the
    // internal vertices along the path are marked
    //
    for (I = 0; I < numVerts; I++)
        for (J = 0; J < numVerts; J++)
           if (I != J)
               if (_TestPath(theGraph, imageVerts[I],
                                       imageVerts[J]) != TRUE)
                   return FALSE;

    // The visited flags should have marked only degree two vertices,
    // so for every marked vertex, we subtract one from the count of
    // the degree two vertices.
    for (I = 0; I < theGraph->N; I++)
        if (theGraph->G[I].visited)
            degrees[2]--;

    /* If every degree 2 vertex is used in a path between image
        vertices, then there are no extra pieces of the graph
        in theGraph.  Specifically, the prior tests identify
        a K_5 and ensure that nothing else could exist in the
        graph except extra degree 2 vertices, which must be
        joined in a cycle so that all are degree 2. */

    return degrees[2] == 0 ? TRUE : FALSE;
}

/********************************************************************
 _TestForK33GraphObstruction()

 theGraph - the graph to test
 degrees -  array of counts of the number of vertices of each degree
            given by the array index.  Only valid up to numVerts-1
 imageVerts - the degree 3 vertices of the K3,3 homeomorph

 This routine tests whether theGraph is a K_{3,3} homeomorph.

 returns TRUE if so, FALSE if not
 ********************************************************************/

int _TestForK33GraphObstruction(graphP theGraph, int *degrees, int *imageVerts)
{
int  I, imageVertPos, temp, success;

	if (degrees[4] != 0)
		return FALSE;

	if (degrees[3] != 6)
    	 return FALSE;

     /* Partition the six image vertices into two sets of 3
            (or report failure) */

     for (imageVertPos = 3; imageVertPos < 6; imageVertPos++)
     {
          I = success = 0;
          do {
             if (_TestPath(theGraph, imageVerts[imageVertPos], imageVerts[0]) == TRUE)
             {
                 success = TRUE;
                 break;
             }

             I++;
             temp = imageVerts[I];
             imageVerts[I] = imageVerts[imageVertPos];
             imageVerts[imageVertPos] = temp;
          }  while (I < 3);

          if (!success)
              return FALSE;
     }

     /* Now test the paths between each of the first three vertices and
            each of the last three vertices */

     for (I = 0; I < theGraph->N; I++)
          theGraph->G[I].visited = 0;

     for (imageVertPos=0; imageVertPos<3; imageVertPos++)
          for (I=3; I<6; I++)
               if (_TestPath(theGraph, imageVerts[imageVertPos],
                                       imageVerts[I]) != TRUE)
                   return FALSE;

     for (I = 0; I < theGraph->N; I++)
          if (theGraph->G[I].visited)
              degrees[2]--;

     /* If every degree 2 vertex is used in a path between image
            vertices, then there are no extra pieces of the graph
            in theGraph.  Specifically, the prior tests identify
            a K_{3,3} and ensure that nothing else could exist in the
            graph except extra degree 2 vertices, which must be
            joined in a cycle so that all are degree 2. */

     return degrees[2] == 0 ? TRUE : FALSE;
}

/********************************************************************
 _CheckKuratowskiSubgraphIntegrity()

 This function checks whether theGraph received as input contains
 either a K_5 or K_{3,3} homeomorph.

 RETURNS:   OK if theGraph contains a K_5 or K_{3,3} homeomorph,
            NOTOK otherwise

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

int  _CheckKuratowskiSubgraphIntegrity(graphP theGraph)
{
int  degrees[5], imageVerts[6];

     if (_getImageVertices(theGraph, degrees, 4, imageVerts, 6) != OK)
         return NOTOK;

     if (_TestForCompleteGraphObstruction(theGraph, 5, degrees, imageVerts) == TRUE)
     {
         return OK;
     }

     if (_TestForK33GraphObstruction(theGraph, degrees, imageVerts) == TRUE)
     {
         return OK;
     }

     return NOTOK;
}

/********************************************************************
 _TestForK23GraphObstruction()

 theGraph - the graph to test
 degrees -  array of counts of the number of vertices of each degree
            given by the array index.  Only valid up to numVerts-1
 imageVerts - the degree 3 vertices of the K2,3 homeomorph

 This routine tests whether theGraph is a K_{2,3} homeomorph.
 This routine operates over the results of _getImageVertices()

 returns TRUE if so, FALSE if not
 ********************************************************************/

int _TestForK23GraphObstruction(graphP theGraph, int *degrees, int *imageVerts)
{
int  I, J, imageVertPos;

     // This function operates over the imageVerts results produced by
     // getImageVertices, which only finds vertices of degree 3 or higher.
     // So, for a K2,3, there must be exactly two degree 3 vertices and
     // no degree 4 vertices.
     if (degrees[3] != 2)
         return FALSE;

     // For K_{2,3}, the three vertices of degree 2 were not
     // detected as image vertices because degree 2 vertices
     // are indistinguishable from the internal path vertices
     // between the image vertices.  So, here we acknowledge
     // that more image vertices need to be selected.
     imageVertPos = 2;

     // Assign the remaining three image vertices to be the
     // neighbors of the first degree 3 image vertex.
     // Ensure that each is distinct from the second
     // degree 3 image vertex. This must be the case because
     // the two degree 3 image vertices are in the same partition
     // and hence must not be adjacent.

     J = gp_GetFirstArc(theGraph, imageVerts[0]);
     while (gp_IsArc(theGraph, J))
     {
         imageVerts[imageVertPos] = theGraph->G[J].v;
         if (imageVerts[imageVertPos] == imageVerts[1])
             return FALSE;
         imageVertPos++;
         J = gp_GetNextArc(theGraph, J);
     }

     /* The paths from imageVerts[0] to each of the new degree 2
          image vertices are the edges we just traversed.
          Now test the paths between each of the degree 2 image
          vertices and imageVerts[1]. */

     for (I = 0; I < theGraph->N; I++)
          theGraph->G[I].visited = 0;

     for (imageVertPos=2; imageVertPos<5; imageVertPos++)
     {
          if (_TestPath(theGraph, imageVerts[imageVertPos],
                                  imageVerts[1]) != TRUE)
                   return FALSE;
          theGraph->G[imageVerts[imageVertPos]].visited = 1;
     }

     for (I = 0; I < theGraph->N; I++)
          if (theGraph->G[I].visited)
              degrees[2]--;

     /* If every degree 2 vertex is used in a path between the
          two degree 3 image vertices, then there are no extra
          pieces of the graph in theGraph.  Specifically, the
          prior tests identify a K_{2,3} and ensure that nothing
          else could exist in the graph... except extra degree 2
          vertices joined in a cycle. We return NOTOK in that case. */

     return degrees[2] == 0 ? TRUE : FALSE;
}

/********************************************************************
 _CheckOuterplanarObstructionIntegrity()

 This function checks whether theGraph received as input contains
 either a K_4 or K_{2,3} homeomorph.

 RETURNS:   OK if theGraph contains a K_4 or K_{2,3} homeomorph,
            NOTOK otherwise

 To be a K_4 homeomorph, there must be exactly 4 vertices of degree 3,
 which are called 'image' vertices, and all other vertices must be
 either degree 2 or degree 0. Furthermore, each of the image vertices
 must be able to reach all of the other image vertices by a single edge
 or a path of degree two vertices.

 To be a K_{2,3} homeomorph, there must be exactly 2 vertices of degree 3.
 All other vertices must be degree 2.  Furthermore, the two degree 3
 vertices must have three internally disjoint paths connecting them,
 and each path must contain at least two edges (i.e. at least one internal
 vertex).  The two degree 3 vertices are image vertices, and an internal
 vertex from each of the three paths contributes the remaining three
 image vertices.

 It is not necessary to check that the paths between the degree three
 vertices are distinct since if the paths had a common vertex, then the
 common vertex would not be degree 2.
 ********************************************************************/

int  _CheckOuterplanarObstructionIntegrity(graphP theGraph)
{
int  degrees[4], imageVerts[5];

     if (_getImageVertices(theGraph, degrees, 3, imageVerts, 5) != OK)
         return NOTOK;

     if (_TestForCompleteGraphObstruction(theGraph, 4, degrees, imageVerts) == TRUE)
     {
         return OK;
     }

     if (_TestForK23GraphObstruction(theGraph, degrees, imageVerts) == TRUE)
     {
         return OK;
     }

/* We get here only if we failed to recognize an outerplanarity
    obstruction, so we return failure */

     return NOTOK;
}

/********************************************************************
 _TestPath()

 This function determines whether there exists a path of degree two
 vertices between two given vertices.  The function marks each
 degree two vertex as visited.  It returns TRUE if it finds the
 path and FALSE otherwise.
 ********************************************************************/

int  _TestPath(graphP theGraph, int U, int V)
{
int  J;

     J = gp_GetFirstArc(theGraph, U);
     while (gp_IsArc(theGraph, J))
     {
         if (_TryPath(theGraph, J, V) == OK)
         {
             _MarkPath(theGraph, J);
             return TRUE;
         }

         J = gp_GetNextArc(theGraph, J);
     }

     return FALSE;
 }

/********************************************************************
 _TryPath()

 This function seeks a given path to a vertex V starting with a
 given edge out of a starting vertex U.  The path is allowed to
 contain zero or more degree two vertices, but we stop as soon as
 a vertex of degree higher than two is encountered.
 The function returns boolean true if that vertex is V, and
 boolean false otherwise.
 ********************************************************************/

int  _TryPath(graphP theGraph, int J, int V)
{
int  Jin, nextVertex;

     nextVertex = theGraph->G[J].v;

     // while nextVertex is strictly degree 2
     while (gp_IsArc(theGraph, gp_GetFirstArc(theGraph, nextVertex)) &&
    		gp_IsArc(theGraph, gp_GetLastArc(theGraph, nextVertex)) &&
    		gp_GetNextArc(theGraph, gp_GetFirstArc(theGraph, nextVertex)) == gp_GetLastArc(theGraph, nextVertex))
     {
         Jin = gp_GetTwinArc(theGraph, J);
         J = gp_GetFirstArc(theGraph, nextVertex);
         if (J == Jin)
             J = gp_GetLastArc(theGraph, nextVertex);

         nextVertex = theGraph->G[J].v;
     }

     return nextVertex == V ? TRUE : FALSE;
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
     // while nextVertex is strictly degree 2
     while (gp_IsArc(theGraph, gp_GetFirstArc(theGraph, nextVertex)) &&
    		gp_IsArc(theGraph, gp_GetLastArc(theGraph, nextVertex)) &&
    		gp_GetNextArc(theGraph, gp_GetFirstArc(theGraph, nextVertex)) == gp_GetLastArc(theGraph, nextVertex))
     {
         theGraph->G[nextVertex].visited = 1;

         Jin = gp_GetTwinArc(theGraph, J);
         J = gp_GetFirstArc(theGraph, nextVertex);
         if (J == Jin)
             J = gp_GetLastArc(theGraph, nextVertex);

         nextVertex = theGraph->G[J].v;
     }
}

/********************************************************************
 _TestSubgraph()
 Checks whether theSubgraph is in fact a subgraph of theGraph.
 For each vertex v in graph G and subgraph H, we iterate the adjacency
 list of H(v) and, for each neighbor w, we mark G(w).  Then, we
 iterate the adjacency list of G(v) and unmark each neighbor.  Then,
 we iterate the adjacency list of H(v) again to ensure that every
 neighbor w was unmarked.  If there exists a marked neighbor, then
 H(v) contains an incident edge that is not incident to G(v).

 Returns TRUE if theSubgraph contains only edges from theGraph,
         FALSE otherwise
 ********************************************************************/

int  _TestSubgraph(graphP theSubgraph, graphP theGraph)
{
int I, J;
int Result = TRUE;
int invokeSortOnGraph = FALSE;
int invokeSortOnSubgraph = FALSE;

    // If the graph is not sorted by DFI, but the alleged subgraph is,
    // then "unsort" the alleged subgraph so both have the same vertex order
    if (!(theGraph->internalFlags & FLAGS_SORTEDBYDFI) &&
         (theSubgraph->internalFlags & FLAGS_SORTEDBYDFI))
    {
        invokeSortOnSubgraph = TRUE;
        gp_SortVertices(theSubgraph);
    }

    // If the graph is not sorted by DFI, but the alleged subgraph is,
    // then "unsort" the alleged subgraph so both have the same vertex order
    if (!(theSubgraph->internalFlags & FLAGS_SORTEDBYDFI) &&
         (theGraph->internalFlags & FLAGS_SORTEDBYDFI))
    {
        invokeSortOnGraph = TRUE;
        gp_SortVertices(theGraph);
    }

/* We clear all visitation flags */

     for (I = 0; I < theGraph->N; I++)
          theGraph->G[I].visited = 0;

/* For each vertex... */

     for (I = 0; I < theSubgraph->N; I++)
     {
          /* For each neighbor w in the adjacency list of vertex I in the
                subgraph, set the visited flag in w in the graph */

          J = gp_GetFirstArc(theSubgraph, I);
          while (gp_IsArc(theSubgraph, J))
          {
        	  if (theSubgraph->G[J].v == NIL)
        	  {
        		  Result = FALSE;
        		  break;
        	  }
              theGraph->G[theSubgraph->G[J].v].visited = 1;
              J = gp_GetNextArc(theSubgraph, J);
          }

          if (Result != TRUE)
        	  break;

          /* For each neighbor w in the adjacency list of vertex I in the graph,
                clear the visited flag in w in the graph */

          J = gp_GetFirstArc(theGraph, I);
          while (gp_IsArc(theGraph, J))
          {
        	  if (theGraph->G[J].v == NIL)
        	  {
        		  Result = FALSE;
        		  break;
        	  }
              theGraph->G[theGraph->G[J].v].visited = 0;
              J = gp_GetNextArc(theGraph, J);
          }

          if (Result != TRUE)
        	  break;

          /* For each neighbor w in the adjacency list of vertex I in the
                subgraph, set the visited flag in w in the graph */

          J = gp_GetFirstArc(theSubgraph, I);
          while (gp_IsArc(theSubgraph, J))
          {
              if (theGraph->G[theSubgraph->G[J].v].visited)
              {
            	  Result = FALSE;
            	  break;
              }
              J = gp_GetNextArc(theSubgraph, J);
          }

          if (Result != TRUE)
        	  break;
     }

    // Restore the DFI sort order of either graph if it had to be reordered at the start
    if (invokeSortOnSubgraph)
        gp_SortVertices(theSubgraph);
    if (invokeSortOnGraph)
        gp_SortVertices(theGraph);

     return Result;
}
