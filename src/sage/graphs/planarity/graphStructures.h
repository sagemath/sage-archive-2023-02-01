#ifndef GRAPHSTRUCTURE_H
#define GRAPHSTRUCTURE_H

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

#include <stdio.h>
#include "appconst.h"
#include "listcoll.h"
#include "stack.h"

#include "graphFunctionTable.h"
#include "graphExtensions.private.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The DEFAULT_EDGE_LIMIT expresses the initial setting for the arcCapacity
 * as a constant factor of N, the number of vertices. We allow 3N edges, but
 * this number can be safely set to a larger integer value.
 */

#define DEFAULT_EDGE_LIMIT      3

/* Simple integer selection macros */

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

#define MIN3(x, y, z) MIN(MIN((x), (y)), MIN((y), (z)))
#define MAX3(x, y, z) MAX(MAX((x), (y)), MAX((y), (z)))

/* Vertex activity categories */

#define VAS_INACTIVE    0
#define VAS_INTERNAL    1
#define VAS_EXTERNAL    2

/* Types:

   TYPE_UNKNOWN - initial assignment

   Edge types: (assigned by depth first search; used throughout algorithms)

   EDGE_DFSCHILD - the arc is an edge to a DFS child; these are embedded first
                        as singleton bicomps.
   EDGE_FORWARD - back edge directed from DFS ancestor to descendant
   EDGE_BACK - DFS tree edge _or_ back edge directed from descendant to
                ancestor.  Embedder ignores these because the ancestors of a
                vertex are only embedded after the vertex.
   EDGE_DFSPARENT - If the arc (u,v) is of type EDGE_DFSCHILD, then the
                        twin arc (v,u) is marked with EDGE_DFSPARENT

   Vertex types: (used when searching paths of interest in a non-planar graph)

   VERTEX_HIGH_RXW - On the external face path between vertices R and X
   VERTEX_LOW_RXW  - X or on the external face path between vertices X and W
   VERTEX_HIGH_RYW - On the external face path between vertices R and Y
   VERTEX_LOW_RYW  - Y or on the external face path between vertices Y and W
*/

#define TYPE_UNKNOWN            0

#define TYPE_VERTEX_VISITED		1

#define EDGE_DFSCHILD           1
#define EDGE_FORWARD            2
#define EDGE_BACK               3
#define EDGE_DFSPARENT          4

#define VERTEX_HIGH_RXW         6
#define VERTEX_LOW_RXW          7
#define VERTEX_HIGH_RYW         8
#define VERTEX_LOW_RYW          9

/* Data members needed by vertices and edges

   Vertices
        v: Carries original vertex number (same as array index)
                DFSNumber then uses it to store DFI.
                SortVertices then restores original vertex numbers when vertices
                are put in DFI order (i.e. not same as array index)
        visited: helps detect vertex visitation during various algorithms
                such as Walkup
        link: array indices that 'point' to the start and end arcs of the adjacency list
        type: Used by Kuratowski subgraph isolator to classify vertices when
                searching for certain paths in a biconnected component.
        flags: Lowest 16 bits a reserved for future expansion of the library.
               Next higher 16 bits can be safely used by consuming applications.
               Currently, no flag bits are used for vertices.

   Edges
        v: The edge record for (u,v) will be in u's list and store the index of
                the neighbour v. Starts out being original vertex number, but
                SortVertices renumbers to DFI so we get constant time access.
        visited: helps detect edge visitation, e.g. during the initial depth
                        first search, during a face reading, and during
                        Kuratowski subgraph isolation
        link: Linkages to other edges in an adjacency list.
        type: Used by DFSNumber to classify edges as DFSCHILD, DFSPARENT,
                FORWARD, BACK. See macro definitions above.
        flags: Lowest 16 bits a reserved for future expansion of the library.
               Next higher 16 bits can be safely used by consuming applications.
               The library uses bits 0 and 1 to indicate the INONLY and OUTONLY
               arcs of a directed edge.
               The planar embedder uses bit 2 on a DFSCHILD edge record of the
               root edge of a bicomp to indicate inverted orientation.
*/

typedef struct
{
     int  v;
     int  visited;
     int  link[2];
     int  type;
     int  flags;
} graphNode;

typedef graphNode * graphNodeP;

#define EDGEFLAG_INVERTED 4
#define GET_EDGEFLAG_INVERTED(theGraph, e) (theGraph->G[e].flags & EDGEFLAG_INVERTED)
#define SET_EDGEFLAG_INVERTED(theGraph, e) (theGraph->G[e].flags |= EDGEFLAG_INVERTED)
#define CLEAR_EDGEFLAG_INVERTED(theGraph, e) (theGraph->G[e].flags &= (~EDGEFLAG_INVERTED))

/* Data members needed by vertices
        DFSParent: The DFI of the DFS tree parent of this vertex
        leastAncestor: min(DFI of neighbors connected by backedge)
        Lowpoint: min(leastAncestor, min(Lowpoint of DFS Children))
        adjacentTo: Used by the embedder; during walk-up, each vertex that is
                directly adjacent via a back edge to the vertex currently
                being embedded will have the forward edge's index stored in
                this field.  During walkdown, each vertex whose AdjacentTo
                field is set will cause a back edge to be embedded.
        pertinentBicompList: used by Walkup to store a list of child bicomps of
                a vertex descendant of the current vertex that are pertinent
                and must be merged by the Walkdown in order to embed the cycle
                edges of the current vertex.  In this implementation,
                externally active pertinent child bicomps are placed at the end
                of the list as an easy way to make sure all internally active
                bicomps are processed first.
        separatedDFSChildList: contains list DFS children of this vertex in
                non-descending order by Lowpoint (sorted in linear time).
                When merging bicomp rooted by edge (r, c) into vertex v (i.e.
                merging root copy r with parent copy v), the vertex c is
                removed from the separatedDFSChildList of v.
                A vertex's status-- inactive, internally active, externally
                active-- is determined by the lesser of its leastAncestor and
                the least lowpoint from among only those DFS children that
                aren't in the same bicomp with the vertex.
        fwdArcList: at the start of embedding, the back edges from a vertex
                to its DFS descendants are separated from the main adjacency
                list and placed in a circular list until they are embedded.
                This member indicates a node in that list.
*/

typedef struct
{
        int DFSParent, leastAncestor, Lowpoint, adjacentTo;
        int pertinentBicompList, separatedDFSChildList, fwdArcList;
} vertexRec;

typedef vertexRec * vertexRecP;

/* This structure defines a pair of links used by each vertex and root copy
    to more efficiently traverse the external face.
    These also help in the creation of a planarity tester that does not need
    to embed the edges, which would be more efficient when one only needs to
    know whether any of a give set of graphs is planar without justifying
    the result with a combinatorial embedding. */

typedef struct
{
    int vertex[2];
    int inversionFlag;
} extFaceLinkRec;

typedef extFaceLinkRec * extFaceLinkRecP;

/* Flags for graph:
        FLAGS_DFSNUMBERED is set if DFSNumber() has succeeded for the graph
        FLAGS_SORTEDBYDFI records whether the graph is in original vertex
                order or sorted by depth first index.  Successive calls to
                SortVertices() toggle this bit.
        FLAGS_OBSTRUCTIONFOUND is set by gp_Embed() if an embedding obstruction
                was isolated in the graph returned.  It is cleared by gp_Embed()
                if an obstruction was not found.  The flag is used by
                gp_TestEmbedResultIntegrity() to decide what integrity tests to run.
*/

#define FLAGS_DFSNUMBERED       1
#define FLAGS_SORTEDBYDFI       2
#define FLAGS_OBSTRUCTIONFOUND  4

/* Variables needed in embedding by Kuratowski subgraph isolator:
        minorType: the type of planarity obstruction found.
        v: the current vertex I
        r: the root of the bicomp on which the Walkdown failed
        x,y: stopping vertices on bicomp rooted by r
        w: pertinent vertex on ext. face path below x and y
        px, py: attachment points of x-y path,
        z: Unused except in minors D and E (not needed in A, B, C).

        ux,dx: endpoints of unembedded edge that helps connext x with
                ancestor of v
        uy,dy: endpoints of unembedded edge that helps connext y with
                ancestor of v
        dw: descendant endpoint in unembedded edge to v
        uz,dz: endpoints of unembedded edge that helps connext z with
                ancestor of v (for minors B and E, not A, C, D).
*/

typedef struct
{
    int minorType;
    int v, r, x, y, w, px, py, z;
    int ux, dx, uy, dy, dw, uz, dz;
} isolatorContext;

typedef isolatorContext * isolatorContextP;

#define MINORTYPE_A         1
#define MINORTYPE_B         2
#define MINORTYPE_C         4
#define MINORTYPE_D         8
#define MINORTYPE_E         16
#define MINORTYPE_E1        32
#define MINORTYPE_E2        64
#define MINORTYPE_E3        128
#define MINORTYPE_E4        256

#define MINORTYPE_E5        512
#define MINORTYPE_E6        1024
#define MINORTYPE_E7        2048

/* Container for graph functions
        G: Vertices stored at 0 to n-1, second vertex buffer at n to 2n-1,
                edges at 2n and above
        V: Additional information about vertices
        N: Number of vertices
        M: Number of edges
        edgeOffset: always 2*N; location of the start of edge records in G
        arcCapacity: the maximum number of edge records allowed in G
        edgeHoles: free locations where edges have been deleted
        theStack: Used by various graph routines needing a stack
        internalFlags: Additional state information about the graph
        embedFlags: controls type of embedding (e.g. planar)

        IC: contains additional useful variables for Kuratowski subgraph isolation.
        BicompLists: storage space for pertinent bicomp lists that develop
                        during embedding
        DFSChildLists: storage space for separated DFS child lists that
                        develop during embedding
        buckets: Used to help bucket sort the separatedDFSChildList elements
                    of all vertices (see _CreateSortedSeparatedDFSChildLists())
        bin: Used to help bucket sort the separatedDFSChildList elements
                    of all vertices (see _CreateSortedSeparatedDFSChildLists())

        extensions: a list of extension data structures
        functions: a table of function pointers that can be overloaded to provide
                   extension behaviors to the graph
*/

typedef struct
{
        graphNodeP G;
        vertexRecP V;
        int N, M, edgeOffset, arcCapacity;
        stackP edgeHoles;
        stackP theStack;
        int internalFlags, embedFlags;

        isolatorContext IC;
        listCollectionP BicompLists, DFSChildLists;
        int *buckets;
        listCollectionP bin;
        extFaceLinkRecP extFace;

        graphExtensionP extensions;
        graphFunctionTable functions;

} baseGraphStructure;

typedef baseGraphStructure * graphP;

/********************************************************************
 _VertexActiveStatus()
 Tells whether a vertex is externally active, internally active
 or inactive.
 ********************************************************************/

#define _VertexActiveStatus(theGraph, theVertex, I) \
        (EXTERNALLYACTIVE(theGraph, theVertex, I) \
         ? VAS_EXTERNAL \
         : PERTINENT(theGraph, theVertex) \
           ? VAS_INTERNAL \
           : VAS_INACTIVE)

/********************************************************************
 PERTINENT()
 A vertex is pertinent in a partially processed graph if there is an
 unprocessed back edge between the vertex I whose edges are currently
 being processed and either the vertex or a DFS descendant D of the
 vertex not in the same bicomp as the vertex.

 The vertex is either directly adjacent to I by an unembedded back edge
 or there is an unembedded back edge (I, D) and the vertex is a cut
 vertex in the partially processed graph along the DFS tree path from
 D to I.

 Pertinence is a dynamic property that can change for a vertex after
 each edge addition.  In other words, a vertex can become non-pertinent
 during step I as more back edges to I are embedded.
 ********************************************************************/

#define PERTINENT(theGraph, theVertex) \
        (theGraph->V[theVertex].adjacentTo != NIL || \
         theGraph->V[theVertex].pertinentBicompList != NIL)

/********************************************************************
 FUTUREPERTINENT()
 A vertex is future-pertinent in a partially processed graph if
 there is an unprocessed back edge between a DFS ancestor A of the
 vertex I whose edges are currently being processed and either the
 vertex or a DFS descendant D of the vertex not in the same bicomp
 as the vertex.

 The vertex is either directly adjacent to A by an unembedded back edge
 or there is an unembedded back edge (A, D) and the vertex is a cut
 vertex in the partially processed graph along the DFS tree path from
 D to A.

 If no more edges are added to the partially processed graph prior to
 processing the edges of A, then the vertex would be pertinent.
 The addition of edges to the partially processed graph can alter
 both the pertinence and future pertinence of a vertex.  For example,
 if the vertex is pertinent due to an unprocessed back edge (I, D1) and
 future pertinent due to an unprocessed back edge (A, D2), then the
 vertex may lose both its pertinence and future pertinence when edge
 (I, D1) is added if D2 is equal to or an ancestor of D1.

 Generally, pertinence and future pertinence are dynamic properties
 that can change for a vertex after each edge addition.
 ********************************************************************/

#define FUTUREPERTINENT(theGraph, theVertex, I) \
        (  theGraph->V[theVertex].leastAncestor < I || \
           (theGraph->V[theVertex].separatedDFSChildList != NIL && \
            theGraph->V[theGraph->V[theVertex].separatedDFSChildList].Lowpoint < I) )

/********************************************************************
 EXTERNALLYACTIVE()
 Tells whether a vertex is still externally active in step I.
 A vertex can become inactive during step I as edges are embedded.

 In planarity-related algorithms, external activity is the same as
 future pertinence.  A vertex must be kept on the external face of
 the partial embedding until its pertinence and future pertinence
 are resolved through edge additions.

 For outerplanarity-related algorithms, all vertices are always
 externally active, since they must always remain on the external face.
 ********************************************************************/

#define EXTERNALLYACTIVE(theGraph, theVertex, I) \
        ( ( theGraph->embedFlags & EMBEDFLAGS_OUTERPLANAR) || \
          FUTUREPERTINENT(theGraph, theVertex, I) )

#ifdef __cplusplus
}
#endif

#endif

