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

#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include "appconst.h"
#include "listcoll.h"
#include "stack.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The EDGELIMIT expresses the maximum number of edges allowed as a constant
        factor of N, the number of vertices. We allow 3N edges, but this
        number can be safely set to a larger integer value. */

#define EDGE_LIMIT      3

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
        link: array indices that 'point' to the start and end of the edge list

        type: Used by Kuratowski subgraph isolator to classify vertices when
                searching for certain paths in a biconnected component.
        sign: Unused

   Edges
        v: The edge record for (u,v) will be in u's list and store the index of
                the neighbour v. Starts out being original vertex number, but
                SortVertices renumbers to DFI so we get constant time access.
        visited: helps detect edge visitation, e.g. during the initial depth
                        first search, during a face reading, and during
                        Kuratowski subgraph isolation
        link: Linkages to other edges in an adjacency list.
        type: Used by DFSNumber to classify edges as DFSCHILD, DFSPARENT,
                FORWARD, BACK or SHORTCIRCUIT. See macro definitions above.
        sign: Initialized to 1, may become -1 during embedding.
                For planar embedder, set to -1 on the DFSCHILD edge record of
                the root edge of a bicomp that is flipped.  Used to recover
                consistent vertex orientation in each bicomp.
*/

typedef struct
{
     int  v;
     int  visited;
     int  link[2];
     int  type;
     int  sign;
} graphNode;

typedef graphNode * graphNodeP;

/* Additional data members needed only by vertices
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
    the result. */

typedef struct
{
    int link[2];
    int inversionFlag;
} extFaceLinkRec;

typedef extFaceLinkRec * extFaceLinkRecP;

/* Flags for graph:
        FLAGS_DFSNUMBERED is set if DFSNumber() has succeeded for the graph
        FLAGS_SORTEDBYDFI records whether the graph is in original vertex
                order or sorted by depth first index.  Successive calls to
                SortVertices() toggle this bit.
*/

#define FLAGS_DFSNUMBERED       1
#define FLAGS_SORTEDBYDFI       2

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

#define FLAGS_MINOR_A         1
#define FLAGS_MINOR_B         2
#define FLAGS_MINOR_C         4
#define FLAGS_MINOR_D         8
#define FLAGS_MINOR_E         16
#define FLAGS_MINOR_E1        32
#define FLAGS_MINOR_E2        64
#define FLAGS_MINOR_E3        128
#define FLAGS_MINOR_E4        256

#define FLAGS_MINOR_E5        512
#define FLAGS_MINOR_E6        1024
#define FLAGS_MINOR_E7        2048

/* Container for graph functions
        G: Vertices stored at 0 to n-1, second vertex buffer at n to 2n-1,
                edges at 2n and above
        V: Additional information about vertices
        N: Number of vertices
        M: Number of edges
        internalFlags: Additional state information about the graph
        embedFlags: controls type of embedding (e.g. planar)
        IC: contains additional useful variables for Kuratowski subgraph isolation.
        BicompLists: storage space for pertinent bicomp lists that develop
                        during embedding
        DFSChildLists: storage space for separated DFS child lists that
                        develop during embedding
        theStack: Used by various routines needing a stack, including
                depth first search, lowpoint, Walkdown, OrientVerticesInBicomp,
                and MarkHighestXYPath in the Kuratowski subgraph isolator
        buckets: Used to help bucket sort the separatedDFSChildList elements
                    of all vertices (see _CreateSortedSeparatedDFSChildLists())
        bin: Used to help bucket sort the separatedDFSChildList elements
                    of all vertices (see _CreateSortedSeparatedDFSChildLists())
*/

typedef struct
{
        graphNodeP G;
        vertexRecP V;
        int N, M, internalFlags, embedFlags;
        isolatorContext IC;
        listCollectionP BicompLists, DFSChildLists;
        stackP theStack;
        int *buckets;
        listCollectionP bin;
        extFaceLinkRecP extFace;
} BM_graph;

typedef BM_graph * graphP;

/* Public functions for graphs */

graphP gp_New(void);
int    gp_InitGraph(graphP theGraph, int N);
void   gp_ReinitializeGraph(graphP theGraph);
int    gp_CopyGraph(graphP dstGraph, graphP srcGraph);
graphP gp_DupGraph(graphP theGraph);

int    gp_CreateRandomGraph(graphP theGraph);
int    gp_CreateRandomGraphEx(graphP theGraph, int numEdges);

void   gp_Free(graphP *pGraph);

int    gp_Read(graphP theGraph, char *FileName);
int    gp_Write(graphP theGraph, char *FileName, int Mode);

int    gp_IsNeighbor(graphP theGraph, int u, int v);
int    gp_GetVertexDegree(graphP theGraph, int v);

int    gp_AddEdge(graphP theGraph, int u, int ulink, int v, int vlink);
void   gp_HideEdge(graphP theGraph, int arcPos);
void   gp_RestoreEdge(graphP theGraph, int arcPos);
int    gp_DeleteEdge(graphP theGraph, int J, int nextLink);

int    gp_CreateDFSTree(graphP theGraph);
int    gp_SortVertices(graphP theGraph);
void   gp_LowpointAndLeastAncestor(graphP theGraph);

int    gp_Embed(graphP theGraph, int embedFlags);

int    gp_CheckEmbeddingIntegrity(graphP theGraph);
int    gp_CheckKuratowskiSubgraphIntegrity(graphP theGraph);

#ifndef SPEED_MACROS

int    gp_GetTwinArc(graphP theGraph, int Arc);

#else
/********************************************************************
 int  gp_GetTwinArc(graphP theGraph, int Arc);
 This macro function returns the calculated twin arc of a given arc.
 If the arc location is even, then the successor is the twin.
 If the arc node is odd, then the predecessor is the twin.

 Logically, we return (Arc & 1) ? Arc-1 : Arc+1
 ********************************************************************/

// The original, first definition appears to be the faster one
#define gp_GetTwinArc(theGraph, Arc) ((Arc) & 1) ? Arc-1 : Arc+1
//#define gp_GetTwinArc(theGraph, Arc) ((Arc)+1-(((Arc)&1)<<1))

#endif

/* Possible Mode parameter of gp_Write */

#define WRITE_ADJLIST   1
#define WRITE_ADJMATRIX 2
#define WRITE_DEBUGINFO 3

/* Possible Flags for gp_Embed */

#define EMBEDFLAGS_PLANAR       1

/* Private functions shared by modules in this implementation */

#define PERTINENT(theEmbedding, theVertex, I) \
        (theEmbedding->V[theVertex].adjacentTo != NIL || \
         theEmbedding->V[theVertex].pertinentBicompList != NIL ? 1 : 0)

#define EXTERNALLYACTIVE(theEmbedding, theVertex, I) \
        (theEmbedding->V[theVertex].leastAncestor < I \
         ? 1 \
         : theEmbedding->V[theVertex].separatedDFSChildList != NIL && \
           theEmbedding->V[theEmbedding->V[theVertex].separatedDFSChildList].Lowpoint < I \
           ? 1 \
           : 0)

#ifndef SPEED_MACROS

int  _VertexActiveStatus(graphP theEmbedding, int theVertex, int I);

#else

#define _VertexActiveStatus(theEmbedding, theVertex, I) \
        (EXTERNALLYACTIVE(theEmbedding, theVertex, I) \
         ? VAS_EXTERNAL \
         : PERTINENT(theEmbedding, theVertex, I) \
           ? VAS_INTERNAL \
           : VAS_INACTIVE)

#endif

#ifdef __cplusplus
}
#endif

#endif

