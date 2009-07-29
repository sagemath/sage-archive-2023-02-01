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

#define GRAPHSTRUCTURE_C

#include <stdlib.h>

#include "graph.h"

/********************************************************************
 Private functions, except exported within library
 ********************************************************************/

void _InitGraphNode(graphP theGraph, int I);
void _ClearIsolatorContext(graphP theGraph);
void _FillVisitedFlags(graphP theGraph, int FillValue);
void _FillVisitedFlagsInBicomp(graphP theGraph, int BicompRoot, int FillValue);
void _FillVisitedFlagsInOtherBicomps(graphP theGraph, int BicompRoot, int FillValue);
void _SetVertexTypeInBicomp(graphP theGraph, int BicompRoot, int theType);

/********************************************************************
 Private functions.
 ********************************************************************/

void _ClearGraph(graphP theGraph);

void _InitVertexRec(graphP theGraph, int I);

int  _GetRandomNumber(int NMin, int NMax);

void _AddArc(graphP theGraph, int u, int v, int arcPos, int link);
void _HideArc(graphP theGraph, int arcPos);

/********************************************************************
 gp_New()
 Constructor for graph object.
 Can create two graphs if restricted to no dynamic memory.
 ********************************************************************/

graphP gp_New()
{
graphP theGraph = (graphP) malloc(sizeof(BM_graph));

     if (theGraph != NULL)
     {
         theGraph->G = NULL;
         theGraph->V = NULL;

         theGraph->BicompLists = NULL;
         theGraph->DFSChildLists = NULL;
         theGraph->theStack = NULL;

         theGraph->buckets = NULL;
         theGraph->bin = NULL;

         theGraph->extFace = NULL;

         _ClearGraph(theGraph);
     }

     return theGraph;
}

/********************************************************************
 gp_InitGraph()
 Allocates memory for graph and vertex records now that N is known.
 For G, we need N vertex nodes, N more vertex nodes for root copies,
        (2 * EDGE_LIMIT * N) edge records.
 For V, we need N vertex records.
 The BicompLists and DFSChildLists are of size N and start out empty.
 The stack, initially empty, is made big enough for a pair of integers per
        edge record, or 2 * 2 * EDGE_LIMIT * N.
 buckets and bin are both O(n) in size.  They are used by
        CreateSortedSeparatedDFSChildLists()
 Returns OK on success, NOTOK on all failures.
 ********************************************************************/

int  gp_InitGraph(graphP theGraph, int N)
{
int I;

     _ClearGraph(theGraph);

/* Allocate memory as described above */

     if ((theGraph->G = (graphNodeP) malloc((2+2*EDGE_LIMIT)*N*sizeof(graphNode))) == NULL ||
         (theGraph->V = (vertexRecP) malloc(N*sizeof(vertexRec))) == NULL ||
         (theGraph->BicompLists = LCNew(N)) == NULL ||
         (theGraph->DFSChildLists = LCNew(N)) == NULL ||
         (theGraph->theStack = sp_New(2 * 2 * EDGE_LIMIT * N)) == NULL ||
         (theGraph->buckets = (int *) malloc(N * sizeof(int))) == NULL ||
         (theGraph->bin = LCNew(N)) == NULL ||
         (theGraph->extFace = (extFaceLinkRecP) malloc(2*N*sizeof(extFaceLinkRec))) == NULL ||
         0)
     {
         _ClearGraph(theGraph);
         return NOTOK;
     }

/* Initialize memory */

     theGraph->N = N;

     for (I = 0; I < (2+2*EDGE_LIMIT)*N; I++)
          _InitGraphNode(theGraph, I);

     for (I = 0; I < N; I++)
          _InitVertexRec(theGraph, I);

     for (I = 0; I < 2*N; I++)
     {
         theGraph->extFace[I].link[0] = theGraph->extFace[I].link[1] = NIL;
         theGraph->extFace[I].inversionFlag = 0;
     }

     return OK;
}

/********************************************************************
 _InitGraphNode()
 Sets the fields in a single graph node structure to initial values
 ********************************************************************/

void _InitGraphNode(graphP theGraph, int I)
{
     theGraph->G[I].v =
     theGraph->G[I].link[0] =
     theGraph->G[I].link[1] = NIL;
     theGraph->G[I].visited = 0;
     theGraph->G[I].type = TYPE_UNKNOWN;
     theGraph->G[I].sign = 1;
}

/********************************************************************
 _InitVertexRec()
 Sets the fields in a single vertex record to initial values
 ********************************************************************/

void _InitVertexRec(graphP theGraph, int I)
{
     theGraph->V[I].leastAncestor =
     theGraph->V[I].Lowpoint = I;
     theGraph->V[I].DFSParent = NIL;
     theGraph->V[I].adjacentTo = NIL;
     theGraph->V[I].pertinentBicompList = NIL;
     theGraph->V[I].separatedDFSChildList = NIL;
     theGraph->V[I].fwdArcList = NIL;
}

/********************************************************************
 _ClearIsolatorContext()
 ********************************************************************/

void _ClearIsolatorContext(graphP theGraph)
{
isolatorContextP IC = &theGraph->IC;

     IC->minorType = 0;
     IC->v = IC->r = IC->x = IC->y = IC->w = IC->px = IC->py = IC->z =
     IC->ux = IC->dx = IC->uy = IC->dy = IC->dw = IC->uz = IC->dz = NIL;
}

/********************************************************************
 _FillVisitedFlags()
 ********************************************************************/

void _FillVisitedFlags(graphP theGraph, int FillValue)
{
int  i, limit;

     for (i=0, limit=2*(theGraph->N + theGraph->M); i<limit; i++)
          theGraph->G[i].visited = FillValue;

}

/********************************************************************
 _FillVisitedFlagsInBicomp()
 ********************************************************************/

void _FillVisitedFlagsInBicomp(graphP theGraph, int BicompRoot, int FillValue)
{
int  V, J;

     sp_ClearStack(theGraph->theStack);
     sp_Push(theGraph->theStack, BicompRoot);
     while (sp_NonEmpty(theGraph->theStack))
     {
          sp_Pop(theGraph->theStack, V);
          theGraph->G[V].visited = FillValue;

          J = theGraph->G[V].link[0];
          while (J >= 2*theGraph->N)
          {
             theGraph->G[J].visited = FillValue;

             if (theGraph->G[J].type == EDGE_DFSCHILD)
                 sp_Push(theGraph->theStack, theGraph->G[J].v);

             J = theGraph->G[J].link[0];
          }
     }
}

/********************************************************************
 _FillVisitedFlagsInOtherBicomps()
 Typically, we want to clear or set all visited flags in the graph
 (see _FillVisitedFlags).  However, in some algorithms this would be
 too costly, so it is necessary to clear or set the visited flags only
 in one bicomp (see _FillVisitedFlagsInBicomp), then do some processing
 that sets some of the flags then performs some tests.  If the tests
 are positive, then we can clear or set all the visited flags in the
 other bicomps (the processing may have set the visited flags in the
 one bicomp in a particular way that we want to retain, so we skip
 the given bicomp).
 ********************************************************************/

void _FillVisitedFlagsInOtherBicomps(graphP theGraph, int BicompRoot, int FillValue)
{
int  R, N;

     N = theGraph->N;

     for (R = N; R < 2*N; R++)
          if (theGraph->G[R].v != NIL && R != BicompRoot)
              _FillVisitedFlagsInBicomp(theGraph, R, FillValue);
}

/********************************************************************
 _SetVertexTypeInBicomp()
 ********************************************************************/

void _SetVertexTypeInBicomp(graphP theGraph, int BicompRoot, int theType)
{
int  V, J;

     sp_ClearStack(theGraph->theStack);
     sp_Push(theGraph->theStack, BicompRoot);
     while (sp_NonEmpty(theGraph->theStack))
     {
          sp_Pop(theGraph->theStack, V);
          theGraph->G[V].type = theType;

          J = theGraph->G[V].link[0];
          while (J >= 2*theGraph->N)
          {
             if (theGraph->G[J].type == EDGE_DFSCHILD)
                 sp_Push(theGraph->theStack, theGraph->G[J].v);

             J = theGraph->G[J].link[0];
          }
     }
}

/********************************************************************
 _ClearGraph()
 Clears all memory used by the graph, restoring it to the state it
 was in immediately after gp_New() created it.
 ********************************************************************/

void _ClearGraph(graphP theGraph)
{
     if (theGraph->G != NULL)
     {
          free(theGraph->G);
          theGraph->G = NULL;
     }
     if (theGraph->V != NULL)
     {
          free(theGraph->V);
          theGraph->V = NULL;
     }

     theGraph->N = theGraph->M = 0;
     theGraph->internalFlags = theGraph->embedFlags = 0;

     _ClearIsolatorContext(theGraph);

     LCFree(&theGraph->BicompLists);
     LCFree(&theGraph->DFSChildLists);

     sp_Free(&theGraph->theStack);

     if (theGraph->buckets != NULL)
     {
         free(theGraph->buckets);
         theGraph->buckets = NULL;
     }

     LCFree(&theGraph->bin);

     if (theGraph->extFace != NULL)
     {
         free(theGraph->extFace);
         theGraph->extFace = NULL;
     }
}

/********************************************************************
 gp_ReinitializeGraph()
 Reinitializes a graph, restoring it to the state it was in immediately
 gp_InitGraph() processed it.
 ********************************************************************/

void gp_ReinitializeGraph(graphP theGraph)
{
int  N = theGraph->N, I;

     theGraph->M = 0;
     theGraph->internalFlags = theGraph->embedFlags = 0;

     for (I = 0; I < (2+2*EDGE_LIMIT)*N; I++)
          _InitGraphNode(theGraph, I);

     for (I = 0; I < N; I++)
          _InitVertexRec(theGraph, I);

     for (I = 0; I < 2*N; I++)
     {
         theGraph->extFace[I].link[0] = theGraph->extFace[I].link[1] = NIL;
         theGraph->extFace[I].inversionFlag = 0;
     }

     _ClearIsolatorContext(theGraph);

     LCReset(theGraph->BicompLists);
     LCReset(theGraph->DFSChildLists);
     sp_ClearStack(theGraph->theStack);
     LCReset(theGraph->bin);
}

/********************************************************************
 gp_Free()
 Frees G and V, then the graph record.  Then sets your pointer to NULL
 (so you must pass the address of your pointer).
 ********************************************************************/

void gp_Free(graphP *pGraph)
{
     if (pGraph == NULL) return;
     if (*pGraph == NULL) return;

     _ClearGraph(*pGraph);

     free(*pGraph);
     *pGraph = NULL;
}

/********************************************************************
 gp_CopyGraph()
 Copies the content of the srcGraph into the dstGraph.  The dstGraph
 must have been previously initialized with the same number of
 vertices as the srcGraph (e.g. gp_InitGraph(dstGraph, srcGraph->N).

 Returns OK for success, NOTOK for failure.
 ********************************************************************/

int  gp_CopyGraph(graphP dstGraph, graphP srcGraph)
{
int I;

     /* Parameter checks */
     if (dstGraph == NULL || srcGraph == NULL)
         return NOTOK;

     if (dstGraph->N != srcGraph->N)
         return NOTOK;

     for (I = 0; I < (2+2*EDGE_LIMIT)*srcGraph->N; I++)
          dstGraph->G[I] = srcGraph->G[I];

     for (I = 0; I < srcGraph->N; I++)
          dstGraph->V[I] = srcGraph->V[I];

     for (I = 0; I < 2*srcGraph->N; I++)
     {
         dstGraph->extFace[I].link[0] = srcGraph->extFace[I].link[0];
         dstGraph->extFace[I].link[1] = srcGraph->extFace[I].link[1];
         dstGraph->extFace[I].inversionFlag = srcGraph->extFace[I].inversionFlag;
     }

     dstGraph->N = srcGraph->N;
     dstGraph->M = srcGraph->M;
     dstGraph->internalFlags = srcGraph->internalFlags;
     dstGraph->embedFlags = srcGraph->embedFlags;

     dstGraph->IC = srcGraph->IC;

     LCCopy(dstGraph->BicompLists, srcGraph->BicompLists);
     LCCopy(dstGraph->DFSChildLists, srcGraph->DFSChildLists);

     sp_Copy(dstGraph->theStack, srcGraph->theStack);

     return OK;
}

/********************************************************************
 gp_DupGraph()
 ********************************************************************/

graphP gp_DupGraph(graphP theGraph)
{
graphP result;

     if ((result = gp_New()) == NULL) return NULL;

     if (gp_InitGraph(result, theGraph->N) != OK ||
         gp_CopyGraph(result, theGraph) != OK)
     {
         gp_Free(&result);
         return NULL;
     }

     return result;
}

/********************************************************************
 gp_CreateRandomGraph()

 Creates a randomly generated graph.  First a tree is created by
 connecting each vertex to some successor.  Then a random number of
 additional random edges are added.  If an edge already exists, then
 we retry until a non-existent edge is picked.

 This function assumes the caller has already called srand().
 ********************************************************************/

int  gp_CreateRandomGraph(graphP theGraph)
{
int N, I, M, u, v;

     N = theGraph->N;

/* Generate a random tree; note that this method virtually guarantees
        that the graph will be renumbered, but it is linear time.
        Also, we are not generating the DFS tree but rather a tree
        that simply ensures the resulting random graph is connected. */

     for (I=1; I < N; I++)
          if (gp_AddEdge(theGraph, _GetRandomNumber(0, I-1), 0, I, 0) != OK)
                return NOTOK;

/* Generate a random number of additional edges
        (actually, leave open a small chance that no
        additional edges will be added). */

     M = _GetRandomNumber(7*N/8, EDGE_LIMIT*N);

     if (M > N*(N-1)/2) M = N*(N-1)/2;

     for (I=N-1; I<M; I++)
     {
          u = _GetRandomNumber(0, N-2);
          v = _GetRandomNumber(u+1, N-1);

          if (gp_IsNeighbor(theGraph, u, v))
              I--;
          else
          {
              if (gp_AddEdge(theGraph, u, 0, v, 0) != OK)
                  return NOTOK;
          }
     }

     return OK;
}

/********************************************************************
 _GetRandomNumber()
 This function generates a random number between NMin and NMax
 inclusive.  It assumes that the caller has called srand().
 It calls rand(), but before truncating to the proper range,
 it adds the high bits of the rand() result into the low bits.
 The result of this is that the randomness appearing in the
 truncated bits also has an affect on the non-truncated bits.
 I used the same technique to improve the spread of hashing functions
 in my Jan.98 Dr. Dobb's Journal article  "Resizable Data Structures".
 ********************************************************************/

int  _GetRandomNumber(int NMin, int NMax)
{
int  N = rand();

     if (NMax < NMin) return NMin;

     N += ((N&0xFFFF0000)>>16);
     N += ((N&0x0000FF00)>>8);
     N %= (NMax-NMin+1);
     return N+NMin;
}

/********************************************************************
 _getUnprocessedChild()
 Support routine for gp_Create RandomGraphEx(), this function
 obtains a child of the given vertex in the randomly generated
 tree that has not yet been processed.  NIL is returned if the
 given vertex has no unprocessed children

 ********************************************************************/

int _getUnprocessedChild(graphP theGraph, int parent)
{
int J = theGraph->G[parent].link[0];
int JTwin = gp_GetTwinArc(theGraph, J);
int child = theGraph->G[J].v;

    /* The tree edges were added to the link[0] side of each vertex,
        and we move processed tree edge records to the link[1] side,
        so if the immediate link[0] edge record is not a tree edge
        then we return NIL because the vertex has no remaining
        unprocessed children */

    if (theGraph->G[J].type == TYPE_UNKNOWN)
        return NIL;

    /* if the child has already been processed, then all children
        have been pushed to the link[1] side and we have just encountered
        the first child we processed, so there are no remaining
        unprocessed children */

    if (theGraph->G[J].visited)
        return NIL;

    /* We have found an edge leading to an unprocessed child, so
        we mark it as processed so that it doesn't get returned
        again in future iterations. */

    theGraph->G[J].visited = 1;
    theGraph->G[JTwin].visited = 1;

    /* Now we move the edge record in the parent vertex to the
        link[1] side of that vertex. */

    if (theGraph->G[J].link[0] != theGraph->G[J].link[1])
    {
        theGraph->G[parent].link[0] = theGraph->G[J].link[0];
        theGraph->G[theGraph->G[J].link[0]].link[1] = parent;
        theGraph->G[J].link[0] = parent;
        theGraph->G[J].link[1] = theGraph->G[parent].link[1];
        theGraph->G[theGraph->G[parent].link[1]].link[0] = J;
        theGraph->G[parent].link[1] = J;
    }

    /* Now we move the edge record in the child vertex to the
        link[1] of the child. */

    if (theGraph->G[J].link[0] != theGraph->G[J].link[1])
    {
        theGraph->G[theGraph->G[JTwin].link[0]].link[1] = theGraph->G[JTwin].link[1];
        theGraph->G[theGraph->G[JTwin].link[1]].link[0] = theGraph->G[JTwin].link[0];
        theGraph->G[JTwin].link[0] = child;
        theGraph->G[JTwin].link[1] = theGraph->G[child].link[1];
        theGraph->G[theGraph->G[child].link[1]].link[0] = JTwin;
        theGraph->G[child].link[1] = JTwin;
    }

    /* Now we set the child's parent and return the child. */

    theGraph->V[child].DFSParent = parent;

    return child;
}

/********************************************************************
 _hasUnprocessedChild()
 Support routine for gp_Create RandomGraphEx(), this function
 obtains a child of the given vertex in the randomly generated
 tree that has not yet been processed.  False (0) is returned
 unless the given vertex has an unprocessed child.
 ********************************************************************/

int _hasUnprocessedChild(graphP theGraph, int parent)
{
int J = theGraph->G[parent].link[0];

    if (theGraph->G[J].type == TYPE_UNKNOWN)
        return 0;

    if (theGraph->G[J].visited)
        return 0;

    return 1;
}

/********************************************************************
 gp_CreateRandomGraphEx()
 Given a graph structure with a pre-specified number of vertices N,
 this function creates a graph with the specified number of edges.

 If numEdges <= 3N-6, then the graph generated is planar.  If
 numEdges is larger, then a maximal planar graph is generated, then
 (numEdges - 3N + 6) additional random edges are added.

 This function assumes the caller has already called srand().
 ********************************************************************/

int  gp_CreateRandomGraphEx(graphP theGraph, int numEdges)
{
#define EDGE_TREE_RANDOMGEN (TYPE_UNKNOWN+1)

int N, I, arc, M, root, v, c, p, last, u, J, e;

     N = theGraph->N;

     if (numEdges > EDGE_LIMIT * N)
         numEdges = EDGE_LIMIT * N;

/* Generate a random tree. */

    for (I=1; I < N; I++)
    {
        v = _GetRandomNumber(0, I-1);
        if (gp_AddEdge(theGraph, v, 0, I, 0) != OK)
            return NOTOK;

        else
	    {
            arc = 2*N + 2*theGraph->M - 2;
		    theGraph->G[arc].type = EDGE_TREE_RANDOMGEN;
		    theGraph->G[gp_GetTwinArc(theGraph, arc)].type = EDGE_TREE_RANDOMGEN;
		    theGraph->G[arc].visited = 0;
		    theGraph->G[gp_GetTwinArc(theGraph, arc)].visited = 0;
	    }
    }

/* Add edges up to the limit or until the graph is maximal planar. */

    M = numEdges <= 3*N - 6 ? numEdges : 3*N - 6;

    root = 0;
    v = last = _getUnprocessedChild(theGraph, root);

    while (v != root && theGraph->M < M)
    {
	     c = _getUnprocessedChild(theGraph, v);

	     if (c != NIL)
	     {
             if (last != v)
             {
		         if (gp_AddEdge(theGraph, last, 1, c, 1) != OK)
			         return NOTOK;
             }

		     if (gp_AddEdge(theGraph, root, 1, c, 1) != OK)
			     return NOTOK;

		     v = last = c;
	     }

	     else
	     {
		     p = theGraph->V[v].DFSParent;
		     while (p != NIL && (c = _getUnprocessedChild(theGraph, p)) == NIL)
		     {
			     v = p;
			     p = theGraph->V[v].DFSParent;
			     if (p != NIL && p != root)
			     {
				     if (gp_AddEdge(theGraph, last, 1, p, 1) != OK)
					     return NOTOK;
			     }
		     }

		     if (p != NIL)
		     {
                 if (p == root)
                 {
                     if (gp_AddEdge(theGraph, v, 1, c, 1) != OK)
				         return NOTOK;

                     if (v != last)
                     {
			             if (gp_AddEdge(theGraph, last, 1, c, 1) != OK)
				             return NOTOK;
                     }
                 }
                 else
                 {
			         if (gp_AddEdge(theGraph, last, 1, c, 1) != OK)
				         return NOTOK;
                 }

                 if (p != root)
                 {
			        if (gp_AddEdge(theGraph, root, 1, c, 1) != OK)
				         return NOTOK;
                    last = c;
                 }

			     v = c;
		     }
	     }
    }

/* Add additional edges if the limit has not yet been reached. */

    while (theGraph->M < numEdges)
    {
        u = _GetRandomNumber(0, N-1);
        v = _GetRandomNumber(0, N-1);

        if (u != v && !gp_IsNeighbor(theGraph, u, v))
            if (gp_AddEdge(theGraph, u, 0, v, 0) != OK)
                return NOTOK;
    }

/* Clear the edge types back to 'unknown' */

    for (e = 0; e < numEdges; e++)
    {
        J = 2*N + 2*e;
        theGraph->G[J].type = TYPE_UNKNOWN;
        theGraph->G[gp_GetTwinArc(theGraph, J)].type = TYPE_UNKNOWN;
        theGraph->G[J].visited = 0;
        theGraph->G[gp_GetTwinArc(theGraph, J)].visited = 0;
    }

/* Put all DFSParent indicators back to NIL */

    for (I = 0; I < N; I++)
        theGraph->V[I].DFSParent = NIL;

    return OK;

#undef EDGE_TREE_RANDOMGEN
}

/********************************************************************
 gp_IsNeighbor()
 Checks whether v is already in u's adjacency list.
 Returns 1 for yes, 0 for no.
 ********************************************************************/

int  gp_IsNeighbor(graphP theGraph, int u, int v)
{
int  J;

     J = theGraph->G[u].link[0];
     while (J >= 2*theGraph->N)
     {
          if (theGraph->G[J].v == v) return 1;
          J = theGraph->G[J].link[0];
     }
     return 0;
}

/********************************************************************
 gp_GetVertexDegree()

 Counts the number of edge records in the adjacency list of a given
 vertex V.  The while loop condition is 2N or higher because our
 data structure keeps records at locations 0 to N-1 for vertices
 AND N to 2N-1 for copies of vertices.  So edge records are stored
 at locations 2N and above.
 ********************************************************************/

int  gp_GetVertexDegree(graphP theGraph, int v)
{
int  J, degree;

     if (theGraph==NULL || v==NIL) return 0;

     degree = 0;

     J = theGraph->G[v].link[0];
     while (J >= 2*theGraph->N)
     {
         degree++;
         J = theGraph->G[J].link[0];
     }

     return degree;
}

/********************************************************************
 _AddArc()
 This routine adds arc (u,v) to u's edge list, storing the record for
 v at position arcPos.  The record is either added to the link[0] or
 link[1] side of vertex u, depending on the link parameter.
 The links of a vertex record can be viewed as previous (link[0]) and
 next (link[1]) pointers.  Thus, an edge record is appended to u's
 list by hooking it to u.link[0], and it is prepended by hooking it
 to u.link[1].  The use of exclusive-or (i.e. 1^link) is simply to get
 the other link (if link is 0 then 1^link is 1, and vice versa).
 ********************************************************************/

void _AddArc(graphP theGraph, int u, int v, int arcPos, int link)
{
     theGraph->G[arcPos].v = v;
     if (theGraph->G[u].link[0] == NIL)
     {
         theGraph->G[u].link[0] = theGraph->G[u].link[1] = arcPos;
         theGraph->G[arcPos].link[0] = theGraph->G[arcPos].link[1] = u;
     }
     else
     {
     int u0 = theGraph->G[u].link[link];

         theGraph->G[arcPos].link[link] = u0;
         theGraph->G[arcPos].link[1^link] = u;

         theGraph->G[u].link[link] = arcPos;

         theGraph->G[u0].link[1^link] = arcPos;
     }
}

/********************************************************************
 gp_AddEdge()
 Adds the undirected edge (u,v) to the graph by placing edge records
 representing u into v's circular edge record list and v into u's
 circular edge record list.

 upos receives the location in G where the u record in v's list will be
        placed, and vpos is the location in G of the v record we placed in
 u's list.  These are used to initialize the short circuit links.

 ulink (0|1) indicates whether the edge record to v in u's list should
        become adjacent to u by its 0 or 1 link, i.e. u[ulink] == vpos.
 vlink (0|1) indicates whether the edge record to u in v's list should
        become adjacent to v by its 0 or 1 link, i.e. v[vlink] == upos.

 ********************************************************************/

int  gp_AddEdge(graphP theGraph, int u, int ulink, int v, int vlink)
{
int  upos, vpos;

     if (theGraph==NULL || u<0 || v<0 || u>=2*theGraph->N || v>=2*theGraph->N)
         return NOTOK;

     /* We enforce the edge limit */

     if (theGraph->M >= EDGE_LIMIT*theGraph->N)
         return NONPLANAR;

     vpos = 2*theGraph->N + 2*theGraph->M;
     upos = gp_GetTwinArc(theGraph, vpos);

     _AddArc(theGraph, u, v, vpos, ulink);
     _AddArc(theGraph, v, u, upos, vlink);

     theGraph->M++;
     return OK;
}

/********************************************************************
 _HideArc()
 This routine removes an arc from an edge list, but does not delete
 it from the data structure.  Many algorithms must temporarily remove
 an edge, perform some calculation, and eventually put the edge back.
 This routine supports that operation.

 The neighboring adjacency list nodes are cross-linked, but the
 link[0] and link[1] fields of the arc are retained so it can
 reinsert itself when _RestoreArc() is called.
 ********************************************************************/

void _HideArc(graphP theGraph, int arcPos)
{
int  link0, link1;

     link0 = theGraph->G[arcPos].link[0];
     link1 = theGraph->G[arcPos].link[1];
     if (link0==NIL || link1==NIL) return;

     theGraph->G[link0].link[1] = link1;
     theGraph->G[link1].link[0] = link0;
}

/********************************************************************
 _RestoreArc()
 This routine reinserts an arc into the edge list from which it
 was previously removed by _HideArc().

 The assumed processing model is that arcs will be restored in reverse
 of the order in which they were hidden, i.e. it is assumed that the
 hidden arcs will be pushed on a stack and the arcs will be popped
 from the stack for restoration.
 ********************************************************************/

void _RestoreArc(graphP theGraph, int arcPos)
{
int  link0, link1;

     link0 = theGraph->G[arcPos].link[0];
     link1 = theGraph->G[arcPos].link[1];
     if (link0==NIL || link1==NIL) return;

     theGraph->G[link0].link[1] = arcPos;
     theGraph->G[link1].link[0] = arcPos;
}

/********************************************************************
 gp_HideEdge()
 This routine removes an arc and its twin arc from its edge list,
 but does not delete them from the data structure.  Many algorithms must
 temporarily remove an edge, perform some calculation, and eventually
 put the edge back. This routine supports that operation.

 For each arc, the neighboring adjacency list nodes are cross-linked,
 but the link[0] and link[1] fields of the arc are retained so it can
 be reinserted by calling gp_RestoreEdge().
 ********************************************************************/

void gp_HideEdge(graphP theGraph, int arcPos)
{
     _HideArc(theGraph, arcPos);
     _HideArc(theGraph, gp_GetTwinArc(theGraph, arcPos));
}

/********************************************************************
 gp_RestoreEdge()
 This routine reinserts an arc and its twin arc into the edge list
 from which it was previously removed by gp_HideEdge().

 The assumed processing model is that edges will be restored in
 reverse of the order in which they were hidden, i.e. it is assumed
 that the hidden edges will be pushed on a stack and the edges will
 be popped from the stack for restoration.

 Note: Since both arcs of an edge are restored, only one arc need
        be pushed on the stack for restoration.  This routine
        restores the two arcs in the opposite order from the order
        in which they are hidden by gp_HideEdge().
 ********************************************************************/

void gp_RestoreEdge(graphP theGraph, int arcPos)
{
     _RestoreArc(theGraph, gp_GetTwinArc(theGraph, arcPos));
     _RestoreArc(theGraph, arcPos);
}

/****************************************************************************
 gp_DeleteEdge()

 This function deletes the given edge record J and its twin, reducing the
 number of edges M in the graph.
 Before the Jth record is deleted, its link[nextLink] is collected as the
 return result.  This is useful because it is the 'next' edge record in the
 adjacency list of a vertex, which is otherwise hard to obtain from record
 J once it is deleted.
 ****************************************************************************/

int  gp_DeleteEdge(graphP theGraph, int J, int nextLink)
{
int  JTwin = gp_GetTwinArc(theGraph, J);
int  N = theGraph->N, M = theGraph->M;
int  nextArc, JPos, MPos, i;

/* Calculate the nextArc after J so that, when J is deleted, the return result
        informs a calling loop of the next edge to be processed. */

     nextArc = theGraph->G[J].link[nextLink];

/* Delete the edge records J and JTwin. */

     theGraph->G[theGraph->G[J].link[0]].link[1] = theGraph->G[J].link[1];
     theGraph->G[theGraph->G[J].link[1]].link[0] = theGraph->G[J].link[0];
     theGraph->G[theGraph->G[JTwin].link[0]].link[1] = theGraph->G[JTwin].link[1];
     theGraph->G[theGraph->G[JTwin].link[1]].link[0] = theGraph->G[JTwin].link[0];

/* If records J and JTwin are not the last in the edge record array, then
        we move the last two edge records to replace J and JTwin.
        Also, if nextArc is moved in this process (i.e. if it is part of
        the last edge that is moved to replace J and JTwin), then we
        change the nextArc variable so that it continues to indicate the
        next edge record after J in the adjacency list containing J. */

     JPos = (J < JTwin ? J : JTwin) - 2*N;
     MPos = 2*(M-1);
     if (JPos < MPos)
     {
         for (i=0; i<=1; i++, JPos++, MPos++)
         {
             if (nextArc == 2*N+MPos) nextArc = 2*N+JPos;
             theGraph->G[2*N+JPos] = theGraph->G[2*N+MPos];
             theGraph->G[theGraph->G[2*N+JPos].link[0]].link[1] = 2*N+JPos;
             theGraph->G[theGraph->G[2*N+JPos].link[1]].link[0] = 2*N+JPos;
         }
     }

/* Now we reduce the number of edges in the data structure, and then
        return the previously calculated successor of J. */

     theGraph->M--;
     return nextArc;
}

#ifndef SPEED_MACROS

/********************************************************************
 gp_GetTwinArc()
 This function returns the calculated twin arc of a given arc.
 If the arc location is even, then the successor is the twin.
 If the arc node is odd, then the predecessor is the twin.
 ********************************************************************/

int  gp_GetTwinArc(graphP theGraph, int Arc)
{
     return (Arc & 1) ? Arc-1 : Arc+1;
}

#endif
